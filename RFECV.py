#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Project: GBM Gene Feature Selection via Ensemble RFE Voting
===============================================================================
"""

import logging
import json
import pickle
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold
from sklearn.base import clone
import xgboost as xgb

# :: Configure Logging
log_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"pipeline_{log_timestamp}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class ConfigManager:
    """Handles loading and validating configuration settings."""
    
    def __init__(self, config_path: str = "config.json"):
        self.script_dir = Path(__file__).parent.resolve()
        self.config_path = self.script_dir / config_path
        self.config = self._load_config()
        self._resolve_paths()
        
    def _load_config(self) -> Dict:
        if not self.config_path.exists():
            logger.warning(f"Config file not found at {self.config_path}. Creating default...")
            self._create_default_config()
            
        with open(self.config_path, 'r', encoding='utf-8') as f:
            return json.load(f)
            
    def _create_default_config(self):
        default = {
            "paths": {
                "input_data": "data/GBM_NMF_group4.csv",
                "cache_dir": "cache",
                "output_dir": "results",
                "output_prefix": "GBM_group4"
            },
            "params": {
                "iteration_list": [10, 20, 30, 40, 50],
                "cv_splits": 5,
                "scoring": "accuracy",
                "vote_threshold": 1.0
            },
            "models": {
                "rf": {"type": "RandomForest", "random_state": 42},
                "svm": {"type": "SVC", "kernel": "linear", "random_state": 42},
                "xgb": {"type": "XGB", "random_state": 42, "use_label_encoder": False}
            }
        }
        with open(self.config_path, 'w', encoding='utf-8') as f:
            json.dump(default, f, indent=4)
        logger.info(f"Default config created at {self.config_path}")
        
    def _resolve_paths(self):
        """Converts relative paths in config to absolute paths based on script location."""
        for key, value in self.config['paths'].items():
            if isinstance(value, str):
                abs_path = self.script_dir / value
                self.config['paths'][key] = abs_path
                
    def get(self, *keys):
        """Nested access to config values."""
        val = self.config
        for k in keys:
            val = val[k]
        return val


class EnsembleFeatureSelector:
    """
    Performs ensemble feature selection using RFE voting across multiple models.
    Supports multiple iteration counts for stability analysis.
    """
    
    def __init__(self,  pd.DataFrame, config: ConfigManager, n_iterations: int):
        self.data = data
        self.config = config
        self.n_iterations = n_iterations
        self.X = None
        self.y = None
        self.voting_matrix = None
        self.feature_names = None
        
    def prepare_data(self):
        """Separate features and target labels."""
        if 'cluster' not in self.data.columns:
            raise ValueError("Target column 'cluster' not found in data.")
        
        self.X = self.data.drop(columns=['cluster'])
        self.y = self.data['cluster']
        self.feature_names = self.X.columns.tolist()
        logger.info(f"Data prepared. Features: {self.X.shape[1]}, Samples: {self.X.shape[0]}")
        
    def _get_estimator(self, model_name: str):
        """Factory method to create model instances based on config."""
        model_cfg = self.config.get('models', model_name)
        model_type = model_cfg.get('type')
        params = {k: v for k, v in model_cfg.items() if k != 'type'}
        
        if model_type == 'RandomForest':
            return RandomForestClassifier(**params)
        elif model_type == 'SVC':
            return SVC(**params)
        elif model_type == 'XGB':
            return xgb.XGBClassifier(**params)
        else:
            raise ValueError(f"Unknown model type: {model_type}")

    def _run_single_model_rfe(self, model_name: str) -> np.ndarray:
        """
        Run iterative RFECV for a single model type.
        Returns a matrix of shape (n_iterations, n_features) with rankings.
        """
        logger.info(f"Starting {self.n_iterations} RFE iterations for model: {model_name}")
        rankings = []
        
        cv_splits = self.config.get('params', 'cv_splits')
        scoring = self.config.get('params', 'scoring')
        cv = StratifiedKFold(n_splits=cv_splits, shuffle=True)
        
        for i in range(1, self.n_iterations + 1):
            estimator = self._get_estimator(model_name)
            if hasattr(estimator, 'random_state'):
                estimator.set_params(random_state=i)
                
            rfecv = RFECV(
                estimator=estimator,
                cv=cv,
                scoring=scoring,
                n_jobs=-1
            )
            
            y_train = self.y
            if model_name == 'xgb':
                y_train = self.y - 1
            
            try:
                rfecv.fit(self.X, y_train)
                rankings.append(rfecv.ranking_)
            except Exception as e:
                logger.warning(f"Iteration {i} failed for {model_name}: {e}")
                rankings.append(np.ones(len(self.feature_names)) * 2)
                
        return np.array(rankings)

    def run_ensemble_selection(self) -> pd.DataFrame:
        """Execute feature selection across all configured models."""
        cache_dir = Path(self.config.get('paths', 'cache_dir'))
        cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Cache file name includes iteration count for differentiation
        cache_file = cache_dir / f"rfe_cache_iter{self.n_iterations}.pkl"
        
        # Try loading from cache
        if cache_file.exists():
            try:
                logger.info(f"Loading cached results from {cache_file}")
                with open(cache_file, 'rb') as f:
                    cache_data = pickle.load(f)
                self.voting_matrix = cache_data['voting_matrix']
                self.feature_names = cache_data['feature_names']
                return self._process_votes()
            except Exception as e:
                logger.warning(f"Cache loading failed: {e}. Recomputing...")
        
        # Compute fresh results
        if self.X is None:
            self.prepare_data()
            
        all_rankings = []
        model_names = list(self.config.config['models'].keys())
        
        for name in model_names:
            rankings = self._run_single_model_rfe(name)
            all_rankings.append(rankings)
            
        self.voting_matrix = np.vstack(all_rankings)
        logger.info(f"Total voting matrix shape: {self.voting_matrix.shape}")
        
        # Save to cache
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump({
                    'voting_matrix': self.voting_matrix,
                    'feature_names': self.feature_names
                }, f)
            logger.info(f"Results cached to {cache_file}")
        except Exception as e:
            logger.error(f"Failed to save cache: {e}")
            
        return self._process_votes()
    
    def _process_votes(self) -> pd.DataFrame:
        """Calculate vote counts and filter features based on threshold."""
        vote_counts = np.sum(self.voting_matrix == 1, axis=0)
        max_possible_votes = self.voting_matrix.shape[0]
        threshold_ratio = self.config.get('params', 'vote_threshold')
        threshold_count = int(max_possible_votes * threshold_ratio)
        
        logger.info(f"Total votes possible: {max_possible_votes}, Threshold: {threshold_count}")
        
        selected_indices = np.where(vote_counts >= threshold_count)[0]
        selected_features = [self.feature_names[i] for i in selected_indices]
        
        logger.info(f"Selected {len(selected_features)} features out of {len(self.feature_names)}")
        
        return self.data[selected_features + ['cluster']]

    def save_results(self, result_df: pd.DataFrame):
        """Save the final selected features and the list of names."""
        output_dir = Path(self.config.get('paths', 'output_dir'))
        output_dir.mkdir(parents=True, exist_ok=True)
        
        prefix = self.config.get('paths', 'output_prefix')
        threshold = int(self.config.get('params', 'vote_threshold') * 100)
        
        # File names include iteration count for differentiation
        file_name = f"{prefix}_iter{self.n_iterations}_selected_{threshold}pct.csv"
        result_path = output_dir / file_name
        result_df.to_csv(result_path, index=True)
        logger.info(f"Dataset saved to {result_path}")
        
        probe_list = [col for col in result_df.columns if col != 'cluster']
        probe_path = output_dir / f"{prefix}_iter{self.n_iterations}_probes_{threshold}pct.csv"
        pd.Series(probe_list).to_csv(probe_path, index=False, header=False)
        logger.info(f"Probe list saved to {probe_path}")
        
        return {
            'n_iterations': self.n_iterations,
            'n_selected': len(probe_list),
            'n_total_features': len(self.feature_names),
            'file_path': str(result_path)
        }


class MultiIterationRunner:
    """
    Orchestrates multiple feature selection runs with different iteration counts.
    Generates summary statistics across all iterations.
    """
    
    def __init__(self,  pd.DataFrame, config: ConfigManager):
        self.data = data
        self.config = config
        self.results_summary = []
        
    def run_all_iterations(self) -> pd.DataFrame:
        """Run feature selection for all configured iteration counts."""
        iteration_list = self.config.get('params', 'iteration_list')
        logger.info(f"Starting batch run for iterations: {iteration_list}")
        
        for n_iter in iteration_list:
            logger.info(f"{'='*60}")
            logger.info(f"Running iteration count: {n_iter}")
            logger.info(f"{'='*60}")
            
            selector = EnsembleFeatureSelector(self.data, self.config, n_iter)
            selected_data = selector.run_ensemble_selection()
            result_info = selector.save_results(selected_data)
            
            self.results_summary.append(result_info)
            
            # Optional: Add delay between runs to prevent resource exhaustion
            # time.sleep(1)
        
        # Generate summary report
        self._generate_summary_report()
        
        return pd.DataFrame(self.results_summary)
    
    def _generate_summary_report(self):
        """Generate and save a summary report of all iteration runs."""
        output_dir = Path(self.config.get('paths', 'output_dir'))
        output_dir.mkdir(parents=True, exist_ok=True)
        
        summary_df = pd.DataFrame(self.results_summary)
        summary_path = output_dir / "iteration_summary.csv"
        summary_df.to_csv(summary_path, index=False)
        logger.info(f"Summary report saved to {summary_path}")
        
        # Log summary statistics
        logger.info("\n" + "="*60)
        logger.info("ITERATION SUMMARY")
        logger.info("="*60)
        logger.info(summary_df.to_string(index=False))
        logger.info("="*60)
        
        # Generate stability analysis (optional)
        self._analyze_feature_stability()
    
    def _analyze_feature_stability(self):
        """
        Analyze how feature selection stabilizes across different iteration counts.
        """
        output_dir = Path(self.config.get('paths', 'output_dir'))
        output_dir.mkdir(parents=True, exist_ok=True)
        
        stability_data = []
        iteration_list = self.config.get('params', 'iteration_list')
        
        for result in self.results_summary:
            stability_data.append({
                'iterations': result['n_iterations'],
                'selected_features': result['n_selected'],
                'selection_ratio': result['n_selected'] / result['n_total_features']
            })
        
        stability_df = pd.DataFrame(stability_data)
        stability_path = output_dir / "feature_stability_analysis.csv"
        stability_df.to_csv(stability_path, index=False)
        logger.info(f"Stability analysis saved to {stability_path}")


def main():
    """Main execution pipeline."""
    logger.info("=== Starting Feature Selection Pipeline ===")
    logger.info(f"Script running from: {Path(__file__).parent.resolve()}")
    
    # 1. Load Configuration
    config = ConfigManager()
    
    # 2. Load Data
    input_path = Path(config.get('paths', 'input_data'))
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        logger.error("Please ensure your data is placed correctly or update config.json")
        return
        
    try:
        data = pd.read_csv(input_path, sep=',', index_col=0, header=0)
        logger.info(f"Data loaded: {data.shape}")
    except Exception as e:
        logger.error(f"Failed to load  {e}")
        return

    # 3. Initialize Multi-Iteration Runner
    runner = MultiIterationRunner(data, config)
    
    # 4. Run All Iterations
    summary = runner.run_all_iterations()
    
    # 5. Final Report
    logger.info("\n=== Pipeline Completed Successfully ===")
    logger.info(f"Total runs completed: {len(summary)}")
    logger.info(f"Results saved to: {config.get('paths', 'output_dir')}")

if __name__ == "__main__":
    main()