#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Project: ERSR-GBM-Subtypes
Module:  mapping/main.py
Purpose: Gene-Probe ID Mapping Entry Point
===============================================================================
"""

import json
import logging
from pathlib import Path
from gene_probe_mapper import GeneProbeMapper

# :: Configure Logging ---------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("mapping.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# :: Load Configuration --------------------------------------------------------
def load_config(config_path: str = "mapping/config.json"):
    """Load configuration from JSON file."""
    script_dir = Path(__file__).parent.resolve()
    config_file = script_dir / config_path
    
    if not config_file.exists():
        logger.warning(f"Config not found: {config_file}. Creating default...")
        create_default_config(config_file)
    
    with open(config_file, 'r', encoding='utf-8') as f:
        return json.load(f)

def create_default_config(config_file: Path):
    """Create default configuration file."""
    default = {
        "paths": {
            "ers_genes": "data/ERSRgenes.csv",
            "id_map": "data/IDmap.csv",
            "output_dir": "results/mapping",
            "output_file": "map_result.csv"
        },
        "params": {
            "gene_column": "Gene",
            "probe_gene_column": "UCSC_RefGene_Name",
            "probe_id_column": "ID",
            "probe_chr_column": "CHR",
            "probe_island_column": "UCSC_CpG_Islands_Name",
            "probe_enhancer_column": "Enhancer"
        },
        "output": {
            "columns": ["gene", "ori_idx", "ID", "CHR", "UCSC_CpG_Islands_Name", "Enhancer"],
            "format": "csv"
        }
    }
    
    config_file.parent.mkdir(parents=True, exist_ok=True)
    with open(config_file, 'w', encoding='utf-8') as f:
        json.dump(default, f, indent=4, ensure_ascii=False)
    
    logger.info(f"Default config created: {config_file}")

# :: Main Execution -----------------------------------------------------------
def main():
    """Main execution pipeline."""
    logger.info("=== Starting Gene-Probe Mapping Pipeline ===")
    
    # Load config
    config = load_config()
    
    # Resolve paths (relative to project root)
    project_root = Path(__file__).parent.parent.resolve()
    for key in config['paths']:
        if not Path(config['paths'][key]).is_absolute():
            config['paths'][key] = str(project_root / config['paths'][key])
    
    # Initialize mapper
    mapper = GeneProbeMapper(config)
    
    # Run mapping
    mapper.load_data()
    mapper.build_mapping()
    mapper.save_results()
    
    logger.info("=== Mapping Pipeline Completed Successfully ===")
    logger.info(f"Results saved to: {config['paths']['output_dir']}")

if __name__ == "__main__":
    main()
