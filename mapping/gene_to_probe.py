#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Module: gene_to_probe.py
Purpose: Core logic for Gene-Probe ID mapping
===============================================================================
"""

import pandas as pd
import csv
import logging
from pathlib import Path
from typing import Dict, List

logger = logging.getLogger(__name__)

class GeneProbeMapper:
    """
    Maps gene names to probe IDs using annotation files.
    Handles multi-gene probes (semicolon-separated).
    """
    
    def __init__(self, config: dict):
        self.config = config
        self.ers_genes = None
        self.id_map = None
        self.mapping_result = {}
        
    def load_data(self):
        """Load ERS genes and ID mapping reference files."""
        logger.info("Step 1: Loading reference data...")
        
        # Load ERS genes
        ers_path = self.config['paths']['ers_genes']
        if not Path(ers_path).exists():
            raise FileNotFoundError(f"ERS genes file not found: {ers_path}")
        
        self.ers_genes = pd.read_csv(
            ers_path, 
            sep=',', 
            header=0, 
            usecols=[0]
        )
        logger.info(f"  Loaded ERS genes: {len(self.ers_genes)} genes")
        
        # Load ID map
        id_path = self.config['paths']['id_map']
        if not Path(id_path).exists():
            raise FileNotFoundError(f"ID map file not found: {id_path}")
        
        self.id_map = pd.read_csv(
            id_path, 
            sep=',', 
            header=0
        )
        logger.info(f"  Loaded ID map: {len(self.id_map)} probes")
        
    def build_mapping(self):
        """Build gene-to-probe mapping."""
        logger.info("Step 2: Building gene-probe mapping...")
        
        # Initialize result dictionary
        gene_col = self.config['params']['gene_column']
        for gene in self.ers_genes[gene_col]:
            self.mapping_result[gene] = []
        
        # Iterate through all probes
        probe_gene_col = self.config['params']['probe_gene_column']
        total_probes = len(self.id_map)
        
        for idx in range(total_probes):
            if idx % 50000 == 0:
                logger.info(f"  Processing probe {idx}/{total_probes}...")
            
            gene_names = self.id_map.loc[idx, probe_gene_col]
            
            # Skip NaN values
            if pd.isna(gene_names) or str(gene_names) == "nan":
                continue
            
            # Handle multi-gene probes (semicolon-separated)
            if ";" in str(gene_names):
                gene_list = str(gene_names).split(';')
            else:
                gene_list = [gene_names]
            
            # Map each gene
            for gene in gene_list:
                gene = gene.strip()  # Remove whitespace
                if gene in self.mapping_result:
                    # Build record
                    record = {
                        'gene': gene,
                        'ori_idx': idx + 1,  # 1-based index
                        'ID': self.id_map.loc[idx, self.config['params']['probe_id_column']],
                        'CHR': self.id_map.loc[idx, self.config['params']['probe_chr_column']],
                        'UCSC_CpG_Islands_Name': self.id_map.loc[idx, self.config['params']['probe_island_column']],
                        'Enhancer': self.id_map.loc[idx, self.config['params']['probe_enhancer_column']]
                    }
                    self.mapping_result[gene].append(record)
        
        # Log statistics
        total_mappings = sum(len(v) for v in self.mapping_result.values())
        genes_with_probes = sum(1 for v in self.mapping_result.values() if len(v) > 0)
        
        logger.info(f"  Total mappings: {total_mappings}")
        logger.info(f"  Genes with probes: {genes_with_probes}/{len(self.mapping_result)}")
        
    def save_results(self):
        """Save mapping results to CSV file."""
        logger.info("Step 3: Saving results...")
        
        output_dir = Path(self.config['paths']['output_dir'])
        output_dir.mkdir(parents=True, exist_ok=True)
        
        output_file = output_dir / self.config['paths']['output_file']
        columns = self.config['output']['columns']
        
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            
            # Write header
            writer.writerow(columns)
            
            # Write data
            row_count = 0
            for gene, records in self.mapping_result.items():
                for record in records:
                    row = [record[col] for col in columns]
                    writer.writerow(row)
                    row_count += 1
            
        logger.info(f"  Saved {row_count} records to {output_file}")
        
        # Also save summary statistics
        summary_file = output_dir / "mapping_summary.json"
        import json
        summary = {
            'total_ers_genes': len(self.ers_genes),
            'total_probes_scanned': len(self.id_map),
            'total_mappings': row_count,
            'genes_with_probes': sum(1 for v in self.mapping_result.values() if len(v) > 0),
            'genes_without_probes': sum(1 for v in self.mapping_result.values() if len(v) == 0)
        }
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=4, ensure_ascii=False)
        
        logger.info(f"  Saved summary to {summary_file}")
