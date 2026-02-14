#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract target genes (excluding TFs) for PFC cell types across different stages for Upset plot visualization.

Input: data/networks/csv/network_data.csv
Output:
- results/networks/analysis/pfc_celltype_stage_targets.csv
- results/networks/analysis/pfc_celltype_stage_targets.json
"""

from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from typing import List, Optional, Tuple

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message

__all__ = ["extract_pfc_targets", "run"]

NETWORK_DATA_CSV = "results/networks/analysis/network_data.csv"
OUTPUT_DIR = "results/networks/analysis"
OUTPUT_CSV = os.path.join(OUTPUT_DIR, "pfc_celltype_stage_targets.csv")
OUTPUT_JSON = os.path.join(OUTPUT_DIR, "pfc_celltype_stage_targets.json")

PFC_REGION_NAMES = ["Prefrontal cortex"]


def extract_pfc_targets(
    network_csv: str = NETWORK_DATA_CSV,
    output_csv: str = OUTPUT_CSV,
    output_json: str = OUTPUT_JSON,
    output_dir: str = OUTPUT_DIR,
    pairs: Optional[List[Tuple[str, str]]] = None,
) -> dict:
    """
    Extract target genes (excluding TFs) for each cell type across different stages in PFC region.

    Parameters:
        network_csv: Path to network data CSV file
        output_csv: Path to output CSV file
        output_json: Path to output JSON file
        output_dir: Output directory

    Returns:
        dict: {cell_type: {stage: [target_genes]}}
    """
    os.makedirs(output_dir, exist_ok=True)

    # Read network data
    log_message(f"Reading network data: {network_csv}", message_type="info")
    if not os.path.exists(network_csv):
        raise FileNotFoundError(f"Network data file not found: {network_csv}")

    network_data = pd.read_csv(network_csv)
    log_message(f"Loaded {len(network_data):,} network edges", message_type="info")

    # Check required columns
    required_cols = ["Region", "CellType", "TF", "Target", "Stage"]
    missing_cols = [col for col in required_cols if col not in network_data.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Filter by pairs if specified, otherwise use PFC region names
    if pairs is not None:
        # Filter by specific (Region, CellType) pairs
        pair_mask = pd.Series([False] * len(network_data))
        for region, celltype in pairs:
            mask = (network_data["Region"] == region) & (network_data["CellType"] == celltype)
            pair_mask = pair_mask | mask
        filtered_data = network_data[pair_mask].copy()
        log_message(
            f"Filtered by {len(pairs)} pairs: {len(filtered_data):,} edges (from {len(network_data):,} total)",
            message_type="info",
        )
    else:
        # Filter PFC region data (original behavior)
        pfc_mask = network_data["Region"].isin(PFC_REGION_NAMES)
        filtered_data = network_data[pfc_mask].copy()
        log_message(
            f"Filtered PFC region data: {len(filtered_data):,} edges (from {len(network_data):,} total)",
            message_type="info",
        )

    if len(filtered_data) == 0:
        log_message("Warning: No data found after filtering", message_type="warning")
        # Check actual Region values in data
        unique_regions = network_data["Region"].unique()
        log_message(f"Region values in data: {unique_regions}", message_type="info")
        return {}

    # Extract all TF list (from TF column)
    all_tfs = set(filtered_data["TF"].unique())
    log_message(f"Extracted {len(all_tfs):,} TFs", message_type="info")

    # Get all stages and cell types (extracted from actual data)
    all_stages = sorted(filtered_data["Stage"].unique())
    all_cell_types = sorted(filtered_data["CellType"].unique())
    log_message(
        f"Found {len(all_stages)} stages: {', '.join(all_stages)}", message_type="info"
    )
    log_message(
        f"Found {len(all_cell_types)} cell types: {', '.join(all_cell_types)}",
        message_type="info",
    )

    # Extract target genes (excluding TFs) for each cell type and stage
    # Use cell types from actual data, not hardcoded list
    celltype_stage_targets = defaultdict(lambda: defaultdict(set))
    csv_rows = []

    for cell_type in all_cell_types:
        cell_data = filtered_data[filtered_data["CellType"] == cell_type]
        if len(cell_data) == 0:
            log_message(
                f"Warning: Cell type '{cell_type}' has no data", message_type="warning"
            )
            continue

        for stage in all_stages:
            stage_data = cell_data[cell_data["Stage"] == stage]
            if len(stage_data) == 0:
                continue

            # Extract Target genes
            targets = set(stage_data["Target"].unique())

            # Exclude TFs
            non_tf_targets = targets - all_tfs

            if len(non_tf_targets) > 0:
                celltype_stage_targets[cell_type][stage] = non_tf_targets
                # Prepare data for CSV
                for target in sorted(non_tf_targets):
                    csv_rows.append(
                        {
                            "CellType": cell_type,
                            "Stage": stage,
                            "Target": target,
                        }
                    )

                log_message(
                    f"  {cell_type} - {stage}: {len(non_tf_targets)} non-TF target genes "
                    f"(total targets: {len(targets)}, TFs: {len(targets & all_tfs)})",
                    message_type="info",
                )

    # Convert to regular dict and sort
    result_dict = {}
    for cell_type in sorted(celltype_stage_targets.keys()):
        result_dict[cell_type] = {}
        for stage in sorted(celltype_stage_targets[cell_type].keys()):
            result_dict[cell_type][stage] = sorted(
                celltype_stage_targets[cell_type][stage]
            )

    # Save CSV file
    if csv_rows:
        csv_df = pd.DataFrame(csv_rows)
        csv_df.to_csv(output_csv, index=False)
        log_message(
            f"Saved CSV file: {output_csv} ({len(csv_rows):,} rows)", message_type="success"
        )
    else:
        log_message("Warning: No data to save", message_type="warning")

    # Save JSON file (for R to read)
    with open(output_json, "w", encoding="utf-8") as f:
        json.dump(result_dict, f, indent=2, ensure_ascii=False)
    log_message(f"Saved JSON file: {output_json}", message_type="success")

    # Summary statistics
    total_targets = sum(
        len(targets) for ct_dict in result_dict.values() for targets in ct_dict.values()
    )
    log_message(
        f"Total: {len(result_dict)} cell types, {total_targets:,} target gene records",
        message_type="success",
    )

    return result_dict


def run(
    network_csv: str = NETWORK_DATA_CSV,
    output_csv: str = OUTPUT_CSV,
    output_json: str = OUTPUT_JSON,
    output_dir: str = OUTPUT_DIR,
    pairs: Optional[List[Tuple[str, str]]] = None,
) -> None:
    """Run the analysis pipeline"""
    extract_pfc_targets(
        network_csv=network_csv,
        output_csv=output_csv,
        output_json=output_json,
        output_dir=output_dir,
        pairs=pairs,
    )


def main():
    # pairs = [
    #     ("Prefrontal cortex", "Excitatory neurons"),
    #     ("Prefrontal cortex", "Astrocytes")
    # ]
    
    # log_message(f"Processing {len(pairs)} (Region, CellType) pairs", message_type="info")
    # for region, celltype in pairs:
    #     log_message(f"  - {region}, {celltype}", message_type="info")

    # Default paths
    network_csv = NETWORK_DATA_CSV
    output_csv = OUTPUT_CSV
    output_json = OUTPUT_JSON
    output_dir = OUTPUT_DIR

    run(
        network_csv=network_csv,
        output_csv=output_csv,
        output_json=output_json,
        output_dir=output_dir,
        # pairs=pairs,
    )


if __name__ == "__main__":
    main()
