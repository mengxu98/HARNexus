#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Process HAR-CSN atlas network data from CSV files

Steps:
1. Load TF-HAR pairs data from results/har_tf/human/har_tf_pairs_scores.csv
2. Process all network CSV files from results/networks/har_csn_atlas/csv/
3. Extract brain region, stage, and cell type from filenames
4. Combine all networks into a unified network_data.csv
"""

import os
import sys
import glob
import pandas as pd
import re

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message

HAR_TF_PAIRS_FILE = "results/har_tf/human/har_tf_pairs_scores.csv"
NETWORK_CSV_DIR = "results/networks/har_csn_atlas/csv"
OUTPUT_DIR = "results/networks/analysis"
TF_HAR_OUTPUT = os.path.join("results/har_tf/human/tf_har_combinations.csv")
NETWORK_OUTPUT = os.path.join(OUTPUT_DIR, "network_data.csv")


def get_components(filename):
    basename = os.path.basename(filename).replace(".csv", "")
    parts = basename.split("_")

    stage_idx = None
    for i, part in enumerate(parts):
        if re.match(r"^S[0-9]+$", part):
            stage_idx = i
            break

    if stage_idx is not None:
        stage = parts[stage_idx]
        brain_region = " ".join(parts[:stage_idx])
        cell_type = " ".join(parts[stage_idx + 1 :])
    else:
        stage = parts[-2] if len(parts) >= 2 else parts[-1]
        brain_region = " ".join(parts[:-2]) if len(parts) > 2 else parts[0]
        cell_type = parts[-1]

    return {"brain_region": brain_region, "stage": stage, "cell_type": cell_type}


log_message("Loading TF-HAR pairs data...", message_type="info")

har_tf_pairs = pd.read_csv(HAR_TF_PAIRS_FILE)

tf_har_df = har_tf_pairs[["TF", "har"]].drop_duplicates()
tf_har_df.columns = ["TF", "HAR"]

log_message(f"TF-HAR combinations: {len(tf_har_df):,} pairs", message_type="info")
log_message(f"Unique TFs: {len(tf_har_df['TF'].unique()):,}", message_type="info")
log_message(f"Unique HARs: {len(tf_har_df['HAR'].unique()):,}", message_type="info")

os.makedirs(OUTPUT_DIR, exist_ok=True)

tf_har_df.to_csv(TF_HAR_OUTPUT, index=False)
log_message(f"Saved TF-HAR combinations to {TF_HAR_OUTPUT}", message_type="success")

csv_files = glob.glob(os.path.join(NETWORK_CSV_DIR, "*.csv"))
log_message(f"Found {len(csv_files):,} network CSV files", message_type="info")

log_message("Processing network files...", message_type="info")

all_data_list = []

for csv_file in csv_files:
    components = get_components(csv_file)
    brain_region = components["brain_region"]
    stage = components["stage"]
    cell_type = components["cell_type"]

    log_message(
        f"Processing {os.path.basename(csv_file)}: "
        f"{brain_region} - {stage} - {cell_type}",
        message_type="running",
    )

    try:
        network = pd.read_csv(csv_file)

        if network.empty:
            log_message(
                f"Warning: {os.path.basename(csv_file)} is empty. Skipping",
                message_type="warning",
            )
            continue

        required_cols = ["regulator", "target", "weight"]
        if not all(col in network.columns for col in required_cols):
            log_message(
                f"Warning: {os.path.basename(csv_file)} does not have expected columns"
                f"Expected: {required_cols}, Found: {list(network.columns)}. Skipping",
                message_type="warning",
            )
            continue

        network_df = pd.DataFrame(
            {
                "Region": brain_region,
                "Stage": stage,
                "CellType": cell_type,
                "TF": network["regulator"],
                "Target": network["target"],
                "Weight": network["weight"],
            }
        )

        all_data_list.append(network_df)

    except Exception as e:
        log_message(
            f"Error loading {os.path.basename(csv_file)}: {str(e)}",
            message_type="error",
        )
        continue

if all_data_list:
    network_df = pd.concat(all_data_list, ignore_index=True)
    log_message(f"Total network edges: {len(network_df):,}", message_type="success")

    network_df.to_csv(NETWORK_OUTPUT, index=False)
    log_message(f"Saved network data to {NETWORK_OUTPUT}", message_type="success")
else:
    log_message(
        "No network data was processed. Please check the CSV files",
        message_type="warning",
    )

log_message("Atlas processing complete!", message_type="success")
