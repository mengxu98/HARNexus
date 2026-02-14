#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message


def check_dir(path: str) -> str:
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def stage_to_num(stage: str):
    try:
        return int(str(stage).replace("S", ""))
    except Exception:
        return None


def main():
    p = argparse.ArgumentParser(description="Compute atlas network statistics.")
    p.add_argument(
        "--no-exclude-tfs-in-genes",
        action="store_true",
        help="Do not remove edges where Target is a TF (default: remove such edges)",
    )
    args = p.parse_args()
    exclude_tfs_in_genes = not args.no_exclude_tfs_in_genes

    log_message("Loading network data...")
    res_dir = check_dir("results/networks/analysis/")
    network_file = os.path.join(res_dir, "network_data.csv")

    if not os.path.exists(network_file):
        log_message(f"Network file not found: {network_file}", message_type="error")
        sys.exit(1)

    network_data = pd.read_csv(network_file)
    required_cols = {"Region", "Stage", "CellType", "TF", "Target"}
    missing = sorted(required_cols - set(network_data.columns))
    if missing:
        log_message(
            f"Missing required columns: {', '.join(missing)}", message_type="error"
        )
        sys.exit(1)

    if exclude_tfs_in_genes:
        tfs = network_data["TF"].unique()
        n_before = len(network_data)
        network_data = network_data[~network_data["Target"].isin(tfs)].copy()
        n_removed = n_before - len(network_data)
        if n_removed > 0:
            log_message(f"Excluded {n_removed:,} edges where Target is a TF", message_type="info")

    n_edges = len(network_data)
    n_regions = network_data["Region"].nunique(dropna=True)
    n_stages = network_data["Stage"].nunique(dropna=True)
    n_celltypes = network_data["CellType"].nunique(dropna=True)

    log_message(f"Edges: {n_edges:,}")
    log_message(f"Regions: {n_regions:,}")
    log_message(f"Stages: {n_stages:,}")
    log_message(f"Cell types: {n_celltypes:,}")

    n_region_stage = network_data.drop_duplicates(["Region", "Stage"]).shape[0]
    n_region_celltype = network_data.drop_duplicates(["Region", "CellType"]).shape[0]
    n_stage_celltype = network_data.drop_duplicates(["Stage", "CellType"]).shape[0]
    n_region_stage_celltype = network_data.drop_duplicates(
        ["Region", "Stage", "CellType"]
    ).shape[0]

    log_message(f"Unique Region×Stage: {n_region_stage:,}")
    log_message(f"Unique Region×CellType: {n_region_celltype:,}")
    log_message(f"Unique Stage×CellType: {n_stage_celltype:,}")
    log_message(f"Unique Region×Stage×CellType: {n_region_stage_celltype:,}")

    log_message("Calculating network statistics...")
    stats = (
        network_data.groupby(["Region", "Stage", "CellType"], dropna=False)
        .agg(
            Edges_count=("Target", "size"),
            TFs_count=("TF", pd.Series.nunique),
            Genes_count=("Target", pd.Series.nunique),
        )
        .reset_index()
        .rename(columns={"CellType": "Cell_type"})
    )

    stats["Stage_num"] = stats["Stage"].map(stage_to_num)
    stats = stats.sort_values(
        ["Region", "Cell_type", "Stage_num", "Stage"], na_position="last"
    ).drop(columns=["Stage_num"])

    out_file = os.path.join(res_dir, "network_statistics.csv")
    stats.to_csv(out_file, index=False)
    log_message(f"Saved network statistics to {out_file}")


if __name__ == "__main__":
    main()

