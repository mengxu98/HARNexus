#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract target genes of all-stages TFs from har_csn_atlas, save subnetwork and gene table.

Uses all_stages_TFs.csv (Region, CellType, TF, stage_list, ...) to filter edges from
har_csn_atlas/csv: only edges where (Region, CellType) and TF match. No all_stages_TF_target_edges.csv.

Inputs: all_stages_TFs.csv; har_csn_atlas/csv/*.csv
Outputs: all_stages_subnetwork_edges.csv, all_stages_network_genes.csv
"""

from __future__ import annotations

import importlib.util
import os
import sys
from collections import defaultdict
from typing import List, Optional

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message

__all__ = [
    "load_subnetwork_from_tfs",
    "extract_genes_from_edges",
    "run",
]

TFS_CSV = "results/networks/analysis/all_stages_TFs.csv"
SUBNETWORK_CSV = "results/networks/analysis/all_stages_subnetwork_edges.csv"
GENES_CSV = "results/networks/analysis/all_stages_network_genes.csv"
NETWORK_CSV_DIR = "results/networks/har_csn_atlas/csv"
OUTPUT_DIR = "results/networks/analysis"
EDGE_COLS = ["Region", "CellType", "TF", "Target", "Stage", "Weight"]


def _load_06():
    """Load 06_all_stages_tf_target_pairs for load_all_edges."""
    _dir = os.path.dirname(os.path.abspath(__file__))
    _path = os.path.join(_dir, "06_all_stages_tf_target_pairs.py")
    if not os.path.isfile(_path):
        raise FileNotFoundError(f"06_all_stages_tf_target_pairs.py not found: {_path}")
    spec = importlib.util.spec_from_file_location("tf_pairs", _path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_subnetwork_from_tfs(tfs_df: pd.DataFrame, csv_dir: str) -> pd.DataFrame:
    """
    根据 all_stages_TFs 的 (Region, CellType, TF) 从 har_csn_atlas/csv 提取子网络边。
    返回 DataFrame，列: Region, CellType, TF, Target, Stage, Weight。
    """
    rc_tfs: dict[tuple[str, str], set[str]] = defaultdict(set)
    for _, r in tfs_df.iterrows():
        k = (str(r["Region"]).strip(), str(r["CellType"]).strip())
        rc_tfs[k].add(str(r["TF"]).strip())

    mod = _load_06()
    all_edges = mod.load_all_edges(csv_dir)
    rows = []
    for e in all_edges:
        k = (e["Region"], e["CellType"])
        if k not in rc_tfs or e["TF"] not in rc_tfs[k]:
            continue
        rows.append({c: e[c] for c in EDGE_COLS})

    return pd.DataFrame(rows)


def extract_genes_from_edges(
    edges_df: pd.DataFrame,
    tfs_df: pd.DataFrame,
    dimension_cols: tuple[str, ...] = ("Region", "CellType"),
    include_stage: bool = True,
    include_n_regulators: bool = True,
    include_n_stages: bool = True,
) -> pd.DataFrame:
    """
    Build a gene table from edges and all-stages TFs.

    - dimension_cols: required dimension columns for grouping (default Region, CellType).
    - include_stage: if True, add Stage to output and aggregate per (Region, CellType, Stage).
    - TFs: from tfs_df; when include_stage, expand stage_list to one row per (Region, CellType, Stage, gene=TF), is_TF=TRUE.
    - Targets: from edges (Region, CellType, Stage, Target); if gene is also in TF set for that (R,C) and Stage, is_TF=TRUE (same as 03).

    Output columns: gene, is_TF, Region, CellType [, Stage] [, n_regulators] [, n_stages].
    """
    needed = ("Region", "CellType", "TF", "Target", "Stage")
    for c in needed:
        if c not in edges_df.columns:
            raise ValueError(f"edges_df missing column: {c}")
    for c in ("Region", "CellType", "TF", "stage_list"):
        if c not in tfs_df.columns:
            raise ValueError(f"tfs_df missing column: {c}")

    # TF set per (Region, CellType, Stage): TFs that have this Stage in stage_list
    tf_set: dict[tuple, set[str]] = defaultdict(set)
    tf_n_stages: dict[tuple, int] = {}
    for _, r in tfs_df.iterrows():
        region, cell, tf = r["Region"], r["CellType"], r["TF"]
        stages = [s.strip() for s in str(r.get("stage_list", "") or "").split(";") if s.strip()]
        n = int(r.get("n_stages") or len(stages) or 0)
        for st in stages:
            k = (region, cell, st)
            tf_set[k].add(tf)
            tf_n_stages[(region, cell, tf)] = n

    # n_regulators: distinct TF count per (Region, CellType, Stage, Target)
    n_reg: dict[tuple, set[str]] = defaultdict(set)
    for _, r in edges_df.iterrows():
        k = (r["Region"], r["CellType"], r["Stage"], r["Target"])
        n_reg[k].add(r["TF"])

    rows: list[dict] = []
    seen: set[tuple] = set()

    # 1) TF rows
    for _, r in tfs_df.iterrows():
        region, cell, tf = r["Region"], r["CellType"], r["TF"]
        stages = [s.strip() for s in str(r.get("stage_list", "") or "").split(";") if s.strip()]
        n_st = int(r.get("n_stages") or len(stages) or 0)
        for st in (stages if include_stage else [None]):
            key = (region, cell, st, tf) if include_stage else (region, cell, tf)
            if key in seen:
                continue
            seen.add(key)
            row = {
                "gene": tf,
                "is_TF": "TRUE",
                "Region": region,
                "CellType": cell,
            }
            if include_stage and st is not None:
                row["Stage"] = st
            if include_n_stages:
                row["n_stages"] = n_st
            if include_n_regulators:
                row["n_regulators"] = ""
            rows.append(row)

    # 2) Target rows (skip if (R,C[,Stage],gene) already as TF)
    for _, r in edges_df.iterrows():
        region, cell, st, tg = r["Region"], r["CellType"], r["Stage"], r["Target"]
        key = (region, cell, st, tg) if include_stage else (region, cell, tg)
        if key in seen:
            continue
        is_tf = tg in tf_set.get((region, cell, st), set())
        seen.add(key)
        row = {
            "gene": tg,
            "is_TF": "TRUE" if is_tf else "FALSE",
            "Region": region,
            "CellType": cell,
        }
        if include_stage:
            row["Stage"] = st
        if include_n_regulators:
            row["n_regulators"] = len(n_reg.get((region, cell, st, tg), set()))
        if include_n_stages:
            row["n_stages"] = ""
        rows.append(row)

    out = pd.DataFrame(rows)
    # column order
    base = ["gene", "is_TF", "Region", "CellType"]
    if include_stage:
        base.append("Stage")
    if include_n_regulators:
        base.append("n_regulators")
    if include_n_stages:
        base.append("n_stages")
    out = out[[c for c in base if c in out.columns]]
    by = ["Region", "CellType"] + (["Stage"] if include_stage else []) + ["gene", "is_TF"]
    asc = [True] * (2 + (1 if include_stage else 0)) + [True, False]
    return out.sort_values(by=by, ascending=asc).reset_index(drop=True)


def run(
    tfs_csv: str = TFS_CSV,
    out_subnetwork: str = SUBNETWORK_CSV,
    out_genes: str = GENES_CSV,
    csv_dir: str = NETWORK_CSV_DIR,
    out_dir: str = OUTPUT_DIR,
    dimension_cols: tuple[str, ...] = ("Region", "CellType"),
    include_stage: bool = True,
    regions: Optional[List[str]] = None,
    celltypes: Optional[List[str]] = None,
) -> None:
    os.makedirs(out_dir, exist_ok=True)
    tfs_df = pd.read_csv(tfs_csv, dtype=str)
    
    # Filter by Region and CellType if specified
    if regions is not None:
        tfs_df = tfs_df[tfs_df["Region"].isin(regions)]
        log_message(f"Filtered by regions: {regions}", message_type="info")
    if celltypes is not None:
        tfs_df = tfs_df[tfs_df["CellType"].isin(celltypes)]
        log_message(f"Filtered by celltypes: {celltypes}", message_type="info")
    
    log_message(f"Loaded all_stages TFs: {len(tfs_df):,} rows from {tfs_csv}", message_type="info")

    edges_df = load_subnetwork_from_tfs(tfs_df, csv_dir)
    if edges_df.empty:
        log_message("No subnetwork edges; cannot build gene table.", message_type="warning")
        return

    log_message(f"Extracted subnetwork: {len(edges_df):,} edges", message_type="info")

    sub_path = os.path.join(out_dir, os.path.basename(out_subnetwork)) if not os.path.isabs(out_subnetwork) else out_subnetwork
    edges_df.to_csv(sub_path, index=False)
    log_message(f"Saved subnetwork to {sub_path}", message_type="success")

    genes_df = extract_genes_from_edges(
        edges_df,
        tfs_df,
        dimension_cols=dimension_cols,
        include_stage=include_stage,
        include_n_regulators=True,
        include_n_stages=True,
    )
    out_path = os.path.join(out_dir, os.path.basename(out_genes)) if not os.path.isabs(out_genes) else out_genes
    genes_df.to_csv(out_path, index=False)
    log_message(f"Saved {len(genes_df):,} gene rows to {out_path}", message_type="success")


def main():
    pairs = [
        ("Prefrontal cortex", "Excitatory neurons"),
        ("Prefrontal cortex", "Astrocytes")
    ]
    
    regions = sorted(list(set([pair[0] for pair in pairs])))
    celltypes = sorted(list(set([pair[1] for pair in pairs])))
    
    log_message(f"Processing {len(pairs)} (Region, CellType) pairs", message_type="info")
    for region, celltype in pairs:
        log_message(f"  - {region}, {celltype}", message_type="info")

    tfs_csv = TFS_CSV
    out_subnetwork = SUBNETWORK_CSV
    out_genes = GENES_CSV
    csv_dir = NETWORK_CSV_DIR
    out_dir = OUTPUT_DIR
    include_stage = True

    run(
        tfs_csv=tfs_csv,
        out_subnetwork=out_subnetwork,
        out_genes=out_genes,
        csv_dir=csv_dir,
        out_dir=out_dir,
        include_stage=include_stage,
        regions=regions,
        celltypes=celltypes,
    )


if __name__ == "__main__":
    main()
