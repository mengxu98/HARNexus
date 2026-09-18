#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Find TFs, target genes, or TF-gene pairs that appear in **every** stage of each (Region, CellType).

"All stages" = stages for which a network CSV exists: {Region}_{Stage}_{CellType}.csv.

Modes (--mode):
- TFs:      TFs with >=1 edge (as regulator) in every stage. Use --with-targets to get their
            target edges by stage (for target gene dynamics).
- genes:    Target genes with >=1 edge (as target) in every stage. Use --with-regulators to get
            their regulator edges by stage (for upstream TF dynamics).
- TF-gene:  (TF, Target) pairs with an edge in every stage (same as original behavior).

Optional --region and --celltype filter to a specific Region and/or CellType; omit to search all.

Inputs: results/networks/har_csn_atlas/csv/*.csv (regulator, target, weight).
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message

# Public API for programmatic use
__all__ = [
    "run",
    "find_all_stages_tfs",
    "find_all_stages_genes",
    "find_all_stages_tf_gene_pairs",
    "get_rc_stages_from_files",
    "load_all_edges",
    "MODE_TFS",
    "MODE_GENES",
    "MODE_TF_GENE",
]

MODE_TFS = "TFs"
MODE_GENES = "genes"
MODE_TF_GENE = "TF-gene"
MODES = (MODE_TFS, MODE_GENES, MODE_TF_GENE)
NETWORK_CSV_DIR = "results/networks/har_csn_atlas/csv"
OUTPUT_DIR = "results/networks/analysis"
REQUIRED_COLS = ["regulator", "target", "weight"]


def get_components(filename: str) -> dict:
    """Parse brain_region, stage, cell_type from CSV filename."""
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


def get_rc_stages_from_files(csv_dir: str) -> dict[tuple[str, str], list[str]]:
    """
    For each (Region, CellType), return the sorted list of stages for which a CSV exists.
    Only (Region, CellType) with >=2 stages are returned.
    """
    files = glob.glob(os.path.join(csv_dir, "*.csv"))
    rc_stages: dict[tuple[str, str], set[str]] = defaultdict(set)
    for f in files:
        comp = get_components(f)
        key = (comp["brain_region"], comp["cell_type"])
        rc_stages[key].add(comp["stage"])
    return {k: sorted(v) for k, v in rc_stages.items() if len(v) >= 2}


def _filter_rc(
    rc_multi: dict[tuple[str, str], list[str]],
    region: str | list[str] | None,
    celltype: str | list[str] | None,
) -> dict[tuple[str, str], list[str]]:
    """Subset rc_multi by optional region and/or celltype. None means no filter.
    region/celltype can be a str, or list[str] to match any of."""
    if region is None and celltype is None:
        return rc_multi

    def _match(val: str | list[str] | None, key: str) -> bool:
        if val is None:
            return True
        if isinstance(val, list):
            return key in val
        return key == val

    return {
        (r, c): stages
        for (r, c), stages in rc_multi.items()
        if _match(region, r) and _match(celltype, c)
    }


def _to_float(x) -> float | None:
    if x is None or (isinstance(x, str) and x.strip() == ""):
        return None
    try:
        v = float(x)
        return v
    except (ValueError, TypeError):
        return None


def load_all_edges(csv_dir: str) -> list[dict]:
    """Load all network CSVs; return list of dicts (Region, Stage, CellType, TF, Target, Weight)."""
    files = glob.glob(os.path.join(csv_dir, "*.csv"))
    rows: list[dict] = []

    for f in files:
        comp = get_components(f)
        region = comp["brain_region"]
        stage = comp["stage"]
        cell = comp["cell_type"]

        try:
            with open(f, newline="", encoding="utf-8") as fp:
                r = csv.DictReader(fp)
                if r.fieldnames is None or not all(
                    c in r.fieldnames for c in REQUIRED_COLS
                ):
                    log_message(
                        f"Skip {os.path.basename(f)}: missing {REQUIRED_COLS}",
                        message_type="warning",
                    )
                    continue
                for rec in r:
                    w = _to_float(rec.get("weight"))
                    if w is None or w == 0.0:
                        continue
                    rows.append(
                        {
                            "Region": region,
                            "Stage": stage,
                            "CellType": cell,
                            "TF": (rec.get("regulator") or "").strip(),
                            "Target": (rec.get("target") or "").strip(),
                            "Weight": w,
                        }
                    )
        except Exception as e:
            log_message(f"Skip {os.path.basename(f)}: {e}", message_type="warning")
            continue

    return rows


def find_all_stages_tf_gene_pairs(
    edges: list[dict],
    rc_multi: dict[tuple[str, str], list[str]],
    region: str | list[str] | None = None,
    celltype: str | list[str] | None = None,
) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    """
    Keep (TF, Target) pairs that have an edge in **every** stage. Optional region/celltype filter.

    Returns: (pairs, summary, edges_out, target_genes)
    """
    rc = _filter_rc(rc_multi, region, celltype)
    pair_stages: dict[tuple[str, str], dict[tuple[str, str], set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    pair_weights: dict[
        tuple[str, str], dict[tuple[str, str], list[tuple[str, float]]]
    ] = defaultdict(lambda: defaultdict(list))

    for r in edges:
        key = (r["Region"], r["CellType"])
        if key not in rc:
            continue
        pt = (r["TF"], r["Target"])
        if not r["TF"] or not r["Target"]:
            continue
        pair_stages[key][pt].add(r["Stage"])
        pair_weights[key][pt].append((r["Stage"], r["Weight"]))

    pairs_rows: list[dict] = []
    edges_out_rows: list[dict] = []

    for (region, cell), stages_all in sorted(rc.items()):
        st_set = set(stages_all)
        for (tf, target), st_present in pair_stages[(region, cell)].items():
            if st_present != st_set:
                continue
            wlist = pair_weights[(region, cell)][(tf, target)]
            abs_weights = [abs(w) for _, w in wlist]
            sum_abs = sum(abs_weights)
            mean_abs = sum_abs / len(abs_weights) if abs_weights else 0.0

            pairs_rows.append(
                {
                    "Region": region,
                    "CellType": cell,
                    "TF": tf,
                    "Target": target,
                    "n_stages": len(stages_all),
                    "stage_list": ";".join(stages_all),
                    "sum_abs_weight": round(sum_abs, 6),
                    "mean_abs_weight": round(mean_abs, 6),
                }
            )
            for st, w in wlist:
                edges_out_rows.append(
                    {
                        "Region": region,
                        "CellType": cell,
                        "Stage": st,
                        "TF": tf,
                        "Target": target,
                        "Weight": w,
                    }
                )

    # Summary: only (Region, CellType) that have >=1 all-stages pair; skip n_pairs==0
    summary_rows: list[dict] = []
    for (region, cell), stages_all in sorted(rc.items()):
        sub = [p for p in pairs_rows if p["Region"] == region and p["CellType"] == cell]
        n_pairs = len(sub)
        if n_pairs == 0:
            continue
        n_tfs = len({p["TF"] for p in sub})
        n_tg = len({p["Target"] for p in sub})
        summary_rows.append(
            {
                "Region": region,
                "CellType": cell,
                "n_stages": len(stages_all),
                "stage_list": ";".join(stages_all),
                "n_pairs": n_pairs,
                "n_unique_TFs": n_tfs,
                "n_unique_targets": n_tg,
            }
        )

    # Target genes: unique Target per (Region, CellType) with n_regulators
    tg_map: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    for p in pairs_rows:
        k = (p["Region"], p["CellType"], p["Target"])
        tg_map[k].append(p["TF"])
    target_genes_rows: list[dict] = []
    for (region, cell, target), tfs in sorted(tg_map.items()):
        first = next(
            p
            for p in pairs_rows
            if p["Region"] == region and p["CellType"] == cell and p["Target"] == target
        )
        target_genes_rows.append(
            {
                "Region": region,
                "CellType": cell,
                "Target": target,
                "n_stages": first["n_stages"],
                "stage_list": first["stage_list"],
                "n_regulators": len(set(tfs)),
            }
        )

    return pairs_rows, summary_rows, edges_out_rows, target_genes_rows


def find_all_stages_tfs(
    edges: list[dict],
    rc_multi: dict[tuple[str, str], list[str]],
    region: str | list[str] | None = None,
    celltype: str | list[str] | None = None,
) -> tuple[list[dict], list[dict]]:
    """
    TFs that have >=1 edge (as regulator) in every stage. Optional region/celltype filter.

    Returns:
      - tfs_rows: [{Region, CellType, TF, n_stages, stage_list}]
      - edges_for_tfs: all edges where TF is in the all-stages TF set (Region, CellType, TF, Target, Stage, Weight)
        for studying target gene dynamics across stages.
    """
    rc = _filter_rc(rc_multi, region, celltype)
    tf_stages: dict[tuple[str, str], dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )

    for r in edges:
        key = (r["Region"], r["CellType"])
        if key not in rc or not r["TF"]:
            continue
        tf_stages[key][r["TF"]].add(r["Stage"])

    tfs_rows: list[dict] = []
    all_stages_tf_set: dict[tuple[str, str], set[str]] = {}
    for (r, c), stages_all in sorted(rc.items()):
        st_set = set(stages_all)
        for tf, st_present in tf_stages.get((r, c), {}).items():
            if st_present != st_set:
                continue
            tfs_rows.append(
                {
                    "Region": r,
                    "CellType": c,
                    "TF": tf,
                    "n_stages": len(stages_all),
                    "stage_list": ";".join(stages_all),
                }
            )
            all_stages_tf_set.setdefault((r, c), set()).add(tf)

    edges_for_tfs: list[dict] = []
    for r in edges:
        key = (r["Region"], r["CellType"])
        if key not in all_stages_tf_set or r["TF"] not in all_stages_tf_set[key]:
            continue
        edges_for_tfs.append(
            {
                "Region": r["Region"],
                "CellType": r["CellType"],
                "TF": r["TF"],
                "Target": r["Target"],
                "Stage": r["Stage"],
                "Weight": r["Weight"],
            }
        )
    return tfs_rows, edges_for_tfs


def find_all_stages_genes(
    edges: list[dict],
    rc_multi: dict[tuple[str, str], list[str]],
    region: str | list[str] | None = None,
    celltype: str | list[str] | None = None,
) -> tuple[list[dict], list[dict]]:
    """
    Target genes that have >=1 edge (as target) in every stage. Optional region/celltype filter.

    Returns:
      - genes_rows: [{Region, CellType, Target, n_stages, stage_list, n_regulators}]
      - edges_for_genes: all edges where Target is in the all-stages gene set (Region, CellType, TF, Target, Stage, Weight)
        for studying how upstream TFs change across stages.
    """
    rc = _filter_rc(rc_multi, region, celltype)
    tg_stages: dict[tuple[str, str], dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    tg_tfs: dict[tuple[str, str, str], set[str]] = defaultdict(set)

    for r in edges:
        key = (r["Region"], r["CellType"])
        if key not in rc or not r["Target"]:
            continue
        tg_stages[key][r["Target"]].add(r["Stage"])
        tg_tfs[(r["Region"], r["CellType"], r["Target"])].add(r["TF"])

    genes_rows: list[dict] = []
    all_stages_gene_set: dict[tuple[str, str], set[str]] = {}
    for (r, c), stages_all in sorted(rc.items()):
        st_set = set(stages_all)
        for tg, st_present in tg_stages.get((r, c), {}).items():
            if st_present != st_set:
                continue
            genes_rows.append(
                {
                    "Region": r,
                    "CellType": c,
                    "Target": tg,
                    "n_stages": len(stages_all),
                    "stage_list": ";".join(stages_all),
                    "n_regulators": len(tg_tfs.get((r, c, tg), set())),
                }
            )
            all_stages_gene_set.setdefault((r, c), set()).add(tg)

    edges_for_genes: list[dict] = []
    for r in edges:
        key = (r["Region"], r["CellType"])
        if (
            key not in all_stages_gene_set
            or r["Target"] not in all_stages_gene_set[key]
        ):
            continue
        edges_for_genes.append(
            {
                "Region": r["Region"],
                "CellType": r["CellType"],
                "TF": r["TF"],
                "Target": r["Target"],
                "Stage": r["Stage"],
                "Weight": r["Weight"],
            }
        )
    return genes_rows, edges_for_genes


def _write_csv(path: str, rows: list[dict], fieldnames: list[str]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def parse_args():
    p = argparse.ArgumentParser(
        description="Find TFs, genes, or TF-gene pairs present in all stages (HAR-CSN atlas)."
    )
    p.add_argument(
        "--mode",
        choices=MODES,
        default=MODE_TFS,
        help="TFs=all-stages TFs; genes=all-stages target genes; TF-gene=all-stages pairs (default).",
    )
    p.add_argument(
        "--region", default=None, help="Filter to this Region (omit for all)."
    )
    p.add_argument(
        "--celltype", default=None, help="Filter to this CellType (omit for all)."
    )
    p.add_argument(
        "--with-targets",
        action="store_true",
        help="[TFs mode] Also output TF->Target edges by stage for target dynamics.",
    )
    p.add_argument(
        "--with-regulators",
        action="store_true",
        help="[genes mode] Also output TF->Target edges by stage for upstream TF dynamics.",
    )
    p.add_argument(
        "--csv_dir", default=NETWORK_CSV_DIR, help="har_csn_atlas csv directory"
    )
    p.add_argument("--out_dir", default=OUTPUT_DIR, help="Output directory")
    return p.parse_args()


def run(
    csv_dir: str,
    out_dir: str,
    mode: str = MODE_TFS,
    region: str | list[str] = None,
    celltype: str | list[str] = None,
    with_targets: bool = True,
    with_regulators: bool = True,
    *,
    _edges: list[dict] | None = None,
    _rc_multi: dict | None = None,
) -> None:
    """
    Main logic: load data (or use provided _edges, _rc_multi), run the chosen mode, write CSVs.
    Can be called programmatically or from main().
    """
    os.makedirs(out_dir, exist_ok=True)
    if _edges is not None and _rc_multi is not None:
        edges, rc_multi = _edges, _rc_multi
    else:
        edges = load_all_edges(csv_dir)
        rc_multi = get_rc_stages_from_files(csv_dir)

    if mode == MODE_TFS:
        tfs_rows, edges_tfs = find_all_stages_tfs(edges, rc_multi, region, celltype)
        out_tfs = os.path.join(out_dir, "all_stages_TFs.csv")
        _write_csv(
            out_tfs, tfs_rows, ["Region", "CellType", "TF", "n_stages", "stage_list"]
        )
        log_message(f"All-stages TFs: {len(tfs_rows):,}", message_type="success")
        log_message(f"Saved: {out_tfs}", message_type="success")
        if with_targets:
            out_te = os.path.join(out_dir, "all_stages_TF_target_edges.csv")
            _write_csv(
                out_te,
                edges_tfs,
                ["Region", "CellType", "TF", "Target", "Stage", "Weight"],
            )
            log_message(
                f"TF target edges (for dynamics): {len(edges_tfs):,} -> {out_te}",
                message_type="success",
            )

    elif mode == MODE_GENES:
        genes_rows, edges_genes = find_all_stages_genes(
            edges, rc_multi, region, celltype
        )
        out_genes = os.path.join(out_dir, "all_stages_genes.csv")
        _write_csv(
            out_genes,
            genes_rows,
            ["Region", "CellType", "Target", "n_stages", "stage_list", "n_regulators"],
        )
        log_message(f"All-stages genes: {len(genes_rows):,}", message_type="success")
        log_message(f"Saved: {out_genes}", message_type="success")
        if with_regulators:
            out_re = os.path.join(out_dir, "all_stages_gene_regulator_edges.csv")
            _write_csv(
                out_re,
                edges_genes,
                ["Region", "CellType", "TF", "Target", "Stage", "Weight"],
            )
            log_message(
                f"Gene regulator edges (for dynamics): {len(edges_genes):,} -> {out_re}",
                message_type="success",
            )

    else:  # MODE_TF_GENE
        pairs_rows, summary_rows, edges_out_rows, target_genes_rows = (
            find_all_stages_tf_gene_pairs(edges, rc_multi, region, celltype)
        )
        _write_csv(
            os.path.join(out_dir, "all_stages_tf_target_pairs.csv"),
            pairs_rows,
            [
                "Region",
                "CellType",
                "TF",
                "Target",
                "n_stages",
                "stage_list",
                "sum_abs_weight",
                "mean_abs_weight",
            ],
        )
        _write_csv(
            os.path.join(out_dir, "all_stages_tf_target_summary.csv"),
            summary_rows,
            [
                "Region",
                "CellType",
                "n_stages",
                "stage_list",
                "n_pairs",
                "n_unique_TFs",
                "n_unique_targets",
            ],
        )
        _write_csv(
            os.path.join(out_dir, "all_stages_tf_target_edges.csv"),
            edges_out_rows,
            ["Region", "CellType", "Stage", "TF", "Target", "Weight"],
        )
        _write_csv(
            os.path.join(out_dir, "all_stages_target_genes.csv"),
            target_genes_rows,
            ["Region", "CellType", "Target", "n_stages", "stage_list", "n_regulators"],
        )
        log_message(
            f"All-stages TF-target pairs: {len(pairs_rows):,}", message_type="success"
        )
        log_message(
            f"All-stages target genes: {len(target_genes_rows):,}", message_type="info"
        )
        log_message(
            f"(Region,CellType) groups: {len(summary_rows):,}", message_type="info"
        )
        log_message(
            "Saved: all_stages_tf_target_*.csv, all_stages_target_genes.csv",
            message_type="success",
        )


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

    csv_dir = "results/networks/har_csn_atlas/csv"
    out_dir = "results/networks/analysis"
    mode = MODE_TF_GENE
    with_targets = True
    with_regulators = True

    log_message("Loading HAR-CSN network CSVs...", message_type="running")
    edges = load_all_edges(csv_dir)
    log_message(f"Total edges (Weight!=0): {len(edges):,}", message_type="info")
    if not edges:
        log_message("No edge data; exiting.", message_type="warning")
        return

    rc_multi = get_rc_stages_from_files(csv_dir)
    rc = _filter_rc(rc_multi, regions, celltypes)
    log_message(
        f"(Region,CellType) with >=2 stage networks: {len(rc_multi):,} (after filter: {len(rc):,})",
        message_type="info",
    )
    log_message(f"Mode: {mode}", message_type="info")

    run(
        csv_dir,
        out_dir,
        mode=mode,
        region=regions,
        celltype=celltypes,
        with_targets=with_targets,
        with_regulators=with_regulators,
        _edges=edges,
        _rc_multi=rc_multi,
    )


if __name__ == "__main__":
    main()
