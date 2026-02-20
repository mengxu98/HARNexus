import numpy as np
import pandas as pd
import plotly.graph_objects as go
import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import sys
import os
import io
import base64
from datetime import datetime
from typing import List, Optional, Dict, Any
import networkx as nx

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from functions.utils import (
    COLOR_CELLTYPES,
    COLOR_STAGES,
    DEFAULT_NODE_COLOR,
    log_message,
)
from functions.utils_network import (
    load_network_data,
    load_har_tf_data,
    generate_colors,
)

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))


def _resolve_path(path: Optional[str]) -> Optional[str]:
    """Resolve path relative to project root if it is not absolute."""
    if path is None or not path:
        return path
    if os.path.isabs(path):
        return path
    return os.path.join(PROJECT_ROOT, path)


class InteractiveNetworkVisualizer:
    def __init__(
        self,
        data_path: Optional[str] = None,
        har_tf_path: Optional[str] = None,
        network_dir: str = "results/networks/analysis",
        regions: Optional[List[str]] = None,
        stages: Optional[List[str]] = None,
        celltypes: Optional[List[str]] = None,
        tfs: Optional[List[str]] = None,
        targets: Optional[List[str]] = None,
        hars: Optional[List[str]] = None,
    ):
        """
        Initialize network visualizer

        Parameters:
        -----------
        data_path : str, optional
            Directly specify data file path (if provided, will use this path)
        har_tf_path : str, optional
            HAR-TF mapping file path
        network_dir : str
            Data directory base path (default: results/networks/analysis)
        regions : List[str], optional
            List of specified brain regions
        stages : List[str], optional
            List of specified developmental stages
        celltypes : List[str], optional
            List of specified cell types
        tfs : List[str], optional
            List of specified TFs (for filtering)
        targets : List[str], optional
            List of specified target genes (for filtering)
        hars : List[str], optional
            List of specified HARs (for filtering)
        """
        self.data_path = _resolve_path(data_path)
        self.har_tf_path = _resolve_path(
            har_tf_path
            if har_tf_path is not None
            else "results/har_tf/human/har_tf_pairs_scores.csv"
        )
        self.network_dir = _resolve_path(network_dir)
        self.filter_regions = regions
        self.filter_stages = stages
        self.filter_celltypes = celltypes
        self.filter_tfs = tfs
        self.filter_targets = targets
        self.filter_hars = hars
        self.network_data = None
        self.har_tf_data = None
        self.subnetwork_data = None
        self._cache = {}
        self.load_data()

    def load_data(self):
        try:
            if self.har_tf_path:
                self.har_tf_data = load_har_tf_data(self.har_tf_path)
            else:
                log_message("No HAR-TF data provided.", message_type="info")
                self.har_tf_data = None

            if self.data_path:
                log_message(
                    f"Loading network data from {self.data_path}...",
                    message_type="info",
                )
                self.network_data = pd.read_csv(self.data_path)
                if "Weight" not in self.network_data.columns:
                    if "weight" in self.network_data.columns:
                        self.network_data = self.network_data.rename(
                            columns={"weight": "Weight"}
                        )
                    else:
                        raise ValueError("Weight column not found in data")
                if "TF" not in self.network_data.columns:
                    if "regulator" in self.network_data.columns:
                        self.network_data = self.network_data.rename(
                            columns={"regulator": "TF"}
                        )
                    else:
                        raise ValueError("TF column not found in data")
                if "Target" not in self.network_data.columns:
                    if "target" in self.network_data.columns:
                        self.network_data = self.network_data.rename(
                            columns={"target": "Target"}
                        )
                    else:
                        raise ValueError("Target column not found in data")
                if "Stage" not in self.network_data.columns:
                    self.network_data["Stage"] = "Unknown_Stage"
                if "Region" not in self.network_data.columns:
                    self.network_data["Region"] = "Unknown_Region"
                if "CellType" not in self.network_data.columns:
                    self.network_data["CellType"] = "Unknown_CellType"
            else:
                self.network_data = load_network_data(
                    network_dir=self.network_dir,
                    regions=self.filter_regions,
                    stages=self.filter_stages,
                    celltypes=self.filter_celltypes,
                    tfs=self.filter_tfs,
                    targets=self.filter_targets,
                    hars=self.filter_hars,
                    har_tf_data=self.har_tf_data,
                )

            log_message("Cleaning and optimizing data...", message_type="info")

            self.network_data["Weight"] = self.network_data["Weight"].fillna(0)
            self.network_data["TF"] = self.network_data["TF"].fillna("Unknown_TF")
            self.network_data["Target"] = self.network_data["Target"].fillna(
                "Unknown_Target"
            )
            self.network_data["Stage"] = self.network_data["Stage"].fillna(
                "Unknown_Stage"
            )
            self.network_data["Region"] = self.network_data["Region"].fillna(
                "Unknown_Region"
            )
            self.network_data["CellType"] = self.network_data["CellType"].fillna(
                "Unknown_CellType"
            )
            self.network_data = self.network_data.dropna(subset=["TF", "Target"])

            self.network_data = self.network_data[
                (self.network_data["Weight"] > 0)
                & (self.network_data["TF"] != "Unknown_TF")
                & (self.network_data["Target"] != "Unknown_Target")
            ].copy()

            self.network_data["Weight"] = self.network_data["Weight"].astype("float32")

            log_message(
                f"Network data loaded successfully. Shape: {self.network_data.shape}",
                message_type="success",
            )

        except Exception as e:
            log_message(f"Error loading data: {str(e)}", message_type="error")
            sys.exit(1)

    def get_available_options(self):
        return {
            "tfs": sorted(self.network_data["TF"].unique().tolist()),
            "targets": sorted(self.network_data["Target"].unique().tolist()),
            "stages": sorted(self.network_data["Stage"].unique().tolist()),
            "regions": sorted(self.network_data["Region"].unique().tolist()),
            "celltypes": sorted(self.network_data["CellType"].unique().tolist()),
            "har_tfs": sorted(self.har_tf_data["TF"].unique().tolist())
            if self.har_tf_data is not None
            else [],
            "hars": sorted(self.har_tf_data["HAR"].unique().tolist())
            if self.har_tf_data is not None
            else [],
        }

    def get_top_nodes(self, node_type: str, top_n: int, filters: Dict[str, Any] = None):
        data = self.network_data.copy()

        if filters:
            for key, value in filters.items():
                if value and key in data.columns:
                    if isinstance(value, list):
                        data = data[data[key].isin(value)]
                    else:
                        data = data[data[key] == value]

        if node_type == "tf":
            top_nodes = (
                data.groupby("TF")["Weight"].sum().nlargest(top_n).index.tolist()
            )
        elif node_type == "target":
            top_nodes = (
                data.groupby("Target")["Weight"].sum().nlargest(top_n).index.tolist()
            )
        else:
            top_nodes = []

        return top_nodes

    def create_sankey_diagram(
        self,
        show_har: bool = True,
        top_n_tfs: int = 10,
        top_n_targets: int = 10,
        specific_tfs: List[str] = None,
        specific_targets: List[str] = None,
        specific_stages: List[str] = None,
        specific_regions: List[str] = None,
        specific_celltypes: List[str] = None,
        specific_hars: List[str] = None,
        only_show_levels: Optional[List[str]] = None,
        height: int = 700,
        width: int = 1200,
        theme: str = "light",
    ):
        flow_data = self.network_data.copy()

        if specific_hars and self.har_tf_data is not None:
            relevant_tfs = self.har_tf_data[
                self.har_tf_data["HAR"].isin(specific_hars)
            ]["TF"].unique()
            flow_data = flow_data[flow_data["TF"].isin(relevant_tfs)]

        if specific_celltypes:
            if isinstance(specific_celltypes, str):
                flow_data = flow_data[flow_data["CellType"] == specific_celltypes]
            else:
                flow_data = flow_data[flow_data["CellType"].isin(specific_celltypes)]
        if specific_regions:
            flow_data = flow_data[flow_data["Region"].isin(specific_regions)]
        if specific_stages:
            flow_data = flow_data[flow_data["Stage"].isin(specific_stages)]

        if len(flow_data) == 0:
            return None, "No data found for the specified filters", [], []

        def get_weighted_data(data, group_by_columns):
            result = data.groupby(group_by_columns)["Weight"].sum().reset_index()
            result["AbsWeight"] = result["Weight"].abs()
            return result

        priority_columns = []
        if specific_celltypes:
            priority_columns.append("CellType")
        if specific_stages:
            priority_columns.append("Stage")
        if specific_regions:
            priority_columns.append("Region")

        if specific_tfs:
            available_tfs = set(flow_data["TF"].unique())
            top_tfs = [tf for tf in specific_tfs if tf in available_tfs]
            if not top_tfs:
                return (
                    None,
                    "None of the specified TFs found in the data.",
                    [],
                    [],
                )
            if len(top_tfs) < len(specific_tfs):
                missing = set(specific_tfs) - set(top_tfs)
                log_message(
                    f"Warning: Some specified TFs not found: {missing}",
                    message_type="warning",
                )
        elif specific_targets:
            target_data = flow_data[flow_data["Target"].isin(specific_targets)]
            if len(target_data) == 0:
                return (
                    None,
                    "No specified target genes found in the data after filtering.",
                    [],
                    [],
                )
            tf_target_strength = (
                target_data.groupby(["TF", "Target"])["Weight"]
                .apply(lambda x: x.abs().sum())
                .reset_index(name="AbsWeight")
            )
            tf_strength = (
                tf_target_strength.groupby("TF")["AbsWeight"]
                .sum()
                .sort_values(ascending=False)
            )
            targets_in_data = set(tf_target_strength["Target"].unique())
            missing_targets = set(specific_targets) & targets_in_data
            cover_map = tf_target_strength.groupby("TF")["Target"].agg(set).to_dict()
            selected_tfs = []
            selected_set = set()
            effective_top_n_tfs = max(top_n_tfs, len(specific_targets))
            while missing_targets and len(selected_tfs) < effective_top_n_tfs:
                best_tf = None
                best_cover_n = 0
                best_strength = -1.0
                for tf, covered in cover_map.items():
                    if tf in selected_set:
                        continue
                    cover_n = len(covered & missing_targets)
                    if cover_n == 0:
                        continue
                    strength = float(tf_strength.get(tf, 0.0))
                    if (cover_n > best_cover_n) or (
                        cover_n == best_cover_n and strength > best_strength
                    ):
                        best_tf = tf
                        best_cover_n = cover_n
                        best_strength = strength
                if best_tf is None:
                    break
                selected_tfs.append(best_tf)
                selected_set.add(best_tf)
                missing_targets -= cover_map.get(best_tf, set())
            top_tfs = list(selected_tfs)
            if len(top_tfs) < top_n_tfs:
                for tf in tf_strength.index.tolist():
                    if tf not in selected_set:
                        top_tfs.append(tf)
                        selected_set.add(tf)
                        if len(top_tfs) >= top_n_tfs:
                            break
        else:
            if priority_columns:
                group_cols = priority_columns + ["TF"]
                weighted_data = get_weighted_data(flow_data, group_cols)
                top_tfs = (
                    weighted_data.groupby("TF")["AbsWeight"]
                    .sum()
                    .sort_values(ascending=False)
                    .head(top_n_tfs)
                    .index.tolist()
                )
            else:
                top_tfs = (
                    flow_data.groupby("TF")["Weight"]
                    .apply(lambda x: x.abs().sum())
                    .sort_values(ascending=False)
                    .head(top_n_tfs)
                    .index.tolist()
                )

        flow_data_filtered = flow_data[flow_data["TF"].isin(top_tfs)]
        if specific_targets:
            available_targets = set(flow_data_filtered["Target"].unique())
            top_targets = [t for t in specific_targets if t in available_targets]
            if not top_targets:
                return (
                    None,
                    "No specified target genes found for selected TFs.",
                    [],
                    [],
                )
        else:
            if priority_columns:
                group_cols = priority_columns + ["Target"]
                weighted_data = get_weighted_data(flow_data_filtered, group_cols)
                top_targets = (
                    weighted_data.groupby("Target")["AbsWeight"]
                    .sum()
                    .sort_values(ascending=False)
                    .head(top_n_targets)
                    .index.tolist()
                )
            else:
                top_targets = (
                    flow_data_filtered.groupby("Target")["Weight"]
                    .apply(lambda x: x.abs().sum())
                    .sort_values(ascending=False)
                    .head(top_n_targets)
                    .index.tolist()
                )

        flow_data = flow_data[
            flow_data["TF"].isin(top_tfs) & flow_data["Target"].isin(top_targets)
        ].copy()

        if len(flow_data) == 0:
            return None, "No connections found for selected nodes", [], []

        def _should_show_level(specific_vals):
            if specific_vals is None:
                return True
            if isinstance(specific_vals, str):
                return False
            return len(specific_vals) >= 2

        def _validate_and_filter_complete_paths(
            flow_data, top_tfs, top_targets, show_har, har_tf_data, levels_to_show
        ):
            filtered_data = flow_data.copy()
            valid_tfs = set(top_tfs)
            valid_targets = set(top_targets)

            if show_har and har_tf_data is not None:
                tfs_with_har = set(har_tf_data["TF"].unique())
                valid_tfs = valid_tfs & tfs_with_har

            tf_target_pairs = filtered_data.groupby("TF")["Target"].apply(set)
            valid_tfs = {
                tf
                for tf in valid_tfs
                if tf in tf_target_pairs
                and len(tf_target_pairs[tf] & valid_targets) > 0
            }

            target_tf_pairs = filtered_data.groupby("Target")["TF"].apply(set)
            valid_targets = {
                t
                for t in valid_targets
                if t in target_tf_pairs and len(target_tf_pairs[t] & valid_tfs) > 0
            }

            filtered_data = filtered_data[
                filtered_data["TF"].isin(valid_tfs)
                & filtered_data["Target"].isin(valid_targets)
            ].copy()

            for level in ["Stage", "Region", "CellType"]:
                if level in levels_to_show:
                    level_data = filtered_data.groupby([level, "TF", "Target"]).size()
                    valid_levels = set()
                    for (level_val, tf, target), _ in level_data.items():
                        if tf in valid_tfs and target in valid_targets:
                            valid_levels.add(level_val)
                    filtered_data = filtered_data[
                        filtered_data[level].isin(valid_levels)
                    ]

            final_valid_tfs = set()
            final_valid_targets = set()

            for tf in valid_tfs:
                tf_data = filtered_data[filtered_data["TF"] == tf]
                if len(tf_data) > 0:
                    tf_targets = set(tf_data["Target"].unique())
                    if len(tf_targets & valid_targets) > 0:
                        final_valid_tfs.add(tf)

            for target in valid_targets:
                target_data = filtered_data[filtered_data["Target"] == target]
                if len(target_data) > 0:
                    target_tfs = set(target_data["TF"].unique())
                    if len(target_tfs & final_valid_tfs) > 0:
                        final_valid_targets.add(target)

            filtered_data = filtered_data[
                filtered_data["TF"].isin(final_valid_tfs)
                & filtered_data["Target"].isin(final_valid_targets)
            ].copy()

            return filtered_data, list(final_valid_tfs), list(final_valid_targets)

        levels_to_show_preview = ["TF"]
        if show_har and self.har_tf_data is not None:
            levels_to_show_preview.insert(0, "HAR")
        if _should_show_level(specific_stages):
            levels_to_show_preview.append("Stage")
        elif not specific_stages:
            levels_to_show_preview.append("Stage")
        if _should_show_level(specific_regions):
            levels_to_show_preview.append("Region")
        elif not specific_regions:
            levels_to_show_preview.append("Region")
        if _should_show_level(specific_celltypes):
            levels_to_show_preview.append("CellType")
        elif not specific_celltypes:
            levels_to_show_preview.append("CellType")
        levels_to_show_preview.append("Target")

        flow_data, valid_tfs, valid_targets = _validate_and_filter_complete_paths(
            flow_data,
            top_tfs,
            top_targets,
            show_har,
            self.har_tf_data,
            levels_to_show_preview,
        )

        if len(flow_data) == 0:
            return (
                None,
                "No complete connections found after data integrity validation",
                [],
                [],
            )

        top_tfs = valid_tfs
        top_targets = valid_targets

        only_show_levels = only_show_levels or []
        stage_region_celltype = ["Stage", "Region", "CellType"]

        def _level_entries(level_key, specific_vals, df_col, only_show_alt=None):
            in_data = set(flow_data[df_col].unique())
            only_show = level_key in only_show_levels or (
                only_show_alt and only_show_alt in only_show_levels
            )
            if specific_vals and only_show:
                normalized = (
                    [specific_vals]
                    if isinstance(specific_vals, str)
                    else list(specific_vals)
                )
                return sorted(set(normalized) & in_data)
            return sorted(in_data)

        levels_to_show = ["TF"]
        level_nodes = {
            "HAR": [],
            "TF": top_tfs,
            "Stage": _level_entries("Stage", specific_stages, "Stage"),
            "Region": _level_entries("Region", specific_regions, "Region"),
            "CellType": _level_entries(
                "CellType",
                specific_celltypes,
                "CellType",
                "specific_celltypes",
            ),
            "Target": top_targets,
        }

        har_to_tf_map = {}

        if show_har and self.har_tf_data is not None:
            levels_to_show.insert(0, "HAR")
            har_tf_filtered = self.har_tf_data[self.har_tf_data["TF"].isin(top_tfs)]
            if specific_hars:
                har_tf_filtered = har_tf_filtered[
                    har_tf_filtered["HAR"].isin(specific_hars)
                ]

            har_tf_groups = (
                har_tf_filtered.groupby("TF")["HAR"]
                .agg(lambda x: sorted(x))
                .reset_index()
            )

            har_nodes = []
            for _, row in har_tf_groups.iterrows():
                tf = row["TF"]
                har_list = row["HAR"]
                first_har = har_list[0]
                har_count = len(har_list)
                har_display = f"HARs: {first_har}...({har_count})"
                har_nodes.append(har_display)
                har_to_tf_map[har_display] = tf

            level_nodes["HAR"] = har_nodes

        if not only_show_levels:
            if _should_show_level(specific_stages):
                levels_to_show.append("Stage")
            elif not specific_stages:
                levels_to_show.append("Stage")
            if _should_show_level(specific_regions):
                levels_to_show.append("Region")
            elif not specific_regions:
                levels_to_show.append("Region")
            if _should_show_level(specific_celltypes):
                levels_to_show.append("CellType")
            elif not specific_celltypes:
                levels_to_show.append("CellType")
        else:
            for L in stage_region_celltype:
                if L in only_show_levels or (
                    L == "CellType" and "specific_celltypes" in only_show_levels
                ):
                    levels_to_show.append(L)
        levels_to_show.append("Target")

        nodes = []
        node_labels = []
        node_to_idx = {}
        current_idx = 0
        italic_levels = {"TF", "Target"}

        for level in levels_to_show:
            level_nodes_list = level_nodes[level]
            for node in level_nodes_list if level_nodes_list else [None]:
                nodes.append(node)
                if level in italic_levels and node is not None:
                    node_labels.append(f"<i>{node}</i>")
                else:
                    node_labels.append(node)
                node_to_idx[node] = current_idx
                current_idx += 1

        final_width = width if width else 1000

        sources = []
        targets = []
        values = []
        link_colors = []
        all_connections = []

        for i in range(len(levels_to_show) - 1):
            current_level = levels_to_show[i]
            next_level = levels_to_show[i + 1]

            if current_level == "HAR" and next_level == "TF":
                for har_group in level_nodes["HAR"]:
                    tf = har_to_tf_map[har_group]
                    if har_group in node_to_idx and tf in node_to_idx:
                        sources.append(node_to_idx[har_group])
                        targets.append(node_to_idx[tf])
                        tf_weight = flow_data[flow_data["TF"] == tf]["Weight"].sum()
                        values.append(float(tf_weight))
                        all_connections.append(tf)
            else:
                grouped_data = (
                    flow_data.groupby([current_level, next_level])["Weight"]
                    .sum()
                    .reset_index()
                )

                for _, row in grouped_data.iterrows():
                    if (
                        row[current_level] in node_to_idx
                        and row[next_level] in node_to_idx
                    ):
                        sources.append(node_to_idx[row[current_level]])
                        targets.append(node_to_idx[row[next_level]])
                        values.append(float(row["Weight"]))

                        if current_level == "TF":
                            all_connections.append(row[current_level])
                        else:
                            mask = flow_data[next_level] == row[next_level]
                            if current_level != "HAR":
                                mask &= flow_data[current_level] == row[current_level]
                            path_tfs = flow_data[mask]["TF"].unique()
                            if len(path_tfs) > 0:
                                tf_weights = {
                                    tf: flow_data[(flow_data["TF"] == tf) & mask][
                                        "Weight"
                                    ].sum()
                                    for tf in path_tfs
                                }
                                dominant_tf = max(
                                    tf_weights.items(), key=lambda x: x[1]
                                )[0]
                                all_connections.append(dominant_tf)

        n_tfs = len(top_tfs)
        colors = self._generate_colors(n_tfs)
        tf_color_map = {tf: colors[i % len(colors)] for i, tf in enumerate(top_tfs)}

        link_colors = [tf_color_map[tf] for tf in all_connections]

        node_colors = []
        for level in levels_to_show:
            level_nodes_list = level_nodes[level]
            if level == "TF":
                for j, node in enumerate(
                    level_nodes_list if level_nodes_list else [None]
                ):
                    if node is not None and j < len(colors):
                        node_colors.append(colors[j % len(colors)])
                    else:
                        node_colors.append(DEFAULT_NODE_COLOR)
            elif level == "HAR":
                for har_group in level_nodes_list if level_nodes_list else [None]:
                    if har_group is not None:
                        tf = har_to_tf_map[har_group]
                        node_colors.append(tf_color_map[tf])
                    else:
                        node_colors.append(DEFAULT_NODE_COLOR)
            elif level == "Stage":
                for node in level_nodes_list if level_nodes_list else [None]:
                    if node is not None:
                        node_colors.append(COLOR_STAGES.get(node, DEFAULT_NODE_COLOR))
                    else:
                        node_colors.append(DEFAULT_NODE_COLOR)
            elif level == "CellType":
                for node in level_nodes_list if level_nodes_list else [None]:
                    if node is not None:
                        node_colors.append(
                            COLOR_CELLTYPES.get(node, DEFAULT_NODE_COLOR)
                        )
                    else:
                        node_colors.append(DEFAULT_NODE_COLOR)
            else:
                for node in level_nodes_list if level_nodes_list else [None]:
                    if node is not None:
                        mask = flow_data[level] == node
                        path_tfs = flow_data[mask]["TF"].unique()
                        if len(path_tfs) > 0:
                            tf_weights = {
                                tf: flow_data[(flow_data["TF"] == tf) & mask][
                                    "Weight"
                                ].sum()
                                for tf in path_tfs
                            }
                            dominant_tf = max(tf_weights.items(), key=lambda x: x[1])[0]
                            node_colors.append(tf_color_map[dominant_tf])
                        else:
                            node_colors.append(DEFAULT_NODE_COLOR)
                    else:
                        node_colors.append(DEFAULT_NODE_COLOR)

        SMALL_GAP = 0.055
        LARGE_GAP = 0.13
        small_gap_pairs = {("TF", "Stage"), ("Stage", "Region")}
        level_x = [0.02]
        for i in range(1, len(levels_to_show)):
            prev, nxt = levels_to_show[i - 1], levels_to_show[i]
            d = SMALL_GAP if (prev, nxt) in small_gap_pairs else LARGE_GAP
            level_x.append(level_x[-1] + d)
        total = level_x[-1] - level_x[0]
        if total > 0:
            level_x = [0.02 + (x - level_x[0]) / total * 0.96 for x in level_x]
        level_to_x = {levels_to_show[i]: level_x[i] for i in range(len(levels_to_show))}
        node_x = []
        node_y = []
        y_lo, y_hi = 0.03, 0.97
        for level in levels_to_show:
            lst = level_nodes[level] if level_nodes[level] else [None]
            n = len(lst)
            x = level_to_x[level]
            for j in range(n):
                node_x.append(x)
                node_y.append(y_lo + (j + 0.5) / max(n, 1) * (y_hi - y_lo))

        fig = go.Figure(
            data=[
                go.Sankey(
                    arrangement="fixed",
                    node=dict(
                        pad=15,
                        thickness=20,
                        line=dict(color="black", width=0.5),
                        label=node_labels,
                        color=node_colors,
                        x=node_x,
                        y=node_y,
                    ),
                    link=dict(
                        source=sources,
                        target=targets,
                        value=values,
                        color=link_colors,
                    ),
                )
            ]
        )

        title_parts = ["Network wiring"]
        if specific_celltypes:
            ct_str = (
                specific_celltypes
                if isinstance(specific_celltypes, str)
                else specific_celltypes[0]
            )
            title_parts.append(f"of {ct_str}")
        if specific_regions:
            title_parts.append(f"in {', '.join(specific_regions)}")
        if specific_stages:
            title_parts.append(f"at {', '.join(specific_stages)}")
        title_text = " ".join(title_parts)
        if specific_tfs:
            title_text += f"<br>Selected TFs: {', '.join(top_tfs)}"
        if specific_targets:
            title_text += f"<br>Selected Targets: {', '.join(top_targets)}"

        level_labels = {
            "HAR": "HARs",
            "CellType": "Cell types",
            "Stage": "Stages",
            "Region": "Regions",
            "TF": "TFs",
            "Target": "Target genes",
        }
        is_dark = theme == "dark" if theme else False
        template = "plotly_dark" if is_dark else "plotly_white"
        bg_color = "#121212" if is_dark else "#ffffff"
        text_color = "#ffffff" if is_dark else "#212529"
        col_annotations = [
            dict(
                x=level_to_x[level],
                y=0.03,
                xref="x",
                yref="paper",
                text=level_labels.get(level, level),
                showarrow=False,
                yanchor="top",
                xanchor="center",
                font=dict(size=10, color=text_color, family="Arial"),
            )
            for level in levels_to_show
        ]

        fig.update_layout(
            title=dict(
                text=title_text,
                x=0.1,
                y=0.9,
                xanchor="left",
                yanchor="top",
                font=dict(size=14, color=text_color, family="Arial"),
            ),
            font=dict(size=10, color=text_color, family="Arial"),
            height=height,
            width=final_width,
            template=template,
            paper_bgcolor=bg_color,
            plot_bgcolor=bg_color,
            margin=dict(b=52, l=80, r=80, t=80),
            autosize=False,
            annotations=col_annotations,
            xaxis=dict(visible=False, range=[0, 1]),
            yaxis=dict(visible=False),
        )

        return fig, "Success", top_tfs, top_targets

    def get_node_subnetwork(self, node_name: str, node_type: str, max_depth: int = 1):
        try:
            if node_type == "tf":
                subnetwork_data = self.network_data[
                    (self.network_data["TF"] == node_name)
                    & (self.network_data["Weight"] > 0)
                ].copy()
            elif node_type == "target":
                subnetwork_data = self.network_data[
                    (self.network_data["Target"] == node_name)
                    & (self.network_data["Weight"] > 0)
                ].copy()
            else:
                return None

            if len(subnetwork_data) == 0:
                return None

            if len(subnetwork_data) > 100:
                subnetwork_data = subnetwork_data.nlargest(100, "Weight")

            G = nx.DiGraph()

            edges_data = []
            for _, row in subnetwork_data.iterrows():
                edges_data.append(
                    (
                        row["TF"],
                        row["Target"],
                        {
                            "weight": row["Weight"],
                            "celltype": row["CellType"],
                            "stage": row["Stage"],
                            "region": row["Region"],
                        },
                    )
                )

            G.add_edges_from(edges_data)

            if node_name in G.nodes():
                try:
                    subgraph = nx.ego_graph(
                        G, node_name, radius=max_depth, undirected=False
                    )
                except Exception:
                    subgraph = G
            else:
                subgraph = G

            if len(subgraph.nodes()) > 50:
                node_weights = {}
                for node in subgraph.nodes():
                    total_weight = sum(
                        data.get("weight", 0)
                        for source, target, data in subgraph.edges(data=True)
                        if source == node or target == node
                    )
                    node_weights[node] = total_weight

                top_nodes = sorted(
                    node_weights.items(), key=lambda x: x[1], reverse=True
                )[:50]
                subgraph = subgraph.subgraph([node for node, _ in top_nodes])

            try:
                pos = nx.spring_layout(subgraph, k=2, iterations=30, seed=42)
            except Exception:
                pos = nx.random_layout(subgraph, seed=42)

            edge_x = []
            edge_y = []
            edge_info = []

            for edge in subgraph.edges():
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])

                edge_data = subgraph[edge[0]][edge[1]]
                edge_info.append(
                    {
                        "source": edge[0],
                        "target": edge[1],
                        "weight": edge_data.get("weight", 0),
                        "celltype": edge_data.get("celltype", "Unknown"),
                        "stage": edge_data.get("stage", "Unknown"),
                        "region": edge_data.get("region", "Unknown"),
                    }
                )

            node_x = []
            node_y = []
            node_text = []
            node_info = []

            for node in subgraph.nodes():
                x, y = pos[node]
                node_x.append(x)
                node_y.append(y)
                node_text.append(node)

                is_tf = any(edge[0] == node for edge in subgraph.edges())
                node_info.append(
                    {
                        "name": node,
                        "type": "TF" if is_tf else "Target",
                        "connections": subgraph.degree(node),
                    }
                )

            edge_trace = go.Scatter(
                x=edge_x,
                y=edge_y,
                line=dict(width=1.5, color="#888"),
                hoverinfo="none",
                mode="lines",
            )

            node_trace = go.Scatter(
                x=node_x,
                y=node_y,
                mode="markers+text",
                hoverinfo="text",
                text=node_text,
                textposition="middle center",
                marker=dict(
                    showscale=True,
                    colorscale="Viridis",
                    reversescale=False,
                    color=[],
                    size=15,
                    colorbar=dict(thickness=10, xanchor="left", titleside="right"),
                    line=dict(width=1, color="black"),
                ),
            )

            node_adjacencies = []
            node_text_info = []
            for node in subgraph.nodes():
                degree = subgraph.degree(node)
                node_adjacencies.append(degree)
                node_text_info.append(f"{node}<br>Edges: {degree}")

            node_trace.marker.color = node_adjacencies
            node_trace.text = node_text_info

            fig = go.Figure(
                data=[edge_trace, node_trace],
                layout=go.Layout(
                    title=f"Subnetwork for {node_name} ({len(subgraph.nodes())} nodes, {len(subgraph.edges())} edges)",
                    titlefont_size=14,
                    showlegend=False,
                    hovermode="closest",
                    margin=dict(b=20, l=5, r=5, t=40),
                    annotations=[
                        dict(
                            text="Drag to explore • Zoom to see details",
                            showarrow=False,
                            xref="paper",
                            yref="paper",
                            x=0.005,
                            y=-0.002,
                            xanchor="left",
                            yanchor="bottom",
                            font=dict(color="gray", size=10),
                        )
                    ],
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                ),
            )

            return fig, edge_info, node_info

        except Exception as e:
            log_message(f"Error in get_node_subnetwork: {str(e)}", message_type="error")
            return None

    def _generate_colors(self, n):
        return generate_colors(n)

    def export_network_data(
        self,
        show_har: bool = True,
        top_n_tfs: int = 10,
        top_n_targets: int = 10,
        specific_tfs: List[str] = None,
        specific_targets: List[str] = None,
        specific_stages: List[str] = None,
        specific_regions: List[str] = None,
        specific_celltypes: List[str] = None,
        specific_hars: List[str] = None,
        output_format: str = "csv",
    ):
        """
        Export filtered network data to CSV or JSON format

        Parameters:
        -----------
        output_format : str
            Export format ('csv' or 'json')

        Returns:
        --------
        str : Exported data as string
        """
        flow_data = self.network_data[
            (self.network_data["Weight"] > 0)
            & (self.network_data["TF"] != "Unknown_TF")
            & (self.network_data["Target"] != "Unknown_Target")
        ].copy()

        if specific_hars and self.har_tf_data is not None:
            relevant_tfs = self.har_tf_data[
                self.har_tf_data["HAR"].isin(specific_hars)
            ]["TF"].unique()
            flow_data = flow_data[flow_data["TF"].isin(relevant_tfs)]

        if specific_celltypes:
            flow_data = flow_data[flow_data["CellType"].isin(specific_celltypes)]
        if specific_regions:
            flow_data = flow_data[flow_data["Region"].isin(specific_regions)]
        if specific_stages:
            flow_data = flow_data[flow_data["Stage"].isin(specific_stages)]

        if specific_tfs:
            available_tfs = set(flow_data["TF"].unique())
            top_tfs = [tf for tf in specific_tfs if tf in available_tfs]
        else:
            top_tfs = self.get_top_nodes(
                "tf",
                top_n_tfs,
                {
                    "CellType": specific_celltypes,
                    "Region": specific_regions,
                    "Stage": specific_stages,
                },
            )

        if specific_targets:
            available_targets = set(flow_data["Target"].unique())
            top_targets = [
                target for target in specific_targets if target in available_targets
            ]
        else:
            top_targets = self.get_top_nodes(
                "target",
                top_n_targets,
                {
                    "CellType": specific_celltypes,
                    "Region": specific_regions,
                    "Stage": specific_stages,
                },
            )

        export_data = flow_data[
            flow_data["TF"].isin(top_tfs) & flow_data["Target"].isin(top_targets)
        ].copy()

        if show_har and self.har_tf_data is not None:
            har_mapping = self.har_tf_data[self.har_tf_data["TF"].isin(top_tfs)]
            export_data = export_data.merge(har_mapping, on="TF", how="left")

        if output_format == "csv":
            output = io.StringIO()
            export_data.to_csv(output, index=False)
            return output.getvalue()
        elif output_format == "json":
            return export_data.to_json(orient="records", indent=2)
        else:
            raise ValueError("output_format must be 'csv' or 'json'")

    def get_network_statistics(
        self,
        specific_tfs: List[str] = None,
        specific_targets: List[str] = None,
        specific_stages: List[str] = None,
        specific_regions: List[str] = None,
        specific_celltypes: List[str] = None,
        actual_tfs: List[str] = None,
        actual_targets: List[str] = None,
    ):
        """
        Get comprehensive network statistics

        Returns:
        --------
        dict : Network statistics
        """
        flow_data = self.network_data[
            (self.network_data["Weight"] > 0)
            & (self.network_data["TF"] != "Unknown_TF")
            & (self.network_data["Target"] != "Unknown_Target")
        ].copy()

        if specific_celltypes:
            flow_data = flow_data[flow_data["CellType"].isin(specific_celltypes)]
        if specific_regions:
            flow_data = flow_data[flow_data["Region"].isin(specific_regions)]
        if specific_stages:
            flow_data = flow_data[flow_data["Stage"].isin(specific_stages)]
        if actual_tfs:
            flow_data = flow_data[flow_data["TF"].isin(actual_tfs)]
        elif specific_tfs:
            flow_data = flow_data[flow_data["TF"].isin(specific_tfs)]
        if actual_targets:
            flow_data = flow_data[flow_data["Target"].isin(actual_targets)]
        elif specific_targets:
            flow_data = flow_data[flow_data["Target"].isin(specific_targets)]

        stats = {
            "total_connections": len(flow_data),
            "unique_tfs": flow_data["TF"].nunique(),
            "unique_targets": flow_data["Target"].nunique(),
            "unique_celltypes": flow_data["CellType"].nunique(),
            "unique_regions": flow_data["Region"].nunique(),
            "unique_stages": flow_data["Stage"].nunique(),
            "total_weight": flow_data["Weight"].sum(),
            "avg_weight": flow_data["Weight"].mean(),
            "max_weight": flow_data["Weight"].max(),
            "min_weight": flow_data["Weight"].min(),
            "top_tfs_by_weight": flow_data.groupby("TF")["Weight"]
            .apply(lambda x: x.abs().sum())
            .nlargest(10)
            .to_dict(),
            "top_targets_by_weight": flow_data.groupby("Target")["Weight"]
            .apply(lambda x: x.abs().sum())
            .nlargest(10)
            .to_dict(),
            "connections_by_celltype": flow_data.groupby("CellType").size().to_dict(),
            "connections_by_region": flow_data.groupby("Region").size().to_dict(),
            "connections_by_stage": flow_data.groupby("Stage").size().to_dict(),
        }

        return stats

    def get_region_stage_metrics(
        self,
        specific_tfs: List[str] = None,
        specific_targets: List[str] = None,
        specific_stages: List[str] = None,
        specific_regions: List[str] = None,
        specific_celltypes: List[str] = None,
        actual_tfs: List[str] = None,
        actual_targets: List[str] = None,
    ) -> pd.DataFrame:
        """
        Get Region x Stage aggregated counts (Edges, TFs, Genes) for current filters.
        Returns DataFrame with columns Region, Stage, Edges_count, TFs_count, Genes_count.
        """
        flow_data = self.network_data[
            (self.network_data["Weight"] > 0)
            & (self.network_data["TF"] != "Unknown_TF")
            & (self.network_data["Target"] != "Unknown_Target")
        ].copy()
        if specific_celltypes:
            flow_data = flow_data[flow_data["CellType"].isin(specific_celltypes)]
        if specific_regions:
            flow_data = flow_data[flow_data["Region"].isin(specific_regions)]
        if specific_stages:
            flow_data = flow_data[flow_data["Stage"].isin(specific_stages)]
        if actual_tfs:
            flow_data = flow_data[flow_data["TF"].isin(actual_tfs)]
        elif specific_tfs:
            flow_data = flow_data[flow_data["TF"].isin(specific_tfs)]
        if actual_targets:
            flow_data = flow_data[flow_data["Target"].isin(actual_targets)]
        elif specific_targets:
            flow_data = flow_data[flow_data["Target"].isin(specific_targets)]

        all_regions = sorted(flow_data["Region"].dropna().unique())
        all_stages = sorted(
            flow_data["Stage"].dropna().unique(),
            key=_stage_order_key,
        )
        if not all_regions or not all_stages:
            empty = pd.DataFrame(
                columns=["Region", "Stage", "Edges_count", "TFs_count", "Genes_count"]
            )
            return empty

        agg = (
            flow_data.groupby(["Region", "Stage"])
            .agg(
                Edges_count=("Weight", "count"),
                TFs_count=("TF", "nunique"),
                Genes_count=("Target", "nunique"),
            )
            .reset_index()
        )
        full = pd.DataFrame(
            [(r, s) for r in all_regions for s in all_stages],
            columns=["Region", "Stage"],
        )
        merged = full.merge(agg, on=["Region", "Stage"], how="left")
        merged["Edges_count"] = merged["Edges_count"].fillna(0).astype(int)
        merged["TFs_count"] = merged["TFs_count"].fillna(0).astype(int)
        merged["Genes_count"] = merged["Genes_count"].fillna(0).astype(int)
        return merged


def estimate_text_width(text: str, font_size: int = 10) -> float:
    """
    Estimate text width in pixels based on character count and font size.

    Parameters:
    -----------
    text : str
        The text to measure
    font_size : int
        Font size in pixels (default: 10)

    Returns:
    --------
    float : Estimated width in pixels
    """
    if not text:
        return 0.0
    avg_char_width = font_size * 0.6
    return len(str(text)) * avg_char_width


def capitalize_first_word(text: str) -> str:
    """Capitalize only the first word of a string"""
    if not text:
        return text
    words = text.split()
    if words:
        words[0] = words[0].capitalize()
    return " ".join(words)


def _stage_order_key(s: str) -> float:
    """Sort key for stage strings (e.g. S13 -> 13)."""
    try:
        return float(str(s).replace("S", ""))
    except (ValueError, TypeError):
        return float("inf")


def load_network_statistics_csv(network_dir: str) -> Optional[pd.DataFrame]:
    """Load network_statistics.csv from network_dir. Returns None if missing."""
    path = os.path.join(network_dir, "network_statistics.csv")
    if not os.path.isfile(path):
        return None
    df = pd.read_csv(path)
    for col in (
        "Region",
        "Stage",
        "Cell_type",
        "Edges_count",
        "TFs_count",
        "Genes_count",
    ):
        if col not in df.columns:
            return None
    return df


def build_single_region_stage_heatmap(
    region_stage_df: pd.DataFrame,
    metric: str,
    title: str,
    theme: str = "light",
) -> go.Figure:
    if region_stage_df is None or region_stage_df.empty:
        empty = go.Figure()
        empty.add_annotation(
            text="No data",
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.5,
            showarrow=False,
        )
        empty.update_layout(template="plotly_white", height=280, title=title)
        return empty

    all_regions = sorted(region_stage_df["Region"].dropna().unique())
    all_stages = sorted(
        region_stage_df["Stage"].dropna().unique(),
        key=_stage_order_key,
    )
    if metric not in region_stage_df.columns:
        metric = "Edges_count"
    pivot = region_stage_df.pivot_table(
        index="Region", columns="Stage", values=metric, aggfunc="first"
    )
    region_rev = list(reversed(all_regions))
    pivot = pivot.reindex(index=region_rev, columns=all_stages)
    z = pivot.values
    if z.size == 0:
        z = np.full((len(region_rev), len(all_stages)), np.nan)
    else:
        z = np.where(np.isnan(z) | (z == 0), np.nan, z)

    non_zero = region_stage_df[region_stage_df[metric] > 0][metric]
    if len(non_zero) > 0:
        zmin, zmax = float(non_zero.min()), float(non_zero.max())
    else:
        zmin, zmax = 0.0, 1.0

    n_regions = len(region_rev)
    n_stages = len(all_stages)
    margin_l = 50
    margin_r = 38
    margin_t = 24
    margin_b = 50
    max_plot_h = 420
    max_plot_w = 380
    plot_h = min(n_regions * 22, max_plot_h)
    plot_w = min(n_stages * 28, max_plot_w)
    fig_height = plot_h + margin_t + margin_b
    fig_width = plot_w + margin_l + margin_r
    fig_height = 550
    fig_width = 550

    fig = go.Figure(
        data=go.Heatmap(
            x=all_stages,
            y=region_rev,
            z=z,
            colorscale="Viridis",
            zmin=zmin,
            zmax=zmax,
            hovertemplate="Region: %{y}<br>Stage: %{x}<br>%{z}<extra></extra>",
            colorbar=dict(
                len=0.35,
                lenmode="fraction",
                y=0.5,
                yanchor="middle",
                thickness=7,
                outlinewidth=0,
            ),
        )
    )
    is_dark = theme == "dark"
    template = "plotly_dark" if is_dark else "plotly_white"
    fig.update_layout(
        title=dict(
            text=title,
            x=0,
            xanchor="left",
            font=dict(size=13),
        ),
        template=template,
        height=fig_height,
        width=fig_width,
        margin=dict(l=margin_l, r=margin_r, t=margin_t, b=margin_b),
        autosize=False,
        showlegend=False,
        xaxis=dict(
            tickangle=270,
            tickvals=list(range(n_stages)),
            ticktext=all_stages,
            title_text="Stage",
            automargin=True,
        ),
        yaxis=dict(
            autorange="reversed",
            tickvals=list(range(n_regions)),
            ticktext=region_rev,
            title_text="Brain region",
        ),
    )
    return fig


def parse_text_input(text: str) -> List[str]:
    """Parse comma-separated text input into a list of strings"""
    if not text:
        return []
    return [item.strip() for item in text.split(",") if item.strip()]


def merge_lists(list1: Optional[List[str]], list2: Optional[List[str]]) -> List[str]:
    """Merge two lists, removing duplicates"""
    result = set()
    if list1:
        result.update(list1)
    if list2:
        result.update(list2)
    return sorted(list(result))


def create_dash_app(visualizer: InteractiveNetworkVisualizer):
    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

    options = visualizer.get_available_options()

    app.index_string = """
    <!DOCTYPE html>
    <html>
        <head>
            {%metas%}
            <title>{%title%}</title>
            {%favicon%}
            {%css%}
            <style>
                :root {
                    --bg-color: #ffffff;
                    --text-color: #212529;
                    --card-bg: #ffffff;
                    --card-border: #dee2e6;
                    --primary-color: #0d6efd;
                    --secondary-color: #6c757d;
                }
                
                [data-theme="dark"] {
                    --bg-color: #121212;
                    --text-color: #ffffff;
                    --card-bg: #1e1e1e;
                    --card-border: #3a3a3a;
                    --primary-color: #4a9eff;
                    --secondary-color: #8e8e8e;
                }
                
                body {
                    background-color: var(--bg-color);
                    color: var(--text-color);
                    transition: background-color 0.2s ease, color 0.2s ease;
                }
                
                .card, .btn, .form-control, .form-select, .dropdown-menu, .alert, .btn-group {
                    box-shadow: none !important;
                }
                
                .btn-group .btn {
                    border-radius: 0;
                }
                
                .btn-group .btn:not(:first-child) {
                    margin-left: -1px;
                }
                
                .btn-group .btn:focus {
                    z-index: 1;
                }
                
                .alert {
                    border-radius: 0;
                    border-width: 1px;
                    border-style: solid;
                }
                
                [data-theme="dark"] body {
                    background-color: var(--bg-color);
                    color: var(--text-color);
                }
                
                [data-theme="dark"] .card, 
                [data-theme="dark"] .btn, 
                [data-theme="dark"] .form-control, 
                [data-theme="dark"] .form-select, 
                [data-theme="dark"] .dropdown-menu, 
                [data-theme="dark"] .alert {
                    box-shadow: none !important;
                }
                
                [data-theme="dark"] .text-muted {
                    color: #b0b0b0 !important;
                }
                
                [data-theme="dark"] label {
                    color: var(--text-color);
                }
                
                [data-theme="dark"] .card-body {
                    color: var(--text-color);
                }
                
                [data-theme="dark"] .card-header {
                    color: var(--text-color);
                }
                
                [data-theme="dark"] .dropdown-menu {
                    background-color: var(--card-bg);
                    border-color: var(--card-border);
                }
                
                [data-theme="dark"] .dropdown-item {
                    color: var(--text-color);
                }
                
                [data-theme="dark"] .dropdown-item:hover {
                    background-color: var(--card-border);
                    color: var(--text-color);
                }
                
                textarea {
                    border-radius: 0;
                    box-shadow: none;
                }
                
                [data-theme="dark"] textarea {
                    background-color: var(--card-bg);
                    color: var(--text-color);
                    border-color: var(--card-border);
                    border-radius: 0;
                    box-shadow: none;
                }
                
                [data-theme="dark"] textarea:focus {
                    background-color: var(--card-bg);
                    color: var(--text-color);
                    border-color: var(--primary-color);
                    box-shadow: none;
                    outline: 2px solid var(--primary-color);
                    outline-offset: -2px;
                }
                
                .card {
                    background-color: var(--card-bg);
                    border-color: var(--card-border);
                    border-radius: 0;
                    border-width: 1px;
                    border-style: solid;
                    box-shadow: none;
                    transition: none;
                }
                
                [data-theme="dark"] .card {
                    box-shadow: none;
                }
                
                .card-header {
                    background-color: var(--card-bg);
                    border-bottom: 1px solid var(--card-border);
                    border-width: 0 0 1px 0;
                    border-style: solid;
                    font-weight: 600;
                    padding: 1rem;
                    display: flex;
                    justify-content: space-between;
                    align-items: center;
                }
                
                .btn {
                    border-radius: 0;
                    font-weight: 500;
                    border-width: 1px;
                    border-style: solid;
                    transition: background-color 0.15s ease, border-color 0.15s ease, color 0.15s ease;
                }
                
                .btn:hover {
                    transform: none;
                    box-shadow: none;
                    opacity: 0.9;
                }
                
                .form-control, .form-select {
                    border-radius: 0;
                    border: 1px solid var(--card-border);
                    border-width: 1px;
                    border-style: solid;
                    transition: border-color 0.15s ease;
                    background-color: var(--card-bg);
                    color: var(--text-color);
                }
                
                .form-control:focus, .form-select:focus {
                    border-color: var(--primary-color);
                    box-shadow: none;
                    outline: 2px solid var(--primary-color);
                    outline-offset: -2px;
                    background-color: var(--card-bg);
                    color: var(--text-color);
                }
                
                .form-control::placeholder, .form-select::placeholder {
                    color: var(--secondary-color);
                    opacity: 0.7;
                }
                
                [data-theme="dark"] .form-control, [data-theme="dark"] .form-select {
                    background-color: var(--card-bg);
                    color: var(--text-color);
                }
                
                [data-theme="dark"] .form-control:focus, [data-theme="dark"] .form-select:focus {
                    background-color: var(--card-bg);
                    color: var(--text-color);
                }
                
                h1 {
                    font-weight: 600;
                    letter-spacing: 0;
                }
                
                .page-title {
                    font-size: 1.35rem;
                }
                
                .theme-switch-container {
                    display: inline-flex;
                    align-items: center;
                }
                
                .theme-switch-container .dropdown-toggle {
                    padding: 0;
                    border: none !important;
                    background: none !important;
                    box-shadow: none !important;
                    text-decoration: none !important;
                }
                
                .theme-switch-container .dropdown-toggle:hover,
                .theme-switch-container .dropdown-toggle:focus {
                    text-decoration: none !important;
                }
                
                .theme-switch-container .dropdown-toggle::after {
                    display: none !important;
                }
                
                .theme-icon-symbol {
                    font-size: 1.15em;
                    opacity: 0.9;
                }
                
                .theme-switch-container .dropdown-menu {
                    min-width: 120px;
                    border-radius: 0;
                    border: 1px solid var(--card-border);
                    border-width: 1px;
                    border-style: solid;
                    padding: 4px 0;
                    box-shadow: none;
                }
                
                .network-statistics-card-header {
                    padding: 0.35rem 1rem;
                }
                
                .network-stats-heatmap .js-plotly-plot,
                .network-stats-heatmap .plotly {
                    margin: 0 !important;
                    padding: 0 !important;
                }
                
                .network-stats-heatmap .svg-container {
                    margin: 0 !important;
                }
                
                .network-stats-heatmaps-row > [class*="col-"] {
                    padding-left: 0 !important;
                    padding-right: 0 !important;
                }
                
                .network-stats-heatmaps-row .network-stats-heatmap {
                    margin-left: -2px;
                    margin-right: -2px;
                }
                
                .network-stats-heatmaps-row .js-plotly-plot {
                    margin-left: 0 !important;
                    margin-right: 0 !important;
                }
                
                .theme-switch-container .dropdown-item {
                    padding: 8px 14px;
                    font-size: 14px;
                    display: flex;
                    align-items: center;
                }
                
                .theme-switch-container .dropdown-item.active {
                    background-color: #0d6efd;
                    color: #fff;
                }
                
                .theme-switch-container .dropdown-item:not(.active):hover {
                    background-color: rgba(0,0,0,0.05);
                }
                
                [data-theme="dark"] .theme-switch-container .dropdown-menu {
                    border-color: var(--card-border);
                }
                
                [data-theme="dark"] .theme-switch-container .dropdown-item:not(.active):hover {
                    background-color: var(--card-border);
                }
                
                /* Dash Dropdown styles for dark mode */
                [data-theme="dark"] .Select-control {
                    background-color: var(--card-bg) !important;
                    border-color: var(--card-border) !important;
                    color: var(--text-color) !important;
                }
                
                [data-theme="dark"] .Select-control:hover {
                    border-color: var(--primary-color) !important;
                }
                
                [data-theme="dark"] .Select-input > input {
                    color: var(--text-color) !important;
                }
                
                [data-theme="dark"] .Select-placeholder {
                    color: var(--secondary-color) !important;
                }
                
                [data-theme="dark"] .Select-value-label {
                    color: var(--text-color) !important;
                }
                
                [data-theme="dark"] .Select-menu-outer {
                    background-color: var(--card-bg) !important;
                    border-color: var(--card-border) !important;
                }
                
                [data-theme="dark"] .Select-option {
                    background-color: var(--card-bg) !important;
                    color: var(--text-color) !important;
                }
                
                [data-theme="dark"] .Select-option:hover,
                [data-theme="dark"] .Select-option:focus {
                    background-color: var(--card-border) !important;
                    color: var(--text-color) !important;
                }
                
                [data-theme="dark"] .Select-option.is-selected {
                    background-color: var(--primary-color) !important;
                    color: var(--text-color) !important;
                }
                
                [data-theme="dark"] .Select-option.is-focused {
                    background-color: var(--card-border) !important;
                    color: var(--text-color) !important;
                }
                
                [data-theme="dark"] .Select-multi-value-wrapper {
                    color: var(--text-color) !important;
                }
                
                [data-theme="dark"] .Select-value {
                    background-color: var(--card-border) !important;
                    border-color: var(--card-border) !important;
                    color: var(--text-color) !important;
                }
                
                [data-theme="dark"] .Select-value-label {
                    color: var(--text-color) !important;
                }
                
                [data-theme="dark"] .Select-value-icon {
                    border-color: var(--text-color) transparent transparent !important;
                }
                
                [data-theme="dark"] .Select-value-icon:hover {
                    background-color: var(--card-border) !important;
                }
                
                .navbar-fixed {
                    position: sticky;
                    top: 0;
                    z-index: 1000;
                    background-color: var(--bg-color);
                    border-bottom: 1px solid var(--card-border);
                    margin-left: -12px;
                    margin-right: -12px;
                    padding-left: 12px;
                    padding-right: 12px;
                }
                
                .nav-icon-btn {
                    transition: opacity 0.2s ease, transform 0.2s ease;
                }
                
                .nav-icon-btn:hover {
                    opacity: 0.8;
                    transform: scale(1.05);
                }
                
                .nav-icon-btn:active {
                    transform: scale(0.95);
                }
                
                #theme-dropdown-label {
                    transition: opacity 0.2s ease, transform 0.2s ease;
                }
                
                #theme-dropdown-label:hover {
                    opacity: 0.8;
                    transform: scale(1.05);
                }
                
                #theme-dropdown-label:active {
                    transform: scale(0.95);
                }
                
                .nav-icon-btn .dropdown-toggle {
                    padding: 0 !important;
                    border: none !important;
                    background: none !important;
                    box-shadow: none !important;
                }
                
                .nav-icon-btn .dropdown-toggle::after {
                    display: none !important;
                }
                
                #success-message-container {
                    position: relative;
                    margin-top: 10px;
                }
                
                #success-message-container .alert {
                    animation: fadeOut 0.5s ease-out 1s forwards;
                }
                
                @keyframes fadeOut {
                    from {
                        opacity: 1;
                    }
                    to {
                        opacity: 0;
                        height: 0;
                        padding: 0;
                        margin: 0;
                        overflow: hidden;
                    }
                }
            </style>
            <script>
                document.addEventListener('DOMContentLoaded', function() {
                    const observer = new MutationObserver(function(mutations) {
                        mutations.forEach(function(mutation) {
                            mutation.addedNodes.forEach(function(node) {
                                if (node.nodeType === 1 && node.id === 'success-message-container') {
                                    setTimeout(function() {
                                        if (node.parentNode) {
                                            node.style.transition = 'opacity 0.5s ease-out';
                                            node.style.opacity = '0';
                                            setTimeout(function() {
                                                if (node.parentNode) {
                                                    node.style.display = 'none';
                                                }
                                            }, 500);
                                        }
                                    }, 1000);
                                }
                            });
                        });
                    });
                    observer.observe(document.body, { childList: true, subtree: true });
                });
            </script>
        </head>
        <body>
            {%app_entry%}
            <footer>
                {%config%}
                {%scripts%}
                {%renderer%}
            </footer>
        </body>
    </html>
    """

    app.layout = dbc.Container(
        [
            dcc.Store(id="page-load", data=True),
            dcc.Store(id="theme-store", data="light"),
            dcc.Store(id="theme-mode-store", data="auto"),
            dcc.Interval(id="theme-interval", interval=60 * 1000, n_intervals=0),
            html.Div(id="_dummy-output", style={"display": "none"}),
            html.Div(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    html.Div(
                                        [
                                            html.H1(
                                                capitalize_first_word(
                                                    "Interactive network visualization"
                                                ),
                                                className="mb-0 page-title",
                                            ),
                                            html.Div(
                                                [
                                                    dbc.DropdownMenu(
                                                        [
                                                            dbc.DropdownMenuItem(
                                                                [
                                                                    html.Span(
                                                                        "⊙",
                                                                        className="me-2 theme-icon-symbol",
                                                                    ),
                                                                    "Light",
                                                                ],
                                                                id="theme-opt-light",
                                                            ),
                                                            dbc.DropdownMenuItem(
                                                                [
                                                                    html.Span(
                                                                        "◑",
                                                                        className="me-2 theme-icon-symbol",
                                                                    ),
                                                                    "Dark",
                                                                ],
                                                                id="theme-opt-dark",
                                                            ),
                                                            dbc.DropdownMenuItem(
                                                                [
                                                                    html.Span(
                                                                        "◐",
                                                                        className="me-2 theme-icon-symbol",
                                                                    ),
                                                                    "Auto",
                                                                ],
                                                                id="theme-opt-auto",
                                                            ),
                                                        ],
                                                        label=html.Div(
                                                            id="theme-dropdown-label",
                                                            style={
                                                                "display": "inline-flex",
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
                                                                "width": "30px",
                                                                "height": "30px",
                                                                "borderRadius": "50%",
                                                                "backgroundColor": "#ffffff",
                                                                "border": "1px solid #d0d7de",
                                                                "cursor": "pointer",
                                                                "userSelect": "none",
                                                                "flexShrink": "0",
                                                            },
                                                            children=html.Div(
                                                                id="theme-icon-inner",
                                                                style={
                                                                    "width": "20px",
                                                                    "height": "20px",
                                                                    "borderRadius": "50%",
                                                                    "background": "linear-gradient(to right, #0969da 50%, #ffffff 50%)",
                                                                    "border": "1px solid #d0d7de",
                                                                },
                                                            ),
                                                        ),
                                                        id="theme-dropdown",
                                                        color="link",
                                                        className="p-0 border-0 nav-icon-btn",
                                                        menu_variant="light",
                                                        align_end=True,
                                                    ),
                                                    html.Span(
                                                        "|",
                                                        className="nav-separator",
                                                        style={
                                                            "margin": "0 12px",
                                                            "color": "#d0d7de",
                                                            "fontWeight": "300",
                                                            "userSelect": "none",
                                                            "fontSize": "18px",
                                                        },
                                                    ),
                                                    html.A(
                                                        href="https://github.com/mengxu98/HARNexus",
                                                        id="nav-github-btn",
                                                        target="_blank",
                                                        rel="noopener noreferrer",
                                                        className="nav-icon-btn",
                                                        style={
                                                            "display": "inline-flex",
                                                            "alignItems": "center",
                                                            "justifyContent": "center",
                                                            "width": "30px",
                                                            "height": "30px",
                                                            "borderRadius": "50%",
                                                            "backgroundColor": "#374151",
                                                            "textDecoration": "none",
                                                            "border": "none",
                                                            "flexShrink": "0",
                                                        },
                                                        title="View on GitHub",
                                                        children=html.Img(
                                                            src="data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='20' height='20' viewBox='0 0 16 16' fill='white'%3E%3Cpath d='M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.012 8.012 0 0 0 16 8c0-4.42-3.58-8-8-8z'/%3E%3C/svg%3E",
                                                            alt="GitHub",
                                                            style={
                                                                "width": "20px",
                                                                "height": "20px",
                                                                "display": "block",
                                                            },
                                                        ),
                                                    ),
                                                ],
                                                className="theme-switch-container d-flex align-items-center",
                                                style={
                                                    "marginLeft": "auto",
                                                    "alignSelf": "center",
                                                },
                                            ),
                                        ],
                                        className="d-flex flex-nowrap align-items-baseline justify-content-between mb-2",
                                        style={"marginTop": "8px"},
                                    ),
                                    html.Hr(className="mb-0"),
                                ]
                            )
                        ]
                    ),
                ],
                id="navbar-container",
                className="navbar-fixed",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Card(
                                [
                                    dbc.CardHeader(
                                        capitalize_first_word("Filter controls")
                                    ),
                                    dbc.CardBody(
                                        [
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                capitalize_first_word(
                                                                    "Top n TFs:"
                                                                )
                                                            ),
                                                            dbc.Input(
                                                                id="top-n-tfs",
                                                                type="number",
                                                                value=10,
                                                                min=1,
                                                                max=50,
                                                                className="mb-3",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                capitalize_first_word(
                                                                    "Top n target genes:"
                                                                )
                                                            ),
                                                            dbc.Input(
                                                                id="top-n-targets",
                                                                type="number",
                                                                value=10,
                                                                min=1,
                                                                max=100,
                                                                className="mb-3",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                ]
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                capitalize_first_word(
                                                                    "Chart height:"
                                                                )
                                                            ),
                                                            dbc.Input(
                                                                id="chart-height",
                                                                type="number",
                                                                value=700,
                                                                min=300,
                                                                max=1000,
                                                                className="mb-3",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                capitalize_first_word(
                                                                    "Chart width:"
                                                                )
                                                            ),
                                                            dbc.Input(
                                                                id="chart-width",
                                                                type="number",
                                                                value=1200,
                                                                min=300,
                                                                max=2000,
                                                                className="mb-3",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                ]
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                capitalize_first_word(
                                                                    "Cell types:"
                                                                )
                                                            ),
                                                            dcc.Dropdown(
                                                                id="celltype-dropdown",
                                                                options=[
                                                                    {
                                                                        "label": ct,
                                                                        "value": ct,
                                                                    }
                                                                    for ct in options[
                                                                        "celltypes"
                                                                    ]
                                                                ],
                                                                multi=True,
                                                                placeholder="Select cell types...",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                capitalize_first_word(
                                                                    "Regions:"
                                                                )
                                                            ),
                                                            dcc.Dropdown(
                                                                id="region-dropdown",
                                                                options=[
                                                                    {
                                                                        "label": r,
                                                                        "value": r,
                                                                    }
                                                                    for r in options[
                                                                        "regions"
                                                                    ]
                                                                ],
                                                                multi=True,
                                                                placeholder="Select regions...",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                ],
                                                className="mt-3",
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                capitalize_first_word(
                                                                    "Stages:"
                                                                )
                                                            ),
                                                            dcc.Dropdown(
                                                                id="stage-dropdown",
                                                                options=[
                                                                    {
                                                                        "label": s,
                                                                        "value": s,
                                                                    }
                                                                    for s in options[
                                                                        "stages"
                                                                    ]
                                                                ],
                                                                multi=True,
                                                                placeholder="Select stages...",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                capitalize_first_word(
                                                                    "Show har level:"
                                                                )
                                                            ),
                                                            dbc.Switch(
                                                                id="show-har-switch",
                                                                value=True,
                                                                className="mb-3",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                ],
                                                className="mt-3",
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                capitalize_first_word(
                                                                    "Specific TFs:"
                                                                )
                                                            ),
                                                            dcc.Dropdown(
                                                                id="tf-dropdown",
                                                                options=[
                                                                    {
                                                                        "label": tf,
                                                                        "value": tf,
                                                                    }
                                                                    for tf in options[
                                                                        "tfs"
                                                                    ]
                                                                ],
                                                                multi=True,
                                                                placeholder="Select specific TFs...",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                capitalize_first_word(
                                                                    "Specific target genes:"
                                                                )
                                                            ),
                                                            dcc.Dropdown(
                                                                id="target-dropdown",
                                                                options=[
                                                                    {
                                                                        "label": t,
                                                                        "value": t,
                                                                    }
                                                                    for t in options[
                                                                        "targets"
                                                                    ]
                                                                ],
                                                                multi=True,
                                                                placeholder="Select specific target genes...",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                ],
                                                className="mt-3",
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                capitalize_first_word(
                                                                    "Specific TFs (manual input):"
                                                                )
                                                            ),
                                                            dcc.Textarea(
                                                                id="tf-text-input",
                                                                placeholder="Enter TFs separated by commas (e.g., IRF1, ZBTB33, ETS2)...",
                                                                style={
                                                                    "width": "100%",
                                                                    "height": "60px",
                                                                },
                                                                className="mb-2",
                                                            ),
                                                        ],
                                                        width=12,
                                                    ),
                                                ],
                                                className="mt-3",
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                capitalize_first_word(
                                                                    "Specific target genes (manual input):"
                                                                )
                                                            ),
                                                            dcc.Textarea(
                                                                id="target-text-input",
                                                                placeholder="Enter target genes separated by commas (e.g., DOCK10, B2M, MCF2L)...",
                                                                style={
                                                                    "width": "100%",
                                                                    "height": "60px",
                                                                },
                                                                className="mb-2",
                                                            ),
                                                        ],
                                                        width=12,
                                                    ),
                                                ],
                                                className="mt-3",
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                capitalize_first_word(
                                                                    "Specific HARs:"
                                                                )
                                                            ),
                                                            dcc.Dropdown(
                                                                id="har-dropdown",
                                                                options=[
                                                                    {
                                                                        "label": har,
                                                                        "value": har,
                                                                    }
                                                                    for har in options[
                                                                        "hars"
                                                                    ]
                                                                ]
                                                                if options.get("hars")
                                                                else [],
                                                                multi=True,
                                                                placeholder="Select specific HARs...",
                                                            ),
                                                        ],
                                                        width=12,
                                                    ),
                                                ],
                                                className="mt-3",
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Button(
                                                                capitalize_first_word(
                                                                    "Generate sankey"
                                                                ),
                                                                id="generate-button",
                                                                color="primary",
                                                                className="mt-3 w-100",
                                                            )
                                                        ],
                                                        width=6,
                                                    ),
                                                    dbc.Col(
                                                        [
                                                            dbc.Button(
                                                                capitalize_first_word(
                                                                    "Export csv"
                                                                ),
                                                                id="export-csv-button",
                                                                color="success",
                                                                className="mt-3 w-100",
                                                            )
                                                        ],
                                                        width=6,
                                                    ),
                                                ]
                                            ),
                                        ]
                                    ),
                                ]
                            )
                        ],
                        width=3,
                    ),
                    dbc.Col(
                        [
                            dbc.Card(
                                [
                                    dbc.CardHeader(
                                        [
                                            html.Span(
                                                capitalize_first_word(
                                                    "Network sankey visualization"
                                                ),
                                                style={"flex": "1"},
                                            ),
                                            dbc.ButtonGroup(
                                                [
                                                    dbc.Button(
                                                        "PDF",
                                                        id="export-pdf-btn",
                                                        color="primary",
                                                        size="sm",
                                                    ),
                                                    dbc.Button(
                                                        "PNG",
                                                        id="export-png-btn",
                                                        color="primary",
                                                        size="sm",
                                                    ),
                                                    dbc.Button(
                                                        "JPG",
                                                        id="export-jpg-btn",
                                                        color="primary",
                                                        size="sm",
                                                    ),
                                                    dbc.Button(
                                                        "SVG",
                                                        id="export-svg-btn",
                                                        color="primary",
                                                        size="sm",
                                                    ),
                                                ],
                                                className="ms-auto",
                                            ),
                                        ],
                                        style={
                                            "display": "flex",
                                            "alignItems": "center",
                                        },
                                    ),
                                    dbc.CardBody(
                                        [
                                            html.Div(
                                                [
                                                    dcc.Graph(
                                                        id="sankey-plot",
                                                        style={
                                                            "width": "100%",
                                                            "minWidth": "0",
                                                        },
                                                        config={
                                                            "responsive": True,
                                                            "displayModeBar": True,
                                                        },
                                                    ),
                                                ],
                                                id="sankey-container",
                                                style={
                                                    "width": "100%",
                                                    "overflowX": "auto",
                                                    "overflowY": "auto",
                                                },
                                            ),
                                            html.Div(id="click-info", className="mt-3"),
                                            html.Div(
                                                id="download-container",
                                                className="mt-3",
                                            ),
                                            html.Div(
                                                id="export-status", className="mt-2"
                                            ),
                                        ]
                                    ),
                                ]
                            )
                        ],
                        width=9,
                    ),
                ],
                className="mb-4",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Card(
                                [
                                    dbc.CardHeader(
                                        capitalize_first_word("Network statistics"),
                                        className="network-statistics-card-header",
                                    ),
                                    dbc.CardBody(
                                        [
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dcc.Graph(
                                                                id="heatmap-edges",
                                                                config={
                                                                    "responsive": False,
                                                                    "displayModeBar": False,
                                                                },
                                                                style={"width": "100%"},
                                                                className="network-stats-heatmap",
                                                            ),
                                                        ],
                                                        width=12,
                                                        md=4,
                                                        className="pe-0",
                                                    ),
                                                    dbc.Col(
                                                        [
                                                            dcc.Graph(
                                                                id="heatmap-tfs",
                                                                config={
                                                                    "responsive": False,
                                                                    "displayModeBar": False,
                                                                },
                                                                style={"width": "100%"},
                                                                className="network-stats-heatmap",
                                                            ),
                                                        ],
                                                        width=12,
                                                        md=4,
                                                        className="ps-0 pe-0",
                                                    ),
                                                    dbc.Col(
                                                        [
                                                            dcc.Graph(
                                                                id="heatmap-genes",
                                                                config={
                                                                    "responsive": False,
                                                                    "displayModeBar": False,
                                                                },
                                                                style={"width": "100%"},
                                                                className="network-stats-heatmap",
                                                            ),
                                                        ],
                                                        width=12,
                                                        md=4,
                                                        className="ps-0",
                                                    ),
                                                ],
                                                className="mb-3 g-0 network-stats-heatmaps-row",
                                            ),
                                            html.Div(id="network-stats"),
                                        ],
                                        style={
                                            "paddingTop": "0.25rem",
                                            "paddingBottom": "0.25rem",
                                        },
                                    ),
                                ]
                            )
                        ],
                        width=12,
                    )
                ],
                className="mt-3",
            ),
        ],
        fluid=True,
    )

    @app.callback(
        [
            Output("sankey-plot", "figure"),
            Output("sankey-plot", "style"),
            Output("click-info", "children"),
            Output("network-stats", "children"),
            Output("heatmap-edges", "figure"),
            Output("heatmap-tfs", "figure"),
            Output("heatmap-genes", "figure"),
        ],
        [Input("generate-button", "n_clicks"), Input("page-load", "data")],
        [
            State("show-har-switch", "value"),
            State("top-n-tfs", "value"),
            State("top-n-targets", "value"),
            State("celltype-dropdown", "value"),
            State("region-dropdown", "value"),
            State("stage-dropdown", "value"),
            State("tf-dropdown", "value"),
            State("tf-text-input", "value"),
            State("target-dropdown", "value"),
            State("target-text-input", "value"),
            State("har-dropdown", "value"),
            State("chart-height", "value"),
            State("chart-width", "value"),
            State("theme-store", "data"),
        ],
    )
    def update_sankey_plot(
        n_clicks,
        page_load,
        show_har,
        top_n_tfs,
        top_n_targets,
        celltypes,
        regions,
        stages,
        tfs_dropdown,
        tfs_text,
        targets_dropdown,
        targets_text,
        hars,
        height,
        width,
        theme,
    ):
        ctx = dash.callback_context
        if not ctx.triggered or ctx.triggered[0]["prop_id"] == "page-load.data":
            if show_har is None:
                show_har = True
            if top_n_tfs is None:
                top_n_tfs = 10
            if top_n_targets is None:
                top_n_targets = 10
            if height is None:
                height = 700
            if width is None:
                width = 1200
        elif n_clicks is None:
            return (
                go.Figure(),
                {"width": "100%", "minWidth": "0"},
                "",
                "",
                go.Figure(),
                go.Figure(),
                go.Figure(),
            )

        tfs_text_list = parse_text_input(tfs_text) if tfs_text else []
        tfs = merge_lists(tfs_dropdown, tfs_text_list)

        targets_text_list = parse_text_input(targets_text) if targets_text else []
        targets = merge_lists(targets_dropdown, targets_text_list)

        result = visualizer.create_sankey_diagram(
            show_har=show_har,
            top_n_tfs=top_n_tfs,
            top_n_targets=top_n_targets,
            specific_tfs=tfs if tfs else None,
            specific_targets=targets if targets else None,
            specific_stages=stages,
            specific_regions=regions,
            specific_celltypes=celltypes,
            specific_hars=hars if hars else None,
            height=height,
            width=width,
            theme=theme if theme else "light",
        )

        if isinstance(result, tuple) and len(result) >= 2:
            fig, message = result[0], result[1]
            actual_tfs = result[2] if len(result) > 2 else None
            actual_targets = result[3] if len(result) > 3 else None
        else:
            fig, message = result, "Success"
            actual_tfs = None
            actual_targets = None

        if fig is None:
            return (
                go.Figure(),
                {"width": "100%", "minWidth": "0"},
                dbc.Alert(message, color="warning"),
                "",
                go.Figure(),
                go.Figure(),
                go.Figure(),
            )

        stats = visualizer.get_network_statistics(
            specific_tfs=tfs if tfs else None,
            specific_targets=targets if targets else None,
            specific_stages=stages,
            specific_regions=regions,
            specific_celltypes=celltypes,
            actual_tfs=actual_tfs,
            actual_targets=actual_targets,
        )

        display_tf_count = top_n_tfs
        if actual_tfs:
            display_tf_count = min(len(actual_tfs), top_n_tfs)
        elif tfs:
            display_tf_count = min(len(tfs), top_n_tfs)

        display_target_count = top_n_targets
        if actual_targets:
            display_target_count = min(len(actual_targets), top_n_targets)
        elif targets:
            display_target_count = min(len(targets), top_n_targets)

        is_dark = theme == "dark" if theme else False
        template = "plotly_dark" if is_dark else "plotly_white"
        text_color = "#ffffff" if is_dark else "#212529"
        paper_color = "#1e1e1e" if is_dark else "#ffffff"
        plot_bg = "#1e1e1e" if is_dark else "#ffffff"
        grid_color = "#3a3a3a" if is_dark else "#e9ecef"

        summary_fig = go.Figure()
        summary_metrics = [
            ("Edges", stats["total_connections"]),
            ("TFs", stats["unique_tfs"]),
            ("Targets", stats["unique_targets"]),
            ("Cell Types", stats["unique_celltypes"]),
            ("Regions", stats["unique_regions"]),
            ("Stages", stats["unique_stages"]),
        ]
        summary_fig.add_trace(
            go.Bar(
                x=[m[0] for m in summary_metrics],
                y=[m[1] for m in summary_metrics],
                marker_color="rgb(25, 50, 200)",
                text=[f"{m[1]:,}" for m in summary_metrics],
                textposition="outside",
                textfont=dict(color=text_color, size=11),
                cliponaxis=False,
            )
        )
        summary_fig.update_layout(
            title=capitalize_first_word("Summary metrics"),
            height=400,
            showlegend=False,
            template=template,
            paper_bgcolor=paper_color,
            plot_bgcolor=plot_bg,
            font=dict(color=text_color),
            xaxis=dict(gridcolor=grid_color, tickangle=-45),
            yaxis=dict(gridcolor=grid_color),
            margin=dict(l=50, r=50, t=60, b=80),
        )

        top_tfs_items = sorted(
            stats["top_tfs_by_weight"].items(), key=lambda x: x[1], reverse=True
        )[:display_tf_count]
        if top_tfs_items:
            top_tfs_fig = go.Figure()
            weights = [weight for _, weight in top_tfs_items]
            tfs_list = [tf for tf, _ in top_tfs_items]
            top_tfs_fig.add_trace(
                go.Bar(
                    x=weights,
                    y=tfs_list,
                    orientation="h",
                    marker_color="rgb(205, 50, 65)",
                    text=[f"{w:.2f}" for w in weights],
                    textposition="outside",
                    textfont=dict(color=text_color, size=10),
                    cliponaxis=False,
                )
            )
            top_tfs_fig.update_layout(
                title=capitalize_first_word("Top TFs by weight"),
                height=400,
                xaxis_title="Weight",
                yaxis_title="TFs",
                showlegend=False,
                template=template,
                paper_bgcolor=paper_color,
                plot_bgcolor=plot_bg,
                font=dict(color=text_color),
                xaxis=dict(gridcolor=grid_color),
                yaxis=dict(gridcolor=grid_color, autorange="reversed"),
                margin=dict(l=100, r=80, t=60, b=50),
            )
        else:
            top_tfs_fig = go.Figure()
            top_tfs_fig.add_annotation(
                text="No TFs available",
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                showarrow=False,
            )
            top_tfs_fig.update_layout(
                height=400,
                showlegend=False,
                template=template,
                paper_bgcolor=paper_color,
                plot_bgcolor=plot_bg,
                font=dict(color=text_color),
            )

        top_targets_items = sorted(
            stats["top_targets_by_weight"].items(), key=lambda x: x[1], reverse=True
        )[:display_target_count]
        if top_targets_items:
            top_targets_fig = go.Figure()
            weights = [weight for _, weight in top_targets_items]
            targets_list = [target for target, _ in top_targets_items]
            top_targets_fig.add_trace(
                go.Bar(
                    x=weights,
                    y=targets_list,
                    orientation="h",
                    marker_color="rgb(44, 160, 44)",
                    text=[f"{w:.2f}" for w in weights],
                    textposition="outside",
                    textfont=dict(color=text_color, size=10),
                    cliponaxis=False,
                )
            )
            top_targets_fig.update_layout(
                title=capitalize_first_word("Top target genes by weight"),
                height=400,
                xaxis_title="Weight",
                yaxis_title="Target genes",
                showlegend=False,
                template=template,
                paper_bgcolor=paper_color,
                plot_bgcolor=plot_bg,
                font=dict(color=text_color),
                xaxis=dict(gridcolor=grid_color),
                yaxis=dict(gridcolor=grid_color, autorange="reversed"),
                margin=dict(l=100, r=80, t=60, b=50),
            )
        else:
            top_targets_fig = go.Figure()
            top_targets_fig.add_annotation(
                text="No target genes available",
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                showarrow=False,
            )
            top_targets_fig.update_layout(
                height=400,
                showlegend=False,
                template=template,
                paper_bgcolor=paper_color,
                plot_bgcolor=plot_bg,
                font=dict(color=text_color),
            )

        stats_html = [
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dcc.Graph(figure=summary_fig),
                        ],
                        width=12,
                        md=4,
                        lg=4,
                        xl=4,
                        className="pe-0",
                    ),
                    dbc.Col(
                        [
                            dcc.Graph(figure=top_tfs_fig),
                        ],
                        width=12,
                        md=4,
                        lg=4,
                        xl=4,
                        className="ps-0 pe-0",
                    ),
                    dbc.Col(
                        [
                            dcc.Graph(figure=top_targets_fig),
                        ],
                        width=12,
                        md=4,
                        lg=4,
                        xl=4,
                        className="ps-0",
                    ),
                ],
                className="mb-3 g-0",
            ),
        ]

        success_message = html.Div(
            [
                dbc.Alert(
                    capitalize_first_word("Sankey diagram generated successfully!"),
                    color="success",
                    className="mb-2",
                    style={
                        "fontSize": "12px",
                        "padding": "8px 12px",
                        "marginTop": "10px",
                    },
                )
            ],
            id="success-message-container",
        )

        chart_height = height if height else 600
        if fig and hasattr(fig, "layout") and hasattr(fig.layout, "height"):
            chart_height = fig.layout.height if fig.layout.height else chart_height

        graph_style = {
            "width": "100%",
            "minWidth": "0",
            "height": f"{chart_height}px",
        }

        region_stage_df = visualizer.get_region_stage_metrics(
            specific_tfs=tfs if tfs else None,
            specific_targets=targets if targets else None,
            specific_stages=stages,
            specific_regions=regions,
            specific_celltypes=celltypes,
            actual_tfs=actual_tfs,
            actual_targets=actual_targets,
        )
        theme_val = theme if theme else "light"
        heatmap_edges = build_single_region_stage_heatmap(
            region_stage_df,
            "Edges_count",
            "Edge counts",
            theme_val,
        )
        heatmap_tfs = build_single_region_stage_heatmap(
            region_stage_df, "TFs_count", "TF counts", theme_val
        )
        heatmap_genes = build_single_region_stage_heatmap(
            region_stage_df,
            "Genes_count",
            "Target gene counts",
            theme_val,
        )

        return (
            fig,
            graph_style,
            success_message,
            stats_html,
            heatmap_edges,
            heatmap_tfs,
            heatmap_genes,
        )

    @app.callback(
        Output("export-status", "children"),
        [
            Input("export-pdf-btn", "n_clicks"),
            Input("export-png-btn", "n_clicks"),
            Input("export-jpg-btn", "n_clicks"),
            Input("export-svg-btn", "n_clicks"),
        ],
        [State("sankey-plot", "figure")],
    )
    def export_chart(pdf_clicks, png_clicks, jpg_clicks, svg_clicks, figure):
        ctx = dash.callback_context
        if not ctx.triggered or not figure:
            return ""

        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        format_map = {
            "export-pdf-btn": "pdf",
            "export-png-btn": "png",
            "export-jpg-btn": "jpg",
            "export-svg-btn": "svg",
        }

        export_format = format_map.get(button_id)
        if not export_format:
            return ""

        try:
            import plotly.io as pio

            if isinstance(figure, dict):
                fig = go.Figure(figure)
            else:
                fig = figure

            img_bytes = pio.to_image(
                fig,
                format=export_format,
                width=fig.layout.width,
                height=fig.layout.height,
            )

            encoded_img = base64.b64encode(img_bytes).decode()
            mime_type_map = {
                "pdf": "application/pdf",
                "png": "image/png",
                "jpg": "image/jpeg",
                "svg": "image/svg+xml",
            }
            mime_type = mime_type_map[export_format]

            filename = f"sankey_diagram.{export_format}"

            return html.Div(
                [
                    html.A(
                        f"Download {export_format.upper()}",
                        id=f"download-{export_format}",
                        download=filename,
                        href=f"data:{mime_type};base64,{encoded_img}",
                        className="btn btn-success btn-sm",
                        style={"margin-right": "5px"},
                    )
                ]
            )
        except Exception as e:
            return dbc.Alert(
                f"Export failed: {str(e)}. Make sure kaleido is installed: pip install kaleido",
                color="danger",
            )

    @app.callback(
        Output("download-container", "children"),
        [Input("export-csv-button", "n_clicks")],
        [
            State("show-har-switch", "value"),
            State("top-n-tfs", "value"),
            State("top-n-targets", "value"),
            State("celltype-dropdown", "value"),
            State("region-dropdown", "value"),
            State("stage-dropdown", "value"),
            State("tf-dropdown", "value"),
            State("tf-text-input", "value"),
            State("target-dropdown", "value"),
            State("target-text-input", "value"),
            State("har-dropdown", "value"),
        ],
    )
    def export_data(
        n_clicks,
        show_har,
        top_n_tfs,
        top_n_targets,
        celltypes,
        regions,
        stages,
        tfs_dropdown,
        tfs_text,
        targets_dropdown,
        targets_text,
        hars,
    ):
        if n_clicks is None:
            return ""

        tfs_text_list = parse_text_input(tfs_text) if tfs_text else []
        tfs = merge_lists(tfs_dropdown, tfs_text_list)

        targets_text_list = parse_text_input(targets_text) if targets_text else []
        targets = merge_lists(targets_dropdown, targets_text_list)

        try:
            csv_data = visualizer.export_network_data(
                show_har=show_har,
                top_n_tfs=top_n_tfs,
                top_n_targets=top_n_targets,
                specific_tfs=tfs if tfs else None,
                specific_targets=targets if targets else None,
                specific_stages=stages,
                specific_regions=regions,
                specific_celltypes=celltypes,
                specific_hars=hars if hars else None,
                output_format="csv",
            )

            encoded_data = base64.b64encode(csv_data.encode()).decode()

            return html.Div(
                [
                    html.A(
                        capitalize_first_word("Download CSV file"),
                        id="download-link",
                        download="network_data.csv",
                        href=f"data:text/csv;base64,{encoded_data}",
                        className="btn btn-success",
                    )
                ]
            )
        except Exception as e:
            return dbc.Alert(f"Export failed: {str(e)}", color="danger")

    @app.callback(
        [
            Output("theme-store", "data"),
            Output("theme-mode-store", "data"),
            Output("theme-dropdown-label", "children"),
            Output("theme-dropdown-label", "style"),
            Output("theme-opt-light", "className"),
            Output("theme-opt-dark", "className"),
            Output("theme-opt-auto", "className"),
        ],
        [
            Input("theme-opt-light", "n_clicks"),
            Input("theme-opt-dark", "n_clicks"),
            Input("theme-opt-auto", "n_clicks"),
            Input("theme-interval", "n_intervals"),
            Input("page-load", "data"),
        ],
        [
            State("theme-store", "data"),
            State("theme-mode-store", "data"),
        ],
    )
    def update_theme(
        _n_light, _n_dark, _n_auto, _n_interval, _page_load, theme, theme_mode
    ):
        from dash import ctx

        theme = theme or "light"
        theme_mode = theme_mode or "auto"
        triggered = getattr(ctx, "triggered_id", None)

        if triggered == "theme-opt-light":
            theme = "light"
            theme_mode = "light"
        elif triggered == "theme-opt-dark":
            theme = "dark"
            theme_mode = "dark"
        elif triggered == "theme-opt-auto":
            theme_mode = "auto"
            hour = datetime.now().hour
            theme = "dark" if hour < 6 or hour >= 18 else "light"
        elif triggered == "page-load.data" and theme_mode == "auto":
            hour = datetime.now().hour
            theme = "dark" if hour < 6 or hour >= 18 else "light"
        elif triggered == "theme-interval" and theme_mode == "auto":
            hour = datetime.now().hour
            theme = "dark" if hour < 6 or hour >= 18 else "light"

        if theme_mode == "light":
            icon_bg = "#ffffff"
            icon_inner_bg = "#0969da"
            icon_border = "#d0d7de"
        elif theme_mode == "dark":
            icon_bg = "#24292f"
            icon_inner_bg = "#0969da"
            icon_border = "#24292f"
        else:
            icon_bg = "#ffffff"
            icon_inner_bg = "linear-gradient(to right, #0969da 50%, #ffffff 50%)"
            icon_border = "#d0d7de"

        label_icon = html.Div(
            id="theme-icon-inner",
            style={
                "width": "20px",
                "height": "20px",
                "borderRadius": "50%",
                "background": icon_inner_bg,
                "border": f"1px solid {icon_border}",
            },
        )

        label_style = {
            "display": "inline-flex",
            "alignItems": "center",
            "justifyContent": "center",
            "width": "30px",
            "height": "30px",
            "borderRadius": "50%",
            "backgroundColor": icon_bg,
            "border": f"1px solid {icon_border}",
            "cursor": "pointer",
            "userSelect": "none",
            "flexShrink": "0",
        }

        active = "dropdown-item active"
        inactive = "dropdown-item"
        cl_light = active if theme_mode == "light" else inactive
        cl_dark = active if theme_mode == "dark" else inactive
        cl_auto = active if theme_mode == "auto" else inactive

        return theme, theme_mode, label_icon, label_style, cl_light, cl_dark, cl_auto

    app.clientside_callback(
        """
        function(theme) {
            if (theme) {
                document.documentElement.setAttribute('data-theme', theme);
            }
            return window.dash_clientside.no_update;
        }
        """,
        Output("_dummy-output", "children"),
        Input("theme-store", "data"),
    )

    return app


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Enhanced Network Visualization Tool")
    parser.add_argument(
        "--data-path",
        type=str,
        default=None,
        help="Direct path to network data CSV file (optional, if not provided, will use smart loading)",
    )
    parser.add_argument(
        "--har-tf-path",
        type=str,
        default="results/har_tf/human/har_tf_pairs_scores.csv",
        help="Path to HAR-TF combinations CSV file (default: results/har_tf/human/har_tf_pairs_scores.csv)",
    )
    parser.add_argument(
        "--base-dir",
        type=str,
        default="results/networks/analysis",
        help="Base directory for network data (default: results/networks/analysis)",
    )
    parser.add_argument(
        "--regions",
        type=str,
        nargs="+",
        default=None,
        help="Brain regions to filter (e.g., 'Cerebral cortex' 'Prefrontal cortex')",
    )
    parser.add_argument(
        "--stages",
        type=str,
        nargs="+",
        default=None,
        help="Developmental stages to filter (e.g., S13 S14 S15)",
    )
    parser.add_argument(
        "--celltypes",
        type=str,
        nargs="+",
        default=None,
        help="Cell types to filter (e.g., 'Astrocytes' 'Excitatory neurons')",
    )
    parser.add_argument(
        "--tfs",
        type=str,
        nargs="+",
        default=None,
        help="Transcription factors to filter",
    )
    parser.add_argument(
        "--targets", type=str, nargs="+", default=None, help="Target genes to filter"
    )
    parser.add_argument(
        "--hars",
        type=str,
        nargs="+",
        default=None,
        help="HARs to filter (requires --har-tf-path)",
    )
    parser.add_argument(
        "--port", type=int, default=8000, help="Port for web server (default: 8000)"
    )

    args = parser.parse_args()

    try:
        log_message("Initializing Enhanced Network Visualizer...", message_type="info")

        visualizer = InteractiveNetworkVisualizer(
            data_path=args.data_path,
            har_tf_path=args.har_tf_path,
            network_dir=args.base_dir,
            regions=args.regions,
            stages=args.stages,
            celltypes=args.celltypes,
            tfs=args.tfs,
            targets=args.targets,
            hars=args.hars,
        )

        log_message("Creating Dash application...", message_type="info")
        app = create_dash_app(visualizer)

        log_message("Starting web server...", message_type="info")
        log_message(
            f"Open your browser and go to: http://127.0.0.1:{args.port}",
            message_type="success",
        )

        app.run(debug=True, host="127.0.0.1", port=args.port)

    except Exception as e:
        log_message(f"Error: {str(e)}", message_type="error")
        import traceback

        traceback.print_exc()
        sys.exit(1)
