import pandas as pd
import plotly.graph_objects as go
from pathlib import Path
import sys
import os
from typing import List, Optional
from tqdm import tqdm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import (
    COLOR_CELLTYPES,
    COLOR_STAGES,
    DEFAULT_NODE_COLOR,
    log_message,
)
from functions.utils_network import load_network_data, load_har_tf_data, generate_colors


def load_network_data_wrapper(
    file_path: Optional[str] = None,
    network_dir: str = "results/networks/har_csn_atlas",
    regions: Optional[List[str]] = None,
    stages: Optional[List[str]] = None,
    celltypes: Optional[List[str]] = None,
    tfs: Optional[List[str]] = None,
    targets: Optional[List[str]] = None,
    hars: Optional[List[str]] = None,
    har_tf_data: Optional[pd.DataFrame] = None,
):
    try:
        if file_path:
            log_message(f"Loading data from {file_path}...", message_type="info")
            data = pd.read_csv(file_path)
        else:
            return load_network_data(
                network_dir=network_dir,
                regions=regions,
                stages=stages,
                celltypes=celltypes,
                tfs=tfs,
                targets=targets,
                hars=hars,
                har_tf_data=har_tf_data,
            )

        log_message("Processing data...", message_type="info")
        with tqdm(total=6, desc="Cleaning data", ascii=True) as pbar:
            if "Weight" not in data.columns:
                if "weight" in data.columns:
                    data = data.rename(columns={"weight": "Weight"})
                else:
                    raise ValueError("Weight column not found in data")

            if "TF" not in data.columns:
                if "regulator" in data.columns:
                    data = data.rename(columns={"regulator": "TF"})
                else:
                    raise ValueError("TF column not found in data")

            if "Target" not in data.columns:
                if "target" in data.columns:
                    data = data.rename(columns={"target": "Target"})
                else:
                    raise ValueError("Target column not found in data")

            data["Weight"] = data["Weight"].fillna(0)
            pbar.update(1)
            data["TF"] = data["TF"].fillna("Unknown_TF")
            pbar.update(1)
            data["Target"] = data["Target"].fillna("Unknown_Target")
            pbar.update(1)

            if "Stage" not in data.columns:
                data["Stage"] = "Unknown_Stage"
            if "Region" not in data.columns:
                data["Region"] = "Unknown_Region"
            if "CellType" not in data.columns:
                data["CellType"] = "Unknown_CellType"

            data["Stage"] = data["Stage"].fillna("Unknown_Stage")
            pbar.update(1)
            data["Region"] = data["Region"].fillna("Unknown_Region")
            pbar.update(1)
            data["CellType"] = data["CellType"].fillna("Unknown_CellType")
            pbar.update(1)

        data = data.dropna()
        log_message("Data loaded successfully.", message_type="success")
        return data
    except Exception as e:
        log_message(f"Error loading data: {str(e)}", message_type="error")
        log_message(
            "Please ensure you have run the network processing script first:",
            message_type="info",
        )
        log_message("Rscript step01-network_processing.R", message_type="info")
        log_message(
            "Or provide your own data in the correct format:", message_type="info"
        )
        log_message("TF, Target, Weight, CellType, Stage, Region", message_type="info")
        log_message("The required columns are:", message_type="info")
        log_message("TF, Target, Weight", message_type="info")
        log_message("The optional columns are:", message_type="info")
        log_message("CellType, Stage, Region", message_type="info")
        sys.exit(1)


def create_network_sankey(
    data_path: Optional[str] = None,
    output_path: str = "figures/networks/sankey/",
    fig_name: Optional[str] = None,
    har_tf_path: Optional[str] = None,
    network_dir: str = "results/networks/analysis",
    top_n_tfs: int = 5,
    top_n_targets: int = 10,
    specific_tfs: Optional[List[str]] = None,
    specific_targets: Optional[List[str]] = None,
    specific_stages: Optional[List[str]] = None,
    specific_regions: Optional[List[str]] = None,
    specific_celltypes: Optional[List[str]] = None,
    specific_hars: Optional[List[str]] = None,
    only_show_levels: Optional[List[str]] = None,
    font_size: int = 10,
    font_family: str = "Arial",
    font_color: str = "black",
    height: int = 600,
    width: int = 900,
    output_formats: Optional[List[str]] = None,
    show_background: bool = False,
    background_color: str = "aliceblue",
):
    har_tf_data = load_har_tf_data(har_tf_path)
    data = load_network_data_wrapper(
        file_path=data_path,
        network_dir=network_dir,
        regions=specific_regions,
        stages=specific_stages,
        celltypes=specific_celltypes,
        tfs=specific_tfs,
        targets=specific_targets,
        hars=specific_hars,
        har_tf_data=har_tf_data,
    )

    log_message("Creating Sankey diagram...", message_type="info")
    with tqdm(total=7, desc="Generating visualization", ascii=True) as pbar:
        flow_data = data[
            (data["Weight"] > 0)
            & (data["TF"] != "Unknown_TF")
            & (data["Target"] != "Unknown_Target")
        ].copy()
        pbar.update(1)

        if only_show_levels is None:
            only_show_levels = []

        if specific_celltypes:
            if isinstance(specific_celltypes, str):
                flow_data = flow_data[flow_data["CellType"] == specific_celltypes]
            else:
                flow_data = flow_data[flow_data["CellType"].isin(specific_celltypes)]

        if specific_regions:
            flow_data = flow_data[flow_data["Region"].isin(specific_regions)]

        if specific_stages:
            flow_data = flow_data[flow_data["Stage"].isin(specific_stages)]

        pbar.update(1)

        if len(flow_data) > 0:

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

            if specific_tfs is not None:
                available_tfs = set(flow_data["TF"].unique())
                top_tfs = [tf for tf in specific_tfs if tf in available_tfs]
                if not top_tfs:
                    log_message(
                        "None of the specified TFs found in the data.",
                        message_type="warning",
                    )
                    return None
                if len(top_tfs) < len(specific_tfs):
                    missing_tfs = set(specific_tfs) - set(top_tfs)
                    log_message(
                        f"Warning: Some specified TFs not found: {missing_tfs}",
                        message_type="warning",
                    )
            else:
                if specific_targets:
                    # Strategy: guarantee each specified target has at least one TF->Target edge
                    # in the selected TF set, then (optionally) fill remaining slots by TF strength.
                    target_data = flow_data[flow_data["Target"].isin(specific_targets)]
                    if len(target_data) == 0:
                        log_message(
                            "No specified target genes found in the data after filtering.",
                            message_type="warning",
                        )
                        return None

                    # Compute TF-target edge strength on the filtered data
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

                    # Greedy set cover to ensure all specified targets are represented
                    targets_in_data = set(tf_target_strength["Target"].unique())
                    missing_targets = set(specific_targets) & targets_in_data
                    missing_no_edges = set(specific_targets) - targets_in_data

                    cover_map = (
                        tf_target_strength.groupby("TF")["Target"].agg(set).to_dict()
                    )

                    selected_tfs: List[str] = []
                    selected_set = set()
                    # Do NOT add a new parameter: derive an effective cap from existing args
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

                    if missing_no_edges:
                        log_message(
                            f"Warning: Some specified targets have no edges after filtering: {missing_no_edges}",
                            message_type="warning",
                        )
                    if missing_targets:
                        log_message(
                            f"Warning: Could not cover all specified targets with <= {effective_top_n_tfs} TFs. Still missing: {missing_targets}",
                            message_type="warning",
                        )

                    # Fill remaining TF slots (up to top_n_tfs) by overall TF strength
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

            har_nodes = []
            har_to_tf_map = {}
            har_tf_details = []

            if har_tf_data is not None:
                har_tf_groups = (
                    har_tf_data[har_tf_data["TF"].isin(top_tfs)]
                    .groupby("TF")["HAR"]
                    .agg(lambda x: sorted(x))
                    .reset_index()
                )

                for _, row in har_tf_groups.iterrows():
                    tf = row["TF"]
                    har_list = row["HAR"]
                    first_har = har_list[0]
                    har_count = len(har_list)
                    har_display = f"HARs: {first_har}...({har_count})"
                    har_full_list = ", ".join(har_list)
                    har_tf_details.append({"TF": tf, "HARs": har_full_list})
                    har_nodes.append(har_display)
                    har_to_tf_map[har_display] = tf

                har_tf_details_df = pd.DataFrame(har_tf_details)

                filename_parts = ["network"]
                if specific_celltypes:
                    cell_type_str = (
                        "_".join(specific_celltypes)
                        if isinstance(specific_celltypes, list)
                        else specific_celltypes
                    )
                    filename_parts.append(f"cell_{cell_type_str}")
                if specific_regions:
                    region_str = "_".join(specific_regions)
                    filename_parts.append(f"region_{region_str}")
                if specific_stages:
                    stage_str = "_".join(specific_stages)
                    filename_parts.append(f"stage_{stage_str}")
                if specific_tfs is not None:
                    filename_parts.append("specific_tfs")
                if specific_targets is not None:
                    filename_parts.append("specific_targets")

                if fig_name is None:
                    fig_name = "_".join(filename_parts)

                har_tf_details_df.to_csv(
                    f"{output_path}{fig_name}_har_tf_mappings.csv", index=False
                )

            flow_data_filtered = flow_data[flow_data["TF"].isin(top_tfs)]
            if specific_targets is not None:
                available_targets = set(flow_data_filtered["Target"].unique())
                top_targets = [
                    target for target in specific_targets if target in available_targets
                ]
                if not top_targets:
                    log_message(
                        "No specified target genes found in the data for selected TFs.",
                        message_type="warning",
                    )
                    return None
                if len(top_targets) < len(specific_targets):
                    missing_targets = set(specific_targets) - set(top_targets)
                    log_message(
                        f"Warning: Some specified targets not found: {missing_targets}",
                        message_type="warning",
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
            pbar.update(1)

            flow_data = flow_data[
                flow_data["TF"].isin(top_tfs) & flow_data["Target"].isin(top_targets)
            ].copy()

            if len(flow_data) == 0:
                log_message(
                    "No valid connections found for the specified parameters.",
                    message_type="warning",
                )
                return None

            levels_to_show = ["TF"]

            if har_tf_data is not None:
                levels_to_show.insert(0, "HAR")

            stage_region_celltype = ["Stage", "Region", "CellType"]
            if only_show_levels:
                levels_to_show.extend(
                    L
                    for L in stage_region_celltype
                    if L in only_show_levels
                    or (L == "CellType" and "specific_celltypes" in only_show_levels)
                )
            else:
                levels_to_show.extend(stage_region_celltype)

            levels_to_show.append("Target")

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

            level_nodes = {
                "HAR": har_nodes if har_tf_data is not None else [],
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

            nodes = []
            node_labels = []  # display labels: TF and Target in italic
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

            final_width = width if width else 900

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
                                    mask &= (
                                        flow_data[current_level] == row[current_level]
                                    )
                                path_tfs = flow_data[mask]["TF"].unique()
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
            colors = generate_colors(n_tfs)
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
                            node_colors.append(
                                COLOR_STAGES.get(node, DEFAULT_NODE_COLOR)
                            )
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
                                dominant_tf = max(
                                    tf_weights.items(), key=lambda x: x[1]
                                )[0]
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
            level_to_x = {
                levels_to_show[i]: level_x[i] for i in range(len(levels_to_show))
            }

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
            pbar.update(1)

            title_parts = ["Network wiring"]
            if specific_celltypes:
                cell_type_str = (
                    specific_celltypes
                    if isinstance(specific_celltypes, str)
                    else specific_celltypes[0]
                )
                title_parts.append(f"of {cell_type_str}")
            if specific_regions:
                title_parts.append(f"in {', '.join(specific_regions)}")
            if specific_stages:
                title_parts.append(f"at {', '.join(specific_stages)}")

            title_text = " ".join(title_parts)

            if specific_tfs is not None:
                title_text += f"<br>Selected TFs: {', '.join(top_tfs)}"
            if specific_targets is not None:
                title_text += f"<br>Selected Targets: {', '.join(top_targets)}"

            level_labels = {
                "HAR": "HARs",
                "CellType": "Celltypes",
                "Stage": "Stages",
                "Region": "Regions",
                "TF": "TFs",
                "Target": "Target genes",
            }
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
                    font=dict(
                        size=font_size,
                        color=font_color,
                        family=font_family if font_family else "Arial",
                    ),
                )
                for level in levels_to_show
            ]

            _plot_bg = background_color if show_background else "white"
            fig.update_layout(
                title=dict(
                    text=title_text,
                    x=0.1,
                    y=0.9,
                    xanchor="left",
                    yanchor="top",
                    font=dict(
                        size=font_size,
                        color=font_color,
                        family=font_family if font_family else "Arial",
                    ),
                ),
                font=dict(
                    size=font_size,
                    color=font_color,
                    family=font_family if font_family else "Arial",
                ),
                height=height,
                width=final_width,
                margin=dict(b=52, l=80, r=80, t=80),
                autosize=False,
                annotations=col_annotations,
                xaxis=dict(visible=False, range=[0, 1]),
                yaxis=dict(visible=False),
                plot_bgcolor=_plot_bg,
                paper_bgcolor="white",
            )

            fig.update_traces(
                textfont=dict(
                    size=font_size,
                    color=font_color,
                    family=font_family or "Arial",
                ),
                selector=dict(type="sankey"),
            )

            Path(output_path).mkdir(parents=True, exist_ok=True)

            log_message("Saving outputs...", message_type="info")
            fig.write_html(f"{output_path}{fig_name}.html")

            if output_formats is None:
                output_formats = ["pdf"]

            for fmt in output_formats:
                try:
                    fig.write_image(f"{output_path}{fig_name}.{fmt}")
                    log_message(
                        f"Saved {fmt.upper()} version: {fig_name}.{fmt}",
                        message_type="info",
                    )
                except Exception as e:
                    log_message(
                        f"Warning: Could not save {fmt.upper()} version: {str(e)}",
                        message_type="warning",
                    )
            pbar.update(1)

            log_message("Visualization completed successfully.", message_type="success")
            return fig
        else:
            log_message(
                "No data found for the specified parameters.", message_type="warning"
            )
            return None


if __name__ == "__main__":
    try:
        Path("figures/networks/sankey").mkdir(parents=True, exist_ok=True)
        network_dir = "results/networks/analysis"
        har_tf_path = "results/har_tf/human/har_tf_pairs_scores.csv"
        output_path = "figures/networks/sankey/"

        create_network_sankey(
            network_dir=network_dir,
            har_tf_path=har_tf_path,
            output_path=output_path,
            fig_name="har_csn_atlas",
            top_n_tfs=10,
            top_n_targets=20,
            height=700,
            width=1200,
            output_formats=["pdf", "svg"],
        )
        create_network_sankey(
            network_dir=network_dir,
            har_tf_path=har_tf_path,
            output_path=output_path,
            fig_name="har_csn_atlas_2",
            top_n_tfs=10,
            top_n_targets=10,
            height=600,
            width=1000,
            output_formats=["pdf", "svg", "png"],
        )
        create_network_sankey(
            network_dir=network_dir,
            har_tf_path=har_tf_path,
            output_path=output_path,
            fig_name="GPR89B",
            top_n_tfs=20,
            specific_targets=["GPR89B"],
            height=420,
            width=1100,
            output_formats=["pdf", "svg"],
        )

        create_network_sankey(
            network_dir=network_dir,
            har_tf_path=har_tf_path,
            output_path=output_path,
            fig_name="GPR89B-2",
            top_n_tfs=10,
            specific_targets=[
                "GPR89B",
                "FRMPD2B",
                "SRGAP2C",
                "NOTCH2NL",
                "ARHGAP11B",
                "TBC1D3",
                "CROCCP2",
                "LRRC37B",
            ],
            height=400,
            width=1100,
            output_formats=["pdf", "svg"],
        )

        create_network_sankey(
            network_dir=network_dir,
            har_tf_path=har_tf_path,
            output_path="figures/pfc_astrocytes/sankey/",
            fig_name="astrocytes_genes_all_stages",
            specific_regions=["Prefrontal cortex"],
            # specific_stages=["S6", "S7", "S11", "S12"],
            specific_celltypes=["Astrocytes"],
            top_n_tfs=10,
            specific_targets=[
                "CREB5",
                "ADAM9",
                "VIM",
                "EMP1",
                "SEMA5A",
                "WEE1",
                "VCAN",
                "NLGN1",
                "ADAMTS6",
                "PTPRZ1",
                # "ADAM9",
                # "CREB5",
                # "EEPD1",
                # "EMP1",
                # "FABP7",
                # "H1-0",
                # "HSPA1B",
                # "IQGAP2",
                # "PTN",
                # "PTPRZ1",
                # "TFPI",
                # "VIM",
                # "WEE1",
                # "QKI",
                # "SEMA5B",
                # # "ADAMTS6",
                # "LRRTM3",
                # # "NLGN1",
                # "GPC6",
                # "PTPRT",
                # # "SEMA5A",
                # "KCND2",
                # "DCLK2",
                # "VCAN",
            ],
            height=300,
            width=800,
            output_formats=["pdf", "svg"],
        )

        create_network_sankey(
            network_dir=network_dir,
            har_tf_path=har_tf_path,
            output_path="figures/pfc_astrocytes/sankey/",
            fig_name="astrocytes_genes",
            specific_regions=["Prefrontal cortex"],
            specific_stages=["S6", "S7", "S11", "S12"],
            specific_celltypes=["Astrocytes"],
            top_n_tfs=10,
            specific_targets=[
                "CREB5",
                "ADAM9",
                "VIM",
                "EMP1",
                "SEMA5A",
                "WEE1",
                "VCAN",
                "NLGN1",
                "ADAMTS6",
                "PTPRZ1",
                # "ADAM9",
                # "CREB5",
                # "EEPD1",
                # "EMP1",
                # "FABP7",
                # "H1-0",
                # "HSPA1B",
                # "IQGAP2",
                # "PTN",
                # "PTPRZ1",
                # "TFPI",
                # "VIM",
                # "WEE1",
                # "QKI",
                # "SEMA5B",
                # # "ADAMTS6",
                # "LRRTM3",
                # # "NLGN1",
                # "GPC6",
                # "PTPRT",
                # # "SEMA5A",
                # "KCND2",
                # "DCLK2",
                # "VCAN",
            ],
            height=300,
            width=800,
            output_formats=["pdf", "svg"],
        )

        create_network_sankey(
            network_dir=network_dir,
            har_tf_path=har_tf_path,
            output_path="figures/pfc_astrocytes/sankey/",
            fig_name="astrocytes_genes2",
            specific_regions=["Prefrontal cortex"],
            specific_stages=["S6", "S7", "S11", "S12"],
            specific_celltypes=["Astrocytes"],
            only_show_levels=["Stage"],
            top_n_tfs=10,
            specific_targets=[
                "CREB5",
                "ADAM9",
                "VIM",
                "EMP1",
                "SEMA5A",
                "WEE1",
                "VCAN",
                "NLGN1",
                "ADAMTS6",
                "PTPRZ1",
            ],
            height=300,
            width=600,
            output_formats=["pdf", "svg"],
        )

        create_network_sankey(
            network_dir=network_dir,
            har_tf_path=har_tf_path,
            output_path="figures/pfc_astrocytes/sankey/",
            fig_name="astrocytes_subcelltype_genes",
            specific_regions=["Prefrontal cortex"],
            specific_celltypes=["Astrocytes"],
            top_n_tfs=10,
            specific_targets=[
                # Early / progenitor-like
                "MKI67",
                "TOP2A",
                "NES",
                # Middle / immature
                "POU3F3",
                "SOX2",
                "LHX2",
                # Late / mature
                "AQP4",
                "GJA1",
                "GRM3",
            ],
            height=300,
            width=800,
            output_formats=["pdf", "svg"],
        )

        log_message(
            f"Visualization files have been generated in {output_path}",
            message_type="success",
        )

    except Exception as e:
        log_message(f"Error creating Sankey diagram: {str(e)}", message_type="error")
        import traceback

        traceback.print_exc()
        sys.exit(1)
