import pandas as pd
import numpy as np
import plotly.graph_objects as go
from pathlib import Path
import sys
import matplotlib.pyplot as plt
from typing import List, Optional
from tqdm import tqdm
from colorama import init, Fore, Style
from datetime import datetime

# Initialize colorama
init()


def get_colored_time():
    """Get current time with color formatting"""
    current_time = datetime.now().strftime("%H:%M:%S")
    return f"{Fore.CYAN}[{current_time}]{Style.RESET_ALL}"


def load_network_data(file_path: str):
    """Load network data from CSV file with error handling"""
    try:
        print(f"{get_colored_time()} Loading data from {file_path}...")
        data = pd.read_csv(file_path)
        print(f"{get_colored_time()} Processing data...")
        with tqdm(total=6, desc="Cleaning data", ascii=True) as pbar:
            data["Weight"] = data["Weight"].fillna(0)
            pbar.update(1)
            data["TF"] = data["TF"].fillna("Unknown_TF")
            pbar.update(1)
            data["Target"] = data["Target"].fillna("Unknown_Target")
            pbar.update(1)
            data["Stage"] = data["Stage"].fillna("Unknown_Stage")
            pbar.update(1)
            data["Region"] = data["Region"].fillna("Unknown_Region")
            pbar.update(1)
            data["CellType"] = data["CellType"].fillna("Unknown_CellType")
            pbar.update(1)

        data = data.dropna()
        print(f"{get_colored_time()} Data loaded successfully.")
        return data
    except Exception as e:
        print(
            f"{get_colored_time()} {Fore.RED}Error loading data from {file_path}: {str(e)}{Style.RESET_ALL}"
        )
        print(
            f"{get_colored_time()} Please ensure you have run the network processing script first:"
        )
        print(f"{get_colored_time()} Rscript step01-network_processing.R")
        print(f"{get_colored_time()} Or provide your own data in the correct format:")
        print(f"{get_colored_time()} TF, Target, Weight, CellType, Stage, Region")
        print(f"{get_colored_time()} The required columns are:")
        print(f"{get_colored_time()} TF, Target, Weight")
        print(f"{get_colored_time()} The optional columns are:")
        print(f"{get_colored_time()} CellType, Stage, Region")
        sys.exit(1)


def load_har_tf_data(
    file_path: Optional[str] = None,
):
    """Load HAR-TF combinations from CSV file"""
    if file_path is None:
        print(
            f"{get_colored_time()} {Fore.YELLOW}No HAR-TF data file provided.{Style.RESET_ALL}"
        )
        print(
            f"{get_colored_time()} {Fore.YELLOW}HAR-TF data will not be included in the Sankey diagram.{Style.RESET_ALL}"
        )
        return None
    try:
        print(f"{get_colored_time()} Loading HAR-TF data from {file_path}...")
        har_tf_data = pd.read_csv(file_path)
        print(f"{get_colored_time()} HAR-TF data loaded successfully.")
        return har_tf_data
    except Exception as e:
        print(
            f"{get_colored_time()} {Fore.RED}Error loading HAR-TF data from {file_path}: {str(e)}{Style.RESET_ALL}"
        )
        sys.exit(1)


def create_network_sankey(
    data_path: str,
    output_path: str,
    fig_name: Optional[str] = None,
    har_tf_path: Optional[str] = None,
    top_n_tfs: int = 5,
    top_n_targets: int = 10,
    specific_tfs: Optional[List[str]] = None,
    specific_targets: Optional[List[str]] = None,
    specific_stages: Optional[List[str]] = None,
    specific_regions: Optional[List[str]] = None,
    specific_celltypes: Optional[List[str]] = None,
    locked_levels: Optional[List[str]] = None,
    font_size: int = 10,
    font_family: str = "Arial",
    font_color: str = "black",
    height: int = 600,
    width: int = 900,
):
    """
    Create a dynamic-level Sankey diagram for gene regulatory networks.
    Now includes HAR (Human Accelerated Region) as the first level if har_tf_path is provided.

    Parameters:
    -------------------------------------------------------------------
    data_path : str
        Path to the CSV file containing network data
    output_path : str
        Path to the output directory for the Sankey diagram
    fig_name : str, optional
        Name of the output figure file (default: None)
    har_tf_path : str, optional
        Path to the CSV file containing HAR-TF combinations. If None, HAR level will be omitted.
    specific_celltypes : List[str], optional
        List of specific cell types to visualize.
    top_n_tfs : int
        Number of top TFs to include if specific_tfs is not provided (default: 5)
    top_n_targets : int
        Number of top targets to include if specific_targets is not provided (default: 10)
    specific_tfs : List[str], optional
        List of specific TFs to visualize. If provided, top_n_tfs is ignored.
    specific_targets : List[str], optional
        List of specific targets to visualize. If provided, top_n_targets is ignored.
    specific_stages : List[str], optional
        List of specific developmental stages to visualize.
    specific_regions : List[str], optional
        List of specific brain regions to visualize.
    locked_levels : List[str], optional
        Levels to lock for visualization. Can include "specific_celltypes", "region", and/or "stage".
        Number of displayed levels will be adjusted based on number of locked levels.
    font_size : int
        Font size for the plot (default: 10)
    font_family : str
        Font family for the plot (default: "Arial")
    font_color : str
        Color of the text (default: "black")
    height : int
        Height of the plot in pixels (default: 600)
    width : int
        Width of the plot in pixels (default: 900)
    -------------------------------------------------------------------

    Returns:
    plotly.graph_objects.Figure or None
        The Sankey diagram figure object if successful, None otherwise.
    """
    # Load the network data and HAR-TF combinations
    data = load_network_data(data_path)
    har_tf_data = load_har_tf_data(har_tf_path)

    print(f"{get_colored_time()} Creating Sankey diagram...")
    with tqdm(total=7, desc="Generating visualization", ascii=True) as pbar:
        # Initialize flow_data with basic filtering
        flow_data = data[
            (data["Weight"] > 0)
            & (data["TF"] != "Unknown_TF")
            & (data["Target"] != "Unknown_Target")
        ].copy()
        pbar.update(1)

        # Initialize locked_levels if None
        if locked_levels is None:
            locked_levels = []

        # Apply filters based on specified parameters
        if specific_celltypes:
            if isinstance(specific_celltypes, str):
                flow_data = flow_data[flow_data["CellType"] == specific_celltypes]
            else:  # specific_celltypes is a list
                flow_data = flow_data[flow_data["CellType"].isin(specific_celltypes)]

        if specific_regions:
            flow_data = flow_data[flow_data["Region"].isin(specific_regions)]

        if specific_stages:
            flow_data = flow_data[flow_data["Stage"].isin(specific_stages)]

        pbar.update(1)

        if len(flow_data) > 0:
            # Priority-based selection for TFs and targets
            def get_weighted_data(data, group_by_columns):
                return data.groupby(group_by_columns)["Weight"].sum().reset_index()

            # Define priority order for grouping
            priority_columns = []
            if specific_celltypes:
                priority_columns.append("CellType")
            if specific_stages:
                priority_columns.append("Stage")
            if specific_regions:
                priority_columns.append("Region")

            # Get TFs and their associated HARs
            if specific_tfs is not None:
                available_tfs = set(flow_data["TF"].unique())
                top_tfs = [tf for tf in specific_tfs if tf in available_tfs]
                if not top_tfs:
                    print("None of the specified TFs found in the data.")
                    return None
                if len(top_tfs) < len(specific_tfs):
                    missing_tfs = set(specific_tfs) - set(top_tfs)
                    print(f"Warning: Some specified TFs not found: {missing_tfs}")
            else:
                if specific_targets:
                    # First, get all cell types in the data
                    cell_types = flow_data["CellType"].unique()

                    # For each cell type, get TFs connected to the specific targets
                    cell_type_tfs = {}
                    for ct in cell_types:
                        ct_data = flow_data[flow_data["CellType"] == ct]
                        ct_target_data = ct_data[
                            ct_data["Target"].isin(specific_targets)
                        ]

                        # Get weighted TFs for this cell type
                        if priority_columns:
                            # Remove CellType from priority columns for this calculation
                            ct_priority_cols = [
                                col for col in priority_columns if col != "CellType"
                            ]
                            if ct_priority_cols:
                                group_cols = ct_priority_cols + ["TF"]
                                weighted_data = get_weighted_data(
                                    ct_target_data, group_cols
                                )
                            else:
                                weighted_data = (
                                    ct_target_data.groupby("TF")["Weight"]
                                    .sum()
                                    .reset_index()
                                )
                        else:
                            weighted_data = (
                                ct_target_data.groupby("TF")["Weight"]
                                .sum()
                                .reset_index()
                            )

                        # Get top TFs for this cell type
                        cell_type_tfs[ct] = set(
                            weighted_data.nlargest(top_n_tfs, "Weight")["TF"].tolist()
                        )

                    # Find TFs common across all cell types
                    common_tfs = (
                        set.intersection(*cell_type_tfs.values())
                        if cell_type_tfs
                        else set()
                    )

                    if common_tfs:
                        # If we have common TFs, use them
                        # Sort by total weight across all conditions
                        common_tf_weights = (
                            flow_data[flow_data["TF"].isin(common_tfs)]
                            .groupby("TF")["Weight"]
                            .sum()
                            .sort_values(ascending=False)
                        )
                        top_tfs = common_tf_weights.index.tolist()[:top_n_tfs]
                    else:
                        # If no common TFs, fall back to overall top TFs for the specific targets
                        print(
                            "No common TFs found across all cell types for the specified targets."
                        )
                        target_data = flow_data[
                            flow_data["Target"].isin(specific_targets)
                        ]
                        if priority_columns:
                            group_cols = priority_columns + ["TF"]
                            weighted_data = get_weighted_data(target_data, group_cols)
                            top_tfs = (
                                weighted_data.groupby("TF")["Weight"]
                                .sum()
                                .sort_values(ascending=False)
                                .head(top_n_tfs)
                                .index.tolist()
                            )
                        else:
                            top_tfs = (
                                target_data.groupby("TF")["Weight"]
                                .sum()
                                .sort_values(ascending=False)
                                .head(top_n_tfs)
                                .index.tolist()
                            )
                else:
                    # If no specific targets, use original logic
                    if priority_columns:
                        group_cols = priority_columns + ["TF"]
                        weighted_data = get_weighted_data(flow_data, group_cols)
                        top_tfs = (
                            weighted_data.groupby("TF")["Weight"]
                            .sum()
                            .sort_values(ascending=False)
                            .head(top_n_tfs)
                            .index.tolist()
                        )
                    else:
                        top_tfs = (
                            flow_data.groupby("TF")["Weight"]
                            .sum()
                            .sort_values(ascending=False)
                            .head(top_n_tfs)
                            .index.tolist()
                        )

            # Get HARs associated with selected TFs and group them
            har_nodes = []
            har_to_tf_map = {}
            har_tf_details = []  # Store complete mapping for file output

            if har_tf_data is not None:
                har_tf_groups = (
                    har_tf_data[har_tf_data["TF"].isin(top_tfs)]
                    .groupby("TF")["HAR"]
                    .agg(lambda x: sorted(x))  # Keep as list instead of joining
                    .reset_index()
                )

                # Create combined HAR nodes (one per TF) with simplified display
                for _, row in har_tf_groups.iterrows():
                    tf = row["TF"]
                    har_list = row["HAR"]
                    # Create simplified display name
                    first_har = har_list[0]
                    har_count = len(har_list)
                    har_display = f"HARs: {first_har}...({har_count})"
                    # Store complete list for file output
                    har_full_list = ", ".join(har_list)
                    har_tf_details.append({"TF": tf, "HARs": har_full_list})
                    # Add to visualization data
                    har_nodes.append(har_display)
                    har_to_tf_map[har_display] = tf

                # Save HAR-TF details to file
                har_tf_details_df = pd.DataFrame(har_tf_details)

                # Generate filename based on parameters
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

                # Save HAR-TF mappings
                har_tf_details_df.to_csv(
                    f"{output_path}{fig_name}_har_tf_mappings.csv", index=False
                )

            # Get targets based on input parameters or total weight
            flow_data_filtered = flow_data[flow_data["TF"].isin(top_tfs)]
            if specific_targets is not None:
                available_targets = set(flow_data_filtered["Target"].unique())
                top_targets = [
                    target for target in specific_targets if target in available_targets
                ]
                if not top_targets:
                    print(
                        "No specified target genes found in the data for selected TFs."
                    )
                    return None
                if len(top_targets) < len(specific_targets):
                    missing_targets = set(specific_targets) - set(top_targets)
                    print(
                        f"Warning: Some specified targets not found: {missing_targets}"
                    )
            else:
                if priority_columns:
                    group_cols = priority_columns + ["Target"]
                    weighted_data = get_weighted_data(flow_data_filtered, group_cols)
                    top_targets = (
                        weighted_data.groupby("Target")["Weight"]
                        .sum()
                        .sort_values(ascending=False)
                        .head(top_n_targets)
                        .index.tolist()
                    )
                else:
                    top_targets = (
                        flow_data_filtered.groupby("Target")["Weight"]
                        .sum()
                        .sort_values(ascending=False)
                        .head(top_n_targets)
                        .index.tolist()
                    )
            pbar.update(1)

            # Filter data for visualization
            flow_data = flow_data[
                flow_data["TF"].isin(top_tfs) & flow_data["Target"].isin(top_targets)
            ].copy()

            if len(flow_data) == 0:
                print("No valid connections found for the specified parameters.")
                return None

            # Determine which levels to show
            # TFs are always shown
            levels_to_show = ["TF"]

            # Add HAR level if har_tf_data is available
            if har_tf_data is not None:
                levels_to_show.insert(0, "HAR")

            # If no specific parameters are provided, show all levels
            if not any([specific_celltypes, specific_stages, specific_regions]):
                levels_to_show.extend(["Stage", "Region", "CellType"])
            else:
                # Add intermediate levels based on what's specified and not locked
                if not specific_stages and "stage" not in locked_levels:
                    levels_to_show.append("Stage")
                elif specific_stages and "stage" not in locked_levels:
                    levels_to_show.append("Stage")

                if not specific_regions and "region" not in locked_levels:
                    levels_to_show.append("Region")
                elif specific_regions and "region" not in locked_levels:
                    levels_to_show.append("Region")

                if not specific_celltypes and "specific_celltypes" not in locked_levels:
                    levels_to_show.append("CellType")
                elif specific_celltypes and "specific_celltypes" not in locked_levels:
                    levels_to_show.append("CellType")

            levels_to_show.append("Target")

            # Update level nodes with grouped HARs
            level_nodes = {
                "HAR": har_nodes if har_tf_data is not None else [],
                "TF": top_tfs,
                "Stage": sorted(flow_data["Stage"].unique()),
                "Region": sorted(flow_data["Region"].unique()),
                "CellType": sorted(flow_data["CellType"].unique()),
                "Target": top_targets,
            }

            # Build nodes list in order
            nodes = []
            node_to_idx = {}
            current_idx = 0

            for level in levels_to_show:
                for node in level_nodes[level]:
                    nodes.append(node)
                    node_to_idx[node] = current_idx
                    current_idx += 1

            # Prepare source, target, and value lists for Sankey diagram
            sources = []
            targets = []
            values = []
            link_colors = []
            all_connections = []

            # Create connections between adjacent levels
            for i in range(len(levels_to_show) - 1):
                current_level = levels_to_show[i]
                next_level = levels_to_show[i + 1]

                if current_level == "HAR" and next_level == "TF":
                    # Handle HAR groups to TF connections
                    for har_group in level_nodes["HAR"]:
                        tf = har_to_tf_map[har_group]
                        if har_group in node_to_idx and tf in node_to_idx:
                            sources.append(node_to_idx[har_group])
                            targets.append(node_to_idx[tf])
                            # Use the weight from the network data for this TF
                            tf_weight = flow_data[flow_data["TF"] == tf]["Weight"].sum()
                            values.append(float(tf_weight))
                            all_connections.append(tf)
                else:
                    # Handle other level connections
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

                            # Track the TF for each connection by following the path back
                            if current_level == "TF":
                                all_connections.append(row[current_level])
                            else:
                                # Get all TFs that connect to this target through the current path
                                mask = flow_data[next_level] == row[next_level]
                                if current_level != "HAR":
                                    mask &= (
                                        flow_data[current_level] == row[current_level]
                                    )
                                path_tfs = flow_data[mask]["TF"].unique()
                                # Use the TF with the highest weight for this path
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

            # Create color scheme for nodes
            n_tfs = len(top_tfs)
            colors = _generate_colors(n_tfs)
            tf_color_map = {tf: colors[i % len(colors)] for i, tf in enumerate(top_tfs)}

            # Update link colors to match TF colors, including HAR groups
            link_colors = [tf_color_map[tf] for tf in all_connections]

            # Generate node colors with HAR nodes matching their TF colors
            node_colors = []
            for level in levels_to_show:
                if level == "TF":
                    node_colors.extend(colors[: len(level_nodes[level])])
                elif level == "HAR":
                    # Color each HAR group according to its corresponding TF
                    for har_group in level_nodes["HAR"]:
                        tf = har_to_tf_map[har_group]
                        node_colors.append(tf_color_map[tf])
                else:
                    # For intermediate nodes, use the color of the dominant TF
                    for node in level_nodes[level]:
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
                            node_colors.append("rgba(200, 200, 200, 0.8)")

            fig = go.Figure(
                data=[
                    go.Sankey(
                        node=dict(
                            pad=15,
                            thickness=20,
                            line=dict(color="black", width=0.5),
                            label=nodes,
                            color=node_colors,
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

            # Update layout with title
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

            # Update the layout with font settings
            fig.update_layout(
                title=dict(
                    text=title_text,
                    x=0.1,  # Left-aligned
                    y=0.9,
                    xanchor="left",
                    yanchor="top",
                    font=dict(
                        size=font_size,
                        color=font_color,
                        family=font_family,
                    ),
                ),
                font=dict(
                    size=font_size,
                    color=font_color,
                    family=font_family,
                ),
                height=height,
                width=width,
            )

            # Update node labels with custom font settings
            fig.update_traces(
                textfont=dict(
                    size=font_size,
                    color=font_color,
                    family=font_family,
                ),
                selector=dict(type="sankey"),
            )

            # Create output directory if it doesn't exist
            Path(output_path).mkdir(parents=True, exist_ok=True)

            # Save the plot
            print(f"{get_colored_time()} Saving outputs...")
            fig.write_html(f"{output_path}{fig_name}.html")
            try:
                fig.write_image(f"{output_path}{fig_name}.pdf")
            except Exception as e:
                print(
                    f"{get_colored_time()} {Fore.YELLOW}Warning: Could not save PDF version: {str(e)}{Style.RESET_ALL}"
                )
                print(f"{get_colored_time()} HTML version is still available.")
            pbar.update(1)

            print(
                f"{get_colored_time()} {Fore.GREEN}Visualization completed successfully.{Style.RESET_ALL}"
            )
            return fig
        else:
            print(
                f"{get_colored_time()} {Fore.RED}No data found for the specified parameters.{Style.RESET_ALL}"
            )
            return None


def _generate_colors(n):
    """Generate n distinct colors"""
    if n <= 0:
        return []

    # Use a combination of qualitative colormaps for better distinction
    colors = []
    base_colors = [
        "rgb(25, 50, 200)",  # blue
        "rgb(205, 50, 65)",  # red
        "rgb(255, 127, 14)",  # orange
        "rgb(44, 160, 44)",  # green
        "rgb(25, 190, 200)",  # cyan
        "rgb(148, 103, 189)",  # purple
        "rgb(227, 119, 194)",  # pink
        "rgb(188, 189, 34)",  # olive
    ]

    # Repeat the base colors with different alpha values if needed
    while len(colors) < n:
        alpha = 0.8 - (0.3 * (len(colors) // len(base_colors)))
        color_idx = len(colors) % len(base_colors)
        base_color = (
            base_colors[color_idx].replace("rgb", "rgba").replace(")", f", {alpha})")
        )
        colors.append(base_color)

    return colors


if __name__ == "__main__":
    try:
        Path("data/networks/csv").mkdir(parents=True, exist_ok=True)
        Path("results/networks_sankey").mkdir(parents=True, exist_ok=True)
        data_path = "data/networks/csv/network_data.csv"
        har_tf_path = "data/networks/csv/tf_har_combinations.csv"
        output_path = "results/networks_sankey/"

        # Example 1: Five-level diagram
        create_network_sankey(
            data_path=data_path,
            har_tf_path=har_tf_path,
            output_path=output_path,
            fig_name="all_levels_testing",
            top_n_tfs=10,
            top_n_targets=30,
            height=600,
            width=900,
        )

        create_network_sankey(
            data_path=data_path,
            output_path=output_path,
            fig_name="all_levels_testing_no_HAR",
            top_n_tfs=10,
            top_n_targets=30,
            height=600,
            width=900,
        )

        create_network_sankey(
            data_path=data_path,
            har_tf_path=har_tf_path,
            output_path=output_path,
            fig_name="all_levels_testing_font",
            top_n_tfs=10,
            top_n_targets=30,
            height=600,
            width=900,
            font_size=12,
            font_color="blue",
            font_family="Arial",
        )

        # Example 2: Four-level diagram (lock PFC region)
        create_network_sankey(
            data_path=data_path,
            har_tf_path=har_tf_path,
            output_path=output_path,
            fig_name="fig.6g_6",
            specific_regions=["PFC"],
            specific_stages=["S8", "S9"],
            top_n_tfs=10,
            specific_targets=["EDA2R", "ELN", "CDCP1", "GFAP", "IL17D", "CXCL14"],
            height=450,
            width=1000,
        )

        create_network_sankey(
            data_path=data_path,
            har_tf_path=har_tf_path,
            output_path=output_path,
            fig_name="fig.6g_4",
            specific_regions=["PFC"],
            specific_stages=["S8", "S9"],
            top_n_tfs=10,
            specific_targets=["IGDCC4", "CDH3", "SFRP4", "CDON"],
            height=450,
            width=1000,
        )

        create_network_sankey(
            data_path=data_path,
            har_tf_path=har_tf_path,
            output_path=output_path,
            fig_name="fig.6g_10",
            specific_regions=["PFC"],
            specific_stages=["S8", "S9"],
            top_n_tfs=10,
            specific_targets=[
                "EDA2R",
                "ELN",
                "CDCP1",
                "GFAP",
                "IL17D",
                "CXCL14",
                "IGDCC4",
                "CDH3",
                "SFRP4",
                "CDON",
            ],
            height=300,
            width=900,
        )

        create_network_sankey(
            data_path=data_path,
            har_tf_path=har_tf_path,
            output_path=output_path,
            fig_name="fig.6g_EDA2R",
            specific_regions=["PFC"],
            specific_stages=["S8", "S9"],
            specific_celltypes=["Micro", "Perc"],
            top_n_tfs=10,
            specific_targets=["EDA2R"],
            height=300,
            width=800,
        )

        # Example 3: Four-level diagram (lock Astro cell type)
        create_network_sankey(
            data_path=data_path,
            har_tf_path=har_tf_path,
            output_path=output_path,
            specific_celltypes=["Astro"],
            locked_levels=["specific_celltypes"],
            top_n_tfs=5,
            top_n_targets=10,
        )

        # Example 4: Three-level diagram (lock cell type and stage)
        create_network_sankey(
            data_path=data_path,
            har_tf_path=har_tf_path,
            output_path=output_path,
            specific_celltypes=["ExN"],
            specific_stages=["S1", "S2"],
            locked_levels=["specific_celltypes", "stage"],
            top_n_tfs=5,
        )

        create_network_sankey(
            data_path=data_path,
            output_path=output_path,
            fig_name="no_HAR",
            specific_celltypes=["ExN"],
            specific_stages=["S1", "S2"],
            locked_levels=["specific_celltypes", "stage"],
            top_n_tfs=5,
        )

        print(
            f"{get_colored_time()} Visualization files have been generated in {output_path}"
        )

    except Exception as e:
        print(
            f"{get_colored_time()} {Fore.RED}\nError creating Sankey diagram: {str(e)}{Style.RESET_ALL}"
        )
        import traceback

        traceback.print_exc()
        sys.exit(1)
