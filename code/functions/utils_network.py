"""
Network utility functions for Sankey diagram visualization
"""
import pandas as pd
from pathlib import Path
from typing import List, Optional
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message


def load_network_data(
    network_dir: str = "results/networks/har_csn_atlas",
    regions: Optional[List[str]] = None,
    stages: Optional[List[str]] = None,
    celltypes: Optional[List[str]] = None,
    tfs: Optional[List[str]] = None,
    targets: Optional[List[str]] = None,
    hars: Optional[List[str]] = None,
    har_tf_data: Optional[pd.DataFrame] = None,
    file_path: Optional[str] = None,
):
    """
    Load network data intelligently, automatically selecting data files based on specified conditions.
    
    Parameters:
    -----------
    network_dir : str
        Data directory base path
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
        List of specified HARs (for filtering, requires har_tf_data)
    har_tf_data : pd.DataFrame, optional
        HAR-TF mapping data
    file_path : str, optional
        Direct path to network data CSV file (if provided, will use this instead of smart loading)
    
    Returns:
    --------
    pd.DataFrame : Loaded network data
    """
    # If direct file path is provided, use it
    if file_path and Path(file_path).exists():
        log_message(f"Loading from direct path: {file_path}", message_type="info")
        combined_data = pd.read_csv(file_path)
    else:
        csv_dir = Path(network_dir) / "csv"
        all_data_path = Path(network_dir) / "network_data.csv"
        
        # If specific Region, Stage, CellType combination is specified, try to read corresponding CSV files
        if regions and stages and celltypes:
            data_frames = []
            matched_files = []
            
            for region in regions:
                for stage in stages:
                    for celltype in celltypes:
                        # Build filename: {Region}_{Stage}_{CellType}.csv
                        filename = f"{region}_{stage}_{celltype}.csv"
                        file_path = csv_dir / filename
                        
                        if file_path.exists():
                            try:
                                log_message(f"Loading {filename}...", message_type="running")
                                df = pd.read_csv(file_path)

                                # Handle column name differences: regulator -> TF, target -> Target, weight -> Weight
                                if "regulator" in df.columns:
                                    df = df.rename(
                                        columns={
                                            "regulator": "TF",
                                            "target": "Target",
                                            "weight": "Weight",
                                        }
                                    )

                                # Add metadata columns (extracted from filename)
                                df["Region"] = region
                                df["Stage"] = stage
                                df["CellType"] = celltype

                                data_frames.append(df)
                                matched_files.append(filename)
                            except Exception as e:
                                log_message(
                                    f"Warning: Failed to load {filename}: {str(e)}",
                                    message_type="warning",
                                )

            if data_frames:
                log_message(
                    f"Found {len(matched_files)} matching files, merging...",
                    message_type="info",
                )
                combined_data = pd.concat(data_frames, ignore_index=True)
                log_message(
                    f"Merged data shape: {combined_data.shape}", message_type="info"
                )
            else:
                log_message(
                    "No matching CSV files found, loading from main data file...",
                    message_type="warning",
                )
                if all_data_path.exists():
                    combined_data = pd.read_csv(all_data_path)
                else:
                    raise FileNotFoundError(f"Main data file not found: {all_data_path}")
        else:
            # If no complete conditions specified, read main data file
            log_message(
                f"Loading from main data file: {all_data_path}", message_type="info"
            )
            if not all_data_path.exists():
                raise FileNotFoundError(f"Main data file not found: {all_data_path}")
            combined_data = pd.read_csv(all_data_path)
    
    # Ensure required columns exist
    required_columns = ["TF", "Target", "Weight"]
    optional_columns = ["Region", "Stage", "CellType"]
    
    for col in required_columns:
        if col not in combined_data.columns:
            raise ValueError(f"Required column '{col}' not found in data")
    
    # Add missing optional columns
    for col in optional_columns:
        if col not in combined_data.columns:
            combined_data[col] = "Unknown"
    
    # Data cleaning
    combined_data = combined_data.fillna({
        "Weight": 0,
        "TF": "Unknown_TF",
        "Target": "Unknown_Target",
        "Stage": "Unknown_Stage",
        "Region": "Unknown_Region",
        "CellType": "Unknown_CellType"
    })
    
    # Filter by HAR (if HARs and har_tf_data are provided)
    if hars and har_tf_data is not None:
        # Find TFs related to specified HARs
        relevant_tfs = har_tf_data[har_tf_data["HAR"].isin(hars)]["TF"].unique()
        combined_data = combined_data[combined_data["TF"].isin(relevant_tfs)]
        log_message(
            f"Filtered by HARs: {len(relevant_tfs)} TFs found", message_type="info"
        )

    # Filter by TF
    if tfs:
        combined_data = combined_data[combined_data["TF"].isin(tfs)]
        log_message(f"Filtered by TFs: {len(tfs)} TFs", message_type="info")

    # Filter by Target
    if targets:
        combined_data = combined_data[combined_data["Target"].isin(targets)]
        n_found = combined_data["Target"].nunique()
        if n_found < len(targets):
            log_message(
                f"Filtered by Targets: {len(targets)} specified, {n_found} found in data",
                message_type="info",
            )
        else:
            log_message(f"Filtered by Targets: {len(targets)} targets", message_type="info")
    
    # Filter invalid data
    combined_data = combined_data[
        (combined_data["Weight"] != 0) &
        (combined_data["TF"] != "Unknown_TF") &
        (combined_data["Target"] != "Unknown_Target")
    ].copy()
    
    # Optimize data types
    combined_data["Weight"] = combined_data["Weight"].astype("float32")

    log_message(
        f"Data loaded successfully. Final shape: {combined_data.shape}",
        message_type="success",
    )

    return combined_data


def load_har_tf_data(
    file_path: Optional[str] = None,
):
    """
    Load HAR-TF combinations from CSV file
    
    Parameters:
    -----------
    file_path : str, optional
        Path to HAR-TF CSV file. If None, uses default path:
        results/har_tf/human/har_tf_pairs_scores.csv
    
    Returns:
    --------
    pd.DataFrame or None : HAR-TF mapping data, or None if file not found
    """
    if file_path is None:
        file_path = "results/har_tf/human/har_tf_pairs_scores.csv"
    
    try:
        log_message(f"Loading HAR-TF data from {file_path}...", message_type="info")
        har_tf_data = pd.read_csv(file_path)
        
        # Handle column names: ensure TF and HAR columns exist (column name might be lowercase har)
        if "HAR" not in har_tf_data.columns and "har" in har_tf_data.columns:
            har_tf_data = har_tf_data.rename(columns={"har": "HAR"})
        
        # Keep only TF and HAR columns, and remove duplicates
        if "TF" in har_tf_data.columns and "HAR" in har_tf_data.columns:
            har_tf_data = har_tf_data[["TF", "HAR"]].drop_duplicates()
        else:
            raise ValueError(
                f"Required columns 'TF' and 'HAR' (or 'har') not found in {file_path}"
            )
        
        log_message(
            f"HAR-TF data loaded successfully. Shape: {har_tf_data.shape}",
            message_type="success",
        )
        return har_tf_data
    except FileNotFoundError:
        log_message(
            f"HAR-TF data file not found: {file_path}",
            message_type="warning",
        )
        log_message(
            "HAR-TF data will not be included in the Sankey diagram.",
            message_type="warning",
        )
        return None
    except Exception as e:
        log_message(
            f"Error loading HAR-TF data from {file_path}: {str(e)}",
            message_type="error",
        )
        return None


def _hex_to_rgba(hex_str: str, alpha: float) -> str:
    """Convert #RRGGBB to rgba(r,g,b,alpha)."""
    h = hex_str.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def generate_colors(n):
    """
    Generate n distinct colors for network visualization.
    Base palette from prepare_env.R colors9.
    """
    if n <= 0:
        return []

    # Base palette (from prepare_env.R colors9)
    base_colors = [
        "#8076A3",
        "#ED5736",
        "#0AA344",
        "#2177B8",
        "#D70440",
        "#F9BD10",
        "#B14B28",
        "#006D87",
        "#5E7987",
    ]

    colors = []
    # Repeat base colors with different alpha when n > len(base_colors).
    # Keep alpha in [0.2, 1] (Plotly rejects negative or >1 alpha).
    while len(colors) < n:
        raw = 0.8 - (0.3 * (len(colors) // len(base_colors)))
        alpha = max(0.2, min(1.0, raw))
        color_idx = len(colors) % len(base_colors)
        colors.append(_hex_to_rgba(base_colors[color_idx], alpha))

    return colors

