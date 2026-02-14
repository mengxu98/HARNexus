import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message


def check_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return path


log_message("Loading network data...")
res_dir = check_dir("results/networks/analysis/")

regions = ["Prefrontal cortex", "Cerebral cortex"]

source_file = os.path.join(res_dir, "network_data.csv")
if not os.path.exists(source_file):
    log_message(f"Source file {source_file} not found.", message_type="error")
    sys.exit(1)

all_network_data = pd.read_csv(source_file)
available_regions = all_network_data["Region"].unique()
log_message(f"Available regions: {', '.join(available_regions)}")

for region in regions:
    if region not in available_regions:
        log_message(
            f"Warning: Region '{region}' not found in data, skipping",
            message_type="warning",
        )
        continue

    log_message(f"Processing region: {region}")
    file_path = os.path.join(res_dir, f"network_{region}.csv")

    network_data = all_network_data[all_network_data["Region"] == region].copy()
    network_data.to_csv(file_path, index=False)
    log_message(f"Filtered {len(network_data):,} edges for {region}")


def calculate_similarity(network_data, type_col):
    log_message(f"Calculating similarity matrix for: {{.val {type_col}}}")

    if type_col == "Stage":
        groups = network_data["Stage"].unique()
        group_numbers = []
        for g in groups:
            try:
                num = int(str(g).replace("S", ""))
                group_numbers.append(num)
            except ValueError:
                group_numbers.append(float("inf"))

        sorted_indices = np.argsort(group_numbers)
        groups_ordered = groups[sorted_indices]
    else:
        groups = network_data[type_col].unique()
        groups_ordered = sorted(groups)

    n_groups = len(groups_ordered)
    similarity_matrix = np.zeros((n_groups, n_groups))

    edge_sets = {}
    for g in groups_ordered:
        group_data = network_data[network_data[type_col] == g]
        edge_sets[g] = set(group_data["TF"] + "-" + group_data["Target"])

    for i in range(n_groups):
        for j in range(n_groups):
            set1 = edge_sets[groups_ordered[i]]
            set2 = edge_sets[groups_ordered[j]]

            intersection = len(set1.intersection(set2))
            union_len = len(set1.union(set2))

            similarity_matrix[i, j] = intersection / union_len if union_len > 0 else 0

    df_similarity = pd.DataFrame(
        similarity_matrix, index=groups_ordered, columns=groups_ordered
    )
    return df_similarity


tfs_file = "results/har_tf/tfs.csv"
if os.path.exists(tfs_file):
    tfs_df = pd.read_csv(tfs_file)
    tfs_set = set(tfs_df.iloc[:, 0].astype(str))
    log_message(f"Loaded {len(tfs_set)} TFs from {tfs_file}")
else:
    log_message(
        f"TF file {tfs_file} not found, skipping TF annotation", message_type="warning"
    )
    tfs_set = set()

for region in regions:
    if region not in available_regions:
        continue

    file_path = os.path.join(res_dir, f"network_{region}.csv")

    if not os.path.exists(file_path):
        continue

    log_message(f"Analyzing region: {region}")
    network_data = pd.read_csv(file_path)

    log_message("Calculating similarity matrices...")
    similarity_matrices = {}

    for type_col in ["Stage", "CellType"]:
        similarity_matrices[type_col] = calculate_similarity(network_data, type_col)

    similarity_matrices["Stage"].to_csv(
        os.path.join(res_dir, f"similarity_matrix_Stage_{region}.csv")
    )
    similarity_matrices["CellType"].to_csv(
        os.path.join(res_dir, f"similarity_matrix_CellType_{region}.csv")
    )

    log_message("Extracting cell type-specific genes...")
    celltype_genes_list = []
    cell_types = network_data["CellType"].unique()

    for cell_type in cell_types:
        cell_type_data = network_data[network_data["CellType"] == cell_type]
        target_genes = cell_type_data["Target"].unique()

        for gene in target_genes:
            is_tf = gene in tfs_set
            celltype_genes_list.append(
                {"gene": gene, "info": cell_type, "is_TF": is_tf}
            )

    if celltype_genes_list:
        celltype_genes_df = pd.DataFrame(celltype_genes_list)
        celltype_genes_df["CellType"] = celltype_genes_df["info"]
        output_file = os.path.join(
            res_dir, f"celltype_specific_genes_{region}.csv"
        )
        celltype_genes_df.to_csv(output_file, index=False)
        log_message(f"Saved {len(celltype_genes_df)} gene entries to {output_file}")
        log_message(
            f"  Unique genes: {celltype_genes_df['gene'].nunique()}, "
            f"Cell types: {celltype_genes_df['info'].nunique()}, "
            f"TFs: {celltype_genes_df['is_TF'].sum()}"
        )
    else:
        log_message("No cell type-specific genes found", message_type="warning")

    log_message("Extracting stage-specific genes...")
    stage_genes_list = []
    stages = network_data["Stage"].unique()

    for stage in stages:
        stage_data = network_data[network_data["Stage"] == stage]
        target_genes = stage_data["Target"].unique()

        for gene in target_genes:
            is_tf = gene in tfs_set
            cell_types_in_stage = stage_data[stage_data["Target"] == gene]["CellType"].unique()
            cell_type_str = ", ".join(sorted(cell_types_in_stage))
            
            stage_genes_list.append(
                {
                    "gene": gene,
                    "info": stage,
                    "is_TF": is_tf,
                    "Stage": stage,
                    "CellType": cell_type_str
                }
            )

    if stage_genes_list:
        stage_genes_df = pd.DataFrame(stage_genes_list)
        output_file = os.path.join(
            res_dir, f"stage_specific_genes_{region}.csv"
        )
        stage_genes_df.to_csv(output_file, index=False)
        log_message(f"Saved {len(stage_genes_df)} gene entries to {output_file}")
        log_message(
            f"  Unique genes: {stage_genes_df['gene'].nunique()}, "
            f"Stages: {stage_genes_df['Stage'].nunique()}, "
            f"TFs: {stage_genes_df['is_TF'].sum()}"
        )
    else:
        log_message("No stage-specific genes found", message_type="warning")
