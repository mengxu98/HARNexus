import os
import sys
import pandas as pd
import numpy as np
from typing import Dict, Set

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message, check_dir

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))


def load_network_from_csv(csv_path: str) -> pd.DataFrame:
    if not os.path.exists(csv_path):
        return pd.DataFrame()
    try:
        network_df = pd.read_csv(csv_path)
        return network_df
    except Exception as e:
        log_message(f"Error loading network from {csv_path}: {e}", message_type="error")
        return pd.DataFrame()


def calculate_jaccard_similarity(set1: Set, set2: Set) -> float:
    if len(set1) == 0 and len(set2) == 0:
        return 0.0
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union > 0 else 0.0


def calculate_network_stats(network_df: pd.DataFrame) -> Dict:
    if network_df.empty or len(network_df) == 0:
        return {"n_edges": 0, "n_tfs": 0, "n_targets": 0, "n_genes": 0}

    n_edges = len(network_df)
    n_tfs = network_df["regulator"].nunique()
    n_targets = network_df["target"].nunique()
    n_genes = len(
        set(network_df["regulator"].unique()) | set(network_df["target"].unique())
    )

    return {
        "n_edges": n_edges,
        "n_tfs": n_tfs,
        "n_targets": n_targets,
        "n_genes": n_genes,
    }


def main():
    sample_pairs = [
        {"human": "h4", "chimp": "c4"},
        {"human": "h3", "chimp": "c2"},
        {"human": "h1", "chimp": "c1"},
    ]

    exclude_tfs_in_genes = True

    for pair in sample_pairs:
        human_sample = pair["human"]
        chimp_sample = pair["chimp"]

        log_message(
            f"Network comparison for Human {human_sample} and Chimpanzee {chimp_sample}"
        )

        res_dir = os.path.join(PROJECT_ROOT, f"results/species_networks/{human_sample}_{chimp_sample}/")
        output_dir = check_dir(res_dir)
        csv_dir = os.path.join(res_dir, "csv/")

        log_message("Loading networks from CSV files...")

        network_files = []
        if os.path.exists(csv_dir):
            for file in os.listdir(csv_dir):
                if file.endswith(".csv"):
                    if file.startswith(f"human_{human_sample}_"):
                        celltype = file.replace(f"human_{human_sample}_", "").replace(
                            ".csv", ""
                        )
                        network_files.append(
                            {
                                "Species": "Human",
                                "CellType": celltype,
                                "CSV_File": os.path.join(csv_dir, file),
                            }
                        )
                    elif file.startswith(f"chimp_{chimp_sample}_"):
                        celltype = file.replace(f"chimp_{chimp_sample}_", "").replace(
                            ".csv", ""
                        )
                        network_files.append(
                            {
                                "Species": "Chimpanzee",
                                "CellType": celltype,
                                "CSV_File": os.path.join(csv_dir, file),
                            }
                        )

        human_files = [f for f in network_files if f["Species"] == "Human"]
        chimp_files = [f for f in network_files if f["Species"] == "Chimpanzee"]

        human_celltypes = sorted(list(set([f["CellType"] for f in human_files])))

        if not human_celltypes:
            log_message(
                f"Warning: No human network CSV files found in {csv_dir}. "
                f"Please ensure network construction completed successfully.",
                message_type="warning"
            )
            continue

        log_message(f"Human cell types: {', '.join(human_celltypes)}")

        human_networks = {}
        chimp_networks = {}

        for file_info in human_files:
            ct = file_info["CellType"]
            csv_path = file_info["CSV_File"]
            net_df = load_network_from_csv(csv_path)
            if exclude_tfs_in_genes:
                net_df = net_df[~net_df["target"].isin(net_df["regulator"])]
            if not net_df.empty and len(net_df) > 0:
                human_networks[ct] = net_df
            else:
                human_networks[ct] = pd.DataFrame()

        for file_info in chimp_files:
            ct = file_info["CellType"]
            csv_path = file_info["CSV_File"]
            net_df = load_network_from_csv(csv_path)
            if exclude_tfs_in_genes:
                net_df = net_df[~net_df["target"].isin(net_df["regulator"])]
            if not net_df.empty and len(net_df) > 0:
                chimp_networks[ct] = net_df
            else:
                chimp_networks[ct] = pd.DataFrame()

        human_celltypes = [
            ct
            for ct in human_celltypes
            if ct in human_networks and not human_networks[ct].empty
        ]

        comparison_results = []
        human_specific_targets_dict = {}

        for ct in human_celltypes:
            log_message(f"Processing cell type: {ct}...")

            net_human = human_networks[ct]

            stats_human = calculate_network_stats(net_human)
            n_edges_human = stats_human["n_edges"]
            n_tfs_human = stats_human["n_tfs"]
            n_targets_human = stats_human["n_targets"]
            n_genes_human = stats_human["n_genes"]

            edges_human = set(net_human["regulator"] + "->" + net_human["target"])
            tfs_human = set(net_human["regulator"].unique())
            targets_human = set(net_human["target"].unique())

            net_chimp = chimp_networks.get(ct, pd.DataFrame())
            chimp_exists = not net_chimp.empty and len(net_chimp) > 0

            if chimp_exists:
                stats_chimp = calculate_network_stats(net_chimp)
                n_edges_chimp = stats_chimp["n_edges"]
                n_tfs_chimp = stats_chimp["n_tfs"]
                n_targets_chimp = stats_chimp["n_targets"]
                n_genes_chimp = stats_chimp["n_genes"]

                edges_chimp = set(net_chimp["regulator"] + "->" + net_chimp["target"])
                tfs_chimp = set(net_chimp["regulator"].unique())
                targets_chimp = set(net_chimp["target"].unique())

                jaccard_tfs = calculate_jaccard_similarity(tfs_human, tfs_chimp)
                jaccard_targets = calculate_jaccard_similarity(
                    targets_human, targets_chimp
                )
                jaccard_edges = calculate_jaccard_similarity(edges_human, edges_chimp)

                human_specific_targets = list(targets_human - targets_chimp)

                log_message(
                    f"  Jaccard similarity - TFs: {jaccard_tfs:.3f}, "
                    f"Targets: {jaccard_targets:.3f}, Edges: {jaccard_edges:.3f}"
                )
                log_message(f"  Human-specific targets: {len(human_specific_targets)}")
            else:
                n_edges_chimp = None
                n_tfs_chimp = None
                n_targets_chimp = None
                n_genes_chimp = None
                jaccard_tfs = None
                jaccard_targets = None
                jaccard_edges = None
                human_specific_targets = list(targets_human)

                log_message(
                    "  Chimpanzee network not available, setting Jaccard similarity to NA"
                )
                log_message(
                    f"  All human targets considered as human-specific: {len(human_specific_targets)}"
                )

            comparison_results.append(
                {
                    "CellType": ct,
                    "Human_Edges": n_edges_human,
                    "Human_TFs": n_tfs_human,
                    "Human_Targets": n_targets_human,
                    "Human_Genes": n_genes_human,
                    "Chimp_Edges": n_edges_chimp,
                    "Chimp_TFs": n_tfs_chimp,
                    "Chimp_Targets": n_targets_chimp,
                    "Chimp_Genes": n_genes_chimp,
                    "Jaccard_TFs": jaccard_tfs,
                    "Jaccard_Targets": jaccard_targets,
                    "Jaccard_Edges": jaccard_edges,
                }
            )

            human_specific_targets_dict[ct] = human_specific_targets

        summary_df = pd.DataFrame(comparison_results)

        if summary_df.empty or "CellType" not in summary_df.columns:
            log_message(
                "Warning: No valid networks found for comparison. Skipping similarity calculations.",
                message_type="warning"
            )
            if summary_df.empty:
                summary_df = pd.DataFrame(columns=[
                    "CellType", "Human_Edges", "Human_TFs", "Human_Targets", 
                    "Human_Genes", "Chimp_Edges", "Chimp_TFs", "Chimp_Targets",
                    "Chimp_Genes", "Jaccard_TFs", "Jaccard_Targets", "Jaccard_Edges"
                ])
            summary_df.to_csv(
                os.path.join(output_dir, "jaccard_similarity.csv"), index=False
            )
            pd.DataFrame(columns=["CellType", "TargetGene"]).to_csv(
                os.path.join(output_dir, "human_specific_targets.csv"), index=False
            )
            log_message(
                f"Network comparison completed for {human_sample} and {chimp_sample} (no valid networks found)!"
            )
            continue

        summary_df.to_csv(
            os.path.join(output_dir, "jaccard_similarity.csv"), index=False
        )

        human_specific_targets_list = []
        for ct, targets in human_specific_targets_dict.items():
            for target in targets:
                human_specific_targets_list.append(
                    {"CellType": ct, "TargetGene": target}
                )

        if human_specific_targets_list:
            human_specific_targets_df = pd.DataFrame(human_specific_targets_list)
            human_specific_targets_df.to_csv(
                os.path.join(output_dir, "human_specific_targets.csv"), index=False
            )
        else:
            pd.DataFrame(columns=["CellType", "TargetGene"]).to_csv(
                os.path.join(output_dir, "human_specific_targets.csv"), index=False
            )

        log_message("Calculating similarity matrices between cell types...")

        chimp_celltypes_exist = {}
        for ct in human_celltypes:
            net_chimp = chimp_networks.get(ct, pd.DataFrame())
            chimp_celltypes_exist[ct] = not net_chimp.empty and len(net_chimp) > 0

        human_celltypes_ordered = sorted(human_celltypes)
        n_ct = len(human_celltypes_ordered)

        edge_sets_human = {}
        for ct in human_celltypes_ordered:
            net_human = human_networks[ct]
            if not net_human.empty and len(net_human) > 0:
                edge_sets_human[ct] = set(
                    net_human["regulator"] + "-" + net_human["target"]
                )
            else:
                edge_sets_human[ct] = set()

        edge_sets_chimp = {}
        for ct in human_celltypes_ordered:
            if chimp_celltypes_exist[ct]:
                net_chimp = chimp_networks[ct]
                edge_sets_chimp[ct] = set(
                    net_chimp["regulator"] + "-" + net_chimp["target"]
                )
            else:
                edge_sets_chimp[ct] = set()

        log_message("Calculating similarity matrix for human networks...")
        similarity_matrix_human = np.full((n_ct, n_ct), np.nan)

        for i, ct_i in enumerate(human_celltypes_ordered):
            edges_i = edge_sets_human[ct_i]

            for j, ct_j in enumerate(human_celltypes_ordered):
                edges_j = edge_sets_human[ct_j]

                similarity_matrix_human[i, j] = calculate_jaccard_similarity(
                    edges_i, edges_j
                )

        log_message("Calculating similarity matrix for chimp networks...")
        similarity_matrix_chimp = np.full((n_ct, n_ct), np.nan)

        for i, ct_i in enumerate(human_celltypes_ordered):
            for j, ct_j in enumerate(human_celltypes_ordered):
                if chimp_celltypes_exist[ct_i] and chimp_celltypes_exist[ct_j]:
                    edges_i = edge_sets_chimp[ct_i]
                    edges_j = edge_sets_chimp[ct_j]

                    similarity_matrix_chimp[i, j] = calculate_jaccard_similarity(
                        edges_i, edges_j
                    )

        log_message("Calculating human-chimp similarity for diagonal...")
        similarity_vector_human_chimp = {}

        for ct in human_celltypes_ordered:
            if chimp_celltypes_exist[ct]:
                edges_human = edge_sets_human[ct]
                edges_chimp = edge_sets_chimp[ct]

                similarity_vector_human_chimp[ct] = calculate_jaccard_similarity(
                    edges_human, edges_chimp
                )
            else:
                similarity_vector_human_chimp[ct] = None

        summary_df["Jaccard_Edges_HumanChimp"] = summary_df["CellType"].map(
            similarity_vector_human_chimp
        )

        similarity_matrix_human_df = pd.DataFrame(
            similarity_matrix_human,
            index=human_celltypes_ordered,
            columns=human_celltypes_ordered,
        )
        similarity_matrix_human_df.to_csv(
            os.path.join(output_dir, "similarity_matrix_human.csv")
        )

        similarity_matrix_chimp_df = pd.DataFrame(
            similarity_matrix_chimp,
            index=human_celltypes_ordered,
            columns=human_celltypes_ordered,
        )
        similarity_matrix_chimp_df.to_csv(
            os.path.join(output_dir, "similarity_matrix_chimp.csv")
        )

        similarity_vector_human_chimp_df = pd.DataFrame(
            {
                "CellType": list(similarity_vector_human_chimp.keys()),
                "Jaccard_Edges_HumanChimp": list(
                    similarity_vector_human_chimp.values()
                ),
            }
        )
        similarity_vector_human_chimp_df.to_csv(
            os.path.join(output_dir, "similarity_vector_human_chimp.csv"), index=False
        )

        summary_df.to_csv(
            os.path.join(output_dir, "jaccard_similarity.csv"), index=False
        )

        log_message(
            f"Network comparison completed for {human_sample} and {chimp_sample}!"
        )
        log_message(f"Processed {len(comparison_results)} cell types")
        log_message("Saved similarity matrix (Jaccard similarity based on edges)")
        log_message(f"Results saved to {output_dir}")


if __name__ == "__main__":
    main()
