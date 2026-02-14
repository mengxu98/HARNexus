import os
import sys
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message

# Get project root directory (two levels up from this script)
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))


def extract_genes_from_networks(networks_dict, species_name):
    gene_records = []

    for celltype, net_df in networks_dict.items():
        if net_df is None or net_df.empty:
            continue

        tfs = set(net_df["regulator"].unique())
        targets = set(net_df["target"].unique())

        for gene in tfs:
            gene_records.append({"gene": gene, "is_TF": True, "CellType": celltype})

        for gene in targets:
            if gene not in tfs:
                gene_records.append(
                    {"gene": gene, "is_TF": False, "CellType": celltype}
                )

    if not gene_records:
        log_message(f"  No genes found for {species_name}", message_type="warning")
        return pd.DataFrame(columns=["gene", "is_TF", "CellType"])

    genes_df = pd.DataFrame(gene_records)

    genes_df = genes_df.sort_values(
        by=["gene", "CellType", "is_TF"], ascending=[True, True, False]
    )

    genes_df = genes_df.drop_duplicates(subset=["gene", "CellType"], keep="first")

    return genes_df


def load_networks_from_csv(res_dir, species, sample):
    csv_dir = os.path.join(res_dir, "csv")
    networks_dict = {}

    if not os.path.exists(csv_dir):
        log_message(f"CSV directory not found: {csv_dir}", message_type="warning")
        return networks_dict

    prefix = f"{species}_{sample}_"
    for filename in os.listdir(csv_dir):
        if filename.startswith(prefix) and filename.endswith(".csv"):
            celltype = filename.replace(prefix, "").replace(".csv", "")
            csv_path = os.path.join(csv_dir, filename)

            try:
                net_df = pd.read_csv(csv_path)
                if (
                    not net_df.empty
                    and "regulator" in net_df.columns
                    and "target" in net_df.columns
                ):
                    networks_dict[celltype] = net_df
                    log_message(f"  Loaded {species} {celltype}: {len(net_df)} edges")
                else:
                    log_message(
                        f"  Skipping {species} {celltype}: invalid format",
                        message_type="warning",
                    )
            except Exception as e:
                log_message(
                    f"  Error loading {species} {celltype}: {e}", message_type="error"
                )

    return networks_dict


def main():
    # Define sample pairs to loop through
    sample_pairs = [
        {"human": "h4", "chimp": "c4"},
        {"human": "h3", "chimp": "c2"},
        {"human": "h1", "chimp": "c1"},
    ]

    # Loop through each sample pair
    for pair in sample_pairs:
        human_sample = pair["human"]
        chimp_sample = pair["chimp"]

        res_dir = os.path.join(PROJECT_ROOT, f"results/species_networks/{human_sample}_{chimp_sample}/")

        log_message(f"Loading network data for {human_sample} and {chimp_sample}...")

        log_message("Loading human networks...")
        human_networks = load_networks_from_csv(res_dir, "human", human_sample)

        log_message("Loading chimp networks...")
        chimp_networks = load_networks_from_csv(res_dir, "chimp", chimp_sample)

        log_message("Extracting genes from networks...")

        log_message("Extracting genes from human networks...")
        human_genes_df = extract_genes_from_networks(human_networks, "human")

        log_message("Extracting genes from chimp networks...")
        chimp_genes_df = extract_genes_from_networks(chimp_networks, "chimp")

        human_genes_df = human_genes_df.copy()
        human_genes_df["is_TF"] = human_genes_df["is_TF"].map(
            {True: "TRUE", False: "FALSE"}
        )
        human_output_file = os.path.join(res_dir, "network_genes_by_celltype_human.csv")
        human_genes_df.to_csv(human_output_file, index=False, quoting=0)
        log_message(
            "Human genes saved to {.file ",
            human_output_file,
            "}: ",
            len(human_genes_df),
            " records",
        )

        chimp_genes_df = chimp_genes_df.copy()
        chimp_genes_df["is_TF"] = chimp_genes_df["is_TF"].map(
            {True: "TRUE", False: "FALSE"}
        )
        chimp_output_file = os.path.join(res_dir, "network_genes_by_celltype_chimp.csv")
        chimp_genes_df.to_csv(chimp_output_file, index=False, quoting=0)
        log_message(
            "Chimp genes saved to {.file ",
            chimp_output_file,
            "}: ",
            len(chimp_genes_df),
            " records",
        )

        log_message(
            f"Gene extraction completed for {human_sample} and {chimp_sample}!",
            message_type="success",
        )


if __name__ == "__main__":
    main()
