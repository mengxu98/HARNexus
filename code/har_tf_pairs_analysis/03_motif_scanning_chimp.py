#!/usr/bin/env python3

import subprocess
import os
import sys
import pandas as pd
import re

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message

# memelite: https://github.com/jmschrei/memesuite-lite
from memelite import fimo


def create_har_mapping():
    """
    Create a mapping from HAR sequence names to HAR names
    """
    har_df = pd.read_csv("results/har_tf/chimp/har_coords.csv")

    har_mapping = {}
    for _, row in har_df.iterrows():
        if "HARsv2_" in str(row["har"]):
            har_name = row["har"]
            har_mapping[har_name] = har_name
        else:
            coord1 = f"{row['chrom']}:{row['start']}-{row['end']}"
            coord2 = f"{row['chrom'].replace('chr', '')}:{row['start']}-{row['end']}"
            coord3 = f"{row['chrom']}:{row['start']},{row['end']}"
            coord4 = f"{row['chrom'].replace('chr', '')}:{row['start']},{row['end']}"

            har_mapping[coord1] = row["har"]
            har_mapping[coord2] = row["har"]
            har_mapping[coord3] = row["har"]
            har_mapping[coord4] = row["har"]

            har_name = row["har"]
            if har_name.startswith("ZOOHAR."):
                har_mapping[har_name] = har_name

    log_message(
        f"Created HAR mapping with {len(har_mapping)} entries", message_type="info"
    )
    return har_mapping


def load_annotation_mapping():
    """
    Load the annotation file and create a mapping from motif IDs to database names
    """
    annotation_file = "data/motifs/fourdatabase_human_mouse.annotation"

    if not os.path.exists(annotation_file):
        log_message(
            f"Annotation file does not exist: {annotation_file}", message_type="warning"
        )
        return {}

    try:
        df = pd.read_csv(annotation_file, sep="\t", header=None)
        df.columns = ["ensg_id", "tf_name", "species", "motif_id", "logo", "database"]

        motif_to_db = {}
        for _, row in df.iterrows():
            motif_id = row["motif_id"]
            database = row["database"]
            motif_to_db[motif_id] = database

        log_message(
            f"Loaded {len(motif_to_db)} motif to database mappings", message_type="info"
        )
        return motif_to_db

    except Exception as e:
        log_message(f"Failed to load annotation file: {str(e)}", message_type="error")
        return {}


def split_fourdatabase_meme():
    """
    Split the fourdatabase_all.meme file into separate meme files for each database
    """
    input_file = "data/motifs/fourdatabase_all.meme"
    output_dir = "data/motifs/fourdatabase_split"

    os.makedirs(output_dir, exist_ok=True)

    motif_to_db = load_annotation_mapping()

    db_files = {}
    db_counts = {}

    with open(input_file, "r") as f:
        lines = f.readlines()

    header_lines = []
    current_motif = None
    current_database = None
    current_motif_lines = []

    for i, line in enumerate(lines):
        line = line.strip()

        if line.startswith("MOTIF"):
            if current_motif and current_database and current_motif_lines:
                if current_database not in db_files:
                    db_file = os.path.join(
                        output_dir, f"{current_database.lower()}.meme"
                    )
                    db_files[current_database] = open(db_file, "w")
                    db_counts[current_database] = 0

                    for header_line in header_lines:
                        db_files[current_database].write(header_line)

                for motif_line in current_motif_lines:
                    db_files[current_database].write(motif_line)
                db_counts[current_database] += 1

            parts = line.split()
            if len(parts) >= 2:
                motif_id = parts[1]
                current_motif = motif_id

                if motif_id.startswith("M0") and len(motif_id) >= 6:
                    current_database = "CISBP"
                elif motif_id.startswith("V_"):
                    current_database = "TRANSFAC"
                elif "H11MO" in motif_id:
                    current_database = "HOCOMOCO"
                elif motif_id.startswith("MA") and "." in motif_id:
                    current_database = "JASPAR"
                else:
                    current_database = motif_to_db.get(motif_id, "UNKNOWN")

                current_motif_lines = [line + "\n"]
            else:
                current_motif = None
                current_database = None
                current_motif_lines = []

        elif not current_motif and not line.startswith("MOTIF"):
            header_lines.append(line + "\n")
        elif current_motif:
            current_motif_lines.append(line + "\n")

    if current_motif and current_database and current_motif_lines:
        if current_database not in db_files:
            db_file = os.path.join(output_dir, f"{current_database.lower()}.meme")
            db_files[current_database] = open(db_file, "w")
            db_counts[current_database] = 0

            for header_line in header_lines:
                db_files[current_database].write(header_line)

        for motif_line in current_motif_lines:
            db_files[current_database].write(motif_line)
        db_counts[current_database] += 1

    for db, file_handle in db_files.items():
        file_handle.close()
        log_message(
            f"Created {db} database file, containing {db_counts[db]} motifs",
            message_type="success",
        )

    return list(db_files.keys())


def run_fimo_for_database(database_name, har_seqs_file):
    """
    Run FIMO scan for a specific database
    """
    if database_name == "JASPAR":
        # meme_file = "data/motifs/fourdatabase_split/jaspar.meme"
        meme_file = "data/motifs/JASPAR2024_CORE_vertebrates.meme"
    else:
        meme_file = f"data/motifs/fourdatabase_split/{database_name.lower()}.meme"

    if not os.path.exists(meme_file):
        log_message(
            f"Database file does not exist: {meme_file}", message_type="warning"
        )
        return []

    log_message(f"Processing {database_name} database...", message_type="running")
    log_message(f"Using meme file: {meme_file}", message_type="info")

    try:
        hits = fimo(meme_file, har_seqs_file)
        total_hits = sum(len(h) for h in hits)
        log_message(f"Found {total_hits} {database_name} hits", message_type="info")
        return hits
    except Exception as e:
        log_message(
            f"FIMO scan failed ({database_name}): {str(e)}", message_type="error"
        )
        return []


def combine_hits_from_databases(database_hits_list):
    """
    Combine hits results from multiple databases
    """
    all_hits = []

    for database_name, hits in database_hits_list:
        for hit_df in hits:
            if len(hit_df) > 0:
                hit_df = hit_df.copy()
                hit_df["database"] = database_name
                all_hits.append(hit_df)

    if all_hits:
        combined_hits = pd.concat(all_hits, ignore_index=True)
        return combined_hits
    else:
        return pd.DataFrame()


def process_database_hits(database_name, hits):
    """
    Process hits results from a single database
    """
    if not hits:
        return None

    result_dir = "results/har_tf/chimp"

    hits_df = combine_hits_from_databases([(database_name, hits)])

    if len(hits_df) == 0:
        return None

    har_mapping = create_har_mapping()

    hits_df["har_name"] = hits_df["sequence_name"].apply(
        lambda x: map_sequence_to_har(x, har_mapping)
    )

    mapped_df = hits_df.dropna(subset=["har_name"])

    if len(mapped_df) == 0:
        log_message(
            f"No {database_name} hits can be mapped to HAR name", message_type="warning"
        )
        return None

    # First, get motif-to-TF mapping
    motif_to_tf = get_motif_to_tf_mapping()

    # Clean motif names
    clean_motif_names = mapped_df["motif_name"].copy()
    mask_duplicate = clean_motif_names.str.contains(" ", regex=False)
    if mask_duplicate.any():
        duplicate_parts = clean_motif_names[mask_duplicate].str.split(" ", expand=True)
        if len(duplicate_parts.columns) >= 2:
            mask_same_parts = (
                duplicate_parts[0] == duplicate_parts[1]
            ) & mask_duplicate
            clean_motif_names[mask_same_parts] = duplicate_parts.loc[mask_same_parts, 0]

    # Update the motif_name column with cleaned names
    mapped_df["motif_name"] = clean_motif_names

    # Map motif names to TF names
    if database_name == "JASPAR":
        jaspar_motif_ids = clean_motif_names.str.split(" ", expand=True)[0]
        mapped_df["tf_name"] = jaspar_motif_ids.map(motif_to_tf)

        unmapped_mask = mapped_df["tf_name"].isna()
        if unmapped_mask.any():
            tf_fallback = (
                clean_motif_names.str.split(" ", n=1, expand=True)
                .iloc[:, 1]
                .fillna("")
                .str.replace(r"[^A-Za-z0-9_+-]", "", regex=True)
                .str.upper()
            )
            tf_fallback = tf_fallback.replace("", pd.NA)
            mapped_df.loc[unmapped_mask, "tf_name"] = tf_fallback[unmapped_mask]
    else:
        mapped_df["tf_name"] = clean_motif_names.map(motif_to_tf)

    if database_name == "TRANSFAC":
        unmapped_mask = mapped_df["tf_name"].isna()
        if unmapped_mask.any():
            transfac_motif_ids = clean_motif_names[unmapped_mask].str.split(
                " ", expand=True
            )[0]
            mapped_df.loc[unmapped_mask, "tf_name"] = transfac_motif_ids.map(
                motif_to_tf
            )

    # Remove rows where tf_name is None
    mapped_df = mapped_df.dropna(subset=["tf_name"])

    if len(mapped_df) == 0:
        log_message(
            f"No {database_name} hits can be mapped to TF names", message_type="warning"
        )
        return None

    # Now aggregate by motif first, then by TF
    # Step 1: Aggregate by motif (keep all motif information)
    motif_summary = (
        mapped_df.sort_values("score", ascending=False)
        .groupby(["har_name", "motif_name"])
        .agg({
            "score": "max",
            "p-value": "min",
            "start": "first",
            "end": "first", 
            "strand": "first",
            "tf_name": "first"
        })
        .reset_index()
    )

    # Step 2: Aggregate by TF (combine multiple motifs for same TF)
    summary_df = (
        motif_summary.groupby(["har_name", "tf_name"])
        .agg({
            "score": "max",  # Take the best score among all motifs for this TF
            "p-value": "min",  # Take the best p-value among all motifs for this TF
            "start": "first",  # Keep start from the best motif
            "end": "first",    # Keep end from the best motif
            "strand": "first", # Keep strand from the best motif
            "motif_name": lambda x: ";".join(x)  # Concatenate all motif names
        })
        .reset_index()
    )

    # Add database column
    summary_df["database"] = database_name

    output_file = os.path.join(result_dir, f"{database_name.lower()}_hits.csv")
    summary_df.to_csv(output_file, index=False)

    log_message(
        f"{database_name}: {len(summary_df)} hits saved to {output_file}",
        message_type="success",
    )

    return summary_df


def map_sequence_to_har(sequence_name, har_mapping):
    """
    Map sequence names to HAR names using exact matching
    """
    if sequence_name in har_mapping:
        return har_mapping[sequence_name]

    return None


# Global variable cache for annotation mapping
_annotation_mapping_cache = None
_motif_to_tf_cache = None
_htf_target_mapping_cache = None


def get_annotation_mapping():
    """
    Get annotation mapping, using cache to avoid duplicate loading
    """
    global _annotation_mapping_cache
    if _annotation_mapping_cache is None:
        _annotation_mapping_cache = load_annotation_mapping()
    return _annotation_mapping_cache


def get_motif_to_tf_mapping():
    """
    Get motif to TF name mapping, using cache to avoid duplicate loading
    """
    global _motif_to_tf_cache
    if _motif_to_tf_cache is None:
        annotation_file = "data/motifs/fourdatabase_human_mouse.annotation"
        if os.path.exists(annotation_file):
            try:
                df = pd.read_csv(annotation_file, sep="\t", header=None)
                df.columns = [
                    "ensg_id",
                    "tf_name",
                    "species",
                    "motif_id",
                    "logo",
                    "database",
                ]
                df_copy = df.copy()
                # Only keep human TFs, exclude mouse TFs
                df_copy = df_copy[df_copy["species"] == "Homo sapiens"]
                df_copy["tf_name"] = df_copy["tf_name"].str.upper()
                _motif_to_tf_cache = dict(zip(df_copy["motif_id"], df_copy["tf_name"]))
            except:
                _motif_to_tf_cache = {}
        else:
            _motif_to_tf_cache = {}
    return _motif_to_tf_cache


def get_htf_target_mapping():
    """
    Get hTFtarget motif to TF name mapping, using cache to avoid duplicate loading
    """
    global _htf_target_mapping_cache
    if _htf_target_mapping_cache is None:
        htf_annotation_file = "data/motifs/hTFtarget.annotation"
        if os.path.exists(htf_annotation_file):
            try:
                df = pd.read_csv(htf_annotation_file, sep="\t", header=None)
                df.columns = ["motif_name", "tf_name"]
                _htf_target_mapping_cache = dict(zip(df["motif_name"], df["tf_name"]))
            except:
                _htf_target_mapping_cache = {}
        else:
            _htf_target_mapping_cache = {}
    return _htf_target_mapping_cache


def standardize_tf_name(motif_name, database_name):
    """
    Standardize TF names
    """
    if database_name == "HTFTARGET":
        htf_mapping = get_htf_target_mapping()

        # Only process human TFs, exclude mouse TFs
        if " Homo sapiens" in motif_name:
            clean_motif_name = motif_name.split(" Homo sapiens")[0]
            if clean_motif_name in htf_mapping:
                return htf_mapping[clean_motif_name]
            else:
                return None
        elif " Mus musculus" in motif_name:
            # Skip mouse TFs
            return None
        else:
            # For motifs without species info, try to map but be cautious
            clean_motif_name = motif_name
            if clean_motif_name in htf_mapping:
                return htf_mapping[clean_motif_name]
            else:
                return None
    else:
        motif_to_tf = get_motif_to_tf_mapping()

        if " " in motif_name:
            parts = motif_name.split()
            if len(parts) == 2 and parts[0] == parts[1]:
                clean_motif_name = parts[0]
            else:
                clean_motif_name = motif_name
        else:
            clean_motif_name = motif_name

        if clean_motif_name in motif_to_tf:
            return motif_to_tf[clean_motif_name]

        if database_name == "JASPAR" and " " in clean_motif_name:
            motif_id = clean_motif_name.split()[0]
            if motif_id in motif_to_tf:
                return motif_to_tf[motif_id]

        if database_name == "TRANSFAC" and " " in clean_motif_name:
            motif_id = clean_motif_name.split()[0]
            if motif_id in motif_to_tf:
                return motif_to_tf[motif_id]

        return None


def process_htf_target_database(har_seqs_file):
    """
    Process hTFtarget database TFBS prediction
    """
    meme_file = "data/motifs/hTFtarget_prediction.motifs_matrix.meme"
    annotation_file = "data/motifs/hTFtarget.annotation"

    if not os.path.exists(meme_file):
        log_message(
            f"hTFtarget meme file does not exist: {meme_file}", message_type="warning"
        )
        return None

    if not os.path.exists(annotation_file):
        log_message(
            f"hTFtarget annotation file does not exist: {annotation_file}",
            message_type="warning",
        )
        return None

    log_message("Processing hTFtarget database...", message_type="running")

    try:
        hits = fimo(meme_file, har_seqs_file)
        if not hits:
            log_message("hTFtarget FIMO scan found no hits", message_type="warning")
            return None

        log_message(f"Found {len(hits)} hTFtarget hits", message_type="info")

        hits_df = combine_hits_from_databases([("HTFTARGET", hits)])

        if len(hits_df) == 0:
            return None

        har_mapping = create_har_mapping()

        hits_df["har_name"] = hits_df["sequence_name"].apply(
            lambda x: map_sequence_to_har(x, har_mapping)
        )

        mapped_df = hits_df.dropna(subset=["har_name"])

        if len(mapped_df) == 0:
            log_message(
                "No hTFtarget hits can be mapped to HAR name", message_type="warning"
            )
            return None

        # Get hTFtarget mapping
        htf_mapping = get_htf_target_mapping()

        # Clean motif names and map to TF names
        clean_motif_names = mapped_df["motif_name"].copy()
        clean_motif_names = clean_motif_names.str.replace(
            " Homo sapiens", "", regex=False
        )
        clean_motif_names = clean_motif_names.str.replace(
            " Mus musculus", "", regex=False
        )

        # Map motif names to TF names
        mapped_df["tf_name"] = clean_motif_names.map(htf_mapping)
        
        # Remove rows where tf_name is None
        mapped_df = mapped_df.dropna(subset=["tf_name"])

        if len(mapped_df) == 0:
            log_message(
                "No hTFtarget hits can be mapped to TF names", message_type="warning"
            )
            return None

        # Now aggregate by motif first, then by TF
        # Step 1: Aggregate by motif (keep all motif information)
        motif_summary = (
            mapped_df.sort_values("score", ascending=False)
            .groupby(["har_name", "motif_name"])
            .agg({
                "score": "max",
                "p-value": "min",
                "start": "first",
                "end": "first", 
                "strand": "first",
                "tf_name": "first"
            })
            .reset_index()
        )

        # Step 2: Aggregate by TF (combine multiple motifs for same TF)
        summary_df = (
            motif_summary.groupby(["har_name", "tf_name"])
            .agg({
                "score": "max",  # Take the best score among all motifs for this TF
                "p-value": "min",  # Take the best p-value among all motifs for this TF
                "start": "first",  # Keep start from the best motif
                "end": "first",    # Keep end from the best motif
                "strand": "first", # Keep strand from the best motif
                "motif_name": lambda x: ";".join(x)  # Concatenate all motif names
            })
            .reset_index()
        )

        # Add database column
        summary_df["database"] = "HTFTARGET"

        result_dir = "results/har_tf/chimp"
        output_file = os.path.join(result_dir, "htf_target_hits.csv")
        summary_df.to_csv(output_file, index=False)

        log_message(
            f"hTFtarget: {len(summary_df)} hits saved to {output_file}",
            message_type="success",
        )

        return summary_df

    except Exception as e:
        log_message(
            f"Processing hTFtarget database failed: {str(e)}", message_type="error"
        )
        return None


def main():
    """
    Main function
    """
    log_message("-" * 60, message_type="info")
    log_message(
        "Separate fourdatabase_all.meme and predict TFBS for each database",
        message_type="info",
    )
    log_message("-" * 60, message_type="info")

    result_dir = "results/har_tf/chimp"
    os.makedirs(result_dir, exist_ok=True)

    har_seqs_file = os.path.join(result_dir, "har_seqs_panTro5.fa")
    if not os.path.exists(har_seqs_file):
        log_message(
            f"HAR sequence file does not exist: {har_seqs_file}", message_type="error"
        )
        return

    log_message("Step 1: Separate fourdatabase_all.meme...", message_type="running")
    databases = split_fourdatabase_meme()

    if not databases:
        log_message("No databases found to separate", message_type="error")
        return

    log_message("Step 2: Run FIMO scan...", message_type="running")
    all_database_hits = []
    processed_results = []

    for database in databases:
        if database == "UNKNOWN":
            log_message(f"Skip processing {database} database", message_type="info")
            continue

        hits = run_fimo_for_database(database, har_seqs_file)
        if hits:
            all_database_hits.append((database, hits))
            result_df = process_database_hits(database, hits)
            if result_df is not None:
                processed_results.append(result_df)

    log_message("Step 3: Process hTFtarget database...", message_type="running")
    htf_target_result = process_htf_target_database(har_seqs_file)
    if htf_target_result is not None:
        processed_results.append(htf_target_result)

    if processed_results:
        log_message("Step 4: Combine all results...", message_type="running")
        combined_df = pd.concat(processed_results, ignore_index=True)

        combined_file = os.path.join(result_dir, "combined_hits.csv")
        combined_df.to_csv(combined_file, index=False)

        log_message(
            f"Combined results saved to: {combined_file}", message_type="success"
        )
        log_message(
            f"Total {len(combined_df)} prediction results", message_type="success"
        )

        db_stats = combined_df["database"].value_counts()
        log_message("Database statistics:", message_type="info")
        for db, count in db_stats.items():
            log_message(f"  {db}: {count} predictions")

    log_message("-" * 60, message_type="success")
    log_message(
        "fourdatabse separation and hTFtarget TFBS prediction completed!",
        message_type="success",
    )
    log_message("-" * 60, message_type="success")


subprocess.run(["bash", "code/datasets/download_motif_files.sh"])

if __name__ == "__main__":
    main()
