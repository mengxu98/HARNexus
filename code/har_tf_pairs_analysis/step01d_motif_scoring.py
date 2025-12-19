import os
import csv
import re
import time
import sys
from Bio import SeqIO
import math

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message


def pwm_to_log_odds(pwm, bg, pseudocount=1e-4, use_log2=False):
    def log2(x):
        return math.log(x, 2)

    transformed = []
    for row in pwm:
        probs = [max(v, pseudocount) for v in row]
        s = sum(probs)
        probs = [v / s for v in probs]
        if use_log2:
            log_row = [log2(probs[i] / bg[i]) for i in range(4)]
        else:
            log_row = [probs[i] / bg[i] for i in range(4)]
        transformed.append(log_row)
    return transformed


# Define 4 databases to analyze
databases = {
    "CISBP": "data/motifs/fourdatabase_split/cisbp.meme",
    "HTFTARGET": "data/motifs/hTFtarget_prediction.motifs_matrix.meme",
    "HOCOMOCO": "data/motifs/fourdatabase_split/hocomoco.meme",
    "JASPAR": "data/motifs/JASPAR2024_CORE_vertebrates.meme",
}

human_fasta = "results/har_tf/human/har_seqs_hg38.fa"
chimp_fasta = "results/har_tf/chimp/har_seqs_panTro5.fa"
human_csv = "results/har_tf/human/har_tf_pairs_scores.csv"
chimp_csv = "results/har_tf/chimp/har_tf_pairs_scores.csv"
output_csv = "results/har_tf/motif_score_comparison.csv"

# Database column mapping
db_column_map = {
    "CISBP": "cisbp_motif_name",
    "JASPAR": "jaspar_motif_name",
    "HOCOMOCO": "hocomoco_motif_name",
    "HTFTARGET": "htftarget_motif_name",
}


def parse_meme_pwm(meme_path):
    motifs = []
    names = []
    with open(meme_path) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("MOTIF"):
            parts = line.split()
            name = parts[1]
            i += 1
            # move to the header line of the matrix
            while i < len(lines) and not lines[i].strip().startswith(
                "letter-probability matrix"
            ):
                i += 1
            if i >= len(lines):
                break
            header = lines[i].strip()
            # parse w (number of rows in the matrix)
            m = re.search(r"w\s*=\s*(\d+)", header)
            if not m:
                i += 1
                continue
            w = int(m.group(1))
            i += 1
            matrix = []
            for _ in range(w):
                if i >= len(lines):
                    break
                row_vals = lines[i].strip().split()
                if len(row_vals) < 4:
                    break
                # only take the first 4 columns
                row = list(map(float, row_vals[:4]))
                matrix.append(row)
                i += 1
            # verify the completeness of the matrix
            if len(matrix) == w and all(len(r) == 4 for r in matrix):
                names.append(name)
                motifs.append(matrix)
        else:
            i += 1
    return names, motifs


def is_finite_matrix(mat):
    for row in mat:
        for v in row:
            if not math.isfinite(v):
                return False
    return True


def load_motif_info_from_csv(csv_path):
    """
    Load HAR-TF-motif mapping from CSV file.
    Returns: {(har, TF, db): [motif_list]}
    """
    motif_map = {}
    
    if not os.path.exists(csv_path):
        log_message(f"CSV file not found: {csv_path}", message_type="warning")
        return motif_map
    
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            har = row["har"]
            tf = row["TF"]
            
            # Extract motifs for each database
            for db_name, column_name in db_column_map.items():
                motif_str = row.get(column_name, "").strip()
                if not motif_str or motif_str == "":
                    continue
                
                # Split by semicolon and filter empty strings
                motif_list = [m.strip() for m in motif_str.split(";") if m.strip()]
                if motif_list:
                    key = (har, tf, db_name)
                    motif_map[key] = motif_list
    
    log_message(f"Loaded {len(motif_map)} HAR-TF-database mappings from {csv_path}", message_type="info")
    return motif_map


def match_motif_name(target_name, meme_name, db_name):
    """
    Match motif name from CSV with name in MEME file.
    Supports prefix matching and fuzzy matching.
    """
    # Exact match
    if target_name == meme_name:
        return True
    
    # Prefix matching: target_name should be a prefix of meme_name
    if meme_name.startswith(target_name):
        # Check if next character is space or end of string
        if len(meme_name) == len(target_name) or meme_name[len(target_name)] in [" ", "\t"]:
            return True
    
    # For HTFTARGET: check if target_name is contained in meme_name
    if db_name == "HTFTARGET":
        if target_name in meme_name:
            return True
    
    # For CISBP: handle cases like "M03329" matching "M03329" or "M03329 ..."
    if db_name == "CISBP":
        # Split meme_name by space and check if first part matches
        parts = meme_name.split()
        if parts and parts[0] == target_name:
            return True
    
    # For JASPAR: handle cases like "MA0486.2" matching "MA0486.2 PAX6"
    if db_name == "JASPAR":
        parts = meme_name.split()
        if parts and parts[0] == target_name:
            return True
    
    return False


def process_database(db_name, meme_file, required_motif_names=None):
    """
    Process a single database: load motifs, filter, and convert to PSSM.
    If required_motif_names is provided, only load those motifs.
    Returns: (motif_name_to_pssm_dict, all_motif_names_in_file)
    """
    if not os.path.exists(meme_file):
        log_message(
            f"{meme_file} not found, skipping {db_name}", message_type="warning"
        )
        return None, None

    if required_motif_names:
        log_message(
            f"Loading {db_name} database (selective: {len(required_motif_names)} motifs)...", 
            message_type="running"
        )
    else:
        log_message(
            f"Loading {db_name} database from {meme_file}...", message_type="running"
        )
    
    motif_names, motif_matrices = parse_meme_pwm(meme_file)
    
    # Create a set of required motif names for fast lookup
    required_set = set(required_motif_names) if required_motif_names else None
    
    # Filter and match motifs
    motif_name_to_pssm = {}
    all_loaded_names = []
    
    for name, matrix in zip(motif_names, motif_matrices):
        # Filter invalid matrices
        if not matrix:
            continue
        if not all(len(row) == 4 for row in matrix):
            continue
        if len(matrix) == 0:
            continue
        if len(matrix) > 2000:
            continue
        
        # If required_motif_names is provided, check if this motif matches
        if required_set:
            matched = False
            for target_name in required_set:
                if match_motif_name(target_name, name, db_name):
                    matched = True
                    break
            
            if not matched:
                continue
        
        # Generate log-odds
        background = [0.25, 0.25, 0.25, 0.25]
        pssm = pwm_to_log_odds(matrix, background)
        
        if not is_finite_matrix(pssm):
            continue
        
        # Store using the MEME file name as key (for lookup by matching)
        motif_name_to_pssm[name] = pssm
        all_loaded_names.append(name)
    
    if not motif_name_to_pssm:
        log_message(f"No valid motifs found in {db_name}", message_type="warning")
        return None, None
    
    log_message(
        f"{db_name}: {len(motif_name_to_pssm)} motifs loaded", message_type="success"
    )
    return motif_name_to_pssm, all_loaded_names


# Load motif information from CSV files
log_message("Loading motif information from CSV files...", message_type="running")
human_motif_map = load_motif_info_from_csv(human_csv)
chimp_motif_map = load_motif_info_from_csv(chimp_csv)

# Collect all unique motif names needed for each database
all_required_motifs = {}
for db_name in databases.keys():
    motif_set = set()
    
    # Collect from human CSV
    for (har, tf, db), motif_list in human_motif_map.items():
        if db == db_name:
            motif_set.update(motif_list)
    
    # Collect from chimp CSV
    for (har, tf, db), motif_list in chimp_motif_map.items():
        if db == db_name:
            motif_set.update(motif_list)
    
    all_required_motifs[db_name] = list(motif_set)
    log_message(f"{db_name}: {len(motif_set)} unique motifs needed", message_type="info")

# Load motifs from databases (selective loading)
all_databases_data = {}
for db_name, meme_file in databases.items():
    required_motifs = all_required_motifs.get(db_name, [])
    motif_name_to_pssm, _ = process_database(db_name, meme_file, required_motifs)
    if motif_name_to_pssm is not None:
        all_databases_data[db_name] = motif_name_to_pssm


# === load FASTA ===
def load_seqs(fasta_path):
    return {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(fasta_path, "fasta")
    }


human_seqs = load_seqs(human_fasta)
chimp_seqs = load_seqs(chimp_fasta)

BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def seq_to_int_array(seq: str):
    mapping = BASE_TO_IDX
    out = []
    for ch in seq:
        out.append(mapping.get(ch, -1))
    return out


def reverse_complement_int_array(seq_int: list):
    # the reverse complement of the integer array (A<->T, C<->G; -1 remains -1)
    comp_map = {0: 3, 1: 2, 2: 1, 3: 0}
    return [comp_map.get(v, -1) for v in reversed(seq_int)]


def best_pwm_score_one_strand_int(seq_int: list, log_pwm):
    w = len(log_pwm)
    n = len(seq_int)
    if n < w:
        return float("-inf")
    best = float("-inf")
    for i in range(n - w + 1):
        score = 0.0
        valid = True
        for j in range(w):
            idx = seq_int[i + j]
            if idx < 0:
                valid = False
                break
            score += log_pwm[j][idx]
        if valid and score > best:
            best = score
    return best


def best_pwm_score_both_strands(seq: str, log_pwm):
    s = seq.upper()
    seq_int = seq_to_int_array(s)
    rc_int = reverse_complement_int_array(seq_int)
    fwd = best_pwm_score_one_strand_int(seq_int, log_pwm)
    rev = best_pwm_score_one_strand_int(rc_int, log_pwm)
    return max(fwd, rev)


rows = []
start_time = time.time()

# Collect all unique HAR-TF pairs from both CSV files
all_har_tf_pairs = set()
for key in human_motif_map.keys():
    all_har_tf_pairs.add((key[0], key[1]))  # (har, TF)
for key in chimp_motif_map.keys():
    all_har_tf_pairs.add((key[0], key[1]))  # (har, TF)

total_pairs = len(all_har_tf_pairs)
log_message(
    f"Starting motif scanning for {total_pairs} HAR-TF pairs across {len(all_databases_data)} databases...",
    message_type="info",
)

done_pairs = 0
for har_id, tf_name in all_har_tf_pairs:
    # Check if sequences exist
    if har_id not in human_seqs or har_id not in chimp_seqs:
        continue
    
    human_seq = human_seqs[har_id]
    chimp_seq = chimp_seqs[har_id]

    # Process each database
    for db_name in databases.keys():
        # Get motif lists from both human and chimp CSV
        human_motifs = human_motif_map.get((har_id, tf_name, db_name), [])
        chimp_motifs = chimp_motif_map.get((har_id, tf_name, db_name), [])
        
        # Combine unique motifs from both sources
        all_motifs_for_pair = set(human_motifs + chimp_motifs)
        
        if not all_motifs_for_pair:
            continue
        
        # Get PSSM dictionary for this database
        db_pssms = all_databases_data.get(db_name, {})
        if not db_pssms:
            continue
        
        # Scan each motif
        for motif_name in all_motifs_for_pair:
            # Find matching PSSM
            log_pwm = None
            # Try exact match first
            if motif_name in db_pssms:
                log_pwm = db_pssms[motif_name]
            else:
                # Try to find by matching
                for pssm_key, pssm_value in db_pssms.items():
                    if match_motif_name(motif_name, pssm_key, db_name):
                        log_pwm = pssm_value
                        break
            
            if log_pwm is None:
                continue
            
            # Scan sequences
            best_human = best_pwm_score_both_strands(human_seq, log_pwm)
            best_chimp = best_pwm_score_both_strands(chimp_seq, log_pwm)
            
            if best_human == float("-inf") and best_chimp == float("-inf"):
                continue
            
            diff = best_human - best_chimp
            direction = (
                "gain_in_human"
                if diff > 0.5
                else ("gain_in_chimp" if diff < -0.5 else "similar")
            )
            rows.append(
                [
                    har_id,
                    tf_name,
                    motif_name,  # Simplified motif name from CSV
                    db_name,  # Database name
                    best_human,
                    best_chimp,
                    round(diff, 3),
                    direction,
                ]
            )

    done_pairs += 1
    if done_pairs % 1000 == 0 or done_pairs == total_pairs:
        elapsed = time.time() - start_time
        log_message(
            f"Processed {done_pairs}/{total_pairs} HAR-TF pairs in {elapsed:.1f}s",
            message_type="info",
        )

os.makedirs(os.path.dirname(output_csv), exist_ok=True)
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(
        [
            "HAR_ID",
            "TF_Name",
            "Motif_Name",
            "Database",
            "Human_Score",
            "Chimp_Score",
            "Score_Diff",
            "Direction",
        ]
    )
    writer.writerows(rows)

log_message(
    f"Motif comparison complete. Results saved to: {output_csv}", message_type="success"
)
log_message(f"Total results: {len(rows)}", message_type="success")
