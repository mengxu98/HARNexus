#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Hi-C validation of HAR target genes using 4DN contact list-combined (pairs)

Steps:
1. Generate gene promoters (TSS ±2kb, strand-aware) from GTF
2. Load HAR regions
3. Scan Hi-C pairs at read-level
4. Identify HAR–gene promoter contacts
5. Output HAR–gene Hi-C support table
"""

import gzip
import os
import sys
import time
import pandas as pd
from intervaltree import IntervalTree, Interval
from tqdm import tqdm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message


# https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
GTF_FILE = "data/genome/Homo_sapiens.GRCh38.115.gtf.gz"
HAR_FILE = "data/genome/HARs_PMID40011774.csv"  # chr start end HAR_ID
# https://data.4dnucleome.org/files-processed/4DNFIQWVV324/
PAIRS_FILE = "data/hic/4DNFIQWVV324.pairs.gz"  # contact list-combined (pairs)

PROMOTER_FLANK = 5000  # ±2 kb
OUT_FILE = "results/hic/HAR_gene_HiC_supported.csv"


script_start_time = time.time()
log_message("Starting Hi-C validation of HAR target genes...", message_type="info")


# Step 1: Generate gene promoters

log_message("Generating gene promoters from GTF...", message_type="running")
start_time = time.time()

gtf = pd.read_csv(GTF_FILE, sep="\t", comment="#", header=None, low_memory=False)

# Keep only gene entries
gtf = gtf[gtf[2] == "gene"].copy()

gtf["gene_name"] = gtf[8].str.extract('gene_name "([^"]+)"')

def get_promoter(row, flank=PROMOTER_FLANK):
    if row[6] == "+":
        tss = row[3]
    else:
        tss = row[4]
    start = max(0, tss - flank)
    end = tss + flank
    return start, end


gtf[["prom_start", "prom_end"]] = gtf.apply(get_promoter, axis=1, result_type="expand")

gtf[0] = gtf[0].astype(str)
gtf[0] = gtf[0].apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)

promoters = gtf[[0, "prom_start", "prom_end", "gene_name"]].copy()
promoters.columns = ["chr", "start", "end", "gene"]

promoters = promoters[promoters["gene"].notna() & (promoters["gene"] != "")].copy()

elapsed = time.time() - start_time
log_message(
    f"Promoters generated: {len(promoters):,} ({elapsed:.1f}s)", message_type="success"
)
log_message(
    f"Promoter chromosome format sample: {promoters['chr'].unique()[:5].tolist()}",
    message_type="info",
)


# Step 2: Load HAR regions

log_message("Loading HAR regions...", message_type="running")
start_time = time.time()

har_df = pd.read_csv(HAR_FILE)
har_df = har_df[["chr_hg38", "start_hg38", "end_hg38", "HAR_ID"]].copy()
har_df.columns = ["chr", "start", "end", "har_id"]

har_coords = {}
for _, r in har_df.iterrows():
    har_coords[r.har_id] = f"{r.chr}:{r.start}-{r.end}"

elapsed = time.time() - start_time
log_message(
    f"HAR regions loaded: {len(har_df):,} ({elapsed:.1f}s)", message_type="success"
)
log_message(
    f"HAR chromosome format sample: {har_df['chr'].unique()[:5].tolist()}",
    message_type="info",
)


# Step 3: Build interval trees

log_message("Building interval trees...", message_type="running")
start_time = time.time()

har_trees = {}
for _, r in har_df.iterrows():
    har_trees.setdefault(r.chr, IntervalTree()).add(Interval(r.start, r.end, r.har_id))

prom_trees = {}
for _, r in promoters.iterrows():
    prom_trees.setdefault(r.chr, IntervalTree()).add(Interval(r.start, r.end, r.gene))

elapsed = time.time() - start_time
log_message(f"Interval trees ready ({elapsed:.1f}s)", message_type="success")
har_chrs = set(har_trees.keys())
prom_chrs = set(prom_trees.keys())
common_chrs = har_chrs & prom_chrs
log_message(
    f"Common chromosomes between HAR and promoters: {len(common_chrs)}",
    message_type="info",
)


# Step 4: Scan Hi-C pairs at read-level

log_message("Scanning Hi-C pairs (read-level)...", message_type="running")
start_time = time.time()

contact_dict = {}
total_pairs = 0

with gzip.open(PAIRS_FILE, "rt") as f:
    for line in tqdm(f, desc="Hi-C pairs"):
        if line.startswith("#"):
            continue

        cols = line.split()
        chr1, pos1 = cols[1], int(cols[2])
        chr2, pos2 = cols[3], int(cols[4])
        total_pairs += 1

        # end1 = HAR, end2 = promoter
        if chr1 in har_trees and chr2 in prom_trees:
            for h in har_trees[chr1][pos1]:
                for g in prom_trees[chr2][pos2]:
                    key = (h.data, g.data)
                    contact_dict[key] = contact_dict.get(key, 0) + 1

        # end2 = HAR, end1 = promoter
        if chr2 in har_trees and chr1 in prom_trees:
            for h in har_trees[chr2][pos2]:
                for g in prom_trees[chr1][pos1]:
                    key = (h.data, g.data)
                    contact_dict[key] = contact_dict.get(key, 0) + 1

elapsed = time.time() - start_time
log_message(
    f"Total Hi-C pairs scanned: {total_pairs:,} ({elapsed:.1f}s)", message_type="info"
)
log_message(f"Unique HAR–gene contacts: {len(contact_dict):,}", message_type="success")

# Step 5: Output results

log_message("Generating output table...", message_type="running")
start_time = time.time()

os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)

out_df = pd.DataFrame(
    [(h, g, c) for (h, g), c in contact_dict.items()],
    columns=["HAR_ID", "Gene", "n_HiC_validation"],
)

out_df = out_df[out_df["Gene"].notna() & (out_df["Gene"] != "")].copy()

out_df["HAR"] = out_df["HAR_ID"].map(har_coords)

out_df = out_df[["HAR", "HAR_ID", "Gene", "n_HiC_validation"]].copy()
out_df = out_df.rename(columns={"Gene": "hic_gene"})
out_df = out_df.dropna()

out_df.to_csv(OUT_FILE, sep=",", index=False)

elapsed = time.time() - start_time
log_message(
    f"HAR–gene pairs with Hi-C support: {len(out_df):,} ({elapsed:.1f}s)",
    message_type="success",
)
log_message(f"Results written to: {OUT_FILE}", message_type="success")

total_elapsed = time.time() - script_start_time
log_message(
    f"Hi-C validation complete. Total time: {total_elapsed:.1f}s",
    message_type="success",
)
