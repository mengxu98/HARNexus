#!/usr/bin/env python3
"""
Download CCDS current human data and extract Public CCDS genes
"""

import pandas as pd
import urllib.request
import os
import sys
import csv

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message

data_dir = "data/ccds_genes"
os.makedirs(data_dir, exist_ok=True)

url = "https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt"
local_file = os.path.join(data_dir, "CCDS.current.txt")
output_file = os.path.join(data_dir, "ccds_genes.csv")

if os.path.exists(local_file):
    log_message(
        "CCDS data file already exists: {.path ", local_file, "}", message_type="info"
    )
else:
    log_message(
        "Downloading CCDS data from {.path ", url, "}...", message_type="running"
    )
    urllib.request.urlretrieve(url, local_file)
    log_message(
        "Downloaded CCDS data to {.path ", local_file, "}", message_type="success"
    )

log_message("Reading CCDS data...", message_type="running")
with open(local_file, "rt") as f:
    header = f.readline().strip().lstrip("#").split("\t")
    df = pd.read_csv(f, sep="\t", names=header)

log_message("CCDS data shape: {.val ", df.shape, "}", message_type="info")

log_message("Extracting Public CCDS genes...", message_type="running")
ccds_genes = df.loc[df["ccds_status"] == "Public", "gene"].dropna().unique()

log_message(
    "Number of Public CCDS genes: {.val ", len(ccds_genes), "}", message_type="info"
)

log_message("Saving CCDS genes to {.path ", output_file, "}...", message_type="running")
ccds_genes_df = pd.DataFrame({"gene": ccds_genes})
ccds_genes_df.to_csv(output_file, index=False, quoting=csv.QUOTE_NONE, escapechar="\\")

log_message(
    "Saved {.val ",
    len(ccds_genes),
    "} genes to {.path ",
    output_file,
    "}",
    message_type="success",
)
