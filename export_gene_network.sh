#!/bin/bash
# Usage: bash export_gene_network.sh GENENAME
# Example: bash export_gene_network.sh BPTF

GENE="$1"
if [ -z "$GENE" ]; then
  echo "Please provide a gene name as an argument, e.g.: ./export_gene_network.sh BPTF"
  exit 1
fi

INPUT="data/networks/csv/network_data.csv"
OUTDIR="results"
OUTFILE="$OUTDIR/${GENE}_in_network_data.csv"
CLEANED_OUTFILE="$OUTDIR/${GENE}_in_network_data_cleaned.csv"

mkdir -p "$OUTDIR"

awk -F',' -v gene="$GENE" '($4=="\""gene"\"" || $5=="\""gene"\"") {print $0}' "$INPUT" >"$OUTFILE"

sort "$OUTFILE" | uniq >"$CLEANED_OUTFILE"

echo "Exported and cleaned all network data containing $GENE to: $CLEANED_OUTFILE"
