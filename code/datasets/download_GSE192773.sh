#!/bin/bash

# Download script for GSE192773 data
# paper:
#   https://doi.org/10.1038/s41586-023-06338-4
#   https://doi.org/10.1038/s41467-025-60665-w
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192773
# code:

set -e

# Source the enhanced download utilities
source "$(dirname "$0")/../functions/utils.sh"

# Create data directory if it doesn't exist
DATA_DIR="../../data/BrainData/raw/GSE192773"

log_message "Starting GSE192773 data download..."

# Define download list
DOWNLOAD_LIST="
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FAll%5FSpecies%5FUMI%5FCount%5FMatrix%5FBarcodes%2Etxt%2Egz|GSE192773_All_Species_UMI_Count_Matrix_Barcodes.txt.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FAll%5FSpecies%5FUMI%5FCount%5FMatrix%5FGenes%2Etxt%2Egz|GSE192773_All_Species_UMI_Count_Matrix_Genes.txt.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FAll%5FSpecies%5FUMI%5FCount%5FMatrix%5FMetaData%2Etxt%2Egz|GSE192773_All_Species_UMI_Count_Matrix_MetaData.txt.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FAll%5FSpecies%5FUMI%5FCount%5FMatrix%5FsnRNAseq%2Emtx%2Egz|GSE192773_All_Species_UMI_Count_Matrix_snRNAseq.mtx.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FSeq%5FTemplate%5FEvolutionPaper%5FRNA%5FUpdated%2Exlsx|GSE192773_Seq_Template_EvolutionPaper_RNA_Updated.xlsx|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FSeurat%5FAstrocyte%5FRNA%2ERDS%2Egz|GSE192773_Seurat_Astrocyte_RNA.RDS.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FSeurat%5FExcitatory%5FRNA%2ERDS%2Egz|GSE192773_Seurat_Excitatory_RNA.RDS.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FSeurat%5FInhibitory%5FRNA%2ERDS%2Egz|GSE192773_Seurat_Inhibitory_RNA.RDS.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FSeurat%5FMicroglia%5FRNA%2ERDS%2Egz|GSE192773_Seurat_Microglia_RNA.RDS.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FSeurat%5FOPCOli%5FRNA%2ERDS%2Egz|GSE192773_Seurat_OPCOli_RNA.RDS.gz|0
"

# Function to decompress gzip files
decompress_gz_file() {
    local gz_file="$1"
    local decompressed_file="${gz_file%.gz}"
    
    # Check if compressed file exists
    if [ ! -f "$gz_file" ]; then
        log_warning "Compressed file ${YELLOW}$(basename "$gz_file")${NC} not found, skipping decompression"
        return 1
    fi
    
    # Check if decompressed file already exists
    if [ -f "$decompressed_file" ]; then
        log_message "Decompressed file ${GREEN}$(basename "$decompressed_file")${NC} already exists, skipping decompression"
        return 0
    fi
    
    # Decompress the file
    log_message "Decompressing ${WHITE}$(basename "$gz_file")${NC}..."
    if gunzip -c "$gz_file" > "$decompressed_file"; then
        local size=$(stat -f%z "$decompressed_file" 2>/dev/null || echo "0")
        local size_mb=$((size / 1024 / 1024))
        log_success "Successfully decompressed ${GREEN}$(basename "$decompressed_file")${NC} (${WHITE}${size_mb}MB${NC})"
        return 0
    else
        log_error "Failed to decompress ${RED}$(basename "$gz_file")${NC}"
        return 1
    fi
}

# Perform parallel batch download (set concurrency as needed)
batch_download_parallel "$DOWNLOAD_LIST" "$DATA_DIR" 5 3

# Decompress all .gz files
log_message "Starting decompression of .gz files..."
compressed_count=0
decompressed_count=0

while IFS='|' read -r url filename expected_size; do
    # Skip empty lines and comments
    [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
    
    if [[ "$filename" == *.gz ]]; then
        compressed_count=$((compressed_count + 1))
        gz_file="$DATA_DIR/$filename"
        if decompress_gz_file "$gz_file"; then
            decompressed_count=$((decompressed_count + 1))
        fi
    fi
done <<< "$DOWNLOAD_LIST"

log_message "Decompression completed: ${GREEN}$decompressed_count${NC}/${CYAN}$compressed_count${NC} files decompressed"

# Clean up temporary files
cleanup_temp_files "$DATA_DIR"

log_success "GSE192773 data download and decompression completed!"