#!/bin/bash

# title: Human prefrontal cortex gene regulatory dynamics from gestation to adulthood at single-cell resolution
# paper: https://doi.org/10.1016/j.cell.2022.09.039
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/36318921
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168408
# code:

set -e

# Source the enhanced download utilities
source "$(dirname "$0")/../functions/utils.sh"

# Create data directory if it doesn't exist
DATA_DIR="../../data/BrainData/raw/GSE168408"

log_message "Starting GSE168408 data download..."

# Define download list
DOWNLOAD_LIST="
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE168408&format=file|GSE168408.tar|0
"

# Perform batch download
batch_download_parallel "$DOWNLOAD_LIST" "$DATA_DIR" 10 10

# Verify downloaded files
log_message "Verifying downloaded files..."
verification_failed=0

while IFS='|' read -r url filename expected_size; do
    # Skip empty lines and comments
    [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
    
    file_path="$DATA_DIR/$filename"
    if [ ! -f "$file_path" ]; then
        log_error "File ${RED}$filename${NC} not found after download"
        verification_failed=$((verification_failed + 1))
        continue
    fi
    
    if ! verify_file_integrity "$file_path" "$expected_size"; then
        verification_failed=$((verification_failed + 1))
    fi
done <<< "$DOWNLOAD_LIST"

if [ $verification_failed -gt 0 ]; then
    log_error "File verification failed for ${RED}$verification_failed${NC} file(s)"
    exit 1
fi

log_success "All files verified successfully!"

# Extract archive files (idempotent)
extract_archive_file() {
    local archive_file="$1"
    local extract_dir="$2"
    local marker_file="$3"
    
    if [ ! -f "$archive_file" ]; then
        log_warning "Archive file ${YELLOW}$(basename "$archive_file")${NC} not found, skipping extraction"
        return 1
    fi
    
    if [ -f "$marker_file" ]; then
        log_message "Archive already extracted (marker exists), skipping extraction: ${CYAN}$(basename "$archive_file")${NC}"
        return 0
    fi
    
    log_message "Extracting ${WHITE}$(basename "$archive_file")${NC}..."
    
    if [[ "$archive_file" == *.tar || "$archive_file" == *.tar.gz ]]; then
        if tar -xf "$archive_file" -C "$extract_dir" 2>/dev/null; then
            touch "$marker_file"
            log_success "Successfully extracted ${GREEN}$(basename "$archive_file")${NC}"
            return 0
        else
            log_error "Failed to extract ${RED}$(basename "$archive_file")${NC}"
            return 1
        fi
    elif [[ "$archive_file" == *.zip ]]; then
        if unzip -q -o "$archive_file" -d "$extract_dir" 2>/dev/null || unzip -o "$archive_file" -d "$extract_dir"; then
            touch "$marker_file"
            log_success "Successfully extracted ${GREEN}$(basename "$archive_file")${NC}"
            return 0
        else
            log_error "Failed to extract ${RED}$(basename "$archive_file")${NC}"
            return 1
        fi
    else
        log_warning "Unsupported archive format: ${YELLOW}$(basename "$archive_file")${NC}"
        return 1
    fi
}

# Extract archive files
log_message "Extracting archive files..."
while IFS='|' read -r url filename expected_size; do
    [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
    
    if [[ "$filename" == *.tar || "$filename" == *.tar.gz || "$filename" == *.zip ]]; then
        archive_file="$DATA_DIR/$filename"
        marker_file="$DATA_DIR/.${filename//\//_}_extracted"
        extract_archive_file "$archive_file" "$DATA_DIR" "$marker_file"
    fi
done <<< "$DOWNLOAD_LIST"

# Check extracted data integrity
log_message "Checking extracted data integrity..."
total_archives=0
successful_extractions=0

while IFS='|' read -r url filename expected_size; do
    [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
    
    if [[ "$filename" == *.tar || "$filename" == *.tar.gz || "$filename" == *.zip ]]; then
        archive_file="$DATA_DIR/$filename"
        marker_file="$DATA_DIR/.${filename//\//_}_extracted"
        
        if [ -f "$marker_file" ]; then
            total_archives=$((total_archives + 1))
            file_count=$(find "$DATA_DIR" -maxdepth 2 -type f ! -name "$filename" ! -name ".*" ! -name "*.tmp" ! -name "*.resume" 2>/dev/null | wc -l | tr -d ' ')
            if [ "$file_count" -gt 0 ]; then
                successful_extractions=$((successful_extractions + 1))
                log_message "Archive ${GREEN}$(basename "$filename")${NC} extracted, found ${CYAN}$file_count${NC} files"
            else
                log_warning "No extracted content found for ${YELLOW}$(basename "$filename")${NC}"
            fi
        fi
    fi
done <<< "$DOWNLOAD_LIST"

if [ $total_archives -gt 0 ]; then
    log_message "Extraction check: ${GREEN}$successful_extractions${NC}/${CYAN}$total_archives${NC} archives extracted successfully"
fi

# Clean up temporary files
cleanup_temp_files "$DATA_DIR"

log_success "GSE168408 data download and verification completed!"
