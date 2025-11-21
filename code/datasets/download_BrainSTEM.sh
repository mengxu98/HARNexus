#!/bin/bash

# Download script for GSE235493 data
# paper: https://doi.org/10.1126/sciadv.adu7944
# data: 
#	fetal whole-brain: https://doi.org/10.5281/zenodo.13879662
#	midbrain atlases: https://doi.org/10.5281/zenodo.13879918
# code: https://doi.org/10.5281/zenodo.15243470

set -e

# Source the enhanced download utilities
source "$(dirname "$0")/../functions/utils.sh"

# Create data directory if it doesn't exist
DATA_DIR="../../data/BrainData/raw/BrainSTEM"

log_message "Starting BrainSTEM data download..."

# Define download list
DOWNLOAD_LIST="
https://zenodo.org/records/17291276/files/fetalMidbrainAtlas.noRNA.rds?download=1|fetalMidbrainAtlas.noRNA.rds|0
https://zenodo.org/records/17291276/files/fetalWholeBrainAtlas.noRNA.rds?download=1|fetalWholeBrainAtlas.noRNA.rds|0
https://zenodo.org/api/records/13879918/files-archive|midbrain_atlases.zip|0
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

# Extract zip files (idempotent)
extract_zip_file() {
    local zip_file="$1"
    local extract_dir="$2"
    local marker_file="$3"
    
    if [ ! -f "$zip_file" ]; then
        log_warning "Zip file ${YELLOW}$(basename "$zip_file")${NC} not found, skipping extraction"
        return 1
    fi
    
    if [ -f "$marker_file" ]; then
        log_message "Archive already extracted (marker exists), skipping extraction: ${CYAN}$(basename "$zip_file")${NC}"
        return 0
    fi
    
    log_message "Extracting ${WHITE}$(basename "$zip_file")${NC}..."
    if unzip -q -o "$zip_file" -d "$extract_dir" 2>/dev/null || unzip -o "$zip_file" -d "$extract_dir"; then
        touch "$marker_file"
        log_success "Successfully extracted ${GREEN}$(basename "$zip_file")${NC}"
        return 0
    else
        log_error "Failed to extract ${RED}$(basename "$zip_file")${NC}"
        return 1
    fi
}

# Extract zip archives
while IFS='|' read -r url filename expected_size; do
    # Skip empty lines and comments
    [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
    
    if [[ "$filename" == *.zip ]]; then
        zip_file="$DATA_DIR/$filename"
        marker_file="$DATA_DIR/.${filename%.zip}_extracted"
        extract_zip_file "$zip_file" "$DATA_DIR" "$marker_file"
    fi
done <<< "$DOWNLOAD_LIST"

# Check extracted data integrity
log_message "Checking extracted data integrity..."
total_files=0
found_files=0

while IFS='|' read -r url filename expected_size; do
    # Skip empty lines and comments
    [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
    
    if [[ "$filename" == *.zip ]]; then
        zip_file="$DATA_DIR/$filename"
        marker_file="$DATA_DIR/.${filename%.zip}_extracted"
        
        if [ -f "$marker_file" ]; then
            total_files=$((total_files + 1))
            # Check if extracted directory exists and has content
            extracted_dir="${filename%.zip}"
            if [ -d "$DATA_DIR/$extracted_dir" ] && [ "$(ls -A "$DATA_DIR/$extracted_dir" 2>/dev/null)" ]; then
                found_files=$((found_files + 1))
                file_count=$(find "$DATA_DIR/$extracted_dir" -type f | wc -l | tr -d ' ')
                log_message "Extracted ${GREEN}$extracted_dir${NC} contains ${CYAN}$file_count${NC} files"
            else
                log_warning "Extracted directory ${YELLOW}$extracted_dir${NC} is empty or missing"
            fi
        fi
    fi
done <<< "$DOWNLOAD_LIST"

if [ $total_files -gt 0 ]; then
    log_message "Data extraction check: ${GREEN}$found_files${NC}/${CYAN}$total_files${NC} archives extracted successfully"
fi

# Clean up temporary files
cleanup_temp_files "$DATA_DIR"

log_success "BrainSTEM data download and verification completed!"
