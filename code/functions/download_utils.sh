#!/bin/bash

# Enhanced download utilities with resume capability, progress display, and integrity checking
# Usage: source this file in other download scripts

set -e

# Logging function similar to R's log_message
log_message() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "ℹ [$timestamp] $*"
}

# Error logging function
log_error() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "✗ [$timestamp] ERROR: $*" >&2
}

# Success logging function
log_success() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "✓ [$timestamp] $*"
}

# Download function with resume capability and progress display
download_with_resume() {
    local url="$1"
    local output_file="$2"
    local max_retries="${3:-3}"
    local retry_delay="${4:-10}"
    
    local filename=$(basename "$output_file")
    local temp_file="${output_file}.tmp"
    local resume_file="${output_file}.resume"
    
    # Create directory if it doesn't exist
    mkdir -p "$(dirname "$output_file")"
    
    # Check if file already exists and is complete
    if [ -f "$output_file" ]; then
        log_message "File $filename already exists, skipping download"
        return 0
    fi
    
    # Check for resume file
    local resume_flag=""
    if [ -f "$resume_file" ]; then
        local resume_size=$(cat "$resume_file")
        if [ "$resume_size" -gt 0 ]; then
            resume_flag="-C -"
            log_message "Resuming download of $filename from byte $resume_size"
        fi
    fi
    
    local attempt=1
    while [ $attempt -le $max_retries ]; do
        log_message "Downloading $filename (attempt $attempt/$max_retries)..."
        
        # Download with progress bar and resume capability
        if curl -L \
            --connect-timeout 30 \
            --max-time 3600 \
            --retry 3 \
            --retry-delay 5 \
            --retry-max-time 300 \
            --fail \
            --show-error \
            --progress-bar \
            $resume_flag \
            -o "$temp_file" \
            "$url" 2>&1 | tee /dev/stderr; then
            
            # Download successful
            mv "$temp_file" "$output_file"
            rm -f "$resume_file"
            log_success "Successfully downloaded $filename"
            return 0
        else
            # Download failed
            local exit_code=$?
            log_error "Download failed for $filename (exit code: $exit_code)"
            
            # Save current file size for resume
            if [ -f "$temp_file" ]; then
                local current_size=$(stat -f%z "$temp_file" 2>/dev/null || echo "0")
                echo "$current_size" > "$resume_file"
                log_message "Saved resume position: $current_size bytes"
            fi
            
            if [ $attempt -lt $max_retries ]; then
                log_message "Waiting $retry_delay seconds before retry..."
                sleep $retry_delay
            fi
        fi
        
        attempt=$((attempt + 1))
    done
    
    log_error "Failed to download $filename after $max_retries attempts"
    rm -f "$temp_file"
    return 1
}

# Verify file integrity using file size and basic checks
verify_file_integrity() {
    local file_path="$1"
    local expected_size="${2:-0}"
    
    if [ ! -f "$file_path" ]; then
        log_error "File $file_path does not exist"
        return 1
    fi
    
    local actual_size=$(stat -f%z "$file_path" 2>/dev/null || echo "0")
    
    if [ "$expected_size" -gt 0 ] && [ "$actual_size" -lt "$expected_size" ]; then
        log_error "File $file_path is incomplete (expected: $expected_size, actual: $actual_size)"
        return 1
    fi
    
    # Basic file type check for common formats
    local filename=$(basename "$file_path")
    case "$filename" in
        *.rds)
            # Check if it's a valid RDS file (starts with RDS magic bytes)
            if ! head -c 4 "$file_path" | grep -q "RDS"; then
                log_error "File $file_path may not be a valid RDS file"
                return 1
            fi
            ;;
        *.gz)
            # Check if it's a valid gzip file
            if ! gzip -t "$file_path" 2>/dev/null; then
                log_error "File $file_path is not a valid gzip file"
                return 1
            fi
            ;;
    esac
    
    log_success "File integrity check passed for $filename"
    return 0
}

# Batch download function
batch_download() {
    local download_list="$1"
    local data_dir="$2"
    local max_retries="${3:-3}"
    
    log_message "Starting batch download of $(echo "$download_list" | wc -l) files..."
    
    local success_count=0
    local total_count=0
    
    while IFS='|' read -r url filename expected_size; do
        # Skip empty lines and comments
        [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
        
        total_count=$((total_count + 1))
        local output_file="$data_dir/$filename"
        
        log_message "Processing file $total_count: $filename"
        
        if download_with_resume "$url" "$output_file" "$max_retries"; then
            if verify_file_integrity "$output_file" "$expected_size"; then
                success_count=$((success_count + 1))
            else
                log_error "Integrity check failed for $filename"
            fi
        else
            log_error "Download failed for $filename"
        fi
        
        echo "" # Add blank line for readability
    done <<< "$download_list"
    
    log_message "Batch download completed: $success_count/$total_count files successful"
    
    if [ $success_count -eq $total_count ]; then
        log_success "All files downloaded successfully!"
        return 0
    else
        log_error "Some files failed to download"
        return 1
    fi
}

# Clean up temporary files
cleanup_temp_files() {
    local data_dir="$1"
    log_message "Cleaning up temporary files..."
    find "$data_dir" -name "*.tmp" -delete 2>/dev/null || true
    find "$data_dir" -name "*.resume" -delete 2>/dev/null || true
    log_success "Cleanup completed"
}
