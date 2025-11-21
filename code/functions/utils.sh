#!/bin/bash

# Enhanced download utilities with resume capability, progress display, and integrity checking
# Usage: source this file in other download scripts

set -e

# Color codes for terminal output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
NC='\033[0m' # No Color

# Logging function
log_message() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo -e "${BLUE}ℹ${NC} [${NC}$timestamp${NC}] $*"
}

# Error logging function
log_error() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo -e "${RED}✗${NC} [${NC}$timestamp${NC}] ${RED}ERROR:${NC} $*" >&2
}

# Success logging function
log_success() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo -e "${GREEN}✓${NC} [${NC}$timestamp${NC}] ${GREEN}$*${NC}"
}

# Warning logging function
log_warning() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo -e "${YELLOW}!${NC} [${NC}$timestamp${NC}] ${YELLOW}WARNING:${NC} $*"
}

# Progress logging function
log_progress() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo -e "${PURPLE}◌${NC} [${NC}$timestamp${NC}] ${PURPLE}$*${NC}"
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
        log_message "File ${GREEN}$filename${NC} already exists, skipping download"
        return 0
    fi
    
    # Check for resume file or existing partial file
    local resume_flag=""
    if [ -f "$resume_file" ]; then
        local resume_size=$(cat "$resume_file")
        if [ "$resume_size" -gt 0 ]; then
            resume_flag="-C -"
            log_warning "Resuming download of ${YELLOW}$filename${NC} from byte ${CYAN}$resume_size${NC}"
        fi
    elif [ -f "$output_file" ]; then
        # If output file exists but no resume file, create resume file
        local existing_size=$(stat -f%z "$output_file" 2>/dev/null || echo "0")
        if [ "$existing_size" -gt 0 ]; then
            echo "$existing_size" > "$resume_file"
            resume_flag="-C -"
            log_warning "Resuming download of ${YELLOW}$filename${NC} from byte ${CYAN}$existing_size${NC}"
        fi
    fi
    
    local attempt=1
    while [ $attempt -le $max_retries ]; do
        log_message "Downloading ${WHITE}$filename${NC} (attempt ${CYAN}$attempt${NC}/${CYAN}$max_retries${NC})..."
        
        # Download with clean output and lightweight progress (every 10s)
        (
            # Obtain total size if available for percentage display
            total_size=$(curl -sI "$url" | awk -F": " 'tolower($1)=="content-length"{print $2}' | tr -d '\r')
            [ -z "$total_size" ] && total_size=0
            connecting_logged=0
            last_logged_mb=-1

            # If resuming, start from the resume position
            if [ -n "$resume_flag" ] && [ -f "$resume_file" ]; then
                resume_size=$(cat "$resume_file")
                if [ "$resume_size" -gt 0 ]; then
                    # Copy existing partial file to temp file for resume
                    if [ -f "$output_file" ]; then
                        cp "$output_file" "$temp_file"
                    fi
                fi
            fi

            curl -L \
                --connect-timeout 30 \
                --max-time 3600 \
                --retry 3 \
                --retry-delay 10 \
                --retry-max-time 300 \
                --fail \
                --show-error \
                --silent \
                $resume_flag \
                -o "$temp_file" \
                "$url" &
            curl_pid=$!
            while kill -0 "$curl_pid" 2>/dev/null; do
                if [ -f "$temp_file" ]; then
                    current_size=$(stat -f%z "$temp_file" 2>/dev/null || echo "0")
                    if [ "$current_size" -eq 0 ]; then
                        if [ "$connecting_logged" -eq 0 ]; then
                            log_progress "Downloading $filename ... ${YELLOW}连接中${NC}"
                            connecting_logged=1
                        fi
                    else
                        current_mb=$((current_size / 1024 / 1024))
                        if [ "$current_mb" -ne "$last_logged_mb" ]; then
                            if [ "$total_size" -gt 0 ]; then
                                percent=$(( current_size * 100 / total_size ))
                                # Color code based on progress
                                if [ "$percent" -lt 25 ]; then
                                    color="${RED}"
                                elif [ "$percent" -lt 50 ]; then
                                    color="${YELLOW}"
                                elif [ "$percent" -lt 75 ]; then
                                    color="${BLUE}"
                                else
                                    color="${GREEN}"
                                fi
                                log_progress "Downloading $filename ... ${WHITE}${current_mb}MB${NC} (${color}${percent}%${NC})"
                            else
                                log_progress "Downloading $filename ... ${WHITE}${current_mb}MB${NC}"
                            fi
                            last_logged_mb=$current_mb
                        fi
                    fi
                fi
                sleep 10
            done
            wait "$curl_pid"
        )
        if [ $? -eq 0 ]; then
            # Download successful
            mv "$temp_file" "$output_file"
            rm -f "$resume_file"
            
            # Display file size information
            local downloaded_size=$(stat -f%z "$output_file" 2>/dev/null || echo "0")
            local size_mb=$((downloaded_size / 1024 / 1024))
            log_success "Successfully downloaded ${GREEN}$filename${NC} (${WHITE}${size_mb}MB${NC})"
            return 0
        else
            # Download failed
            local exit_code=$?
            log_error "Download failed for ${RED}$filename${NC} (exit code: ${RED}$exit_code${NC})"
            
            # Save current file size for resume and move temp file to output for next attempt
            if [ -f "$temp_file" ]; then
                local current_size=$(stat -f%z "$temp_file" 2>/dev/null || echo "0")
                if [ "$current_size" -gt 0 ]; then
                    mv "$temp_file" "$output_file"
                    echo "$current_size" > "$resume_file"
                    log_message "Saved resume position: ${CYAN}$current_size${NC} bytes"
                else
                    rm -f "$temp_file"
                fi
            fi
            
            if [ $attempt -lt $max_retries ]; then
                log_warning "Waiting ${YELLOW}$retry_delay${NC} seconds before retry..."
                sleep $retry_delay
            fi
        fi
        
        attempt=$((attempt + 1))
    done
    
    log_error "Failed to download ${RED}$filename${NC} after ${RED}$max_retries${NC} attempts"
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
        
        log_message "Processing file ${CYAN}$total_count${NC}: ${WHITE}$filename${NC}"
        
        # If file exists, verify first; on failure, remove and re-download once
        if [ -f "$output_file" ]; then
            if verify_file_integrity "$output_file" "$expected_size"; then
                success_count=$((success_count + 1))
                echo ""
                continue
            else
                log_error "Integrity check failed for ${RED}$filename${NC}, re-downloading..."
                rm -f "$output_file"
            fi
        fi

        if download_with_resume "$url" "$output_file" "$max_retries"; then
            if verify_file_integrity "$output_file" "$expected_size"; then
                success_count=$((success_count + 1))
            else
                log_error "Integrity check failed for ${RED}$filename${NC}"
            fi
        else
            log_error "Download failed for ${RED}$filename${NC}"
        fi
        
        echo "" # Add blank line for readability
    done <<< "$download_list"
    
    log_message "Batch download completed: ${GREEN}$success_count${NC}/${CYAN}$total_count${NC} files successful"
    
    if [ $success_count -eq $total_count ]; then
        log_success "All files downloaded successfully!"
        return 0
    else
        log_error "Some files failed to download"
        return 1
    fi
}

# Parallel batch download with concurrency control
batch_download_parallel() {
    local download_list="$1"
    local data_dir="$2"
    local max_retries="${3:-3}"
    local concurrency="${4:-3}"

    log_message "Starting parallel batch download (concurrency=${CYAN}$concurrency${NC}) of $(echo "$download_list" | wc -l) files..."

    local temp_dir="$data_dir/.dl_tmp_$$"
    mkdir -p "$temp_dir"

    # Helper to process a single file (runs in background)
    _perform_single_download() {
        local url="$1"
        local filename="$2"
        local expected_size="$3"
        local data_dir="$4"
        local max_retries="$5"

        local output_file="$data_dir/$filename"
        local status_file="$temp_dir/${filename}.status"

        # If file exists, verify first; on failure, remove and re-download once
        if [ -f "$output_file" ]; then
            if verify_file_integrity "$output_file" "$expected_size"; then
                echo ok > "$status_file"
                return 0
            else
                rm -f "$output_file"
            fi
        fi

        if download_with_resume "$url" "$output_file" "$max_retries"; then
            if verify_file_integrity "$output_file" "$expected_size"; then
                echo ok > "$status_file"
            else
                echo fail > "$status_file"
            fi
        else
            echo fail > "$status_file"
        fi
    }

    # Launch jobs with concurrency control
    local -a pids=()
    local total_count=0

    while IFS='|' read -r url filename expected_size; do
        [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
        total_count=$((total_count + 1))

        # Throttle to concurrency
        while [ ${#pids[@]} -ge $concurrency ]; do
            wait "${pids[0]}" 2>/dev/null || true
            pids=("${pids[@]:1}")
        done

        log_message "Queueing file ${CYAN}$total_count${NC}: ${WHITE}$filename${NC}"
        _perform_single_download "$url" "$filename" "$expected_size" "$data_dir" "$max_retries" &
        pids+=("$!")
    done <<< "$download_list"

    # Wait for remaining jobs
    for pid in "${pids[@]}"; do
        wait "$pid" 2>/dev/null || true
    done

    # Summarize results
    local success_count=$(grep -l "^ok$" "$temp_dir"/*.status 2>/dev/null | wc -l | tr -d ' ')
    local total_status=$(ls "$temp_dir"/*.status 2>/dev/null | wc -l | tr -d ' ')
    [ -z "$total_status" ] && total_status=0

    log_message "Parallel batch download completed: ${GREEN}$success_count${NC}/${CYAN}$total_status${NC} files successful"

    rm -rf "$temp_dir"

    if [ "$success_count" -eq "$total_status" ]; then
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
