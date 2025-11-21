#!/bin/bash
set -e

# Source the enhanced download utilities
source "$(dirname "$0")/../functions/utils.sh"

log_message "=============================="
log_message "Setting up environment for HAR–TF motif analysis"
log_message "=============================="

# ----------- 1. Create directory structure -----------
BASE_DIR=$(pwd)
DATA_DIR="${BASE_DIR}/data"
GENOME_DIR="${DATA_DIR}/genome"
MOTIF_DIR="${DATA_DIR}/motifs"

mkdir -p "$GENOME_DIR" "$MOTIF_DIR"

log_success "[1/6] Created directories:"
log_message "   ${CYAN}${GENOME_DIR}${NC}"
log_message "   ${CYAN}${MOTIF_DIR}${NC}"

# ----------- 2. Download GRCh38 genome -----------
cd "$GENOME_DIR"

# Check if final processed file already exists
if [ -f "GRCh38.fa" ]; then
    log_success "[2/6] GRCh38.fa already exists, skipping download and processing"
else
    log_message "[2/6] Downloading GRCh38 (UCSC version, with chr prefix)..."

    # Define download list for genome
    GENOME_DOWNLOAD_LIST="
http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz|GRCh38.fa.gz|0
"

    # Download genome file
    batch_download "$GENOME_DOWNLOAD_LIST" "$GENOME_DIR" 3

    # Extract and rename genome file
    if [ -f "GRCh38.fa.gz" ]; then
        log_message "Extracting GRCh38.fa.gz..."
        gunzip -f GRCh38.fa.gz
        if [ -f "hg38.fa" ]; then
            mv hg38.fa GRCh38.fa
            log_success "Renamed hg38.fa to GRCh38.fa"
        fi
    fi
fi

# Create index
if [ ! -f "GRCh38.fa.fai" ]; then
  log_message "[3/6] Indexing GRCh38.fa ..."
  if command -v samtools &> /dev/null; then
    samtools faidx GRCh38.fa
    log_success "Genome index created: GRCh38.fa.fai"
  else
    log_warning "samtools not found. Please install samtools to create genome index."
    log_message "You can install it with: conda install -c bioconda samtools"
    log_message "Or: brew install samtools"
  fi
else
  log_message "[3/6] GRCh38.fa.fai already exists, skipping indexing."
fi

# Add README file
cat > README.txt <<EOF
Human reference genome GRCh38 (UCSC hg38 version).
Source: UCSC Genome Browser (http://hgdownload.soe.ucsc.edu)
Date: $(date)
File: GRCh38.fa
Index: GRCh38.fa.fai
EOF

log_success "Genome ready: $(ls -lh GRCh38.fa)"

# ----------- 3. Download standard TF motif databases -----------
cd "$MOTIF_DIR"

# Check if motif files already exist
jaspar_exists=false
hocomoco_exists=false

if [ -f "JASPAR2024_CORE_vertebrates.meme" ] || [ -f "JASPAR2022_CORE_vertebrates.meme" ] || [ -f "JASPAR2020_CORE_vertebrates.meme" ] || [ -f "JASPAR_CORE_vertebrates.meme" ]; then
    jaspar_exists=true
    log_message "Found existing JASPAR motif file"
fi

if [ -f "HOCOMOCOv11_core_mono_meme_format.meme" ]; then
    hocomoco_exists=true
    log_message "Found existing HOCOMOCO motif file"
fi

# Only download files that don't exist
MOTIF_DOWNLOAD_LIST=""

if [ "$jaspar_exists" = false ]; then
    log_message "[4/6] JASPAR motif database not found, will attempt download"
    MOTIF_DOWNLOAD_LIST="${MOTIF_DOWNLOAD_LIST}
https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt|JASPAR2024_CORE_vertebrates.meme|0"
else
    log_success "[4/6] JASPAR motif database already exists, skipping download"
fi

if [ "$hocomoco_exists" = false ]; then
    log_message "[4/6] HOCOMOCO motif database not found, will attempt download"
    MOTIF_DOWNLOAD_LIST="${MOTIF_DOWNLOAD_LIST}
https://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme|HOCOMOCOv11_core_mono_meme_format.meme|0"
else
    log_success "[4/6] HOCOMOCO motif database already exists, skipping download"
fi

# Download motif files (continue even if some fail)
if [ -n "$MOTIF_DOWNLOAD_LIST" ]; then
    log_message "[4/6] Downloading missing motif databases..."
    batch_download "$MOTIF_DOWNLOAD_LIST" "$MOTIF_DIR" 3 || true
else
    log_success "[4/6] All motif databases already exist, skipping download"
fi

# If JASPAR download failed and no processed file exists, provide manual download instructions
if [ ! -f "JASPAR2024_CORE_vertebrates.meme" ] && [ ! -f "JASPAR2022_CORE_vertebrates.meme" ] && [ ! -f "JASPAR2020_CORE_vertebrates.meme" ] && [ ! -f "JASPAR_CORE_vertebrates.meme" ]; then
    log_warning "JASPAR download failed. Please download manually:"
    log_message "1. Visit: https://jaspar.elixir.no/downloads/"
    log_message "2. Select 'Vertebrates' → 'Single batch file (txt)' → 'PFMs (non-redundant)' → 'MEME'"
    log_message "3. Save as: ${MOTIF_DIR}/JASPAR2024_CORE_vertebrates.meme"
    log_message "4. Or use the HOCOMOCO database which was downloaded successfully"
fi

# Extract JASPAR file if needed
if [ -f "JASPAR2024_CORE_vertebrates.meme.gz" ]; then
    log_message "Extracting JASPAR2024_CORE_vertebrates.meme.gz..."
    gunzip -f JASPAR2024_CORE_vertebrates.meme.gz
    log_success "JASPAR 2024 motifs extracted"
elif [ -f "JASPAR2022_CORE_vertebrates.meme.gz" ]; then
    log_message "Extracting JASPAR2022_CORE_vertebrates.meme.gz..."
    gunzip -f JASPAR2022_CORE_vertebrates.meme.gz
    log_success "JASPAR 2022 motifs extracted"
elif [ -f "JASPAR2020_CORE_vertebrates.meme.gz" ]; then
    log_message "Extracting JASPAR2020_CORE_vertebrates.meme.gz..."
    gunzip -f JASPAR2020_CORE_vertebrates.meme.gz
    log_success "JASPAR 2020 motifs extracted"
fi

# ----------- 4. Download AnimalTFDB TFBS prediction files -----------
log_message "[5/6] Downloading AnimalTFDB TFBS prediction files..."

# Function to download AnimalTFDB files with browser headers
download_animaltfdb_file() {
    local url="$1"
    local filename="$2"
    local description="$3"
    
    if [ -f "$filename" ] && [ -s "$filename" ]; then
        log_success "$description already exists"
        return 0
    fi
    
    log_message "Downloading $description..."
    log_message "URL: $url"
    
    # Method 1: Try curl with browser headers
    log_message "Trying curl with browser headers..."
    if curl -L \
        -H "User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36" \
        -H "Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8" \
        -H "Accept-Language: en-US,en;q=0.5" \
        -H "Accept-Encoding: gzip, deflate" \
        -H "Connection: keep-alive" \
        -H "Upgrade-Insecure-Requests: 1" \
        --retry 3 \
        --retry-delay 5 \
        --connect-timeout 30 \
        --max-time 300 \
        -o "$filename" "$url" 2>/dev/null; then
        
        if [ -s "$filename" ]; then
            log_success "Successfully downloaded $description using curl"
            return 0
        else
            log_warning "Downloaded file is empty, removing..."
            rm -f "$filename"
        fi
    fi
    
    # Method 2: Try wget with browser headers
    if command -v wget &> /dev/null; then
        log_message "Trying wget with browser headers..."
        if wget \
            --user-agent="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36" \
            --header="Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8" \
            --header="Accept-Language: en-US,en;q=0.5" \
            --header="Accept-Encoding: gzip, deflate" \
            --header="Connection: keep-alive" \
            --header="Upgrade-Insecure-Requests: 1" \
            --retry-connrefused \
            --waitretry=5 \
            --read-timeout=30 \
            --timeout=300 \
            -t 3 \
            -O "$filename" "$url" 2>/dev/null; then
            
            if [ -s "$filename" ]; then
                log_success "Successfully downloaded $description using wget"
                return 0
            else
                log_warning "Downloaded file is empty, removing..."
                rm -f "$filename"
            fi
        fi
    fi
    
    # Method 3: Try to get the actual content by following redirects
    log_message "Trying to get direct content..."
    local temp_file="${filename}.tmp"
    
    # Get the final URL after redirects
    local final_url=$(curl -s -I -L "$url" | grep -i "location:" | tail -1 | cut -d' ' -f2 | tr -d '\r\n')
    if [ -z "$final_url" ]; then
        final_url="$url"
    fi
    
    log_message "Final URL: $final_url"
    
    # Try to download the content directly
    if curl -L -s "$final_url" -o "$temp_file" 2>/dev/null; then
        if [ -s "$temp_file" ]; then
            mv "$temp_file" "$filename"
            log_success "Successfully downloaded $description using direct content method"
            return 0
        else
            rm -f "$temp_file"
        fi
    fi
    
    # Method 4: Try with different curl options
    log_message "Trying curl with different options..."
    if curl -L -s \
        --compressed \
        --location-trusted \
        --retry 3 \
        --retry-delay 5 \
        -o "$filename" "$url" 2>/dev/null; then
        
        if [ -s "$filename" ]; then
            log_success "Successfully downloaded $description using compressed curl"
            return 0
        else
            rm -f "$filename"
        fi
    fi
    
    log_error "Failed to download $description after trying multiple methods"
    log_message "Please try downloading manually:"
    log_message "1. Open: $url"
    log_message "2. Copy the content"
    log_message "3. Save as: $MOTIF_DIR/$filename"
    return 1
}

# Download AnimalTFDB files
download_animaltfdb_file \
    "https://guolab.wchscu.cn/AnimalTFDB4_static/download/TFBS_file/fourdatabase_all.meme" \
    "fourdatabase_all.meme" \
    "fourdatabase_all.meme (Combined motif database)"

download_animaltfdb_file \
    "https://guolab.wchscu.cn/AnimalTFDB4_static/download/TFBS_file/hTFtarget_prediction.motifs_matrix.meme" \
    "hTFtarget_prediction.motifs_matrix.meme" \
    "hTFtarget_prediction.motifs_matrix.meme (hTFtarget motifs)"

download_animaltfdb_file \
    "https://guolab.wchscu.cn/AnimalTFDB4_static/download/TFBS_file/fourdatabase_human_mouse.annotation" \
    "fourdatabase_human_mouse.annotation" \
    "fourdatabase_human_mouse.annotation (Annotation file)"

download_animaltfdb_file \
    "https://guolab.wchscu.cn/AnimalTFDB4_static/download/TFBS_file/hTFtarget.annotation" \
    "hTFtarget.annotation" \
    "hTFtarget.annotation (hTFtarget annotation)"

# ----------- 5. Create combined motif file -----------
log_message "[6/6] Creating combined motif file..."

if [ -f "combined_motifs.meme" ]; then
    log_success "combined_motifs.meme already exists, skipping creation"
else
    log_message "Creating combined_motifs.meme from available motif files..."
    
    # Start with MEME header
    cat > combined_motifs.meme <<EOF
MEME version 5

ALPHABET= ACGT

strands: + -

Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 T 0.25000

EOF

    # Add JASPAR motifs if available
    if [ -f "JASPAR2024_CORE_vertebrates.meme" ]; then
        log_message "Adding JASPAR2024 motifs to combined file..."
        # Skip the header and add motifs
        tail -n +$(grep -n "^MOTIF" JASPAR2024_CORE_vertebrates.meme | head -1 | cut -d: -f1) JASPAR2024_CORE_vertebrates.meme >> combined_motifs.meme
    elif [ -f "JASPAR2022_CORE_vertebrates.meme" ]; then
        log_message "Adding JASPAR2022 motifs to combined file..."
        tail -n +$(grep -n "^MOTIF" JASPAR2022_CORE_vertebrates.meme | head -1 | cut -d: -f1) JASPAR2022_CORE_vertebrates.meme >> combined_motifs.meme
    fi

    # Add HOCOMOCO motifs if available
    if [ -f "HOCOMOCOv11_core_mono_meme_format.meme" ]; then
        log_message "Adding HOCOMOCO motifs to combined file..."
        tail -n +$(grep -n "^MOTIF" HOCOMOCOv11_core_mono_meme_format.meme | head -1 | cut -d: -f1) HOCOMOCOv11_core_mono_meme_format.meme >> combined_motifs.meme
    fi

    log_success "Created combined_motifs.meme"
fi

# ----------- 6. Create documentation -----------
cat > README_motifs.txt <<EOF
Motif Databases for TFBS prediction:
====================================

Standard Databases:
- JASPAR2024_CORE_vertebrates.meme  (JASPAR, 2024)
- HOCOMOCOv11_core_mono_meme_format.meme  (Autosome Lab, 2023)
- combined_motifs.meme  (Combined JASPAR + HOCOMOCO)

AnimalTFDB Databases:
- fourdatabase_all.meme: Combined motif database from four databases (JASPAR, HOCOMOCO, TRANSFAC, CISBP)
- hTFtarget_prediction.motifs_matrix.meme: hTFtarget prediction motifs
- fourdatabase_human_mouse.annotation: Annotation file for fourdatabse motifs
- hTFtarget.annotation: Annotation file for hTFtarget motifs

Usage:
1. Use JASPAR2024_CORE_vertebrates.meme with FIMO for standard TFBS prediction
2. Use fourdatabase_all.meme with FIMO for comprehensive TFBS prediction
3. Use hTFtarget_prediction.motifs_matrix.meme with FIMO for hTFtarget prediction
4. Use annotation files to map motif IDs to TF names and species

All files are in MEME motif format and compatible with FIMO, MOODS, and motifmatchr.

Sources:
- JASPAR: https://jaspar.elixir.no/
- HOCOMOCO: https://hocomoco11.autosome.ru/
- AnimalTFDB: https://guolab.wchscu.cn/AnimalTFDB4_static/
Date: $(date)
EOF

log_success "Motif databases ready:"
log_message "Standard motif files:"
ls -lh *.meme 2>/dev/null | grep -E "(JASPAR|HOCOMOCO|combined)" || true
log_message "AnimalTFDB motif files:"
ls -lh fourdatabase_all.meme hTFtarget_prediction.motifs_matrix.meme 2>/dev/null || true
log_message "Annotation files:"
ls -lh *.annotation 2>/dev/null || true

log_message ""
log_message "=============================="
log_success "✅ Environment setup completed."
log_message "Genome:     ${CYAN}${GENOME_DIR}/GRCh38.fa${NC}"
log_message "Standard Motifs: ${CYAN}${MOTIF_DIR}/JASPAR2024_CORE_vertebrates.meme${NC}, ${CYAN}HOCOMOCOv11_core_mono_meme_format.meme${NC}"
log_message "AnimalTFDB Motifs: ${CYAN}${MOTIF_DIR}/fourdatabase_all.meme${NC}, ${CYAN}hTFtarget_prediction.motifs_matrix.meme${NC}"
log_message "=============================="
log_message ""

# Clean up temporary files
cleanup_temp_files "$GENOME_DIR"
cleanup_temp_files "$MOTIF_DIR"

cd "$BASE_DIR"