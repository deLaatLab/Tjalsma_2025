#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# ==============================
# Configuration
# ==============================

# URLs for downloading reference files
GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"

# Filenames
GENOME_FA_GZ="GRCh38.primary_assembly.genome.fa.gz"
GTF_FILE_GZ="gencode.v44.annotation.gtf.gz"

# Directory to store reference files
REF_DIR="./reference"
mkdir -p "$REF_DIR"

# Path variables
GENOME_FA_PATH="$REF_DIR/$GENOME_FA_GZ"
GTF_FILE_PATH="$REF_DIR/$GTF_FILE_GZ"
FILTERED_GTF="$REF_DIR/gencode.v44.annotation.filtered.gtf"
FILTERED_GTF_GZ="$FILTERED_GTF.gz"

# Transcript IDs to exclude
EXCLUDE_TRANSCRIPTS="ENST00000642908.1|ENST00000647543.1|ENST00000648735.1|ENST00000380252.6|ENST00000646569.1|ENST00000643199.1|ENST00000380259.7|ENST00000292896.3|ENST00000380237.5|ENST00000380315.2|ENST00000475226.1|ENST00000643122.1|ENST00000429817.1"

# STAR index parameters
GENOME_DIR="/hpc/hub_laat/users/pkrijger/projects/Sjoerd/384/reference/STAR_101/"
THREADS=8
SJDB_OVERHANG=100

# Conda environment
ENV_NAME="STARsolo"

# ==============================
# Functions
# ==============================

# Function to download a file if it doesn't exist
download_file() {
    local url=$1
    local output_path=$2

    if [[ -f "$output_path" ]]; then
        echo "File '$output_path' already exists. Skipping download."
    else
        echo "Downloading $url..."
        wget -O "$output_path" "$url"
        echo "Downloaded '$output_path'."
    fi
}

# Function to filter GTF file
filter_gtf() {
    echo "Filtering GTF file to remove specified globin transcripts..."
    gunzip -c "$GTF_FILE_PATH" | grep -vE "transcript_id \"($EXCLUDE_TRANSCRIPTS)\"" > "$FILTERED_GTF"
    echo "Filtering complete. Compressing the filtered GTF file..."
    gzip "$FILTERED_GTF"
    echo "Filtered GTF file compressed to '$FILTERED_GTF_GZ'."
}

# Function to set up STAR environment using mamba
setup_star_env() {
    echo "Creating mamba environment '$ENV_NAME' with STAR..."
    mamba create -y -n "$ENV_NAME" -c bioconda star=2.7.11a
    echo "Environment '$ENV_NAME' created."
}

# Function to activate the environment and run STAR
run_star() {
    echo "Activating environment '$ENV_NAME'..."
    source "$(mamba info --base)/etc/profile.d/mamba.sh"
    mamba activate "$ENV_NAME"

    echo "Generating STAR genome index..."
    STAR --runThreadN "$THREADS" \
         --runMode genomeGenerate \
         --genomeDir "$GENOME_DIR" \
         --genomeFastaFiles "$GENOME_FA_PATH" \
         --sjdbGTFfile "$FILTERED_GTF_GZ" \
         --sjdbOverhang "$SJDB_OVERHANG"
    echo "STAR genome index generation complete."
}

# ==============================
# Main Script
# ==============================

echo "Starting STAR index setup..."

# Step 1: Download reference files
download_file "$GENOME_URL" "$GENOME_FA_PATH"
download_file "$GTF_URL" "$GTF_FILE_PATH"

# Step 2: Filter the GTF annotation file
filter_gtf

# Step 3: Set up STAR environment
setup_star_env

# Step 4: Run STAR to generate genome index
run_star

echo "STAR index setup completed successfully."
