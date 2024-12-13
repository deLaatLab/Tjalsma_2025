#!/bin/bash

# Script: run_mageck_mle.sh
# Description: Reproducible MAGeCK MLE workflow for CRISPR screen analysis
# Author: [Your Name]
# Date: [Today’s Date]

# Exit on errors
set -e

# Parameters
ENV_NAME="mageck_env"
COUNT_TABLE="ST_CRISPRi_combined_240814.tsv"
DESIGN_MATRIX="design_matrix.txt"
CONTROL_SGRNA="nontargeting.txt"
OUTPUT_PREFIX="ST_CRISPRi_240814_Matrix_simple_mTSS.mle"
WORK_DIR="/home/p.krijger_cbs-niob.local/projects/Sjoerd/screen/MAGeCK/wo_permutation"

# Step 1: Setup Environment
echo "Setting up MAGeCK environment..."
mamba create -n "$ENV_NAME" -y
mamba activate "$ENV_NAME"
mamba install bioconda::mageck -y

# Step 2: Verify Installation
echo "Verifying MAGeCK installation..."
mageck --version

# Step 3: Run MAGeCK MLE
echo "Running MAGeCK MLE..."
cd "$WORK_DIR"

mageck mle \
  -k "$COUNT_TABLE" \
  --design-matrix "$DESIGN_MATRIX" \
  --norm-method control \
  --control-sgrna "$CONTROL_SGRNA" \
  --output-prefix "$OUTPUT_PREFIX"

echo "MAGeCK MLE analysis complete. Output files prefixed with: $OUTPUT_PREFIX"
