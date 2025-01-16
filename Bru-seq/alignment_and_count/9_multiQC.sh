#!/bin/bash

# Activate the mamba environment for multiqc
mamba activate qc

# Create a directory for multiqc output if it doesn't exist
[ -d "multiqc" ] || mkdir multiqc

# Prepare the multiqc command
multiqc_cmd="multiqc ./stats ./preseq ./rseqc ./counts --outdir ./multiqc"

# Submit the multiqc command to Slurm using sbatch
sbatch --account=hub_laat --job-name=multiqc_analysis --time=10:00 --mem=10GB --cpus-per-task=1 ~/bin/run_slurm.sh "$multiqc_cmd"

# Deactivate the environment
conda deactivate
