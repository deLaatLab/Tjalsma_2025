#!/bin/bash

# Activate environment
mamba activate rseqc

# Create directories if they don't exist
[ -d "logs" ] || mkdir logs
[ -d "rseqc" ] || mkdir rseqc

for bamfile in ./bams_MQ15/*.bam
do
    # Extract the base name of the file for naming the output
    base=$(basename $bamfile .bam)

    # Run RSeQC read_distribution.py
    sbatch --account=hub_laat --job-name=RSeQC_MQ15_${base} --time=1:00:00 --mem=5GB --cpus-per-task=1 \
    ~/bin/run_slurm.sh "read_distribution.py -i $bamfile -r /hpc/hub_laat/users/pkrijger/projects/TACL/Bruseq/reference/hg38_RefSeq.bed > ./rseqc/bwa_${base}_read_distribution.txt"
done
