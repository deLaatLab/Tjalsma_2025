#!/bin/bash

mamba activate align

[ -d "sorted_bams" ] || mkdir sorted_bams
[ -d "logs" ] || mkdir logs

for bn in ./bams/*.bam; do
    en=$(basename ${bn} .bam)
    sorted_bam="./sorted_bams/${en}.coordsorted.bam"

    if [ ! -f "$sorted_bam" ]; then
        echo "Sorting ${en}"
        sbatch --account=hub_laat --job-name=sort_${en} --time=30:00 --mem=30GB --cpus-per-task=6 ~/bin/run_slurm.sh "samtools sort -m4G --threads 6 ${bn} > $sorted_bam"
    else
        echo "Skipping sorting for ${en}, sorted BAM file already exists."
    fi
done
