#!/bin/bash

conda activate align

[ -d "logs" ] || mkdir logs

for bn in ./sorted_bams/*.coordsorted.bam; do
    en=$(basename ${bn} .coordsorted.bam)
    index_file="${bn}.bai"

    if [ ! -f "$index_file" ]; then
        echo "Indexing ${en}"
        sbatch --account=hub_laat --job-name=index_${en} --time=30:00 --mem=20GB --cpus-per-task=1 ~/bin/run_slurm.sh "samtools index ${bn}"
    else
        echo "Skipping indexing for ${en}, index file already exists."
    fi
done
