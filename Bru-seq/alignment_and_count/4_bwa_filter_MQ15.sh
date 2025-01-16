#!/bin/bash

conda activate align

[ -d "bams_MQ15" ] || mkdir bams_MQ15

for bn in ./sorted_bams/*.coordsorted.bam; do
    en=$(basename ${bn} .coordsorted.bam)
    mq15_bam="./bams_MQ15/${en}.MQ15.bam"

    if [ ! -f "$mq15_bam" ]; then
        echo "Filtering ${en} for MQ15"
        sbatch --account=hub_laat --job-name=filter_${en} --time=30:00 --mem=10GB --cpus-per-task=1 ~/bin/run_slurm.sh "samtools view -h -q 15 -Sb ${bn} > $mq15_bam"
    else
        echo "Skipping filtering for ${en}, MQ15 BAM file already exists."
    fi
done
