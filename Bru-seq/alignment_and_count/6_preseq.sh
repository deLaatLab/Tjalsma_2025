#!/bin/bash

mamba activate qc

[ -d "preseq" ] || mkdir preseq

#for PE use the -P option

for bn in ./bams_MQ15/*.MQ15.bam; do
    en=$(basename ${bn} .MQ15.bam)
    ccurve="./preseq/${en}.ccurve"

    if [ ! -f "$ccurve" ]; then
        echo "Sorting ${en}"
        sbatch --account=hub_laat --job-name=preseq_${en} --time=30:00 --mem=10GB --cpus-per-task=1 ~/bin/run_slurm.sh "preseq lc_extrap -o $ccurve -B ${bn}"
    else
        echo "Skipping complexity curve for ${en}, file already exists."
    fi
done

mamba deactivate