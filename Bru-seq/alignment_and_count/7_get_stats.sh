#!/bin/bash

# Activate the mamba environment
mamba activate align

# Create directories for stats files if they don't exist
[ -d "stats" ] || mkdir stats

# Function to submit jobs for generating stats
submit_stats_job() {
    local bn=$1
    local prefix=$2
    local en=$(basename ${bn} .bam)
#    en=${en%_Aligned.sortedByCoord.out}  # Remove '_Aligned.sortedByCoord.out' for STAR bams
    en=${en%.coordsorted}               # Remove '.coordsorted' for BWA bams
    en=${en%.MQ15}               # Remove '.MQ15' for BWA 15 bams


    local stats_file="./stats/${en}.${prefix}.stats"
    local flagstat_file="./stats/${en}.${prefix}.flagstat"
    local idxstats_file="./stats/${en}.${prefix}.idxstats"

    # Submit job for samtools stats if the file does not exist
    if [ ! -f "$stats_file" ]; then
        sbatch --account=hub_laat --job-name=stats_${prefix}_${en} --time=10:00 --mem=10GB --cpus-per-task=1 ~/bin/run_slurm.sh "samtools stats ${bn} > $stats_file"
    else
        echo "Stats file for $en already exists, skipping."
    fi

    # Submit job for samtools flagstat if the file does not exist
    if [ ! -f "$flagstat_file" ]; then
        sbatch --account=hub_laat --job-name=flagstat_${prefix}_${en} --time=10:00 --mem=5GB --cpus-per-task=1 ~/bin/run_slurm.sh "samtools flagstat ${bn} > $flagstat_file"
    else
        echo "Flagstat file for $en already exists, skipping."
    fi

    # Submit job for samtools idxstats if the file does not exist
    if [ ! -f "$idxstats_file" ]; then
        sbatch --account=hub_laat --job-name=idxstats_${prefix}_${en} --time=10:00 --mem=5GB --cpus-per-task=1 ~/bin/run_slurm.sh "samtools idxstats ${bn} > $idxstats_file"
    else
        echo "Idxstats file for $en already exists, skipping."
    fi
}


# Generate stats for STAR bams using sbatch
#for bn in ./STAR_bams/*_Aligned.sortedByCoord.out.bam; do
#    submit_stats_job $bn "STAR"
#done

# Generate stats for BWA bams using sbatch
for bn in ./sorted_bams/*.coordsorted.bam; do
    submit_stats_job $bn "BWA"
done

# Generate stats for BWA bams using sbatch
for bn in ./bams_MQ15/*.MQ15.bam; do
    submit_stats_job $bn "BWA_MQ15"
done


# Deactivate the environment
conda deactivate



