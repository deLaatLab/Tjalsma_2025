#!/bin/bash

conda activate align

[ -d "bams" ] || mkdir bams
[ -d "logs" ] || mkdir logs

while IFS=$'\t' read -r name fastq_R1 fastq_R2; do
    if [[ $name != "name" ]]; then
        r1="./fastqs/${fastq_R1}"
        r2="./fastqs/${fastq_R2}"
        bam_file="./bams/${name}.bam"

        if [ ! -f "$bam_file" ]; then
            echo "Mapping ${name}"
            if [ "$fastq_R2" != "NA" ]; then
                sbatch --account=hub_laat --job-name=map_${name} --time=3:00:00 --mem=30GB --cpus-per-task=12 ~/bin/run_slurm.sh "bwa mem -t 12 ../reference/hg38_dEGFP_SV40polyA/BWA/hg38_d2EGFPSV40polyA.fasta ${r1} ${r2} | samtools view -h -b - > $bam_file"
            else
                sbatch --account=hub_laat --job-name=map_${name} --time=3:00:00 --mem=30GB --cpus-per-task=12 ~/bin/run_slurm.sh "bwa mem -t 12 ../reference/hg38_dEGFP_SV40polyA/BWA/hg38_d2EGFPSV40polyA.fasta ${r1} | samtools view -h -b - > $bam_file"
            fi
        else
            echo "Skipping ${name}, output BAM file already exists."
        fi
    fi
done < info.tsv
