#!/bin/bash

# Default values
MQ=15
STRAND_SEPARATED=0

# Function to parse parameter value
parse_param() {
    echo "$1" | cut -d '=' -f 2
}

# Check and parse MQ parameter if provided
if [[ $1 == MQ=* ]]; then
    MQ=$(parse_param "$1")
fi

# Check and parse strand parameter if provided
if [[ $2 == strand=* ]]; then
    STRAND=$(parse_param "$2")
    if [ "$STRAND" == "TRUE" ]; then
        STRAND_SEPARATED=1
    elif [ "$STRAND" == "FALSE" ]; then
        STRAND_SEPARATED=0
    fi
fi
mamba activate deeptools

[ -d "logs" ] || mkdir logs
[ -d "bigwigs" ] || mkdir bigwigs

for bn in ./sorted_bams/*.coordsorted.bam; do
    en=$(basename ${bn} .coordsorted.bam)

    # BigWig file paths including MQ in the name
    bw_file="./bigwigs/${en}.BWA.MQ${MQ}.bw"
    fw_bw_file="./bigwigs/${en}.BWA.fw.MQ${MQ}.bw"
    rv_bw_file="./bigwigs/${en}.BWA.rv.MQ${MQ}.bw"

    # Generate unstranded bigWig
    if [ ! -f "$bw_file" ]; then
        sbatch --account=hub_laat --job-name=bw_${en} --time=30:00 --mem=20GB --cpus-per-task=6 ~/bin/run_slurm.sh "bamCoverage -b ${bn} -o ${bw_file} -p 6 --ignoreForNormalization chrM --binSize 10 --minMappingQuality ${MQ} --normalizeUsing RPGC --effectiveGenomeSize 2913022398"
    else
        echo "Skipping bigwig creation for ${en}, file ${bw_file} already exists."
    fi

    # Generate stranded bigWigs if flag is set
    if [ "$STRAND_SEPARATED" -eq 1 ]; then
        if [ ! -f "$fw_bw_file" ]; then
            sbatch --account=hub_laat --job-name=bw_fw_${en} --time=30:00 --mem=20GB --cpus-per-task=6 ~/bin/run_slurm.sh "bamCoverage -b ${bn} -o ${fw_bw_file} -p 6 --ignoreForNormalization chrM --binSize 10 --minMappingQuality ${MQ} --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --filterRNAstrand forward"
        else
            echo "Skipping forward strand bigwig creation for ${en}, file ${fw_bw_file} already exists."
        fi

        if [ ! -f "$rv_bw_file" ]; then
            sbatch --account=hub_laat --job-name=bw_rv_${en} --time=30:00 --mem=20GB --cpus-per-task=6 ~/bin/run_slurm.sh "bamCoverage -b ${bn} -o ${rv_bw_file} -p 6 --ignoreForNormalization chrM --binSize 10 --minMappingQuality ${MQ} --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --filterRNAstrand reverse"
        else
            echo "Skipping reverse strand bigwig creation for ${en}, file ${rv_bw_file} already exists."
        fi
    fi
done