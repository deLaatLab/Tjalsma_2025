#!/bin/bash

#counting can be done directly on the bam files w/o filtering and sorting. 
#Read sorting is implemented in featurecounts on the fly and it only incurs minimal time cost

mamba activate featurecounts

[ -d "logs" ] || mkdir logs
[ -d "counts" ] || mkdir counts

sbatch --account=hub_laat --job-name=count_SE --time=3:00:00 --mem=20GB --cpus-per-task=12 ~/bin/run_slurm.sh \
"featureCounts \
-Q 15 \
-T 12 \
-a ../reference/hg38_dEGFP_SV40polyA/d2EGFPSV40polyA_gencode.v44.annotation.gtf \
-o ./counts/bwa_mq15_featureCounts_genecounts_SE_stranded.txt \
-s 2 \
-t gene \
./sorted_bams/*.bam"

mamba deactivate


#sbatch --account=hub_laat --job-name=count_SE --time=3:00:00 --mem=20GB --cpus-per-task=12 ~/bin/run_slurm.sh \
#"featureCounts \
#-Q 15 \
#-T 12 \
#-p \
#--countReadPairs \
#-a ./reference/hg38_dEGFP_SV40polyA/d2EGFPSV40polyA_gencode.v44.annotation.gtf \
#-o ./counts/bwa_mq15_featureCounts_genecounts_pairs_stranded.txt \
#-s 2 \
#-t gene \
#./bams/*.bam"
