#!/bin/bash

# Download hg38 from Gencode
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz

# Unzip hg38
gunzip GRCh38.primary_assembly.genome.fa.gz

# Combine hg38 and d2EGFP_SV40polyA to create a hybrid genome
cat GRCh38.primary_assembly.genome.fa d2EGFP_SV40polyA.fasta > hg38_d2EGFPSV40polyA.fasta

#double check whether cat worked and all chroms are present
grep "^>" hg38_d2EGFPSV40polyA.fasta > cat_check.log