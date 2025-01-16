#!/bin/bash


# Comprehensive gene annotation (the Comprehensive and Basic annotations contain the same number of genes, but comprehensive more transcripts/exons..)
#STAR manual: The use of the most comprehensive annotations for a given species is strongly recommended
#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
#gunzip gencode.v44.annotation.gtf.gz


# Input GTF and output GTF
input_gtf="gencode.v44.annotation.gtf"
output_gtf="d2EGFPSV40polyA_gencode.v44.annotation.gtf"

# Copy existing GTF to new file
cp $input_gtf $output_gtf

# Add annotations for d2EGFP_SV40polyA
echo -e "d2EGFP_SV40polyA\tsynthetic\tgene\t1\t1050\t.\t+\t.\tgene_id \"d2EGFP_SV40polyA\"; gene_type \"synthetic_gene\"; gene_name \"d2EGFP_SV40polyA\";" >> $output_gtf
echo -e "d2EGFP_SV40polyA\tsynthetic\ttranscript\t1\t1050\t.\t+\t.\tgene_id \"d2EGFP_SV40polyA\"; transcript_id \"d2EGFP_SV40polyA.1\"; gene_name \"d2EGFP_SV40polyA\";" >> $output_gtf
echo -e "d2EGFP_SV40polyA\tsynthetic\texon\t1\t1050\t.\t+\t.\tgene_id \"d2EGFP_SV40polyA\"; transcript_id \"d2EGFP_SV40polyA.1\"; exon_number \"1\"; gene_name \"d2EGFP_SV40polyA\";" >> $output_gtf
