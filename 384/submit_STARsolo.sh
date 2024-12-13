
#!/bin/bash

mamba activate STARsolo
ulimit -n 4096

# Set paths
DEMUX_FILE="/hpc/hub_laat/users/pkrijger/projects/Sjoerd/384/STARsolo/Demux_info_ST_384.tsv"
GENOME_DIR="/hpc/hub_laat/users/pkrijger/projects/Sjoerd/384/reference/STAR_101/"
FASTQ_DIR="/hpc/hub_laat/users/pkrijger/projects/Sjoerd/384/fastq/WDLA105/combined/"
OUTPUT_DIR="/hpc/hub_laat/users/pkrijger/projects/Sjoerd/384/STARsolo/output_NoDedup/"
WHITELIST_DIR="/hpc/hub_laat/users/pkrijger/projects/Sjoerd/384/STARsolo/whitelists/"

# Create the whitelist directory if it doesn't exist
mkdir -p "$WHITELIST_DIR"
mkdir -p "$OUTPUT_DIR"

# Loop through each unique 'Name' in the Demux file
tail -n +2 "$DEMUX_FILE" | cut -f9 | sort | uniq | while read -r NAME; do
  # Extract all barcodes for the current Name and create a whitelist file
  grep -P "\t$NAME" "$DEMUX_FILE" | awk '{print $8}' > "${WHITELIST_DIR}${NAME}_whitelist.txt"

  # Construct the FASTQ file paths
  R1_FASTQ="${FASTQ_DIR}${NAME}_R1.fastq.gz"
  R2_FASTQ="${FASTQ_DIR}${NAME}_R2.fastq.gz"

  # Create a SLURM script for each job
  SLURM_SCRIPT="star_${NAME}.slurm"

  cat <<EOT > "$SLURM_SCRIPT"
#!/bin/bash
#SBATCH --job-name=STARsolo_${NAME}     # Job name
#SBATCH --output=${OUTPUT_DIR}${NAME}_%j.out   # Output file
#SBATCH --error=${OUTPUT_DIR}${NAME}_%j.err    # Error file
#SBATCH --time=6:00:00                # Time limit
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=12             # Number of CPU cores per task
#SBATCH --mem=72G                     # Memory limit


# Activate conda environment
mamba activate STARsolo

#In the --readFilesIn option, the 1st file has to be cDNA read, and the 2nd file has to be the barcode (cell+UMI) read, i.e.
#https://teichlab.github.io/scg_lib_structs/methods_html/CEL-seq_family.html


# Run STARSolo
STAR \
--runThreadN 24 \
--genomeDir $GENOME_DIR \
--readFilesIn $R2_FASTQ $R1_FASTQ \
--readFilesCommand zcat \
--soloType CB_UMI_Simple \
--soloCBwhitelist ${WHITELIST_DIR}${NAME}_whitelist.txt \
--soloCBstart 7 --soloCBlen 8 \
--soloUMIstart 1 --soloUMIlen 6 --soloBarcodeReadLength 0 \
--soloStrand Forward \
--soloFeatures Gene \
--outFileNamePrefix ${OUTPUT_DIR}${NAME}_ \
--soloCBmatchWLtype 1MM_multi \
--soloUMIdedup NoDedup \
--soloMultiMappers Unique EM \
--soloUMIfiltering - \
--soloCellFilter None \
--soloCellReadStats Standard \
--outSAMtype BAM SortedByCoordinate \
--outBAMsortingThreadN 6 \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CB UB GX GN
EOT

  # Submit the job
  sbatch "$SLURM_SCRIPT"
done



