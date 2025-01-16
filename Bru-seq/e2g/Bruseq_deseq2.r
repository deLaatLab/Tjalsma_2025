# deseq2 ------------------------------------------------------------------

library(DESeq2)

annoData <- read.table("d2EGFPSV40polyA_gencode.v44.annotation.tsv", header = TRUE, sep = "\t")
countData <- read.table("Bruseq_Tjalsma_2025_counts.tsv", header = TRUE, sep = "\t")


# metadata
meta_data <- data.frame(
    row.names = colnames(countData),
    condition = colnames(countData),
    prot = colnames(countData)
)

meta_data$prot <- sub(pattern = "LDB1.*", replacement = "LDB1", x = meta_data$prot)
meta_data$prot <- sub(pattern = "MED23.*", replacement = "MED23", x = meta_data$prot)
meta_data$prot <- sub(pattern = "RAD21.*", replacement = "RAD21", x = meta_data$prot)


meta_data$treatment <- ifelse(grepl("DMSO", meta_data$condition), "DMSO",
    ifelse(grepl("500NM", meta_data$condition), "500",
        ifelse(grepl("50NM", meta_data$condition), "50",
            ifelse(grepl("IAA", meta_data$condition), "IAA", NA)
        )
    )
)

meta_data$treatProt <- paste0(meta_data$treatment, "_", meta_data$prot)
meta_data$condition <- factor(meta_data$condition)
meta_data$prot <- factor(meta_data$prot)
meta_data$treatment <- factor(meta_data$treatment)
meta_data$treatProt <- factor(meta_data$treatProt)
meta_data$treatment <- relevel(meta_data$treatment, ref = "DMSO")

print(meta_data)


# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(
        countData = as.matrix(countData),
        colData = meta_data,
        design = ~treatProt
)

# Filter such that each control sample has a least 5 counts
mincount_control <- 5
control_samples <- which(meta_data$treatment == "DMSO")
dds_filtered <- dds[apply(counts(dds)[, control_samples], 1, function(x) all(x >= mincount_control)), ]


# Perform DESeq2
dds_filtered <- DESeq(dds_filtered)

# extract results
countFolder <- "Bruseq_Tjalsma_2025"
dir.create(file.path(countFolder, "deseq2"))

treatProts <- unique(meta_data$treatProt)

for (treatProt in treatProts) {
    #treatProt="IAA_RAD21"
    message(treatProt)
    for (treatProt_control in treatProts[treatProts != treatProt]) {
        #treatProt_control="DMSO_RAD21"
        message(treatProt_control)

        results <- results(dds_filtered, contrast = c("treatProt", treatProt, treatProt_control))
        results <- merge(annoData, as.data.frame(results), by.x = "GeneID", by.y = "row.names")
        results <- results[order(results$Chromosome, results$Start, results$End), ]

        write.table(results, file = paste0(countFolder, "/deseq2/", treatProt, "vs", treatProt_control, "_mincount", mincount_control, "_deseq2.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
    }
}
