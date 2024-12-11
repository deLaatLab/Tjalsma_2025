if (!require("ShortRead", character.only = TRUE)) stop("Package not found: ShortRead")


gRNA<-read.table('/home/p.krijger_cbs-niob.local/projects/Sjoerd/screen/CRISPRi_LaatCustom_librarytable.txt'
                 , header = TRUE)
info<-read.table('/home/p.krijger_cbs-niob.local/projects/Sjoerd/screen/info.tsv'
                 , header = TRUE
                 , sep = '\t'
                 , stringsAsFactors = FALSE)



mergeDF<-gRNA
fastqF <- "/home/p.krijger_cbs-niob.local/projects/Sjoerd/screen/fastq"


for(i in 1:nrow(info)){
  message(i,"/",nrow(info))
  barcode<-info$barcode[i]
  fastq<-info$fastq[i]
  name<-paste0(info$expr_id[i], "_", info$rep_id[i])
  outfastq<-file.path(fastqF, paste0(name, ".fastq.gz"))


  message(name)
  reads <- readFastq(fastq)
  
  message("   Extract reads that contain barcodes")
  demultiplex = srFilter(function(x) {
    substr(sread(x), 1, nchar(barcode)) == barcode
  }, name = "demultiplex")
  
  demux_reads <- reads[demultiplex(reads)]
  
  message("   Write demultiplexed reads to fastq file")
  writeFastq(demux_reads, outfastq, mode = "a")

  message("   Get the barcodes, gRNA starts depending on barcode position")
  start_gRNA=7+nchar(barcode)
  
  #remove reads shorter than 19nt
  demux_reads <- demux_reads[width(demux_reads) >= start_gRNA + 19]
  
  gRNA_bc <- subseq(sread(demux_reads), start = start_gRNA, width = 20)
  
  countDF<-as.data.frame(sort(table(as.character(gRNA_bc))))
  mergeDF<-merge(x = mergeDF, y = countDF, by.x='sequence',by.y='Var1',all.x=TRUE)
  mergeDF$Freq[is.na(mergeDF$Freq)]<-0
  colnames(mergeDF) <- gsub('Freq',name,names(mergeDF))

}

head(mergeDF)

write.table(mergeDF, file = '/home/p.krijger_cbs-niob.local/projects/Sjoerd/screen/ST_CRISPRi_individual_experiments_240809.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

#Make MAKeCK input table
mergeDF <- read.table('/home/p.krijger_cbs-niob.local/projects/Sjoerd/screen/ST_CRISPRi_individual_experiments_240809.tsv', sep = '\t', header = TRUE)
mergeDF$gene_TSS <- paste0(mergeDF$gene, "_", mergeDF$transcripts)
combined_DF <- mergeDF[, c("sgID", "gene_TSS")]

# Loop through each condition to sum relevant columns
conditions <- c("E0_Unsorted","E0_Low", "E0_High", "E50_Low", "E50_High", "EC100_Low", "EC100_High", "E100_Low", "E100_High")

for (condition in conditions) {
  cols <- grep(condition, names(mergeDF), value = TRUE)
  combined_DF[[condition]] <- rowSums(mergeDF[, cols, drop = FALSE])
}

head(combined_DF)


write.table(combined_DF, file = '/home/p.krijger_cbs-niob.local/projects/Sjoerd/screen/ST_CRISPRi_combined_240814.tsv', sep = '\t', quote = FALSE, row.names = FALSE)