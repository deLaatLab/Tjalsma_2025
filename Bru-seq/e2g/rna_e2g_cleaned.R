library(dplyr)
library(tibble)
library(ggplot2)
library(ggExtra)
library(plyr)
library(reshape) 
library(tidyverse)

###load in ABC table
abc<-read.table("ENCFF202FID.bed.gz", header=FALSE)
colnames(abc)<-colnames(read.table("ENCFF202FID.bed.gz", header=TRUE, nrows = 1, comment.char = ''))
colnames(abc)

###Add BrU-seq DEseq data
RAD21_deseq <-read.delim("IAA_RAD21vsDMSO_RAD21_mincount5_deseq2.tsv")
MED23_deseq <-read.delim("500_MED23vsDMSO_MED23_mincount5_deseq2.tsv")
LDB1_500NM_deseq <-read.delim("500_LDB1vsDMSO_LDB1_mincount5_deseq2.tsv")

rownames(RAD21_deseq) <- RAD21_deseq[, "GeneID"]
rownames(MED23_deseq) <- MED23_deseq[, "GeneID"]
rownames(LDB1_500NM_deseq) <- LDB1_500NM_deseq[, "GeneID"]

#Take columns
RAD21_deseq <- RAD21_deseq[, c("baseMean", "log2FoldChange", "padj")]
MED23_deseq <- MED23_deseq[, c("baseMean", "log2FoldChange",  "padj")]
LDB1_500NM_deseq <- LDB1_500NM_deseq[, c("baseMean", "log2FoldChange",  "padj")]

normcounts <- RAD21_deseq

normcounts$GeneID2 <- gsub("\\..*$", "", rownames(normcounts))

colnames(normcounts) <- c("RAD21_baseMean", "RAD21_log2FoldChange",  "RAD21_padj", "GeneID2")

##Merge different tables
data<-normcounts
data_abc <- merge(data, abc, by.x = "GeneID2", by.y = "TargetGeneEnsemblID", all.x = FALSE)

MED23_deseq$GeneID2 <- gsub("\\..*$", "", rownames(MED23_deseq))
colnames(MED23_deseq) <- c("MED23_baseMean", "MED23_log2FoldChange", "MED23_padj", "GeneID2")
data_abc <- merge(data_abc, MED23_deseq, by = "GeneID2")

LDB1_500NM_deseq$GeneID2 <- gsub("\\..*$", "", rownames(LDB1_500NM_deseq))
colnames(LDB1_500NM_deseq) <- c("LDB1_500NM_baseMean", "LDB1_500NM_log2FoldChange",  "LDB1_500NM_padj", "GeneID2")
data_abc <- merge(data_abc, LDB1_500NM_deseq, by = "GeneID2")

data_selection_500_merged <- data_abc


###LONG RANGE GENES
data_selection_500_merged_filter <- data_selection_500_merged %>%
  filter(distanceToTSS.Feature>=50000 & distanceToTSS.Feature<500000 &Score >=0.8)  

###SHORT RANGE GENES
data_selection_500_merged_short <- data_selection_500_merged %>%
  filter(distanceToTSS.Feature>2000 & distanceToTSS.Feature<10000 & Score>=0.8)  

data_selection_500_merged_stringent_nomid <- data_selection_500_merged %>%
  filter(distanceToTSS.Feature>=10000 & Score>=0.3)  

#filter short range genes not to have a long-range enhancer
data_selection_500_merged_short <- data_selection_500_merged_short %>%
  filter(!(TargetGene %in% data_selection_500_merged_stringent_nomid$TargetGene))

###MID RANGE GENES
data_selection_500_merged_mid <- data_selection_500_merged %>%
  filter(distanceToTSS.Feature>10000 & distanceToTSS.Feature<40000 & Score>=0.8) 

data_selection_500_merged_stringent_nofar <- data_selection_500_merged %>%
  filter(distanceToTSS.Feature>=40000 & Score>=0.3)  

#filter mid range genes not to have a long-range enhancer
data_selection_500_merged_mid <- data_selection_500_merged_mid %>%
  filter(!(TargetGene %in% data_selection_500_merged_filter$TargetGene)&!(TargetGene %in% data_selection_500_merged_stringent_nofar$TargetGene))

###NO ENHANCER GENES
data_selection_500_merged_all_enhancers <- data_selection_500_merged %>%
  filter(distanceToTSS.Feature>2000 & Score>=0.3)  

data_selection_500_merged_none_enhancers <- data_selection_500_merged %>%
  filter(!(TargetGene %in% data_selection_500_merged_all_enhancers$TargetGene))  


##Select genes from MA table
##Define gene categories
long_genes <- data_selection_500_merged[data_selection_500_merged$TargetGene %in% data_selection_500_merged_filter$TargetGene ,]
mid_genes <- data_selection_500_merged[data_selection_500_merged$TargetGene %in% data_selection_500_merged_mid$TargetGene ,]
short_genes <- data_selection_500_merged[data_selection_500_merged$TargetGene %in% data_selection_500_merged_short$TargetGene ,]
no_enhancers <- data_selection_500_merged[data_selection_500_merged$TargetGene %in% data_selection_500_merged_none_enhancers$TargetGene ,]

length(unique(short_genes$GeneID2))
length(unique(mid_genes$GeneID2))
length(unique(long_genes$GeneID2))
length(unique(no_enhancers$GeneID2))


#Here, if there are multiple enhancers > 0.8, pick one with highest ABC score
data_selection_500_topenhancers <- NULL


for (gene in 1:length(unique(long_genes$GeneID2))) {
  gene_name <- unique(long_genes$GeneID2)[gene]
  highest_element <- max(data_selection_500_merged_filter[data_selection_500_merged_filter$GeneID2 == gene_name, "Score"])
  to_add <- data_selection_500_merged_filter[(data_selection_500_merged_filter$GeneID2 == gene_name)&(data_selection_500_merged_filter$Score == highest_element),]
  to_add$Category <- "long"
  data_selection_500_topenhancers <- rbind(data_selection_500_topenhancers, to_add)
}

for (gene in 1:length(unique(short_genes$GeneID2))) {
  gene_name <- unique(short_genes$GeneID2)[gene]
  highest_element <- max(data_selection_500_merged_short[data_selection_500_merged_short$GeneID2 == gene_name, "Score"])
  to_add <- data_selection_500_merged_short[(data_selection_500_merged_short$GeneID2 == gene_name)&(data_selection_500_merged_short$Score == highest_element),]
  to_add$Category <- "short"
  data_selection_500_topenhancers <- rbind(data_selection_500_topenhancers, to_add)
}

for (gene in 1:length(unique(mid_genes$GeneID2))) {
  gene_name <- unique(mid_genes$GeneID2)[gene]
  highest_element <- max(data_selection_500_merged_mid[data_selection_500_merged_mid$GeneID2 == gene_name, "Score"])
  to_add <- data_selection_500_merged_mid[(data_selection_500_merged_mid$GeneID2 == gene_name)&(data_selection_500_merged_mid$Score == highest_element),]
  to_add$Category <- "mid"
  data_selection_500_topenhancers <- rbind(data_selection_500_topenhancers, to_add)
}

for (gene in 1:length(unique(no_enhancers$GeneID2))) {
  gene_name <- unique(no_enhancers$GeneID2)[gene]
  to_add <- data_selection_500_merged_none_enhancers[(data_selection_500_merged_none_enhancers$GeneID2 == gene_name),][1,]
  to_add$Category <- "no_enhancer"
  data_selection_500_topenhancers <- rbind(data_selection_500_topenhancers, to_add)
}

data_selection_500_topenhancers <- data_selection_500_topenhancers[!duplicated(data_selection_500_topenhancers), ]


##set no enhancer genes to not inCCD (because they do not have an EP link)
##this is just for ggplot later
data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="no_enhancer","inCCD.Feature"] <- 0


data_selection_500_topenhancers <- data_selection_500_topenhancers %>%
  filter(RAD21_baseMean > 50)

####plot ubiquitously expressed ratio per group


colnames(data_selection_500_topenhancers)
feature_tocheck <- "ubiquitousExpressedGene.Feature"
summary_data <- data_selection_500_topenhancers %>%
  group_by(Category) %>%
  dplyr::summarise(
    Ubiquitous = sum(eval(parse(text=feature_tocheck)) == 1),
    Nonubiquitous = sum(eval(parse(text=feature_tocheck)) == 0)
  )
expression_data <- as.data.frame(summary_data)
expression_data <- melt(expression_data)

# Create stacked bar plot 
type_order <- c("no_enhancer", "short","mid","long")
colors_bars <- c("#1abc9c","#e74c3c")
p <- expression_data %>%
  group_by(variable) %>%
  group_by(Category) %>%
  ggplot(aes(fill=factor(variable, c("Ubiquitous", "Nonubiquitous")), x=factor(Category, type_order), y=value))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=colors_bars)+
  labs(x= "Category",y= "Proportion of genes ubiquitously expressed", fill="Expression status", title="Cell-type specific expression")+theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position = "none",
        plot.title = element_text(size = 10),   # Title text size
        axis.title.x = element_text(size = 6), # X-axis label text size
        axis.title.y = element_text(size = 6), # Y-axis label text size
        axis.text = element_text(size = 6),
        panel.background = element_rect(fill = "#ffffff", colour="black", linewidth=0.25),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth=0.25),
        axis.line = element_blank())
p



#plotting MAplot, per category no enhancer, short long
sample_tocheck <- "RAD21"

colors_dots <- c("#afafaf", "#f78028" ,"#1abc9c","#0339A5")
alpha <- 1
size_dots <- 3
categories <- c("no_enhancer", "long", "mid","short")
data_selection_500_topenhancers$plotFC <- data_selection_500_topenhancers[,paste0(sample_tocheck, "_log2FoldChange")]
data_selection_500_topenhancers$expression <- log2(data_selection_500_topenhancers[,paste0(sample_tocheck, "_baseMean")])

data_selection_500_topenhancers <- data_selection_500_topenhancers %>%
  arrange(factor(Category, levels= categories))

plot <- ggplot(data_selection_500_topenhancers, aes(x = expression, y = plotFC))+
  geom_hline(yintercept = 0)+
  geom_point(aes(color=factor(Category, categories)),alpha = alpha, show.legend = FALSE, size = size_dots)+
  labs(x = "log2 expression level", y = paste0("log2FC ",sample_tocheck, " depletion"), color="Category")+
  ggtitle(paste0(sample_tocheck," per enhancer category"))+
  scale_colour_manual(values = colors_dots)

plot_colors <- theme(
  panel.background = element_rect(fill = "#ffffff", colour="black", linewidth=0.25),
  panel.grid.major = element_blank(),   # Remove major grid lines
  panel.grid.minor = element_blank()
)
p <- plot+plot_colors

p <- ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)
p


#Plot cumulative distributions
colors_dots <- c("#afafaf","#0339A5","#1abc9c", "#f78028" )
categories <- c("no_enhancer", "short", "mid","long")
data_selection_500_topenhancers <- data_selection_500_topenhancers %>%
  arrange(factor(Category, levels =categories))

#Draw the cdf plot (choose with x= which factor to plot)
x_min <- -2  # set your desired minimum
x_max <- 2   # set your desired maximum
plot <- ggplot(data = data_selection_500_topenhancers, aes(x = RAD21_log2FoldChange, 
                                                           group = factor(Category, categories), 
                                                           col = factor(Category, categories)))+
  geom_vline(xintercept = 0,linetype = "dashed")+
  geom_hline(yintercept = 0)+
  stat_ecdf(linewidth=0.5)+
  labs(title="Comparison EP genes",y="Cumulative distribution fraction", color="EP category")+
  scale_colour_manual(values = colors_dots)+
  scale_y_continuous(expand=c(0, 0), limits=c(0, 1))+
  scale_fill_manual("Category")+
  theme(legend.position = "none",
        plot.title = element_text(size = 10),   # Title text size
        axis.title.x = element_text(size = 6), # X-axis label text size
        axis.title.y = element_text(size = 6), # Y-axis label text size
        axis.text = element_text(size = 6),
        panel.background = element_rect(fill = "#ffffff", colour="black", linewidth=0.25),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth=0.25),
        axis.line = element_blank())+
  coord_cartesian(xlim = c(x_min, x_max)) 

plot

###Plot log2FCs per category
colors_dots <- c("#f78028" ,"#1abc9c","#afafaf","#0339A5")
categories <- c("no_enhancer", "short", "mid","long")

p <- data_selection_500_topenhancers %>%
  group_by(Category) %>%
  ggplot(aes(x=factor(Category, categories), y=RAD21_log2FoldChange, fill=Category))+
  geom_text(aes(label=after_stat(count)), y=2, stat='count', colour="black", size=2)+
  geom_hline(yintercept = 0)+
  geom_violin(show.legend = FALSE)+geom_boxplot(width=0.2, fill = "white",show.legend = FALSE, size=0.25,lwd=0.1,outlier.size=0.25)+
  labs(x = "Enhancer-promoter category", title="log2FCs")+
  scale_fill_manual(values = colors_dots)+ylim(-3,3)

plot <- p + theme(legend.position = "none",
           plot.title = element_text(size = 10),   # Title text size
           axis.title.x = element_text(size = 6), # X-axis label text size
           axis.title.y = element_text(size = 6), # Y-axis label text size
           axis.text = element_text(size = 6),
           panel.background = element_rect(fill = "#ffffff", colour="black", linewidth=0.25),
           panel.grid.major = element_blank(),   # Remove major grid lines
           panel.grid.minor = element_blank(),
           axis.ticks = element_line(linewidth=0.25),
           axis.line = element_blank())
plot

#Statistics
feature <- "RAD21_log2FoldChange"

##No enhancer vs short
wilcox.test(data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="no_enhancer",feature],
        data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="short",feature])
##No enhancer vs mid
wilcox.test(data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="no_enhancer",feature],
            data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="mid",feature])
##No enhancer vs long
wilcox.test(data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="no_enhancer",feature],
            data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="long",feature])
##Short vs long
wilcox.test(data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="short",feature],
        data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="long",feature])


##Plot ratio of genes in down vs not
##Choose with factor what to plot

table_tocheck <- data_selection_500_topenhancers
factor <- "RAD21"

table_tocheck$Type <- paste0(table_tocheck$inCCD.Feature,table_tocheck$Category)

down <- table_tocheck[(table_tocheck[,paste0(factor,"_padj")] < 0.05) & (table_tocheck[,paste0(factor,"_log2FoldChange")] < -0.8),]
stable <- table_tocheck[((table_tocheck[,paste0(factor,"_padj")] > 0.05) &(table_tocheck[,paste0(factor,"_log2FoldChange")] > -0.2))|
                          ((table_tocheck[,paste0(factor,"_padj")] > 0.05) &(table_tocheck[,paste0(factor,"_log2FoldChange")] < 0.2)),]

effect <- data.frame(matrix(nrow = 4, ncol = 0))
effect$down <- as.vector(table(down$Category))

effect$stable <-as.vector(table(stable$Category))
effect$Category <- names(table(down$Category))
effect$Category <- c("Long","Mid","No enhancer", "Short")
effect <- melt(effect)
colnames(effect) <- c("Category","affected", "genes")
affected_order <- c("stable", "down")
category_order <- c("No enhancer", "Short","Mid","Long")
effect$Category <- factor(effect$Category, levels=category_order)

colors_dots <- c("#afafafaf" ,"#0339A5","#1abc9c","#f78028")
effect <- effect %>%
  arrange(factor(affected, levels = affected_order))
effect <- effect %>%
  arrange(factor(Category, levels = category_order))
p<- effect %>%
  group_by(affected) %>%
  group_by(Category) %>%
  ggplot(aes(x=factor(affected, affected_order), fill=factor(Category), y=genes)) + 
  geom_bar(position="fill", stat="identity", width = 0.8, show.legend = FALSE)+
  scale_fill_manual(values = colors_dots)+
  labs(x= "Effect after depletion",y= "Proportion of genes present")+
  ggtitle(paste0(factor, ". Stable N=", nrow(stable), " Down N=", nrow(down)))+theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position = "none",
        plot.title = element_text(size = 10),   # Title text size
        axis.title.x = element_text(size = 6), # X-axis label text size
        axis.title.y = element_text(size = 6), # Y-axis label text size
        axis.text = element_text(size = 6),
        panel.background = element_rect(fill = "#ffffff", colour="black", linewidth=0.25),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth=0.25),
        axis.line = element_blank())
p


##Look at CTCF-Contact Domains (CCD)
###See if long-range genes contact any inCCD enhancer
data_selection_500_CCD <- NULL
features_toselect <- c("TargetGene","inCCD.Feature", "normalizedEP300_enhActivity.Feature",
                       "RAD21_log2FoldChange", "LDB1_500NM_log2FoldChange", 
                       "RAD21_baseMean",  "LDB1_500NM_baseMean", 
                       "RAD21_padj",  "LDB1_500NM_padj")

for (gene in 1:length(unique(long_genes$TargetGene))) {
  gene_name <- unique(long_genes$TargetGene)[gene]
  enhancers <- (data_selection_500_merged_filter[data_selection_500_merged_filter$TargetGene == gene_name, ])
  
  enhancers_CCD <- enhancers[enhancers$inCCD.Feature ==1,]
  if (sum(enhancers_CCD$inCCD.Feature) >= 1){
    highest_element <- max(enhancers_CCD[enhancers_CCD$TargetGene == gene_name, "Score"])
    to_add <- data_selection_500_merged_filter[(data_selection_500_merged_filter$TargetGene == gene_name)&(data_selection_500_merged_filter$Score == highest_element),features_toselect]
    to_add$Category <- "long"
  } else {
    highest_element <- max(enhancers[enhancers$TargetGene == gene_name, "Score"])
    to_add <- data_selection_500_merged_filter[(data_selection_500_merged_filter$TargetGene == gene_name)&(data_selection_500_merged_filter$Score == highest_element),features_toselect]
    to_add$Category <- "long"
  }
  data_selection_500_CCD <- rbind(data_selection_500_CCD, to_add)
}
##same for mid-range genes
for (gene in 1:length(unique(mid_genes$TargetGene))) {
  gene_name <- unique(mid_genes$TargetGene)[gene]
  enhancers <- (data_selection_500_merged_mid[data_selection_500_merged_mid$TargetGene == gene_name, ])
  
  enhancers_CCD <- enhancers[enhancers$inCCD.Feature ==1,]
  if (sum(enhancers_CCD$inCCD.Feature) >= 1){
    highest_element <- max(enhancers_CCD[enhancers_CCD$TargetGene == gene_name, "Score"])
    to_add <- data_selection_500_merged_mid[(data_selection_500_merged_mid$TargetGene == gene_name)&(data_selection_500_merged_mid$Score == highest_element),features_toselect]
    to_add$Category <- "mid"
  } else {
    highest_element <- max(enhancers[enhancers$TargetGene == gene_name, "Score"])
    to_add <- data_selection_500_merged_mid[(data_selection_500_merged_mid$TargetGene == gene_name)&(data_selection_500_merged_mid$Score == highest_element),features_toselect]
    to_add$Category <- "mid"
  }
  data_selection_500_CCD <- rbind(data_selection_500_CCD, to_add)
}
##same for short-range genes
for (gene in 1:length(unique(short_genes$TargetGene))) {
  gene_name <- unique(short_genes$TargetGene)[gene]
  enhancers <- (data_selection_500_merged_short[data_selection_500_merged_short$TargetGene == gene_name, ])
  
  enhancers_CCD <- enhancers[enhancers$inCCD.Feature ==1,]
  if (sum(enhancers_CCD$inCCD.Feature) >= 1){
    highest_element <- max(enhancers_CCD[enhancers_CCD$TargetGene == gene_name, "Score"])
    to_add <- data_selection_500_merged_short[(data_selection_500_merged_short$TargetGene == gene_name)&(data_selection_500_merged_short$Score == highest_element),features_toselect]
    to_add$Category <- "short"
  } else {
    highest_element <- max(enhancers[enhancers$TargetGene == gene_name, "Score"])
    to_add <- data_selection_500_merged_short[(data_selection_500_merged_short$TargetGene == gene_name)&(data_selection_500_merged_short$Score == highest_element),features_toselect]
    to_add$Category <- "short"
  }
  data_selection_500_CCD <- rbind(data_selection_500_CCD, to_add)
}

###Same expression filter as above
data_selection_500_CCD <- data_selection_500_CCD %>%
    filter(RAD21_baseMean > 50)

##Compare RAD21 log2FC in or out CCD 
table_tocheck <- data_selection_500_CCD 

colors_dots <- c("#9e9ac8","#6a51a3")
categories <- c("short", "mid","long")
plot <- table_tocheck %>%
  group_by(Category) %>%
  ggplot(aes(x=factor(Category, categories), y=RAD21_log2FoldChange, fill=as.factor(inCCD.Feature)))+
  geom_hline(yintercept = 0)+
  geom_text(aes(label=after_stat(count)), y=2, colour="black", size=2, stat="count",position=position_dodge2(width=0.8))+
  geom_boxplot(position = position_dodge(width = 0.8), width=0.6,lwd=0.1,outlier.size=0.25)+
  labs(x = "Enhancer category", y="log2FC RAD21 depletion", title="inCCD filter")+
  scale_fill_manual(values = colors_dots)+
  scale_x_discrete(expand = expansion(mult = c(0.25, 0.25))) +
  theme(legend.position = "none",
        plot.title = element_text(size = 10),   # Title text size
        axis.title.x = element_text(size = 6), # X-axis label text size
        axis.title.y = element_text(size = 6), # Y-axis label text size
        axis.text = element_text(size = 6),
        panel.background = element_rect(fill = "#ffffff", colour="black", linewidth=0.25),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth=0.25),
        axis.line = element_blank())

plot


###Compare features between no enhancer, short and long-range
table_tocheck <- data_selection_500_topenhancers

colnames(table_tocheck)

colors_dots <-  c("#f78028","#1abc9c","#afafafaf","#0339A5")
plot <- table_tocheck %>%
  group_by(Category) %>%
  ggplot(aes(x=factor(Category, levels = c("no_enhancer","short", "mid","long")), y=normalizedEP300_enhActivity.Feature, fill=as.factor(Category)))+
  geom_hline(yintercept = 0)+
  geom_boxplot(width=0.8, show.legend = FALSE,lwd=0.1,outlier.size=0.25)+
  labs(x = "Enhancer-promoter category")+
  scale_fill_manual(values = colors_dots)+
  scale_y_continuous(expand = expansion(add = c(0, 0)))+
  theme(legend.position = "none",
        plot.title = element_text(size = 10),   # Title text size
        axis.title.x = element_text(size = 6), # X-axis label text size
        axis.title.y = element_text(size = 6), # Y-axis label text size
        axis.text = element_text(size = 6),
        panel.background = element_rect(fill = "#ffffff", colour="black", linewidth=0.25),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth=0.25),
        axis.line = element_blank())
plot


##only enhancer feature
###features without checking effect, diff between no enhancer, short and long-range
table_tocheck <- data_selection_500_topenhancers
table_tocheck <- table_tocheck %>%
  filter(!(Category =="no_enhancer"))
colnames(table_tocheck)

colors_dots <-  c("#f78028","#1abc9c","#0339A5")
plot <- table_tocheck %>%
  group_by(Category) %>%
  ggplot(aes(x=factor(Category, levels = c("short", "mid","long")), y=normalizedEP300_enhActivity.Feature, fill=as.factor(Category)))+
  geom_hline(yintercept = 0)+
  geom_boxplot(width=0.8, show.legend = FALSE,lwd=0.1,outlier.size=0.25)+
  labs(x = "Enhancer-promoter category")+
  scale_fill_manual(values = colors_dots)+
  scale_y_continuous(expand = expansion(add = c(0, 0.01)))+
  theme(legend.position = "none",
      plot.title = element_text(size = 10),   # Title text size
      axis.title.x = element_text(size = 6), # X-axis label text size
      axis.title.y = element_text(size = 6), # Y-axis label text size
      axis.text = element_text(size = 6),
      panel.background = element_rect(fill = "#ffffff", colour="black", linewidth=0.25),
      panel.grid.major = element_blank(),   # Remove major grid lines
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(linewidth=0.25),
      axis.line = element_blank())
plot


##Plot features on depletion effect
table_tocheck <- data_selection_500_topenhancers
colnames(table_tocheck)
factor <- "RAD21"
feature <- "normalizedEP300_enhActivity.Feature"

down <- table_tocheck[(table_tocheck[,paste0(factor,"_padj")] < 0.05) & (table_tocheck[,paste0(factor,"_log2FoldChange")] <  -0.8),]
stable <- table_tocheck[((table_tocheck[,paste0(factor,"_padj")] > 0.05) &(table_tocheck[,paste0(factor,"_log2FoldChange")] > -0.2))|
                          ((table_tocheck[,paste0(factor,"_padj")] > 0.05) &(table_tocheck[,paste0(factor,"_log2FoldChange")] < 0.2)),]

long_down <- down[down$Category == "long", "TargetGene"]
long_stable <- stable[stable$Category == "long", "TargetGene"]

long_down_feature <- table_tocheck[table_tocheck$TargetGene %in% long_down,c(feature, "TargetGene")]
long_stable_feature <- table_tocheck[table_tocheck$TargetGene %in% long_stable,c(feature, "TargetGene")]

long_down_feature$Category <- "down_long"
long_stable_feature$Category <- "stable_long"

colnames(long_down_feature) <- c("feature", "TargetGene", "Category")
colnames(long_stable_feature) <- c("feature", "TargetGene", "Category")

data_compare <- rbind(long_down_feature, long_stable_feature)
data_compare$Category <- as.factor(data_compare$Category)

##Plot  levels on long-range genes stable vs long
colors_dots <- c("#fd8d3c","#fee08b")
plot <- data_compare %>%
  ggplot(aes(fill=Category, y=feature, x=factor(Category, levels = c("stable_long", "down_long"))))+
  geom_violin(show.legend = TRUE)+geom_boxplot(width=0.1, fill = "white",show.legend = FALSE,lwd=0.1,outlier.size=0.25)+
  scale_fill_manual(values = colors_dots)+
  ggtitle(paste0(factor))+
  geom_text(aes(label=after_stat(count)), y=85, colour="black", size=4, stat="count",position=position_dodge2(width=0.9))+
  labs(x="Category and effect", y=feature)+theme(legend.position = "none",
  axis.title.x = element_text(size = 16), # X-axis label text size
  axis.title.y = element_text(size = 16), # Y-axis label text size
  axis.text = element_text(size = 12),
  panel.background = element_rect(fill = "#ffffff", colour="black", linewidth=0.25),
  panel.grid.major = element_blank(),   # Remove major grid lines
  panel.grid.minor = element_blank(),
  axis.ticks = element_line(linewidth=0.25),
  axis.line = element_blank())
plot



##Make quartiles and plot RAD21 effect
#Select genes with defined enhancer
table_tocheck <- data_selection_500_topenhancers
colnames(table_tocheck)
table_tocheck <- table_tocheck %>%
  filter(!(Category =="no_enhancer"))

##Quartile the level on these enhancers by a specific feature
quantiles <- 4
feature <- "normalizedEP300_enhActivity.Feature"
table_tocheck <- table_tocheck %>%
  mutate(quantile=as.factor(ntile(eval(parse(text=feature)), quantiles)))

#compare lowest with highest quartile
table_tocheck <- table_tocheck %>%
  filter(quantile==1|quantile==quantiles)

##add back no enhancer
data_selection_500_topenhancers_noenhancer <- data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="no_enhancer",]
data_selection_500_topenhancers_noenhancer$quantile <- 0
table_tocheck <- rbind(table_tocheck,data_selection_500_topenhancers_noenhancer)

categories <- c("no_enhancer","short","mid","long")

##Filter only long-range genes
colors_dots <-  c("#9e9ac8","#6a51a3")

###Plot RAD21 log2FC on EP300 quartiles
plot <- table_tocheck %>%
  group_by(Category) %>%
    ggplot(aes(x=factor(Category, categories), y=RAD21_log2FoldChange, fill=as.factor(quantile)))+
  geom_hline(yintercept = 0, linewidth=0.25)+
  geom_text(aes(label=after_stat(count)), y=1.2, colour="black", size=4, stat="count",position=position_dodge2(width=0.9))+
  geom_boxplot(position = position_dodge(width = 0.8), width=0.7,lwd=0.1,outlier.size=0.25)+
  labs(x = "Enhancer-promoter category", y="log2FC RAD21 depletion", title="EP300 levels")+
  scale_fill_manual(values = colors_dots)+
  scale_x_discrete(expand = expansion(mult = c(0.25, 0.25))) +
  theme(legend.position = "none",
        plot.title = element_text(size = 10),   # Title text size
        axis.title.x = element_text(size = 6), # X-axis label text size
        axis.title.y = element_text(size = 6), # Y-axis label text size
        axis.text = element_text(size = 6),
        panel.background = element_rect(fill = "#ffffff", colour="black", linewidth=0.25),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth=0.25),
        axis.line = element_blank())

plot

#plot quartiles but grouped other way around
colors_dots <-  c("#0339A5","#1abc9c","#f78028")

table_tocheck <- table_tocheck %>%
  filter(!(Category =="no_enhancer"))


###Plot RAD21 log2FC on quartiles
plot <- table_tocheck %>%
  group_by(Category) %>%
  ggplot(aes(fill=factor(Category, categories), y=RAD21_log2FoldChange, x=as.factor(quantile)))+
  geom_hline(yintercept = 0, linewidth=0.25)+
  geom_text(aes(label=after_stat(count)), y=1.2, colour="black", size=4, stat="count",position=position_dodge2(width=0.9))+
  geom_boxplot(position = position_dodge(width = 0.8), width=0.6, lwd=0.05, outlier.size=0.25)+
  labs(x = "Enhancer-promoter category", y="log2FC depletion", title=feature)+
  scale_fill_manual(values = colors_dots)+
  scale_x_discrete(expand = expansion(mult = c(0.25, 0.25))) +
  theme(legend.position = "none",
        plot.title = element_text(size = 12),   # Title text size
        axis.title.x = element_text(size = 12), # X-axis label text size
        axis.title.y = element_text(size = 12), # Y-axis label text size
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "#ffffff", colour="black", linewidth=0.25),
        panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth=0.25),
        axis.line = element_blank())

plot


###Statistics
factor <- "RAD21_log2FoldChange"

wilcox.test(data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="no_enhancer", factor],
            table_tocheck[table_tocheck$Category=="short"&table_tocheck$quantile=="1", factor])
wilcox.test(data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="no_enhancer", factor],
            table_tocheck[table_tocheck$Category=="mid"&table_tocheck$quantile=="1", factor])
wilcox.test(data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="no_enhancer", factor],
            table_tocheck[table_tocheck$Category=="long"&table_tocheck$quantile=="1", factor])

wilcox.test(data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="no_enhancer", factor],
            table_tocheck[table_tocheck$Category=="short"&table_tocheck$quantile=="4", factor])
wilcox.test(data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="no_enhancer", factor],
            table_tocheck[table_tocheck$Category=="mid"&table_tocheck$quantile=="4", factor])
wilcox.test(data_selection_500_topenhancers[data_selection_500_topenhancers$Category=="no_enhancer", factor],
            table_tocheck[table_tocheck$Category=="long"&table_tocheck$quantile=="4", factor])



