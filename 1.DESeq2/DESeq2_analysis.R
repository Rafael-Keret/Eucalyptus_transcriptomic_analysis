
# DESeq2 pipeline for differential gene expression analysis of transcriptomic data from Eucalyptus grandis.
# Author: Rafael Keret

# LOAD PACKAGES
 
library("DESeq2")
library("dplyr")
library("ggrepel")
library("stringr")

# IMPORT AND PREPARE FEATURECOUNTS DATA FILE (i.e. counts matrix)

fc_output <- read.csv("./Data/input/FC_NCBI_GTF", head = TRUE, sep = "\t", skip = 1)
fc_output <- fc_output[, c(1, 7:14)]
colnames(fc_output) <- c("qseqid", "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4")
fc_output$qseqid <- gsub("LOC", "", fc_output$qseqid)

length(unique(fc_output$qseqid))

fc_output <- data.frame(fc_output, row.names = 1)

# (1) Set factor or treatment level(s) from pheno data file (i.e. what you are testing against in your experiment)

DESeq_pheno_data <- read.csv("./Data/input/DESeq_pheno_data.csv", head = TRUE, sep = ",", row.names = 1)

DESeq_pheno_data$treatment <- factor(DESeq_pheno_data$treatment)

# (2) Create a DESeq object with read counts and pheno data
 # If you have multiple factors add these to "design =" below
 # Then you need to specify your factors in the order of increasing importance
 # For instance, if you have species and drought, but you are more interested in the overall impact of drought
 # Then specify species as your first factor, and drought as your second (i.e. "design = ~ species + drought")

dds <- DESeqDataSetFromMatrix(countData = fc_output, colData = DESeq_pheno_data, design = ~treatment)

# (3) Specify the reference for treatment factor
 # Set control (well-watered) as a reference (i.e. droughted divided by control to calculate fold change)
 # The one that is specified first becomes the reference which to compare to 

dds$treatment <- factor(dds$treatment, levels = c("C", "D"))

# (4) Filter for genes that have read counts (counts per million) greater than or equal to 5
 # Soewarto et.al, 2019 M&M provide a good explanation of this

keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

# STATISTICS

# (1) Identify differentially expressed genes (DEGs)
 # The DESeq function automatically normalizes for sequencing depth and composition before performing DEA
 # The method is called "median of ratios"

dds <- DESeq(dds)
deseq_result <- results(dds)
deseq_result

# (2) Look at a summary of our results

summary(deseq_result)  # This summarizes the data based on a p-value < 0.1

deseq_result_0.01 <- results(dds, alpha = 0.01) # can change it to a lower p-value = 0.01
summary(deseq_result_0.01)

# (3) Change deseq_result into dataframe, and order on p-value

deseq_result <- as.data.frame(deseq_result)
deseq_result <- deseq_result[order(deseq_result$pvalue),]
write.csv(deseq_result, "./Data/deseq_result.csv")

# (4) Create background gene list to retrieve FASTA sequences or gene models from NCBI

gene_list <- rownames(deseq_result)

write.table(gene_list, "./Data/ncbi_gene_list.txt", sep="\t", row.names = FALSE, 
            quote = FALSE, col.names = FALSE)

# (5) Some queries that can be used
 # Ask if a specific gene is DE and down or up regulated? 
 # Pick an arbitrary gene for this example (LOC104456901), and call it from deseq_result

deseq_result["104456901",]

# (6) Extract the most differentially expressed genes
 # Select genes with a q-value (i.e. adjusted p-value) < 0.05 and a |log2FoldChange| > 1 (i.e. 2x expression)

deseq_filtered <- deseq_result %>% filter(deseq_result$padj < 0.05)
deseq_filtered <- deseq_filtered %>% filter(abs(deseq_filtered$log2FoldChange) > 1)

write.csv(deseq_filtered, "./Data/deseq_filtered.csv")

# (7) Create a differentially expressed gene list

DEG_list <- rownames(deseq_filtered)

write.table(DEG_list, "./Data/ncbi_DEG_list.txt", sep="\t", row.names = FALSE, 
            quote = FALSE, col.names = FALSE)

# (8) Create a normalized dataframe based on library depth and composition

normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts <- as.data.frame(normalized_counts)
write.csv(normalized_counts, "./Data/normalized_counts.csv")

# DISPERSION PLOT

plotDispEsts(dds)

# PCA PLOT

# (1) Variance stabilizing transformation

vsd <- vst(dds, blind = FALSE)

# (2) Use transformed values to generate PCA plot

plotPCA(vsd, intgroup = "treatment")

# HEATMAPS

# Heatmap of sample-to-sample distance matrix (with clustering) based on normalised counts
# To determine which samples as a whole are most related to each other, based on gene expression patterns

# (1) Load packages

library("pheatmap")
library("RColorBrewer")

# (2) Create distance matrix

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# (3) Colour scheme 

colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# (4) Plot 

pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, col = colours)

# (5) Heatmap of log2 transformed normalized counts (i.e. library size & comp)
 # Top 30 differentially expressed genes

top_deg <- deseq_result[order(deseq_result$padj), ][1:30, ]
top_deg <- row.names(top_deg)

rld <- rlog(dds, blind = FALSE)

treatment_col = list(treatment = c(C = "skyblue", D = "orange2"))

pheatmap(assay(rld)[top_deg, ], cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE, 
         annotation_col = DESeq_pheno_data, annotation_colors = treatment_col, 
         color = colorRampPalette(c("black", "green", "red"))(50))

# (6) Heatmap of Z-scores
 # Top 30 differentially expressed genes

cal_z_score <- function(x) {(x - mean(x)) / sd(x)}

z_score_all <- t(apply(normalized_counts, 1, cal_z_score))  

z_score_subset <- z_score_all[top_deg,]

pheatmap(z_score_subset, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE, 
         annotation_col = DESeq_pheno_data, annotation_colors = treatment_col, 
         color=colorRampPalette(c("green", "black", "red"))(50))

# MA PLOT

# (1) Load package and plot raw untransformed MA plot

library("apeglm")

plotMA(dds, ylim = c(-8, 8))

# Remove noise of low counts with high dispersion (background noise)
# This finds instances where the log fold change was exaggerated and shrinks them down
# Previously, we got higher dispersions (i.e. variances) at a lower "mean of normalised counts"
# So the dispersion plot tells us that the fold change will be exaggerated at lower mean counts
# because the variance was higher in genes of lower expression

# (2) Transform via "apeglm" model

res_apeglm_MA <- lfcShrink(dds, coef = "treatment_D_vs_C", type = "apeglm")
plotMA(res_apeglm_MA, ylim = c(-5, 5))

# (3) Transform via "normal" model

res_normal_MA <- lfcShrink(dds, coef = "treatment_D_vs_C", type = "normal")
plotMA(res_normal_MA, ylim = c(-5, 5))

# VOLCANO PLOT APEGLM

# (1) Load packages

library("data.table")
library("ggplot2")
library("ggrepel")
library("dplyr")

# (2) Change the normalized res_apeglm_MA data, to a dataframe

resLFC_ape <- as.data.frame(res_apeglm_MA) %>% setDT(keep.rownames = "gene_id")

# (3) Label genes as up, down or not regulated

resLFC_ape$diffexpressed <- "NO"
resLFC_ape$diffexpressed[resLFC_ape$log2FoldChange > 1 & resLFC_ape$padj < 0.05] <- "UP"
resLFC_ape$diffexpressed[resLFC_ape$log2FoldChange < -1 & resLFC_ape$padj < 0.05] <- "DOWN"

# (4) Create a labeling column for genes of interest

resLFC_ape$delabel <- NA

# (5) Set the threshold p-value at which labels should be included

threshold <- head(arrange(resLFC_ape, pvalue), 10)$pvalue[10]

# (6) Add gene id labels based on values that are lower than or equal to the threshold
 # Start with a query and end with what we want to perform should the query be fulfilled 

resLFC_ape$delabel[resLFC_ape$pvalue <= threshold & !is.na(resLFC_ape$pvalue)] <- (resLFC_ape$gene_id[resLFC_ape$pvalue <= threshold & !is.na(resLFC_ape$pvalue)])

# (7) Plot volcano
 # Noise reduced data

ggplot(data = resLFC_ape, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel)) + 
  geom_point() + theme_minimal() + geom_text_repel() + scale_colour_manual(values = c("green3", "grey", "red3")) +
  theme(text = element_text(size = 12)) + geom_vline(xintercept = c(-1, 1), col = "black", linetype = 2) + 
  geom_hline(yintercept = -log10(0.01), col = "black", linetype = 2)

# VOLCANO PLOT NORMAL

# (1) Change the normalized res_normal_MA data, to a dataframe

resLFC_norm <- as.data.frame(res_normal_MA) %>% setDT(keep.rownames = "gene_id")

# (2) Label genes as up, down or not regulated

resLFC_norm$diffexpressed <- "NO"
resLFC_norm$diffexpressed[resLFC_norm$log2FoldChange > 1 & resLFC_norm$padj < 0.05] <- "UP"
resLFC_norm$diffexpressed[resLFC_norm$log2FoldChange < -1 & resLFC_norm$padj < 0.05] <- "DOWN"

# (3) Label genes of interest

resLFC_norm$delabel <- NA

# (4) Set the threshold for p-value at which labels should be included

threshold <- head(arrange(resLFC_norm, pvalue), 10)$pvalue[10]

# (5) Add gene id labels based on values that are lower than or equal to the threshold
 # start with a query and end with what we want to perform should the query be fulfilled 

resLFC_norm$delabel[resLFC_norm$pvalue <= threshold & !is.na(resLFC_norm$pvalue)] <- (resLFC_norm$gene_id[resLFC_norm$pvalue <= threshold & !is.na(resLFC_norm$pvalue)])

# (6) Plot volcano
 # Noise reduced data

ggplot(data = resLFC_norm, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel)) + 
  geom_point() + theme_minimal() + geom_text_repel() + scale_colour_manual(values = c("green3", "grey", "red3")) +
  theme(text = element_text(size = 12)) + geom_vline(xintercept = c(-1, 1), col = "black", linetype = 2) + 
  geom_hline(yintercept = -log10(0.01), col = "black", linetype = 2)

# RAW VOLCANO PLOT

# (1) Convert to dataframe and change first column name

deseq_result <- as.data.frame(deseq_result) %>% setDT(keep.rownames = "gene_id")

# (2) Label genes as up, down or not regulated

deseq_result$diffexpressed <- "NO"
deseq_result$diffexpressed[deseq_result$log2FoldChange > 1 & deseq_result$padj < 0.05] <- "UP"
deseq_result$diffexpressed[deseq_result$log2FoldChange < -1 & deseq_result$padj < 0.05] <- "DOWN"

deseq_result$delabel <- NA

# (3) Check and set the threshold for p-value and fold change at which labels should be included

threshold <- head(arrange(deseq_result, pvalue), 10)$pvalue[10]

# (4) Add gene labels based on values that are lower than or equal to the threshold
 # Start with a query and end with what we want to perform should the query be fulfilled 

deseq_result$delabel[deseq_result$pvalue <= threshold & !is.na(deseq_result$pvalue)] <- (deseq_result$gene_id[deseq_result$pvalue <= threshold & !is.na(deseq_result$pvalue)])

# (5) Plot volcano

ggplot(data = deseq_result, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel)) + 
  geom_point() + theme_minimal() + geom_text_repel() + scale_colour_manual(values = c("green3", "grey", "red3")) +
  theme(text = element_text(size = 12)) + geom_vline(xintercept = c(-1, 1), col = "black", linetype = 2) + 
  geom_hline(yintercept = -log10(0.01), col = "black", linetype = 2)
