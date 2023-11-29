
# Functional analysis of differentially expressed genes in Eucalyptus grandis subject to control and droughted conditions
# Author: Rafael Keret

# REFERENCES

# (1) ClusterProfiler
 #  Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, Feng T, Zhou L, Tang W, Zhan L, Fu x, Liu S, Bo X, Yu G (2021). “clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.” The Innovation, 2(3), 100141. doi: 10.1016/j.xinn.2021.100141.

# (2) Arabidopsis thaliana database
 #  Carlson M (2019). org.At.tair.db: Genome wide annotation for Arabidopsis. R package version 3.8.2.

# LOAD PACKAGES

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("Rgraphviz")
BiocManager::install("enrichplot")

library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("dplyr")
library("tidyverse")
library("stringi")
library("pheatmap")

# IMPORT AND PREPARE DATA

# (1) Import DESeq2 output data

deseq_out <- read.csv("./Data/input/deseq_result.csv", head = TRUE)
colnames(deseq_out)[1]  <- "qseqid"

# (2) Rename gene symbols to ENTREZ ID 

deseq_out$qseqid <- stri_replace_all_regex(deseq_out$qseqid, 
                                           pattern=c("EucgrC_r007", "ISU1", "TUB2", "TUB5", "MIOX", "TUB4", "TUB3", "TUA1", 
                                                     "UGP", "NFU4", "EGM2", "WOX13.2", "UXS1", "EGM3"), 
                                           replacement=c("9845767", "104418091", "104432674", "104433897", "104436537", 
                                                         "104441567", "104448849", "104449214", "104449759", "104454102", 
                                                         "104444074", "104419574", "104420918", "104426035"), 
                                           vectorize=FALSE)

length(unique(deseq_out$qseqid))

# (3) Import BLASTX orthologs (i.e. Arabidopsis thaliana orthologs)

blast <- read.table("./Data/input/blastx_top_hits.txt", sep="\t", header = TRUE)
blast$qseqid <- as.character(blast$qseqid)

# (4) Left join by blast, to retain all ENTREZ ID's that have a corresponding A.thaliana ortholog
 # Only export "qseqid" and "log2FoldChange" that are relevant to MapMan software

deseq_orth <- left_join(blast, deseq_out, by = "qseqid") 
deseq_orth <- deseq_orth[which(!duplicated(deseq_orth$qseqid)), ]
length(unique(deseq_orth$qseqid))

write.csv(deseq_orth, file = "./Data/output/deseq_orth.csv", row.names = FALSE)
write.table(na.omit(deseq_orth[, c("qseqid", "log2FoldChange")]), 
            file = "./Data/output/deseq_orth_mapman.txt", row.names = FALSE, sep = "\t")

# (5) Create a filtered dataframe for investigation of xylogenesis related genes
 # Select differentially expressed genes with a "padj < 0.05", and a "|L2FC| > 1.0" for data mining

deseq_filtered <- deseq_orth %>% filter(deseq_orth$padj < 0.05)
deseq_filtered <- deseq_filtered %>% filter(abs(deseq_filtered$log2FoldChange) > 1)

# (6) Create gene ontology (GO) ID and terms table

GO <- read.table("./Data/input/GO_At.txt", sep="\t", header = TRUE, quote = "")
GO <- GO[, c(1,2,5)]

# (7) Left join to assign GO terms to E. grandis ENTREZ IDs

GO_assignment <- left_join(deseq_orth, GO, by = "sseqid", relationship = "many-to-many")

length(unique(GO_assignment$qseqid))

# CLUSTERPROFILER GSEA

# (1) Create term2gene and term2name tables

term2gene <- GO_assignment[, c(14, 1)] %>% select(ID:qseqid)

term2name <- GO_assignment[, c(14, 15)]

length(unique(term2gene$qseqid))

# (2) Create a ranked list for GSEA, ClusterProfiler uses Log2foldchange as a ranking metric
 # POSSILBE ALTERNATIVE RANKING METRIC: Deseq_orth <- mutate(Deseq_orth, rank = sign(log2FoldChange) * -log10(padj))

# (3) Extract Log2FoldChange values and corresponding query sequence ID ("qseqid")

all_genes <- deseq_orth$log2FoldChange

names(all_genes) <- deseq_orth$qseqid

# (4) Remove NAs

genes <- na.omit(all_genes)

# (5) Sort the ranked values in decreasing order

genes <- sort(genes, decreasing = TRUE)

### (3) Gene Set Enrichment Analysis

# (6) Gene Set Enrichment Analysis
 # Perform gene set enrichment analysis using GO terms
 # The software (warning message) recommended to NOT set the permutation level
 # Remember to set the seed to TRUE for consistent results as this is a randomization method

gseGO <- GSEA(genes, 
              minGSSize = 10, 
              maxGSSize = 480, 
              pvalueCutoff = 0.05, 
              pAdjustMethod = "BH", 
              TERM2GENE = term2gene, 
              TERM2NAME = term2name, 
              verbose = TRUE, 
              seed = TRUE, 
              eps = 0)

gse_DF <- as.data.frame(gseGO)
write.csv(gse_DF, file = "./Data/output/gsea_Eucalyptus.csv", row.names = FALSE)

# HISTOGRAM OF RELEVANT / NON-REDUNDANT ENRICHED CATEGORIES

# (1) Load packages

library("forcats")
require(DOSE)

# (2) Histogram of categories of interest
 # Choose biological processes (BP) from the gene set enrichment analysis (GSEA) list that are relevant to the research questions
 # Additionally remove potentially redundant categories

selected_ont <- c(257, 212, 161, 253, 145, 107, 185, 77, 98, 122, 197, 238, 12, 24, 22, 8, 171)

gseGO_histogram <- data.frame(Description = gseGO$Description[selected_ont],
                              ID = gseGO$ID[selected_ont], 
                              NES = gseGO$NES[selected_ont], 
                              qvalue = gseGO$qvalue[selected_ont])

# (3) Function to assign a generalized process to a BP based on partial matches

assign_process <- function(x) {
  query_terms <- c("GO:0009744", "GO:0009746", "GO:0009750", "GO:0009804", "GO:0006749", "GO:0010038",
                   "GO:0006012", "GO:0031407", "GO:0010876", "GO:0019915", "GO:0046503", "GO:0006639",
                   "GO:0042542", "GO:0071456", "GO:0036294", "GO:0000302", "GO:0009813")
  set_terms <- c("osmotic", "osmotic", "osmotic", "osmotic", "osmotic", "osmotic", 
                 "Nucleation", "Nucleation", "Nucleation", "Nucleation", "Nucleation", "Nucleation",
                 "Stress", "Stress", "Stress", "Stress", "Stress")
  
  for (i in seq_along(query_terms)) {
    if (grepl(query_terms[i], x, fixed = TRUE)) {
      return(set_terms[i])
    }
  }
  return(NA)  # If no partial match is found, return NA or any other default value
}

gseGO_histogram$process <- sapply(gseGO_histogram$ID, assign_process)

# (4) Arrange based on process and NES column, then rename "qvalue"

gseGO_histogram <- arrange(gseGO_histogram, process, NES)
gseGO_histogram <- rename(gseGO_histogram, p.adj = qvalue)

# (5) Select order of bars, and plot

desired_order <- gseGO_histogram$Description

histo <- ggplot(gseGO_histogram, aes(NES, factor(Description, levels = desired_order), fill = p.adj)) +
  geom_col() +
  scale_fill_gradientn(colours = c("#b3eebe", "#46bac2", "#371ea3"),
                       guide = guide_colorbar(reverse = TRUE)) +
  theme_minimal() +
  ylab(NULL) +
  xlab("NES") + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 90), expand = c(0.005, 0)) + 
  theme(axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 16), 
        axis.text.x = element_text(size = 16)) +
  theme(plot.margin = unit(c(1, 0, 3.0, 3.0), "lines")) + 
  theme( legend.text = element_text(size = 16), legend.title = element_text(size = 16), 
         legend.key.size = unit(1, "lines"))

# CATEGORY NET PLOT

# (1) Category net plot, Including CATEGORY labels

set.seed(4)
cnet <- cnetplot(gseGO, categorySize = "pvalue", color.params = list(foldChange = genes), 
               showCategory = gseGO$Description[c(257, 212, 161, 253, 145, 107, 185, 
                                                  77, 98, 122, 197, 238, 12, 24, 
                                                  22, 8, 171, 94)],
               layout = 'fr',
               node_label = "category", 
               cex.params = list(category_node = 0.8, gene_node = 0.2, category_label = 0.85), 
               color_category = "firebrick")

# (2) Join histogram and Category net plot

cowplot::plot_grid(histo, cnet, ncol = 1, labels = LETTERS[1:2], rel_widths = c(.8, .8, 1.2), label_size = 22)

# EXTRACT AND LABEL DIFFERENTIALLY EXPRESSED GENES

# (1) Read in gene symbol and functions table

gene_symbols_function <- read.csv("./Data/input/gene_symbols_function.csv", head = TRUE)

# (2) Left join to add gene names and descriptions to filtered E. grandis ENTREZ ID's

gene_symbol <- left_join(deseq_filtered, gene_symbols_function, by = "sseqid", relationship = "many-to-many")

gene_symbol <- select(gene_symbol, c(1, 2, 9, 10, 13, 15, 16))

gene_symbol <- distinct(gene_symbol)

length(unique(gene_symbol$qseqid))

write.csv(gene_symbol, "./Data/output/gene_symbol.csv", row.names = FALSE)

# QUERY FILTERING OF DEGs 

# (1) Import MapMan classification list to help identify genes relevant to embolism prevention and recovery

mapman_classifications <- read.delim("./Data/input/MapMan_classifications.txt", header = TRUE)
mapman_classifications <- rename(mapman_classifications, qseqid = id, L2FC = deseq_orth_mapman.txt)
mapman_classifications$qseqid <- as.character(mapman_classifications$qseqid)
mapman_classifications <- left_join(deseq_filtered, mapman_classifications, by = "qseqid") %>% 
  select(BinCode, BinName, qseqid, type, description, L2FC)

# (2) Identify genes relevant to water transport or cavitation prevention
 # Aquaporins / Lipids

prevention_ID <- c("104446186", "104423212", "104432876", "104435641", "104432875", 
                   "104454378", "104420926", "104422516", "104414444", "104449893", 
                   "104441673", "104450070", "104426097", "104436271", "104450169", 
                   "104456523", "104446880", "104446890", "104422632", "104436828", 
                   "104432538", "104443704", "104414777", "104430096")

# (3) Extract relevant ENTREZ ID's from differentially expressed gene set

prevention_query <- gene_symbol %>% 
  filter(grepl(paste(prevention_ID, collapse = "|"), qseqid, ignore.case = TRUE),
         abs(log2FoldChange) >= 1.0, padj <= 0.05)

# (4) Replace TAIR ID, with gene symbol, and create gene id/symbol column

prevention_query$external_gene_name <- stri_replace_all_regex(prevention_query$external_gene_name, 
                                                           pattern=c("AT1G12100", "AT2G10940", "AT2G25890", 
                                                                     "AT2G37870", "AT3G01570", "AT3G16175", 
                                                                     "AT3G18570", "AT3G22600", "AT3G23470"), 
                                                           replacement=c("LTP", "LTP", "OLE", "LTP4", "OLE", 
                                                                         "TE", "OLE", "LTP19", "CFA"), 
                                                           vectorize=FALSE)

prevention_query$ID_symbol <- paste(prevention_query$qseqid, "(", prevention_query$external_gene_name, ")", sep = "")
prevention_query <- rename(prevention_query, p.adj = padj)

# (5) Arrange based on prevention_ID

prevention_query <- prevention_query[order(match(prevention_query$qseqid, prevention_ID)), ]

# (6) Plot data

prevention_order <- prevention_query$ID_symbol

pre <- ggplot(prevention_query, 
             aes(x = factor(ID_symbol, levels = prevention_order), y = log2FoldChange, fill = p.adj)) +
  geom_bar(stat = "identity", orientation = "x") +
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3")) +
  labs(x = NULL,
       y = expression("Log"[2]~"fold change"), 
       fill = "p.adj") + theme_minimal() +
  geom_errorbar(aes(y = log2FoldChange, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), 
                width = 0.4, colour = "black", alpha = 0.9, size = 0.02) +
  theme(panel.grid.major.x = element_line(color = "grey90", size = 0.5),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines")) + 
  theme(axis.title.y = element_text(size = 18, margin = margin(r = 15)), 
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18, margin = margin(r = 15)), 
        axis.text.x = element_text(size = 18, angle = 90, vjust = 0.3, hjust = 1)) + 
  theme( legend.text = element_text(size = 12), legend.title = element_text(size = 12), 
         legend.key.size = unit(1.0, "lines")) + 
  scale_y_continuous(breaks = seq(-6, 6, by = 2))

# OSMOTIC REGULATION

# (1) Identify genes relevant to osmotic regulation

osmotic_ID <- c("104440059", "104437966", "104440058", "104437963", "104436385", "104453398", "104434467", "104449979",
                "104448846", "104454126", "104422020", "104415036", "104449331", "104456017", "104415812", "104432350",
                "104441268", "104419047", "104456866", "104444849", "104414952", "104416472", "104448611", "104426131", 
                "104449640", "104441763", "104414892", "104414895", "104448272", "104439566", "104419638", "104433084", 
                "104448775", "104457372", "104425170", "104452949", "104437424", "104456105", "104450952")

# (2) Extract relevant ENTREZ ID's from differentially expressed gene set

osmotic_query <- gene_symbol %>% 
  filter(grepl(paste(osmotic_ID, collapse = "|"), qseqid, ignore.case = TRUE),
         abs(log2FoldChange) >= 1.0, padj <= 0.05)

# (3) Replace TAIR ID, with gene symbol, and create gene id/symbol column

osmotic_query$external_gene_name <- stri_replace_all_regex(osmotic_query$external_gene_name, 
                                                pattern=c("AT1G06030", "AT1G54730", "AT2G36950", 
                                                          "AT5G18840", "AT1G01490", "AT3G06130", 
                                                          "AT5G23760", "AT5G27690", "ESL1", 
                                                          "ATBFRUCT1"), 
                                                replacement=c("FRK2/6", "ERD6-like5", "ATHMP20", 
                                                              "ERD6-like16", "ATHMP1", "ATHMP25", 
                                                              "CTP", "ATHMP49", "ERD6-like1", 
                                                              "CWINV1"), 
                                                vectorize=FALSE)

osmotic_query$ID_symbol <- paste(osmotic_query$qseqid, "(", osmotic_query$external_gene_name, ")", sep = "")
osmotic_query <- rename(osmotic_query, p.adj = padj)

# (4) Arrange based on osmotic_ID

osmotic_query <- osmotic_query[order(match(osmotic_query$qseqid, osmotic_ID)), ]

# (5) Plot data

osmotic_order <- osmotic_query$ID_symbol

osm <- ggplot(osmotic_query, 
              aes(x = factor(ID_symbol, levels = osmotic_order), y = log2FoldChange, fill = p.adj)) +
  geom_bar(stat = "identity", orientation = "x") +
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3")) +
  labs(x = NULL,
       y = expression("Log"[2]~"fold change"), 
       fill = "p.adj") + theme_minimal() +
  geom_errorbar(aes(y = log2FoldChange, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), 
                width = 0.4, colour = "black", alpha = 0.9, size = 0.02) +
  theme(panel.grid.major.x = element_line(color = "grey90", size = 0.5),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines")) + 
  theme(axis.title.y = element_text(size = 18, margin = margin(r = 15)), 
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18, margin = margin(r = 15)), 
        axis.text.x = element_text(size = 18, angle = 90, vjust = 0.3, hjust = 1)) + 
  theme( legend.text = element_text(size = 12), legend.title = element_text(size = 12), 
         legend.key.size = unit(1.0, "lines")) + 
  scale_y_continuous(breaks = seq(-6, 8, by = 2))

# (6) Combine and export dataframes

prevention_osmotic <- rbind(prevention_query, osmotic_query)

write.csv(prevention_osmotic, "./Data/output/prevention_osmotic_genes.csv", row.names = FALSE)
