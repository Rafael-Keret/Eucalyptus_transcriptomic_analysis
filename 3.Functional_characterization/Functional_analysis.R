
# Functional analysis of the differentially expressed genes occurring in E. grandis xylem during drought.
# Author: Rafael Keret

# REFERENCES

# (1) ClusterProfiler
 # Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, Feng T, Zhou L, Tang W, Zhan L, Fu x, Liu S, Bo X, Yu G (2021). “clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.” The Innovation, 2(3), 100141. doi: 10.1016/j.xinn.2021.100141.

# (2) Arabidopsis thaliana database
 # Carlson M (2019). org.At.tair.db: Genome wide annotation for Arabidopsis. R package version 3.8.2.

# INSTALL AND LOAD PACKAGES

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

# (6) Gene Set Enrichment Analysis
 # Perform gene set enrichment analysis using GO terms
 # The software (warning message) recommended to NOT set the permutation level
 # Remember to set the seed to TRUE for consistent results as this is a randomization method

gseGO <- GSEA(genes, 
             minGSSize = 50, 
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

# HISTOGRAM OF RELEVANT/NON-REDUNDANT ENRICHED CATEGORIES

# (1) Load packages

library("forcats")
require(DOSE)

# (2) Histogram of categories of interest
 # Choose biological processes (BP) from the gene set enrichment analysis (GSEA) list that are relevant to the research questions
 # Additionally remove potentially redundant categories

selected_ont <- c(31, 107, 43, 20, 19, 133, 191, 184, 192, 168, 139, 
                  126, 134, 199, 200, 37, 177, 25, 153, 81, 79, 76, 9, 1, 
                  175, 12, 173)

gseGO_histogram <- data.frame(Description = gseGO$Description[selected_ont],
                              ID = gseGO$ID[selected_ont], 
                              NES = gseGO$NES[selected_ont], 
                              qvalue = gseGO$qvalue[selected_ont])

# (3) Function to assign a generalized process to a BP based on partial matches

assign_process <- function(x) {
  query_terms <- c("GO:0006260", "GO:0044786", "GO:0009826", "GO:0007018", "GO:0000910", 
                   "GO:0016051", "GO:1901659", "GO:1901361", "GO:0072329", "GO:0015980", 
                   "GO:0009698", "GO:0009813", "GO:0016102", "GO:0016114", "GO:0016144", 
                   "GO:0045491", "GO:0010411", "GO:0010410", "GO:0009825", "GO:0070592", 
                   "GO:0042545", "GO:0009664", "GO:0000302", "GO:0009408", "GO:0034050", 
                   "GO:0042542", "GO:0009626")
  set_terms <- c("Development", "Development", "Development", "Development", "Development", 
                 "Primary", "Primary", "Primary", "Primary", "Primary", 
                 "Secondary", "Secondary", "Secondary", "Secondary", "Secondary", 
                 "Cell wall", "Cell wall", "Cell wall", "Cell wall", "Cell wall", 
                 "Cell wall", "Cell wall", "Stress", "Stress", "Stress", "Stress", 
                 "Stress")
  
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
  labs(x = "NES", y = NULL) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 90), expand = c(0.005, 0)) + 
  theme(axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 15), 
        axis.text.x = element_text(size = 15)) +
  theme(plot.margin = unit(c(1, 0, 3.0, 3.0), "lines")) + 
  theme( legend.text = element_text(size = 15), legend.title = element_text(size = 15), 
         legend.key.size = unit(1, "lines"))

# ENRICHMENT NETWORK MAP OF SIGNIFICANTLY ENRICHED ONTOLOGIES

# (1) Enrichment map showing category interactions, functional association or shared genes
 # To remove redundant GO terms
 # node_label = "category" or "none" 

gseGO_2 <- pairwise_termsim(gseGO)

set.seed(24)
emapp <- emapplot(gseGO_2, showCategory = gseGO_2$Description[c(1:200)], color = "NES", 
         cex.params = list(category_node = 0.5, category_label = 0.1, 
                           line = 0.3, label_group = 0.5), node_label = "category", 
         cluster.params = list(cluster = FALSE, legend = TRUE, n = 2, 
                               label_style = "shadowtext", label_words_n = 1, 
                               label_format = 5, repel = TRUE, group_legend = TRUE, 
                               method = stats::kmeans), 
         edge.params = list(show = TRUE, min = 0.08))

# (2) Join histogram and network plot

cowplot::plot_grid(emapp, histo, ncol = 1, labels = LETTERS[1:2], rel_widths = c(.8, .8, 1.2), label_size = 22)

# EXTRACT AND LABEL DIFFERENTIALLY EXPRESSED GENES

# (1) Read in gene symbol and functions table

gene_symbols_function <- read.csv("./Data/input/gene_symbols_function.csv", head = TRUE)

# (2) Left join to add gene names and descriptions to filtered E. grandis ENTREZ IDs

gene_symbol <- left_join(deseq_filtered, gene_symbols_function, by = "sseqid", relationship = "many-to-many")

gene_symbol <- select(gene_symbol, c(1, 2, 9, 10, 13, 15, 16))

gene_symbol <- distinct(gene_symbol)

length(unique(gene_symbol$qseqid))

write.csv(gene_symbol, "./Data/output/gene_symbol.csv", row.names = FALSE)

# QUERY FILTERING OF DEGs 

# (1) Import MapMan classification list to help identify genes relevant to xylogenesis

mapman_classifications <- read.delim("./Data/input/MapMan_classifications.txt", header = TRUE)
mapman_classifications <- rename(mapman_classifications, qseqid = id, L2FC = deseq_orth_mapman.txt)
mapman_classifications$qseqid <- as.character(mapman_classifications$qseqid)
mapman_classifications <- left_join(deseq_filtered, mapman_classifications, by = "qseqid") %>% 
  select(BinCode, BinName, qseqid, type, description, L2FC)

# (2) Identify genes relevant to wood formation or xylogenesis
 # Genes related to cellular expansion and cell wall component biosynthesis

xylem_development_ID <- c("104433767", "104433819", "104433834", "104433865", "104433885", "104454356", "104442333", "104436371", 
                          "104444630", "104444640", "104444649", "104449794", "104452753", "104452752", "104421716", "104421040", 
                          "104427839", "104421714", "104450069", "104448602", "104434058", "104433548", "104443623", "104443626", 
                          "104428622", "104432727", "104432728", "104424356", "104440602", "104441874", "104419526", "104443687", 
                          "104414336", "104456622", "104453686", "104426781", "104424340", "104424341", "104456165", "108953967", 
                          "104454145", "104425949", "104417282", "104414384", "104454440", "104430693", "104419915", "104420951", 
                          "104449520", "104446186", "104423212", "104449672", "104420130", "104436967", "104450842", "104450844", 
                          "104454116", "104442048", "104441676", "104441677", "104420942", "104415214", "104432862", "104421673", 
                          "104441055", "104416110", "104425811", "104449510", "104435787", "104457357", "104415154", "104450362")

# (3) Extract relevant ENTREZ IDs from differentially expressed gene set

xylogenesis_query <- gene_symbol %>% 
  filter(grepl(paste(xylem_development_ID, collapse = "|"), qseqid, ignore.case = TRUE),
         abs(log2FoldChange) >= 1.0, padj <= 0.05)

# (4) Replace TAIR ID, with gene symbol, and create gene id/symbol column

xylogenesis_query$external_gene_name <- stri_replace_all_regex(xylogenesis_query$external_gene_name, 
                                                  pattern=c("AT1G02460", "AT1G04680", "AT1G14890", "AT1G55770", "AT2G38150", "AT3G18180", 
                                                            "AT3G24130", "AT3G50990", "AT3G53190", "AT3G55700", "AT3G57380", "AT3G59850", 
                                                            "AT4G19420", "AT4G24780", "AT4G30380", "AT5G19730", "AT5G39580", "AT5G62360", 
                                                            "AT1G62660"), 
                                                  replacement=c("PLL", "PLL26", "PMEI", "PMEI", "1,4-GTF", "GT61", "PME29", "PRX36", "PLL17", 
                                                                "UGT76F1", "GT61", "PLL", "PAE", "PLL19", "EXLB2", "PME53", "PRX62", "PMEI13", 
                                                                "VI1"), 
                                                  vectorize=FALSE)

xylogenesis_query$ID_symbol <- paste(xylogenesis_query$qseqid, "(", xylogenesis_query$external_gene_name, ")", sep="")

# (5) Arrange based on xylem_development_ID

xylogenesis_query <- xylogenesis_query[order(match(xylogenesis_query$qseqid, xylem_development_ID)), ]
xylogenesis_query <- rename(xylogenesis_query, p.adj = padj)

# (6) Break data up for 2 plots

xylogenesis_query_1 <- xylogenesis_query[1:36, ]
xylogenesis_query_2 <- xylogenesis_query[37:72, ]

# (7) Invert xylogenesis_query_2 for plotting purposes

xylogenesis_query_2 <- xylogenesis_query_2[rev(seq_len(nrow(xylogenesis_query_2))), ]

# (8) Plot data

xylogeneis_order <- xylogenesis_query$ID_symbol

x1 <- ggplot(xylogenesis_query_1, 
             aes(x = log2FoldChange, y = factor(ID_symbol, levels = xylogeneis_order), fill = p.adj)) +
  geom_bar(stat = "identity", orientation = "y") +
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3")) +
  labs(x = expression("Log"[2]~"fold change"), y = NULL, fill = "p.adj") + theme_minimal() +
  geom_errorbar(aes(x = log2FoldChange, xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), 
                width = 0.4, colour = "black", alpha = 0.9, size = 0.02) +
  theme(panel.grid.major.x = element_line(color = "grey90", size = 0.5),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5)) +
  theme(plot.margin = unit(c(2, 0, 0, 1.8), "lines")) + 
  theme(axis.text.y = element_text(size = 18), axis.title.x = element_text(size = 18), 
        axis.text.x = element_text(size = 18)) + 
  theme( legend.text = element_text(size = 12), legend.title = element_text(size = 12), 
         legend.key.size = unit(1.0, "lines")) +
  scale_y_discrete(limits=rev) + scale_x_continuous(breaks = seq(-6, 6, by = 3))

x2 <- ggplot(xylogenesis_query_2, 
             aes(x = log2FoldChange, y = factor(ID_symbol, levels = xylogeneis_order), fill = p.adj)) +
  geom_bar(stat = "identity", orientation = "y") +
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3")) +
  labs(x = expression("Log"[2]~"fold change"), y = NULL, fill = "p.adj") + theme_minimal() + 
  geom_errorbar(aes(x = log2FoldChange, xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), 
                width = 0.4, colour = "black", alpha = 0.9, size = 0.02) +
  theme(panel.grid.major.x = element_line(color = "grey90", size = 0.5),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5)) +
  theme(plot.margin = unit(c(2, 0, 0, 1.8), "lines")) + 
  theme(axis.text.y = element_text(size = 18), axis.title.x = element_text(size = 18), 
        axis.text.x = element_text(size = 18)) + 
  theme( legend.text = element_text(size = 12), legend.title = element_text(size = 12), 
         legend.key.size = unit(1.0, "lines")) +
  scale_y_discrete(limits=rev) + scale_x_continuous(breaks = seq(-6, 6, by = 3))

cowplot::plot_grid(x1, x2, ncol = 2, labels = LETTERS[1:2], rel_widths = c(1, 1, 1.2), label_size = 26)

# TRANSCRIPTION FACTORS GOVERNING XYLEM DEVELOPMENT

# (1) Enter TF ENTREZ IDs relevant to wood formation

TF_ID <- c("104419933", "104448414", "104418025", "104450926", "104452408", "104436535", 
           "104424951", "104433270", "104436346", "104442190", "104450012", "104453968", 
           "104448388", "104426573", "104414409", "104415262", "104434101", "104426233", 
           "104437240", "104438144", "104418936", "104453423", "104443670", "104424339", 
           "104442593", "104417786", "104433840", "104425859", "104424612", "104431568", 
           "104432590", "104454661", "104449656", "104455407", "104419371")

# (2) Extract relevant ENTREZ IDs from differentially expressed gene set

TF <- gene_symbol %>% 
  filter(grepl(paste(TF_ID, collapse = "|"), qseqid, ignore.case = TRUE),
         abs(log2FoldChange) >= 1.0, padj <= 0.05)

# (3) Replace TAIR ID, with gene symbol, and create gene id/symbol column

TF$external_gene_name <- stri_replace_all_regex(TF$external_gene_name, 
                                                               pattern=c("MYB102", "MYB305", "MYB57", "ATAF1"), 
                                                               replacement=c("MYB41", "MYB13", "MYB4", "NAC2"), 
                                                               vectorize=FALSE)

TF$ID_symbol <- paste(TF$qseqid, "(", TF$external_gene_name, ";", round(TF$log2FoldChange, 2), ")", sep = " ")

# (4) Order TFs based on common function

TF <- subset(TF, grepl(paste(TF_ID, collapse = "|"), 
                       qseqid))[match(TF_ID, 
                                      subset(TF, grepl(paste(TF_ID, collapse = "|"), 
                                                       qseqid))$qseqid), ]

# (5) Create matrix for heatmap

TF_matrix <- matrix(TF$log2FoldChange, nrow = nrow(TF))

# (6) Set the row and column names for the heatmap

rownames(TF_matrix) <- TF$ID_symbol
colnames(TF_matrix) <- "MYB and NAC TFs"

# (7) Specify breaks for colour palette, to set zero as "white"

breaks_TF <- c(seq(min(TF_matrix), -1e-10, length.out = 50), 0, seq(1e-10, max(TF_matrix), length.out = 50))

pheatmap(TF_matrix, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_row = 16, fontsize_col = 11, 
         breaks = breaks_TF, color = colorRampPalette(c("blue", "white", "red"))(length(breaks_TF) - 1), 
         cellwidth = 11, cellheight = 11, show_colnames = FALSE, gaps_row = c(1:35))

# PHENYLPROPANOID PATHWAY PEROXIDASES

# (1) Identify peroxidases relevant to lignin polymerisation

Per_ID <- c("104450831", "104415764", "104434128", "104434134", "104421005", 
                "104425008", "104422800", "104450451", "104449271")

# (2) Extract relevant ENTREZ IDs from differentially expressed gene set

Per <- gene_symbol %>% 
  filter(grepl(paste(Per_ID, collapse = "|"), qseqid, ignore.case = TRUE),
         abs(log2FoldChange) >= 1.0, padj <= 0.05)

# (3) Replace TAIR ID, with gene symbol

Per$external_gene_name <- stri_replace_all_regex(Per$external_gene_name, 
                                                 pattern=c("AT5G51890", "AT2G22420", "AT1G71695", "PER64"), 
                                                 replacement=c("PRX66", "PRX17", "PRX12", "PRX64"), 
                                                 vectorize=FALSE)

# (4) Arrange based on Per_ID

Per <- Per[order(match(Per$qseqid, Per_ID)), ]
Per <- rename(Per, p.adj = padj)

# (5) Merge Gene ID and Symbol

Per$ID_symbol <- paste(Per$qseqid, "(", Per$external_gene_name, ")", sep="")

# (6) Plot data

Per_order <- Per$ID_symbol

ggplot(Per, aes(x = log2FoldChange, y = factor(ID_symbol, levels = Per_order), fill = p.adj)) +
  geom_bar(stat = "identity", orientation = "y") +
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3")) +
  labs(x = expression("Log"[2]~"fold change"),
       y = "Peroxidase", 
       fill = "p.adj") + theme_minimal() + 
  geom_errorbar(aes(x = log2FoldChange, xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), 
                width = 0.4, colour = "black", alpha = 0.9, size = 0.02) +
  theme(panel.grid.major.x = element_line(color = "grey90", size = 0.5),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines")) + 
  theme(axis.title.y = element_text(size = 22), axis.text.y = element_text(size = 22), 
        axis.title.x = element_text(size = 22), axis.text.x = element_text(size = 22)) + 
  theme( legend.text = element_text(size = 18), legend.title = element_text(size = 12), 
         legend.key.size = unit(1.5, "lines")) + 
  scale_x_continuous(breaks = seq(-2, 6, by = 2))
