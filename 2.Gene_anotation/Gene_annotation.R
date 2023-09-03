
# Annotating Eucalyptus grandis ENTREZ IDs with GO, KEGG terms and gene symbols, using A. thaliana orthologs identified by blastX.
# Author: Rafael Keret

# INSTALL AND LOAD PACKAGES
  
BiocManager::install("AnnotationDbi")
BiocManager::install("org.At.tair.db")
BiocManager::install("GO.db")
BiocManager::install("clusterProfiler")
BiocManager::install("limma")

library("AnnotationDbi")
library("org.At.tair.db")
library("GO.db")
library("dplyr")
library("stringi")
library("clusterProfiler")
library("splitstackshape")

# IMPORT AND PREPARE DATA

# (1) BLASTX results table

blast <- read.table("./Data/input/blastx_orthologs_ncbi.txt", sep="\t", header = FALSE)
colnames(blast)  <- c("qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore")
blast$sseqid <- gsub("\\..*", "", blast$sseqid)

# (2) Select the top A.thaliana TAIR ID hit for each unique E.grandis ENTREZ gene ID 

blast <- blast %>% 
  group_by(qseqid) %>%
  slice_min(evalue) %>%
  slice_max(bitscore) %>%
  slice_max(length) %>%
  slice_min(mismatch) %>%
  ungroup

length(unique(blast$qseqid))

# (3) Retain unique E. grandis ENTREZ IDs with a corresponding TAIR ortholog

blast <- blast[which(!duplicated(blast$qseqid)), ]

# (4) Rename chromosomal regions to ENTREZ ID

blast$qseqid <- stri_replace_all_regex(blast$qseqid, 
                                       pattern=c("NC_014570.1:c139556-136747", "NC_052613.1:3272908-3275532", 
                                                 "NC_052613.1:38624984-38627143", "NC_052613.1:c55328170-55325438", 
                                                 "NC_052614.1:11229172-11231700", "NC_052616.1:c28054277-28051374", 
                                                 "NC_052617.1:c20483764-20479776", "NC_052617.1:c34031899-34029393", 
                                                 "NC_052617.1:c40716758-40710803", "NC_052618.1:51481974-51505466", 
                                                 "NC_052620.1:c30982539-30978927", "NC_052620.1:c35825863-35822327", 
                                                 "NC_052621.1:753837-759217", "NC_052622.1:32327277-32331798"), 
                                       replacement=c("9845767", "104418091", "104432674", "104433897", "104436537", 
                                                     "104441567", "104448849", "104449214", "104449759", "104454102", 
                                                     "104444074", "104419574", "104420918", "104426035"), 
                                       vectorize=FALSE)

# (5) Export final BLASTX results table

write.table(blast, "./Data/output/blastx_top_hits.txt", sep="\t", row.names = FALSE, 
            quote = FALSE, col.names = TRUE)

# CLUSTER PROFILER GO ANNOTATION

# (1) Select levels 1 to 10 for GO term searches, and only select biological processes ("BP")

ggo <- groupGO(gene     = blast$sseqid,
               OrgDb    = org.At.tair.db,
               ont      = "BP",
               level    = c(1:10),
               readable = FALSE, 
               keyType = "TAIR")

ggo_DF <- as.data.frame(ggo) 

# According to the "ggo_DF" table ("Biological process") 13078 / 13577 TAIR IDs have been annotated

# (2) Split the individual TAIR IDs in the GeneID column into individual rows

GO_At <- concat.split.multiple(ggo_DF, "geneID", seps="/", "long")
colnames(GO_At)[5] <- "sseqid"  
GO_At <- distinct(GO_At)

write.table(GO_At, "./Data/output/GO_At.txt", sep="\t", row.names = FALSE, 
            quote = FALSE, col.names = TRUE)

# KEGG annotation

library("limma")

# (1) Assign KEGG pathway IDs to A.thaliana TAIR IDs

path_id <- bitr_kegg(blast$sseqid, fromType = "kegg", toType = "Path", organism = "ath")

colnames(path_id)  <- c("keggID", "PathwayID")

# (2) Add descriptions to the KEGG pathway IDs

path_name <- getKEGGPathwayNames(species = "ath")

kegg_At <- left_join(path_id, path_name, "PathwayID")

kegg_At$Description <- gsub("-.*", "", kegg_At$Description)
kegg_At <- distinct(kegg_At)
colnames(kegg_At)[1]  <- "sseqid"

write.table(kegg_At, "./Data/output/kegg_At.txt", sep="\t", row.names = FALSE, 
            quote = FALSE, col.names = TRUE)

# GENE SYMBOL, DESCRIPTION AND ONTOLOGY

# (1) Load biomaRt package then retrieve A.thaliana gene symbol (TAIR) and function

library("biomaRt")

tair_mart <- useMart(biomart = "plants_mart",
                     host = "plants.ensembl.org", dataset = "athaliana_eg_gene")

At_symbols <- getBM( values = blast$sseqid, 
                     mart = tair_mart, 
                     attributes = c("ensembl_gene_id", "entrezgene_id", 
                                    "description", "external_gene_name"), 
                     filters = "ensembl_gene_id")

colnames(At_symbols)[c(1, 3)] <- c("sseqid", "function")

# (2) Export gene symbol and function file

write.csv(At_symbols, "./Data/output/gene_symbols_function.csv", row.names = FALSE)
