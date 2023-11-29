
# Creating a custom mapping file for MapMan
# Author: Rafael Keret

# LOAD PACKAGES

library("dplyr")
library("tidyverse")
library("stringi")

# Building custom MapMan mapping file

# (1) Read in the MapMan arabidopsis thaliana mapping file

mapping_file <- read.csv("./Data/input/Ath_AGI_ISOFORM_MODEL_TAIR10_Aug2012.csv", skip = 1)
colnames(mapping_file)[3] <- "sseqid"
mapping_file$sseqid <- sub("\\..*", "", sub("'", "", mapping_file$sseqid))

# (2) Read in the deseq_orth file

deseq_orth <- read.csv("./Data/input/deseq_orth.csv", head = TRUE)
deseq_orth$sseqid <- tolower(deseq_orth$sseqid)

# (3) Merge deseq_orth to mapping_file to create a custom mapping file for E. grandis based on ENTREZ ID's
 # For example, the TAIR ID will be replaced by its corresponding E. grandis ENTREZ ID

# (4) Merge dataframes but keep all rows for both

custom_mapping_E_grandis <- merge(deseq_orth[, c("qseqid", "sseqid")], mapping_file, by = "sseqid", all.y = TRUE)
custom_mapping_E_grandis$qseqid <- as.character(custom_mapping_E_grandis$qseqid)

# (5) Fill empty rows in the qseqid column with TAIR ID's 
 # Not all TAIR ID's had a corresponding ENTREZ ID, hence the blanks can be refilled with the original TAIR ID's

custom_mapping_E_grandis <- custom_mapping_E_grandis %>% mutate(qseqid = coalesce(qseqid, sseqid))

# (6) The original file had two single quotation marks ('') for the blank spaces and surrounding the gene ID's ('LOC123')

custom_mapping_E_grandis$qseqid <- ifelse(grepl("'", custom_mapping_E_grandis$qseqid), 
                                          paste0("'", custom_mapping_E_grandis$qseqid), 
                                          paste0("'", custom_mapping_E_grandis$qseqid, "'"))

# (7) Reorder custom mapping file

custom_mapping_E_grandis <- custom_mapping_E_grandis[, -1]  
colnames(custom_mapping_E_grandis)[1] <- "IDENTIFIER"  
custom_mapping_E_grandis <- custom_mapping_E_grandis[, c(2, 3, 1, 4, 5)]
custom_mapping_E_grandis <- custom_mapping_E_grandis[order(custom_mapping_E_grandis[, 1]), ]

# (8) Removing duplicate ID's
 # A single TAIR ID can have multiple matching ENTREZ ID's
 # But since the ENTREZ ID's are unique, these wont be captured as duplicates in the mapping file

check_duplicates <- custom_mapping_E_grandis[duplicated(custom_mapping_E_grandis$IDENTIFIER) 
                                             | duplicated(custom_mapping_E_grandis$IDENTIFIER, 
                                                          fromLast=TRUE), ]

custom_mapping_E_grandis <- distinct(custom_mapping_E_grandis)

length(unique(custom_mapping_E_grandis$IDENTIFIER))

# (9) Write new custom mapping file as a tab delimited text file

write.table(custom_mapping_E_grandis, file = "./Data/output/custom_mapping_E_grandis.txt", sep = "\t", na = "", quote = FALSE, row.names = FALSE)
