# Transcriptomic analysis pipeline

Experiments were conducted to determine how drought induced transcriptomic remodeling influences the xylem anatomy and physiology of _Eucalyptus grandis_. During a prolonged drought period, total RNA was harvested from _E. grandis_ stem samples, comprising mainly of developing xylem and cambium tissues. Truseq stranded mRNA libraries were created and subject to paired-end next generation RNA sequencing using the Novaseq 6000 Illumina platform (NCBI: ). The resulting sequencing libraries were processed/analysed with FastQC v0.11.9, Trimmomatic v0.3.2, Hisat2 v2.2.1 and FeatureCounts v2.0.5 (Subread). Next, the gene counts were subject to differential expression analysis using the DESeq2 package in R (v4.3.1). Lastly, functional annotations and gene ontologies were assigned to the _E. grandis_ ENTREZ IDs prior to gene set enrichment analysis using clusterProfiler v4.0 within R.

This repository contains the following folders: 

1. DESeq2 - Annotated R code and relevant data to perform differential expression analysis on the _E. grandis_ gene set.
2. Gene_annotation - Relevant R code required to annotate the _E. grandis_ ENTREZ IDs with gene ontologies and symbols using _Arabidopsis thaliana_ orthologs.
3. Functional_characterization - Performing gene set enrichment analysis and identifying expression candidates responsible for controlling the wood anatomical traits in _E. grandis_.
4. MapMan - Creating a custom MapMan mapping file using _E. grandis_ ENTREZ IDs.

All folders represent R projects that contain all the necessary code, and datasets required to repeat the analysis. 
