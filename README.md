# Analysis of the mRNA transcripts expressed in the xylem of _Eucalyptus grandis_ under droughted and control conditions.

Experiments were conducted to determine how drought induced transcriptomic remodelling influences the xylem anatomy and physiology of _Eucalyptus grandis_. After a prolonged drought period, total RNA was extracted from _E. grandis_ stem samples, comprising mainly of developing xylem and cambium tissues. TruSeq Stranded mRNA libraries were created and subject to paired-end Next Generation Sequencing (i.e. RNAseq) using the Novaseq 6000 Illumina platform (NCBI: PRJNA1012834). The resulting sequencing libraries were processed/analysed with FastQC v0.11.9, Trimmomatic v0.3.2, Hisat2 v2.2.1 and FeatureCounts v2.0.5 (Subread). Next, the read counts were subject to differential expression analysis using the DESeq2 package in R (v4.3.1). The closest _Arabidopsis thaliana_ orthologs, identified by blastx, were used to assign gene ontology, functional descriptions and symbols to the _E. grandis_ ENTREZ IDs using the clusterProfiler v4.0 and biomaRt packages within R. Finally, gene set enrichment analysis (GSEA) was performed with clusterProfilers non-model organism feature, using R. 

This repository contains the following folders: 

1. DESeq2 - Annotated R code and relevant RNAseq counts data to perform differential expression analysis on _E. grandis_ subject to droughted versus well-watered (control) conditions.
2. Gene_annotation - Relevant R code required to annotate the _E. grandis_ ENTREZ IDs with gene ontology, functional descriptions and symbols using _Arabidopsis thaliana_ orthologs.
3. Functional_characterization_structural - Performing gene set enrichment analysis and identifying the gene expression candidates responsible for controlling the developmental or structural xylem anatomical traits in _E. grandis_.
4. Functional_characterization_physiological - Performing gene set enrichment analysis and identifying the gene expression candidates responsible for controlling the physiological processes occuring in the xylem of drought stressed _E. grandis_.
5. MapMan - Creating a custom MapMan mapping file using _E. grandis_ ENTREZ IDs.

All folders represent R projects that contain all the necessary code, and datasets required to repeat the analysis. 
