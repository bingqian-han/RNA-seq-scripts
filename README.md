# Geographical variation in gene regulatory responses to heat stress in a marine snail

Code and data for the manuscript by Bingqian Han et al.

Each folder contains one analysis: the R script, its input files and its output files.

- `01counts_ind_filter`: filtering of genes and individuals
- `02PCA`: PCA of gene expression
- `03VariancePartition`: variance partitioning
- `04Permutation`: permutation tests of treatment groups
- `05DEseq2`: differential expression analysis
- `06GOKEGG`: GO and KEGG enrichment
- `07WGCNA`: co-expression network analysis
- `08Frontloading`: frontloading of heat response genes
- `09Hsp70`: *Hsp* copy number
- `10TF2`: Tf2-like retrotransposons
- `11heartrate`: raw heart-rate data and analysis
- `12temperature_loggers`: raw field rock temperature logger data and analysis

Raw RNA-seq reads are available at NCBI SRA under BioProject PRJNA1167455, and raw pool-seq reads under BioProject PRJNA1399097.
