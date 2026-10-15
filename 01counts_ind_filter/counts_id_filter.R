setwd("/data/01counts_ind_filter/")

library(DESeq2)
library(ggplot2)
library(ggrepel)

# Load raw gene counts from HTSeq-count.
counts <- read.table("counts_htseqcount.txt", header = TRUE)

# Retain genes with at least 10 counts in at least 5 samples in each population (North = NO, Central = SW, South = SP).
keep_North <- rowSums(counts[, grepl("^NO", colnames(counts))] >= 10) >= 5
keep_Central <- rowSums(counts[, grepl("^SW", colnames(counts))] >= 10) >= 5
keep_South <- rowSums(counts[, grepl("^SP", colnames(counts))] >= 10) >= 5
filter <- counts[keep_North & keep_Central & keep_South, ]

# Define PCA using all retained genes.
plot_pca <- function(filter, title) {
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(filter), colData = data.frame(row.names = colnames(filter)), design = ~1)
  vsd <- vst(dds, blind = TRUE)
  pca <- prcomp(t(assay(vsd)), center = TRUE, scale. = FALSE)
  percentVar <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pcaData <- data.frame(sample = rownames(pca$x), pca$x[, 1:2])
  ggplot(pcaData, aes(PC1, PC2, label = sample)) + geom_point(size = 2.5) + geom_text_repel(size = 3, max.overlaps = Inf, seed = 123) +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", percentVar[1]), y = sprintf("PC2 (%.1f%%)", percentVar[2])) + theme_classic()
}

# Inspect PCA across all samples and within each population.
print(plot_pca(filter, "All samples"))
print(plot_pca(filter[, grepl("^NO", colnames(filter)), drop = FALSE], "North"))
print(plot_pca(filter[, grepl("^SW", colnames(filter)), drop = FALSE], "Central"))
print(plot_pca(filter[, grepl("^SP", colnames(filter)), drop = FALSE], "South"))

# Exclude SPL47, identified as an outlier in the South PCA.
filter <- filter[, colnames(filter) != "SPL47", drop = FALSE]

# Recalculate PCA using the final filtered count matrix.
print(plot_pca(filter, "All samples (SPL47 excluded)"))
print(plot_pca(filter[, grepl("^SP", colnames(filter)), drop = FALSE], "South (SPL47 excluded)"))

# Export the final gene- and sample-filtered count matrix.

write.table(filter, "counts_filter.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)


# Remove SPL47 and align sample metadata with the final count matrix.
group <- read.table("group_102ids.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
group <- group[group$id != "SPL47", , drop = FALSE]
rownames(group) <- group$id
group <- group[colnames(filter), , drop = FALSE]

write.table(group, "group.txt", sep = "\t", quote = FALSE, row.names = FALSE)


