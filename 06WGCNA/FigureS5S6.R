# =============================================================================
# Co-expression networks (WGCNA) of each population, module-trait correlations and GO enrichment of heat stress modules
#
# Part 1: for each population separately, expression values are rlog-transformed (DESeq2, blind to the experimental design),
# samples are checked by hierarchical clustering, a soft-thresholding power is chosen from the scale-free topology fit,
# and a signed co-expression network is built with blockwiseModules (power 8 in all three populations).
# Part 2: module eigengenes are correlated with temperature, time, shore level and long heat stress (Figure S5).
# Part 3: GO enrichment of the modules correlated with long heat stress (|r| > 0.5), and the GO terms of the up-regulated
# modules compared across populations (Figure S6, Table S6).
#
# Input:  counts_filter.txt, group.txt (from 01counts_ind_filter), gene_go_description_all.txt (from 05GO)
# Output: net_North_coexpression.rds, net_Central_coexpression.rds, net_South_coexpression.rds (networks),
#         FigureS5.pdf, FigureS6.pdf, Table_S6.csv
# =============================================================================

setwd("/data/06WGCNA/")
library(DESeq2)
library(WGCNA)
library(dendextend)
library(colorspace)
library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(ggtree)
library(patchwork)

# Counts and sample information.
count <- read.table("../01counts_ind_filter/counts_filter.txt", header = TRUE, row.names = 1, check.names = FALSE)
group <- read.delim("../01counts_ind_filter/group.txt", stringsAsFactors = FALSE, check.names = FALSE)
rownames(group) <- group$id
group <- group[colnames(count), ]

# Heat stress treatment: long = 31°C for 3 h; no/short = all other treatments (17°C for 15 mins; 17°C for 3 h; 31°C for 15 mins).
group$treatment <- factor(ifelse(group$temperature == "31degree" & group$time == "3h", "long", "no/short"), levels = c("no/short", "long"))
table(group$population, group$treatment)

# Candidate soft-thresholding powers, and multithreading / WGCNA correlation for network construction.
powers <- c(1:10, seq(from = 12, to = 20, by = 2))
enableWGCNAThreads()
cor <- WGCNA::cor

############################################################
#Part 1: co-expression networks

#################### North

# rlog-transformed expression of North samples.
group_North <- group[group$population == "North", ]
dds_North <- DESeqDataSetFromMatrix(countData = count[, group_North$id], colData = group_North, design = ~ 1)
datExpr_North <- t(assay(rlog(dds_North, blind = TRUE)))

# Check for genes and samples with too many missing values or zero variance.
gsg_North <- goodSamplesGenes(datExpr_North, verbose = 3)
gsg_North$allOK

# Sample tree (average linkage on Euclidean distances).
sample_dend_North <- as.dendrogram(hclust(dist(datExpr_North), method = "average"))
group_dend_North <- group_North[labels(sample_dend_North), ]
sample_colors_North <- rep("grey", nrow(group_dend_North))
sample_colors_North[group_dend_North$treatment == "no/short" & group_dend_North$shore_level == "high"] <- darken("#00A1D5", 0.3)
sample_colors_North[group_dend_North$treatment == "no/short" & group_dend_North$shore_level == "low"] <- lighten("#00A1D5", 0.3)
sample_colors_North[group_dend_North$treatment == "long" & group_dend_North$shore_level == "high"] <- darken("#DC5E83", 0.3)
sample_colors_North[group_dend_North$treatment == "long" & group_dend_North$shore_level == "low"] <- lighten("#DC5E83", 0.3)
labels_colors(sample_dend_North) <- sample_colors_North
plot(hang.dendrogram(sample_dend_North, hang = 0.1), sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
legend("bottomright", inset = 0.05, title = "", legend = c("no/short, High shore", "no/short, Low shore", "long, High shore", "long, Low shore"),
       fill = c(darken("#00A1D5", 0.3), lighten("#00A1D5", 0.3), darken("#DC5E83", 0.3), lighten("#DC5E83", 0.3)), horiz = FALSE, cex = 1.5, bty = "n")

# Soft-thresholding power: scale-free topology fit (signed R^2, target 0.90) and mean connectivity for each candidate power.
sft_North <- pickSoftThreshold(datExpr_North, powerVector = powers, verbose = 5, networkType = "signed")
sft_North$fitIndices[, c("Power", "SFT.R.sq", "mean.k.", "slope")]
par(mfrow = c(1, 2))
plot(sft_North$fitIndices[, 1], -sign(sft_North$fitIndices[, 3]) * sft_North$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n", main = "North Scale independence")
text(sft_North$fitIndices[, 1], -sign(sft_North$fitIndices[, 3]) * sft_North$fitIndices[, 2], labels = powers, col = "red")
abline(h = 0.90, col = "red")
plot(sft_North$fitIndices[, 1], sft_North$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "North Mean connectivity")
text(sft_North$fitIndices[, 1], sft_North$fitIndices[, 5], labels = powers, col = "red")
abline(h = 100, col = "red")
abline(h = 500, col = "red")
par(mfrow = c(1, 1))
#save(datExpr_North, file = "datExpr_North.RData")

# Signed co-expression network with power 8 (all genes in one block; modules of at least 50 genes; similar modules merged at 0.25).
net_North <- blockwiseModules(datExpr_North, power = 8, maxBlockSize = 25000, corType = "pearson", networkType = "signed", TOMType = "signed",
                              minModuleSize = 50, mergeCutHeight = 0.25, numericLabels = TRUE, saveTOMs = FALSE, verbose = 3, randomSeed = 12345)
#saveRDS(net_North, "net_North_coexpression.rds")


#################### Central

# rlog-transformed expression of Central samples.
group_Central <- group[group$population == "Central", ]
dds_Central <- DESeqDataSetFromMatrix(countData = count[, group_Central$id], colData = group_Central, design = ~ 1)
datExpr_Central <- t(assay(rlog(dds_Central, blind = TRUE)))

# Check for genes and samples with too many missing values or zero variance.
gsg_Central <- goodSamplesGenes(datExpr_Central, verbose = 3)
gsg_Central$allOK

# Sample tree (average linkage on Euclidean distances).
sample_dend_Central <- as.dendrogram(hclust(dist(datExpr_Central), method = "average"))
group_dend_Central <- group_Central[labels(sample_dend_Central), ]
sample_colors_Central <- rep("grey", nrow(group_dend_Central))
sample_colors_Central[group_dend_Central$treatment == "no/short" & group_dend_Central$shore_level == "high"] <- darken("#00A1D5", 0.3)
sample_colors_Central[group_dend_Central$treatment == "no/short" & group_dend_Central$shore_level == "low"] <- lighten("#00A1D5", 0.3)
sample_colors_Central[group_dend_Central$treatment == "long" & group_dend_Central$shore_level == "high"] <- darken("#DC5E83", 0.3)
sample_colors_Central[group_dend_Central$treatment == "long" & group_dend_Central$shore_level == "low"] <- lighten("#DC5E83", 0.3)
labels_colors(sample_dend_Central) <- sample_colors_Central
plot(hang.dendrogram(sample_dend_Central, hang = 0.1), sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
legend("bottomright", inset = 0.05, title = "", legend = c("no/short, High shore", "no/short, Low shore", "long, High shore", "long, Low shore"),
       fill = c(darken("#00A1D5", 0.3), lighten("#00A1D5", 0.3), darken("#DC5E83", 0.3), lighten("#DC5E83", 0.3)), horiz = FALSE, cex = 1.5, bty = "n")

# Soft-thresholding power: scale-free topology fit (signed R^2, target 0.90) and mean connectivity for each candidate power.
sft_Central <- pickSoftThreshold(datExpr_Central, powerVector = powers, verbose = 5, networkType = "signed")
sft_Central$fitIndices[, c("Power", "SFT.R.sq", "mean.k.", "slope")]
par(mfrow = c(1, 2))
plot(sft_Central$fitIndices[, 1], -sign(sft_Central$fitIndices[, 3]) * sft_Central$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n", main = "Central Scale independence")
text(sft_Central$fitIndices[, 1], -sign(sft_Central$fitIndices[, 3]) * sft_Central$fitIndices[, 2], labels = powers, col = "red")
abline(h = 0.90, col = "red")
plot(sft_Central$fitIndices[, 1], sft_Central$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Central Mean connectivity")
text(sft_Central$fitIndices[, 1], sft_Central$fitIndices[, 5], labels = powers, col = "red")
abline(h = 100, col = "red")
abline(h = 500, col = "red")
par(mfrow = c(1, 1))
#save(datExpr_Central, file = "datExpr_Central.RData")

# Signed co-expression network with power 8 (all genes in one block; modules of at least 50 genes; similar modules merged at 0.25).
net_Central <- blockwiseModules(datExpr_Central, power = 8, maxBlockSize = 25000, corType = "pearson", networkType = "signed", TOMType = "signed",
                                minModuleSize = 50, mergeCutHeight = 0.25, numericLabels = TRUE, saveTOMs = FALSE, verbose = 3, randomSeed = 12345)
#saveRDS(net_Central, "net_Central_coexpression.rds")


#################### South

# rlog-transformed expression of South samples.
group_South <- group[group$population == "South", ]
dds_South <- DESeqDataSetFromMatrix(countData = count[, group_South$id], colData = group_South, design = ~ 1)
datExpr_South <- t(assay(rlog(dds_South, blind = TRUE)))

# Check for genes and samples with too many missing values or zero variance.
gsg_South <- goodSamplesGenes(datExpr_South, verbose = 3)
gsg_South$allOK

# Sample tree (average linkage on Euclidean distances), coloured by heat stress treatment and shore level, to check for outliers.
sample_dend_South <- as.dendrogram(hclust(dist(datExpr_South), method = "average"))
group_dend_South <- group_South[labels(sample_dend_South), ]
sample_colors_South <- rep("grey", nrow(group_dend_South))
sample_colors_South[group_dend_South$treatment == "no/short" & group_dend_South$shore_level == "high"] <- darken("#00A1D5", 0.3)
sample_colors_South[group_dend_South$treatment == "no/short" & group_dend_South$shore_level == "low"] <- lighten("#00A1D5", 0.3)
sample_colors_South[group_dend_South$treatment == "long" & group_dend_South$shore_level == "high"] <- darken("#DC5E83", 0.3)
sample_colors_South[group_dend_South$treatment == "long" & group_dend_South$shore_level == "low"] <- lighten("#DC5E83", 0.3)
labels_colors(sample_dend_South) <- sample_colors_South
plot(hang.dendrogram(sample_dend_South, hang = 0.1), sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
legend("bottomright", inset = 0.05, title = "", legend = c("no/short, High shore", "no/short, Low shore", "long, High shore", "long, Low shore"),
       fill = c(darken("#00A1D5", 0.3), lighten("#00A1D5", 0.3), darken("#DC5E83", 0.3), lighten("#DC5E83", 0.3)), horiz = FALSE, cex = 1.5, bty = "n")

# Soft-thresholding power: scale-free topology fit (signed R^2, target 0.90) and mean connectivity for each candidate power.
sft_South <- pickSoftThreshold(datExpr_South, powerVector = powers, verbose = 5, networkType = "signed")
sft_South$fitIndices[, c("Power", "SFT.R.sq", "mean.k.", "slope")]
par(mfrow = c(1, 2))
plot(sft_South$fitIndices[, 1], -sign(sft_South$fitIndices[, 3]) * sft_South$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n", main = "South Scale independence")
text(sft_South$fitIndices[, 1], -sign(sft_South$fitIndices[, 3]) * sft_South$fitIndices[, 2], labels = powers, col = "red")
abline(h = 0.90, col = "red")
plot(sft_South$fitIndices[, 1], sft_South$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "South Mean connectivity")
text(sft_South$fitIndices[, 1], sft_South$fitIndices[, 5], labels = powers, col = "red")
abline(h = 100, col = "red")
abline(h = 500, col = "red")
par(mfrow = c(1, 1))
#save(datExpr_South, file = "datExpr_South.RData")

# Signed co-expression network with power 8 (all genes in one block; modules of at least 50 genes; similar modules merged at 0.25).
net_South <- blockwiseModules(datExpr_South, power = 8, maxBlockSize = 25000, corType = "pearson", networkType = "signed", TOMType = "signed",
                              minModuleSize = 50, mergeCutHeight = 0.25, numericLabels = TRUE, saveTOMs = FALSE, verbose = 3, randomSeed = 12345)
#saveRDS(net_South, "net_South_coexpression.rds")


############################################################
#Part 2: module-trait correlations (Figure S5)

# Networks read back here.
net_North <- readRDS("net_North_coexpression.rds")
net_Central <- readRDS("net_Central_coexpression.rds")
net_South <- readRDS("net_South_coexpression.rds")

# Traits correlated with the module eigengenes (column names in the trait tables; 1 = yes, 0 = no) and their labels in Figure S5.
trait_labels <- c("Temperature (31°C)", "Time (3 hours)", "Shore (high)", "Long heat stress")
trait_names <- c("temperature_31", "time_3h", "shore_high", "treatment_long")

#################### North

# Module colours from the network, renamed Module_1 to Module_14 in the order used in the paper (grey = genes not assigned to a module).
group_North <- group[rownames(datExpr_North), ]
stopifnot(identical(names(net_North$colors), colnames(datExpr_North)))
module_colors_North <- labels2colors(net_North$colors)
module_name_map_North <- setNames(paste0("Module_", 1:14), c("brown", "blue", "turquoise", "green", "red", "yellow", "black", "pink", "tan", "cyan", "magenta", "greenyellow", "purple", "salmon"))
ordered_colors_North <- names(module_name_map_North)
stopifnot(setequal(setdiff(unique(module_colors_North), "grey"), ordered_colors_North))

# Module eigengenes (first principal component of each module) and their Pearson correlation with each trait.
MEs_North <- moduleEigengenes(datExpr_North, colors = module_colors_North, excludeGrey = TRUE, align = "along average")$eigengenes
MEs_North <- MEs_North[, paste0("ME", ordered_colors_North)]
traits_North <- data.frame(temperature_31 = as.integer(group_North$temperature == "31degree"), time_3h = as.integer(group_North$time == "3h"),
                           shore_high = as.integer(group_North$shore_level == "high"), treatment_long = as.integer(group_North$treatment == "long"), row.names = group_North$id)
moduleTraitCor_North <- WGCNA::cor(MEs_North, traits_North, method = "pearson", use = "all.obs")
moduleTraitPvalue_North <- corPvalueStudent(moduleTraitCor_North, nSamples = nrow(MEs_North))
module_size_North <- as.integer(table(factor(module_colors_North, levels = ordered_colors_North)))

# Heatmap: correlation (P value) for each module x trait; module colour shown as a stripe on the left; Module 1 at the top.
heatmap_North <- expand.grid(module_number = 1:14, trait = trait_names, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
heatmap_North$module_color <- ordered_colors_North[heatmap_North$module_number]
heatmap_North$r <- as.vector(moduleTraitCor_North)
heatmap_North$pvalue <- as.vector(moduleTraitPvalue_North)
heatmap_North$x <- match(heatmap_North$trait, trait_names)
heatmap_North$y <- 15 - heatmap_North$module_number
heatmap_North$text <- paste0(signif(heatmap_North$r, 2), " (", signif(heatmap_North$pvalue, 2), ")")
stripe_North <- heatmap_North[heatmap_North$trait == "treatment_long", ]
pr_North <- ggplot(heatmap_North, aes(x, y)) + geom_tile(aes(fill = r)) +
  geom_tile(data = stripe_North, aes(x = 0.34, y = y), fill = stripe_North$module_color, width = 0.18, height = 1) + geom_text(aes(label = text), size = 2.8) +
  scale_fill_gradientn(colors = blueWhiteRed(50), limits = c(-1, 1), name = "r") +
  scale_x_continuous(breaks = 1:4, labels = trait_labels, limits = c(0.2, 4.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = 14:1, labels = paste0("Module ", 1:14, " (", module_size_North, ")"), expand = c(0, 0.5)) +
  labs(title = "North", x = NULL, y = NULL) + theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"), axis.text.x = element_text(face = "bold"), axis.text.y = element_text(hjust = 0),
        plot.title = element_text(hjust = 0.5, face = "bold"))
pr_North

# Genes of each module, and the modules correlated with long heat stress (|r| > 0.5).
genes_North <- data.frame(gene = colnames(datExpr_North), module_name = unname(module_name_map_North[module_colors_North]), stringsAsFactors = FALSE)
table(genes_North$module_name, useNA = "ifany")
long_North <- data.frame(population = "North", module_name = unname(module_name_map_North), module_size = module_size_North,
                         r = moduleTraitCor_North[, "treatment_long"], pvalue = moduleTraitPvalue_North[, "treatment_long"], row.names = NULL)
long_North[abs(long_North$r) > 0.5, ]

#################### Central

# Module colours from the network, renamed Module_1 to Module_15 in the order used in the paper (grey = genes not assigned to a module).
group_Central <- group[rownames(datExpr_Central), ]
stopifnot(identical(names(net_Central$colors), colnames(datExpr_Central)))
module_colors_Central <- labels2colors(net_Central$colors)
module_name_map_Central <- setNames(paste0("Module_", 1:15), c("black", "blue", "turquoise", "brown", "green", "magenta", "red", "midnightblue", "cyan", "salmon", "purple", "tan", "yellow", "greenyellow", "pink"))
ordered_colors_Central <- names(module_name_map_Central)
stopifnot(setequal(setdiff(unique(module_colors_Central), "grey"), ordered_colors_Central))

# Module eigengenes (first principal component of each module) and their Pearson correlation with each trait.
MEs_Central <- moduleEigengenes(datExpr_Central, colors = module_colors_Central, excludeGrey = TRUE, align = "along average")$eigengenes
MEs_Central <- MEs_Central[, paste0("ME", ordered_colors_Central)]
traits_Central <- data.frame(temperature_31 = as.integer(group_Central$temperature == "31degree"), time_3h = as.integer(group_Central$time == "3h"),
                             shore_high = as.integer(group_Central$shore_level == "high"), treatment_long = as.integer(group_Central$treatment == "long"), row.names = group_Central$id)
moduleTraitCor_Central <- WGCNA::cor(MEs_Central, traits_Central, method = "pearson", use = "all.obs")
moduleTraitPvalue_Central <- corPvalueStudent(moduleTraitCor_Central, nSamples = nrow(MEs_Central))
module_size_Central <- as.integer(table(factor(module_colors_Central, levels = ordered_colors_Central)))

# Heatmap: correlation (P value) for each module x trait; module colour shown as a stripe on the left; Module 1 at the top.
heatmap_Central <- expand.grid(module_number = 1:15, trait = trait_names, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
heatmap_Central$module_color <- ordered_colors_Central[heatmap_Central$module_number]
heatmap_Central$r <- as.vector(moduleTraitCor_Central)
heatmap_Central$pvalue <- as.vector(moduleTraitPvalue_Central)
heatmap_Central$x <- match(heatmap_Central$trait, trait_names)
heatmap_Central$y <- 16 - heatmap_Central$module_number
heatmap_Central$text <- paste0(signif(heatmap_Central$r, 2), " (", signif(heatmap_Central$pvalue, 2), ")")
stripe_Central <- heatmap_Central[heatmap_Central$trait == "treatment_long", ]
pr_Central <- ggplot(heatmap_Central, aes(x, y)) + geom_tile(aes(fill = r)) +
  geom_tile(data = stripe_Central, aes(x = 0.34, y = y), fill = stripe_Central$module_color, width = 0.18, height = 1) + geom_text(aes(label = text), size = 2.8) +
  scale_fill_gradientn(colors = blueWhiteRed(50), limits = c(-1, 1), name = "r") +
  scale_x_continuous(breaks = 1:4, labels = trait_labels, limits = c(0.2, 4.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = 15:1, labels = paste0("Module ", 1:15, " (", module_size_Central, ")"), expand = c(0, 0.5)) +
  labs(title = "Central", x = NULL, y = NULL) + theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"), axis.text.x = element_text(face = "bold"), axis.text.y = element_text(hjust = 0),
        plot.title = element_text(hjust = 0.5, face = "bold"))
pr_Central

# Genes of each module, and the modules correlated with long heat stress (|r| > 0.5).
genes_Central <- data.frame(gene = colnames(datExpr_Central), module_name = unname(module_name_map_Central[module_colors_Central]), stringsAsFactors = FALSE)
table(genes_Central$module_name, useNA = "ifany")
long_Central <- data.frame(population = "Central", module_name = unname(module_name_map_Central), module_size = module_size_Central,
                           r = moduleTraitCor_Central[, "treatment_long"], pvalue = moduleTraitPvalue_Central[, "treatment_long"], row.names = NULL)
long_Central[abs(long_Central$r) > 0.5, ]

#################### South

# Module colours from the network, renamed Module_1 to Module_16 in the order used in the paper (grey = genes not assigned to a module).
group_South <- group[rownames(datExpr_South), ]
stopifnot(identical(names(net_South$colors), colnames(datExpr_South)))
module_colors_South <- labels2colors(net_South$colors)
module_name_map_South <- setNames(paste0("Module_", 1:16), c("red", "turquoise", "blue", "brown", "green", "black", "purple", "cyan", "greenyellow", "lightcyan", "magenta", "midnightblue", "salmon", "pink", "tan", "yellow"))
ordered_colors_South <- names(module_name_map_South)
stopifnot(setequal(setdiff(unique(module_colors_South), "grey"), ordered_colors_South))

# Module eigengenes (first principal component of each module) and their Pearson correlation with each trait.
MEs_South <- moduleEigengenes(datExpr_South, colors = module_colors_South, excludeGrey = TRUE, align = "along average")$eigengenes
MEs_South <- MEs_South[, paste0("ME", ordered_colors_South)]
traits_South <- data.frame(temperature_31 = as.integer(group_South$temperature == "31degree"), time_3h = as.integer(group_South$time == "3h"),
                           shore_high = as.integer(group_South$shore_level == "high"), treatment_long = as.integer(group_South$treatment == "long"), row.names = group_South$id)
moduleTraitCor_South <- WGCNA::cor(MEs_South, traits_South, method = "pearson", use = "all.obs")
moduleTraitPvalue_South <- corPvalueStudent(moduleTraitCor_South, nSamples = nrow(MEs_South))
module_size_South <- as.integer(table(factor(module_colors_South, levels = ordered_colors_South)))

# Heatmap: correlation (P value) for each module x trait; module colour shown as a stripe on the left; Module 1 at the top.
heatmap_South <- expand.grid(module_number = 1:16, trait = trait_names, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
heatmap_South$module_color <- ordered_colors_South[heatmap_South$module_number]
heatmap_South$r <- as.vector(moduleTraitCor_South)
heatmap_South$pvalue <- as.vector(moduleTraitPvalue_South)
heatmap_South$x <- match(heatmap_South$trait, trait_names)
heatmap_South$y <- 17 - heatmap_South$module_number
heatmap_South$text <- paste0(signif(heatmap_South$r, 2), " (", signif(heatmap_South$pvalue, 2), ")")
stripe_South <- heatmap_South[heatmap_South$trait == "treatment_long", ]
pr_South <- ggplot(heatmap_South, aes(x, y)) + geom_tile(aes(fill = r)) +
  geom_tile(data = stripe_South, aes(x = 0.34, y = y), fill = stripe_South$module_color, width = 0.18, height = 1) + geom_text(aes(label = text), size = 2.8) +
  scale_fill_gradientn(colors = blueWhiteRed(50), limits = c(-1, 1), name = "r") +
  scale_x_continuous(breaks = 1:4, labels = trait_labels, limits = c(0.2, 4.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = 16:1, labels = paste0("Module ", 1:16, " (", module_size_South, ")"), expand = c(0, 0.5)) +
  labs(title = "South", x = NULL, y = NULL) + theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"), axis.text.x = element_text(face = "bold"), axis.text.y = element_text(hjust = 0),
        plot.title = element_text(hjust = 0.5, face = "bold"))
pr_South

# Genes of each module, and the modules correlated with long heat stress (|r| > 0.5).
genes_South <- data.frame(gene = colnames(datExpr_South), module_name = unname(module_name_map_South[module_colors_South]), stringsAsFactors = FALSE)
table(genes_South$module_name, useNA = "ifany")
long_South <- data.frame(population = "South", module_name = unname(module_name_map_South), module_size = module_size_South,
                         r = moduleTraitCor_South[, "treatment_long"], pvalue = moduleTraitPvalue_South[, "treatment_long"], row.names = NULL)
long_South[abs(long_South$r) > 0.5, ]

#################### Figure S5

figS5 <- pr_North / pr_Central / pr_South + plot_annotation(tag_levels = "A")
figS5
ggsave("FigureS5.pdf", figS5, width = 12, height = 15)
ggsave("FigureS5.png", figS5, width = 12, height = 15)

############################################################
#Part 3: GO enrichment of the modules correlated with long heat stress (Figure S6, Table S6)

# Modules correlated with long heat stress (|r| > 0.5):
# North Modules 1, 2 and 13, Central Module 1, South Modules 1 and 16.

North_M1_genes <- genes_North$gene[genes_North$module_name %in% "Module_1"]
North_M2_genes <- genes_North$gene[genes_North$module_name %in% "Module_2"]
North_M13_genes <- genes_North$gene[genes_North$module_name %in% "Module_13"]
Central_M1_genes <- genes_Central$gene[genes_Central$module_name %in% "Module_1"]
South_M1_genes <- genes_South$gene[genes_South$module_name %in% "Module_1"]
South_M16_genes <- genes_South$gene[genes_South$module_name %in% "Module_16"]

# Gene-GO annotation (from 05GO); background = retained genes with at least one GO term.
gene_GO <- read.delim("../05GO/gene_go_description_all.txt", header = TRUE, stringsAsFactors = FALSE, quote = "")
gene_GO$gene <- trimws(gene_GO$gene)
gene_GO$go <- trimws(gene_GO$go)
gene_GO$description <- trimws(gene_GO$description)
gene_GO <- gene_GO[!is.na(gene_GO$gene) & gene_GO$gene != "" & grepl("^GO:[0-9]{7}$", gene_GO$go), ]
TERM2GENE <- distinct(gene_GO, go, gene)
TERM2NAME <- distinct(gene_GO, go, description)
GO_universe <- intersect(rownames(count), unique(gene_GO$gene))
length(GO_universe)

# GO enrichment of each module; a GO term is significant when p.adjust <= 0.01 and qvalue <= 0.01.


# North Module 1
length(North_M1_genes)
go_North_M1 <- enricher(gene = North_M1_genes, universe = GO_universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
North_M1_all <- mutate(go_North_M1@result, Group = "North_M1")
sig_North_M1 <- North_M1_all[North_M1_all$p.adjust <= 0.01 & !is.na(North_M1_all$qvalue) & North_M1_all$qvalue <= 0.01, ]
nrow(sig_North_M1)

# North Module 2
length(North_M2_genes)
go_North_M2 <- enricher(gene = North_M2_genes, universe = GO_universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
North_M2_all <- mutate(go_North_M2@result, Group = "North_M2")
sig_North_M2 <- North_M2_all[North_M2_all$p.adjust <= 0.01 & !is.na(North_M2_all$qvalue) & North_M2_all$qvalue <= 0.01, ]
nrow(sig_North_M2)

# North Module 13
length(North_M13_genes)
go_North_M13 <- enricher(gene = North_M13_genes, universe = GO_universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
North_M13_all <- mutate(go_North_M13@result, Group = "North_M13")
sig_North_M13 <- North_M13_all[North_M13_all$p.adjust <= 0.01 & !is.na(North_M13_all$qvalue) & North_M13_all$qvalue <= 0.01, ]
nrow(sig_North_M13)

# Central Module 1
length(Central_M1_genes)
go_Central_M1 <- enricher(gene = Central_M1_genes, universe = GO_universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
Central_M1_all <- mutate(go_Central_M1@result, Group = "Central_M1")
sig_Central_M1 <- Central_M1_all[Central_M1_all$p.adjust <= 0.01 & !is.na(Central_M1_all$qvalue) & Central_M1_all$qvalue <= 0.01, ]
nrow(sig_Central_M1)

# South Module 1
length(South_M1_genes)
go_South_M1 <- enricher(gene = South_M1_genes, universe = GO_universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
South_M1_all <- mutate(go_South_M1@result, Group = "South_M1")
sig_South_M1 <- South_M1_all[South_M1_all$p.adjust <= 0.01 & !is.na(South_M1_all$qvalue) & South_M1_all$qvalue <= 0.01, ]
nrow(sig_South_M1)

# South Module 16
length(South_M16_genes)
go_South_M16 <- enricher(gene = South_M16_genes, universe = GO_universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
South_M16_all <- mutate(go_South_M16@result, Group = "South_M16")
sig_South_M16 <- South_M16_all[South_M16_all$p.adjust <= 0.01 & !is.na(South_M16_all$qvalue) & South_M16_all$qvalue <= 0.01, ]
nrow(sig_South_M16)

# North Module 13 and South Module 16 have no significantly enriched GO terms, so the following analyses focus only on
# the modules whose genes are up-regulated under long heat stress (North Modules 1 and 2, Central Module 1, South Module 1).

#################### Figure S6: GO terms enriched in the up-regulated modules, ordered by a GO tree

# All GO terms significant in at least one of the four modules.
sig_GO_up <- bind_rows(sig_North_M1, sig_North_M2, sig_Central_M1, sig_South_M1)
go_WGCNA <- distinct(sig_GO_up, ID, Description)
GO_ids <- go_WGCNA$ID
nrow(go_WGCNA)

# Enrichment object containing exactly these GO terms (all cutoffs set to 1, no term is dropped), only to cluster them by
# shared annotated genes (Jaccard similarity); its P values are not used.
gene_GO_WGCNA <- gene_GO[gene_GO$go %in% GO_ids, ]
go_WGCNA_all <- enricher(gene = unique(gene_GO_WGCNA$gene), TERM2GENE = gene_GO_WGCNA[c("go", "gene")], TERM2NAME = gene_GO_WGCNA[c("go", "description")],
                         minGSSize = 1, maxGSSize = Inf, pvalueCutoff = 1, qvalueCutoff = 1)
stopifnot(setequal(go_WGCNA_all@result$ID, GO_ids), !anyDuplicated(go_WGCNA_all@result$Description))
N_GO <- nrow(go_WGCNA_all@result)
go_WGCNA_sim <- pairwise_termsim(go_WGCNA_all, method = "JC", showCategory = N_GO)
tree_WGCNA <- treeplot(go_WGCNA_sim, showCategory = N_GO, cluster.params = list(method = "average", n = 4, color = c("#FDBF6F", "#B2DF8A", "#CAB2D6", "#a6cee3"), label_words_n = 0),
                       hilight.params = list(hilight = FALSE))

# Four functional clusters, each defined by the common ancestor of two GO terms at its ends:
# heat response, ribosome-related process, Tf2-like and mitochondrial respiration.
clade_anchor_IDs <- c("GO:0042594", "GO:0051085", "GO:0031597", "GO:0022626", "GO:0002040", "GO:0039660", "GO:0006123", "GO:0045271")
clade_anchor_labels <- go_WGCNA_all@result$Description[match(clade_anchor_IDs, go_WGCNA_all@result$ID)]
node_heat <- ggtree::MRCA(tree_WGCNA, clade_anchor_labels[1], clade_anchor_labels[2])
node_ribosome <- ggtree::MRCA(tree_WGCNA, clade_anchor_labels[3], clade_anchor_labels[4])
node_Tf2 <- ggtree::MRCA(tree_WGCNA, clade_anchor_labels[5], clade_anchor_labels[6])
node_mito <- ggtree::MRCA(tree_WGCNA, clade_anchor_labels[7], clade_anchor_labels[8])
tree_WGCNA <- ggtree::flip(tree_WGCNA, node_heat, node_ribosome)
tree_WGCNA <- groupClade(tree_WGCNA, c(node_heat, node_ribosome, node_Tf2, node_mito), group_name = "group")

# GO terms in tree order (bottom to top) with their cluster; check that each cluster contains its representative GO term.
tips_WGCNA <- arrange(filter(tree_WGCNA$data, isTip), y)
tips_WGCNA$ID <- go_WGCNA_all@result$ID[match(tips_WGCNA$label, go_WGCNA_all@result$Description)]
stopifnot(!anyNA(tips_WGCNA$ID), setequal(tips_WGCNA$ID, GO_ids),
          identical(as.character(tips_WGCNA$group[match(c("GO:0051085", "GO:0022626", "GO:0015074", "GO:0006120"), tips_WGCNA$ID)]), c("1", "2", "3", "4")))

# Bubbles: gene count and FDR of every GO term in each module, taken from the full results (grey = not significant;
# count 0 and no bubble = no module gene annotated to the term, so the term was not tested).
GO_plot_all <- filter(bind_rows(North_M1_all, North_M2_all, Central_M1_all, South_M1_all), ID %in% GO_ids)
GO_plot_all <- left_join(expand_grid(ID = tips_WGCNA$ID, Group = c("North_M1", "North_M2", "Central_M1", "South_M1")), GO_plot_all, by = c("ID", "Group"))
GO_plot_all$Description <- go_WGCNA$Description[match(GO_plot_all$ID, go_WGCNA$ID)]
GO_plot_all$Count[is.na(GO_plot_all$Count)] <- 0
GO_plot_all$significant <- !is.na(GO_plot_all$p.adjust) & GO_plot_all$p.adjust <= 0.01 & !is.na(GO_plot_all$qvalue) & GO_plot_all$qvalue <= 0.01
GO_plot_all$ColorCol <- ifelse(GO_plot_all$significant, -log10(pmax(GO_plot_all$p.adjust, .Machine$double.xmin)), NA_real_)
GO_plot_all$y <- tips_WGCNA$y[match(GO_plot_all$ID, tips_WGCNA$ID)]
GO_plot_all$x <- unname(c(North_M1 = 7.52, North_M2 = 8.82, Central_M1 = 10.42, South_M1 = 12.02)[GO_plot_all$Group])

# Tree branches coloured by cluster, with the tree width rescaled to 5.57.
nodes_WGCNA <- as.data.frame(tree_WGCNA$data)
nodes_WGCNA$x <- nodes_WGCNA$x / max(nodes_WGCNA$x) * 5.57
branch_colors <- c("0" = "black", "1" = "#FDBF6F", "2" = "#B2DF8A", "3" = "#CAB2D6", "4" = "#a6cee3")
edges_WGCNA <- nodes_WGCNA[nodes_WGCNA$node != nodes_WGCNA$parent, ]
edges_WGCNA$x_parent <- nodes_WGCNA$x[match(edges_WGCNA$parent, nodes_WGCNA$node)]
edges_WGCNA$y_parent <- nodes_WGCNA$y[match(edges_WGCNA$parent, nodes_WGCNA$node)]
edges_WGCNA$color <- unname(branch_colors[as.character(edges_WGCNA$group)])

# Coloured background behind the GO IDs of each cluster, and cluster names on the right.
bands_WGCNA <- summarise(group_by(filter(tips_WGCNA, as.character(group) != "0"), group), ymin = min(y) - 0.5, ymax = max(y) + 0.5, .groups = "drop")
bands_WGCNA$fill <- unname(branch_colors[as.character(bands_WGCNA$group)])
bands_WGCNA$label <- unname(c("1" = "Heat response genes", "2" = "Ribosome-related process", "3" = "Tf2-like", "4" = "Mitochondrial\nrespiration")[as.character(bands_WGCNA$group)])
bands_WGCNA$text_color <- unname(c("1" = "#918009", "2" = "#3A6831", "3" = "#CF92ED", "4" = "#82B6D1")[as.character(bands_WGCNA$group)])
bands_WGCNA$mid <- (bands_WGCNA$ymin + bands_WGCNA$ymax) / 2
bands_WGCNA$x_text <- ifelse(as.character(bands_WGCNA$group) == "3", 13.37, 13.48)

# Population headers (top), module names (bottom), and a further group of GO terms outside the four clusters (labelled only).
headers_WGCNA <- data.frame(population = c("North", "Central", "South"), xmin = c(6.87, 9.67, 11.27), xmax = c(9.47, 11.17, 12.77), x = c(8.17, 10.42, 12.02))
module_headers <- data.frame(x = c(7.52, 8.82, 10.42, 12.02), label = c("Module 1", "Module 2", "Module 1", "Module 1"))
viral_y <- range(tips_WGCNA$y[tips_WGCNA$ID %in% c("GO:0019068", "GO:0042575")])

# Tree, GO IDs, bubbles and labels share one y axis, so every row lines up.
figS6 <- ggplot() +
  geom_segment(data = edges_WGCNA, aes(x = x_parent, xend = x, y = y, yend = y), color = edges_WGCNA$color, linewidth = 0.85, lineend = "butt") +
  geom_segment(data = edges_WGCNA, aes(x = x_parent, xend = x_parent, y = y_parent, yend = y), color = edges_WGCNA$color, linewidth = 0.85, lineend = "butt") +
  geom_rect(data = bands_WGCNA, aes(ymin = ymin, ymax = ymax), xmin = 5.62, xmax = 6.82, fill = bands_WGCNA$fill, alpha = 0.45, color = NA) +
  geom_text(data = tips_WGCNA, aes(y = y, label = ID), x = (5.62 + 6.82) / 2, hjust = 0.5, size = 2.7) +
  geom_segment(data = tips_WGCNA, aes(y = y, yend = y), x = 6.87, xend = 12.77, color = "grey90", linewidth = 0.2) +
  geom_rect(data = headers_WGCNA, aes(xmin = xmin, xmax = xmax), ymin = N_GO + 0.7, ymax = N_GO + 2.7, fill = "grey90", color = "black", linewidth = 0.25) +
  geom_text(data = headers_WGCNA, aes(x = x, label = population), y = N_GO + 1.7, size = 3.8) +
  geom_point(data = GO_plot_all[GO_plot_all$Count > 0, ], aes(x = x, y = y, color = ColorCol, size = Count)) +
  scale_size_continuous(range = c(1, 4), breaks = c(50, 100), name = "Gene count",
                        guide = guide_legend(order = 2, nrow = 1, title.position = "left", override.aes = list(shape = 21, fill = NA, color = "black", stroke = 0.4))) +
  scale_color_gradient(low = "#97D1A0", high = "#1E803D", na.value = "grey80", breaks = c(50, 100), name = expression(plain("-log")[10]~FDR),
                       guide = guide_colorbar(order = 1, direction = "horizontal", title.position = "left", barwidth = grid::unit(1.4, "cm"), barheight = grid::unit(0.3, "cm"))) +
  geom_text(data = module_headers, aes(x = x, label = label), y = -0.6, size = 3.6) +
  geom_segment(data = bands_WGCNA, aes(y = ymin + 0.2, yend = ymax - 0.2), x = 12.95, xend = 12.95, color = bands_WGCNA$fill, linewidth = 0.6) +
  geom_text(data = bands_WGCNA, aes(x = x_text, y = mid, label = label), angle = 90, fontface = "bold", size = 3.7, color = bands_WGCNA$text_color) +
  annotate("segment", x = 13.50, xend = 13.50, y = viral_y[1], yend = viral_y[2], linewidth = 0.5) +
  annotate("text", x = 13.88, y = mean(viral_y), label = "Viral processes and DNA metabolism", angle = 90, size = 3.5, fontface = "bold") +
  scale_x_continuous(limits = c(-0.1, 14.2), expand = c(0, 0)) + scale_y_continuous(limits = c(-6.5, N_GO + 4), expand = c(0, 0)) +
  theme_void() + theme(legend.position = "inside", legend.position.inside = c(0.65, 0.015), legend.justification = c(0.5, 0), legend.box = "horizontal",
                       legend.title = element_text(size = 9), legend.text = element_text(size = 8), plot.margin = margin(21, 12, 11, 20),
                       plot.background = element_rect(fill = "white", color = NA))
figS6
ggsave("FigureS6.pdf", figS6, width = 11, height = max(7, N_GO * 0.14 + 0.74), limitsize = FALSE)
ggsave("FigureS6.png", figS6, width = 11, height = max(7, N_GO * 0.14 + 0.74), limitsize = FALSE)


#################### Table S6: gene count and FDR of each GO term in the four modules, in the order of Figure S6 (top to bottom)

table_S6 <- pivot_wider(GO_plot_all[, c("ID", "Description", "Group", "Count", "p.adjust")], names_from = Group, values_from = c(Count, p.adjust), names_glue = "{.value}_{Group}")
table_S6 <- table_S6[match(rev(tips_WGCNA$ID), table_S6$ID), c("ID", "Description", "Count_North_M1", "p.adjust_North_M1", "Count_North_M2", "p.adjust_North_M2",
                                                               "Count_Central_M1", "p.adjust_Central_M1", "Count_South_M1", "p.adjust_South_M1")]
colnames(table_S6) <- c("GO_ID", "GO_term_name", "Gene_count_North_Module1", "FDR_North_Module1", "Gene_count_North_Module2", "FDR_North_Module2",
                        "Gene_count_Central_Module1", "FDR_Central_Module1", "Gene_count_South_Module1", "FDR_South_Module1")

# FDR with 6 decimals; values below 1E-06 in scientific notation (3 significant digits), as in Table S5.
table_S6$FDR_North_Module1 <- ifelse(is.na(table_S6$FDR_North_Module1), NA, ifelse(table_S6$FDR_North_Module1 < 1e-6, toupper(formatC(table_S6$FDR_North_Module1, format = "e", digits = 2)), formatC(table_S6$FDR_North_Module1, format = "f", digits = 6)))
table_S6$FDR_North_Module2 <- ifelse(is.na(table_S6$FDR_North_Module2), NA, ifelse(table_S6$FDR_North_Module2 < 1e-6, toupper(formatC(table_S6$FDR_North_Module2, format = "e", digits = 2)), formatC(table_S6$FDR_North_Module2, format = "f", digits = 6)))
table_S6$FDR_Central_Module1 <- ifelse(is.na(table_S6$FDR_Central_Module1), NA, ifelse(table_S6$FDR_Central_Module1 < 1e-6, toupper(formatC(table_S6$FDR_Central_Module1, format = "e", digits = 2)), formatC(table_S6$FDR_Central_Module1, format = "f", digits = 6)))
table_S6$FDR_South_Module1 <- ifelse(is.na(table_S6$FDR_South_Module1), NA, ifelse(table_S6$FDR_South_Module1 < 1e-6, toupper(formatC(table_S6$FDR_South_Module1, format = "e", digits = 2)), formatC(table_S6$FDR_South_Module1, format = "f", digits = 6)))
dim(table_S6)
write.csv(table_S6, "Table_S6.csv", row.names = FALSE, na = "NA")


