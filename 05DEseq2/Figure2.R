setwd("C:/path/to/RNAgithub/05DEseq2/")

# =============================================================================
# Differential expression within each population (DESeq2) and Figure 2.
# Four groups per population: shore level (high, low) x heat stress (no/short, long). Two contrasts define the general temperature
# response genes and the shore-level-specific response genes (padj < 0.05, |log2FC| >= 0.5).
# Also: overlap of these genes among populations (SuperExactTest) and their distribution across linkage groups.
# Input: counts_filter.txt, group.txt (from 01counts_ind_filter), LG_gene_position_20176.txt
# Output: results of all genes for each population and contrast (*_generaltemp_allgenes.txt, *_shorelevel_allgenes.txt,
#         *_hightemp_allgenes.txt, *_lowtemp_allgenes.txt), Figure2.pdf
# =============================================================================

library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(eulerr)
library(grid)
library(ggpubr)
library(SuperExactTest)


# Load the count matrix and sample information.
counts <- read.table("../01counts_ind_filter/counts_filter.txt", header = TRUE)
group <- read.table("../01counts_ind_filter/group.txt", header = TRUE, stringsAsFactors = FALSE)
rownames(group) <- group$id
group <- group[colnames(counts), ]

# Heat stress treatment: snails kept at 31°C for 3 hours form the long heat stress group;
# all other treatments (17°C for 15 mins, 17°C for 3 hours, 31°C for 15 mins) are pooled as no/short heat stress.
group$stress <- ifelse(group$temperature == "31degree" & group$time == "3h", "Long", "NoShort")
group$group <- paste0(ifelse(group$shore_level == "high", "High", "Low"), group$stress)
group$group <- factor(group$group, levels = c("HighNoShort", "HighLong", "LowNoShort", "LowLong"))


#################### North
# Differential expression among the four groups (shore level x heat stress) in the North population.
counts_North <- counts[, group$population == "North"]
group_North <- group[group$population == "North", ]
dds_North <- DESeqDataSetFromMatrix(countData = counts_North, colData = group_North, design = ~ 0 + group)
dds_North <- DESeq(dds_North)

# General response: 0.5 * ((HighLong - HighNoShort) + (LowLong - LowNoShort)); mean response across shore levels.
res_North_general <- results(dds_North, contrast = list(c("groupHighLong", "groupLowLong"), c("groupHighNoShort", "groupLowNoShort")), listValues = c(0.5, -0.5))
# Shore-level-specific response: (HighLong - HighNoShort) - (LowLong - LowNoShort).
res_North_shorelevel <- results(dds_North, contrast = list(c("groupHighLong", "groupLowNoShort"), c("groupHighNoShort", "groupLowLong")))
# Heat stress response within each shore level.
res_North_high <- results(dds_North, contrast = c("group", "HighLong", "HighNoShort"))
res_North_low <- results(dds_North, contrast = c("group", "LowLong", "LowNoShort"))

# Save the full results for all genes.
North_generaltemp <- data.frame(gene = rownames(res_North_general), as.data.frame(res_North_general))
North_shorelevel <- data.frame(gene = rownames(res_North_shorelevel), as.data.frame(res_North_shorelevel))
North_hightemp <- data.frame(gene = rownames(res_North_high), as.data.frame(res_North_high))
North_lowtemp <- data.frame(gene = rownames(res_North_low), as.data.frame(res_North_low))
write.table(North_generaltemp, file = "North_generaltemp_allgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(North_shorelevel, file = "North_shorelevel_allgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(North_hightemp, file = "North_hightemp_allgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(North_lowtemp, file = "North_lowtemp_allgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Heat stress response of each gene in high- vs low-shore snails, used for Figure 2A.
highlowtemp_North <- dplyr::left_join(North_hightemp, North_lowtemp, by = "gene")
fchighlow_North <- data.frame(gene_id = highlowtemp_North$gene, highfc = highlowtemp_North$log2FoldChange.x, highpadj = highlowtemp_North$padj.x,
                              lowfc = highlowtemp_North$log2FoldChange.y, lowpadj = highlowtemp_North$padj.y)

# Significant genes: padj < 0.05 and |log2FC| >= 0.5 for both general and shore-level-specific responses. Genes significant in both tests form the intersection group.
general_genes_North <- North_generaltemp$gene[which(North_generaltemp$padj < 0.05 & abs(North_generaltemp$log2FoldChange) >= 0.5)]
shore_genes_North <- North_shorelevel$gene[which(North_shorelevel$padj < 0.05 & abs(North_shorelevel$log2FoldChange) >= 0.5)]
general_only_genes_North <- setdiff(general_genes_North, shore_genes_North)
shore_only_genes_North <- setdiff(shore_genes_North, general_genes_North)
intersection_genes_North <- intersect(general_genes_North, shore_genes_North)

# Label each gene by the category.
fchighlow_North$category <- "Not significant"
fchighlow_North$category[fchighlow_North$gene_id %in% general_only_genes_North] <- "0.5*((HH-HL)+(LH-LL)): general temperature response genes"
fchighlow_North$category[fchighlow_North$gene_id %in% shore_only_genes_North] <- "(HH-HL)-(LH-LL): shore-level-specific response genes"
fchighlow_North$category[fchighlow_North$gene_id %in% intersection_genes_North] <- "Intersection"

# Subsets of plotting data for counting genes in each category in Figure 2A.
# FC and padj columns are from the high- and low-shore comparisons; categories are defined by the general and shore-level-specific tests.
general_only_result_North <- fchighlow_North[fchighlow_North$gene_id %in% general_only_genes_North, ]
shore_only_result_North <- fchighlow_North[fchighlow_North$gene_id %in% shore_only_genes_North, ]
intersection_result_North <- fchighlow_North[fchighlow_North$gene_id %in% intersection_genes_North, ]


#################### Central
# Differential expression among the four groups (shore level x heat stress) in the Central population.
counts_Central <- counts[, group$population == "Central"]
group_Central <- group[group$population == "Central", ]
dds_Central <- DESeqDataSetFromMatrix(countData = counts_Central, colData = group_Central, design = ~ 0 + group)
dds_Central <- DESeq(dds_Central)

# General response: 0.5 * ((HighLong - HighNoShort) + (LowLong - LowNoShort)); mean response across shore levels.
res_Central_general <- results(dds_Central, contrast = list(c("groupHighLong", "groupLowLong"), c("groupHighNoShort", "groupLowNoShort")), listValues = c(0.5, -0.5))
# Shore-level-specific response: (HighLong - HighNoShort) - (LowLong - LowNoShort).
res_Central_shorelevel <- results(dds_Central, contrast = list(c("groupHighLong", "groupLowNoShort"), c("groupHighNoShort", "groupLowLong")))
# Heat stress response within each shore level.
res_Central_high <- results(dds_Central, contrast = c("group", "HighLong", "HighNoShort"))
res_Central_low <- results(dds_Central, contrast = c("group", "LowLong", "LowNoShort"))

# Save the full results for all genes.
Central_generaltemp <- data.frame(gene = rownames(res_Central_general), as.data.frame(res_Central_general))
Central_shorelevel <- data.frame(gene = rownames(res_Central_shorelevel), as.data.frame(res_Central_shorelevel))
Central_hightemp <- data.frame(gene = rownames(res_Central_high), as.data.frame(res_Central_high))
Central_lowtemp <- data.frame(gene = rownames(res_Central_low), as.data.frame(res_Central_low))
write.table(Central_generaltemp, file = "Central_generaltemp_allgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Central_shorelevel, file = "Central_shorelevel_allgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Central_hightemp, file = "Central_hightemp_allgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Central_lowtemp, file = "Central_lowtemp_allgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Heat stress response of each gene in high- vs low-shore snails, used for Figure 2B.
highlowtemp_Central <- dplyr::left_join(Central_hightemp, Central_lowtemp, by = "gene")
fchighlow_Central <- data.frame(gene_id = highlowtemp_Central$gene, highfc = highlowtemp_Central$log2FoldChange.x, highpadj = highlowtemp_Central$padj.x,
                                lowfc = highlowtemp_Central$log2FoldChange.y, lowpadj = highlowtemp_Central$padj.y)

# Significant genes: padj < 0.05 and |log2FC| >= 0.5 for both general and shore-level-specific responses. Genes significant in both tests form the intersection group.
general_genes_Central <- Central_generaltemp$gene[which(Central_generaltemp$padj < 0.05 & abs(Central_generaltemp$log2FoldChange) >= 0.5)]
shore_genes_Central <- Central_shorelevel$gene[which(Central_shorelevel$padj < 0.05 & abs(Central_shorelevel$log2FoldChange) >= 0.5)]
general_only_genes_Central <- setdiff(general_genes_Central, shore_genes_Central)
shore_only_genes_Central <- setdiff(shore_genes_Central, general_genes_Central)
intersection_genes_Central <- intersect(general_genes_Central, shore_genes_Central)

# Label each gene by the category.
fchighlow_Central$category <- "Not significant"
fchighlow_Central$category[fchighlow_Central$gene_id %in% general_only_genes_Central] <- "0.5*((HH-HL)+(LH-LL)): general temperature response genes"
fchighlow_Central$category[fchighlow_Central$gene_id %in% shore_only_genes_Central] <- "(HH-HL)-(LH-LL): shore-level-specific response genes"
fchighlow_Central$category[fchighlow_Central$gene_id %in% intersection_genes_Central] <- "Intersection"

# Subsets of plotting data for counting genes in each category in Figure 2B.
# FC and padj columns are from the high- and low-shore comparisons; categories are defined by the general and shore-level-specific tests.
general_only_result_Central <- fchighlow_Central[fchighlow_Central$gene_id %in% general_only_genes_Central, ]
shore_only_result_Central <- fchighlow_Central[fchighlow_Central$gene_id %in% shore_only_genes_Central, ]
intersection_result_Central <- fchighlow_Central[fchighlow_Central$gene_id %in% intersection_genes_Central, ]


#################### South
# Differential expression among the four groups (shore level x heat stress) in the South population.
counts_South <- counts[, group$population == "South"]
group_South <- group[group$population == "South", ]
dds_South <- DESeqDataSetFromMatrix(countData = counts_South, colData = group_South, design = ~ 0 + group)
dds_South <- DESeq(dds_South)

# General response: 0.5 * ((HighLong - HighNoShort) + (LowLong - LowNoShort)); mean response across shore levels.
res_South_general <- results(dds_South, contrast = list(c("groupHighLong", "groupLowLong"), c("groupHighNoShort", "groupLowNoShort")), listValues = c(0.5, -0.5))
# Shore-level-specific response: (HighLong - HighNoShort) - (LowLong - LowNoShort).
res_South_shorelevel <- results(dds_South, contrast = list(c("groupHighLong", "groupLowNoShort"), c("groupHighNoShort", "groupLowLong")))
# Heat stress response within each shore level.
res_South_high <- results(dds_South, contrast = c("group", "HighLong", "HighNoShort"))
res_South_low <- results(dds_South, contrast = c("group", "LowLong", "LowNoShort"))

# Save the full results for all genes.
South_generaltemp <- data.frame(gene = rownames(res_South_general), as.data.frame(res_South_general))
South_shorelevel <- data.frame(gene = rownames(res_South_shorelevel), as.data.frame(res_South_shorelevel))
South_hightemp <- data.frame(gene = rownames(res_South_high), as.data.frame(res_South_high))
South_lowtemp <- data.frame(gene = rownames(res_South_low), as.data.frame(res_South_low))
write.table(South_generaltemp, file = "South_generaltemp_allgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(South_shorelevel, file = "South_shorelevel_allgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(South_hightemp, file = "South_hightemp_allgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(South_lowtemp, file = "South_lowtemp_allgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Heat stress response of each gene in high- vs low-shore snails, used for Figure 2C.
highlowtemp_South <- dplyr::left_join(South_hightemp, South_lowtemp, by = "gene")
fchighlow_South <- data.frame(gene_id = highlowtemp_South$gene, highfc = highlowtemp_South$log2FoldChange.x, highpadj = highlowtemp_South$padj.x,
                              lowfc = highlowtemp_South$log2FoldChange.y, lowpadj = highlowtemp_South$padj.y)

# Significant genes: padj < 0.05 and |log2FC| >= 0.5 for both general and shore-level-specific responses. Genes significant in both tests form the intersection group.
general_genes_South <- South_generaltemp$gene[which(South_generaltemp$padj < 0.05 & abs(South_generaltemp$log2FoldChange) >= 0.5)]
shore_genes_South <- South_shorelevel$gene[which(South_shorelevel$padj < 0.05 & abs(South_shorelevel$log2FoldChange) >= 0.5)]
general_only_genes_South <- setdiff(general_genes_South, shore_genes_South)
shore_only_genes_South <- setdiff(shore_genes_South, general_genes_South)
intersection_genes_South <- intersect(general_genes_South, shore_genes_South)

# Label each gene by the category.
fchighlow_South$category <- "Not significant"
fchighlow_South$category[fchighlow_South$gene_id %in% general_only_genes_South] <- "0.5*((HH-HL)+(LH-LL)): general temperature response genes"
fchighlow_South$category[fchighlow_South$gene_id %in% shore_only_genes_South] <- "(HH-HL)-(LH-LL): shore-level-specific response genes"
fchighlow_South$category[fchighlow_South$gene_id %in% intersection_genes_South] <- "Intersection"

# Subsets of plotting data for counting genes in each category in Figure 2C.
# FC and padj columns are from the high- and low-shore comparisons; categories are defined by the general and shore-level-specific tests.
general_only_result_South <- fchighlow_South[fchighlow_South$gene_id %in% general_only_genes_South, ]
shore_only_result_South <- fchighlow_South[fchighlow_South$gene_id %in% shore_only_genes_South, ]
intersection_result_South <- fchighlow_South[fchighlow_South$gene_id %in% intersection_genes_South, ]


#################### Figure 2A-C

# Panel A: North.
p2A <- ggplot(fchighlow_North, aes(x = lowfc, y = highfc)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, linetype = "dashed") + geom_vline(xintercept = 0, colour = "black", linewidth = 0.4, linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.4, linetype = "dashed") +
  geom_point(data = fchighlow_North[fchighlow_North$category == "Not significant", ], aes(colour = category), size = 1, na.rm = TRUE) +
  geom_point(data = fchighlow_North[fchighlow_North$category != "Not significant", ], aes(colour = category), size = 1, na.rm = TRUE) +
  scale_colour_manual(values = c("0.5*((HH-HL)+(LH-LL)): general temperature response genes" = "#029F73", "(HH-HL)-(LH-LL): shore-level-specific response genes" = "#ff7f00",
                                 "Intersection" = "#631879", "Not significant" = "#D3D3D3"),
                      breaks = c("0.5*((HH-HL)+(LH-LL)): general temperature response genes", "(HH-HL)-(LH-LL): shore-level-specific response genes", "Intersection", "Not significant"),
                      labels = c("General temperature response genes", "Shore-level-specific response genes", "Intersection", "Not significant"), drop = FALSE,
                      guide = guide_legend(nrow = 1, override.aes = list(size = 2.5))) +
  scale_x_continuous(breaks = seq(-20, 10, by = 10)) + scale_y_continuous(breaks = seq(-30, 10, by = 10)) + coord_cartesian(xlim = c(-29, 11), ylim = c(-30, 10), expand = FALSE) +
  labs(x = "", y = expression(Log[2]*FC~"(" * italic(High[long]) - italic(High["no/short"]) * ")"), colour = "Category", title = "North") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.justification = "center", legend.title = element_text(size = 13), legend.text = element_text(size = 12),
        legend.key.spacing.x = unit(14, "pt"), legend.box.margin = margin(6, 0, 0, 0), axis.title.y = element_text(size = 13, margin = margin(r = 10)),
        axis.title.x = element_text(size = 13, margin = margin(t = 20)), axis.text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 14),
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.25), panel.grid.minor = element_line(colour = "grey95", linewidth = 0.2))

# Panel B: Central.
p2B <- ggplot(fchighlow_Central, aes(x = lowfc, y = highfc)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, linetype = "dashed") + geom_vline(xintercept = 0, colour = "black", linewidth = 0.4, linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.4, linetype = "dashed") +
  geom_point(data = fchighlow_Central[fchighlow_Central$category == "Not significant", ], aes(colour = category), size = 1, na.rm = TRUE) +
  geom_point(data = fchighlow_Central[fchighlow_Central$category != "Not significant", ], aes(colour = category), size = 1, na.rm = TRUE) +
  scale_colour_manual(values = c("0.5*((HH-HL)+(LH-LL)): general temperature response genes" = "#029F73", "(HH-HL)-(LH-LL): shore-level-specific response genes" = "#ff7f00",
                                 "Intersection" = "#631879", "Not significant" = "#D3D3D3"),
                      breaks = c("0.5*((HH-HL)+(LH-LL)): general temperature response genes", "(HH-HL)-(LH-LL): shore-level-specific response genes", "Intersection", "Not significant"),
                      labels = c("General temperature response genes", "Shore-level-specific response genes", "Intersection", "Not significant"), drop = FALSE,
                      guide = guide_legend(nrow = 1, override.aes = list(size = 2.5))) +
  scale_x_continuous(breaks = seq(-20, 10, by = 10)) + scale_y_continuous(breaks = seq(-30, 10, by = 10)) + coord_cartesian(xlim = c(-29, 11), ylim = c(-30, 10), expand = FALSE) +
  labs(x = expression(Log[2]*FC~"(" * italic(Low[long]) - italic(Low["no/short"]) * ")"), y = "", colour = "Category", title = "Central") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.justification = "center", legend.title = element_text(size = 13), legend.text = element_text(size = 12),
        legend.key.spacing.x = unit(14, "pt"), legend.box.margin = margin(6, 0, 0, 0), axis.title.y = element_text(size = 13, margin = margin(r = 10)),
        axis.title.x = element_text(size = 13, margin = margin(t = 20)), axis.text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 14),
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.25), panel.grid.minor = element_line(colour = "grey95", linewidth = 0.2))

# Panel C: South.
p2C <- ggplot(fchighlow_South, aes(x = lowfc, y = highfc)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, linetype = "dashed") + geom_vline(xintercept = 0, colour = "black", linewidth = 0.4, linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.4, linetype = "dashed") +
  geom_point(data = fchighlow_South[fchighlow_South$category == "Not significant", ], aes(colour = category), size = 1, na.rm = TRUE) +
  geom_point(data = fchighlow_South[fchighlow_South$category != "Not significant", ], aes(colour = category), size = 1, na.rm = TRUE) +
  scale_colour_manual(values = c("0.5*((HH-HL)+(LH-LL)): general temperature response genes" = "#029F73", "(HH-HL)-(LH-LL): shore-level-specific response genes" = "#ff7f00",
                                 "Intersection" = "#631879", "Not significant" = "#D3D3D3"),
                      breaks = c("0.5*((HH-HL)+(LH-LL)): general temperature response genes", "(HH-HL)-(LH-LL): shore-level-specific response genes", "Intersection", "Not significant"),
                      labels = c("General temperature response genes", "Shore-level-specific response genes", "Intersection", "Not significant"), drop = FALSE,
                      guide = guide_legend(nrow = 1, override.aes = list(size = 2.5))) +
  scale_x_continuous(breaks = seq(-20, 10, by = 10)) + scale_y_continuous(breaks = seq(-30, 10, by = 10)) + coord_cartesian(xlim = c(-29, 11), ylim = c(-30, 10), expand = FALSE) +
  labs(x = "", y = "", colour = "Category", title = "South") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.justification = "center", legend.title = element_text(size = 13), legend.text = element_text(size = 12),
        legend.key.spacing.x = unit(14, "pt"), legend.box.margin = margin(6, 0, 0, 0), axis.title.y = element_text(size = 13, margin = margin(r = 10)),
        axis.title.x = element_text(size = 13, margin = margin(t = 20)), axis.text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 14),
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.25), panel.grid.minor = element_line(colour = "grey95", linewidth = 0.2))


#################### Figure 2D: number of genes in each category

category_counts_df <- data.frame(Category = c("0.5*((HH−HL)+(LH−LL))", "(HH−HL)−(LH−LL)", "Intersection"),
                                 North = c(nrow(general_only_result_North), nrow(shore_only_result_North), nrow(intersection_result_North)),
                                 Central = c(nrow(general_only_result_Central), nrow(shore_only_result_Central), nrow(intersection_result_Central)),
                                 South = c(nrow(general_only_result_South), nrow(shore_only_result_South), nrow(intersection_result_South)))
category_counts_long <- tidyr::pivot_longer(category_counts_df, cols = -Category, names_to = "Population", values_to = "GeneNumber")
category_counts_long$Category <- factor(category_counts_long$Category, levels = c("0.5*((HH−HL)+(LH−LL))", "(HH−HL)−(LH−LL)", "Intersection"))
category_counts_long$Population <- factor(category_counts_long$Population, levels = c("North", "Central", "South"))

p2D <- ggplot(category_counts_long, aes(x = Population, y = GeneNumber, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = GeneNumber), position = position_dodge(width = 0.9), vjust = -0.3, size = 12, size.unit = "pt") +
  scale_fill_manual(values = c("0.5*((HH−HL)+(LH−LL))" = "#029F73", "(HH−HL)−(LH−LL)" = "#ff7f00", "Intersection" = "#631879")) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(title = "", x = "", y = "Number of significant genes", fill = "Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 13, colour = "black"), axis.text.y = element_text(size = 12), axis.title = element_text(size = 13), legend.title = element_text(size = 11),
        legend.text = element_text(size = 9), legend.position = "none",
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.3), panel.grid.minor = element_line(colour = "grey95", linewidth = 0.2))


#################### Figure 2E: overlap of general temperature response genes among populations

# How many general temperature response genes are unique to one population, shared by two, or shared by all three.
general_union <- union(union(general_genes_North, general_genes_Central), general_genes_South)
venn_counts <- c("North" = length(setdiff(general_genes_North, union(general_genes_Central, general_genes_South))),
                 "Central" = length(setdiff(general_genes_Central, union(general_genes_North, general_genes_South))),
                 "South" = length(setdiff(general_genes_South, union(general_genes_North, general_genes_Central))),
                 "North&Central" = length(setdiff(intersect(general_genes_North, general_genes_Central), general_genes_South)),
                 "North&South" = length(setdiff(intersect(general_genes_North, general_genes_South), general_genes_Central)),
                 "Central&South" = length(setdiff(intersect(general_genes_Central, general_genes_South), general_genes_North)),
                 "North&Central&South" = length(intersect(intersect(general_genes_North, general_genes_Central), general_genes_South)))
venn_percent <- round(100 * venn_counts / length(general_union))
# Count above percentage (as in Figure 3B), so that the labels fit the small overlap regions.
labels_with_percent <- paste0(venn_counts, "\n(", venn_percent, "%)")
names(labels_with_percent) <- names(venn_counts)
venn_region_table <- data.frame(region = names(venn_counts), n = as.integer(venn_counts), percent = as.numeric(venn_percent), label = labels_with_percent)

# Color of Euler diagram
euler_fit <- euler(venn_counts)
region_colours <- c("North" = "#ddeed9", "Central" = "#ddeed9", "South" = "#ddeed9", "North&Central" = "#A8E6A1", "North&South" = "#A8E6A1",
                    "Central&South" = "#A8E6A1", "North&Central&South" = "#1b9e77")
euler_region_order <- names(euler_fit$original.values)
euler_diagram <- plot(euler_fit, fills = list(fill = region_colours[euler_region_order]), labels = FALSE,
                      quantities = list(labels = labels_with_percent[euler_region_order], cex = 12 / 12, font = 1, col = "#333333"), edges = list(col = "black", lwd = 2))

# Add population names to the diagram.
euler_grob <- grid.grabExpr({
  print(euler_diagram)
  grid.text("North", x = unit(0.12, "npc"), y = unit(0.73, "npc"), gp = gpar(fontsize = 13))
  grid.text("South", x = unit(0.59, "npc"), y = unit(0.975, "npc"), gp = gpar(fontsize = 13))
  grid.text("Central", x = unit(0.84, "npc"), y = unit(0.025, "npc"), gp = gpar(fontsize = 13))
}, wrap = TRUE, width = 3.5, height = 3.5)
p2E <- ggdraw() + draw_grob(euler_grob, x = 0.04, y = 0.04, width = 0.92, height = 0.92)


#################### Figure 2F: response magnitude of general response genes shared by all three populations

# For general response genes shared by all three populations, compare how strongly they respond in each population.
common_genes <- intersect(intersect(general_genes_North, general_genes_Central), general_genes_South)
North_common <- North_generaltemp[match(common_genes, North_generaltemp$gene), ]
Central_common <- Central_generaltemp[match(common_genes, Central_generaltemp$gene), ]
South_common <- South_generaltemp[match(common_genes, South_generaltemp$gene), ]

# Direction of the response of the shared general response genes in the three populations.
direction_common <- data.frame(gene = common_genes, North = sign(North_common$log2FoldChange),
                               Central = sign(Central_common$log2FoldChange), South = sign(South_common$log2FoldChange))
same_direction <- direction_common$North == direction_common$Central & direction_common$North == direction_common$South
sum(same_direction)                                  # shared genes changing in the same direction in all three populations
sum(same_direction & direction_common$North == 1)    # up-regulated in all three populations
sum(same_direction & direction_common$North == -1)   # down-regulated in all three populations

log2FC_df <- data.frame(gene = rep(common_genes, 3), Population = rep(c("North", "Central", "South"), each = length(common_genes)),
                        log2FC = c(North_common$log2FoldChange, Central_common$log2FoldChange, South_common$log2FoldChange))
log2FC_df$Population <- factor(log2FC_df$Population, levels = c("North", "Central", "South"))
log2FC_df$abs_log2FC <- abs(log2FC_df$log2FC)

# Mean and median |log2FC| per population.
mean_response <- aggregate(abs_log2FC ~ Population, data = log2FC_df, FUN = mean)
median_response <- aggregate(abs_log2FC ~ Population, data = log2FC_df, FUN = median)
log2fc_data <- data.frame(Population = c("North", "Central", "South"),
                          Mean = mean_response$abs_log2FC[match(c("North", "Central", "South"), mean_response$Population)],
                          Median = median_response$abs_log2FC[match(c("North", "Central", "South"), median_response$Population)])

# Labels for the two lines ("Mean", "Median").
label_data <- data.frame(Population = c("Central", "Central"), Statistic = c("Mean", "Median"),
                         log2FC = c(log2fc_data$Mean[log2fc_data$Population == "Central"], log2fc_data$Median[log2fc_data$Population == "Central"]), label = c("Mean", "Median"))
label_data$Population <- factor(label_data$Population, levels = c("North", "Central", "South"))
log2fc_data$Population <- factor(log2fc_data$Population, levels = c("North", "Central", "South"))
log2fc_long <- tidyr::pivot_longer(log2fc_data, cols = c(Mean, Median), names_to = "Statistic", values_to = "log2FC")
log2fc_long$value_label <- format(round(log2fc_long$log2FC, 3), nsmall = 3)

p2F <- ggplot(log2fc_long, aes(x = Population, y = log2FC, group = Statistic)) +
  geom_line(aes(linetype = Statistic, colour = Statistic), linewidth = 1) + geom_point(aes(shape = Statistic), colour = "black", size = 3) +
  geom_text(aes(label = value_label), hjust = 0, nudge_x = 0.04, vjust = -1.2, size = 12, size.unit = "pt", colour = "black") +
  geom_text(data = label_data[label_data$Statistic == "Mean", ], aes(label = label), colour = "black", vjust = 0.5, nudge_x = -0.25, nudge_y = -0.08, size = 12, size.unit = "pt", fontface = "italic") +
  geom_text(data = label_data[label_data$Statistic == "Median", ], aes(label = label), colour = "black", vjust = 0.5, nudge_x = -0.25, nudge_y = -0.08, size = 12, size.unit = "pt", fontface = "italic") +
  scale_colour_manual(values = c("Mean" = "#029F73", "Median" = "#029F73")) + scale_linetype_manual(values = c("Mean" = "dashed", "Median" = "solid")) +
  scale_shape_manual(values = c("Mean" = 15, "Median" = 16)) + scale_x_discrete(expand = c(0.2, 0)) + coord_cartesian(clip = "off") +
  labs(y = expression("|log"[2]*plain("FC")*"|"), x = "", title = "") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(size = 13, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"), axis.title.y = element_text(size = 13), legend.position = "none",
        plot.title = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 20),
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.3), panel.grid.minor = element_line(colour = "grey95", linewidth = 0.2))


#################### Combine Figure 2A-F

# Use the panel A legend as the shared legend.
p2A_grob <- ggplotGrob(p2A)
legend <- p2A_grob$grobs[[which(p2A_grob$layout$name == "guide-box-bottom")]]

p_figure2 <- ggarrange(p2A, p2B, p2C, p2D, p2E, p2F, labels = c("A", "B", "C", "D", "E", "F"), ncol = 3, nrow = 2,
                       common.legend = TRUE, legend = "bottom", legend.grob = legend, font.label = list(size = 17, face = "bold"))
p_figure2
ggsave("Figure2.pdf", p_figure2, width = 12, height = 8, device = grDevices::pdf, useDingbats = FALSE)



#################### Overlap test of general response genes among populations (SuperExactTest)

# Test whether the three populations share more general response genes than expected by chance.
temperature_gene_sets <- list(North = general_genes_North, South = general_genes_South, Central = general_genes_Central)
superexact_result <- supertest(temperature_gene_sets, n = nrow(counts))
summary(superexact_result)


#################### Shore-level-specific overlap among populations

# All significant shore-level genes (orange + purple), orange genes only (Table S3) and purple genes only (Table S4).
shore_gene_sets <- list(North = shore_genes_North, Central = shore_genes_Central, South = shore_genes_South)
shore_only_gene_sets <- list(North = shore_only_genes_North, Central = shore_only_genes_Central, South = shore_only_genes_South)
intersection_gene_sets <- list(North = intersection_genes_North, Central = intersection_genes_Central, South = intersection_genes_South)

# Genes shared by each pair of populations (including those also shared by the third) and by all three.
shared_shore_genes <- function(gene_sets) {
  list("North&Central" = sort(intersect(gene_sets$North, gene_sets$Central)),
       "North&South" = sort(intersect(gene_sets$North, gene_sets$South)),
       "Central&South" = sort(intersect(gene_sets$Central, gene_sets$South)),
       "North&Central&South" = sort(Reduce(intersect, gene_sets)))
}
shore_shared_genes <- shared_shore_genes(shore_gene_sets)
lengths(shore_shared_genes)
shore_shared_genes
shore_only_shared_genes <- shared_shore_genes(shore_only_gene_sets)
lengths(shore_only_shared_genes)
shore_only_shared_genes

# Genes in each of the seven non-overlapping regions (unique to one population, shared by two only, shared by all three).
overlap_regions <- function(gene_sets) {
  combinations <- list(North = "North", Central = "Central", South = "South",
                       "North&Central" = c("North", "Central"), "North&South" = c("North", "South"),
                       "Central&South" = c("Central", "South"), "North&Central&South" = names(gene_sets))
  lapply(combinations, function(included) {
    excluded <- setdiff(names(gene_sets), included)
    shared <- Reduce(intersect, gene_sets[included])
    others <- Reduce(union, gene_sets[excluded], init = character())
    sort(setdiff(shared, others))
  })
}
lengths(overlap_regions(shore_gene_sets))
lengths(overlap_regions(shore_only_gene_sets))
shore_superexact_result <- supertest(shore_gene_sets, n = nrow(counts))
summary(shore_superexact_result)


#################### Linkage group (LG) enrichment of general response genes

# Number of retained genes per LG.
pos <- read.table("LG_gene_position_20176.txt", header = TRUE)
pos <- pos[pos$gene_id %in% rownames(counts), ]
gene_lg <- dplyr::summarise(dplyr::group_by(pos, chr), genenum = dplyr::n(), .groups = "drop")
gene_lg <- dplyr::filter(gene_lg, grepl("^LG", chr))

# Test whether any LG carries more or fewer genes of a gene set than expected from its total number of genes
# (two-sided Fisher's exact test, BH correction across LGs).
test_lg <- function(test_genes) {
  gene_lg2 <- dplyr::summarise(dplyr::group_by(pos[pos$gene_id %in% test_genes, ], chr), genenum = dplyr::n(), .groups = "drop")
  genechr <- dplyr::left_join(gene_lg, gene_lg2, by = "chr")
  colnames(genechr) <- c("chr", "allgenenum", "DEnum")
  genechr$DEnum[is.na(genechr$DEnum)] <- 0
  fisher_test <- do.call(rbind, lapply(seq_len(nrow(genechr)), function(k) {
    a <- genechr$DEnum[k]                                   # genes of the set on this LG
    b <- genechr$allgenenum[k] - a                          # other genes on this LG
    c <- length(test_genes) - a                             # genes of the set elsewhere
    d <- nrow(counts) - genechr$allgenenum[k] - c           # other genes elsewhere
    data.frame(chr = genechr$chr[k], DEnum = a, allgenenum = genechr$allgenenum[k],
               pvalue = fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "two.sided")$p.value)
  }))
  fisher_test$padj <- p.adjust(fisher_test$pvalue, method = "BH")
  fisher_test
}

# General response genes shared by all three populations, and those of each population.
lg_shared <- test_lg(common_genes)
lg_North <- test_lg(general_genes_North)
lg_Central <- test_lg(general_genes_Central)
lg_South <- test_lg(general_genes_South)
sapply(list(shared = lg_shared, North = lg_North, Central = lg_Central, South = lg_South), function(x) min(x$padj))


