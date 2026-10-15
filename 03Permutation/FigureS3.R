# =============================================================================
# Permutation test of the number of differentially expressed genes (DEGs) between treatments (Figure S3)
#
# For each pair of the four treatments (17°C 15 min, 17°C 3 h, 31°C 15 min, 31°C 3 h), the samples of the two
# treatments are pooled within each population x shore level, and the treatment labels are reassigned in every
# possible way (21 relabellings for 2 vs 5 samples, 126 for 5 vs 5 samples). For each relabelling, DESeq2 is refitted
# on all samples and the DEGs (padj < 0.05, |log2FC| >= 1) between the two treatments are counted in each
# population x shore level. The observed number of DEGs is then compared with this permutation distribution.
#
# Part 1 runs the permutations. Part 2 draws Figure S3 from the saved results.
#
# Input:  counts_filter.txt, group.txt (from 01counts_ind_filter)
# Output: DEG_1715vs173.txt, DEG_1715vs3115.txt, DEG_1715vs313.txt, DEG_173vs3115.txt, DEG_173vs313.txt, DEG_3115vs313.txt, FigureS3.pdf
# =============================================================================

setwd("/data/03Permutation/")
library(DESeq2)
library(ggplot2)
library(ggpubr)


############################################################
#Part 1: permutations
# Running all permutations takes many hours (one DESeq2 fit per relabelling); the saved DEG_*.txt files are read in Part 2.

# Counts and sample information; condition = population x shore level x treatment (24 groups).
counts <- read.table("../01counts_ind_filter/counts_filter.txt", header = TRUE)
group <- read.table("../01counts_ind_filter/group.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
group$treatment <- paste(group$temperature, group$time, sep = "_")
group$popshore <- paste(group$population, group$shore_level, sep = "_")
group$condition <- paste(group$popshore, group$treatment, sep = "_")
popshores <- c("North_high", "North_low", "Central_high", "Central_low", "South_high", "South_low")

# Relabellings of the pooled samples (one column each). Within each population x shore level the pooled samples are in
# file order (first treatment, then second treatment); order[k] is the sample whose label goes to position k.
# 2 vs 5 samples: every choice of the 2 positions of the smaller treatment (21). 5 vs 5 samples: every choice of the
# 5 positions of the first treatment, with position 1 fixed so that each split is counted once (126).
orders_2first <- combn(7, 2, function(x) { o <- rep(NA, 7); o[x] <- 1:2; o[is.na(o)] <- 3:7; o })
orders_2last <- combn(7, 2, function(x) { o <- rep(NA, 7); o[x] <- 6:7; o[is.na(o)] <- 1:5; o })
combs_5of10 <- combn(10, 5)
orders_5of10 <- apply(combs_5of10[, combs_5of10[1, ] == 1], 2, function(x) { o <- rep(NA, 10); o[x] <- 1:5; o[is.na(o)] <- 6:10; o })

# One relabelling: reassign the two treatment labels within every population x shore level, refit DESeq2 on all samples
# and count the DEGs between the two treatments in each population x shore level.
# South low has one 31°C 3 h sample fewer (SPL47 removed); there the relabelling is shortened as order[order != 1] - 1.
count_DEGs <- function(trt1, trt2, order) {
  grp <- group
  for (ps in popshores) {
    idx <- which(grp$popshore == ps & grp$treatment %in% c(trt1, trt2))
    ord <- if (length(idx) == length(order)) order else order[order != 1] - 1
    grp$condition[idx] <- grp$condition[idx][ord]
  }
  rownames(grp) <- grp$id
  grp <- grp[colnames(counts), ]
  dds <- DESeq(DESeqDataSetFromMatrix(countData = counts, colData = grp, design = ~ condition))
  sapply(popshores, function(ps) {
    res <- results(dds, contrast = c("condition", paste(ps, trt1, sep = "_"), paste(ps, trt2, sep = "_")), alpha = 0.05)
    sum(res$padj < 0.05 & abs(res$log2FoldChange) >= 1, na.rm = TRUE)
  })
}

# All relabellings for one pair of treatments; the observed labelling is moved to the first row before saving.
run_permutation <- function(trt1, trt2, orders, file) {
  degs <- t(apply(orders, 2, function(o) count_DEGs(trt1, trt2, o)))
  observed <- apply(orders, 2, function(o) all(o == seq_along(o)))
  degs <- rbind(degs[observed, , drop = FALSE], degs[!observed, , drop = FALSE])
  write.table(degs, file = file, quote = F, row.names = F, col.names = F)
  degs
}

perm_1715vs173 <- run_permutation("17degree_15min", "17degree_3h", orders_2last, "DEG_1715vs173.txt")
perm_1715vs3115 <- run_permutation("17degree_15min", "31degree_15min", orders_5of10, "DEG_1715vs3115.txt")
perm_1715vs313 <- run_permutation("17degree_15min", "31degree_3h", orders_5of10, "DEG_1715vs313.txt")
perm_173vs3115 <- run_permutation("17degree_3h", "31degree_15min", orders_2first, "DEG_173vs3115.txt")
perm_173vs313 <- run_permutation("17degree_3h", "31degree_3h", orders_2first, "DEG_173vs313.txt")
perm_3115vs313 <- run_permutation("31degree_15min", "31degree_3h", orders_5of10, "DEG_3115vs313.txt")


############################################################
#Part 2: Figure S3

# Permutation results.
d_1715vs173 <- read.table("DEG_1715vs173.txt", header = FALSE)
d_1715vs3115 <- read.table("DEG_1715vs3115.txt", header = FALSE)
d_1715vs313 <- read.table("DEG_1715vs313.txt", header = FALSE)
d_173vs3115 <- read.table("DEG_173vs3115.txt", header = FALSE)
d_173vs313 <- read.table("DEG_173vs313.txt", header = FALSE)
d_3115vs313 <- read.table("DEG_3115vs313.txt", header = FALSE)

x_labels <- c("North_high", "North_low", "Central_high", "Central_low", "South_high", "South_low")
group_colors <- c("North_high" = "#2a4fa2", "North_low" = "#a0b4e6", "Central_high" = "#c9a800", "Central_low" = "#fff2a6", "South_high" = "#b02a2a", "South_low" = "#ffadad")

# Permutation p-value: how often a random relabelling gives at least as many DEGs as the observed labelling
# (the observed labelling is counted as one of the possible labellings), shown as stars.
get_p_stars <- function(obs, null) {
  p <- (sum(null >= obs) + 1) / (length(null) + 1)
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
}

# Violin plot of the permutation distribution per population x shore level (all labellings, including the observed one,
# as in the p-value): black points = relabellings, red point (centred) = observed, stars = permutation p-value.
plot_violin <- function(data, title, ylabel = FALSE) {
  # The first row of each file is the observed labelling ("obs", red point); all other rows are relabellings ("perm").
  df <- data.frame(value = unlist(data), group = factor(rep(x_labels, each = nrow(data)), levels = x_labels),
                   type = rep(c("obs", rep("perm", nrow(data) - 1)), times = ncol(data)))
  stars <- sapply(1:ncol(data), function(i) get_p_stars(data[1, i], data[-1, i]))
  ggplot() + geom_violin(data = df, aes(x = group, y = value, fill = group), trim = FALSE) +
    geom_jitter(data = df[df$type == "perm", ], aes(x = group, y = value), color = "black", size = 0.2, width = 0.2, height = 0) +
    geom_point(data = df[df$type == "obs", ], aes(x = group, y = value), color = "red", size = 1.5) +
    geom_text(data = data.frame(x = x_labels, y = apply(data, 2, max) + 100, label = stars), aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 10) +
    scale_fill_manual(values = group_colors) + labs(title = title, x = NULL, y = if (ylabel) "Number of DEGs" else NULL) + theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5))
}

# Fixed seed so the points are spread horizontally the same way every time (no vertical jitter, so heights are the true DEG numbers).
set.seed(123)
pS3A <- plot_violin(d_1715vs173, "17°C15min vs 17°C3h", ylabel = TRUE)
pS3B <- plot_violin(d_1715vs3115, "17°C15min vs 31°C15min")
pS3C <- plot_violin(d_1715vs313, "17°C15min vs 31°C3h")
pS3D <- plot_violin(d_173vs3115, "17°C3h vs 31°C15min", ylabel = TRUE)
pS3E <- plot_violin(d_173vs313, "17°C3h vs 31°C3h")
pS3F <- plot_violin(d_3115vs313, "31°C15min vs 31°C3h")

figS3 <- ggarrange(pS3A, pS3B, pS3C, pS3D, pS3E, pS3F, labels = c("A", "B", "C", "D", "E", "F"), ncol = 3, nrow = 2, common.legend = FALSE, font.label = list(size = 14, face = "bold"))
figS3

ggsave("FigureS3.pdf", plot = figS3, width = 12, height = 10)
ggsave("FigureS3.png", plot = figS3, width = 12, height = 10)

