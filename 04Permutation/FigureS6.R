setwd("C:/path/to/RNAgithub/04Permutation/")

# =============================================================================
# Figure S6: permutation test of the number of DEGs between treatments
# For each pair of treatments, the treatment labels are shuffled within each population x shore level in every possible way
# (21 permutations for 2 vs 5 samples, 126 for 5 vs 5), DESeq2 is refitted on all samples and the DEGs (padj < 0.05, |log2FC| >= 0.5)
# are counted in each population x shore level. The observed number of DEGs is compared with these permutations.
# Input: counts_filter.txt, group.txt, DEGs_*.txt and DEGs_*_order.txt
# Output: FigureS6.pdf
# =============================================================================

library(DESeq2)
library(ggplot2)
library(ggpubr)


############################################################
#Part 1: permutations

# Counts and sample information
counts <- read.table("../01counts_ind_filter/counts_filter.txt", header = TRUE)
group <- read.table("../01counts_ind_filter/group.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
group$treatment <- paste(group$temperature, group$time, sep = "_")
group$popshore <- paste(group$population, group$shore_level, sep = "_")
group$condition <- paste(group$popshore, group$treatment, sep = "_")
popshores <- c("North_high", "North_low", "Central_high", "Central_low", "South_high", "South_low")

# Permutations: all ways to swap the treatment labels among the pooled samples of one population x shore level (one column per permutation)
# 2 vs 5 samples, smaller treatment first or last in the file: 21 permutations each
orders_2first <- combn(7, 2, function(x) { o <- rep(NA, 7); o[x] <- 1:2; o[is.na(o)] <- 3:7; o })
orders_2last <- combn(7, 2, function(x) { o <- rep(NA, 7); o[x] <- 6:7; o[is.na(o)] <- 1:5; o })
# 5 vs 5 samples: 126 permutations (sample 1 always stays in the first group, so each split is counted only once)
combs_5of10 <- combn(10, 5)
orders_5of10 <- apply(combs_5of10[, combs_5of10[1, ] == 1], 2, function(x) { o <- rep(NA, 10); o[x] <- 1:5; o[is.na(o)] <- 6:10; o })

# One permutation: shuffle the two treatment labels within each population x shore level, refit DESeq2 and count the DEGs in each population x shore level
# South low has one 31°C 3 h sample fewer (SPL47 removed)
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
    sum(res$padj < 0.05 & abs(res$log2FoldChange) >= 0.5, na.rm = TRUE)
  })
}

# Run all permutations of one treatment comparison and save the DEG counts and the permutation orders
run_permutation <- function(trt1, trt2, orders, file) {
  degs <- t(apply(orders, 2, function(o) count_DEGs(trt1, trt2, o)))
  write.table(as.vector(t(degs)), file = file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  order_file <- sub("\\.txt$", "_order.txt", file)
  write.table(t(orders), file = order_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  degs
}

#perm_1715vs173 <- run_permutation("17degree_15min", "17degree_3h", orders_2last, "DEGs_1715vs173.txt")
#perm_1715vs3115 <- run_permutation("17degree_15min", "31degree_15min", orders_5of10, "DEGs_1715vs3115.txt")
#perm_1715vs313 <- run_permutation("17degree_15min", "31degree_3h", orders_5of10, "DEGs_1715vs313.txt")
#perm_173vs3115 <- run_permutation("17degree_3h", "31degree_15min", orders_2first, "DEGs_173vs3115.txt")
#perm_173vs313 <- run_permutation("17degree_3h", "31degree_3h", orders_2first, "DEGs_173vs313.txt")
#perm_3115vs313 <- run_permutation("31degree_15min", "31degree_3h", orders_5of10, "DEGs_3115vs313.txt")


############################################################
# Part 2: Figure S6

# Read the DEG counts of all permutations of one comparison (one row per permutation, one column per population x shore level),
# with the observed labelling moved to the first row
read_permutation <- function(comparison) {
  degs <- matrix(scan(paste0("DEGs_", comparison, ".txt"), quiet = TRUE), ncol = length(popshores), byrow = TRUE)
  orders <- as.matrix(read.table(paste0("DEGs_", comparison, "_order.txt")))
  obs <- which(apply(orders, 1, function(o) all(o == seq_along(o))))
  data <- as.data.frame(degs[c(obs, setdiff(seq_len(nrow(degs)), obs)), ])
  names(data) <- popshores
  data
}
d_1715vs173 <- read_permutation("1715vs173")
d_1715vs3115 <- read_permutation("1715vs3115")
d_1715vs313 <- read_permutation("1715vs313")
d_173vs3115 <- read_permutation("173vs3115")
d_173vs313 <- read_permutation("173vs313")
d_3115vs313 <- read_permutation("3115vs313")

# Permutation p-value: proportion of all labellings (observed included) with at least as many DEGs as observed
get_p <- function(x) mean(x >= x[1])
get_stars <- function(p) if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""

# Observed DEG counts and permutation p-values of all comparisons
d_all <- list("1715vs173" = d_1715vs173, "1715vs3115" = d_1715vs3115, "1715vs313" = d_1715vs313, "173vs3115" = d_173vs3115, "173vs313" = d_173vs313, "3115vs313" = d_3115vs313)
sapply(d_all, function(d) unlist(d[1, ]))
sapply(d_all, function(d) sapply(d, get_p))

# Violin plot of the permutation distribution: black points = permutations, red point = observed, stars = permutation p-value
x_labels <- c("North_high", "North_low", "Central_high", "Central_low", "South_high", "South_low")
group_colors <- c("North_high" = "#2a4fa2", "North_low" = "#a0b4e6", "Central_high" = "#c9a800", "Central_low" = "#fff2a6", "South_high" = "#b02a2a", "South_low" = "#ffadad")
plot_violin <- function(data, title, ylabel = FALSE) {
  df <- data.frame(value = unlist(data), group = factor(rep(x_labels, each = nrow(data)), levels = x_labels), type = rep(c("obs", rep("perm", nrow(data) - 1)), ncol(data)))
  stars <- sapply(data, function(x) get_stars(get_p(x)))
  ggplot() + geom_violin(data = df, aes(x = group, y = value, fill = group), trim = FALSE) +
    geom_jitter(data = df[df$type == "perm", ], aes(x = group, y = value), color = "black", size = 0.2, width = 0.2, height = 0) +
    geom_point(data = df[df$type == "obs", ], aes(x = group, y = value), color = "red", size = 1.5) +
    geom_text(data = data.frame(x = x_labels, y = sapply(data, max) + 100, label = stars), aes(x = x, y = y, label = label), size = 14, size.unit = "pt") +
    scale_fill_manual(values = group_colors) + labs(title = title, x = NULL, y = if (ylabel) "Number of DEGs" else NULL) + theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 12), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14))
}

# Fixed seed so the jittered points look the same every time
set.seed(123)
pS6A <- plot_violin(d_1715vs173, "17°C15min vs 17°C3h", ylabel = TRUE)
pS6B <- plot_violin(d_1715vs3115, "17°C15min vs 31°C15min")
pS6C <- plot_violin(d_1715vs313, "17°C15min vs 31°C3h")
pS6D <- plot_violin(d_173vs3115, "17°C3h vs 31°C15min", ylabel = TRUE)
pS6E <- plot_violin(d_173vs313, "17°C3h vs 31°C3h")
pS6F <- plot_violin(d_3115vs313, "31°C15min vs 31°C3h")

figS6 <- ggarrange(pS6A, pS6B, pS6C, pS6D, pS6E, pS6F, labels = LETTERS[1:6], ncol = 3, nrow = 2, font.label = list(size = 17, face = "bold"))
figS6
ggsave("FigureS6.pdf", plot = figS6, width = 12, height = 10)

