# =============================================================================
# WGCNA: co-expression modules of each population, their response to treatment and shore level, and GO enrichment of the
# heat-responsive modules.
# The rlog expression data and the network of each population (signed, power 8, minimum module size 50, merge cut height 0.25,
# randomSeed 12345) were saved in this folder beforehand.
# Output: FigureS8.pdf (module correspondence among populations), FigureS9.pdf (module eigengenes),
#         FigureS10.pdf (ANOVA and Tukey tests), FigureS11.pdf (GO enrichment of the heat-responsive modules)
# =============================================================================

setwd("C:/path/to/RNAgithub/07WGCNA/")

library(WGCNA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(ggtree)
library(patchwork)


# Colours of the three populations
population_colors <- c(North = "#87AADE", Central = "#FFE57F", South = "#FF9999")
eigengene_colors <- c(North = "#87AADE", Central = "#C99A20", South = "#FF9999")
eigengene_dark <- c(North = "#336DA9", Central = "#956F0A", South = "#C54F56")
eigengene_header <- c(North = "#C4D5EF", Central = "#F1D47A", South = "#F5C0C0")

############################################################
#Part 1: sample information, expression data and networks
count <- read.table("../01counts_ind_filter/counts_filter.txt", header = TRUE, row.names = 1, check.names = FALSE)
group <- read.delim("../01counts_ind_filter/group.txt", stringsAsFactors = FALSE, check.names = FALSE)
rownames(group) <- group$id
group <- group[colnames(count), ]
group$treatment_group <- factor(paste0(ifelse(group$temperature == "31degree", "31C", "17C"), "_", ifelse(group$time == "3h", "3h", "15min")), levels = c("17C_15min", "17C_3h", "31C_15min", "31C_3h"))
group$shore_level <- factor(group$shore_level, levels = c("high", "low"))
load("datExpr_North.RData")
load("datExpr_Central.RData")
load("datExpr_South.RData")
net_North <- readRDS("net_North_coexpression.rds")
net_Central <- readRDS("net_Central_coexpression.rds")
net_South <- readRDS("net_South_coexpression.rds")

############################################################
#Part 2: module eigengenes, two-way ANOVA and Tukey tests in each population (Figures S9 and S10)

#################### North

# Module of each gene, renamed Module_1 to Module_14 (grey = not assigned to any module).
group_North <- group[rownames(datExpr_North), ]
module_colors_North <- labels2colors(net_North$colors)
module_name_map_North <- setNames(paste0("Module_", 1:14), c("brown", "blue", "turquoise", "green", "red", "yellow", "black", "pink", "tan", "cyan", "magenta", "greenyellow", "purple", "salmon"))
ordered_colors_North <- names(module_name_map_North)

# Module eigengenes (first principal component of each module).
MEs_North <- moduleEigengenes(datExpr_North, colors = module_colors_North, excludeGrey = TRUE, align = "along average")$eigengenes
MEs_North <- MEs_North[, paste0("ME", ordered_colors_North)]
module_size_North <- as.integer(table(factor(module_colors_North, levels = ordered_colors_North)))

# Two-way ANOVA of each module eigengene: treatment, shore level and their interaction (type II sums of squares, as the
# design is unbalanced). P values are BH-adjusted across modules.
MEs_matrix_North <- as.matrix(setNames(MEs_North, unname(module_name_map_North)))
fit_full_North <- lm(MEs_matrix_North ~ treatment_group * shore_level, data = group_North)
rss_full_North <- colSums(residuals(fit_full_North)^2)
rss_additive_North <- colSums(residuals(lm(MEs_matrix_North ~ treatment_group + shore_level, data = group_North))^2)
rss_treatment_North <- colSums(residuals(lm(MEs_matrix_North ~ treatment_group, data = group_North))^2)
rss_shore_North <- colSums(residuals(lm(MEs_matrix_North ~ shore_level, data = group_North))^2)
anova_North <- data.frame(population = "North", module_name = unname(module_name_map_North), module_size = module_size_North, module_number = 1:14, df_residual = fit_full_North$df.residual,
                          F_treatment = unname((rss_shore_North - rss_additive_North) / 3 / (rss_full_North / fit_full_North$df.residual)),
                          F_shore = unname((rss_treatment_North - rss_additive_North) / (rss_full_North / fit_full_North$df.residual)),
                          F_interaction = unname((rss_additive_North - rss_full_North) / 3 / (rss_full_North / fit_full_North$df.residual)))
anova_North$FDR_treatment <- p.adjust(pf(anova_North$F_treatment, 3, anova_North$df_residual, lower.tail = FALSE), method = "BH")
anova_North$FDR_shore <- p.adjust(pf(anova_North$F_shore, 1, anova_North$df_residual, lower.tail = FALSE), method = "BH")
anova_North$FDR_interaction <- p.adjust(pf(anova_North$F_interaction, 3, anova_North$df_residual, lower.tail = FALSE), method = "BH")
anova_North$heat_responsive <- anova_North$FDR_treatment < 0.05

# Tukey tests between the four treatments in each module (two shore levels combined).
emm_North <- suppressMessages(emmeans::emmeans(fit_full_North, ~ treatment_group | module_name, mult.levs = list(module_name = unname(module_name_map_North))))
tukey_North <- as.data.frame(summary(emmeans::contrast(emm_North, method = "pairwise", adjust = "tukey")))
tukey_North$comparison <- paste0("P_", sub("-", "_vs_", gsub("[() ]", "", tukey_North$contrast)))
tukey_North$module_name <- as.character(tukey_North$module_name)
anova_North <- dplyr::left_join(anova_North, tidyr::pivot_wider(tukey_North[, c("module_name", "comparison", "p.value")], names_from = comparison, values_from = p.value), by = "module_name")

# Heat-responsive modules (significant treatment effect).
anova_North[anova_North$heat_responsive, c("module_name", "FDR_treatment", "FDR_interaction", "P_17C_15min_vs_31C_15min", "P_17C_15min_vs_31C_3h", "P_17C_3h_vs_31C_3h", "P_31C_15min_vs_31C_3h")]

# Genes of each module.
genes_North <- data.frame(gene = colnames(datExpr_North), module_name = unname(module_name_map_North[module_colors_North]), stringsAsFactors = FALSE)
table(genes_North$module_name, useNA = "ifany")

#################### Central

# Module of each gene, renamed Module_1 to Module_15 (grey = not assigned to any module).
group_Central <- group[rownames(datExpr_Central), ]
module_colors_Central <- labels2colors(net_Central$colors)
module_name_map_Central <- setNames(paste0("Module_", 1:15), c("black", "blue", "turquoise", "brown", "green", "magenta", "red", "midnightblue", "cyan", "salmon", "purple", "tan", "yellow", "greenyellow", "pink"))
ordered_colors_Central <- names(module_name_map_Central)

# Module eigengenes (first principal component of each module).
MEs_Central <- moduleEigengenes(datExpr_Central, colors = module_colors_Central, excludeGrey = TRUE, align = "along average")$eigengenes
MEs_Central <- MEs_Central[, paste0("ME", ordered_colors_Central)]
module_size_Central <- as.integer(table(factor(module_colors_Central, levels = ordered_colors_Central)))

# Two-way ANOVA of each module eigengene: treatment, shore level and their interaction (type II sums of squares, as the
# design is unbalanced). P values are BH-adjusted across modules.
MEs_matrix_Central <- as.matrix(setNames(MEs_Central, unname(module_name_map_Central)))
fit_full_Central <- lm(MEs_matrix_Central ~ treatment_group * shore_level, data = group_Central)
rss_full_Central <- colSums(residuals(fit_full_Central)^2)
rss_additive_Central <- colSums(residuals(lm(MEs_matrix_Central ~ treatment_group + shore_level, data = group_Central))^2)
rss_treatment_Central <- colSums(residuals(lm(MEs_matrix_Central ~ treatment_group, data = group_Central))^2)
rss_shore_Central <- colSums(residuals(lm(MEs_matrix_Central ~ shore_level, data = group_Central))^2)
anova_Central <- data.frame(population = "Central", module_name = unname(module_name_map_Central), module_size = module_size_Central, module_number = 1:15, df_residual = fit_full_Central$df.residual,
                            F_treatment = unname((rss_shore_Central - rss_additive_Central) / 3 / (rss_full_Central / fit_full_Central$df.residual)),
                            F_shore = unname((rss_treatment_Central - rss_additive_Central) / (rss_full_Central / fit_full_Central$df.residual)),
                            F_interaction = unname((rss_additive_Central - rss_full_Central) / 3 / (rss_full_Central / fit_full_Central$df.residual)))
anova_Central$FDR_treatment <- p.adjust(pf(anova_Central$F_treatment, 3, anova_Central$df_residual, lower.tail = FALSE), method = "BH")
anova_Central$FDR_shore <- p.adjust(pf(anova_Central$F_shore, 1, anova_Central$df_residual, lower.tail = FALSE), method = "BH")
anova_Central$FDR_interaction <- p.adjust(pf(anova_Central$F_interaction, 3, anova_Central$df_residual, lower.tail = FALSE), method = "BH")
anova_Central$heat_responsive <- anova_Central$FDR_treatment < 0.05

# Tukey tests between the four treatments in each module (two shore levels combined).
emm_Central <- suppressMessages(emmeans::emmeans(fit_full_Central, ~ treatment_group | module_name, mult.levs = list(module_name = unname(module_name_map_Central))))
tukey_Central <- as.data.frame(summary(emmeans::contrast(emm_Central, method = "pairwise", adjust = "tukey")))
tukey_Central$comparison <- paste0("P_", sub("-", "_vs_", gsub("[() ]", "", tukey_Central$contrast)))
tukey_Central$module_name <- as.character(tukey_Central$module_name)
anova_Central <- dplyr::left_join(anova_Central, tidyr::pivot_wider(tukey_Central[, c("module_name", "comparison", "p.value")], names_from = comparison, values_from = p.value), by = "module_name")

# Heat-responsive modules (significant treatment effect).
anova_Central[anova_Central$heat_responsive, c("module_name", "FDR_treatment", "FDR_interaction", "P_17C_15min_vs_31C_15min", "P_17C_15min_vs_31C_3h", "P_17C_3h_vs_31C_3h", "P_31C_15min_vs_31C_3h")]

# Genes of each module.
genes_Central <- data.frame(gene = colnames(datExpr_Central), module_name = unname(module_name_map_Central[module_colors_Central]), stringsAsFactors = FALSE)
table(genes_Central$module_name, useNA = "ifany")

#################### South

# Module of each gene, renamed Module_1 to Module_16 (grey = not assigned to any module).
group_South <- group[rownames(datExpr_South), ]
module_colors_South <- labels2colors(net_South$colors)
module_name_map_South <- setNames(paste0("Module_", 1:16), c("red", "turquoise", "blue", "brown", "green", "black", "purple", "cyan", "greenyellow", "lightcyan", "magenta", "midnightblue", "salmon", "pink", "tan", "yellow"))
ordered_colors_South <- names(module_name_map_South)

# Module eigengenes (first principal component of each module).
MEs_South <- moduleEigengenes(datExpr_South, colors = module_colors_South, excludeGrey = TRUE, align = "along average")$eigengenes
MEs_South <- MEs_South[, paste0("ME", ordered_colors_South)]
module_size_South <- as.integer(table(factor(module_colors_South, levels = ordered_colors_South)))

# Two-way ANOVA of each module eigengene: treatment, shore level and their interaction (type II sums of squares, as the
# design is unbalanced). P values are BH-adjusted across modules.
MEs_matrix_South <- as.matrix(setNames(MEs_South, unname(module_name_map_South)))
fit_full_South <- lm(MEs_matrix_South ~ treatment_group * shore_level, data = group_South)
rss_full_South <- colSums(residuals(fit_full_South)^2)
rss_additive_South <- colSums(residuals(lm(MEs_matrix_South ~ treatment_group + shore_level, data = group_South))^2)
rss_treatment_South <- colSums(residuals(lm(MEs_matrix_South ~ treatment_group, data = group_South))^2)
rss_shore_South <- colSums(residuals(lm(MEs_matrix_South ~ shore_level, data = group_South))^2)
anova_South <- data.frame(population = "South", module_name = unname(module_name_map_South), module_size = module_size_South, module_number = 1:16, df_residual = fit_full_South$df.residual,
                          F_treatment = unname((rss_shore_South - rss_additive_South) / 3 / (rss_full_South / fit_full_South$df.residual)),
                          F_shore = unname((rss_treatment_South - rss_additive_South) / (rss_full_South / fit_full_South$df.residual)),
                          F_interaction = unname((rss_additive_South - rss_full_South) / 3 / (rss_full_South / fit_full_South$df.residual)))
anova_South$FDR_treatment <- p.adjust(pf(anova_South$F_treatment, 3, anova_South$df_residual, lower.tail = FALSE), method = "BH")
anova_South$FDR_shore <- p.adjust(pf(anova_South$F_shore, 1, anova_South$df_residual, lower.tail = FALSE), method = "BH")
anova_South$FDR_interaction <- p.adjust(pf(anova_South$F_interaction, 3, anova_South$df_residual, lower.tail = FALSE), method = "BH")
anova_South$heat_responsive <- anova_South$FDR_treatment < 0.05

# Tukey tests between the four treatments in each module (two shore levels combined).
emm_South <- suppressMessages(emmeans::emmeans(fit_full_South, ~ treatment_group | module_name, mult.levs = list(module_name = unname(module_name_map_South))))
tukey_South <- as.data.frame(summary(emmeans::contrast(emm_South, method = "pairwise", adjust = "tukey")))
tukey_South$comparison <- paste0("P_", sub("-", "_vs_", gsub("[() ]", "", tukey_South$contrast)))
tukey_South$module_name <- as.character(tukey_South$module_name)
anova_South <- dplyr::left_join(anova_South, tidyr::pivot_wider(tukey_South[, c("module_name", "comparison", "p.value")], names_from = comparison, values_from = p.value), by = "module_name")

# Heat-responsive modules (significant treatment effect).
anova_South[anova_South$heat_responsive, c("module_name", "FDR_treatment", "FDR_interaction", "P_17C_15min_vs_31C_15min", "P_17C_15min_vs_31C_3h", "P_17C_3h_vs_31C_3h", "P_31C_15min_vs_31C_3h")]

# Genes of each module.
genes_South <- data.frame(gene = colnames(datExpr_South), module_name = unname(module_name_map_South[module_colors_South]), stringsAsFactors = FALSE)
table(genes_South$module_name, useNA = "ifany")

#################### Figure S10: ANOVA FDR of every module (A-C) and Tukey tests of the heat-responsive modules (D)

# Arrows in Figure S10D: mean eigengene at 31C 3h higher (up) or lower (down) than in all three other treatments.
emmeans_North <- as.data.frame(summary(emm_North))
direction_North <- dplyr::summarise(dplyr::group_by(emmeans_North, module_name), direction = ifelse(emmean[treatment_group == "31C_3h"] > max(emmean[treatment_group != "31C_3h"]), "\u2191", ifelse(emmean[treatment_group == "31C_3h"] < min(emmean[treatment_group != "31C_3h"]), "\u2193", "")))
anova_North$direction <- direction_North$direction[match(anova_North$module_name, direction_North$module_name)]
emmeans_Central <- as.data.frame(summary(emm_Central))
direction_Central <- dplyr::summarise(dplyr::group_by(emmeans_Central, module_name), direction = ifelse(emmean[treatment_group == "31C_3h"] > max(emmean[treatment_group != "31C_3h"]), "\u2191", ifelse(emmean[treatment_group == "31C_3h"] < min(emmean[treatment_group != "31C_3h"]), "\u2193", "")))
anova_Central$direction <- direction_Central$direction[match(anova_Central$module_name, direction_Central$module_name)]
emmeans_South <- as.data.frame(summary(emm_South))
direction_South <- dplyr::summarise(dplyr::group_by(emmeans_South, module_name), direction = ifelse(emmean[treatment_group == "31C_3h"] > max(emmean[treatment_group != "31C_3h"]), "\u2191", ifelse(emmean[treatment_group == "31C_3h"] < min(emmean[treatment_group != "31C_3h"]), "\u2193", "")))
anova_South$direction <- direction_South$direction[match(anova_South$module_name, direction_South$module_name)]

# Colour: -log10 adjusted P (FDR in A-C, Tukey-adjusted P in D); white = not significant.
scale_fill_S10 <- scale_fill_gradient(low = "#D8EEDC", high = "#1E803D", limits = c(-log10(0.05), 10), breaks = c(2, 5, 10), labels = c("2", "5", "\u226510"), na.value = "white", name = expression(plain("-log")[10]~italic(P)[adj]),
                                      guide = guide_colorbar(direction = "horizontal", title.position = "left", barwidth = grid::unit(3, "cm"), barheight = grid::unit(2.5, "mm")))
theme_S10 <- theme_minimal(base_size = 10, base_family = "sans") +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_text(colour = "black"), axis.text.y = element_text(size = 10, margin = margin(r = 2)),
        axis.text.x.bottom = element_text(size = 9.5, lineheight = 0.95, margin = margin(t = 5, b = 3)), plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 8)),
        plot.tag = element_text(size = 14, face = "bold"), plot.margin = margin(3, 4, 3, 4), legend.title = element_text(size = 10), legend.text = element_text(size = 9.5))

# A: North (heat-responsive modules in bold).
s10_North <- tidyr::pivot_longer(anova_North[, c("module_number", "FDR_treatment", "FDR_shore", "FDR_interaction")], cols = -module_number, names_to = "term", values_to = "FDR")
s10_North$x <- match(s10_North$term, c("FDR_treatment", "FDR_shore", "FDR_interaction"))
s10_North$fill <- ifelse(s10_North$FDR < 0.05, pmin(-log10(s10_North$FDR), 10), NA)
s10_North$label <- ifelse(s10_North$FDR < 0.001, sprintf("%.1e", s10_North$FDR), sprintf("%.2g", s10_North$FDR))
p_S10_North <- ggplot(s10_North, aes(x, 16 - module_number)) + geom_tile(aes(fill = fill), colour = "grey80", linewidth = 0.25) + geom_text(aes(label = label, colour = !is.na(fill) & fill > 5), size = 9.5, size.unit = "pt") +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") + scale_fill_S10 +
  scale_x_continuous(breaks = 1:3, labels = c("Treatment", "Shore level", "Interaction"), position = "bottom", limits = c(0.5, 3.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = 16 - anova_North$module_number, labels = parse(text = paste0(ifelse(anova_North$heat_responsive, "bold", "plain"), "(\"Module ", anova_North$module_number, " (n = ", anova_North$module_size, ")\")")), limits = c(1.5, 15.5), expand = c(0, 0)) +
  labs(title = "North", tag = "A", x = NULL, y = NULL) + theme_S10

# B: Central.
s10_Central <- tidyr::pivot_longer(anova_Central[, c("module_number", "FDR_treatment", "FDR_shore", "FDR_interaction")], cols = -module_number, names_to = "term", values_to = "FDR")
s10_Central$x <- match(s10_Central$term, c("FDR_treatment", "FDR_shore", "FDR_interaction"))
s10_Central$fill <- ifelse(s10_Central$FDR < 0.05, pmin(-log10(s10_Central$FDR), 10), NA)
s10_Central$label <- ifelse(s10_Central$FDR < 0.001, sprintf("%.1e", s10_Central$FDR), sprintf("%.2g", s10_Central$FDR))
p_S10_Central <- ggplot(s10_Central, aes(x, 16 - module_number)) + geom_tile(aes(fill = fill), colour = "grey80", linewidth = 0.25) + geom_text(aes(label = label, colour = !is.na(fill) & fill > 5), size = 9.5, size.unit = "pt") +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") + scale_fill_S10 +
  scale_x_continuous(breaks = 1:3, labels = c("Treatment", "Shore level", "Interaction"), position = "bottom", limits = c(0.5, 3.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = 16 - anova_Central$module_number, labels = parse(text = paste0(ifelse(anova_Central$heat_responsive, "bold", "plain"), "(\"Module ", anova_Central$module_number, " (n = ", anova_Central$module_size, ")\")")), limits = c(0.5, 15.5), expand = c(0, 0)) +
  labs(title = "Central", tag = "B", x = NULL, y = NULL) + theme_S10

# C: South.
s10_South <- tidyr::pivot_longer(anova_South[, c("module_number", "FDR_treatment", "FDR_shore", "FDR_interaction")], cols = -module_number, names_to = "term", values_to = "FDR")
s10_South$x <- match(s10_South$term, c("FDR_treatment", "FDR_shore", "FDR_interaction"))
s10_South$fill <- ifelse(s10_South$FDR < 0.05, pmin(-log10(s10_South$FDR), 10), NA)
s10_South$label <- ifelse(s10_South$FDR < 0.001, sprintf("%.1e", s10_South$FDR), sprintf("%.2g", s10_South$FDR))
p_S10_South <- ggplot(s10_South, aes(x, 17 - module_number)) + geom_tile(aes(fill = fill), colour = "grey80", linewidth = 0.25) + geom_text(aes(label = label, colour = !is.na(fill) & fill > 5), size = 9.5, size.unit = "pt") +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") + scale_fill_S10 +
  scale_x_continuous(breaks = 1:3, labels = c("Treatment", "Shore level", "Interaction"), position = "bottom", limits = c(0.5, 3.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = 17 - anova_South$module_number, labels = parse(text = paste0(ifelse(anova_South$heat_responsive, "bold", "plain"), "(\"Module ", anova_South$module_number, " (n = ", anova_South$module_size, ")\")")), limits = c(0.5, 16.5), expand = c(0, 0)) +
  labs(title = "South", tag = "C", x = NULL, y = NULL) + theme_S10

# D: Tukey tests of the heat-responsive modules, among the three no/short heat stress treatments (upper) and against long heat stress (lower).
heat_modules <- dplyr::bind_rows(anova_North, anova_Central, anova_South)
heat_modules <- heat_modules[heat_modules$heat_responsive, ]
heat_modules$row <- rev(seq_len(nrow(heat_modules)))
heat_modules$label <- paste0(heat_modules$population, ": Module ", heat_modules$module_number, " ", heat_modules$direction)

s10_tukey_short <- tidyr::pivot_longer(heat_modules[, c("row", "P_17C_15min_vs_31C_15min", "P_17C_15min_vs_17C_3h", "P_17C_3h_vs_31C_15min")], cols = -row, names_to = "comparison", values_to = "P")
s10_tukey_short$x <- match(s10_tukey_short$comparison, c("P_17C_15min_vs_31C_15min", "P_17C_15min_vs_17C_3h", "P_17C_3h_vs_31C_15min"))
s10_tukey_short$fill <- ifelse(s10_tukey_short$P < 0.05, pmin(-log10(s10_tukey_short$P), 10), NA)
s10_tukey_short$label <- ifelse(s10_tukey_short$P < 0.001, sprintf("%.1e", s10_tukey_short$P), sprintf("%.2g", s10_tukey_short$P))
p_S10_tukey_short <- ggplot(s10_tukey_short, aes(x, row)) + geom_tile(aes(fill = fill), colour = "grey80", linewidth = 0.25) + geom_text(aes(label = label, colour = !is.na(fill) & fill > 5), size = 9.5, size.unit = "pt") +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") + scale_fill_S10 +
  scale_x_continuous(breaks = 1:3, labels = c("17\u00b0C 15 min\nvs 31\u00b0C 15 min", "17\u00b0C 15 min\nvs 17\u00b0C 3 h", "17\u00b0C 3 h\nvs 31\u00b0C 15 min"), position = "bottom", limits = c(0.5, 3.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = heat_modules$row, labels = heat_modules$label, limits = c(0.5, nrow(heat_modules) + 0.5), expand = c(0, 0)) +
  labs(title = "Tukey comparisons", tag = "D", x = NULL, y = NULL) + theme_S10 + theme(axis.text.y = element_text(face = "bold"))

s10_tukey_long <- tidyr::pivot_longer(heat_modules[, c("row", "P_17C_15min_vs_31C_3h", "P_17C_3h_vs_31C_3h", "P_31C_15min_vs_31C_3h")], cols = -row, names_to = "comparison", values_to = "P")
s10_tukey_long$x <- match(s10_tukey_long$comparison, c("P_17C_15min_vs_31C_3h", "P_17C_3h_vs_31C_3h", "P_31C_15min_vs_31C_3h"))
s10_tukey_long$fill <- ifelse(s10_tukey_long$P < 0.05, pmin(-log10(s10_tukey_long$P), 10), NA)
s10_tukey_long$label <- ifelse(s10_tukey_long$P < 0.001, sprintf("%.1e", s10_tukey_long$P), sprintf("%.2g", s10_tukey_long$P))
p_S10_tukey_long <- ggplot(s10_tukey_long, aes(x, row)) + geom_tile(aes(fill = fill), colour = "grey80", linewidth = 0.25) + geom_text(aes(label = label, colour = !is.na(fill) & fill > 5), size = 9.5, size.unit = "pt") +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") + scale_fill_S10 +
  scale_x_continuous(breaks = 1:3, labels = c("17\u00b0C 15 min\nvs 31\u00b0C 3 h", "17\u00b0C 3 h\nvs 31\u00b0C 3 h", "31\u00b0C 15 min\nvs 31\u00b0C 3 h"), position = "bottom", limits = c(0.5, 3.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = heat_modules$row, labels = heat_modules$label, limits = c(0.5, nrow(heat_modules) + 0.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) + theme_S10 + theme(axis.text.y = element_text(face = "bold"), plot.margin = margin(14, 4, 3, 4))

# Combine A-D with one colour key.
figure_S10 <- p_S10_North + p_S10_Central + p_S10_South + p_S10_tukey_short + p_S10_tukey_long + plot_layout(design = "AB\nCD\nCE", heights = c(15, 7, 7), guides = "collect") & theme(legend.position = "bottom")
# Shade the labels of the heat-responsive modules in the colour of their population, as in Figure S9.
figure_S10 <- patchwork::patchworkGrob(figure_S10)
label_fills_S10 <- list(ifelse(anova_North$heat_responsive, eigengene_header[["North"]], NA),
                        ifelse(anova_Central$heat_responsive, eigengene_header[["Central"]], NA),
                        ifelse(anova_South$heat_responsive, eigengene_header[["South"]], NA),
                        unname(eigengene_header[heat_modules$population]), unname(eigengene_header[heat_modules$population]))
for (panel_S10 in seq_along(label_fills_S10)) {
  axis_id_S10 <- match(paste0("axis-l-", panel_S10), figure_S10$layout$name)
  axis_S10 <- figure_S10$grobs[[axis_id_S10]]$children$axis
  title_id_S10 <- which(vapply(axis_S10$grobs, inherits, logical(1), "titleGrob"))
  title_S10 <- axis_S10$grobs[[title_id_S10]]
  text_S10 <- title_S10$children[[1]]
  background_S10 <- grid::rectGrob(x = grid::unit(-3, "pt"), y = text_S10$y,
                                   width = grid::unit(1, "npc") + grid::unit(1, "pt"), height = grid::unit(0.88 / length(text_S10$label), "npc"),
                                   just = "left", gp = grid::gpar(fill = label_fills_S10[[panel_S10]], col = NA), name = "module_label_background")
  title_S10$children <- do.call(grid::gList, c(list(module_label_background = background_S10), as.list(title_S10$children)))
  title_S10$childrenOrder <- names(title_S10$children)
  axis_S10$grobs[[title_id_S10]] <- title_S10
  figure_S10$grobs[[axis_id_S10]]$children$axis <- axis_S10
}
grid::grid.newpage()
grid::grid.draw(figure_S10)
ggsave("FigureS10.pdf", figure_S10, width = 26, height = 30, units = "cm", device = cairo_pdf, bg = "white")

#################### Figure S9: module eigengenes under the four treatments
# Triangles = high shore, circles = low shore, diamonds = treatment means (two shore levels combined).
make_S9_population <- function(sample_info, eigengenes, anova_table, population, tag) {
  eigengene_data <- tidyr::pivot_longer(
    data.frame(sample_id = rownames(sample_info), treatment_group = sample_info$treatment_group,
               shore_level = sample_info$shore_level, eigengenes, check.names = FALSE),
    cols = starts_with("Module_"), names_to = "module_name", values_to = "ME")
  panel_labels <- paste0(ifelse(anova_table$heat_responsive, "bold", "plain"),
                         "(\"Module ", anova_table$module_number, " (n = ", anova_table$module_size, ")\")")
  eigengene_data$panel <- factor(panel_labels[match(eigengene_data$module_name, anova_table$module_name)], levels = panel_labels)
  eigengene_data$x <- as.integer(eigengene_data$treatment_group)
  means <- eigengene_data |>
    dplyr::group_by(module_name, panel, treatment_group, x) |>
    dplyr::summarise(ME = mean(ME), .groups = "drop")
  me_range <- range(eigengene_data$ME)
  me_span <- diff(me_range)
  # Brackets: Tukey test between 17C 15min and 31C 15min, shown for heat-responsive modules when significant
  p_column <- "P_17C_15min_vs_31C_15min"
  brackets <- anova_table[anova_table$heat_responsive & anova_table[[p_column]] < 0.05, c("module_name", p_column), drop = FALSE]
  brackets$panel <- factor(panel_labels[match(brackets$module_name, anova_table$module_name)], levels = panel_labels)
  brackets$y <- as.numeric(tapply(eigengene_data$ME, eigengene_data$module_name, max)[brackets$module_name]) + me_span * 0.10
  brackets$tip <- brackets$y - me_span * 0.035
  brackets$label_y <- brackets$y + me_span * 0.06
  brackets$label <- if (nrow(brackets)) paste0("italic(P)[adj] == '",
                                               ifelse(brackets[[p_column]] < 0.001, sprintf("%.1e", brackets[[p_column]]), sprintf("%.2g", brackets[[p_column]])), "'") else character(0)
  plot <- ggplot(eigengene_data, aes(x, ME)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
    geom_point(aes(shape = shore_level), position = position_jitter(width = 0.08, height = 0, seed = 1),
               colour = eigengene_colors[[population]], alpha = 1, size = 1.3, stroke = 0.25) +
    geom_line(data = means, aes(group = 1), colour = eigengene_dark[[population]], linewidth = 0.45) +
    geom_point(data = means, colour = eigengene_dark[[population]], size = 2.8, shape = 18) +
    geom_segment(data = brackets, aes(x = 1, xend = 3, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_segment(data = brackets, aes(x = 1, xend = 1, y = y, yend = tip), inherit.aes = FALSE, linewidth = 0.3) +
    geom_segment(data = brackets, aes(x = 3, xend = 3, y = y, yend = tip), inherit.aes = FALSE, linewidth = 0.3) +
    geom_text(data = brackets, aes(x = 2, y = label_y, label = label), inherit.aes = FALSE,
              size = 11.5, size.unit = "pt", vjust = 0, parse = TRUE) +
    ggh4x::facet_wrap2(~ panel, ncol = 4, labeller = label_parsed, axes = "x", remove_labels = "x",
                       strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(
                         fill = ifelse(anova_table$heat_responsive, eigengene_header[[population]], "grey94"),
                         colour = ifelse(anova_table$heat_responsive, eigengene_dark[[population]], NA),
                         linewidth = ifelse(anova_table$heat_responsive, 0.8, 0)))) +
    scale_x_continuous(breaks = 1:4, labels = c("17\u00b0C\n15 min", "17\u00b0C\n3 h", "31\u00b0C\n15 min", "31\u00b0C\n3 h"),
                       limits = c(0.65, 4.35), expand = c(0, 0)) +
    scale_y_continuous(limits = c(me_range[1] - me_span * 0.07, me_range[2] + me_span * 0.38),
                       breaks = scales::breaks_pretty(n = 4), expand = c(0, 0)) +
    scale_shape_manual(values = c(high = 17, low = 16), name = "Shore level",
                       labels = c(high = "High shore", low = "Low shore"),
                       guide = guide_legend(override.aes = list(size = 2.6, alpha = 1, colour = "black"))) +
    labs(title = population, tag = tag, x = NULL, y = "Module eigengene") +
    theme_bw(base_size = 11, base_family = "sans") +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
          panel.border = element_rect(colour = "grey35", linewidth = 0.35),
          axis.ticks = element_line(colour = "grey35", linewidth = 0.25),
          axis.text = element_text(colour = "black", size = 10), axis.text.x = element_text(lineheight = 0.95),
          axis.title.y = element_text(size = 11), strip.background = element_rect(fill = "grey94", colour = NA),
          strip.text = element_text(size = 11, margin = margin(3, 2, 3, 2)),
          panel.spacing = grid::unit(9, "pt"), plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.tag = element_text(size = 14, face = "bold"), plot.margin = margin(4, 6, 4, 4),
          legend.position = "bottom", legend.title = element_text(size = 11), legend.text = element_text(size = 11))
  list(plot = plot, eigengenes = eigengene_data, means = means, brackets = brackets)
}
s9_North <- make_S9_population(group_North, MEs_matrix_North, anova_North, "North", "A")
s9_Central <- make_S9_population(group_Central, MEs_matrix_Central, anova_Central, "Central", "B")
s9_South <- make_S9_population(group_South, MEs_matrix_South, anova_South, "South", "C")
p_S9_North <- s9_North$plot
p_S9_Central <- s9_Central$plot
p_S9_South <- s9_South$plot
figure_S9_page1 <- (p_S9_North / p_S9_Central) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
figure_S9_page2 <- p_S9_South
figure_S9_page1
figure_S9_page2
# Save Figure S9 as two pages (A-B and C)
pages_S9 <- c(tempfile(fileext = ".pdf"), tempfile(fileext = ".pdf"))
ggsave(pages_S9[1], figure_S9_page1, width = 26, height = 36, units = "cm", device = cairo_pdf, bg = "white")
ggsave(pages_S9[2], figure_S9_page2, width = 26, height = 18, units = "cm", device = cairo_pdf, bg = "white")
qpdf::pdf_combine(pages_S9, output = "FigureS9.pdf")
unlink(pages_S9)

############################################################
#Part 3: GO enrichment of the heat-responsive modules (Figure S11)

# Genes of the heat-responsive modules: North Modules 1, 2 and 13, Central Module 1, South Modules 1 and 16.

North_M1_genes <- genes_North$gene[genes_North$module_name %in% "Module_1"]
North_M2_genes <- genes_North$gene[genes_North$module_name %in% "Module_2"]
North_M13_genes <- genes_North$gene[genes_North$module_name %in% "Module_13"]
Central_M1_genes <- genes_Central$gene[genes_Central$module_name %in% "Module_1"]
South_M1_genes <- genes_South$gene[genes_South$module_name %in% "Module_1"]
South_M16_genes <- genes_South$gene[genes_South$module_name %in% "Module_16"]

# Gene-GO annotation (from 06GOKEGG); background = retained genes with at least one GO term.
gene_GO <- read.delim("../06GOKEGG/gene_go_description_all.txt", header = TRUE, stringsAsFactors = FALSE, quote = "")
gene_GO$gene <- trimws(gene_GO$gene)
gene_GO$go <- trimws(gene_GO$go)
gene_GO$description <- trimws(gene_GO$description)
gene_GO <- gene_GO[!is.na(gene_GO$gene) & gene_GO$gene != "" & grepl("^GO:[0-9]{7}$", gene_GO$go), ]
TERM2GENE <- dplyr::distinct(gene_GO, go, gene)
TERM2NAME <- dplyr::distinct(gene_GO, go, description)
GO_universe <- intersect(rownames(count), unique(gene_GO$gene))
length(GO_universe)

# GO enrichment of each module; a GO term is significant when p.adjust <= 0.01 and qvalue <= 0.01.


# North Module 1
length(North_M1_genes)
go_North_M1 <- enricher(gene = North_M1_genes, universe = GO_universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
North_M1_all <- dplyr::mutate(go_North_M1@result, Group = "North_M1")
sig_North_M1 <- North_M1_all[North_M1_all$p.adjust <= 0.01 & !is.na(North_M1_all$qvalue) & North_M1_all$qvalue <= 0.01, ]
nrow(sig_North_M1)

# North Module 2
length(North_M2_genes)
go_North_M2 <- enricher(gene = North_M2_genes, universe = GO_universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
North_M2_all <- dplyr::mutate(go_North_M2@result, Group = "North_M2")
sig_North_M2 <- North_M2_all[North_M2_all$p.adjust <= 0.01 & !is.na(North_M2_all$qvalue) & North_M2_all$qvalue <= 0.01, ]
nrow(sig_North_M2)

# North Module 13
length(North_M13_genes)
go_North_M13 <- enricher(gene = North_M13_genes, universe = GO_universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
North_M13_all <- dplyr::mutate(go_North_M13@result, Group = "North_M13")
sig_North_M13 <- North_M13_all[North_M13_all$p.adjust <= 0.01 & !is.na(North_M13_all$qvalue) & North_M13_all$qvalue <= 0.01, ]
nrow(sig_North_M13)

# Central Module 1
length(Central_M1_genes)
go_Central_M1 <- enricher(gene = Central_M1_genes, universe = GO_universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
Central_M1_all <- dplyr::mutate(go_Central_M1@result, Group = "Central_M1")
sig_Central_M1 <- Central_M1_all[Central_M1_all$p.adjust <= 0.01 & !is.na(Central_M1_all$qvalue) & Central_M1_all$qvalue <= 0.01, ]
nrow(sig_Central_M1)

# South Module 1
length(South_M1_genes)
go_South_M1 <- enricher(gene = South_M1_genes, universe = GO_universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
South_M1_all <- dplyr::mutate(go_South_M1@result, Group = "South_M1")
sig_South_M1 <- South_M1_all[South_M1_all$p.adjust <= 0.01 & !is.na(South_M1_all$qvalue) & South_M1_all$qvalue <= 0.01, ]
nrow(sig_South_M1)

# South Module 16
length(South_M16_genes)
go_South_M16 <- enricher(gene = South_M16_genes, universe = GO_universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
South_M16_all <- dplyr::mutate(go_South_M16@result, Group = "South_M16")
sig_South_M16 <- South_M16_all[South_M16_all$p.adjust <= 0.01 & !is.na(South_M16_all$qvalue) & South_M16_all$qvalue <= 0.01, ]
nrow(sig_South_M16)

# North Module 13 and South Module 16 have no significantly enriched GO terms, so the following analyses focus only on
# the modules whose genes are up-regulated under long heat stress (North Modules 1 and 2, Central Module 1, South Module 1).

#################### Figure S11: GO terms enriched in the up-regulated modules, ordered by a GO tree

# All GO terms significant in at least one of the four modules.
sig_GO_up <- dplyr::bind_rows(sig_North_M1, sig_North_M2, sig_Central_M1, sig_South_M1)
go_WGCNA <- dplyr::distinct(sig_GO_up, ID, Description)
GO_ids <- go_WGCNA$ID
nrow(go_WGCNA)

# Cluster these GO terms by their shared genes (Jaccard similarity, average linkage), as in Figure 3A.
gene_GO_WGCNA <- gene_GO[gene_GO$go %in% GO_ids, ]
go_WGCNA_all <- enricher(gene = unique(gene_GO_WGCNA$gene), TERM2GENE = gene_GO_WGCNA[c("go", "gene")], TERM2NAME = gene_GO_WGCNA[c("go", "description")],
                         minGSSize = 1, maxGSSize = Inf, pvalueCutoff = 1, qvalueCutoff = 1)
N_GO <- nrow(go_WGCNA_all@result)
go_WGCNA_sim <- pairwise_termsim(go_WGCNA_all, method = "JC", showCategory = N_GO)
GO_hclust <- hclust(as.dist(1 - go_WGCNA_sim@termsim), method = "average")
tree_WGCNA <- ggtree(ape::as.phylo(GO_hclust))

# The four functional categories (heat response genes, ribosome-related processes, Tf2-like retrotransposons and
# mitochondrial respiration), each defined by two GO terms at its ends.
clade_anchor_IDs <- c("GO:0042594", "GO:0051085", "GO:0031597", "GO:0022626", "GO:0002040", "GO:0039660", "GO:0006123", "GO:0045271")
clade_anchor_labels <- go_WGCNA_all@result$Description[match(clade_anchor_IDs, go_WGCNA_all@result$ID)]
node_heat <- ggtree::MRCA(tree_WGCNA, clade_anchor_labels[1], clade_anchor_labels[2])
node_ribosome <- ggtree::MRCA(tree_WGCNA, clade_anchor_labels[3], clade_anchor_labels[4])
node_Tf2 <- ggtree::MRCA(tree_WGCNA, clade_anchor_labels[5], clade_anchor_labels[6])
node_mito <- ggtree::MRCA(tree_WGCNA, clade_anchor_labels[7], clade_anchor_labels[8])
tree_WGCNA <- ggtree::flip(tree_WGCNA, node_heat, node_ribosome)
tree_WGCNA <- groupClade(tree_WGCNA, c(node_heat, node_ribosome, node_Tf2, node_mito), group_name = "group")

# GO terms in tree order (bottom to top) with their cluster.
tips_WGCNA <- dplyr::arrange(dplyr::filter(tree_WGCNA$data, isTip), y)
tips_WGCNA$ID <- go_WGCNA_all@result$ID[match(tips_WGCNA$label, go_WGCNA_all@result$Description)]

tips_WGCNA$Description <- go_WGCNA_all@result$Description[match(tips_WGCNA$ID, go_WGCNA_all@result$ID)]

# Wrap long GO names over two lines.
tips_WGCNA$Description_plot <- stringr::str_wrap(tips_WGCNA$Description, width = 58)
tips_WGCNA$row_height <- stringr::str_count(tips_WGCNA$Description_plot, "\n") + 1
tips_WGCNA$y_original <- tips_WGCNA$y
tips_WGCNA$y <- cumsum(tips_WGCNA$row_height) - tips_WGCNA$row_height / 2 + 0.5
tree_WGCNA$data$y <- approx(tips_WGCNA$y_original, tips_WGCNA$y, xout = tree_WGCNA$data$y, rule = 2)$y
Y_GO_S11 <- sum(tips_WGCNA$row_height)

# Width of the longest GO name, so that the bubbles start right after the names.
label_size_S11 <- 11
metric_file_S11 <- tempfile(fileext = ".pdf")
grDevices::cairo_pdf(metric_file_S11, family = "sans")
grid::pushViewport(grid::viewport(gp = grid::gpar(fontfamily = "sans", fontsize = label_size_S11)))
label_width_in_S11 <- max(grid::convertWidth(grid::stringWidth(tips_WGCNA$Description_plot), "inches", valueOnly = TRUE)) + 0.12
grid::popViewport()
grDevices::dev.off()
unlink(metric_file_S11)
x_unit_in_S11 <- (11 - grid::convertWidth(grid::unit(32, "pt"), "inches", valueOnly = TRUE)) / 14.3
label_xmin_S11 <- 2.7
label_xmax_S11 <- label_xmin_S11 + label_width_in_S11 / x_unit_in_S11
x_shift_S11 <- label_xmax_S11 - 3.90
figure_width_S11 <- 11 + (x_shift_S11 - 4.82) * x_unit_in_S11

# Bubbles: gene count and FDR of each GO term in each module (grey = not significant; no bubble = not tested).
GO_plot_all <- dplyr::filter(dplyr::bind_rows(North_M1_all, North_M2_all, Central_M1_all, South_M1_all), ID %in% GO_ids)
GO_plot_all <- dplyr::left_join(tidyr::expand_grid(ID = tips_WGCNA$ID, Group = c("North_M1", "North_M2", "Central_M1", "South_M1")), GO_plot_all, by = c("ID", "Group"))
GO_plot_all$Description <- go_WGCNA$Description[match(GO_plot_all$ID, go_WGCNA$ID)]
GO_plot_all$Count[is.na(GO_plot_all$Count)] <- 0
GO_plot_all$significant <- !is.na(GO_plot_all$p.adjust) & GO_plot_all$p.adjust <= 0.01 & !is.na(GO_plot_all$qvalue) & GO_plot_all$qvalue <= 0.01
GO_plot_all$ColorCol <- ifelse(GO_plot_all$significant, -log10(pmax(GO_plot_all$p.adjust, .Machine$double.xmin)), NA_real_)
GO_plot_all$y <- tips_WGCNA$y[match(GO_plot_all$ID, tips_WGCNA$ID)]
GO_plot_all$x <- unname(c(North_M1 = 4.425, North_M2 = 5.375, Central_M1 = 6.425, South_M1 = 7.475)[GO_plot_all$Group]) + x_shift_S11

# Tree branches coloured by category.
nodes_WGCNA <- as.data.frame(tree_WGCNA$data)
nodes_WGCNA$x <- nodes_WGCNA$x / max(nodes_WGCNA$x) * 2.6
branch_colors <- c("0" = "black", "1" = "#FDBF6F", "2" = "#B2DF8A", "3" = "#CAB2D6", "4" = "#a6cee3")
edges_WGCNA <- nodes_WGCNA[nodes_WGCNA$node != nodes_WGCNA$parent, ]
edges_WGCNA$x_parent <- nodes_WGCNA$x[match(edges_WGCNA$parent, nodes_WGCNA$node)]
edges_WGCNA$y_parent <- nodes_WGCNA$y[match(edges_WGCNA$parent, nodes_WGCNA$node)]
edges_WGCNA$color <- unname(branch_colors[as.character(edges_WGCNA$group)])

# Coloured background behind the GO names of each category, and category names on the right.
bands_WGCNA <- dplyr::summarise(dplyr::group_by(dplyr::filter(tips_WGCNA, as.character(group) != "0"), group), ymin = min(y - row_height / 2), ymax = max(y + row_height / 2), .groups = "drop")
bands_WGCNA$fill <- unname(branch_colors[as.character(bands_WGCNA$group)])
bands_WGCNA$label <- unname(c("1" = "Heat response genes", "2" = "Ribosome-related processes", "3" = "Tf2-like retrotransposons", "4" = "Mitochondrial\nrespiration")[as.character(bands_WGCNA$group)])
bands_WGCNA$text_color <- unname(c("1" = "#C46A00", "2" = "#3A6831", "3" = "#6A3D9A", "4" = "#1F78B4")[as.character(bands_WGCNA$group)])
bands_WGCNA$mid <- (bands_WGCNA$ymin + bands_WGCNA$ymax) / 2
bands_WGCNA$x_text <- ifelse(as.character(bands_WGCNA$group) == "3", 8.38, 8.66) + x_shift_S11

# Population headers (top), module names (bottom), and the viral and DNA metabolism GO terms outside the four categories.
headers_WGCNA <- data.frame(population = c("North", "Central", "South"), xmin = c(3.95, 5.95, 7.00), xmax = c(5.85, 6.90, 7.95), x = c(4.90, 6.425, 7.475))
headers_WGCNA[c("xmin", "xmax", "x")] <- headers_WGCNA[c("xmin", "xmax", "x")] + x_shift_S11
module_headers <- data.frame(x = c(4.425, 5.375, 6.425, 7.475), label = c("Module 1", "Module 2", "Module 1", "Module 1"))
module_headers$x <- module_headers$x + x_shift_S11
viral_y <- range(tips_WGCNA$y[tips_WGCNA$ID %in% c("GO:0019068", "GO:0042575")])

# Figure S11 is split into two pages after the ribosome-related processes.
cut_S11 <- min(tips_WGCNA$y[as.character(tips_WGCNA$group) == "2"] - tips_WGCNA$row_height[as.character(tips_WGCNA$group) == "2"] / 2)

# Page 1: heat response genes and ribosome-related processes.
ymin_page1_S11 <- cut_S11
ymax_page1_S11 <- Y_GO_S11 + 0.5
figure_height_page1_S11 <- max(9, (ymax_page1_S11 - ymin_page1_S11 + 10.5) * 0.18 + 0.45)
tips_page1_S11 <- tips_WGCNA[tips_WGCNA$y > ymin_page1_S11 & tips_WGCNA$y < ymax_page1_S11, ]
bubbles_page1_S11 <- GO_plot_all[GO_plot_all$ID %in% tips_page1_S11$ID & GO_plot_all$Count > 0, ]
edges_h_page1_S11 <- edges_WGCNA[edges_WGCNA$y > ymin_page1_S11 & edges_WGCNA$y < ymax_page1_S11, ]
edges_v_page1_S11 <- edges_WGCNA[pmax(edges_WGCNA$y, edges_WGCNA$y_parent) >= ymin_page1_S11 & pmin(edges_WGCNA$y, edges_WGCNA$y_parent) <= ymax_page1_S11, ]
edges_v_page1_S11$y <- pmax(ymin_page1_S11, pmin(edges_v_page1_S11$y, ymax_page1_S11))
edges_v_page1_S11$y_parent <- pmax(ymin_page1_S11, pmin(edges_v_page1_S11$y_parent, ymax_page1_S11))
bands_page1_S11 <- bands_WGCNA[as.character(bands_WGCNA$group) %in% c("1", "2"), ]
# Tree, GO descriptions, bubbles and labels share one y axis, so every row lines up.
figS11_page1 <- ggplot() +
  geom_segment(data = edges_h_page1_S11, aes(x = x_parent, xend = x, y = y, yend = y), color = edges_h_page1_S11$color, linewidth = 0.85, lineend = "butt") +
  geom_segment(data = edges_v_page1_S11, aes(x = x_parent, xend = x_parent, y = y_parent, yend = y), color = edges_v_page1_S11$color, linewidth = 0.85, lineend = "butt") +
  geom_rect(data = bands_page1_S11, aes(ymin = ymin, ymax = ymax), xmin = label_xmin_S11, xmax = label_xmax_S11, fill = bands_page1_S11$fill, alpha = 0.45, color = NA) +
  geom_text(data = tips_page1_S11, aes(y = y, label = Description_plot), x = label_xmin_S11, hjust = 0, size = label_size_S11, size.unit = "pt", family = "sans", lineheight = 0.9) +
  geom_segment(data = tips_page1_S11, aes(y = y, yend = y), x = 3.95 + x_shift_S11, xend = 7.95 + x_shift_S11, color = "grey90", linewidth = 0.2) +
  geom_rect(data = headers_WGCNA, aes(xmin = xmin, xmax = xmax), ymin = ymax_page1_S11 + 0.7, ymax = ymax_page1_S11 + 2.7, fill = "grey90", color = "black", linewidth = 0.25) +
  geom_text(data = headers_WGCNA, aes(x = x, label = population), y = ymax_page1_S11 + 1.7, size = 13, size.unit = "pt") +
  geom_point(data = bubbles_page1_S11, aes(x = x, y = y, color = ColorCol, size = Count)) +
  scale_size_continuous(range = c(1, 4), limits = range(GO_plot_all$Count[GO_plot_all$Count > 0]), breaks = c(50, 100), name = "Gene count",
                        guide = guide_legend(order = 2, nrow = 1, title.position = "left", override.aes = list(shape = 21, fill = NA, color = "black", stroke = 0.4))) +
  scale_color_gradient(low = "#97D1A0", high = "#1E803D", limits = range(GO_plot_all$ColorCol, na.rm = TRUE), na.value = "grey80", breaks = c(50, 100), name = expression(plain("-log")[10]~FDR),
                       guide = guide_colorbar(order = 1, direction = "horizontal", title.position = "left", barwidth = grid::unit(1.8, "cm"), barheight = grid::unit(0.35, "cm"))) +
  geom_text(data = module_headers, aes(x = x, label = label), y = ymin_page1_S11 - 1.1, size = 11, size.unit = "pt") +
  geom_segment(data = bands_page1_S11, aes(y = ymin + 0.2, yend = ymax - 0.2), x = 8.13 + x_shift_S11, xend = 8.13 + x_shift_S11, color = bands_page1_S11$fill, linewidth = 0.6) +
  geom_text(data = bands_page1_S11, aes(x = x_text, y = mid, label = label), angle = 90, fontface = "bold", size = 13, size.unit = "pt", color = bands_page1_S11$text_color) +
  scale_x_continuous(limits = c(-0.1, 9.38 + x_shift_S11), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + coord_cartesian(ylim = c(ymin_page1_S11 - 7, ymax_page1_S11 + 3.5), expand = FALSE) +
  annotate("text", x = 0, y = ymax_page1_S11 + 1.7, label = "Figure S11 (1/2)", hjust = 0, size = 11, size.unit = "pt") +
  theme_void(base_family = "sans") + theme(legend.position = "inside", legend.position.inside = c(0.64, 0.015), legend.justification = c(0.5, 0), legend.box = "horizontal",
                                           legend.title = element_text(size = 12), legend.text = element_text(size = 11), plot.margin = margin(21, 12, 11, 20),
                                           plot.background = element_rect(fill = "white", color = NA))
figS11_page1
pages_S11 <- c(tempfile(fileext = ".pdf"), tempfile(fileext = ".pdf"))
ggsave(pages_S11[1], figS11_page1, width = figure_width_S11, height = figure_height_page1_S11, device = cairo_pdf, limitsize = FALSE, bg = "white")

# Page 2: Tf2-like retrotransposons, mitochondrial respiration and the remaining GO terms.
ymin_page2_S11 <- 0.5
ymax_page2_S11 <- cut_S11
figure_height_page2_S11 <- max(9, (ymax_page2_S11 - ymin_page2_S11 + 10.5) * 0.18 + 0.45)
tips_page2_S11 <- tips_WGCNA[tips_WGCNA$y > ymin_page2_S11 & tips_WGCNA$y < ymax_page2_S11, ]
bubbles_page2_S11 <- GO_plot_all[GO_plot_all$ID %in% tips_page2_S11$ID & GO_plot_all$Count > 0, ]
edges_h_page2_S11 <- edges_WGCNA[edges_WGCNA$y > ymin_page2_S11 & edges_WGCNA$y < ymax_page2_S11, ]
edges_v_page2_S11 <- edges_WGCNA[pmax(edges_WGCNA$y, edges_WGCNA$y_parent) >= ymin_page2_S11 & pmin(edges_WGCNA$y, edges_WGCNA$y_parent) <= ymax_page2_S11, ]
edges_v_page2_S11$y <- pmax(ymin_page2_S11, pmin(edges_v_page2_S11$y, ymax_page2_S11))
edges_v_page2_S11$y_parent <- pmax(ymin_page2_S11, pmin(edges_v_page2_S11$y_parent, ymax_page2_S11))
bands_page2_S11 <- bands_WGCNA[as.character(bands_WGCNA$group) %in% c("3", "4"), ]
figS11_page2 <- ggplot() +
  geom_segment(data = edges_h_page2_S11, aes(x = x_parent, xend = x, y = y, yend = y), color = edges_h_page2_S11$color, linewidth = 0.85, lineend = "butt") +
  geom_segment(data = edges_v_page2_S11, aes(x = x_parent, xend = x_parent, y = y_parent, yend = y), color = edges_v_page2_S11$color, linewidth = 0.85, lineend = "butt") +
  geom_rect(data = bands_page2_S11, aes(ymin = ymin, ymax = ymax), xmin = label_xmin_S11, xmax = label_xmax_S11, fill = bands_page2_S11$fill, alpha = 0.45, color = NA) +
  geom_text(data = tips_page2_S11, aes(y = y, label = Description_plot), x = label_xmin_S11, hjust = 0, size = label_size_S11, size.unit = "pt", family = "sans", lineheight = 0.9) +
  geom_segment(data = tips_page2_S11, aes(y = y, yend = y), x = 3.95 + x_shift_S11, xend = 7.95 + x_shift_S11, color = "grey90", linewidth = 0.2) +
  geom_rect(data = headers_WGCNA, aes(xmin = xmin, xmax = xmax), ymin = ymax_page2_S11 + 0.7, ymax = ymax_page2_S11 + 2.7, fill = "grey90", color = "black", linewidth = 0.25) +
  geom_text(data = headers_WGCNA, aes(x = x, label = population), y = ymax_page2_S11 + 1.7, size = 13, size.unit = "pt") +
  geom_point(data = bubbles_page2_S11, aes(x = x, y = y, color = ColorCol, size = Count)) +
  scale_size_continuous(range = c(1, 4), limits = range(GO_plot_all$Count[GO_plot_all$Count > 0]), breaks = c(50, 100), name = "Gene count",
                        guide = guide_legend(order = 2, nrow = 1, title.position = "left", override.aes = list(shape = 21, fill = NA, color = "black", stroke = 0.4))) +
  scale_color_gradient(low = "#97D1A0", high = "#1E803D", limits = range(GO_plot_all$ColorCol, na.rm = TRUE), na.value = "grey80", breaks = c(50, 100), name = expression(plain("-log")[10]~FDR),
                       guide = guide_colorbar(order = 1, direction = "horizontal", title.position = "left", barwidth = grid::unit(1.8, "cm"), barheight = grid::unit(0.35, "cm"))) +
  geom_text(data = module_headers, aes(x = x, label = label), y = ymin_page2_S11 - 1.1, size = 11, size.unit = "pt") +
  geom_segment(data = bands_page2_S11, aes(y = ymin + 0.2, yend = ymax - 0.2), x = 8.13 + x_shift_S11, xend = 8.13 + x_shift_S11, color = bands_page2_S11$fill, linewidth = 0.6) +
  geom_text(data = bands_page2_S11, aes(x = x_text, y = mid, label = label), angle = 90, fontface = "bold", size = 13, size.unit = "pt", color = bands_page2_S11$text_color) +
  annotate("segment", x = 8.68 + x_shift_S11, xend = 8.68 + x_shift_S11, y = viral_y[1], yend = viral_y[2], linewidth = 0.5) +
  annotate("text", x = 9.06 + x_shift_S11, y = mean(viral_y), label = "Viral processes and DNA metabolism", angle = 90, size = 13, size.unit = "pt", fontface = "bold") +
  scale_x_continuous(limits = c(-0.1, 9.38 + x_shift_S11), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + coord_cartesian(ylim = c(ymin_page2_S11 - 7, ymax_page2_S11 + 3.5), expand = FALSE) +
  annotate("text", x = 0, y = ymax_page2_S11 + 1.7, label = "Figure S11 (2/2)", hjust = 0, size = 11, size.unit = "pt") +
  theme_void(base_family = "sans") + theme(legend.position = "inside", legend.position.inside = c(0.64, 0.015), legend.justification = c(0.5, 0), legend.box = "horizontal",
                                           legend.title = element_text(size = 12), legend.text = element_text(size = 11), plot.margin = margin(21, 12, 11, 20),
                                           plot.background = element_rect(fill = "white", color = NA))
figS11_page2
ggsave(pages_S11[2], figS11_page2, width = figure_width_S11, height = figure_height_page2_S11, device = cairo_pdf, limitsize = FALSE, bg = "white")
qpdf::pdf_combine(pages_S11, output = "FigureS11.pdf")
unlink(pages_S11)


############################################################
#Part 4: correspondence of modules among the three populations (Figure S8)

# Module of each gene in each population (Unassigned = not in any module).
networks <- list(North = net_North, Central = net_Central, South = net_South)
module_maps <- list(North = module_name_map_North, Central = module_name_map_Central, South = module_name_map_South)
common_genes <- Reduce(intersect, lapply(networks, function(net) names(net$colors)))
positions <- list()
module_labels <- list()
for (pop in names(networks)) {
  gene_colors <- WGCNA::labels2colors(networks[[pop]]$colors)
  names(gene_colors) <- names(networks[[pop]]$colors)
  modules <- unname(module_maps[[pop]][gene_colors[common_genes]])
  modules[is.na(modules)] <- "Unassigned"
  positions[[pop]] <- data.frame(gene = common_genes, module = modules)
  module_order <- c(unname(module_maps[[pop]]), "Unassigned")
  module_labels[[pop]] <- data.frame(population = pop, module = module_order, label = sub("_", " ", module_order), genes = as.integer(table(factor(modules, levels = module_order))))
}
population_gene_flow <- list(positions = positions, labels = bind_rows(module_labels), common_genes = common_genes)

# Figure S8B: pairs of modules from two populations that share more genes than expected (only adjacent populations are shown).
make_module_correspondence <- function(gene_flow, module_maps) {
  populations <- names(gene_flow$positions)
  universe <- length(gene_flow$common_genes)
  labels <- gene_flow$labels[gene_flow$labels$module != "Unassigned", ]
  labels$number <- as.integer(sub("Module_", "", labels$module))
  labels$key <- paste(labels$population, labels$module, sep = ":")
  nodes <- labels[order(match(labels$population, populations), labels$number), ]
  pairs <- combn(populations, 2, simplify = FALSE)
  comparisons <- do.call(rbind, lapply(pairs, function(pair) {
    a <- gene_flow$positions[[pair[1]]]
    b <- gene_flow$positions[[pair[2]]]
    b <- b[match(a$gene, b$gene), ]
    modules_a <- unname(module_maps[[pair[1]]])
    modules_b <- unname(module_maps[[pair[2]]])
    tab <- table(factor(a$module, levels = modules_a), factor(b$module, levels = modules_b))
    d <- as.data.frame(tab, stringsAsFactors = FALSE)
    names(d) <- c("module_a", "module_b", "overlap")
    d$population_a <- pair[1]
    d$population_b <- pair[2]
    d$size_a <- nodes$genes[match(paste(pair[1], d$module_a, sep = ":"), nodes$key)]
    d$size_b <- nodes$genes[match(paste(pair[2], d$module_b, sep = ":"), nodes$key)]
    d$log_p <- phyper(d$overlap - 1, d$size_a, universe - d$size_a, d$size_b,
                      lower.tail = FALSE, log.p = TRUE)
    d$jaccard <- d$overlap / (d$size_a + d$size_b - d$overlap)
    d
  }))
  # Hypergeometric test for every pair of modules, BH-adjusted over all pairs (in log scale, so that very small P values are kept).
  o <- order(comparisons$log_p)
  log_q_sorted <- comparisons$log_p[o] + log(nrow(comparisons)) - log(seq_len(nrow(comparisons)))
  comparisons$log_fdr <- NA_real_
  comparisons$log_fdr[o] <- pmin(0, rev(cummin(rev(log_q_sorted))))
  comparisons$minus_log10_fdr <- -comparisons$log_fdr / log(10)
  comparisons$significant <- comparisons$log_fdr < log(0.05)
  
  # Stack the modules of each population, with gaps between them.
  stack_height <- max(tapply(nodes$genes, nodes$population, sum)) + 7000
  nodes$ymin <- nodes$ymax <- nodes$y <- rep(NA_real_, nrow(nodes))
  for (pop in populations) {
    idx <- which(nodes$population == pop)
    padding <- (stack_height - sum(nodes$genes[idx])) / (length(idx) - 1)
    used <- c(0, head(cumsum(nodes$genes[idx] + padding), -1))
    nodes$ymax[idx] <- stack_height - used
    nodes$ymin[idx] <- nodes$ymax[idx] - nodes$genes[idx]
    nodes$y[idx] <- (nodes$ymin[idx] + nodes$ymax[idx]) / 2
  }
  nodes$x <- match(nodes$population, populations)
  shown <- comparisons[comparisons$significant &
                         match(comparisons$population_b, populations) - match(comparisons$population_a, populations) == 1, ]
  shown$from <- paste(shown$population_a, shown$module_a, sep = ":")
  shown$to <- paste(shown$population_b, shown$module_b, sep = ":")
  shown$x <- nodes$x[match(shown$from, nodes$key)] + 0.045
  shown$xend <- nodes$x[match(shown$to, nodes$key)] - 0.045
  shown$y <- nodes$y[match(shown$from, nodes$key)]
  shown$yend <- nodes$y[match(shown$to, nodes$key)]
  shown$score <- pmin(shown$minus_log10_fdr, 200)
  shown <- shown[order(shown$overlap), ]
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = shown,
                          ggplot2::aes(x = x, xend = xend, y = y, yend = yend,
                                       linewidth = overlap, colour = score), alpha = 0.8, lineend = "round") +
    ggplot2::scale_linewidth_continuous(name = "Shared genes", range = c(0.15, 2.5),
                                        breaks = c(100, 500, 1000, 2000)) +
    ggplot2::scale_colour_gradient(name = expression(plain("-log")[10]~FDR),
                                   low = "#D8EEDC", high = "#1E803D", limits = c(0, 200),
                                   breaks = c(0, 50, 100, 200), labels = c("0", "50", "100", "200+"),
                                   guide = ggplot2::guide_colourbar(direction = "horizontal", title.position = "top",
                                                                    barwidth = grid::unit(3, "cm"), barheight = grid::unit(3, "mm"))) +
    ggplot2::geom_rect(data = nodes,
                       ggplot2::aes(xmin = x - 0.045, xmax = x + 0.045, ymin = ymin, ymax = ymax, fill = population),
                       colour = "grey40", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = population_colors, guide = "none") +
    ggplot2::scale_x_continuous(breaks = 1:3, labels = populations, position = "top",
                                limits = c(0.45, 3.7), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, stack_height), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 10, base_family = "sans") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 13, face = "bold", colour = "black",
                                                       margin = ggplot2::margin(b = 8)),
                   legend.position = "bottom", legend.box = "horizontal",
                   legend.title = ggplot2::element_text(size = 10),
                   legend.text = ggplot2::element_text(size = 8),
                   legend.margin = ggplot2::margin(4, 4, 4, 4),
                   plot.margin = ggplot2::margin(12, 8, 6, 8))
  for (pop in populations) {
    d <- nodes[nodes$population == pop, ]
    left_side <- pop == "North"
    p <- p + ggrepel::geom_label_repel(data = d, ggplot2::aes(x = x, y = y, label = label),
                                       nudge_x = if (left_side) -0.1 else 0.12, hjust = if (left_side) 1 else 0,
                                       direction = "y", size = 10 / ggplot2::.pt, family = "sans", fill = "white", colour = "black",
                                       label.size = 0, label.padding = grid::unit(1, "pt"), box.padding = grid::unit(2, "pt"),
                                       point.padding = grid::unit(2, "pt"), segment.colour = "grey45", segment.size = 0.25,
                                       min.segment.length = 0, max.overlaps = Inf, max.iter = 10000, max.time = Inf, seed = 1,
                                       ylim = c(0, stack_height))
  }
  list(plot = p, nodes = nodes, comparisons = comparisons, displayed_links = shown,
       universe = universe)
}

population_module_flow <- make_module_correspondence(population_gene_flow,
                                                     list(North = module_name_map_North, Central = module_name_map_Central, South = module_name_map_South))

# Figure S8A: distribution of the genes of each module among the modules of the next population.
make_complete_module_flow <- function(gene_flow, module_maps) {
  populations <- names(gene_flow$positions)
  nodes <- gene_flow$labels
  nodes$key <- paste(nodes$population, nodes$module, sep = ":")
  total <- length(gene_flow$common_genes)
  stack_height <- total + 7000
  nodes$ymin <- nodes$ymax <- nodes$y <- rep(NA_real_, nrow(nodes))
  for (pop in populations) {
    levels <- c(unname(module_maps[[pop]]), "Unassigned")
    idx <- which(nodes$population == pop)
    idx <- idx[order(match(nodes$module[idx], levels))]
    padding <- (stack_height - total) / (length(idx) - 1)
    used <- c(0, head(cumsum(nodes$genes[idx] + padding), -1))
    nodes$ymax[idx] <- stack_height - used
    nodes$ymin[idx] <- nodes$ymax[idx] - nodes$genes[idx]
    nodes$y[idx] <- (nodes$ymin[idx] + nodes$ymax[idx]) / 2
  }
  nodes$x <- match(nodes$population, populations)
  links <- do.call(rbind, lapply(1:2, function(gap) {
    pop_a <- populations[gap]
    pop_b <- populations[gap + 1]
    a <- gene_flow$positions[[pop_a]]
    b <- gene_flow$positions[[pop_b]]
    b <- b[match(a$gene, b$gene), ]
    levels_a <- c(unname(module_maps[[pop_a]]), "Unassigned")
    levels_b <- c(unname(module_maps[[pop_b]]), "Unassigned")
    d <- as.data.frame(table(factor(a$module, levels = levels_a), factor(b$module, levels = levels_b)), stringsAsFactors = FALSE)
    names(d) <- c("module_a", "module_b", "genes")
    d <- d[d$genes > 0, ]
    d$gap <- gap
    d$population_a <- pop_a
    d$population_b <- pop_b
    d$from <- paste(pop_a, d$module_a, sep = ":")
    d$to <- paste(pop_b, d$module_b, sep = ":")
    d$top_a <- d$top_b <- rep(NA_real_, nrow(d))
    for (module in levels_a) {
      idx <- which(d$module_a == module)
      idx <- idx[order(match(d$module_b[idx], levels_b))]
      node <- nodes[nodes$key == paste(pop_a, module, sep = ":"), ]
      d$top_a[idx] <- node$ymax - c(0, head(cumsum(d$genes[idx]), -1))
    }
    for (module in levels_b) {
      idx <- which(d$module_b == module)
      idx <- idx[order(match(d$module_a[idx], levels_a))]
      node <- nodes[nodes$key == paste(pop_b, module, sep = ":"), ]
      d$top_b[idx] <- node$ymax - c(0, head(cumsum(d$genes[idx]), -1))
    }
    d$bottom_a <- d$top_a - d$genes
    d$bottom_b <- d$top_b - d$genes
    d$assigned <- d$module_a != "Unassigned" & d$module_b != "Unassigned"
    d
  }))
  links$id <- seq_len(nrow(links))
  t <- seq(0, 1, length.out = 32)
  ease <- 3 * t^2 - 2 * t^3
  polygons <- do.call(rbind, lapply(seq_len(nrow(links)), function(i) {
    d <- links[i, ]
    x <- d$gap + 0.045 + t * 0.91
    data.frame(x = c(x, rev(x)),
               y = c(d$top_a + ease * (d$top_b - d$top_a), rev(d$bottom_a + ease * (d$bottom_b - d$bottom_a))),
               id = i, assigned = d$assigned)
  }))
  nodes$colour <- ifelse(nodes$module == "Unassigned", "#DEDEDE", population_colors[nodes$population])
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = polygons[!polygons$assigned, ],
                          ggplot2::aes(x, y, group = id), fill = "grey78", alpha = 0.35, colour = NA) +
    ggplot2::geom_polygon(data = polygons[polygons$assigned, ],
                          ggplot2::aes(x, y, group = id), fill = "#1E803D", alpha = 0.35, colour = NA) +
    ggplot2::geom_rect(data = nodes,
                       ggplot2::aes(xmin = x - 0.045, xmax = x + 0.045, ymin = ymin, ymax = ymax),
                       fill = nodes$colour, colour = "grey40", linewidth = 0.25) +
    ggplot2::scale_x_continuous(breaks = 1:3, labels = populations, position = "top",
                                limits = c(0.35, 3.85), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, stack_height), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 10, base_family = "sans") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 13, face = "bold", colour = "black", margin = ggplot2::margin(b = 8)),
                   legend.position = "bottom", legend.title = ggplot2::element_text(size = 10), legend.text = ggplot2::element_text(size = 8),
                   plot.caption = ggplot2::element_text(size = 10, hjust = 0.5, margin = ggplot2::margin(t = 4)),
                   plot.margin = ggplot2::margin(12, 10, 6, 8))
  for (pop in populations) {
    d <- nodes[nodes$population == pop, ]
    left_side <- pop == "North"
    p <- p + ggrepel::geom_label_repel(data = d, ggplot2::aes(x = x, y = y, label = label),
                                       nudge_x = if (left_side) -0.1 else 0.12, hjust = if (left_side) 1 else 0,
                                       direction = "y", size = 10 / ggplot2::.pt, family = "sans", colour = "black", fill = "white",
                                       label.size = 0, label.padding = grid::unit(1, "pt"), box.padding = grid::unit(1, "pt"),
                                       point.padding = grid::unit(1, "pt"), segment.colour = "grey45", segment.size = 0.2,
                                       min.segment.length = 0, max.overlaps = Inf, max.iter = 10000, max.time = Inf, seed = 1,
                                       ylim = c(0, stack_height))
  }
  list(plot = p, nodes = nodes, links = links, total_genes = total)
}

population_complete_flow <- make_complete_module_flow(population_gene_flow,
                                                      list(North = module_name_map_North, Central = module_name_map_Central, South = module_name_map_South))

# Combine A and B.
module_correspondence_A <- population_complete_flow$plot + ggplot2::labs(tag = "A") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 14, face = "bold"))
module_correspondence_B <- population_module_flow$plot + ggplot2::labs(tag = "B") +
  ggplot2::guides(linewidth = ggplot2::guide_legend(nrow = 1, title.position = "top")) +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 14, face = "bold"),
                 legend.box.spacing = grid::unit(20, "pt"), legend.spacing.x = grid::unit(12, "pt"),
                 legend.title = ggplot2::element_text(size = 10, margin = ggplot2::margin(b = 5)))
figure_module_correspondence_AB <- (patchwork::free(module_correspondence_A, side = "b") | module_correspondence_B) +
  patchwork::plot_layout(widths = c(1, 1), guides = "keep") &
  ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")
figure_module_correspondence_AB
ggsave("FigureS8.pdf", figure_module_correspondence_AB, width = 26, height = 26, units = "cm", device = cairo_pdf, bg = "white")
