setwd("C:/path/to/RNAgithub/08Frontloading/")

# =============================================================================
# Figure 4: baseline expression and plasticity of the 43 heat response genes shared by all three populations
# Control: 17°C for 15 min; long heat stress: 31°C for 3 h
# Input: counts_filter.txt, group.txt, heat_intersect_43genes_withannotation_Manual_curation.txt
# Output: Figure4.pdf
# =============================================================================
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(cowplot)
library(edgeR)
library(grid)
counts <- read.table("../01counts_ind_filter/counts_filter.txt", header = TRUE)
heat_curated <- read.delim("../06GOKEGG/heat_intersect_43genes_withannotation_Manual_curation.txt", stringsAsFactors = FALSE)
hsp70_genes <- heat_curated$gene[heat_curated$family == "Hsp70"]

# Expression level of each gene: log2 CPM (edgeR) from the final filtered count matrix.
logcpm <- cpm(as.matrix(counts), log = TRUE)

# Compare control (17°C, 15 min) with long heat stress (31°C, 3 h).
group <- read.delim("../01counts_ind_filter/group.txt", stringsAsFactors = FALSE)
group$treatment <- ifelse(group$temperature == "17degree" & group$time == "15min", "control", ifelse(group$temperature == "31degree" & group$time == "3h", "heat", NA))
group_ctrl_heat <- group[!is.na(group$treatment), ]
logcpm_heat_shared <- logcpm[rownames(logcpm) %in% heat_curated$gene, colnames(logcpm) %in% group_ctrl_heat$id, drop = FALSE]

# Mean expression of each gene in each population under control and heat stress.
expr_long <- tidyr::pivot_longer(tibble::rownames_to_column(as.data.frame(logcpm_heat_shared), "gene"), -gene, names_to = "id", values_to = "expression")
expr_long <- dplyr::left_join(expr_long, group_ctrl_heat[, c("id", "population", "treatment")], by = "id")
expr_mean <- tidyr::pivot_wider(dplyr::summarise(dplyr::group_by(expr_long, gene, population, treatment), mean_expr = mean(expression), .groups = "drop"), names_from = treatment, values_from = mean_expr)
expr_pop <- tidyr::pivot_wider(expr_mean, names_from = population, values_from = c(control, heat))


#################### Figure 4A: front-loading of shared heat response genes

# For each pair of populations: baseline difference = difference in control expression; plasticity difference = difference in heat-induced change.
# Front-loading predicts a negative relationship: the population with higher baseline expression shows a smaller induction.
fl_NS <- data.frame(gene = expr_pop$gene, comp = "North vs South", baseline = expr_pop$control_North - expr_pop$control_South,
                    plasticity = (expr_pop$heat_North - expr_pop$control_North) - (expr_pop$heat_South - expr_pop$control_South))
fl_CN <- data.frame(gene = expr_pop$gene, comp = "Central vs North", baseline = expr_pop$control_Central - expr_pop$control_North,
                    plasticity = (expr_pop$heat_Central - expr_pop$control_Central) - (expr_pop$heat_North - expr_pop$control_North))
fl_CS <- data.frame(gene = expr_pop$gene, comp = "Central vs South", baseline = expr_pop$control_Central - expr_pop$control_South,
                    plasticity = (expr_pop$heat_Central - expr_pop$control_Central) - (expr_pop$heat_South - expr_pop$control_South))
frontload <- rbind(fl_NS, fl_CS, fl_CN)
frontload <- frontload[order(match(frontload$gene, expr_pop$gene), match(frontload$comp, c("North vs South", "Central vs South", "Central vs North"))), ]
frontload$comp <- factor(frontload$comp, levels = c("North vs South", "Central vs South", "Central vs North"))
frontload$type <- ifelse(frontload$gene %in% hsp70_genes, "Hsp70", "Other genes")

# Number of genes in each quadrant (Q2 and Q4: higher baseline goes with lower plasticity).
frontload$quadrant <- dplyr::case_when(frontload$baseline > 0 & frontload$plasticity > 0 ~ "Q1", frontload$baseline < 0 & frontload$plasticity > 0 ~ "Q2",
                                       frontload$baseline < 0 & frontload$plasticity < 0 ~ "Q3", frontload$baseline > 0 & frontload$plasticity < 0 ~ "Q4")
table(frontload$comp, frontload$quadrant)

fig4A <- ggplot(frontload) + geom_vline(xintercept = 0, color = "black", linetype = "dashed") + geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(aes(x = baseline, y = plasticity, shape = type, size = type, fill = comp, color = comp), stroke = 0.6) +
  geom_smooth(aes(x = baseline, y = plasticity, color = comp, fill = comp), method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.1) +
  labs(x = "Baseline expression difference", y = "Plasticity difference", color = "Population comparison", fill = "Population comparison", shape = "Gene category") +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"), axis.title = element_text(color = "black", size = 16),
        axis.text = element_text(color = "black", size = 14), legend.title = element_text(size = 15), legend.text = element_text(size = 14),
        legend.position = "inside", legend.position.inside = c(0.2, 0.185), legend.spacing.y = unit(0.05, "cm"), plot.title = element_text(size = 11, hjust = 0.5)) +
  scale_color_manual(values = c("North vs South" = "#1B9E77", "Central vs South" = "#7570B3", "Central vs North" = "#E7298A"), labels = c("North vs South" = "North vs. South", "Central vs South" = "Central vs. South", "Central vs North" = "Central vs. North")) +
  scale_fill_manual(values = c("North vs South" = scales::alpha("#1B9E77", 0.65), "Central vs South" = scales::alpha("#7570B3", 0.65), "Central vs North" = scales::alpha("#E7298A", 0.65)), labels = c("North vs South" = "North vs. South", "Central vs South" = "Central vs. South", "Central vs North" = "Central vs. North")) +
  scale_size_manual(values = c("Hsp70" = 5, "Other genes" = 3)) + scale_shape_manual(values = c("Hsp70" = 24, "Other genes" = 21)) +
  guides(size = "none", color = guide_legend(order = 2), fill = guide_legend(order = 2), shape = guide_legend(order = 1, override.aes = list(size = c("Hsp70" = 3, "Other genes" = 2))))
fig4A


#################### Figure 4B: baseline expression and plasticity in each population
# Baseline expression = mean log2 CPM under control; plasticity = mean log2 CPM under long heat stress - baseline expression
comparison_colors <- c("North vs South" = "#1B9E77", "Central vs South" = "#7570B3", "Central vs North" = "#E7298A")
population_colors <- c("North" = "#87aade", "Central" = "#ffe57f", "South" = "#ff9999")
expr_classes <- expr_mean
expr_classes$plasticity <- expr_classes$heat - expr_classes$control
expr_classes$population <- factor(expr_classes$population, levels = c("North", "Central", "South"))
expr_classes$type <- ifelse(expr_classes$gene %in% hsp70_genes, "Hsp70", "Other genes")
hsp70_expr <- expr_classes[expr_classes$type == "Hsp70", ]
other_expr <- expr_classes[expr_classes$type == "Other genes", ]

# Pairwise Wilcoxon tests between populations (BH correction over the three pairs), separately for each measure and gene class
population_tests <- function(dat, value, gene_type) {
  test <- pairwise.wilcox.test(x = dat[[value]], g = dat$population, paired = FALSE, p.adjust.method = "BH")
  result <- as.data.frame(as.table(test$p.value))
  names(result) <- c("group2", "group1", "p.adj")
  result <- result[!is.na(result$p.adj), ]
  result$group1 <- as.character(result$group1)
  result$group2 <- as.character(result$group2)
  result$comparison <- case_when(
    result$group1 == "North" & result$group2 == "South" ~ "North vs South",
    result$group1 == "Central" & result$group2 == "South" ~ "Central vs South",
    result$group1 == "North" & result$group2 == "Central" ~ "Central vs North")
  result$bracket_color <- unname(comparison_colors[result$comparison])
  result$label <- case_when(result$p.adj < 0.0001 ~ "****", result$p.adj < 0.001 ~ "***",
                            result$p.adj < 0.01 ~ "**", result$p.adj < 0.05 ~ "*", TRUE ~ "ns")
  result$gene_class <- gene_type
  result$metric <- value
  result$n_genes <- length(unique(dat$gene))
  list(test = test, table = result)
}
test_hsp_baseline <- population_tests(hsp70_expr, "control", "Hsp70")
test_hsp_plasticity <- population_tests(hsp70_expr, "plasticity", "Hsp70")
test_other_baseline <- population_tests(other_expr, "control", "Other genes")
test_other_plasticity <- population_tests(other_expr, "plasticity", "Other genes")
all_tests <- bind_rows(test_hsp_baseline$table, test_hsp_plasticity$table,
                       test_other_baseline$table, test_other_plasticity$table)

# Same y range for Hsp70 and other genes, with room for the significance brackets
common_y_limits <- setNames(lapply(c("control", "plasticity"), function(value) {
  brackets <- unlist(lapply(c("Hsp70", "Other genes"), function(gene_type) {
    n_sig <- sum(all_tests$metric == value & all_tests$gene_class == gene_type & all_tests$p.adj < 0.05)
    if (!n_sig) return(numeric(0))
    max(expr_classes[[value]][expr_classes$type == gene_type]) + 0.65 + (seq_len(n_sig) - 1) * 1.8
  }))
  c(floor(min(c(0, expr_classes[[value]]))), ceiling(max(c(expr_classes[[value]], brackets))) + 1)
}), c("control", "plasticity"))

# Boxplots: baseline expression and plasticity in rows, Hsp70 and other genes in columns
box_data <- tidyr::pivot_longer(expr_classes, cols = c(control, plasticity),
                                names_to = "metric", values_to = "value")
box_data$metric <- factor(box_data$metric, levels = c("control", "plasticity"))
box_data$type <- factor(box_data$type, levels = c("Hsp70", "Other genes"))
box_sig <- all_tests %>% dplyr::filter(p.adj < 0.05) %>%
  dplyr::group_by(gene_class, metric) %>%
  dplyr::mutate(y.position = max(box_data$value[box_data$type == dplyr::first(gene_class) &
                                                  box_data$metric == dplyr::first(metric)]) + 0.65 + (dplyr::row_number() - 1) * 1.8) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(type = factor(gene_class, levels = c("Hsp70", "Other genes")),
                metric = factor(metric, levels = c("control", "plasticity")))
# Raise the brackets and stars slightly
box_sig$y.position <- box_sig$y.position + 0.6
box_limits <- expand.grid(metric = c("control", "plasticity"), type = c("Hsp70", "Other genes"),
                          bound = 1:2, stringsAsFactors = FALSE)
box_limits$value <- mapply(function(m, b) common_y_limits[[m]][b], box_limits$metric, box_limits$bound)
box_limits$population <- factor("North", levels = levels(expr_classes$population))
box_limits$metric <- factor(box_limits$metric, levels = c("control", "plasticity"))
box_limits$type <- factor(box_limits$type, levels = c("Hsp70", "Other genes"))
metric_strip_labels <- c(control = "Baseline expression", plasticity = "Plasticity")
fig4B_core <- ggplot(box_data, aes(x = population, y = value)) +
  geom_blank(data = box_limits) +
  geom_boxplot(aes(fill = population), width = 0.7, color = "grey20", alpha = 1,
               outlier.shape = NA, linewidth = 0.5) +
  geom_point(aes(fill = population, shape = type), color = "grey35",
             position = position_jitter(width = 0.10, height = 0, seed = 20261022),
             size = 1.6, stroke = 0.4, alpha = 0.8) +
  facet_grid(metric ~ type, scales = "free_y",
             labeller = labeller(type = c(Hsp70 = "Hsp70", "Other genes" = "Other genes"),
                                 metric = metric_strip_labels)) +
  stat_pvalue_manual(box_sig, label = "label", y.position = "y.position", color = "black",
                     tip.length = 0.01, bracket.size = 0.45, size = 7, fontface = "bold", hide.ns = TRUE, inherit.aes = FALSE) +
  scale_y_continuous(breaks = scales::breaks_width(5), minor_breaks = NULL,
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_fill_manual(values = population_colors) + scale_shape_manual(values = c(Hsp70 = 24, "Other genes" = 21)) +
  labs(x = NULL, y = NULL) + theme_bw(base_size = 16, base_family = "sans") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey96"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
        panel.spacing = grid::unit(0, "pt"), strip.placement = "outside",
        strip.background.x = element_rect(fill = "#D3D3D3", color = "black", linewidth = 0.5),
        strip.background.y = element_rect(fill = "#D3D3D3", color = "black", linewidth = 0.5),
        strip.text.x = element_text(size = 16, face = "bold"),
        strip.text.y.right = element_text(size = 16, face = "bold", angle = 270),
        axis.text.x = element_text(size = 16, face = "bold", color = "black"),
        axis.text.y = element_text(size = 14, color = "black"), legend.position = "none", plot.margin = margin(6, 7, 5, 4))
# Y-axis unit of each row on the left
fig4B_table <- ggplotGrob(fig4B_core)
fig4B_panel_rows <- sort(unique(fig4B_table$layout$t[grepl("^panel($|-)", fig4B_table$layout$name)]))
fig4B_table <- gtable::gtable_add_cols(fig4B_table, grid::unit(24, "pt"), pos = 0)
row_y_units <- list(expression(log[2]~CPM), expression(log[2]~FC))
for (i in seq_along(fig4B_panel_rows)) {
  fig4B_table <- gtable::gtable_add_grob(fig4B_table,
                                         grid::textGrob(row_y_units[[i]], rot = 90, gp = grid::gpar(fontsize = 16, fontfamily = "sans")),
                                         t = fig4B_panel_rows[i], b = fig4B_panel_rows[i], l = 1, r = 1, clip = "off", name = paste0("row_y_unit_", i))
}
fig4B <- cowplot::ggdraw() + cowplot::draw_grob(fig4B_table, y = 0.01, height = 0.96)
fig4B

#################### Combine Figure 4A-B
fig4 <- ggarrange(fig4A, fig4B, labels = c("A", "B"), ncol = 2,
                  widths = c(7, 6.5), font.label = list(size = 20, face = "bold"))
fig4
ggsave("Figure4.pdf", fig4, width = 13.5, height = 7.4, device = cairo_pdf, bg = "white")
