setwd("C:/path/to/RNAgithub/09Hsp70")

# =============================================================================
# Figure S12: relative copy number of the 21 Hsp genes among the 43 shared heat response genes, from pool-seq depth
# Input: Hsp_depth.txt, Housekeeper_depth.txt, housekeeper_original.gtf
# Output: FigureS12.pdf
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# Pool-seq depth at each exon position of the 21 Hsp genes in the six pools (positions without reads have depth 0).
depth <- read.delim("Hsp_depth.txt", check.names = FALSE)
# Panel labels (gene ID and Hsp family), in the order of the input file.
facet_order <- unique(paste0("<span style='font-size:10pt'><b>", depth$gene_id, "</b></span><br><i>", depth$Hsp_type, "</i>"))
samples <- c("Norway.high", "Norway.low", "Sweden.high", "Sweden.low", "Spain.high", "Spain.low")

# Remove exon positions with unusually high depth in Hsp90 (V1264_003735) and MSTRG.27436.
depth$max_depth <- apply(depth[samples], 1, max)
depth$median_depth <- apply(depth[samples], 1, median)
depth <- depth %>% dplyr::filter(ref_gene_id != "V1264_003735" | max_depth <= 220)
depth <- depth %>% dplyr::filter(gene_id != "MSTRG.27436" | (median_depth <= 180 & max_depth <= 250))

# Mean exon depth of each gene in each pool.
sample_summary <- depth %>%
  dplyr::group_by(gene_id, Hsp_type) %>%
  dplyr::summarise(dplyr::across(all_of(samples), mean), .groups = "drop") %>%
  tidyr::pivot_longer(cols = all_of(samples), names_to = "sample", values_to = "mean_exon_depth")

# Depth of the two single-copy housekeeping genes (RPL13A g15582, HPRT1 g13887), as in Figure5_S13-18.R Part 8.
hk_depth <- read.delim("Housekeeper_depth.txt", check.names = FALSE)
pools <- c("Northhigh", "Northlow", "Centralhigh", "Centrallow", "Southhigh", "Southlow")
hk <- read.delim("housekeeper_original.gtf", header = FALSE, quote = "", comment.char = "#")
colnames(hk) <- c("Chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Convert genomic exon positions to positions on the two reference gene sequences.
hk_pos <- hk %>%
  dplyr::filter(feature == "exon") %>%
  dplyr::mutate(origin = dplyr::case_when(Chr == "LG3_SUPER_4" ~ 48883704, Chr == "LG12_SUPER_3" ~ 54491592), Chr = dplyr::recode(Chr, "LG3_SUPER_4" = "g15582.t1", "LG12_SUPER_3" = "g13887.t1"), start = start - origin + 1, end = end - origin + 1) %>%
  dplyr::rowwise() %>% dplyr::mutate(pos = list(seq.int(start, end))) %>%
  tidyr::unnest(pos) %>% dplyr::ungroup() %>% dplyr::distinct(Chr, pos)

# Single-copy depth of each pool = mean exon depth of each housekeeping gene, averaged over the two genes.
hk_gene_depth <- hk_depth %>%
  dplyr::semi_join(hk_pos, by = c("Chr", "pos")) %>%
  tidyr::pivot_longer(cols = all_of(pools), names_to = "pool", values_to = "depth") %>%
  dplyr::group_by(Chr, pool) %>% dplyr::summarise(mean_depth = mean(depth), .groups = "drop")
norm_factor <- hk_gene_depth %>%
  dplyr::group_by(pool) %>% dplyr::summarise(factor = mean(mean_depth), .groups = "drop") %>%
  dplyr::mutate(sample = dplyr::recode(pool, "Northhigh" = "Norway.high", "Northlow" = "Norway.low", "Centralhigh" = "Sweden.high", "Centrallow" = "Sweden.low", "Southhigh" = "Spain.high", "Southlow" = "Spain.low")) %>%
  dplyr::select(sample, factor)
norm_factor

# Relative copy number = Hsp mean exon depth / housekeeping reference depth.
sample_summary <- sample_summary %>%
  dplyr::left_join(norm_factor, by = "sample") %>%
  dplyr::mutate(relative_copy_number = mean_exon_depth / factor, country = sub("\\..*", "", sample))

# Points: high-shore (triangle) and low-shore (circle) pools.
point_data <- sample_summary %>%
  dplyr::mutate(country = dplyr::recode(country, "Norway" = "North", "Sweden" = "Central", "Spain" = "South"), country = factor(country, levels = c("North", "Central", "South")), shore = factor(sub(".*\\.", "", sample), levels = c("high", "low")), facet_label = paste0("<span style='font-size:10pt'><b>", gene_id, "</b></span><br><i>", Hsp_type, "</i>"))
point_data$facet_label <- factor(point_data$facet_label, levels = facet_order)

# Bars: mean of the high- and low-shore pools.
plot_data <- point_data %>%
  dplyr::group_by(gene_id, Hsp_type, country, facet_label) %>%
  dplyr::summarise(relative_copy_number = mean(relative_copy_number), .groups = "drop")

# Figure S12: one panel per gene; dashed line = relative copy number of 1.
country_colors <- c("North" = "#87aade", "Central" = "#ffe57f", "South" = "#ff9999")
y_max <- ceiling(max(point_data$relative_copy_number) * 5) / 5

p <- ggplot(plot_data, aes(x = country, y = relative_copy_number, fill = country)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey35", linewidth = 0.4) +
  geom_point(data = point_data, aes(shape = shore, group = shore), position = position_dodge(width = 0.3), size = 2, stroke = 0.5, colour = "black") +
  facet_wrap(~ facet_label, ncol = 5) +
  scale_fill_manual(values = country_colors, guide = "none") +
  scale_shape_manual(values = c("high" = 24, "low" = 21), labels = c("High shore", "Low shore"), name = NULL) +
  scale_y_continuous(limits = c(0, y_max), breaks = seq(0, y_max, 0.5), expand = expansion(mult = c(0, 0))) +
  guides(shape = guide_legend(ncol = 1, override.aes = list(fill = "white", size = 2.5))) +
  labs(x = "", y = "Relative copy number") +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(axis.title.y = element_text(size = 11.5, face = "plain"), axis.text = element_text(size = 10, colour = "black"), axis.text.x = element_text(angle = 0, hjust = 0.5), strip.text = ggtext::element_markdown(size = 10, lineheight = 0.85, margin = margin(3.5, 2, 2, 2)), panel.grid.minor = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.30, 0.08), legend.justification = c(0.5, 0.5), legend.direction = "vertical", legend.text = element_text(size = 10), legend.key.size = grid::unit(0.45, "cm"), legend.background = element_blank())

p
ggsave("FigureS12.pdf", p, width = 21, height = 21, units = "cm", device = cairo_pdf)
