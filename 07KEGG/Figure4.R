# =============================================================================
# KEGG annotation and KEGG pathway enrichment of general temperature response genes (Figure 4)
#
# Part 1 assigns KEGG orthologs (K numbers) to every gene. Transcripts were annotated with EggNOG-mapper and KofamKOALA;
# K numbers from both sources are pooled and transferred from transcripts to genes using the StringTie assembly, and
# each K number is linked to its KEGG pathways.
#
# Part 2 tests which KEGG pathways are enriched among the general temperature response genes (up- and down-regulated
# separately) in each population.
#
# Part 3 shows the enrichment of the four significant pathways in the three populations (Figure 4A) and how the up-regulated
# genes of each pathway are shared among populations, coloured by the GO clusters of Figure 3 (Figure 4B).
#
# Input:  ko_orthology_Lsaxegg.txt (EggNOG-mapper), ko_orhtology_Lsaxkofam.txt (KofamKOALA), ko_list_db.txt (K number descriptions,
#         KofamKOALA), gene_to_transcript_mapping_StringTie.txt, gene_go_description_all.txt, GO_cluster_terms.txt (from 05GO),
#         counts_filter.txt, North/Central/South_generaltemp_allgenes.txt (from 04DESeq2_Figure2.R)
# Output: gene_ko_description_all.txt   K numbers and KEGG pathways of each gene
#         Figure4.pdf
# =============================================================================

setwd("C:/Users/hanbq/Desktop/RNA20161011/07KEGG/")
library(dplyr)
library(tidyr)
library(KEGGREST)
library(clusterProfiler)
library(grid)


############################################################
#Part 1: KEGG annotation of genes

# EggNOG-mapper: K numbers of each transcript (several K numbers separated by commas).
egg <- read.table("ko_orthology_Lsaxegg.txt", sep = " ", header = FALSE, stringsAsFactors = FALSE, col.names = c("transcript", "K"))
egg$transcript <- sub("^gnl\\|WGS:JBAMIC\\|", "", egg$transcript)
egg$transcript <- sub("\\.p[0-9]*$", "", egg$transcript)
egg$K <- gsub("ko:", "", egg$K)

# KofamKOALA: K number of each transcript.
kofam <- read.table("ko_orhtology_Lsaxkofam.txt", sep = " ", header = FALSE, stringsAsFactors = FALSE, col.names = c("transcript", "K"))
kofam$transcript <- sub("^gnl\\|WGS:JBAMIC\\|", "", kofam$transcript)
kofam$transcript <- sub("\\.p[0-9]*$", "", kofam$transcript)


# K numbers at the gene level: all unique K numbers of a gene's transcripts from both sources, one gene-K pair per row.
stringtie <- read.table("../05GO/gene_to_transcript_mapping_StringTie.txt", sep = "\t", header = TRUE)
eggkofam <- dplyr::left_join(rbind(egg, kofam), stringtie, by = "transcript")
dd <- tidyr::separate_rows(eggkofam[!is.na(eggkofam$gene), c("gene", "K")], K, sep = ",")
dd$K <- gsub(" ", "", dd$K)
dd <- unique(dd[!is.na(dd$K) & dd$K != "NA" & dd$K != "", ])

# Name of each K number.
k_description <- read.table("ko_list_db.txt", sep = "\t", header = TRUE, fill = TRUE, quote = "")
k_description <- k_description[, c(1, 12)]
colnames(k_description) <- c("K", "K_description")
gene_K <- dplyr::left_join(dd, k_description, by = "K")

# KEGG pathways of every K number and the pathway names (KEGG database, queried with KEGGREST).
k_pathway <- keggLink("pathway", "ko")
k_to_pathway <- data.frame(K = sub("ko:", "", names(k_pathway)), pathway = sub("path:", "", k_pathway))
k_to_pathway <- k_to_pathway[grepl("^ko", k_to_pathway$pathway), ]
pathway_names <- keggList("pathway")
pathway_description <- data.frame(pathway = sub("^(path:)?map", "ko", names(pathway_names)), pathway_description = as.vector(pathway_names))
k_to_pathway <- dplyr::left_join(k_to_pathway, pathway_description, by = "pathway")

# Pathways that do not apply to a snail are not used: Human Diseases (ko05), Drug Development (ko07), and a manually selected
# list of Organismal Systems pathways (mainly vertebrate immune, endocrine, circulatory, digestive, nervous and sensory systems)
# and drug resistance pathways.
excluded_pathways <- c("ko04640", "ko04610", "ko04611", "ko04613", "ko04620", "ko04624", "ko04621", "ko04622", "ko04623", "ko04625", "ko04650", "ko04612", "ko04660",
                       "ko04658", "ko04659", "ko04657", "ko04662", "ko04664", "ko04666", "ko04670", "ko04672", "ko04062", "ko04911", "ko04910", "ko04922", "ko04923",
                       "ko04920", "ko03320", "ko04929", "ko04912", "ko04913", "ko04915", "ko04914", "ko04917", "ko04921", "ko04926", "ko04935", "ko04918", "ko04919",
                       "ko04928", "ko04916", "ko04924", "ko04614", "ko04925", "ko04927", "ko04260", "ko04261", "ko04270", "ko04970", "ko04971", "ko04972", "ko04976",
                       "ko04973", "ko04974", "ko04975", "ko04979", "ko04977", "ko04981", "ko04980", "ko04978", "ko04962", "ko04960", "ko04961", "ko04964", "ko04966",
                       "ko04724", "ko04727", "ko04725", "ko04728", "ko04726", "ko04720", "ko04730", "ko04723", "ko04721", "ko04722", "ko04744", "ko04745", "ko04740",
                       "ko04742", "ko04750", "ko04320", "ko04360", "ko04361", "ko04380", "ko04211", "ko04212", "ko04213", "ko04710", "ko04713", "ko04711", "ko04712",
                       "ko04714", "ko04626", "ko04930", "ko04940", "ko04950", "ko04936", "ko04932", "ko04931", "ko04933", "ko04934", "ko01501", "ko01502", "ko01503",
                       "ko01521", "ko01524", "ko01523", "ko01522")
k_to_pathway <- k_to_pathway[!grepl("^ko05|^ko07", k_to_pathway$pathway) & !k_to_pathway$pathway %in% excluded_pathways, ]

# Save the gene-K-pathway table used for KEGG enrichment (genes whose K numbers have no retained pathway keep NA as pathway).
gene_ko <- merge(gene_K, k_to_pathway, by = "K", all.x = TRUE)[, c("gene", "K", "K_description", "pathway", "pathway_description")]
write.table(gene_ko, file = "gene_ko_description_all.txt", sep = "\t", quote = F, row.names = F, col.names = T)


############################################################
#Part 2: KEGG enrichment of general temperature response genes

# Gene-pathway annotation from Part 1 (genes with at least one retained KEGG pathway).
gene_ko <- read.delim("gene_ko_description_all.txt", header = TRUE, stringsAsFactors = FALSE, quote = "")
gene_pathway <- dplyr::distinct(gene_ko[!is.na(gene_ko$pathway), c("gene", "pathway", "pathway_description")])
TERM2GENE <- dplyr::distinct(gene_pathway[, c("pathway", "gene")])
TERM2NAME <- dplyr::distinct(gene_pathway[, c("pathway", "pathway_description")])

# General temperature response of all retained genes in each population (from 04DESeq2_Figure2.R).
deg_North <- read.delim("../04DEseq2/North_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
deg_Central <- read.delim("../04DEseq2/Central_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
deg_South <- read.delim("../04DEseq2/South_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)

# Background for enrichment: all retained genes (final filtered count matrix) that belong to at least one KEGG pathway.
counts <- read.table("../01counts_ind_filter/counts_filter.txt", header = TRUE)
universe <- intersect(rownames(counts), TERM2GENE$gene)
length(universe)


#################### KEGG enrichment of up- and down-regulated general temperature response genes
# A pathway is significantly enriched when p.adjust <= 0.01 and qvalue <= 0.01.

#################### North
up_North <- deg_North$gene[which(deg_North$padj < 0.05 & deg_North$log2FoldChange >= 1)]
kegg_up_North <- enricher(gene = up_North, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
                          minGSSize = 10, maxGSSize = 3000, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_up_North <- kegg_up_North@result
sig_up_North <- res_up_North$ID[res_up_North$p.adjust <= 0.01 & res_up_North$qvalue <= 0.01]

down_North <- deg_North$gene[which(deg_North$padj < 0.05 & deg_North$log2FoldChange <= -1)]
kegg_down_North <- enricher(gene = down_North, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
                            minGSSize = 10, maxGSSize = 3000, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_down_North <- kegg_down_North@result
sig_down_North <- res_down_North$ID[res_down_North$p.adjust <= 0.01 & res_down_North$qvalue <= 0.01]

#################### Central
up_Central <- deg_Central$gene[which(deg_Central$padj < 0.05 & deg_Central$log2FoldChange >= 1)]
kegg_up_Central <- enricher(gene = up_Central, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
                            minGSSize = 10, maxGSSize = 3000, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_up_Central <- kegg_up_Central@result
sig_up_Central <- res_up_Central$ID[res_up_Central$p.adjust <= 0.01 & res_up_Central$qvalue <= 0.01]

down_Central <- deg_Central$gene[which(deg_Central$padj < 0.05 & deg_Central$log2FoldChange <= -1)]
kegg_down_Central <- enricher(gene = down_Central, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
                              minGSSize = 10, maxGSSize = 3000, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_down_Central <- kegg_down_Central@result
sig_down_Central <- res_down_Central$ID[res_down_Central$p.adjust <= 0.01 & res_down_Central$qvalue <= 0.01]

#################### South
up_South <- deg_South$gene[which(deg_South$padj < 0.05 & deg_South$log2FoldChange >= 1)]
kegg_up_South <- enricher(gene = up_South, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
                          minGSSize = 10, maxGSSize = 3000, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_up_South <- kegg_up_South@result
sig_up_South <- res_up_South$ID[res_up_South$p.adjust <= 0.01 & res_up_South$qvalue <= 0.01]

down_South <- deg_South$gene[which(deg_South$padj < 0.05 & deg_South$log2FoldChange <= -1)]
kegg_down_South <- enricher(gene = down_South, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
                            minGSSize = 10, maxGSSize = 3000, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_down_South <- kegg_down_South@result
sig_down_South <- res_down_South$ID[res_down_South$p.adjust <= 0.01 & res_down_South$qvalue <= 0.01]

# Pathways significantly enriched in at least one population: four among up-regulated genes (ko04141 Protein processing in
# endoplasmic reticulum, ko04010 MAPK signaling pathway, ko03010 Ribosome, ko00190 Oxidative phosphorylation); none among
# down-regulated genes.
unique(c(sig_up_North, sig_up_Central, sig_up_South))
unique(c(sig_down_North, sig_down_Central, sig_down_South))


############################################################
#Part 3: Figure 4

# The four pathways significantly enriched among up-regulated genes, in the order of the figure.
pathways_up <- c("ko04141", "ko04010", "ko03010", "ko00190")
pathway_labels <- c("Protein processing in\nendoplasmic reticulum", "MAPK signaling pathway", "Ribosome", "Oxidative phosphorylation")

#################### Data of Figure 4A: gene count and FDR of each pathway in each population
fig4A <- data.frame(pathway = rep(1:4, 3), population = rep(c("North", "Central", "South"), each = 4),
                    Count = c(res_up_North$Count[match(pathways_up, res_up_North$ID)], res_up_Central$Count[match(pathways_up, res_up_Central$ID)],
                              res_up_South$Count[match(pathways_up, res_up_South$ID)]),
                    p.adjust = c(res_up_North$p.adjust[match(pathways_up, res_up_North$ID)], res_up_Central$p.adjust[match(pathways_up, res_up_Central$ID)],
                                 res_up_South$p.adjust[match(pathways_up, res_up_South$ID)]),
                    qvalue = c(res_up_North$qvalue[match(pathways_up, res_up_North$ID)], res_up_Central$qvalue[match(pathways_up, res_up_Central$ID)],
                               res_up_South$qvalue[match(pathways_up, res_up_South$ID)]))
fig4A$significant <- !is.na(fig4A$p.adjust) & fig4A$p.adjust <= 0.01 & fig4A$qvalue <= 0.01
fig4A$log10FDR <- -log10(fig4A$p.adjust)
fig4A
# Pathways not tested in a population (no up-regulated gene in the pathway) have no bubble.
fig4A <- fig4A[!is.na(fig4A$Count), ]

#################### GO cluster of the genes in each pathway (clusters of Figure 3, from 05GO)
gene_GO <- read.delim("../05GO/gene_go_description_all.txt", header = TRUE, stringsAsFactors = FALSE)
GO_cluster_terms <- read.delim("../05GO/GO_cluster_terms.txt", header = TRUE, stringsAsFactors = FALSE)
gene_cluster <- unique(merge(gene_GO[, c("gene", "go")], GO_cluster_terms[, c("GO_ID", "cluster")], by.x = "go", by.y = "GO_ID")[, c("gene", "cluster")])
cluster_colours <- c("Heat response genes" = "#FDBF6F", "Ribosome-related processes" = "#B2DF8A", "Mitochondrial respiration" = "#80B1D3",
                     "Tf2-like" = "#CAB2D6", "Unclustered" = "#B3B3B3")

# Up-regulated genes (in any population) of each pathway, and the GO clusters they belong to. Each pathway is coloured by its
# most common GO cluster: Heat response genes (Protein processing in ER, MAPK), Ribosome-related processes (Ribosome),
# Mitochondrial respiration (Oxidative phosphorylation).
up_any <- unique(c(up_North, up_Central, up_South))
table(gene_cluster$cluster[gene_cluster$gene %in% intersect(TERM2GENE$gene[TERM2GENE$pathway == "ko04141"], up_any)])
table(gene_cluster$cluster[gene_cluster$gene %in% intersect(TERM2GENE$gene[TERM2GENE$pathway == "ko04010"], up_any)])
table(gene_cluster$cluster[gene_cluster$gene %in% intersect(TERM2GENE$gene[TERM2GENE$pathway == "ko03010"], up_any)])
table(gene_cluster$cluster[gene_cluster$gene %in% intersect(TERM2GENE$gene[TERM2GENE$pathway == "ko00190"], up_any)])
main_cluster <- c("Heat response genes", "Heat response genes", "Ribosome-related processes", "Mitochondrial respiration")

#################### Data of Figure 4B: up-regulated genes of each pathway shared among populations (UpSet)
# Up-regulated genes of each pathway (in any population). Each gene is assigned to exactly one combination of populations
# (e.g. "110" = up-regulated in North and Central only); genes of the pathway's main GO cluster and other genes are counted
# separately. Only combinations with genes are shown, ordered by size.
genes_ER <- unique(c(intersect(up_North, TERM2GENE$gene[TERM2GENE$pathway == "ko04141"]), intersect(up_Central, TERM2GENE$gene[TERM2GENE$pathway == "ko04141"]),
                     intersect(up_South, TERM2GENE$gene[TERM2GENE$pathway == "ko04141"])))
genes_MAPK <- unique(c(intersect(up_North, TERM2GENE$gene[TERM2GENE$pathway == "ko04010"]), intersect(up_Central, TERM2GENE$gene[TERM2GENE$pathway == "ko04010"]),
                       intersect(up_South, TERM2GENE$gene[TERM2GENE$pathway == "ko04010"])))
genes_Ribo <- unique(c(intersect(up_North, TERM2GENE$gene[TERM2GENE$pathway == "ko03010"]), intersect(up_Central, TERM2GENE$gene[TERM2GENE$pathway == "ko03010"]),
                       intersect(up_South, TERM2GENE$gene[TERM2GENE$pathway == "ko03010"])))
genes_OXPHOS <- unique(c(intersect(up_North, TERM2GENE$gene[TERM2GENE$pathway == "ko00190"]), intersect(up_Central, TERM2GENE$gene[TERM2GENE$pathway == "ko00190"]),
                         intersect(up_South, TERM2GENE$gene[TERM2GENE$pathway == "ko00190"])))
fig4B_genes <- data.frame(panel = rep(1:4, c(length(genes_ER), length(genes_MAPK), length(genes_Ribo), length(genes_OXPHOS))),
                          gene = c(genes_ER, genes_MAPK, genes_Ribo, genes_OXPHOS))
fig4B_genes$combination <- paste0(as.integer(fig4B_genes$gene %in% up_North), as.integer(fig4B_genes$gene %in% up_Central), as.integer(fig4B_genes$gene %in% up_South))
fig4B_genes$in_main <- paste(fig4B_genes$gene, main_cluster[fig4B_genes$panel]) %in% paste(gene_cluster$gene, gene_cluster$cluster)
fig4B <- fig4B_genes %>% dplyr::group_by(panel, combination) %>% dplyr::summarise(Count = dplyr::n(), n_main = sum(in_main), .groups = "drop") %>%
  dplyr::arrange(panel, -Count, match(combination, c("111", "110", "101", "011", "100", "010", "001"))) %>% as.data.frame()
fig4B


#################### Figure 4 (drawn with grid on a 170 x 97 mm page; positions in mm, font sizes in points)
cairo_pdf("Figure4.pdf", width = 170 / 25.4, height = 97 / 25.4, family = "Arial", bg = "white")
grid.newpage()
# All positions below are given on a 104 mm high layout; it is shifted down by 7 mm to remove the empty space at the bottom.
pushViewport(viewport(y = unit(-7, "mm"), just = "bottom", width = unit(170, "mm"), height = unit(104, "mm")))

#### A: bubbles of the four pathways (rows) in the three populations (columns)
pop_x <- c(19.33, 40, 60.67)                    # centres of the North / Central / South columns (from x = 9 to 71 mm)
pathway_y <- c(90, 84, 78, 72)                  # rows of the four pathways
fdr_palette <- colorRampPalette(c("#AADBB4", "#16813D"))(256)
fdr_limit <- max(10, ceiling(max(fig4A$log10FDR[fig4A$significant]) / 10) * 10)
fig4A$x <- pop_x[match(fig4A$population, c("North", "Central", "South"))]
fig4A$y <- pathway_y[fig4A$pathway]
fig4A$radius <- 1.8 * sqrt(fig4A$Count / 80)    # bubble area proportional to the gene count
fig4A$fill <- ifelse(fig4A$significant, fdr_palette[1 + round(255 * pmin(fig4A$log10FDR / fdr_limit, 1))], "#CCCCCC")

grid.text("A", x = unit(2.5, "mm"), y = unit(100, "mm"), gp = gpar(fontsize = 9, fontface = "bold"))
grid.rect(x = unit(pop_x, "mm"), y = unit(96, "mm"), width = unit(20.07, "mm"), height = unit(4.6, "mm"), gp = gpar(fill = "#E6E6E6", col = "black", lwd = 0.35))
grid.text(c("North", "Central", "South"), x = unit(pop_x, "mm"), y = unit(96, "mm"), gp = gpar(fontsize = 7))
grid.segments(x0 = unit(9, "mm"), x1 = unit(71, "mm"), y0 = unit(pathway_y, "mm"), y1 = unit(pathway_y, "mm"), gp = gpar(col = "#DEDEDE", lwd = 0.5))
grid.text(sub("\n", " ", pathway_labels), x = unit(75, "mm"), y = unit(pathway_y, "mm"), just = "left", gp = gpar(fontsize = 7, fontface = "bold"))
grid.circle(x = unit(fig4A$x, "mm"), y = unit(fig4A$y, "mm"), r = unit(fig4A$radius, "mm"), gp = gpar(fill = fig4A$fill, col = NA))

# Legend of A: gene count (circles) and -log10(FDR) (colour bar), at the top right. Same format as the GO cluster legend:
# titles 6 pt (aligned by their top, because the subscript of -log10 makes that title taller), labels 5 pt, first item 5 mm below the title, items 4.5 mm apart.
legend_counts <- c(20, 40, 60, 80)
legend_y <- c(91.5, 87, 82.5, 78)
grid.text("Gene count", x = unit(134, "mm"), y = unit(97.5, "mm"), just = c("left", "top"), gp = gpar(fontsize = 6))
grid.circle(x = unit(135.8, "mm"), y = unit(legend_y, "mm"), r = unit(1.8 * sqrt(legend_counts / 80), "mm"), gp = gpar(fill = NA, col = "black", lwd = 0.4))
grid.text(legend_counts, x = unit(138.8, "mm"), y = unit(legend_y, "mm"), just = "left", gp = gpar(fontsize = 5))
grid.text(expression(-log[10]*"(FDR)"), x = unit(148, "mm"), y = unit(97.5, "mm"), just = c("left", "top"), gp = gpar(fontsize = 6))
bar_bottom <- 78
bar_top <- 91.5
grid.rect(x = unit(150, "mm"), y = unit(bar_bottom + (1:256 - 0.5) * (bar_top - bar_bottom) / 256, "mm"), width = unit(1.8, "mm"),
          height = unit((bar_top - bar_bottom) / 256 + 0.02, "mm"), gp = gpar(fill = fdr_palette, col = NA))
tick_y <- bar_bottom + c(0, 0.5, 1) * (bar_top - bar_bottom)
grid.segments(x0 = unit(150.9, "mm"), x1 = unit(151.4, "mm"), y0 = unit(tick_y, "mm"), y1 = unit(tick_y, "mm"), gp = gpar(lwd = 0.3))
grid.text(c(0, fdr_limit / 2, fdr_limit), x = unit(151.9, "mm"), y = unit(tick_y, "mm"), just = "left", gp = gpar(fontsize = 5))

#### B: one UpSet panel per pathway (bars = genes in each combination of populations; dots below = the combination)
upset_ymax <- max(85, ceiling(max(fig4B$Count) / 10) * 10)   # common bar scale for the four panels
panel_width <- 34.5
panel_left <- 18 + (0:3) * (panel_width + 3)
set_y <- c(27.5, 23.3, 19.1)                    # rows of North / Central / South in the dot matrix
bar_baseline <- 32
fig4B$k <- ave(fig4B$Count, fig4B$panel, FUN = seq_along)
fig4B$n <- ave(fig4B$Count, fig4B$panel, FUN = length)
fig4B$spacing <- pmin(7.8, (panel_width - 2) / fig4B$n)
fig4B$x <- panel_left[fig4B$panel] + (panel_width - fig4B$n * fig4B$spacing) / 2 + (fig4B$k - 0.5) * fig4B$spacing
fig4B$h_main <- fig4B$n_main / upset_ymax * 33
fig4B$h_other <- (fig4B$Count - fig4B$n_main) / upset_ymax * 33
fig4B$in_N <- substr(fig4B$combination, 1, 1) == "1"
fig4B$in_C <- substr(fig4B$combination, 2, 2) == "1"
fig4B$in_S <- substr(fig4B$combination, 3, 3) == "1"

grid.text("B", x = unit(2.5, "mm"), y = unit(63, "mm"), gp = gpar(fontsize = 9, fontface = "bold"))
grid.text(c("North", "Central", "South"), x = unit(16.5, "mm"), y = unit(set_y, "mm"), just = "right", gp = gpar(fontsize = 7))
grid.rect(x = unit(rep(panel_left + panel_width / 2, 2), "mm"), y = unit(rep(set_y[c(1, 3)], each = 4), "mm"), width = unit(panel_width, "mm"),
          height = unit(3.8, "mm"), gp = gpar(fill = "#EEEEEE", col = NA))

# Bars: genes of the main GO cluster (bottom, cluster colour) and other genes (top, grey), with the total on top.
bars_main <- fig4B[fig4B$h_main > 0, ]
bars_other <- fig4B[fig4B$h_other > 0, ]
grid.rect(x = unit(bars_main$x, "mm"), y = unit(bar_baseline, "mm"), width = unit(bars_main$spacing * 0.84, "mm"), height = unit(bars_main$h_main, "mm"),
          just = "bottom", gp = gpar(fill = cluster_colours[main_cluster[bars_main$panel]], col = "black", lwd = 0.5))
grid.rect(x = unit(bars_other$x, "mm"), y = unit(bar_baseline + bars_other$h_main, "mm"), width = unit(bars_other$spacing * 0.84, "mm"),
          height = unit(bars_other$h_other, "mm"), just = "bottom", gp = gpar(fill = "#D9D9D9", col = "black", lwd = 0.5))
grid.text(fig4B$Count, x = unit(fig4B$x, "mm"), y = unit(bar_baseline + fig4B$h_main + fig4B$h_other + 1.8, "mm"), gp = gpar(fontsize = 7))

# Dot matrix: black = population in the combination; a vertical line joins the populations of a shared combination.
grid.circle(x = unit(rep(fig4B$x, 3), "mm"), y = unit(rep(set_y, each = nrow(fig4B)), "mm"), r = unit(0.6, "mm"),
            gp = gpar(fill = ifelse(c(fig4B$in_N, fig4B$in_C, fig4B$in_S), "#282828", "#C8C8C8"), col = NA))
shared <- fig4B[fig4B$in_N + fig4B$in_C + fig4B$in_S > 1, ]
grid.segments(x0 = unit(shared$x, "mm"), x1 = unit(shared$x, "mm"), y0 = unit(ifelse(shared$in_N, set_y[1], set_y[2]), "mm"),
              y1 = unit(ifelse(shared$in_S, set_y[3], set_y[2]), "mm"), gp = gpar(lwd = 1))
grid.text(pathway_labels, x = unit(panel_left + panel_width / 2, "mm"), y = unit(14.6, "mm"), just = "top", gp = gpar(fontsize = 7, fontface = "bold", lineheight = 0.98))

# Hsp70: the genes shared by all three populations in Protein processing in ER and MAPK signaling include Hsp70 genes
# (label and arrows in the colour of the Heat response genes cluster).
hsp70 <- fig4B[fig4B$panel %in% 1:2 & fig4B$combination == "111", ]
hsp70$target_x <- hsp70$x + ifelse(hsp70$panel == 1, 1, -1) * (hsp70$spacing * 0.42 + 0.8)
hsp70$target_y <- bar_baseline + hsp70$h_main + hsp70$h_other - pmin(0.9, (hsp70$h_main + hsp70$h_other) / 2)
root_x <- mean(hsp70$target_x)
root_y <- max(hsp70$target_y) + 8
grid.text("Hsp70", x = unit(root_x, "mm"), y = unit(root_y + 2.8, "mm"), gp = gpar(fontsize = 8, fontface = "bold", col = cluster_colours["Heat response genes"]))
grid.segments(x0 = unit(root_x, "mm"), y0 = unit(root_y, "mm"), x1 = unit(hsp70$target_x, "mm"), y1 = unit(hsp70$target_y, "mm"),
              arrow = arrow(length = unit(1.4, "mm"), type = "closed"), gp = gpar(col = cluster_colours["Heat response genes"], fill = cluster_colours["Heat response genes"], lty = 2, lwd = 0.65))

#### Legend of the GO clusters (right side, below the legend of A and above the Oxidative phosphorylation panel)
go_legend_y <- c(64.5, 60, 55.5, 51)
grid.text("GO cluster", x = unit(134, "mm"), y = unit(70.5, "mm"), just = c("left", "top"), gp = gpar(fontsize = 6))
grid.rect(x = unit(134.9, "mm"), y = unit(go_legend_y, "mm"), width = unit(1.8, "mm"), height = unit(1.8, "mm"),
          gp = gpar(fill = c(cluster_colours[c("Heat response genes", "Ribosome-related processes", "Mitochondrial respiration")], "#D9D9D9"), col = "black", lwd = 0.3))
grid.text(c("Heat response genes", "Ribosome-related processes", "Mitochondrial respiration", "Other genes"), x = unit(137, "mm"), y = unit(go_legend_y, "mm"),
          just = "left", gp = gpar(fontsize = 5))
dev.off()

