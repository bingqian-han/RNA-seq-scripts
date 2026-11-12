setwd("C:/path/to/RNAgithub/06GOKEGG/")

# =============================================================================
# GO annotation and GO enrichment of general temperature response genes
#
# Part 1 builds a functional annotation for every gene. Transcripts were annotated against SwissProt and NR
# (DIAMOND), InterProScan and FANTASIA; GO terms from all four sources are pooled and transferred from transcripts
# to genes using the StringTie assembly. Protein descriptions from NR, SwissProt and InterProScan are combined
# into one description per gene.
#
# Part 2 tests which GO terms are enriched among the general temperature response genes (up- and down-regulated
# separately) in each population, and combines the up-regulated GO terms of the three populations.
#
# Part 3 clusters the up-regulated GO terms by shared genes into four functional categories and compares these categories
# and their genes among populations (Figure 3A-B), and tests KEGG enrichment of the general temperature response genes
# (Figure 3C-D).
#
# Part 4 shows the gene count and FDR of every up-regulated GO term in each population, ordered by the GO tree
# (Figure S7).
#
# Input:  SwissProt_GO.xlsx, SwissProt_annotation.tab, NR_GO.xlsx, NR_annotation.tab, InterProScan_GO.tsv,
#         InterProScan_annotation.tsv, FANTASIA_results.csv, Lsax_stringtie_merged.gtf, go.db.txt,
#         North/Central/South_generaltemp_allgenes.txt (from 04DESeq2_Figure2.R), counts_filter.txt
#         ko_orthology_Lsaxegg.txt, ko_orhtology_Lsaxkofam.txt, ko_list_db.txt,
#         heat_intersect_43genes_withannotation_Manual_curation.txt (manual curation of gene families)
# Output: gene_to_transcript_mapping_StringTie.txt  gene-transcript relationships
#         gene_go_description_all.txt              GO terms of each gene
#         gene_annotation.txt                      functional description of each gene
#         heat_intersect_annotation_for_curation.txt  shared heat response genes, for manual curation
#         GO_cluster_terms.txt                     functional cluster of each up-regulated GO term
#         gene_ko_description_all.txt, Figure3.pdf, FigureS7.pdf
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(enrichplot)
library(ggtree)
library(ggplot2)
library(tibble)
library(patchwork)
library(ggpubr)
library(ggVennDiagram)
library(edgeR)
library(grid)

############################################################ 
#Part 1: functional annotation of genes

#################### GO annotation of genes, combining four annotation sources

# SwissProt: GO terms of the SwissProt proteins (from UniProt ID mapping).
swissgo <- read_excel("SwissProt_GO.xlsx")
swissgo <- swissgo[, c(1, 6)]
colnames(swissgo) <- c("swissid", "GO")
swissgo$GO <- gsub(" ", "", swissgo$GO)

# SwissProt: DIAMOND search results.
swiss_matches <- read.delim("SwissProt_annotation.tab", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
swissgene <- swiss_matches
swissgene$SubjectID <- gsub(".*\\|(.*)\\|.*", "\\1", swissgene$SubjectID)
swissgene$QueryID <- sub(".*\\|", "", swissgene$QueryID)
swissgene <- swissgene[, 1:2]
swissgene$QueryID <- sub("\\.p[0-9]*$", "", swissgene$QueryID)
colnames(swissgene) <- c("transcript", "swissid")
swissgene$transcript <- gsub(" ", "", swissgene$transcript)

# SwissProt: GO terms of each transcript, through its matched protein.
nrow(swissgene)
nrow(swissgo)
merged_swiss <- dplyr::left_join(swissgene, swissgo, by = "swissid")
nrow(merged_swiss)
merged_swiss <- merged_swiss[, c(1, 3)]
merged_swiss <- merged_swiss[!is.na(merged_swiss$GO), ]
nrow(merged_swiss)
merged_swiss <- merged_swiss[merged_swiss$GO != "", ]

# NR: GO terms of the NR proteins (from UniProt ID mapping).
nrgo <- read_excel("NR_GO.xlsx")
nrow(nrgo)
nrgo <- nrgo[, c(1, 6)]
colnames(nrgo) <- c("nrid", "GO")
nrgo$GO <- gsub(" ", "", nrgo$GO)
nrgo <- nrgo[!is.na(nrgo$GO), ]
nrgo <- nrgo[nrgo$GO != "", ]
length(nrgo$nrid)
length(unique(nrgo$nrid))

# NR: DIAMOND search results.
nr_matches <- read.delim("NR_annotation.tab", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
nrgene <- nr_matches[, c("QueryID", "SubjectID")]
colnames(nrgene) <- c("transcript", "nrid")
nrgene$transcript <- sub(".*\\|", "", nrgene$transcript)
nrgene$transcript <- sub("\\.p[0-9]*$", "", nrgene$transcript)

# NR: GO terms of each transcript, through its matched protein.
nrow(nrgene)
nrow(nrgo)
merged_nr <- dplyr::left_join(nrgene, nrgo, by = "nrid")
merged_nr <- merged_nr[, c(1, 3)]
merged_nr <- merged_nr[!is.na(merged_nr$GO), ]
merged_nr <- merged_nr[merged_nr$GO != "", ]

# InterProScan: GO terms assigned to each transcript from its protein domains.
inter <- read.table("InterProScan_GO.tsv", sep = "\t", header = FALSE)
inter$V1 <- sub("\\.p[0-9]*$", "", inter$V1)
colnames(inter) <- c("transcript", "GO")

# FANTASIA: GO terms predicted for each transcript.
fantasia <- read.csv("FANTASIA_results.csv", sep = ",", header = TRUE)
fantasia <- fantasia[, c("query_accession", "go_id")]
colnames(fantasia) <- c("transcript", "GO")
fantasia$transcript <- sub("\\.p[0-9]*$", "", fantasia$transcript)
fantasia$transcript <- gsub("gnl\\|WGS:JBAMIC\\|", "", fantasia$transcript)
fantasia$transcript <- gsub(" ", "", fantasia$transcript)
fantasia$GO <- gsub(" ", "", fantasia$GO)
fantasia <- fantasia[!is.na(fantasia$GO), ]
fantasia <- fantasia[fantasia$GO != "", ]
fantasia <- fantasia[fantasia$GO != "-", ]

# Pool the GO terms of each transcript from all four sources.
innrsw <- rbind(merged_nr, merged_swiss, inter, fantasia)

# Gene-to-transcript relationships from the merged StringTie assembly (also used for the functional annotation below).
stringtie <- read.table("Lsax_stringtie_merged.gtf", sep = "\t", header = FALSE)
stringtie <- stringtie[stringtie$V3 == "transcript", ]
stringtie$V9 <- gsub("gene_id ", "", stringtie$V9)
stringtie$V9 <- gsub("transcript_id ", "", stringtie$V9)
stringtie$V9 <- gsub("gnl\\|WGS:JBAMIC\\|", "", stringtie$V9)
stringtie$V9 <- gsub("ref.* ", "", stringtie$V9)
stringtie$V9 <- gsub(" ", "", stringtie$V9)
stringtie$V10 <- unlist(lapply(stringtie$V9, function(x) unlist(strsplit(x, ";"))[1]))
stringtie$V11 <- unlist(lapply(stringtie$V9, function(x) unlist(strsplit(x, ";"))[2]))
stringtie <- stringtie[, 10:11]
colnames(stringtie) <- c("gene", "transcript")
write.table(stringtie, file = "gene_to_transcript_mapping_StringTie.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# GO annotation at the gene level: all unique GO terms of a gene's transcripts, one gene-GO pair per row.
innrsw <- dplyr::left_join(innrsw, stringtie, by = "transcript")
innrsw <- innrsw[, c(3, 1, 2)]
dd <- data.frame(gene = sort(na.omit(unique(innrsw$gene))))
dd$go <- tapply(innrsw$GO, innrsw$gene, function(x) {
  y <- unique(unlist(strsplit(x, ";")))
  y <- y[y != "NA"]
  paste(y, collapse = ";")
})
dd <- dd[!is.na(dd$go), ]
dd <- tidyr::separate_rows(dd, go, sep = ";")
dd$go <- gsub(" ", "", dd$go)

# Add the name of each GO term and save the gene-GO table used for GO enrichment.
go.description <- read.table("go.db.txt", sep = "\t", header = FALSE)
colnames(go.description) <- c("go", "description")
go <- dplyr::left_join(dd, go.description, by = "go")
write.table(go, file = "gene_go_description_all.txt", sep = "\t", quote = F, row.names = F, col.names = T)


#################### Functional description of genes, combining three annotation sources

# NR: protein descriptions of each gene's transcripts.
nr.anno <- nr_matches[, c(1, 14)]
colnames(nr.anno) <- c("transcript", "description")
nr.anno$transcript <- sub("^.*\\|", "", nr.anno$transcript)
nr.anno$transcript <- sub("\\.p[0-9]+$", "", nr.anno$transcript)
nr.anno$source <- "NR"
nr.anno <- unique(nr.anno)
nr.gene <- merge(stringtie, nr.anno, by = "transcript")
nr.gene <- nr.gene[, c("gene", "transcript", "description", "source")]
nr.gene.merge <- aggregate(description ~ gene, data = nr.gene, FUN = function(x) paste(unique(x), collapse = "; "))
nr.gene.merge$source <- "NR"

# SwissProt: protein descriptions of each gene's transcripts.
swiss.anno <- swiss_matches[, c(1, 14)]
colnames(swiss.anno) <- c("transcript", "description")
swiss.anno$transcript <- sub("^.*\\|", "", swiss.anno$transcript)
swiss.anno$transcript <- sub("\\.p[0-9]+$", "", swiss.anno$transcript)
swiss.anno$source <- "SwissProt"
swiss.anno <- unique(swiss.anno)
swiss.gene <- merge(stringtie, swiss.anno, by = "transcript")
swiss.gene <- swiss.gene[, c("gene", "transcript", "description", "source")]
swiss.gene <- unique(swiss.gene)
swiss.gene.merge <- aggregate(description ~ gene, data = swiss.gene, FUN = function(x) paste(unique(x), collapse = "; "))
swiss.gene.merge$source <- "SwissProt"

# InterProScan: protein domain descriptions of each gene's transcripts.
interpro.anno <- read.delim("InterProScan_annotation.tsv", sep = "\t", header = FALSE)
colnames(interpro.anno) <- c("transcript", "description")
interpro.anno$transcript <- sub("^.*\\|", "", interpro.anno$transcript)
interpro.anno$transcript <- sub("\\.p[0-9]+$", "", interpro.anno$transcript)
interpro.anno$source <- "InterProScan"
interpro.anno <- interpro.anno[!is.na(interpro.anno$description) & interpro.anno$description != "" & interpro.anno$description != "-", ]
interpro.anno <- unique(interpro.anno)
interpro.gene <- merge(stringtie, interpro.anno, by = "transcript")
interpro.gene <- interpro.gene[, c("gene", "transcript", "description", "source")]
interpro.gene <- unique(interpro.gene)
interpro.gene.merge <- aggregate(description ~ gene, data = interpro.gene, FUN = function(x) paste(unique(x), collapse = "; "))
interpro.gene.merge$source <- "InterProScan"

# One row per gene: all descriptions and the sources they came from.
all.anno <- rbind(nr.gene.merge, swiss.gene.merge, interpro.gene.merge)
all.gene.annotation <- aggregate(description ~ gene, data = all.anno, FUN = function(x) paste(unique(x), collapse = "; "))
all.gene.source <- aggregate(source ~ gene, data = all.anno, FUN = function(x) paste(unique(x), collapse = "; "))
all.gene.annotation <- merge(all.gene.annotation, all.gene.source, by = "gene")
all.gene.annotation <- all.gene.annotation[, c("gene", "description", "source")]
write.table(all.gene.annotation, file = "gene_annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)


############################################################ 
#Part 2: GO enrichment of general temperature response genes
# Up- and down-regulated genes: general temperature response genes (padj < 0.05, |log2FC| >= 0.5) with positive or
# negative log2FC. The same genes are used for the KEGG enrichment in Part 3.

# Gene-GO annotation from Part 1.
gene_GO <- read.delim("gene_go_description_all.txt", header = TRUE, stringsAsFactors = FALSE)
gene_GO$description <- gsub("^ ", "", gene_GO$description)

# General temperature response of all retained genes in each population (from 04DESeq2_Figure2.R).
deg_North <- read.delim("../05DEseq2/North_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
deg_Central <- read.delim("../05DEseq2/Central_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
deg_South <- read.delim("../05DEseq2/South_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)

# Background for enrichment: all retained genes (final filtered count matrix) that have at least one GO term.
counts <- read.table("../01counts_ind_filter/counts_filter.txt", header = TRUE)
universe <- intersect(rownames(counts), gene_GO$gene)


#################### GO enrichment of up- and down-regulated general temperature response genes
# A GO term is significantly enriched when p.adjust <= 0.01 and qvalue <= 0.01.

#################### North

# Up-regulated genes, and how many of them have GO annotation.
up_North <- deg_North$gene[which(deg_North$padj < 0.05 & deg_North$log2FoldChange >= 0.5)]
length(up_North)
sum(up_North %in% gene_GO$gene)
go_up_North <- enricher(gene = up_North, universe = universe, TERM2GENE = gene_GO[c("go", "gene")], TERM2NAME = gene_GO[c("go", "description")],
                        minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_up_North <- go_up_North@result
sig_up_North <- res_up_North$ID[res_up_North$p.adjust <= 0.01 & res_up_North$qvalue <= 0.01]
length(sig_up_North)

# Down-regulated genes.
down_North <- deg_North$gene[which(deg_North$padj < 0.05 & deg_North$log2FoldChange <= -0.5)]
length(down_North)
sum(down_North %in% gene_GO$gene)
go_down_North <- enricher(gene = down_North, universe = universe, TERM2GENE = gene_GO[c("go", "gene")], TERM2NAME = gene_GO[c("go", "description")],
                          minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_down_North <- go_down_North@result
sum(res_down_North$p.adjust <= 0.01 & res_down_North$qvalue <= 0.01)



#################### Central

# Up-regulated genes, and how many of them have GO annotation.
up_Central <- deg_Central$gene[which(deg_Central$padj < 0.05 & deg_Central$log2FoldChange >= 0.5)]
length(up_Central)
sum(up_Central %in% gene_GO$gene)
go_up_Central <- enricher(gene = up_Central, universe = universe, TERM2GENE = gene_GO[c("go", "gene")], TERM2NAME = gene_GO[c("go", "description")],
                          minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_up_Central <- go_up_Central@result
sig_up_Central <- res_up_Central$ID[res_up_Central$p.adjust <= 0.01 & res_up_Central$qvalue <= 0.01]
length(sig_up_Central)

# Down-regulated genes.
down_Central <- deg_Central$gene[which(deg_Central$padj < 0.05 & deg_Central$log2FoldChange <= -0.5)]
length(down_Central)
sum(down_Central %in% gene_GO$gene)
go_down_Central <- enricher(gene = down_Central, universe = universe, TERM2GENE = gene_GO[c("go", "gene")], TERM2NAME = gene_GO[c("go", "description")],
                            minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_down_Central <- go_down_Central@result
sum(res_down_Central$p.adjust <= 0.01 & res_down_Central$qvalue <= 0.01)

#################### South

# Up-regulated genes, and how many of them have GO annotation.
up_South <- deg_South$gene[which(deg_South$padj < 0.05 & deg_South$log2FoldChange >= 0.5)]
length(up_South)
sum(up_South %in% gene_GO$gene)
go_up_South <- enricher(gene = up_South, universe = universe, TERM2GENE = gene_GO[c("go", "gene")], TERM2NAME = gene_GO[c("go", "description")],
                        minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_up_South <- go_up_South@result
sig_up_South <- res_up_South$ID[res_up_South$p.adjust <= 0.01 & res_up_South$qvalue <= 0.01]
length(sig_up_South)

# Down-regulated genes.
down_South <- deg_South$gene[which(deg_South$padj < 0.05 & deg_South$log2FoldChange <= -0.5)]
length(down_South)
sum(down_South %in% gene_GO$gene)
go_down_South <- enricher(gene = down_South, universe = universe, TERM2GENE = gene_GO[c("go", "gene")], TERM2NAME = gene_GO[c("go", "description")],
                          minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_down_South <- go_down_South@result
sum(res_down_South$p.adjust <= 0.01 & res_down_South$qvalue <= 0.01)


#################### Merge up-regulated GO terms across populations

# All GO terms significantly enriched among up-regulated genes in at least one population.
go_up_union <- union(union(sig_up_North, sig_up_Central), sig_up_South)
length(go_up_union)

# For each of these GO terms, gene ratio and adjusted p-value in every population, taken from the full enrichment results,
# so a term that is significant in one population still shows its (non-significant) values in the others.
# NA means the term was not tested in that population (no up-regulated gene annotated to it).
cmp_North <- res_up_North[res_up_North$ID %in% go_up_union, c("ID", "GeneRatio", "p.adjust")]
cmp_Central <- res_up_Central[res_up_Central$ID %in% go_up_union, c("ID", "GeneRatio", "p.adjust")]
cmp_South <- res_up_South[res_up_South$ID %in% go_up_union, c("ID", "GeneRatio", "p.adjust")]
colnames(cmp_North) <- c("ID", "GeneRatio_North", "padj_North")
colnames(cmp_Central) <- c("ID", "GeneRatio_Central", "padj_Central")
colnames(cmp_South) <- c("ID", "GeneRatio_South", "padj_South")

# One row per GO term, with its name.
go_up_names <- data.frame(ID = sort(go_up_union), Description = gene_GO$description[match(sort(go_up_union), gene_GO$go)])
go_up_compare <- dplyr::left_join(dplyr::left_join(dplyr::left_join(go_up_names, cmp_North, by = "ID"), cmp_Central, by = "ID"), cmp_South, by = "ID")


############################################################
#Part 3: GO functional clusters (Figure 3A-B) and KEGG enrichment/overlap (Figure 3C-D)

#################### Figure 3A: clustering of the up-regulated GO terms and their significance across the three populations

# Build an enrichment object containing the GO terms significantly enriched among up-regulated genes in at least one
# population (all cutoffs set to 1 so no term is dropped), so their similarity (shared genes) can be used to cluster
# them into functional groups.
genes_go_up <- unique(gene_GO$gene[gene_GO$go %in% go_up_compare$ID])
gene_GO_up <- gene_GO[gene_GO$go %in% go_up_compare$ID, ]
go_up_all <- enricher(gene = genes_go_up, TERM2GENE = gene_GO_up[c("go", "gene")], TERM2NAME = gene_GO_up[c("go", "description")],
                      minGSSize = 1, maxGSSize = Inf, pvalueCutoff = 1, qvalueCutoff = 1)
nrow(go_up_all@result)
go_up_sim <- pairwise_termsim(go_up_all)

# Tree of GO terms (average linkage), keeping all enriched terms.
N_GO <- nrow(go_up_all@result)
tree_base <- treeplot(go_up_sim, showCategory = N_GO, cluster.params = list(method = "average", n = 3, label_words_n = 0),
                      hilight.params = list(hilight = FALSE))
tree_base$layers <- tree_base$layers[1:2]

# Branches of the four functional clusters in this 79-GO tree (node numbers refer to this particular tree):
# rRNA = 87, heat response = 86, Tf2-like = 88, mitochondrial respiration = 83.
cluster_nodes <- c(rRNA = 87, Heat = 86, TF2 = 88, Mito = 83)
cluster_definitions <- data.frame(group = as.character(1:4), cluster = c("rRNA", "Heat", "TF2", "Mito"), fill = c("#B2DF8A", "#FDBF6F", "#CAB2D6", "#A6CEE3"),
                                  text_color = c("#3A6831", "#C46A00", "#6A3D9A", "#1F78B4"),
                                  title = c("Ribosome-related processes", "Heat response genes", "Tf2-like\nretrotransposons", "Mitochondrial\nrespiration"), stringsAsFactors = FALSE)
branch_colors <- c("0" = "black", setNames(cluster_definitions$fill, cluster_definitions$group))

# Colour the four clusters and rotate two branches so the unclustered GO terms sit between the clusters (topology unchanged).
tree_ordered <- groupClade(tree_base, unname(cluster_nodes), group_name = "group")
tree_ordered <- ggtree::rotate(tree_ordered, 84)
tree_ordered <- ggtree::rotate(tree_ordered, 81)
tree_go <- tree_ordered
tree_go$data$y <- min(tree_go$data$y) + max(tree_go$data$y) - tree_go$data$y

# GO terms of each cluster, from left to right in Figure 3A.
tree_tips <- dplyr::arrange(dplyr::filter(tree_go$data, isTip), y)
tree_tips$GO <- go_up_all@result$ID[match(tree_tips$label, go_up_all@result$Description)]
go_rRNA <- tree_tips$GO[tree_tips$group == "1"]
go_heat <- tree_tips$GO[tree_tips$group == "2"]
go_TF2 <- tree_tips$GO[tree_tips$group == "3"]
go_mito <- tree_tips$GO[tree_tips$group == "4"]

# Cluster of every up-regulated GO term (Unclustered = outside the four clusters), in the order of Figure 3A.
GO_cluster_terms <- data.frame(GO_ID = tree_tips$GO, GO_term_name = tree_tips$label,
                               cluster = c("0" = "Unclustered", "1" = "Ribosome-related processes", "2" = "Heat response genes", "3" = "Tf2-like retrotransposons",
                                           "4" = "Mitochondrial respiration")[as.character(tree_tips$group)])
write.table(GO_cluster_terms, "GO_cluster_terms.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Number of GO terms of each cluster that are significantly enriched (p.adjust <= 0.01 and qvalue <= 0.01) in each population.
sig_count <- data.frame(population = c("North", "Central", "South"),
                        rRNA = c(sum(go_rRNA %in% sig_up_North), sum(go_rRNA %in% sig_up_Central), sum(go_rRNA %in% sig_up_South)),
                        Heat = c(sum(go_heat %in% sig_up_North), sum(go_heat %in% sig_up_Central), sum(go_heat %in% sig_up_South)),
                        TF2 = c(sum(go_TF2 %in% sig_up_North), sum(go_TF2 %in% sig_up_Central), sum(go_TF2 %in% sig_up_South)),
                        Mito = c(sum(go_mito %in% sig_up_North), sum(go_mito %in% sig_up_Central), sum(go_mito %in% sig_up_South)))
sig_count$population <- factor(sig_count$population, levels = c("South", "Central", "North"))
sig_count

# Layout of Figure 3: figure width, left margin of the population labels (panels A and D) and spacing of the population rows
figure3_width_in <- 13.5
population_label_left_pt <- 10
population_order <- c("North", "Central", "South")
population_row_spacing_A_pt <- 34
population_row_spacing_D_pt <- 26
population_band_height_pt <- 22

# Vertical layout of panel A (in points)
panel_padding_A_pt <- 6
panel_header_height_A_pt <- 46
population_y_A <- setNames(panel_padding_A_pt + population_band_height_pt / 2 +
                             c(2, 1, 0) * population_row_spacing_A_pt, population_order)
panel_body_top_A <- max(population_y_A) + population_band_height_pt / 2
panel_top_A <- panel_body_top_A + panel_header_height_A_pt
panel_title_y_A <- panel_body_top_A + panel_header_height_A_pt / 2
tree_base_y_A <- panel_top_A + 10
tree_height_A <- 140
y_limits_A <- c(-40, tree_base_y_A + tree_height_A + 10)
figure3_height_A <- grid::convertHeight(grid::unit(diff(y_limits_A) + 8 + 5, "pt"), "inches", valueOnly = TRUE)


# Horizontal position of each GO term in Figure 3A; small clusters are widened so that their count panel and title fit
runs <- rle(as.character(tree_tips$group))
min_run_width <- setNames(rep(10, nrow(cluster_definitions)), cluster_definitions$group)
min_run_width[cluster_definitions$group[cluster_definitions$cluster == "Mito"]] <- 13
run_width <- ifelse(runs$values == "0", runs$lengths, pmax(runs$lengths, min_run_width[runs$values]))
run_right <- 0.5 + cumsum(run_width)
run_left <- run_right - run_width
cursor <- max(run_right)
x_limits_A <- c(-6, cursor + 0.2)
tree_tips$run_left <- rep(run_left, runs$lengths)
tree_tips$run_right <- rep(run_right, runs$lengths)
tree_tips$x_A <- tree_tips$run_left + (sequence(runs$lengths) - 0.5) / rep(runs$lengths, runs$lengths) * rep(run_width, runs$lengths)

# Internal nodes of the tree are placed above the middle of their children (filled from the tips towards the root).
nodes_A <- as.data.frame(tree_go$data)
nodes_A$x_plot <- NA_real_
nodes_A$x_plot[match(tree_tips$node, nodes_A$node)] <- tree_tips$x_A
while (anyNA(nodes_A$x_plot)) {
  for (i in which(is.na(nodes_A$x_plot))) {
    children <- which(nodes_A$parent == nodes_A$node[i] & nodes_A$node != nodes_A$node[i])
    if (all(!is.na(nodes_A$x_plot[children]))) nodes_A$x_plot[i] <- mean(nodes_A$x_plot[children])
  }
}
nodes_A$y_plot <- tree_base_y_A + (max(nodes_A$x) - nodes_A$x) / max(nodes_A$x) * tree_height_A
edges_A <- nodes_A[nodes_A$node != nodes_A$parent, ]
edges_A$x_parent <- nodes_A$x_plot[match(edges_A$parent, nodes_A$node)]
edges_A$y_parent <- nodes_A$y_plot[match(edges_A$parent, nodes_A$node)]
edges_A$branch_color <- unname(branch_colors[as.character(edges_A$group)])

# Count panel of each cluster, spanning the GO terms of that cluster
panels_A <- dplyr::summarise(dplyr::group_by(dplyr::filter(tree_tips, as.character(group) != "0"), group), xmin = min(run_left) + 0.35, xmax = max(run_right) - 0.35, .groups = "drop")
panels_A <- dplyr::left_join(dplyr::mutate(panels_A, group = as.character(group)), cluster_definitions, by = "group")
panels_A$bar_xmin <- panels_A$xmin
panels_A$bar_xmax <- panels_A$xmax - 1.2

# Same bar scale in all four count panels
requested_limit_A <- c(rRNA = max(sig_count$rRNA) + 3, Heat = max(sig_count$Heat) + 3, TF2 = max(sig_count$TF2) + 3, Mito = max(sig_count$Mito) + 3)[panels_A$cluster]
count_unit_A <- min((panels_A$bar_xmax - panels_A$bar_xmin) / requested_limit_A)
panels_A$limit <- floor((panels_A$bar_xmax - panels_A$bar_xmin) / count_unit_A + 1e-9)
panels_A$mid <- (panels_A$xmin + panels_A$xmax) / 2

# Bars and axis ticks of the count panels.
bars_A <- tidyr::pivot_longer(sig_count, -population, names_to = "cluster", values_to = "Count")
bars_A <- dplyr::left_join(bars_A, panels_A[, c("cluster", "bar_xmin", "bar_xmax", "fill", "limit")], by = "cluster")
bars_A$y <- unname(population_y_A[as.character(bars_A$population)])
bars_A$xend <- bars_A$bar_xmin + bars_A$Count * count_unit_A
limit_rRNA <- panels_A$limit[panels_A$cluster == "rRNA"]
limit_heat <- panels_A$limit[panels_A$cluster == "Heat"]
limit_TF2 <- panels_A$limit[panels_A$cluster == "TF2"]
limit_mito <- panels_A$limit[panels_A$cluster == "Mito"]
tick_rRNA <- pretty(c(0, limit_rRNA), n = 3)
tick_heat <- pretty(c(0, limit_heat), n = 3)
tick_TF2 <- pretty(c(0, limit_TF2), n = 3)
tick_mito <- pretty(c(0, limit_mito), n = 3)
ticks_A <- rbind(data.frame(cluster = "rRNA", label = tick_rRNA[tick_rRNA <= limit_rRNA]), data.frame(cluster = "Heat", label = tick_heat[tick_heat <= limit_heat]),
                 data.frame(cluster = "TF2", label = tick_TF2[tick_TF2 <= limit_TF2]), data.frame(cluster = "Mito", label = tick_mito[tick_mito <= limit_mito]))
ticks_A <- dplyr::left_join(ticks_A, panels_A[, c("cluster", "bar_xmin")], by = "cluster")
ticks_A$x <- ticks_A$bar_xmin + ticks_A$label * count_unit_A
population_labels_A <- data.frame(population = population_order, y = unname(population_y_A))

fig3A <- ggplot() +
  geom_segment(data = edges_A, aes(x = x_parent, xend = x_plot, y = y_parent, yend = y_parent), color = edges_A$branch_color, linewidth = 0.85, lineend = "butt") +
  geom_segment(data = edges_A, aes(x = x_plot, xend = x_plot, y = y_parent, yend = y_plot), color = edges_A$branch_color, linewidth = 0.85, lineend = "butt") +
  geom_rect(data = panels_A, aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = panel_top_A, fill = "white", color = panels_A$fill, linewidth = 1.1) +
  geom_rect(data = bars_A, aes(xmin = bar_xmin, xmax = xend, ymin = y - population_band_height_pt / 2, ymax = y + population_band_height_pt / 2), fill = bars_A$fill, color = NA) +
  geom_text(data = bars_A, aes(x = xend, y = y, label = Count), hjust = -0.25, size = 14, size.unit = "pt") +
  geom_text(data = panels_A, aes(x = mid, label = title), y = panel_title_y_A, color = panels_A$text_color, fontface = "bold", size = 17,
            size.unit = "pt", lineheight = 0.95) +
  geom_text(data = population_labels_A, aes(y = y, label = population), x = x_limits_A[1], hjust = 0, fontface = "bold", size = 17, size.unit = "pt") +
  geom_segment(data = ticks_A, aes(x = x, xend = x), y = 0, yend = -3, color = "grey40", linewidth = 0.3) +
  geom_text(data = ticks_A, aes(x = x, label = label), y = -12, size = 14, size.unit = "pt", color = "grey30") +
  annotate("text", x = (min(panels_A$xmin) + max(panels_A$xmax)) / 2, y = -30, label = "Number of significantly enriched GO terms", size = 16, size.unit = "pt") +
  scale_x_continuous(limits = x_limits_A, expand = c(0, 0)) + scale_y_continuous(limits = y_limits_A, expand = c(0, 0)) +
  coord_cartesian(clip = "off") + theme_void(base_family = "sans") + theme(plot.margin = margin(8, 10, 5, population_label_left_pt))
fig3A


#################### Figure 3B: overlap among populations of the up-regulated genes behind each GO cluster

# Up-regulated genes annotated to any GO term of each cluster, within the enrichment background.
# Figure 3A counts significant GO terms; Figure 3B counts unique genes, including those in non-significant terms of the same cluster.
genes_heat <- unique(gene_GO$gene[gene_GO$go %in% go_heat & gene_GO$gene %in% universe])
genes_rRNA <- unique(gene_GO$gene[gene_GO$go %in% go_rRNA & gene_GO$gene %in% universe])
genes_mito <- unique(gene_GO$gene[gene_GO$go %in% go_mito & gene_GO$gene %in% universe])
genes_TF2 <- unique(gene_GO$gene[gene_GO$go %in% go_TF2 & gene_GO$gene %in% universe])
genes_heat_North <- intersect(up_North, genes_heat)
genes_heat_Central <- intersect(up_Central, genes_heat)
genes_heat_South <- intersect(up_South, genes_heat)
genes_rRNA_North <- intersect(up_North, genes_rRNA)
genes_rRNA_Central <- intersect(up_Central, genes_rRNA)
genes_rRNA_South <- intersect(up_South, genes_rRNA)
genes_mito_North <- intersect(up_North, genes_mito)
genes_mito_Central <- intersect(up_Central, genes_mito)
genes_mito_South <- intersect(up_South, genes_mito)
genes_TF2_North <- intersect(up_North, genes_TF2)
genes_TF2_Central <- intersect(up_Central, genes_TF2)
genes_TF2_South <- intersect(up_South, genes_TF2)

# Venn diagrams of the four clusters (font sizes in points are divided by .pt because ggVennDiagram uses mm)
venn_horizontal_padding <- 1 / 25
venn_plot_margin <- margin(16, 2, 6, 2)

# Venn diagram: heat response genes.
venn_heat <- ggVennDiagram(list(North = genes_heat_North, Central = genes_heat_Central, South = genes_heat_South), label_alpha = 0, label_size = 14 / .pt, set_size = 16 / .pt, edge_size = 0.65) +
  scale_fill_gradientn(colours = c("#FEF8F1", "#FDEBD4", "#FDBF6F")) + scale_x_continuous(expand = expansion(mult = venn_horizontal_padding)) +
  theme_void(base_family = "sans") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17, color = "#C46A00", lineheight = 0.95), legend.position = "none",
                                           plot.margin = venn_plot_margin) + coord_fixed(clip = "off")
venn_heat$layers[[3]]$aes_params$fontface <- "bold"
venn_heat$layers[[3]]$data$X <- c(-1.5, 5.5, 2)

# Venn diagram: ribosome-related processes.
venn_rRNA <- ggVennDiagram(list(North = genes_rRNA_North, Central = genes_rRNA_Central, South = genes_rRNA_South), label_alpha = 0, label_size = 14 / .pt, set_size = 16 / .pt, edge_size = 0.65) +
  scale_fill_gradientn(colours = c("#F4F9ED", "#E4F1D6", "#B2DF8A")) + scale_x_continuous(expand = expansion(mult = venn_horizontal_padding)) +
  theme_void(base_family = "sans") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17, color = "#3A6831", lineheight = 0.95), legend.position = "none",
                                           plot.margin = venn_plot_margin) + coord_fixed(clip = "off")
venn_rRNA$layers[[3]]$aes_params$fontface <- "bold"
venn_rRNA$layers[[3]]$data$X <- c(-1.5, 5.5, 2)

# Venn diagram: mitochondrial respiration.
venn_mito <- ggVennDiagram(list(North = genes_mito_North, Central = genes_mito_Central, South = genes_mito_South), label_alpha = 0, label_size = 14 / .pt, set_size = 16 / .pt, edge_size = 0.65) +
  scale_fill_gradientn(colours = c("#F4F9FC", "#E3EFF6", "#A6CEE3")) + scale_x_continuous(expand = expansion(mult = venn_horizontal_padding)) +
  theme_void(base_family = "sans") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17, color = "#1F78B4", lineheight = 0.95), legend.position = "none",
                                           plot.margin = venn_plot_margin) + coord_fixed(clip = "off")
venn_mito$layers[[3]]$aes_params$fontface <- "bold"
venn_mito$layers[[3]]$data$X <- c(-1.5, 5.5, 2)

# Venn diagram: Tf2-like retrotransposons.
venn_TF2 <- ggVennDiagram(list(North = genes_TF2_North, Central = genes_TF2_Central, South = genes_TF2_South), label_alpha = 0, label_size = 14 / .pt, set_size = 16 / .pt, edge_size = 0.65) +
  scale_fill_gradientn(colours = c("#FAF7FB", "#EFE8F2", "#CAB2D6")) + scale_x_continuous(expand = expansion(mult = venn_horizontal_padding)) +
  theme_void(base_family = "sans") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17, color = "#6A3D9A", lineheight = 0.95), legend.position = "none",
                                           plot.margin = venn_plot_margin) + coord_fixed(clip = "off")
venn_TF2$layers[[3]]$aes_params$fontface <- "bold"
venn_TF2$layers[[3]]$data$X <- c(-1.5, 5.5, 2)

# Move the population labels of the Venn diagrams (in points); the original labels are kept invisible so the plot limits do not change
venn_south_bottom_space_pt <- 4
venn_label_moves_pt <- data.frame(population = c("North", "Central", "South"),
                                  right = c(-13, 4, 0), up = c(3, 3, -7))
place_venn_labels <- function(p) {
  labels <- p$layers[[3]]$data
  p$layers[[3]]$aes_params$alpha <- 0
  for (i in seq_len(nrow(labels))) {
    move <- venn_label_moves_pt[match(labels$name[i], venn_label_moves_pt$population), ]
    label <- grid::textGrob(labels$name[i],
                            x = grid::unit(1 / 2, "npc") + grid::unit(move$right, "pt"),
                            y = grid::unit(1 / 2, "npc") + grid::unit(move$up, "pt"),
                            gp = grid::gpar(fontfamily = "sans", fontsize = 16, fontface = "bold"))
    p <- p + annotation_custom(label, xmin = labels$X[i], xmax = labels$X[i],
                               ymin = labels$Y[i], ymax = labels$Y[i])
  }
  p
}
venn_heat <- place_venn_labels(venn_heat)
venn_rRNA <- place_venn_labels(venn_rRNA)
venn_mito <- place_venn_labels(venn_mito)
venn_TF2 <- place_venn_labels(venn_TF2)

fig3B <- venn_heat + venn_rRNA + venn_mito + venn_TF2 + plot_layout(ncol = 4)
fig3B

# How strongly the genes behind each cluster are shared among populations: Szymkiewicz-Simpson overlap coefficient,
# extended to three sets = genes shared by all three populations / gene number of the smallest population set.
overlap_result <- data.frame(cluster = c("Heat", "rRNA", "Mito", "TF2"),
                             n_North = c(length(genes_heat_North), length(genes_rRNA_North), length(genes_mito_North), length(genes_TF2_North)),
                             n_Central = c(length(genes_heat_Central), length(genes_rRNA_Central), length(genes_mito_Central), length(genes_TF2_Central)),
                             n_South = c(length(genes_heat_South), length(genes_rRNA_South), length(genes_mito_South), length(genes_TF2_South)),
                             n_shared = c(length(intersect(intersect(genes_heat_North, genes_heat_Central), genes_heat_South)),
                                          length(intersect(intersect(genes_rRNA_North, genes_rRNA_Central), genes_rRNA_South)),
                                          length(intersect(intersect(genes_mito_North, genes_mito_Central), genes_mito_South)),
                                          length(intersect(intersect(genes_TF2_North, genes_TF2_Central), genes_TF2_South))))
overlap_result$n_min <- pmin(overlap_result$n_North, overlap_result$n_Central, overlap_result$n_South)
overlap_result$overlap_coefficient <- ifelse(overlap_result$n_min > 0, overlap_result$n_shared / overlap_result$n_min, NA_real_)
overlap_result

#################### Heat response genes shared by all three populations

# Shared heat response genes with their functional annotation, exported for manual curation of gene families.
genes_heat_shared <- intersect(intersect(genes_heat_North, genes_heat_Central), genes_heat_South)
length(genes_heat_shared)
heat_shared_annotation <- merge(data.frame(gene = genes_heat_shared, stringsAsFactors = FALSE), all.gene.annotation, by = "gene", all.x = TRUE, sort = FALSE)
write.table(heat_shared_annotation, "heat_intersect_annotation_for_curation.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Manually curated annotation of the shared heat response genes, including their gene family (e.g. Hsp70).
heat_curated <- read.delim("heat_intersect_43genes_withannotation_Manual_curation.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
hsp70_genes <- heat_curated$gene[heat_curated$family == "Hsp70"]


#################### Figure 3C-D: KEGG enrichment of general temperature response genes
# The KEGG annotation of genes is built once and saved as gene_ko_description_all.txt; later runs read this file.
build_figure3_kegg <- function(width_in = figure3_width_in) {
  if (!file.exists("gene_ko_description_all.txt")) {
    library(KEGGREST)
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
    stringtie <- read.table("gene_to_transcript_mapping_StringTie.txt", sep = "\t", header = TRUE)
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
    
    
    # Pathways that do not apply to a snail: Human Diseases pathways outside ko05, drug resistance and plant pathways,
    # Organismal Systems pathways describing vertebrate-specific cell types, organs or signaling systems, and the vertebrate
    # versions of pathways for which KEGG also provides an invertebrate version.
    excluded_pathways <- c(
      # Human Diseases outside ko05 (diabetes, liver disease, Cushing syndrome, AGE-RAGE in diabetic complications)
      "ko04930", "ko04940", "ko04950", "ko04936", "ko04932", "ko04931", "ko04933", "ko04934",
      # Drug resistance
      "ko01501", "ko01502", "ko01503", "ko01521", "ko01524", "ko01523", "ko01522",
      # Plant pathways
      "ko04712", "ko04626",
      # Vertebrate immune cells and molecules (blood cells, lymphocytes, antibodies, MHC, chemokines, coagulation)
      "ko04640", "ko04610", "ko04611", "ko04613", "ko04625", "ko04650", "ko04612", "ko04660", "ko04658", "ko04659",
      "ko04662", "ko04664", "ko04666", "ko04670", "ko04672", "ko04062",
      # Vertebrate endocrine glands and hormones
      "ko04911", "ko04922", "ko04923", "ko04920", "ko03320", "ko04929", "ko04912", "ko04913", "ko04915", "ko04914",
      "ko04917", "ko04921", "ko04926", "ko04935", "ko04918", "ko04919", "ko04928", "ko04916", "ko04924", "ko04614",
      "ko04925", "ko04927",
      # Vertebrate circulatory, digestive and excretory organs
      "ko04261", "ko04270", "ko04970", "ko04971", "ko04972", "ko04976", "ko04973", "ko04974", "ko04975", "ko04979",
      "ko04977", "ko04981", "ko04980", "ko04978", "ko04962", "ko04960", "ko04961", "ko04964", "ko04966",
      # Vertebrate nervous, sensory and skeletal systems, and mammalian thermogenesis
      "ko04720", "ko04730", "ko04723", "ko04722", "ko04740", "ko04742", "ko04750", "ko04380", "ko04714",
      # Vertebrate versions of pathways with an invertebrate version in KEGG (Toll-like receptor, phototransduction,
      # circadian rhythm and entrainment, longevity)
      "ko04620", "ko04744", "ko04710", "ko04713", 
      # Longevity regulating pathways, which combine parts of other pathways with downstream stress-response genes
      "ko04211", "ko04212", "ko04213")
    
    
    k_to_pathway <- k_to_pathway[!grepl("^ko05|^ko07", k_to_pathway$pathway) & !k_to_pathway$pathway %in% excluded_pathways, ]
    
    # Save the gene-K-pathway table used for KEGG enrichment (genes whose K numbers have no retained pathway keep NA as pathway).
    gene_ko <- merge(gene_K, k_to_pathway, by = "K", all.x = TRUE)[, c("gene", "K", "K_description", "pathway", "pathway_description")]
    write.table(gene_ko, file = "gene_ko_description_all.txt", sep = "\t", quote = F, row.names = F, col.names = T)
    
    
    
  }
  # Gene-pathway annotation (genes with at least one retained KEGG pathway).
  gene_ko <- read.delim("gene_ko_description_all.txt", header = TRUE, stringsAsFactors = FALSE, quote = "")
  gene_pathway <- dplyr::distinct(gene_ko[!is.na(gene_ko$pathway), c("gene", "pathway", "pathway_description")])
  TERM2GENE <- dplyr::distinct(gene_pathway[, c("pathway", "gene")])
  TERM2NAME <- dplyr::distinct(gene_pathway[, c("pathway", "pathway_description")])
  
  # General temperature response of all retained genes in each population (from 04DESeq2_Figure2.R).
  deg_North <- read.delim("../05DEseq2/North_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
  deg_Central <- read.delim("../05DEseq2/Central_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
  deg_South <- read.delim("../05DEseq2/South_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
  
  # Background for enrichment: all retained genes (final filtered count matrix) that belong to at least one KEGG pathway.
  counts <- read.table("../01counts_ind_filter/counts_filter.txt", header = TRUE)
  universe <- intersect(rownames(counts), TERM2GENE$gene)
  length(universe)
  
  
  #################### KEGG enrichment of up- and down-regulated general temperature response genes
  # A pathway is significantly enriched when p.adjust <= 0.01 and qvalue <= 0.01.
  
  #################### North
  up_North <- deg_North$gene[which(deg_North$padj < 0.05 & deg_North$log2FoldChange >= 0.5)]
  kegg_up_North <- enricher(gene = up_North, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
                            minGSSize = 10, maxGSSize = 3000, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
  res_up_North <- kegg_up_North@result
  sig_up_North <- res_up_North$ID[res_up_North$p.adjust <= 0.01 & res_up_North$qvalue <= 0.01]
  
  down_North <- deg_North$gene[which(deg_North$padj < 0.05 & deg_North$log2FoldChange <= -0.5)]
  kegg_down_North <- enricher(gene = down_North, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
                              minGSSize = 10, maxGSSize = 3000, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
  res_down_North <- kegg_down_North@result
  sig_down_North <- res_down_North$ID[res_down_North$p.adjust <= 0.01 & res_down_North$qvalue <= 0.01]
  
  #################### Central
  up_Central <- deg_Central$gene[which(deg_Central$padj < 0.05 & deg_Central$log2FoldChange >= 0.5)]
  kegg_up_Central <- enricher(gene = up_Central, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
                              minGSSize = 10, maxGSSize = 3000, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
  res_up_Central <- kegg_up_Central@result
  sig_up_Central <- res_up_Central$ID[res_up_Central$p.adjust <= 0.01 & res_up_Central$qvalue <= 0.01]
  
  down_Central <- deg_Central$gene[which(deg_Central$padj < 0.05 & deg_Central$log2FoldChange <= -0.5)]
  kegg_down_Central <- enricher(gene = down_Central, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
                                minGSSize = 10, maxGSSize = 3000, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
  res_down_Central <- kegg_down_Central@result
  sig_down_Central <- res_down_Central$ID[res_down_Central$p.adjust <= 0.01 & res_down_Central$qvalue <= 0.01]
  
  #################### South
  up_South <- deg_South$gene[which(deg_South$padj < 0.05 & deg_South$log2FoldChange >= 0.5)]
  kegg_up_South <- enricher(gene = up_South, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
                            minGSSize = 10, maxGSSize = 3000, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
  res_up_South <- kegg_up_South@result
  sig_up_South <- res_up_South$ID[res_up_South$p.adjust <= 0.01 & res_up_South$qvalue <= 0.01]
  
  down_South <- deg_South$gene[which(deg_South$padj < 0.05 & deg_South$log2FoldChange <= -0.5)]
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
  # KEGG data for Figure 3C-D
  
  # The four pathways significantly enriched among up-regulated genes, in the order of the figure.
  pathways_up <- c("ko04141", "ko04010", "ko03010", "ko00190")
  pathway_labels <- c("Protein processing in\nendoplasmic reticulum", "MAPK signaling pathway", "Ribosome", "Oxidative phosphorylation")
  
  #################### Data of Figure 3C: gene count and FDR of each pathway in each population
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
  
  #################### GO cluster of the genes in each pathway (clusters of Figure 3A)
  gene_GO <- read.delim("gene_go_description_all.txt", header = TRUE, stringsAsFactors = FALSE)
  GO_cluster_terms <- read.delim("GO_cluster_terms.txt", header = TRUE, stringsAsFactors = FALSE)
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
  
  #################### Data of Figure 3D: up-regulated genes of each pathway shared among populations (UpSet)
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
  
  
  
  # Panels C and D are drawn in mm on a 170 mm wide layout and scaled to the Figure 3 width; text sizes are not scaled
  drawing_scale <- width_in * 25.4 / 170
  drawing_font_scale <- 14 * 25.4 / 170
  unit <- function(x, units, ...) grid::unit(if (identical(units, "mm")) x * drawing_scale else x, units, ...)
  gpar <- function(...) {
    a <- list(...)
    if (!is.null(a$fontsize)) a$fontsize <- a$fontsize * drawing_font_scale
    if (!is.null(a$lwd)) a$lwd <- a$lwd * drawing_scale
    do.call(grid::gpar, a)
  }
  # Vertical range (mm) of panels C and D in the layout
  panel_windows <- list(C = c(bottom = 70, top = 102), D = c(bottom = 8, top = 65))
  draw_kegg <- function(panel) {
    grid.newpage()
    window <- panel_windows[[panel]]
    pushViewport(viewport(y = unit(-window["bottom"], "mm"), just = "bottom",
                          width = unit(170, "mm"), height = unit(104, "mm")))
    
    if (panel == "C") {
      #### C: bubbles of the four pathways (rows) in the three populations (columns)
      pop_x <- c(19.33, 40, 60.67)                    # centres of the North / Central / South columns (from x = 9 to 71 mm)
      pathway_y <- c(90, 84, 78, 72)                  # rows of the four pathways
      fdr_palette <- colorRampPalette(c("#AADBB4", "#16813D"))(256)
      fdr_limit <- max(10, ceiling(max(fig4A$log10FDR[fig4A$significant]) / 10) * 10)
      fig4A$x <- pop_x[match(fig4A$population, c("North", "Central", "South"))]
      fig4A$y <- pathway_y[fig4A$pathway]
      fig4A$radius <- 1.8 * sqrt(fig4A$Count / 80)    # bubble area proportional to the gene count
      fig4A$fill <- ifelse(fig4A$significant, fdr_palette[1 + round(255 * pmin(fig4A$log10FDR / fdr_limit, 1))], "#CCCCCC")
      
      grid.text("C", x = unit(2.5, "mm"), y = unit(100, "mm"), gp = gpar(fontsize = 9, fontface = "bold"))
      # Population header
      header_edges <- seq(9, 71, length.out = 4)
      header_centres <- (header_edges[-1] + header_edges[-4]) / 2
      grid.rect(x = unit(40, "mm"), y = unit(96, "mm"), width = unit(62, "mm"), height = unit(4.6, "mm"),
                gp = gpar(fill = "#E6E6E6", col = "black", lwd = 0.35))
      grid.segments(x0 = unit(header_edges[2:3], "mm"), x1 = unit(header_edges[2:3], "mm"),
                    y0 = unit(93.7, "mm"), y1 = unit(98.3, "mm"), gp = gpar(col = "black", lwd = 0.35))
      grid.text(c("North", "Central", "South"), x = unit(header_centres, "mm"), y = unit(96, "mm"),
                gp = gpar(fontsize = 17 / drawing_font_scale, fontface = "bold"))
      
      grid.segments(x0 = unit(9, "mm"), x1 = unit(71, "mm"), y0 = unit(pathway_y, "mm"), y1 = unit(pathway_y, "mm"), gp = gpar(col = "#DEDEDE", lwd = 0.5))
      grid.text(sub("\n", " ", pathway_labels), x = unit(73, "mm"), y = unit(pathway_y, "mm"), just = "left", gp = gpar(fontsize = 16 / drawing_font_scale, fontface = "bold"))
      grid.circle(x = unit(fig4A$x, "mm"), y = unit(fig4A$y, "mm"), r = unit(fig4A$radius, "mm"), gp = gpar(fill = fig4A$fill, col = NA))
      
      # Legend of C: gene count (circles) and -log10(FDR) (colour bar)
      legend_counts <- c(20, 40, 60, 80)
      legend_move_pt <- c(right = 9, down = 5)
      pushViewport(viewport(x = grid::unit(1 / 2, "npc") + grid::unit(legend_move_pt["right"], "pt"),
                            y = grid::unit(1 / 2, "npc") - grid::unit(3, "mm") - grid::unit(legend_move_pt["down"], "pt")))
      legend_y <- c(91.5, 87, 82.5, 78)
      # Legend frame
      legend_frame_left <- unit(134, "mm") - grid::unit(9, "pt")
      legend_frame_right <- grid::unit(1, "npc") - grid::unit(16 + legend_move_pt["right"], "pt")
      grid.rect(x = legend_frame_left, y = unit(87, "mm"), width = legend_frame_right - legend_frame_left,
                height = unit(25, "mm"), just = c("left", "centre"),
                gp = gpar(fill = NA, col = "#999999", lwd = 0.4))
      grid.text("Gene count", x = unit(134, "mm"), y = unit(97.5, "mm"), just = c("left", "top"), gp = gpar(fontsize = 7, fontface = "bold"))
      grid.circle(x = unit(135.8, "mm"), y = unit(legend_y, "mm"), r = unit(1.8 * sqrt(legend_counts / 80), "mm"), gp = gpar(fill = NA, col = "black", lwd = 0.4))
      grid.text(legend_counts, x = unit(138.8, "mm"), y = unit(legend_y, "mm"), just = "left", gp = gpar(fontsize = 7))
      grid.text(expression(bold("-"*log[10]*"(FDR)")), x = unit(151, "mm"), y = unit(97.5, "mm"), just = c("left", "top"), gp = gpar(fontsize = 7))
      bar_bottom <- 78
      bar_top <- 91.5
      grid.rect(x = unit(153, "mm"), y = unit(bar_bottom + (1:256 - 0.5) * (bar_top - bar_bottom) / 256, "mm"), width = unit(1.8, "mm"),
                height = unit((bar_top - bar_bottom) / 256 + 0.02, "mm"), gp = gpar(fill = fdr_palette, col = NA))
      tick_y <- bar_bottom + c(0, 0.5, 1) * (bar_top - bar_bottom)
      grid.segments(x0 = unit(153.9, "mm"), x1 = unit(154.4, "mm"), y0 = unit(tick_y, "mm"), y1 = unit(tick_y, "mm"), gp = gpar(lwd = 0.3))
      grid.text(c(0, fdr_limit / 2, fdr_limit), x = unit(154.9, "mm"), y = unit(tick_y, "mm"), just = "left", gp = gpar(fontsize = 7))
      popViewport()
      
    }  # panel C
    
    if (panel == "D") {
      #### D: one UpSet panel per pathway (bars = genes in each combination of populations; dots below = the combination)
      upset_ymax <- max(85, ceiling(max(fig4B$Count) / 10) * 10)   # common bar scale for the four panels
      panel_width <- 34.5
      panel_left <- 18 + (0:3) * (panel_width + 3)
      population_spacing_D <- grid::convertHeight(grid::unit(population_row_spacing_D_pt, "pt"), "mm", valueOnly = TRUE) / drawing_scale
      set_y <- 23.3 + c(1, 0, -1) * population_spacing_D  # North / Central / South; shared by labels, dots and backgrounds
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
      
      grid.text("D", x = unit(2.5, "mm"), y = unit(63, "mm"), gp = gpar(fontsize = 9, fontface = "bold"))
      # Align the population labels of D with those of A
      population_label_origin_D <- unit(6.5, "mm")
      pushViewport(viewport(x = grid::unit(1 / 2, "npc") + grid::unit(population_label_left_pt, "pt") - population_label_origin_D))
      
      grid.text(c("North", "Central", "South"), x = population_label_origin_D, y = unit(set_y, "mm"), just = "left", gp = gpar(fontsize = 17 / drawing_font_scale, fontface = "bold"))
      grid.rect(x = unit(rep(panel_left + panel_width / 2, 2), "mm"), y = unit(rep(set_y[c(1, 3)], each = 4), "mm"), width = unit(panel_width, "mm"),
                height = grid::unit(population_band_height_pt, "pt"), gp = gpar(fill = "#EEEEEE", col = NA))
      
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
      # Pathway names below the panels
      grid.text(pathway_labels[1], x = unit(panel_left[1], "mm"), y = unit(14.6, "mm"), just = c("left", "top"),
                gp = gpar(fontsize = 16 / drawing_font_scale, fontface = "bold", lineheight = 0.98))
      grid.text(pathway_labels[-1], x = unit(panel_left[-1] + panel_width / 2, "mm"), y = unit(14.6, "mm"), just = "top",
                gp = gpar(fontsize = 16 / drawing_font_scale, fontface = "bold", lineheight = 0.98))
      
      # Mark the Hsp70 genes among the genes shared by all three populations in Protein processing in ER and MAPK signaling
      hsp70 <- fig4B[fig4B$panel %in% 1:2 & fig4B$combination == "111", ]
      hsp70$target_x <- hsp70$x + ifelse(hsp70$panel == 1, 1, -1) * (hsp70$spacing * 0.42 + 0.8)
      hsp70$target_y <- bar_baseline + hsp70$h_main + hsp70$h_other - pmin(0.9, (hsp70$h_main + hsp70$h_other) / 2)
      root_x <- mean(hsp70$target_x)
      root_y <- max(hsp70$target_y) + 8
      grid.text("Hsp70", x = unit(root_x, "mm"), y = unit(root_y + 2.8, "mm"), gp = gpar(fontsize = 8, fontface = "bold", col = "#C46A00"))
      grid.segments(x0 = unit(root_x, "mm"), y0 = unit(root_y, "mm"), x1 = unit(hsp70$target_x, "mm"), y1 = unit(hsp70$target_y, "mm"),
                    arrow = arrow(length = unit(1.4, "mm"), type = "closed"), gp = gpar(col = "#C46A00", fill = "#C46A00", lty = 2, lwd = 0.65))
      
      popViewport()  # complete D drawing
    }  # panel D
    popViewport()  # original layout coordinate system
  }
  panel_heights_in <- vapply(panel_windows, function(window) width_in * diff(window) / 170, numeric(1))
  panels <- list(C = grid::grid.grabExpr(draw_kegg("C"), width = width_in, height = panel_heights_in["C"], wrap = TRUE),
                 D = grid::grid.grabExpr(draw_kegg("D"), width = width_in, height = panel_heights_in["D"], wrap = TRUE),
                 heights_in = panel_heights_in)
  attr(panels, "panel_C_data") <- fig4A
  attr(panels, "panel_D_data") <- fig4B
  panels
}

figure3_kegg_panels <- build_figure3_kegg(width_in = figure3_width_in)

#################### Combine Figure 3A-D
panel_gap_pt <- 20
figure3_outer_margin_pt <- 6
points_to_inches <- function(x) grid::convertHeight(grid::unit(x, "pt"), "inches", valueOnly = TRUE)

# Add the panel label and trim the blank margin below panels A and B
frame_go_panel <- function(plot, label, height_in, trim_bottom_pt, clip = "on") {
  labelled <- ggarrange(plot, labels = label, font.label = list(size = 20, face = "bold"))
  original <- ggplotGrob(labelled)
  inside <- grid::grobTree(original,
                           vp = grid::viewport(y = -grid::unit(trim_bottom_pt, "pt"), just = "bottom",
                                               width = grid::unit(1, "npc"), height = grid::unit(height_in, "inches")))
  list(grob = grid::grobTree(inside, vp = grid::viewport(clip = clip)),
       height_in = height_in - points_to_inches(4 + trim_bottom_pt))
}
panel_A <- frame_go_panel(fig3A, "A", figure3_height_A, trim_bottom_pt = 2)
panel_B <- frame_go_panel(fig3B, "B", 4, trim_bottom_pt = 8, clip = "off")
figure3_heights_in <- c(points_to_inches(figure3_outer_margin_pt),
                        panel_A$height_in, points_to_inches(panel_gap_pt),
                        panel_B$height_in, points_to_inches(panel_gap_pt + venn_south_bottom_space_pt),
                        figure3_kegg_panels$heights_in["C"], points_to_inches(panel_gap_pt),
                        figure3_kegg_panels$heights_in["D"], points_to_inches(figure3_outer_margin_pt))
fig3 <- ggarrange(NULL, panel_A$grob, NULL, panel_B$grob, NULL,
                  figure3_kegg_panels$C, NULL, figure3_kegg_panels$D, NULL,
                  ncol = 1, nrow = 9, heights = figure3_heights_in)
figure3_height_in <- sum(figure3_heights_in)
ggsave("Figure3.pdf", fig3, width = figure3_width_in, height = figure3_height_in, device = cairo_pdf, bg = "white")

############################################################
#Part 4: up-regulated GO terms in each population, ordered by the GO tree (Figure S7)

#################### Figure S7

# Tree coordinates of the same GO tree and clusters as Figure 3A
nodes_S7 <- as.data.frame(tree_ordered$data)
nodes_S7$x <- nodes_S7$x / max(nodes_S7$x) * 2.6
tips_S7 <- nodes_S7[nodes_S7$isTip, ]
tips_S7 <- tips_S7[order(tips_S7$y), ]
tips_S7$GO <- go_up_all@result$ID[match(tips_S7$label, go_up_all@result$Description)]
tips_S7$Description <- go_up_all@result$Description[match(tips_S7$GO, go_up_all@result$ID)]
N_GO_S7 <- nrow(tips_S7)

# Width of the longest GO description, so that the bubbles start right after the descriptions
label_size_S7 <- 12
metric_file_S7 <- tempfile(fileext = ".pdf")
grDevices::cairo_pdf(metric_file_S7, family = "sans")
grid::pushViewport(grid::viewport(gp = grid::gpar(fontfamily = "sans", fontsize = label_size_S7)))
label_width_in_S7 <- max(grid::convertWidth(grid::stringWidth(tips_S7$Description), "inches", valueOnly = TRUE)) + 0.12
grid::popViewport()
grDevices::dev.off()
unlink(metric_file_S7)
x_unit_in_S7 <- (10 - grid::convertWidth(grid::unit(22, "pt"), "inches", valueOnly = TRUE)) / 13.7
label_xmin_S7 <- 2.7
label_xmax_S7 <- label_xmin_S7 + label_width_in_S7 / x_unit_in_S7
x_shift_S7 <- label_xmax_S7 - 3.95
figure_width_S7 <- 10 + (x_shift_S7 - 4.5) * x_unit_in_S7

# Tree branches (horizontal and vertical segments), coloured by cluster.
edges_S7 <- nodes_S7[nodes_S7$node != nodes_S7$parent, ]
edges_S7$x_parent <- nodes_S7$x[match(edges_S7$parent, nodes_S7$node)]
edges_S7$y_parent <- nodes_S7$y[match(edges_S7$parent, nodes_S7$node)]
edges_S7$branch_color <- unname(branch_colors[as.character(edges_S7$group)])

# Bubbles: number of up-regulated genes in each GO term (numerator of GeneRatio) and FDR (p.adjust), per population.
go_data_S7 <- go_up_compare[, c("ID", "Description")]
go_data_S7$Count_North <- as.numeric(sub("/.*", "", go_up_compare$GeneRatio_North))
go_data_S7$FDR_North <- go_up_compare$padj_North
go_data_S7$Count_Central <- as.numeric(sub("/.*", "", go_up_compare$GeneRatio_Central))
go_data_S7$FDR_Central <- go_up_compare$padj_Central
go_data_S7$Count_South <- as.numeric(sub("/.*", "", go_up_compare$GeneRatio_South))
go_data_S7$FDR_South <- go_up_compare$padj_South
bubble_S7 <- dplyr::bind_rows(dplyr::transmute(go_data_S7, ID, population = "North", Count = Count_North, FDR = FDR_North),
                              dplyr::transmute(go_data_S7, ID, population = "Central", Count = Count_Central, FDR = FDR_Central),
                              dplyr::transmute(go_data_S7, ID, population = "South", Count = Count_South, FDR = FDR_South))
bubble_S7$y <- tips_S7$y[match(bubble_S7$ID, tips_S7$GO)]
bubble_S7$x <- unname(c(North = 4.425, Central = 5.375, South = 6.325)[bubble_S7$population]) + x_shift_S7

# Bubble colour shows -log10 FDR only for significant terms (FDR < 0.01); non-significant terms are grey; untested terms are not drawn.
bubble_S7$ColorCol <- ifelse(!is.na(bubble_S7$FDR) & bubble_S7$FDR < 0.01, -log10(pmax(bubble_S7$FDR, .Machine$double.xmin)), NA_real_)
bubble_S7 <- bubble_S7[!is.na(bubble_S7$Count), ]

# Coloured background behind the GO descriptions of each cluster (the tree itself stays on white).
label_bands_S7 <- dplyr::summarise(dplyr::group_by(dplyr::filter(tips_S7, as.character(group) != "0"), group), ymin = min(y) - 0.5, ymax = max(y) + 0.5, .groups = "drop")
label_bands_S7$fill <- unname(branch_colors[as.character(label_bands_S7$group)])

# Column headers for the three populations.
headers_S7 <- data.frame(population = c("North", "Central", "South"), xmin = c(3.95, 4.90, 5.85), xmax = c(4.90, 5.85, 6.80), x = c(4.425, 5.375, 6.325))
headers_S7[c("xmin", "xmax", "x")] <- headers_S7[c("xmin", "xmax", "x")] + x_shift_S7

# Cluster names on the right.
annotations_S7 <- label_bands_S7
annotations_S7$label <- cluster_definitions$title[match(as.character(annotations_S7$group), cluster_definitions$group)]
annotations_S7$label[as.character(annotations_S7$group) == "3"] <- "Tf2-like retrotransposons"
annotations_S7$text_color <- cluster_definitions$text_color[match(as.character(annotations_S7$group), cluster_definitions$group)]
annotations_S7$mid <- (annotations_S7$ymin + annotations_S7$ymax) / 2
annotations_S7$x_text <- 7.67 + x_shift_S7

# A further group of GO terms outside the four clusters, labelled only (not a separate cluster).
viral_y_S7 <- range(tips_S7$y[tips_S7$GO %in% c("GO:0044185", "GO:0042575")])

# Tree, GO descriptions, bubbles and labels share one y axis, so every row lines up.
figS7 <- ggplot() +
  geom_segment(data = edges_S7, aes(x = x_parent, xend = x, y = y, yend = y), color = edges_S7$branch_color, linewidth = 0.85, lineend = "butt") +
  geom_segment(data = edges_S7, aes(x = x_parent, xend = x_parent, y = y_parent, yend = y), color = edges_S7$branch_color, linewidth = 0.85, lineend = "butt") +
  geom_rect(data = label_bands_S7, aes(ymin = ymin, ymax = ymax), xmin = label_xmin_S7, xmax = label_xmax_S7, fill = label_bands_S7$fill, alpha = 0.45, color = NA) +
  geom_segment(data = tips_S7, aes(y = y, yend = y), x = 3.95 + x_shift_S7, xend = 6.80 + x_shift_S7, color = "grey90", linewidth = 0.25) +
  geom_text(data = tips_S7, aes(y = y, label = Description), x = label_xmin_S7, hjust = 0, size = label_size_S7, size.unit = "pt", family = "sans", color = "black") +
  geom_rect(data = headers_S7, aes(xmin = xmin, xmax = xmax), ymin = N_GO_S7 + 0.65, ymax = N_GO_S7 + 2.65, fill = "grey90", color = "grey35", linewidth = 0.3) +
  geom_text(data = headers_S7, aes(x = x, label = population), y = N_GO_S7 + 1.65, size = 14, size.unit = "pt") +
  geom_point(data = bubble_S7, aes(x = x, y = y, size = Count, color = ColorCol)) +
  scale_size_continuous(range = c(1, 3.6), breaks = c(25, 50, 75, 100, 125), name = "Gene count",
                        guide = guide_legend(order = 2, nrow = 1, title.position = "left", override.aes = list(shape = 21, fill = NA, color = "black", stroke = 0.4))) +
  scale_color_gradient(low = "#97D1A0", high = "#1E803D", na.value = "grey80", name = expression(plain("-log")[10]~FDR), breaks = c(50, 100),
                       guide = guide_colorbar(order = 1, direction = "horizontal", title.position = "left", barwidth = grid::unit(1.8, "cm"), barheight = grid::unit(0.35, "cm"))) +
  geom_segment(data = annotations_S7, aes(y = ymin + 0.2, yend = ymax - 0.2), x = 7.38 + x_shift_S7, xend = 7.38 + x_shift_S7, color = annotations_S7$fill, linewidth = 0.6) +
  geom_text(data = annotations_S7, aes(x = x_text, y = mid, label = label), angle = 90, size = 14, size.unit = "pt", fontface = "bold", lineheight = 0.95, color = annotations_S7$text_color) +
  annotate("segment", x = 7.96 + x_shift_S7, xend = 7.96 + x_shift_S7, y = viral_y_S7[1], yend = viral_y_S7[2], color = "black", linewidth = 0.5) +
  annotate("text", x = 8.30 + x_shift_S7, y = mean(viral_y_S7), label = "Viral processes and DNA metabolism", angle = 90, size = 14, size.unit = "pt", fontface = "bold", color = "black") +
  scale_x_continuous(limits = c(-0.1, 9.1 + x_shift_S7), expand = c(0, 0)) + scale_y_continuous(limits = c(-5.7, N_GO_S7 + 4), expand = c(0, 0)) +
  theme_void(base_family = "sans") +
  theme(legend.position = "inside", legend.position.inside = c(0.64, 0.015), legend.justification = c(0.5, 0), legend.box = "horizontal",
        legend.title = element_text(size = 12), legend.text = element_text(size = 11), legend.key.width = grid::unit(0.4, "cm"),
        legend.spacing.x = grid::unit(1, "cm"), plot.margin = margin(13, 12, 5, 10), plot.background = element_rect(fill = "white", color = NA))
figS7
ggsave("FigureS7.pdf", figS7, width = figure_width_S7, height = 15.5, device = cairo_pdf, bg = "white")


#################### Number of retained genes with GO terms, KO identifiers and KEGG pathways
counts <- read.table("../01counts_ind_filter/counts_filter.txt", header = TRUE)
gene_GO <- read.delim("gene_go_description_all.txt", header = TRUE, stringsAsFactors = FALSE)
gene_ko <- read.delim("gene_ko_description_all.txt", header = TRUE, stringsAsFactors = FALSE, quote = "")

n_genes <- nrow(counts)
n_GO <- sum(rownames(counts) %in% gene_GO$gene)
n_KO <- sum(rownames(counts) %in% gene_ko$gene)
n_pathway <- sum(rownames(counts) %in% gene_ko$gene[!is.na(gene_ko$pathway)])
c(genes = n_genes, GO = n_GO, GO_percent = round(100 * n_GO / n_genes, 1),
  KO = n_KO, KO_percent = round(100 * n_KO / n_genes, 1), pathway = n_pathway, pathway_percent = round(100 * n_pathway / n_genes, 1))


#################### GO enrichment of shore-level-specific response genes
# All genes significant in the shore-level-specific contrast (orange and purple points in Figure 2A-C), with the same
# background as for the general temperature response genes
shore_GO_universe <- intersect(rownames(counts), gene_GO$gene)
shore_North <- read.delim("../05DEseq2/North_shorelevel_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
shore_Central <- read.delim("../05DEseq2/Central_shorelevel_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
shore_South <- read.delim("../05DEseq2/South_shorelevel_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
shore_genes_North <- shore_North$gene[which(shore_North$padj < 0.05 & abs(shore_North$log2FoldChange) >= 0.5)]
shore_genes_Central <- shore_Central$gene[which(shore_Central$padj < 0.05 & abs(shore_Central$log2FoldChange) >= 0.5)]
shore_genes_South <- shore_South$gene[which(shore_South$padj < 0.05 & abs(shore_South$log2FoldChange) >= 0.5)]

go_shore_North <- enricher(gene = shore_genes_North, universe = shore_GO_universe, TERM2GENE = gene_GO[c("go", "gene")], TERM2NAME = gene_GO[c("go", "description")], minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
go_shore_Central <- enricher(gene = shore_genes_Central, universe = shore_GO_universe, TERM2GENE = gene_GO[c("go", "gene")], TERM2NAME = gene_GO[c("go", "description")], minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
go_shore_South <- enricher(gene = shore_genes_South, universe = shore_GO_universe, TERM2GENE = gene_GO[c("go", "gene")], TERM2NAME = gene_GO[c("go", "description")], minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)

# Significantly enriched GO terms (p.adjust <= 0.01 and qvalue <= 0.01)
significant_shore_go <- function(enrichment) {
  columns <- c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "qvalue", "geneID")
  if (is.null(enrichment)) return(as.data.frame(setNames(rep(list(character()), length(columns)), columns)))
  result <- enrichment@result
  selected <- which(!is.na(result$p.adjust) & !is.na(result$qvalue) & result$p.adjust <= 0.01 & result$qvalue <= 0.01)
  result[selected, columns, drop = FALSE]
}
sig_go_shore_North <- significant_shore_go(go_shore_North)
sig_go_shore_Central <- significant_shore_go(go_shore_Central)
sig_go_shore_South <- significant_shore_go(go_shore_South)
print(sig_go_shore_North)
print(sig_go_shore_Central)
print(sig_go_shore_South)
shore_go_summary <- data.frame(Population = c("North", "Central", "South"),
                               Shore_genes = c(length(shore_genes_North), length(shore_genes_Central), length(shore_genes_South)),
                               GO_annotated_genes = c(sum(shore_genes_North %in% shore_GO_universe), sum(shore_genes_Central %in% shore_GO_universe), sum(shore_genes_South %in% shore_GO_universe)),
                               Significant_GO_terms = c(nrow(sig_go_shore_North), nrow(sig_go_shore_Central), nrow(sig_go_shore_South)))
print(shore_go_summary)
