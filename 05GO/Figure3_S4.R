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
# Part 3 clusters the up-regulated GO terms by shared genes into four functional groups, compares these groups and their
# genes among populations, and tests front-loading of heat response genes shared by all three populations (Figure 3).
#
# Part 4 shows the gene count and FDR of every up-regulated GO term in each population, ordered by the GO tree
# (Figure S4), and lists the same values as a table (Table S5).
#
# Input:  SwissProt_GO.xlsx, SwissProt_annotation.tab, NR_GO.xlsx, NR_annotation.tab, InterProScan_GO.tsv,
#         InterProScan_annotation.tsv, FANTASIA_results.csv, Lsax_stringtie_merged.gtf, go.db.txt,
#         North/Central/South_generaltemp_allgenes.txt (from 04DESeq2_Figure2.R), counts_filter.txt
#         group.txt, heat_intersect_43genes_withannotation_Manual_curation.txt (manual curation of gene families)
# Output: gene_to_transcript_mapping_StringTie.txt  gene-transcript relationships
#         gene_go_description_all.txt              GO terms of each gene
#         gene_annotation.txt                      functional description of each gene
#         heat_intersect_annotation_for_curation.txt  shared heat response genes, for manual curation
#         GO_up_North_full.txt, GO_up_Central_full.txt, GO_up_South_full.txt  full GO enrichment results of up-regulated genes
#         GO_up76go.txt                            up-regulated GO terms compared across populations
#         GO_cluster_terms.txt                     functional cluster of each up-regulated GO term
#         Figure3.pdf, FigureS4.pdf, Table_S5.csv
# =============================================================================

setwd("/data/05GO/")
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

# Gene-GO annotation from Part 1.
gene_GO <- read.delim("gene_go_description_all.txt", header = TRUE, stringsAsFactors = FALSE)
gene_GO$description <- gsub("^ ", "", gene_GO$description)

# General temperature response of all retained genes in each population (from 04DESeq2_Figure2.R).
deg_North <- read.delim("../04DEseq2/North_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
deg_Central <- read.delim("../04DEseq2/Central_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
deg_South <- read.delim("../04DEseq2/South_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)

# Background for enrichment: all retained genes (final filtered count matrix) that have at least one GO term.
counts <- read.table("../01counts_ind_filter/counts_filter.txt", header = TRUE)
universe <- intersect(rownames(counts), gene_GO$gene)


#################### GO enrichment of up- and down-regulated general temperature response genes
# A GO term is significantly enriched when p.adjust <= 0.01 and qvalue <= 0.01.

#################### North

# Up-regulated genes, and how many of them have GO annotation.
up_North <- deg_North$gene[which(deg_North$padj < 0.05 & deg_North$log2FoldChange >= 1)]
length(up_North)
sum(up_North %in% gene_GO$gene)
go_up_North <- enricher(gene = up_North, universe = universe, TERM2GENE = gene_GO[c("go", "gene")], TERM2NAME = gene_GO[c("go", "description")],
                        minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_up_North <- go_up_North@result
sig_up_North <- res_up_North$ID[res_up_North$p.adjust <= 0.01 & res_up_North$qvalue <= 0.01]
length(sig_up_North)
write.table(res_up_North, "GO_up_North_full.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Down-regulated genes.
down_North <- deg_North$gene[which(deg_North$padj < 0.05 & deg_North$log2FoldChange <= -1)]
length(down_North)
sum(down_North %in% gene_GO$gene)
go_down_North <- enricher(gene = down_North, universe = universe, TERM2GENE = gene_GO[c("go", "gene")], TERM2NAME = gene_GO[c("go", "description")],
                          minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_down_North <- go_down_North@result
sum(res_down_North$p.adjust <= 0.01 & res_down_North$qvalue <= 0.01)



#################### Central

# Up-regulated genes, and how many of them have GO annotation.
up_Central <- deg_Central$gene[which(deg_Central$padj < 0.05 & deg_Central$log2FoldChange >= 1)]
length(up_Central)
sum(up_Central %in% gene_GO$gene)
go_up_Central <- enricher(gene = up_Central, universe = universe, TERM2GENE = gene_GO[c("go", "gene")], TERM2NAME = gene_GO[c("go", "description")],
                          minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_up_Central <- go_up_Central@result
sig_up_Central <- res_up_Central$ID[res_up_Central$p.adjust <= 0.01 & res_up_Central$qvalue <= 0.01]
length(sig_up_Central)
write.table(res_up_Central, "GO_up_Central_full.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Down-regulated genes.
down_Central <- deg_Central$gene[which(deg_Central$padj < 0.05 & deg_Central$log2FoldChange <= -1)]
length(down_Central)
sum(down_Central %in% gene_GO$gene)
go_down_Central <- enricher(gene = down_Central, universe = universe, TERM2GENE = gene_GO[c("go", "gene")], TERM2NAME = gene_GO[c("go", "description")],
                            minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_down_Central <- go_down_Central@result
sum(res_down_Central$p.adjust <= 0.01 & res_down_Central$qvalue <= 0.01)

#################### South

# Up-regulated genes, and how many of them have GO annotation.
up_South <- deg_South$gene[which(deg_South$padj < 0.05 & deg_South$log2FoldChange >= 1)]
length(up_South)
sum(up_South %in% gene_GO$gene)
go_up_South <- enricher(gene = up_South, universe = universe, TERM2GENE = gene_GO[c("go", "gene")], TERM2NAME = gene_GO[c("go", "description")],
                        minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
res_up_South <- go_up_South@result
sig_up_South <- res_up_South$ID[res_up_South$p.adjust <= 0.01 & res_up_South$qvalue <= 0.01]
length(sig_up_South)
write.table(res_up_South, "GO_up_South_full.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Down-regulated genes.
down_South <- deg_South$gene[which(deg_South$padj < 0.05 & deg_South$log2FoldChange <= -1)]
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
write.table(go_up_compare, file = "GO_up76go.txt", sep = "\t", quote = F, row.names = F, col.names = T)


############################################################
#Part 3: functional clusters of up-regulated GO terms and front-loading of heat response genes (Figure 3)

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
                                  text_color = c("#3A6831", "#918009", "#CF92ED", "#82B6D1"),
                                  title = c("Ribosome-related processes", "Heat response genes", "Tf2-like", "Mitochondrial\nrespiration"), stringsAsFactors = FALSE)
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
                               cluster = c("0" = "Unclustered", "1" = "Ribosome-related processes", "2" = "Heat response genes", "3" = "Tf2-like",
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

# Horizontal position of each GO term in Figure 3A: consecutive terms of the same cluster form one run; a small cluster gets
# at least 10 units of width so its count panel and title fit (topology and GO membership do not change).
runs <- rle(as.character(tree_tips$group))
run_width <- ifelse(runs$values == "0", runs$lengths, pmax(runs$lengths, 10))
run_right <- 0.5 + cumsum(run_width)
run_left <- run_right - run_width
cursor <- max(run_right)
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
nodes_A$y_plot <- 4.25 + (max(nodes_A$x) - nodes_A$x) / max(nodes_A$x) * 3.3
edges_A <- nodes_A[nodes_A$node != nodes_A$parent, ]
edges_A$x_parent <- nodes_A$x_plot[match(edges_A$parent, nodes_A$node)]
edges_A$y_parent <- nodes_A$y_plot[match(edges_A$parent, nodes_A$node)]
edges_A$branch_color <- unname(branch_colors[as.character(edges_A$group)])

# One count panel under each cluster; the panel of a small cluster also uses the blank area of unclustered terms to its right.
panels_A <- dplyr::summarise(dplyr::group_by(dplyr::filter(tree_tips, as.character(group) != "0"), group), xmin = min(run_left) + 0.35, xmax = max(run_right) - 0.35, .groups = "drop")
panels_A <- dplyr::left_join(dplyr::mutate(panels_A, group = as.character(group)), cluster_definitions, by = "group")
next_run_group <- c(runs$values[-1], NA)
next_run_right <- c(run_right[-1], NA)
expand_run <- runs$values != "0" & runs$lengths <= 6 & next_run_group %in% "0"
panels_A$xmax[match(runs$values[expand_run], panels_A$group)] <- next_run_right[expand_run]
panels_A$bar_xmin <- panels_A$xmin + 0.35
panels_A$bar_xmax <- panels_A$xmax - 1.2

# One count has the same physical width in every panel (set by the panel that needs the most space per count).
requested_limit_A <- c(rRNA = max(sig_count$rRNA) + 3, Heat = max(sig_count$Heat) + 3, TF2 = max(sig_count$TF2) + 3, Mito = max(sig_count$Mito) + 3)[panels_A$cluster]
count_unit_A <- min((panels_A$bar_xmax - panels_A$bar_xmin) / requested_limit_A)
panels_A$limit <- floor((panels_A$bar_xmax - panels_A$bar_xmin) / count_unit_A + 1e-9)
panels_A$mid <- (panels_A$xmin + panels_A$xmax) / 2

# Bars and axis ticks of the count panels.
bars_A <- tidyr::pivot_longer(sig_count, -population, names_to = "cluster", values_to = "Count")
bars_A <- dplyr::left_join(bars_A, panels_A[, c("cluster", "bar_xmin", "bar_xmax", "fill", "limit")], by = "cluster")
bars_A$y <- unname(c(North = 2.55, Central = 1.55, South = 0.55)[as.character(bars_A$population)])
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
population_labels_A <- data.frame(population = c("North", "Central", "South"), y = c(2.55, 1.55, 0.55))

fig3A <- ggplot() +
  geom_segment(data = edges_A, aes(x = x_parent, xend = x_plot, y = y_parent, yend = y_parent), color = edges_A$branch_color, linewidth = 0.85, lineend = "butt") +
  geom_segment(data = edges_A, aes(x = x_plot, xend = x_plot, y = y_parent, yend = y_plot), color = edges_A$branch_color, linewidth = 0.85, lineend = "butt") +
  geom_rect(data = panels_A, aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = 4, fill = "white", color = panels_A$fill, linewidth = 1.1) +
  geom_rect(data = bars_A, aes(xmin = bar_xmin, xmax = xend, ymin = y - 0.30, ymax = y + 0.30), fill = bars_A$fill, color = NA) +
  geom_text(data = bars_A, aes(x = xend, y = y, label = Count), hjust = -0.25, size = 13, size.unit = "pt") +
  geom_text(data = panels_A, aes(x = mid, label = title), y = 3.55, color = panels_A$text_color, fontface = "bold", size = ifelse(panels_A$cluster == "Mito", 16, 17),
            size.unit = "pt", lineheight = 0.95) +
  geom_text(data = population_labels_A, aes(y = y, label = population), x = -0.8, hjust = 1, fontface = "bold", size = 17, size.unit = "pt") +
  geom_segment(data = ticks_A, aes(x = x, xend = x), y = 0, yend = -0.07, color = "grey40", linewidth = 0.3) +
  geom_text(data = ticks_A, aes(x = x, label = label), y = -0.23, size = 10, size.unit = "pt", color = "grey30") +
  annotate("text", x = (min(panels_A$xmin) + max(panels_A$xmax)) / 2, y = -0.62, label = "Number of significantly enriched GO terms", size = 12, size.unit = "pt") +
  scale_x_continuous(limits = c(-6.7, cursor + 0.2), expand = c(0, 0)) + scale_y_continuous(limits = c(-0.9, 7.8), expand = c(0, 0)) +
  coord_cartesian(clip = "off") + theme_void(base_family = "sans") + theme(plot.margin = margin(8, 10, 5, 10))
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

# Venn diagrams; font sizes in points (ggVennDiagram expects mm, so points are divided by .pt).

# Venn diagram: heat response genes.
venn_heat <- ggVennDiagram(list(North = genes_heat_North, Central = genes_heat_Central, South = genes_heat_South), label_alpha = 0, label_size = 12 / .pt, set_size = 16 / .pt, edge_size = 0.65) +
  scale_fill_gradientn(colours = c("#FEF8F1", "#FDEBD4", "#FDBF6F")) + scale_x_continuous(expand = expansion(mult = c(0.09, 0.09))) + ggtitle("Heat response genes") +
  theme_void(base_family = "sans") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17, color = "#918009", lineheight = 0.95), legend.position = "none",
                                           plot.margin = margin(8, 4, 6, 4)) + coord_fixed()
venn_heat$layers[[3]]$aes_params$fontface <- "bold"

# Venn diagram: ribosome-related processes.
venn_rRNA <- ggVennDiagram(list(North = genes_rRNA_North, Central = genes_rRNA_Central, South = genes_rRNA_South), label_alpha = 0, label_size = 12 / .pt, set_size = 16 / .pt, edge_size = 0.65) +
  scale_fill_gradientn(colours = c("#F4F9ED", "#E4F1D6", "#B2DF8A")) + scale_x_continuous(expand = expansion(mult = c(0.09, 0.09))) + ggtitle("Ribosome-related processes") +
  theme_void(base_family = "sans") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17, color = "#3A6831", lineheight = 0.95), legend.position = "none",
                                           plot.margin = margin(8, 4, 6, 4)) + coord_fixed()
venn_rRNA$layers[[3]]$aes_params$fontface <- "bold"

# Venn diagram: mitochondrial respiration.
venn_mito <- ggVennDiagram(list(North = genes_mito_North, Central = genes_mito_Central, South = genes_mito_South), label_alpha = 0, label_size = 12 / .pt, set_size = 16 / .pt, edge_size = 0.65) +
  scale_fill_gradientn(colours = c("#F4F9FC", "#E3EFF6", "#A6CEE3")) + scale_x_continuous(expand = expansion(mult = c(0.09, 0.09))) + ggtitle("Mitochondrial respiration") +
  theme_void(base_family = "sans") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17, color = "#82B6D1", lineheight = 0.95), legend.position = "none",
                                           plot.margin = margin(8, 4, 6, 4)) + coord_fixed()
venn_mito$layers[[3]]$aes_params$fontface <- "bold"

# Venn diagram: Tf2-like.
venn_TF2 <- ggVennDiagram(list(North = genes_TF2_North, Central = genes_TF2_Central, South = genes_TF2_South), label_alpha = 0, label_size = 12 / .pt, set_size = 16 / .pt, edge_size = 0.65) +
  scale_fill_gradientn(colours = c("#FAF7FB", "#EFE8F2", "#CAB2D6")) + scale_x_continuous(expand = expansion(mult = c(0.09, 0.09))) + ggtitle("Tf2-like") +
  theme_void(base_family = "sans") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17, color = "#CF92ED", lineheight = 0.95), legend.position = "none",
                                           plot.margin = margin(8, 4, 6, 4)) + coord_fixed()
venn_TF2$layers[[3]]$aes_params$fontface <- "bold"

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


#################### Figure 3C: front-loading of shared heat response genes

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

fig3C <- ggplot(frontload) + geom_vline(xintercept = 0, color = "black", linetype = "dashed") + geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(aes(x = baseline, y = plasticity, shape = type, size = type, fill = comp, color = comp), stroke = 0.6) +
  geom_smooth(aes(x = baseline, y = plasticity, color = comp, fill = comp), method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.1) +
  labs(x = "Baseline expression difference", y = "Plasticity difference", color = "Population comparison", fill = "Population comparison", shape = "Gene category") +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"), axis.title = element_text(color = "black", size = 16),
        axis.text = element_text(color = "black", size = 11), legend.title = element_text(size = 15), legend.text = element_text(size = 13),
        legend.position = "inside", legend.position.inside = c(0.2, 0.2), legend.spacing.y = unit(0.05, "cm"), plot.title = element_text(size = 11, hjust = 0.5)) +
  scale_color_manual(values = c("North vs South" = "#1B9E77", "Central vs South" = "#7570B3", "Central vs North" = "#E7298A")) +
  scale_fill_manual(values = c("North vs South" = scales::alpha("#1B9E77", 0.65), "Central vs South" = scales::alpha("#7570B3", 0.65), "Central vs North" = scales::alpha("#E7298A", 0.65))) +
  scale_size_manual(values = c("Hsp70" = 5, "Other genes" = 3)) + scale_shape_manual(values = c("Hsp70" = 24, "Other genes" = 21)) +
  guides(size = "none", shape = guide_legend(override.aes = list(size = c("Hsp70" = 3, "Other genes" = 2))))
fig3C


#################### Figure 3D: baseline expression and plasticity of Hsp70 genes in each population

# Per Hsp70 gene and population: baseline = mean control expression; plasticity = heat minus control (log2FC).
hsp70_expr <- expr_mean[expr_mean$gene %in% hsp70_genes, ]
hsp70_expr$plasticity <- hsp70_expr$heat - hsp70_expr$control
hsp70_expr$population <- factor(hsp70_expr$population, levels = c("North", "Central", "South"))

# Pairwise Wilcoxon tests between populations (BH correction); significant pairs (padj < 0.05) are marked on the boxplots.
test_baseline <- pairwise.wilcox.test(x = hsp70_expr$control, g = hsp70_expr$population, paired = FALSE, p.adjust.method = "BH")
test_plasticity <- pairwise.wilcox.test(x = hsp70_expr$plasticity, g = hsp70_expr$population, paired = FALSE, p.adjust.method = "BH")
test_baseline
test_plasticity

# Significant pairs and their star labels, with brackets stacked above the highest box.
sig_baseline <- as.data.frame(as.table(test_baseline$p.value))
colnames(sig_baseline) <- c("group2", "group1", "p.adj")
sig_baseline <- dplyr::mutate(dplyr::filter(sig_baseline, !is.na(p.adj), p.adj < 0.05), group1 = as.character(group1), group2 = as.character(group2),
                              label = dplyr::case_when(p.adj < 0.0001 ~ "****", p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**", TRUE ~ "*"))
sig_baseline$y.position <- max(hsp70_expr$control, na.rm = TRUE) + seq_len(nrow(sig_baseline)) * 0.5

sig_plasticity <- as.data.frame(as.table(test_plasticity$p.value))
colnames(sig_plasticity) <- c("group2", "group1", "p.adj")
sig_plasticity <- dplyr::mutate(dplyr::filter(sig_plasticity, !is.na(p.adj), p.adj < 0.05), group1 = as.character(group1), group2 = as.character(group2),
                                label = dplyr::case_when(p.adj < 0.0001 ~ "****", p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**", TRUE ~ "*"))
sig_plasticity$y.position <- max(hsp70_expr$plasticity, na.rm = TRUE) + seq_len(nrow(sig_plasticity)) * 0.8

# Shared style of the two boxplots.
theme_grid <- theme(panel.background = element_rect(fill = "white", color = NA), panel.grid.major = element_line(color = "grey96"), panel.grid.minor = element_line(color = "grey96"))
theme_box <- theme(legend.position = "none", strip.background = element_rect(fill = "#D3D3D3", color = "black"), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
                   strip.text = element_text(size = 14, face = "bold"), axis.title = element_text(size = 16, face = "bold"),
                   axis.text.x = element_text(size = 16, face = "bold", color = "black"), axis.text.y = element_text(size = 11))

# Boxplot: Hsp70 baseline expression.
box_baseline <- ggplot(hsp70_expr, aes(x = population, y = control)) + geom_boxplot(width = 0.7, fill = "#fdbf6f", color = "black", alpha = 0.9, outlier.shape = NA) +
  facet_wrap(~ "Hsp70 baseline expression", nrow = 1) + labs(x = "", y = expression(log[2]*CPM)) + theme_bw(base_size = 16) + theme_grid + theme_box +
  expand_limits(y = 0) + scale_y_continuous(breaks = scales::breaks_width(2), expand = expansion(mult = c(0.02, 0.08)))
if (nrow(sig_baseline) > 0) box_baseline <- box_baseline + stat_pvalue_manual(sig_baseline, label = "label", y.position = "y.position", tip.length = 0.01, bracket.size = 0.4, size = 14 / .pt, hide.ns = TRUE)

# Boxplot: Hsp70 plasticity.
box_plasticity <- ggplot(hsp70_expr, aes(x = population, y = plasticity)) + geom_boxplot(width = 0.7, fill = "#fdbf6f", color = "black", alpha = 0.9, outlier.shape = NA) +
  facet_wrap(~ "Hsp70 plasticity", nrow = 1) + labs(x = "", y = expression(log[2]*FC)) + theme_bw(base_size = 16) + theme_grid + theme_box +
  scale_y_continuous(breaks = scales::breaks_width(2), expand = expansion(mult = c(0.02, 0.08)))
if (nrow(sig_plasticity) > 0) box_plasticity <- box_plasticity + stat_pvalue_manual(sig_plasticity, label = "label", y.position = "y.position", tip.length = 0.01, bracket.size = 0.4, size = 14 / .pt, hide.ns = TRUE)

fig3D <- (box_baseline | box_plasticity) + plot_layout(widths = c(1, 1))
fig3D


#################### Combine Figure 3A-D (heights in inches: A 5.4, B 4.6, C+D 6.0)

fig3 <- ggarrange(fig3A, fig3B, ggarrange(fig3C, fig3D, labels = c("C", "D"), ncol = 2, font.label = list(size = 20, face = "bold")),
                  labels = c("A", "B", ""), ncol = 1, nrow = 3, heights = c(5.4, 4.6, 6.0), font.label = list(size = 20, face = "bold"))
ggsave("Figure3.pdf", fig3, width = 14, height = 16, device = cairo_pdf, bg = "white")


############################################################
#Part 4: up-regulated GO terms in each population, ordered by the GO tree (Figure S4, Table S5)

#################### Figure S4

# Tree coordinates of the same GO tree and clusters as Figure 3A, with the tree width rescaled to 4.7 to fit the 10 x 12 inch figure.
nodes_S4 <- as.data.frame(tree_ordered$data)
nodes_S4$x <- nodes_S4$x / max(nodes_S4$x) * 4.7
tips_S4 <- nodes_S4[nodes_S4$isTip, ]
tips_S4 <- tips_S4[order(tips_S4$y), ]
tips_S4$GO <- go_up_all@result$ID[match(tips_S4$label, go_up_all@result$Description)]
N_GO_S4 <- nrow(tips_S4)

# Tree branches (horizontal and vertical segments), coloured by cluster.
edges_S4 <- nodes_S4[nodes_S4$node != nodes_S4$parent, ]
edges_S4$x_parent <- nodes_S4$x[match(edges_S4$parent, nodes_S4$node)]
edges_S4$y_parent <- nodes_S4$y[match(edges_S4$parent, nodes_S4$node)]
edges_S4$branch_color <- unname(branch_colors[as.character(edges_S4$group)])

# Bubbles: number of up-regulated genes in each GO term (numerator of GeneRatio) and FDR (p.adjust), per population.
go_data_S4 <- go_up_compare[, c("ID", "Description")]
go_data_S4$Count_North <- as.numeric(sub("/.*", "", go_up_compare$GeneRatio_North))
go_data_S4$FDR_North <- go_up_compare$padj_North
go_data_S4$Count_Central <- as.numeric(sub("/.*", "", go_up_compare$GeneRatio_Central))
go_data_S4$FDR_Central <- go_up_compare$padj_Central
go_data_S4$Count_South <- as.numeric(sub("/.*", "", go_up_compare$GeneRatio_South))
go_data_S4$FDR_South <- go_up_compare$padj_South
bubble_S4 <- dplyr::bind_rows(dplyr::transmute(go_data_S4, ID, population = "North", Count = Count_North, FDR = FDR_North),
                              dplyr::transmute(go_data_S4, ID, population = "Central", Count = Count_Central, FDR = FDR_Central),
                              dplyr::transmute(go_data_S4, ID, population = "South", Count = Count_South, FDR = FDR_South))
bubble_S4$y <- tips_S4$y[match(bubble_S4$ID, tips_S4$GO)]
bubble_S4$x <- unname(c(North = 6.925, Central = 8.675, South = 10.425)[bubble_S4$population])

# Bubble colour shows -log10 FDR only for significant terms (FDR < 0.01); non-significant terms are grey; untested terms are not drawn.
bubble_S4$ColorCol <- ifelse(!is.na(bubble_S4$FDR) & bubble_S4$FDR < 0.01, -log10(pmax(bubble_S4$FDR, .Machine$double.xmin)), NA_real_)
bubble_S4 <- bubble_S4[!is.na(bubble_S4$Count), ]

# Coloured background behind the GO IDs of each cluster (the tree itself stays on white).
label_bands_S4 <- dplyr::summarise(dplyr::group_by(dplyr::filter(tips_S4, as.character(group) != "0"), group), ymin = min(y) - 0.5, ymax = max(y) + 0.5, .groups = "drop")
label_bands_S4$fill <- unname(branch_colors[as.character(label_bands_S4$group)])

# Column headers for the three populations.
headers_S4 <- data.frame(population = c("North", "Central", "South"), xmin = c(6.05, 7.80, 9.55), xmax = c(7.80, 9.55, 11.30), x = c(6.925, 8.675, 10.425))

# Cluster names on the right.
annotations_S4 <- label_bands_S4
annotations_S4$label <- cluster_definitions$title[match(as.character(annotations_S4$group), cluster_definitions$group)]
annotations_S4$text_color <- cluster_definitions$text_color[match(as.character(annotations_S4$group), cluster_definitions$group)]
annotations_S4$mid <- (annotations_S4$ymin + annotations_S4$ymax) / 2
annotations_S4$x_text <- 12.17

# A further group of GO terms outside the four clusters, labelled only (not a separate cluster).
viral_y_S4 <- range(tips_S4$y[tips_S4$GO %in% c("GO:0044185", "GO:0042575")])

# Tree, GO IDs, bubbles and labels share one y axis, so every row lines up.
figS4 <- ggplot() +
  geom_segment(data = edges_S4, aes(x = x_parent, xend = x, y = y, yend = y), color = edges_S4$branch_color, linewidth = 0.85, lineend = "butt") +
  geom_segment(data = edges_S4, aes(x = x_parent, xend = x_parent, y = y_parent, yend = y), color = edges_S4$branch_color, linewidth = 0.85, lineend = "butt") +
  geom_rect(data = label_bands_S4, aes(ymin = ymin, ymax = ymax), xmin = 4.8, xmax = 6.05, fill = label_bands_S4$fill, alpha = 0.45, color = NA) +
  geom_segment(data = tips_S4, aes(y = y, yend = y), x = 6.05, xend = 11.30, color = "grey90", linewidth = 0.25) +
  geom_text(data = tips_S4, aes(y = y, label = GO), x = (4.8 + 6.05) / 2, hjust = 0.5, size = 9, size.unit = "pt", color = "black") +
  geom_rect(data = headers_S4, aes(xmin = xmin, xmax = xmax), ymin = N_GO_S4 + 0.65, ymax = N_GO_S4 + 2.65, fill = "grey90", color = "grey35", linewidth = 0.3) +
  geom_text(data = headers_S4, aes(x = x, label = population), y = N_GO_S4 + 1.65, size = 11, size.unit = "pt") +
  geom_point(data = bubble_S4, aes(x = x, y = y, size = Count, color = ColorCol)) +
  scale_size_continuous(range = c(1, 3.6), breaks = c(25, 50, 75, 100, 125), name = "Gene count",
                        guide = guide_legend(order = 2, nrow = 1, title.position = "left", override.aes = list(shape = 21, fill = NA, color = "black", stroke = 0.4))) +
  scale_color_gradient(low = "#97D1A0", high = "#1E803D", na.value = "grey80", name = expression(plain("-log")[10]~FDR), breaks = c(50, 100),
                       guide = guide_colorbar(order = 1, direction = "horizontal", title.position = "left", barwidth = grid::unit(1.3, "cm"), barheight = grid::unit(0.3, "cm"))) +
  geom_segment(data = annotations_S4, aes(y = ymin + 0.2, yend = ymax - 0.2), x = 11.88, xend = 11.88, color = annotations_S4$fill, linewidth = 0.6) +
  geom_text(data = annotations_S4, aes(x = x_text, y = mid, label = label), angle = 90, size = 11, size.unit = "pt", fontface = "bold", lineheight = 0.95, color = annotations_S4$text_color) +
  annotate("segment", x = 12.46, xend = 12.46, y = viral_y_S4[1], yend = viral_y_S4[2], color = "black", linewidth = 0.5) +
  annotate("text", x = 12.80, y = mean(viral_y_S4), label = "Viral processes and DNA metabolism", angle = 90, size = 11, size.unit = "pt", fontface = "bold", color = "black") +
  scale_x_continuous(limits = c(-0.1, 13.6), expand = c(0, 0)) + scale_y_continuous(limits = c(-5.7, N_GO_S4 + 4), expand = c(0, 0)) +
  theme_void(base_family = "sans") +
  theme(legend.position = "inside", legend.position.inside = c(0.63, 0.015), legend.justification = c(0.5, 0), legend.box = "horizontal",
        legend.title = element_text(size = 9), legend.text = element_text(size = 8), legend.key.width = grid::unit(0.4, "cm"),
        legend.spacing.x = grid::unit(1, "cm"), plot.margin = margin(13, 12, 5, 10), plot.background = element_rect(fill = "white", color = NA))
figS4
ggsave("FigureS4.pdf", figS4, width = 10, height = 12, device = cairo_pdf, bg = "white")


#################### Table S5: gene count and FDR of each GO term per population, in the order of Figure S4 (top to bottom)

table_S5 <- go_data_S4[match(rev(tips_S4$GO), go_data_S4$ID), c("ID", "Description", "Count_North", "FDR_North", "Count_Central", "FDR_Central", "Count_South", "FDR_South")]
colnames(table_S5) <- c("GO_ID", "GO_term_name", "Gene_count_North", "FDR_North", "Gene_count_Central", "FDR_Central", "Gene_count_South", "FDR_South")

# FDR below 1E-06 in scientific notation (3 significant digits), otherwise 6 decimals.
table_S5$FDR_North <- ifelse(is.na(table_S5$FDR_North), NA, ifelse(table_S5$FDR_North < 1e-6, toupper(formatC(table_S5$FDR_North, format = "e", digits = 2)), formatC(table_S5$FDR_North, format = "f", digits = 6)))
table_S5$FDR_Central <- ifelse(is.na(table_S5$FDR_Central), NA, ifelse(table_S5$FDR_Central < 1e-6, toupper(formatC(table_S5$FDR_Central, format = "e", digits = 2)), formatC(table_S5$FDR_Central, format = "f", digits = 6)))
table_S5$FDR_South <- ifelse(is.na(table_S5$FDR_South), NA, ifelse(table_S5$FDR_South < 1e-6, toupper(formatC(table_S5$FDR_South, format = "e", digits = 2)), formatC(table_S5$FDR_South, format = "f", digits = 6)))
head(table_S5[, c("GO_ID", "GO_term_name")], 10)
dim(table_S5)
write.csv(table_S5, "Table_S5.csv", row.names = FALSE, na = "NA")

