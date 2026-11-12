setwd("C:/path/to/RNAgithub/10TF2/")

# =============================================================================
# Heat-induced Tf2-like retrotransposons (Figure 5 and Figures S13, S14, S16-S18) and expression of the HSF1-like gene (Figure S15).
# Steps run outside R (LTR_retriever, bedtools, BLASTN, MAFFT, TEsorter, FIMO, read mapping) are given as comments; their output files are read here.
# Output: Figure5.pdf, FigureS13.pdf-FigureS18.pdf, and Figure5_HSF1_MA04862_positions.tsv (HSF1 site and TATA box in Figure 5B).
# =============================================================================

library(rtracklayer)
library(GenomicRanges)
library(Biostrings)
library(DECIPHER)
library(ape)
library(phangorn)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


############################################################
#Part 1: heat-induced Tf2-like genes and the intact LTR elements that contain them

# Up-regulated general temperature response genes (padj < 0.05, log2FC >= 0.5) in at least one population.
deg_North <- read.delim("../05DEseq2/North_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
deg_Central <- read.delim("../05DEseq2/Central_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
deg_South <- read.delim("../05DEseq2/South_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
up_North <- deg_North$gene[which(deg_North$padj < 0.05 & deg_North$log2FoldChange >= 0.5)]
up_Central <- deg_Central$gene[which(deg_Central$padj < 0.05 & deg_Central$log2FoldChange >= 0.5)]
up_South <- deg_South$gene[which(deg_South$padj < 0.05 & deg_South$log2FoldChange >= 0.5)]
up_union <- unique(c(up_North, up_Central, up_South))

# Up-regulated genes annotated as Tf2-like retrotransposons (125 genes).
annotation <- read.delim("../06GOKEGG/gene_annotation.txt", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
tf2_genes <- unique(annotation$gene[annotation$gene %in% up_union & grepl("Tf2-", annotation$description, ignore.case = FALSE)])
length(tf2_genes)

# Genomic position of each gene = span of its exons in the StringTie annotation.
gtf <- import("../06GOKEGG/Lsax_stringtie_merged.gtf")
gtf_exon <- gtf[gtf$type == "exon" & !is.na(gtf$gene_id) & gtf$gene_id != ""]
gene_ranges <- unlist(range(split(gtf_exon, gtf_exon$gene_id)), use.names = FALSE)
gene_ranges$gene_id <- names(split(gtf_exon, gtf_exon$gene_id))
tf2_gene_ranges <- gene_ranges[gene_ranges$gene_id %in% tf2_genes]
length(tf2_gene_ranges)

# Intact LTR retrotransposons in the reference genome (GCA_037325665.1), identified beforehand with LTRharvest and LTR_retriever.
# repeat_region = whole element; long_terminal_repeat = its two LTRs.
gff <- read.delim("GCA_037325665.1_US_GU_Lsax_2.0_genomic.fna.mod.pass.list.gff3", header = FALSE, quote = "", comment.char = "#", stringsAsFactors = FALSE)
colnames(gff) <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attribute")
repeat_region <- gff %>% dplyr::filter(type == "repeat_region") %>%
  dplyr::mutate(intact_id = paste0(chr, ":", start, "-", end), parent = sub(".*ID=([^;]+).*", "\\1", attribute))
ltr_region <- gff %>% dplyr::filter(type == "long_terminal_repeat") %>%
  dplyr::transmute(parent = sub(".*Parent=([^;]+).*", "\\1", attribute), ltr_id = sub(".*ID=([^;]+).*", "\\1", attribute), ltr_start = start, ltr_end = end)
intact_gr <- GRanges(repeat_region$chr, IRanges(repeat_region$start, repeat_region$end), intact_id = repeat_region$intact_id)

# Overlap between the Tf2 genes and the intact elements.
hits_gene_intact <- findOverlaps(tf2_gene_ranges, intact_gr, ignore.strand = TRUE)
gene_intact <- data.frame(gene = tf2_gene_ranges$gene_id[queryHits(hits_gene_intact)], intact_id = intact_gr$intact_id[subjectHits(hits_gene_intact)],
                          intact_length = width(intact_gr)[subjectHits(hits_gene_intact)],
                          overlap_bp = width(pintersect(tf2_gene_ranges[queryHits(hits_gene_intact)], intact_gr[subjectHits(hits_gene_intact)])))

# Intact elements of which at least 80% is covered by a Tf2 gene: 47 elements.
intact_ids_80 <- unique(gene_intact$intact_id[gene_intact$overlap_bp / gene_intact$intact_length >= 0.8])
length(intact_ids_80)
intact_gr_80 <- intact_gr[intact_gr$intact_id %in% intact_ids_80]
intact_gr_80$name <- intact_gr_80$intact_id
intact_gr_80$score <- 0
#export(intact_gr_80, "TF2_intact_LTR_80percent_overlap.bed", format = "BED")
# Run on Linux: bedtools getfasta with this BED file -> TF2_intact_LTR_80percent_overlap.fa

#Part 1B: overlap of each of the 125 Tf2 genes with the intact elements (covers >= 80% of an intact element,
# lies inside an intact element, other partial overlap, or no overlap).
tf2_range_index <- match(tf2_genes, tf2_gene_ranges$gene_id)
tf2_ltr_gene_summary <- data.frame(
  gene = tf2_genes, chr = as.character(seqnames(tf2_gene_ranges))[tf2_range_index],
  gene_start = start(tf2_gene_ranges)[tf2_range_index], gene_end = end(tf2_gene_ranges)[tf2_range_index],
  gene_strand = as.character(strand(tf2_gene_ranges))[tf2_range_index],
  gene_span_bp = width(tf2_gene_ranges)[tf2_range_index], stringsAsFactors = FALSE)
tf2_ltr_overlap_detail <- gene_intact %>%
  dplyr::left_join(tf2_ltr_gene_summary, by = "gene") %>%
  dplyr::mutate(element_coverage_pct = 100 * overlap_bp / intact_length,
                gene_span_coverage_pct = 100 * overlap_bp / gene_span_bp,
                intact_strand = repeat_region$strand[match(intact_id, repeat_region$intact_id)],
                selected_seed_element = intact_id %in% intact_ids_80)
tf2_ltr_best_hit <- tf2_ltr_overlap_detail %>%
  dplyr::group_by(gene) %>%
  dplyr::arrange(dplyr::desc(element_coverage_pct), dplyr::desc(overlap_bp), intact_id, .by_group = TRUE) %>%
  dplyr::mutate(n_intact_overlaps = dplyr::n()) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>%
  dplyr::select(gene, n_intact_overlaps, best_intact_id = intact_id, intact_length,
                overlap_bp, element_coverage_pct, gene_span_coverage_pct, selected_seed_element)
tf2_ltr_gene_summary <- tf2_ltr_gene_summary %>%
  dplyr::left_join(tf2_ltr_best_hit, by = "gene") %>%
  dplyr::mutate(n_intact_overlaps = tidyr::replace_na(n_intact_overlaps, 0L),
                category = dplyr::case_when(
                  n_intact_overlaps == 0L ~ "no_overlap_with_identified_intact_LTR",
                  element_coverage_pct >= 80 ~ "covers_at_least_80pct_of_intact_LTR",
                  gene_span_coverage_pct >= 100 - 1e-8 ~ "partial_gene_model_contained_in_intact_LTR",
                  TRUE ~ "other_partial_overlap_with_intact_LTR"))
tf2_seed_gene_ids <- tf2_ltr_gene_summary$gene[tf2_ltr_gene_summary$category == "covers_at_least_80pct_of_intact_LTR"]
tf2_remaining_gene_ids <- setdiff(tf2_genes, tf2_seed_gene_ids)
tf2_partial_gene_models <- tf2_ltr_gene_summary %>%
  dplyr::filter(category == "partial_gene_model_contained_in_intact_LTR")
tf2_fragment_summary <- data.frame(
  measure = c("all_upregulated_Tf2_genes", "selected_seed_intact_elements", "genes_covering_ge80pct_of_an_intact_element",
              "remaining_genes", "remaining_genes_contained_as_partial_models_in_intact_elements",
              "distinct_intact_elements_containing_these_partial_models", "of_those_elements_also_in_47_seed_set",
              "remaining_genes_with_other_partial_overlap", "remaining_genes_without_any_intact_overlap"),
  count = c(length(tf2_genes), length(intact_ids_80), length(tf2_seed_gene_ids), length(tf2_remaining_gene_ids),
            nrow(tf2_partial_gene_models), length(unique(tf2_partial_gene_models$best_intact_id)),
            length(intersect(unique(tf2_partial_gene_models$best_intact_id), intact_ids_80)),
            sum(tf2_ltr_gene_summary$category == "other_partial_overlap_with_intact_LTR"),
            sum(tf2_ltr_gene_summary$category == "no_overlap_with_identified_intact_LTR")))

print(tf2_fragment_summary, row.names = FALSE)
print(tf2_partial_gene_models, row.names = FALSE)

#Part 1C: BLASTN of the 125 Tf2 genes (genomic span from the first to the last exon) against the two final
# consensus sequences built in Part 3 (Table S8). Run on Linux (BLAST+ 2.17.0):
# blastn -task blastn -query Tf2_125_gene_spans.normalized.fa -subject consensus_cladeA_cladeB.fa \
#   -strand both -evalue 1e-10 -dust yes -soft_masking true -max_hsps 100 \
#   -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident gaps sstrand" \
#   -out Tf2_125_vs_consensus.blastn.tsv
# Hits are kept with E-value <= 1e-10, identity >= 80% and alignment length >= 80 bp. Overlapping hits are merged for
# coverage; identity is the mean over the hits, weighted by hit length. Best consensus = highest gene coverage, then highest bit score.

tf2_consensus_blast_file <- "Tf2_125_vs_consensus.blastn.tsv"
tf2_consensus_blast_segments <- read.delim(tf2_consensus_blast_file, header = FALSE, stringsAsFactors = FALSE,
                                           col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "nident", "gaps", "sstrand"))
tf2_consensus_qualified <- tf2_consensus_blast_segments %>%
  dplyr::filter(evalue <= 1e-10, pident >= 80, length >= 80)

tf2_consensus_pair_metrics <- tf2_consensus_qualified %>% dplyr::group_by(qseqid, sseqid) %>%
  dplyr::summarise(qualified_blast_segments = dplyr::n(),
                   gene_span_bp = dplyr::first(qlen), consensus_bp = dplyr::first(slen),
                   query_covered_bp = sum(width(IRanges::reduce(IRanges(pmin(qstart, qend), pmax(qstart, qend))))),
                   consensus_covered_bp = sum(width(IRanges::reduce(IRanges(pmin(sstart, send), pmax(sstart, send))))),
                   query_coverage_pct = 100 * query_covered_bp / gene_span_bp,
                   consensus_coverage_pct = 100 * consensus_covered_bp / consensus_bp,
                   blast_segment_length_weighted_identity_pct = 100 * sum(nident) / sum(length),
                   best_blast_segment_bitscore = max(bitscore), min_evalue = min(evalue), .groups = "drop")
tf2_consensus_pooled_coverage <- tf2_consensus_qualified %>% dplyr::group_by(qseqid) %>%
  dplyr::summarise(pooled_query_coverage_pct = 100 * sum(width(IRanges::reduce(IRanges(pmin(qstart, qend), pmax(qstart, qend))))) / dplyr::first(qlen), .groups = "drop")

# Matches of each gene to the two consensus sequences, and its best match.
tf2_consensus_per_gene <- tf2_consensus_pair_metrics %>% dplyr::group_by(qseqid) %>%
  dplyr::summarise(matched_consensuses = paste(sort(sseqid), collapse = ","), n_consensuses = dplyr::n(),
                   max_single_consensus_query_coverage_pct = max(query_coverage_pct),
                   any_query_coverage_ge80pct = any(query_coverage_pct >= 80), .groups = "drop") %>%
  dplyr::left_join(tf2_consensus_pooled_coverage, by = "qseqid")
tf2_consensus_best_match <- tf2_consensus_pair_metrics %>%
  dplyr::arrange(qseqid, dplyr::desc(query_coverage_pct), dplyr::desc(best_blast_segment_bitscore), sseqid) %>%
  dplyr::group_by(qseqid) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>%
  dplyr::transmute(qseqid, best_consensus = sseqid,
                   best_query_covered_bp = query_covered_bp,
                   best_consensus_covered_bp = consensus_covered_bp,
                   best_query_coverage_pct = query_coverage_pct,
                   best_consensus_coverage_pct = consensus_coverage_pct,
                   best_blast_segment_length_weighted_identity_pct = blast_segment_length_weighted_identity_pct,
                   best_qualified_blast_segments = qualified_blast_segments, best_blast_segment_bitscore,
                   best_min_evalue = min_evalue)
tf2_consensus_gene_homology <- tf2_ltr_gene_summary %>%
  dplyr::left_join(tf2_consensus_per_gene, by = c("gene" = "qseqid")) %>%
  dplyr::left_join(tf2_consensus_best_match, by = c("gene" = "qseqid")) %>%
  dplyr::mutate(n_consensuses = tidyr::replace_na(n_consensuses, 0L),
                any_query_coverage_ge80pct = tidyr::replace_na(any_query_coverage_ge80pct, FALSE),
                homology_class = dplyr::case_when(n_consensuses == 0L ~ "none", n_consensuses > 1L ~ "both", TRUE ~ best_consensus))

tf2_consensus_gene_sets <- list(all125 = tf2_consensus_gene_homology,
                                seed47 = tf2_consensus_gene_homology[tf2_consensus_gene_homology$category == "covers_at_least_80pct_of_intact_LTR", ],
                                partial_models5 = tf2_consensus_gene_homology[tf2_consensus_gene_homology$category == "partial_gene_model_contained_in_intact_LTR", ],
                                remaining78 = tf2_consensus_gene_homology[tf2_consensus_gene_homology$gene %in% tf2_remaining_gene_ids, ],
                                no_intact_overlap73 = tf2_consensus_gene_homology[tf2_consensus_gene_homology$category == "no_overlap_with_identified_intact_LTR", ])
tf2_consensus_summary_R <- do.call(rbind, lapply(names(tf2_consensus_gene_sets), function(nm) {
  dat <- tf2_consensus_gene_sets[[nm]]
  data.frame(group = nm, n_genes = nrow(dat), n_qualified_match = sum(dat$n_consensuses > 0),
             n_query_coverage_ge80pct = sum(dat$any_query_coverage_ge80pct),
             n_pooled_query_coverage_ge80pct = sum(dat$pooled_query_coverage_pct >= 80, na.rm = TRUE),
             n_matches_both_consensuses = sum(dat$n_consensuses > 1))
}))

# Identity and aligned length of the best matches, as reported in the Results.
tf2_best_matched_genes <- tf2_consensus_gene_homology %>% dplyr::filter(n_consensuses > 0L)
tf2_consensus_identity_length_summary_R <- data.frame(
  n_genes = nrow(tf2_consensus_gene_homology),
  n_matched_genes = nrow(tf2_best_matched_genes),
  n_matches_both_consensuses = sum(tf2_consensus_gene_homology$n_consensuses > 1L),
  identity_min_pct = min(tf2_best_matched_genes$best_blast_segment_length_weighted_identity_pct),
  identity_max_pct = max(tf2_best_matched_genes$best_blast_segment_length_weighted_identity_pct),
  identity_median_pct = median(tf2_best_matched_genes$best_blast_segment_length_weighted_identity_pct),
  aligned_query_bp_min = min(tf2_best_matched_genes$best_query_covered_bp),
  aligned_query_bp_max = max(tf2_best_matched_genes$best_query_covered_bp),
  n_best_query_coverage_ge80pct = sum(tf2_consensus_gene_homology$best_query_coverage_pct >= 80, na.rm = TRUE),
  pct_best_query_coverage_ge80pct = 100 * sum(tf2_consensus_gene_homology$best_query_coverage_pct >= 80, na.rm = TRUE) / nrow(tf2_consensus_gene_homology))
print(tf2_consensus_summary_R, row.names = FALSE)
print(tf2_consensus_identity_length_summary_R, row.names = FALSE)


############################################################
#Part 2: tree of the 47 elements, the two clades and their first consensus sequences (Figure S13)

# Elements on the minus strand are reverse-complemented so that all sequences have the same orientation.
seed_seq <- readDNAStringSet("TF2_intact_LTR_80percent_overlap.fa")
names(seed_seq) <- sub("::.*", "", names(seed_seq))
seed_strand <- repeat_region$strand[match(names(seed_seq), repeat_region$intact_id)]
table(seed_strand)
seed_seq[seed_strand == "-"] <- reverseComplement(seed_seq[seed_strand == "-"])

# Neighbour-joining tree from the alignment, rooted at the midpoint; branches are rotated (topology unchanged) so that the
# long elements (> 7 kb) are drawn at the top.
seed_aligned <- AlignSeqs(seed_seq)
seed_distance <- DistanceMatrix(seed_aligned, type = "dist")
tree_mid <- ladderize(midpoint(nj(as.dist(seed_distance))))
long_ids <- names(seed_seq)[width(seed_seq) > 7000]
tree_rot <- rotateConstr(tree_mid, c(long_ids, setdiff(tree_mid$tip.label, long_ids)))

# The two main clades are defined by two elements at either end of each clade: all tips below their most recent common
# ancestor belong to the clade (clade A: 20 elements, clade B: 26 elements).
anchor_A <- c("CM074560.1:68560024-68566712", "CM074569.1:37319988-37325259")
anchor_B <- c("CM074567.1:11604713-11608647", "CM074569.1:1253935-1259037")
seed_clade_A <- extract.clade(tree_rot, getMRCA(tree_rot, anchor_A))$tip.label
seed_clade_B <- extract.clade(tree_rot, getMRCA(tree_rot, anchor_B))$tip.label
length(seed_clade_A)
length(seed_clade_B)

# Figure S13
clade_colors <- c(cladeA_cons = "#D55E00", cladeB_cons = "#0072B2")
tip_colors <- ifelse(tree_rot$tip.label %in% seed_clade_A, clade_colors[["cladeA_cons"]], ifelse(tree_rot$tip.label %in% seed_clade_B, clade_colors[["cladeB_cons"]], "grey40"))
cairo_pdf("FigureS13.pdf", width = 11, height = 13)
par(mar = c(2, 1, 1, 1), family = "sans")
plot(tree_rot, type = "phylogram", cex = 0.85, font = 1, no.margin = FALSE, tip.color = tip_colors)
add.scale.bar(x = 0, y = 0, cex = 0.9)
legend("topleft", legend = c("Clade A", "Clade B", "Other"), col = c(clade_colors, "grey40"), pch = 16, bty = "n", cex = 1)
dev.off()

# First consensus of each clade (majority base at each alignment column, threshold 0.5).
consensus1_A <- gsub("-", "", consensusString(AlignSeqs(seed_seq[seed_clade_A]), ambiguityMap = "N", threshold = 0.5))
consensus1_B <- gsub("-", "", consensusString(AlignSeqs(seed_seq[seed_clade_B]), ambiguityMap = "N", threshold = 0.5))
consensus1 <- DNAStringSet(c(cladeA_cons = consensus1_A, cladeB_cons = consensus1_B))
width(consensus1)
#writeXStringSet(consensus1, "TF2_clade_consensus.fa")


############################################################
#Part 3: all intact copies of clade A and clade B in the genome, and the final consensus sequences

# All intact elements of the genome, saved as BED to extract their sequences.
all_intact_gr <- GRanges(repeat_region$chr, IRanges(repeat_region$start, repeat_region$end), strand = ifelse(repeat_region$strand %in% c("+", "-"), repeat_region$strand, "*"))
all_intact_gr$name <- repeat_region$intact_id
all_intact_gr$score <- 0
#export(all_intact_gr, "all_intact_LTR.bed", format = "BED")

# Run on Linux:
# bedtools getfasta -fi GCA_037325665.1_US_GU_Lsax_2.0_genomic.fna -bed all_intact_LTR.bed -s -nameOnly > all_intact_LTR.fa
# makeblastdb -in TF2_clade_consensus.fa -dbtype nucl -out TF2_cons
# blastn -query all_intact_LTR.fa -db TF2_cons -evalue 1e-10 -max_target_seqs 10 -max_hsps 100 -num_threads 8 \
#        -outfmt "6 qseqid sseqid pident length qstart qend sstart send qlen slen bitscore" > all_intact_vs_cons.tsv

blast1 <- read.table("all_intact_vs_cons.tsv", col.names = c("qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "qlen", "slen", "bitscore"))
blast1$qseqid <- sub("\\(.*\\)$", "", blast1$qseqid)

# Hits of one element on the same consensus (e.g. one for each LTR) are merged before the coverage is calculated;
# strand = orientation of the element relative to the consensus.
merged_length <- function(s, e) sum(width(IRanges::reduce(IRanges(pmin(s, e), pmax(s, e)))))
blast1_summary <- blast1 %>% dplyr::group_by(qseqid, sseqid) %>%
  dplyr::summarise(aln_len = sum(length), identity = sum(pident * length) / sum(length), cov_intact = merged_length(qstart, qend) / dplyr::first(qlen) * 100,
                   strand = ifelse(sum(length[send > sstart]) >= sum(length) / 2, "+", "-"), bitscore = sum(bitscore), .groups = "drop")

# Copies of each clade: identity >= 80%, alignment >= 80 bp and >= 80% of the element covered (clade A: 138 copies, clade B: 64).
hits <- blast1_summary %>% dplyr::filter(identity >= 80, aln_len >= 80, cov_intact >= 80)
table(hits$sseqid)

# Genomic strand on which each copy has the same orientation as its consensus (strand "?" treated as +).
members <- hits %>% dplyr::left_join(repeat_region %>% dplyr::transmute(qseqid = intact_id, chr, start, end, strand_extract = ifelse(strand %in% c("+", "-"), strand, "+")), by = "qseqid") %>%
  dplyr::mutate(strand_cons = ifelse(strand_extract == strand, "+", "-")) %>% dplyr::arrange(chr, start)
member_A <- members %>% dplyr::filter(sseqid == "cladeA_cons")
member_B <- members %>% dplyr::filter(sseqid == "cladeB_cons")
#export(GRanges(member_A$chr, IRanges(member_A$start, member_A$end), strand = member_A$strand_cons, name = member_A$qseqid, score = 0), "cladeA_intact.bed", format = "BED")
#export(GRanges(member_B$chr, IRanges(member_B$start, member_B$end), strand = member_B$strand_cons, name = member_B$qseqid, score = 0), "cladeB_intact.bed", format = "BED")

# Sequences of all copies, in the orientation of the consensus.
all_seq <- readDNAStringSet("all_intact_LTR.fa")
names(all_seq) <- sub("\\(.*\\)$", "", names(all_seq))
copies_A <- all_seq[member_A$qseqid]
copies_B <- all_seq[member_B$qseqid]
copies_A[member_A$strand == "-"] <- reverseComplement(copies_A[member_A$strand == "-"])
copies_B[member_B$strand == "-"] <- reverseComplement(copies_B[member_B$strand == "-"])

# Final consensus from all copies: columns covered by fewer than 30% of the copies are dropped (gaps are not counted as a base),
# and the majority base must reach 50%, otherwise N.
make_consensus <- function(aln, min_coverage = 0.3, base_threshold = 0.5) {
  mat <- as.matrix(aln)
  mat <- mat[, colMeans(mat != "-") >= min_coverage, drop = FALSE]
  cons <- apply(mat, 2, function(x) {
    x <- x[x %in% c("A", "C", "G", "T")]
    if (length(x) == 0) return("N")
    freq <- table(x) / length(x)
    if (max(freq) >= base_threshold) names(freq)[which.max(freq)] else "N"
  })
  paste0(cons, collapse = "")
}
consensus2 <- DNAStringSet(c(cladeA_cons = make_consensus(AlignSeqs(copies_A)), cladeB_cons = make_consensus(AlignSeqs(copies_B))))
width(consensus2)
#writeXStringSet(consensus2, "consensus_cladeA_cladeB.fa")

# Run on Linux with the final consensus sequences:
# - TEsorter consensus_cladeA_cladeB.fa -db gydb (and -db rexdb-metazoa); GyDB domain positions plus the RExDB clade A GAG -> domain.txt
# - bedtools getfasta -s -nameOnly with cladeA_intact.bed / cladeB_intact.bed -> cladeA_intact.fa / cladeB_intact.fa
# - BLASTN of the copies of each clade against the consensus -> intactLTR_vs_consensus.blastn.tsv
# - mafft --add cladeA_intact.fa --keeplength cladeA_cons.fa > cladeA_to_cons.aln.fa (same for clade B)
# - RNA-seq reads and pool-seq reads mapped to the two consensus sequences plus two housekeeping genes (g15582.t1, g13887.t1):
#   per-base pool-seq depth -> all_depth.txt; RNA-seq read counts in the RT domain -> tf2_RT_counts.txt; mapped fragments -> library_size.txt


############################################################
#Part 4: structure of the two consensus sequences (LTRs and protein domains)

clade_order <- c("cladeA_cons", "cladeB_cons")
cons <- readDNAStringSet("consensus_cladeA_cladeB.fa")
names(cons) <- sub("[#[:space:]].*$", "", names(cons))
cons <- cons[clade_order]
cons_len <- setNames(width(cons), names(cons))

# 5′ and 3′ LTR of every copy (LTR_retriever): on the plus strand the left LTR is the 5′ LTR, on the minus strand the right LTR.
# Positions within the copy sequence, which is in the orientation of the consensus.
ltr_copy <- members %>% dplyr::select(qseqid, sseqid, strand_cons, start, end) %>% dplyr::left_join(repeat_region %>% dplyr::select(qseqid = intact_id, parent), by = "qseqid") %>%
  dplyr::left_join(ltr_region, by = "parent") %>%
  dplyr::mutate(label = dplyr::if_else((strand_cons == "+" & grepl("^lLTR_", ltr_id)) | (strand_cons == "-" & grepl("^rLTR_", ltr_id)), "5′ LTR", "3′ LTR"),
                LTR_length = ltr_end - ltr_start + 1, relative_start = dplyr::if_else(strand_cons == "+", ltr_start - start + 1, end - ltr_end + 1),
                relative_end = dplyr::if_else(strand_cons == "+", ltr_end - start + 1, end - ltr_start + 1))
ltr5_copy_A <- ltr_copy %>% dplyr::filter(sseqid == "cladeA_cons", label == "5′ LTR")
ltr3_copy_A <- ltr_copy %>% dplyr::filter(sseqid == "cladeA_cons", label == "3′ LTR")
ltr5_copy_B <- ltr_copy %>% dplyr::filter(sseqid == "cladeB_cons", label == "5′ LTR")
ltr3_copy_B <- ltr_copy %>% dplyr::filter(sseqid == "cladeB_cons", label == "3′ LTR")

# Position of the LTRs on the consensus: the LTR of each copy is aligned (whole LTR, global-local alignment) to the first half
# of its consensus (5′ LTR) or to the second half (3′ LTR); the LTR on the consensus is the start-end pair found in most copies.
half_A <- floor(cons_len[["cladeA_cons"]] / 2)
half_B <- floor(cons_len[["cladeB_cons"]] / 2)
aln5_A <- pairwiseAlignment(subseq(copies_A[ltr5_copy_A$qseqid], ltr5_copy_A$relative_start, ltr5_copy_A$relative_end), subseq(cons[["cladeA_cons"]], 1, half_A), type = "global-local")
aln3_A <- pairwiseAlignment(subseq(copies_A[ltr3_copy_A$qseqid], ltr3_copy_A$relative_start, ltr3_copy_A$relative_end), subseq(cons[["cladeA_cons"]], half_A + 1, cons_len[["cladeA_cons"]]), type = "global-local")
aln5_B <- pairwiseAlignment(subseq(copies_B[ltr5_copy_B$qseqid], ltr5_copy_B$relative_start, ltr5_copy_B$relative_end), subseq(cons[["cladeB_cons"]], 1, half_B), type = "global-local")
aln3_B <- pairwiseAlignment(subseq(copies_B[ltr3_copy_B$qseqid], ltr3_copy_B$relative_start, ltr3_copy_B$relative_end), subseq(cons[["cladeB_cons"]], half_B + 1, cons_len[["cladeB_cons"]]), type = "global-local")
pos5_A <- sort(table(paste(start(subject(aln5_A)), end(subject(aln5_A)))), decreasing = TRUE)
pos3_A <- sort(table(paste(start(subject(aln3_A)) + half_A, end(subject(aln3_A)) + half_A)), decreasing = TRUE)
pos5_B <- sort(table(paste(start(subject(aln5_B)), end(subject(aln5_B)))), decreasing = TRUE)
pos3_B <- sort(table(paste(start(subject(aln3_B)) + half_B, end(subject(aln3_B)) + half_B)), decreasing = TRUE)
# Most frequent start-end pairs and the number of copies supporting them.
head(pos5_A); head(pos3_A); head(pos5_B); head(pos3_B)

ltr_cons <- data.frame(sseqid = rep(clade_order, each = 2), label = rep(c("5′ LTR", "3′ LTR"), 2),
                       position = c(names(pos5_A)[1], names(pos3_A)[1], names(pos5_B)[1], names(pos3_B)[1])) %>%
  tidyr::separate(position, into = c("start", "end"), sep = " ", convert = TRUE) %>% dplyr::mutate(LTR_length = end - start + 1, feature = "LTR")
ltr_cons

# 5′ LTR of the two consensus sequences (used for Figure 5B).
ltr5_consensus <- DNAStringSet(c(cladeA_5LTR = as.character(subseq(cons[["cladeA_cons"]], ltr_cons$start[1], ltr_cons$end[1])),
                                 cladeB_5LTR = as.character(subseq(cons[["cladeB_cons"]], ltr_cons$start[3], ltr_cons$end[3]))))

# LTRs and protein domains (TEsorter) of the two consensus sequences, used in Figure 5A and Figure S14.
domain <- read.table("domain.txt", header = FALSE, col.names = c("sseqid", "start", "end", "feature"))
features <- dplyr::bind_rows(ltr_cons %>% dplyr::select(sseqid, start, end, feature, label), domain %>% dplyr::mutate(label = feature)) %>%
  dplyr::mutate(mid = (start + end) / 2, sseqid = factor(sseqid, levels = clade_order))


############################################################
#Part 4B: similarity between the two consensus sequences
# BLASTN between the two consensus sequences (E-value <= 1e-5); test of their best local alignment score against 1,000
# shuffled clade B sequences with the same dinucleotide frequencies; amino acid identity of their RT domains (GyDB, global alignment).

compare_tf2_consensuses <- function(
    consensus_file = "consensus_cladeA_cladeB.fa",
    gydb_protein_file = "consensus_cladeA_cladeB.fa.gydb.dom.faa",
    n_permutations = 1000L, seed = 1L,
    blastn_exe = "C:/Users/hanbq/Documents/Codex/2026-10-22/xia/work/tools/blast/bin/blastn.exe") {
  
  seqs <- Biostrings::readDNAStringSet(consensus_file)
  names(seqs) <- sub("[#[:space:]].*$", "", names(seqs))
  ids <- c("cladeA_cons", "cladeB_cons")
  seqs <- seqs[ids]
  seq_lengths <- setNames(BiocGenerics::width(seqs), ids)
  
  tmp <- tempfile("tf2_consensus_AB_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  query_file <- file.path(tmp, "cladeA_cons.fa")
  subject_file <- file.path(tmp, "cladeB_cons.fa")
  blast_segments_file <- file.path(tmp, "cladeA_vs_cladeB.blastn.tsv")
  Biostrings::writeXStringSet(seqs[1], query_file)
  Biostrings::writeXStringSet(seqs[2], subject_file)
  fields <- c("qseqid", "sseqid", "pident", "length", "qstart", "qend",
              "sstart", "send", "evalue", "bitscore")
  args <- c("-task", "blastn", "-query", shQuote(query_file),
            "-subject", shQuote(subject_file), "-evalue", "1e-5",
            "-outfmt", shQuote(paste("6", paste(fields, collapse = " "))),
            "-out", shQuote(blast_segments_file))
  system2(blastn_exe, args = args, stdout = TRUE, stderr = TRUE)
  # The BLAST output is empty when there is no hit.
  blast_segments <- data.frame(qseqid = character(), sseqid = character(), pident = numeric(),
                               length = integer(), qstart = integer(), qend = integer(),
                               sstart = integer(), send = integer(), evalue = numeric(), bitscore = numeric())
  if (file.info(blast_segments_file)$size > 0)
    blast_segments <- read.table(blast_segments_file, col.names = fields, stringsAsFactors = FALSE)
  union_bp <- function(start, end) {
    if (!length(start)) return(0L)
    sum(BiocGenerics::width(IRanges::reduce(IRanges::IRanges(pmin(start, end), pmax(start, end)))))
  }
  cov_bp <- c(union_bp(blast_segments$qstart, blast_segments$qend), union_bp(blast_segments$sstart, blast_segments$send))
  blast_summary <- data.frame(
    consensus = ids, length_bp = unname(seq_lengths), n_blast_segments = nrow(blast_segments),
    covered_bp = cov_bp, coverage_pct = 100 * cov_bp / unname(seq_lengths))
  
  mat <- Biostrings::nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = FALSE)
  local_score <- function(x, y) max(
    Biostrings::pairwiseAlignment(x, y, type = "local", substitutionMatrix = mat,
                                  gapOpening = 5, gapExtension = 2, scoreOnly = TRUE),
    Biostrings::pairwiseAlignment(x, Biostrings::reverseComplement(y), type = "local",
                                  substitutionMatrix = mat, gapOpening = 5, gapExtension = 2, scoreOnly = TRUE))
  observed <- local_score(seqs[[1]], seqs[[2]])
  set.seed(seed)
  null_scores <- replicate(n_permutations,
                           local_score(seqs[[1]], universalmotif::shuffle_sequences(seqs[2], k = 2)[[1]]))
  nucleotide_summary <- data.frame(
    cladeA_bp = BiocGenerics::width(seqs)[1], cladeB_bp = BiocGenerics::width(seqs)[2],
    observed_score = observed, n_permutations = n_permutations,
    n_null_ge_observed = sum(null_scores >= observed),
    empirical_p = mean(null_scores >= observed),
    null_median = median(null_scores), null_min = min(null_scores), null_max = max(null_scores))
  
  domains <- Biostrings::readAAStringSet(gydb_protein_file)
  rt_index <- vapply(ids, function(id) {
    idx <- which(startsWith(names(domains), paste0(id, "|")) &
                   grepl(";gene=RT(;|$)", names(domains)))
    idx
  }, integer(1))
  rt <- domains[rt_index]
  names(rt) <- ids
  rt_alignment <- Biostrings::pairwiseAlignment(
    rt[[1]], rt[[2]], type = "global", substitutionMatrix = "BLOSUM62",
    gapOpening = 10, gapExtension = 0.5)
  # Identity = identical positions / aligned positions including gaps (PID1).
  rt_summary <- data.frame(cladeA_RT_aa = BiocGenerics::width(rt)[1], cladeB_RT_aa = BiocGenerics::width(rt)[2],
                           identity_pct = Biostrings::pid(rt_alignment, type = "PID1"))
  list(blast_segments = blast_segments, blast_summary = blast_summary,
       nucleotide_summary = nucleotide_summary, null_scores = null_scores,
       rt_proteins = rt, rt_alignment = rt_alignment, rt_summary = rt_summary)
}

tf2_consensus_comparison <- compare_tf2_consensuses()
consensus_AB_blast_segments <- tf2_consensus_comparison$blast_segments
consensus_AB_coverage <- tf2_consensus_comparison$blast_summary
consensus_AB_shuffle_test <- tf2_consensus_comparison$nucleotide_summary
consensus_AB_null_scores <- tf2_consensus_comparison$null_scores
consensus_RT_identity <- tf2_consensus_comparison$rt_summary
consensus_RT_alignment <- tf2_consensus_comparison$rt_alignment
print(consensus_AB_blast_segments)
print(consensus_AB_coverage)
print(consensus_AB_shuffle_test)
print(consensus_RT_identity)
# Results: no BLASTN hit; shuffle test p = 0.381; RT amino acid identity 30.5%.


############################################################
#Part 4C: GAG, INT, RT and RNase H in each copy of clade A and clade B
# Each domain of the consensus (domain.txt, as in Figure 5A) is searched by BLASTN in the copies of its own clade
# (E-value <= 1e-10, identity >= 80%, alignment length >= 80 bp); domain coverage = aligned length / domain length.

check_tf2_copy_domains <- function(
    copy_members, consensus_file = "consensus_cladeA_cladeB.fa",
    domain_file = "domain.txt", copy_file = "all_intact_LTR.fa",
    blastn_exe = "C:/Users/hanbq/Documents/Codex/2026-10-22/xia/work/tools/blast/bin/blastn.exe") {
  clades <- c("cladeA_cons", "cladeB_cons")
  features <- c("GAG", "INT", "RT", "RNaseH")
  cm <- unique(as.data.frame(copy_members[, c("qseqid", "sseqid")]))
  names(cm) <- c("copy_id", "clade")
  cm$copy_id <- sub("\\(.*\\)$", "", cm$copy_id)
  cons <- Biostrings::readDNAStringSet(consensus_file)
  names(cons) <- sub("[#[:space:]].*$", "", names(cons))
  copies <- Biostrings::readDNAStringSet(copy_file)
  names(copies) <- sub("\\(.*\\)$", "", sub("[[:space:]].*$", "", names(copies)))
  dom <- read.table(domain_file, col.names = c("clade", "start", "end", "domain"),
                    stringsAsFactors = FALSE)
  dom <- dom[dom$clade %in% clades & dom$domain %in% features, ]
  dom$query_id <- paste(dom$clade, dom$domain, sep = "__")
  dom$domain_bp <- dom$end - dom$start + 1L
  queries <- Biostrings::DNAStringSet(setNames(vapply(seq_len(nrow(dom)), function(i)
    as.character(Biostrings::subseq(cons[[dom$clade[i]]], start = dom$start[i], end = dom$end[i])),
    character(1)), dom$query_id))
  
  tmp <- tempfile("tf2_copy_domains_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  query_file <- file.path(tmp, "consensus_domains.fa")
  subject_file <- file.path(tmp, "assigned_copies.fa")
  blast_segments_file <- file.path(tmp, "domains_vs_copies.tsv")
  Biostrings::writeXStringSet(queries, query_file)
  Biostrings::writeXStringSet(copies[cm$copy_id], subject_file)
  fields <- c("qseqid", "sseqid", "pident", "length", "qstart", "qend",
              "sstart", "send", "evalue", "bitscore", "qlen", "slen")
  system2(blastn_exe, args = c("-task", "blastn", "-query", shQuote(query_file),
                               "-subject", shQuote(subject_file), "-evalue", "1e-10",
                               "-outfmt", shQuote(paste("6", paste(fields, collapse = " "))),
                               "-out", shQuote(blast_segments_file)), stdout = TRUE, stderr = TRUE)
  blast_segments <- read.table(blast_segments_file, col.names = fields, stringsAsFactors = FALSE)
  blast_segments <- blast_segments[blast_segments$evalue <= 1e-10 & blast_segments$pident >= 80 & blast_segments$length >= 80, ]
  blast_segments <- blast_segments[dom$clade[match(blast_segments$qseqid, dom$query_id)] == cm$clade[match(blast_segments$sseqid, cm$copy_id)], ]
  by_pair <- split(blast_segments, paste(blast_segments$qseqid, blast_segments$sseqid, sep = "\t"))
  per_copy <- merge(cm, dom, by = "clade", sort = FALSE)
  per_copy$aligned_domain_bp <- vapply(seq_len(nrow(per_copy)), function(i) {
    hh <- by_pair[[paste(per_copy$query_id[i], per_copy$copy_id[i], sep = "\t")]]
    if (is.null(hh)) return(0L)
    as.integer(sum(BiocGenerics::width(IRanges::reduce(IRanges::IRanges(
      pmin(hh$qstart, hh$qend), pmax(hh$qstart, hh$qend))))))
  }, integer(1))
  per_copy$coverage_pct <- 100 * per_copy$aligned_domain_bp / per_copy$domain_bp
  per_copy$detected <- per_copy$aligned_domain_bp > 0L
  per_copy$coverage_ge80 <- per_copy$coverage_pct >= 80
  per_copy$status <- ifelse(!per_copy$detected, "No qualified match",
                            ifelse(per_copy$coverage_ge80, ">=80% domain coverage", "Partial: <80% domain coverage"))
  per_copy <- per_copy[order(match(per_copy$clade, clades), per_copy$copy_id,
                             match(per_copy$domain, features)), ]
  rownames(per_copy) <- NULL
  summary <- per_copy %>% dplyr::group_by(clade, domain) %>%
    dplyr::summarise(n_copies = dplyr::n(), n_detected = sum(detected),
                     n_coverage_ge80 = sum(coverage_ge80), n_partial_lt80 = sum(detected & !coverage_ge80),
                     n_no_match = sum(!detected), detected_pct = 100 * mean(detected),
                     coverage_ge80_pct = 100 * mean(coverage_ge80), .groups = "drop") %>%
    dplyr::arrange(match(clade, clades), match(domain, features))
  list(per_copy = per_copy, summary = summary, qualified_blast_segments = blast_segments,
       annotation = dom)
}

tf2_copy_domains <- check_tf2_copy_domains(members)
copy_domain_retention <- tf2_copy_domains$per_copy  # 202 copies x 4 domains = 808 rows.
copy_domain_summary <- tf2_copy_domains$summary   # Eight clade/domain combinations.
copy_domain_exceptions <- copy_domain_retention[!copy_domain_retention$coverage_ge80, ]
copy_domain_completeness <- copy_domain_retention %>% dplyr::group_by(clade, copy_id) %>%
  dplyr::summarise(n_domains_detected = sum(detected), n_domains_coverage_ge80 = sum(coverage_ge80),
                   all_four_detected = all(detected), all_four_coverage_ge80 = all(coverage_ge80), .groups = "drop")
print(copy_domain_summary)


############################################################
#Part 5: HSF1 binding sites in the 5′ LTRs (Figure S16)

# 5′ LTR of every copy, cut from the copy sequence (already in the orientation of the consensus).
ltr5 <- members %>% dplyr::left_join(repeat_region %>% dplyr::select(qseqid = intact_id, parent), by = "qseqid") %>% dplyr::left_join(ltr_region, by = "parent") %>%
  dplyr::filter((strand_cons == "+" & grepl("^lLTR_", ltr_id)) | (strand_cons == "-" & grepl("^rLTR_", ltr_id))) %>%
  dplyr::mutate(ltr_length = ltr_end - ltr_start + 1, relative_start = dplyr::if_else(strand_cons == "+", ltr_start - start + 1, end - ltr_end + 1),
                relative_end = dplyr::if_else(strand_cons == "+", ltr_end - start + 1, end - ltr_start + 1))
ltr5_A <- ltr5 %>% dplyr::filter(sseqid == "cladeA_cons") %>% dplyr::arrange(chr, start) %>% dplyr::mutate(fasta_id = sprintf("cladeA_5LTR_%03d", dplyr::row_number()))
ltr5_B <- ltr5 %>% dplyr::filter(sseqid == "cladeB_cons") %>% dplyr::arrange(chr, start) %>% dplyr::mutate(fasta_id = sprintf("cladeB_5LTR_%03d", dplyr::row_number()))
seq5_A <- subseq(copies_A[ltr5_A$qseqid], start = ltr5_A$relative_start, end = ltr5_A$relative_end)
seq5_B <- subseq(copies_B[ltr5_B$qseqid], start = ltr5_B$relative_start, end = ltr5_B$relative_end)
names(seq5_A) <- ltr5_A$fasta_id
names(seq5_B) <- ltr5_B$fasta_id
#writeXStringSet(seq5_A, "cladeA_5LTR.fa")
#writeXStringSet(seq5_B, "cladeB_5LTR.fa")

# Run on Linux: FIMO with the HSF1 motif (JASPAR MA0486.2), p < 1e-4.
# fimo --o fimo_cladeA_5LTR_MA04862 --thresh 1e-4 MA0486.2.meme cladeA_5LTR.fa -> cladeA_5LTR_MA04862_fimo.tsv
# fimo --o fimo_cladeB_5LTR_MA04862 --thresh 1e-4 MA0486.2.meme cladeB_5LTR.fa -> cladeB_5LTR_MA04862_fimo.tsv
# Positive control: 500 bp upstream of the start codon of the eight Hsp70 genes -> 8HSP70_8ATG_candidates_upstream500_MA04862_fimo.tsv
# Background: oriented 5′ LTRs of other intact elements -> otherLTR_MA04862_fimo.tsv
fimo_A <- read.delim("cladeA_5LTR_MA04862_fimo.tsv", comment.char = "#", quote = "")
fimo_B <- read.delim("cladeB_5LTR_MA04862_fimo.tsv", comment.char = "#", quote = "")
fimo_HSP70 <- read.delim("8HSP70_8ATG_candidates_upstream500_MA04862_fimo.tsv", comment.char = "#", quote = "")
seq_HSP70 <- readDNAStringSet("8HSP70_8ATG_candidates_upstream500.fa")
fimo_other <- read.delim("otherLTR_MA04862_fimo.tsv", comment.char = "#", quote = "")
seq_other <- readDNAStringSet("otherLTR_5LTR.fa")

# Keep only the 5′ LTRs of the other intact Ty3/gypsy elements, the same elements as in the insertion-time analysis (Part 7).
other_ltr_background_bed <- read.delim("otherLTR_5LTR.bed", header = FALSE, quote = "",
                                       col.names = c("chr", "start0", "end", "fasta_id", "score", "strand"))
other_ltr_background_map <- other_ltr_background_bed %>%
  dplyr::mutate(sequence_name = paste0(fasta_id, "(", strand, ")"),
                ltr_start = start0 + 1L, ltr_end = end) %>%
  dplyr::left_join(gff %>% dplyr::filter(type == "long_terminal_repeat") %>%
                     dplyr::transmute(chr, ltr_start = start, ltr_end = end,
                                      parent = sub(".*Parent=([^;]+).*", "\\1", attribute)),
                   by = c("chr", "ltr_start", "ltr_end")) %>%
  dplyr::left_join(repeat_region %>% dplyr::transmute(parent, intact_id,
                                                      SuperFamily = sub(".*classification=LTR/([^;]+).*", "\\1", attribute)), by = "parent")
other_ltr_gypsy_map <- other_ltr_background_map %>% dplyr::filter(SuperFamily == "Gypsy")
seq_other <- seq_other[names(seq_other) %in% other_ltr_gypsy_map$sequence_name]
fimo_other <- fimo_other %>% dplyr::filter(sequence_name %in% names(seq_other))

# HSF1 motif used in all analyses.
hsf1_motif_id <- "MA0486.2"

# One line per sequence (copies labelled with their genomic position).
line_A <- data.frame(sequence_name = ltr5_A$fasta_id, LTR_length = ltr5_A$ltr_length, display_id = ltr5_A$qseqid) %>% dplyr::arrange(sequence_name) %>% dplyr::mutate(ypos = rev(seq_len(dplyr::n())))
line_B <- data.frame(sequence_name = ltr5_B$fasta_id, LTR_length = ltr5_B$ltr_length, display_id = ltr5_B$qseqid) %>% dplyr::arrange(sequence_name) %>% dplyr::mutate(ypos = rev(seq_len(dplyr::n())))
line_HSP70 <- data.frame(sequence_name = names(seq_HSP70), upstream_length = width(seq_HSP70)) %>% dplyr::arrange(sequence_name) %>%
  dplyr::mutate(ypos = rev(seq_len(dplyr::n())), display_id = sub("\\([+-]\\)$", "", sequence_name))

# Motif hits with p < 1e-4; overlapping hits (same site found on both strands or shifted by a few bases) count as one site,
# represented by the hit with the highest score.
point_A <- fimo_A %>% dplyr::filter(p.value < 1e-4) %>% dplyr::group_by(sequence_name, motif_id) %>% dplyr::arrange(start, stop, .by_group = TRUE) %>%
  dplyr::mutate(hit_group = cumsum(start > dplyr::lag(cummax(stop), default = 0L))) %>% dplyr::group_by(sequence_name, motif_id, hit_group) %>% dplyr::slice_max(score, n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>%
  dplyr::left_join(line_A %>% dplyr::select(sequence_name, ypos), by = "sequence_name") %>% dplyr::mutate(mid = (start + stop) / 2)
point_B <- fimo_B %>% dplyr::filter(p.value < 1e-4) %>% dplyr::group_by(sequence_name, motif_id) %>% dplyr::arrange(start, stop, .by_group = TRUE) %>%
  dplyr::mutate(hit_group = cumsum(start > dplyr::lag(cummax(stop), default = 0L))) %>% dplyr::group_by(sequence_name, motif_id, hit_group) %>% dplyr::slice_max(score, n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>%
  dplyr::left_join(line_B %>% dplyr::select(sequence_name, ypos), by = "sequence_name") %>% dplyr::mutate(mid = (start + stop) / 2)
# Hsp70 positions relative to the start codon (last upstream base = -1).
point_HSP70 <- fimo_HSP70 %>% dplyr::filter(p.value < 1e-4) %>% dplyr::group_by(sequence_name, motif_id) %>% dplyr::arrange(start, stop, .by_group = TRUE) %>%
  dplyr::mutate(hit_group = cumsum(start > dplyr::lag(cummax(stop), default = 0L))) %>% dplyr::group_by(sequence_name, motif_id, hit_group) %>% dplyr::slice_max(score, n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>%
  dplyr::left_join(line_HSP70 %>% dplyr::select(sequence_name, ypos, upstream_length), by = "sequence_name") %>% dplyr::mutate(mid_relative = (start + stop) / 2 - upstream_length - 1)

# 5′ LTRs of the other intact LTR elements: hits merged in the same way.
point_other <- fimo_other %>% dplyr::filter(p.value < 1e-4) %>% dplyr::group_by(sequence_name, motif_id) %>% dplyr::arrange(start, stop, .by_group = TRUE) %>%
  dplyr::mutate(hit_group = cumsum(start > dplyr::lag(cummax(stop), default = 0L))) %>% dplyr::group_by(sequence_name, motif_id, hit_group) %>% dplyr::slice_max(score, n = 1, with_ties = FALSE) %>% dplyr::ungroup()

# Number and percentage of sequences with at least one HSF1 site (p < 1e-4).
hsf1_summary <- data.frame(group = c("Hsp70 upstream", "Clade A 5′ LTR", "Clade B 5′ LTR", "Other intact Ty3/gypsy 5′ LTR"),
                           n_sequences = c(length(seq_HSP70), nrow(ltr5_A), nrow(ltr5_B), length(seq_other)),
                           n_with_site = c(dplyr::n_distinct(point_HSP70$sequence_name), dplyr::n_distinct(point_A$sequence_name), dplyr::n_distinct(point_B$sequence_name),
                                           dplyr::n_distinct(point_other$sequence_name)))
hsf1_summary$percent_with_site <- round(100 * hsf1_summary$n_with_site / hsf1_summary$n_sequences, 1)

# Fisher's exact tests: 5′ LTRs with at least one HSF1 site in clade A or clade B vs the other Ty3/gypsy elements (BH over the two tests).
fisher_A <- fisher.test(matrix(c(hsf1_summary$n_with_site[2], hsf1_summary$n_sequences[2] - hsf1_summary$n_with_site[2],
                                 hsf1_summary$n_with_site[4], hsf1_summary$n_sequences[4] - hsf1_summary$n_with_site[4]), nrow = 2))
fisher_B <- fisher.test(matrix(c(hsf1_summary$n_with_site[3], hsf1_summary$n_sequences[3] - hsf1_summary$n_with_site[3],
                                 hsf1_summary$n_with_site[4], hsf1_summary$n_sequences[4] - hsf1_summary$n_with_site[4]), nrow = 2))
hsf1_summary$odds_ratio_vs_other <- c(NA, fisher_A$estimate, fisher_B$estimate, NA)
hsf1_summary$p_vs_other <- c(NA, fisher_A$p.value, fisher_B$p.value, NA)
hsf1_summary$q_vs_other <- c(NA, p.adjust(c(fisher_A$p.value, fisher_B$p.value), "BH"), NA)
hsf1_summary

# The most frequent motif sequence in each clade (in consensus orientation) is located in the 5′ LTR of the consensus;
# this position is marked in Figure 5A.
dominant_A <- point_A %>% dplyr::mutate(cons_sequence = toupper(dplyr::if_else(strand == "+", matched_sequence, as.character(reverseComplement(DNAStringSet(matched_sequence)))))) %>%
  dplyr::count(cons_sequence, sort = TRUE) %>% dplyr::slice_head(n = 1)
dominant_B <- point_B %>% dplyr::mutate(cons_sequence = toupper(dplyr::if_else(strand == "+", matched_sequence, as.character(reverseComplement(DNAStringSet(matched_sequence)))))) %>%
  dplyr::count(cons_sequence, sort = TRUE) %>% dplyr::slice_head(n = 1)
hsf_A <- matchPattern(dominant_A$cons_sequence, cons[["cladeA_cons"]])
hsf_B <- matchPattern(dominant_B$cons_sequence, cons[["cladeB_cons"]])
hsf_A <- hsf_A[start(hsf_A) >= ltr_cons$start[1] & end(hsf_A) <= ltr_cons$end[1]]
hsf_B <- hsf_B[start(hsf_B) >= ltr_cons$start[3] & end(hsf_B) <= ltr_cons$end[3]]
# Centre of the HSF1 site on each consensus, marked in Figure 5A.
hsf_sites <- data.frame(sseqid = factor(clade_order, levels = clade_order), mid = c((start(hsf_A) + end(hsf_A)) / 2, (start(hsf_B) + end(hsf_B)) / 2))
hsf_sites

# Figure S16: one row per sequence, HSF1 sites as triangles and candidate TATA boxes as boxes.
# A: 5′ LTRs of the other intact Ty3/gypsy elements (longest at the top); B: 500 bp upstream of the start codon of the Hsp70 genes;
# C, D: 5′ LTRs of the clade A and clade B copies.
line_other <- data.frame(sequence_name = names(seq_other), LTR_length = width(seq_other)) %>% dplyr::arrange(LTR_length, sequence_name) %>% dplyr::mutate(ypos = seq_len(dplyr::n()))
plot_other <- point_other %>% dplyr::left_join(line_other %>% dplyr::select(sequence_name, ypos), by = "sequence_name") %>% dplyr::mutate(mid = (start + stop) / 2)
length_HSP70 <- max(line_HSP70$upstream_length)
length_other <- max(line_other$LTR_length)
length_A <- max(line_A$LTR_length)
length_B <- max(line_B$LTR_length)
c(Hsp70 = length_HSP70, other = length_other, clade_A = length_A, clade_B = length_B)

# Candidate TATA boxes (TATAWAW) in all sequences; overlapping matches (e.g. TATATAT) count as one. Hsp70 positions are relative to the start codon.
tata_HSP70 <- as.data.frame(vmatchPattern("TATAWAW", seq_HSP70, fixed = "subject"))
tata_HSP70$sequence_name <- names(seq_HSP70)[tata_HSP70$group]
tata_HSP70 <- tata_HSP70 %>% dplyr::group_by(sequence_name) %>% dplyr::arrange(start, .by_group = TRUE) %>%
  dplyr::mutate(hit_group = cumsum(start > dplyr::lag(cummax(end), default = 0L))) %>% dplyr::group_by(sequence_name, hit_group) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>% dplyr::left_join(line_HSP70 %>% dplyr::select(sequence_name, ypos, upstream_length), by = "sequence_name") %>%
  dplyr::mutate(start_relative = start - upstream_length - 1, end_relative = end - upstream_length - 1)
tata_copy_A <- as.data.frame(vmatchPattern("TATAWAW", seq5_A, fixed = "subject"))
tata_copy_A$sequence_name <- names(seq5_A)[tata_copy_A$group]
tata_copy_A <- tata_copy_A %>% dplyr::group_by(sequence_name) %>% dplyr::arrange(start, .by_group = TRUE) %>%
  dplyr::mutate(hit_group = cumsum(start > dplyr::lag(cummax(end), default = 0L))) %>% dplyr::group_by(sequence_name, hit_group) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>% dplyr::left_join(line_A %>% dplyr::select(sequence_name, ypos), by = "sequence_name")
tata_copy_B <- as.data.frame(vmatchPattern("TATAWAW", seq5_B, fixed = "subject"))
tata_copy_B$sequence_name <- names(seq5_B)[tata_copy_B$group]
tata_copy_B <- tata_copy_B %>% dplyr::group_by(sequence_name) %>% dplyr::arrange(start, .by_group = TRUE) %>%
  dplyr::mutate(hit_group = cumsum(start > dplyr::lag(cummax(end), default = 0L))) %>% dplyr::group_by(sequence_name, hit_group) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>% dplyr::left_join(line_B %>% dplyr::select(sequence_name, ypos), by = "sequence_name")
tata_other <- as.data.frame(vmatchPattern("TATAWAW", seq_other, fixed = "subject"))
tata_other$sequence_name <- names(seq_other)[tata_other$group]
tata_other <- tata_other %>% dplyr::group_by(sequence_name) %>% dplyr::arrange(start, .by_group = TRUE) %>%
  dplyr::mutate(hit_group = cumsum(start > dplyr::lag(cummax(end), default = 0L))) %>% dplyr::group_by(sequence_name, hit_group) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>% dplyr::left_join(line_other %>% dplyr::select(sequence_name, ypos), by = "sequence_name")
# Number of sequences with at least one candidate TATA box.
c(other = dplyr::n_distinct(tata_other$sequence_name), Hsp70 = dplyr::n_distinct(tata_HSP70$sequence_name), clade_A = dplyr::n_distinct(tata_copy_A$sequence_name),
  clade_B = dplyr::n_distinct(tata_copy_B$sequence_name))

# Common theme of Figure S16
theme_fimo <- theme_minimal(base_size = 12, base_family = "Arial") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 13, margin = margin(t = 7)),
        axis.text.y = element_blank(), plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(t = 0, b = 3)), plot.tag = element_text(size = 17, face = "bold"),
        plot.tag.position = "topleft", plot.margin = margin(2, 10, 2, 8))


# Figure S16A: boxes around the seven other elements similar to clade A or B (not assigned to the clades by the 80% coverage rule).
s16_targets <- data.frame(
  sequence_name = c("other_5LTR_0034(-)", "other_5LTR_0052(-)",
                    "other_5LTR_0220(-)", "other_5LTR_0411(+)",
                    "other_5LTR_0419(+)", "other_5LTR_0438(-)",
                    "other_5LTR_0064(-)"),
  clade = c(rep("A", 6), "B"), stringsAsFactors = FALSE)
s16_marked <- dplyr::inner_join(s16_targets, line_other, by = "sequence_name") %>%
  dplyr::arrange(clade, ypos) %>% dplyr::group_by(clade) %>%
  dplyr::mutate(box_run = cumsum(ypos != dplyr::lag(ypos, default = dplyr::first(ypos) - 2L) + 1L)) %>%
  dplyr::ungroup()
s16_site_ranges <- dplyr::bind_rows(
  plot_other %>% dplyr::select(sequence_name, start, end = stop),
  tata_other %>% dplyr::select(sequence_name, start, end)) %>%
  dplyr::inner_join(s16_marked %>% dplyr::select(sequence_name, clade, box_run, ypos), by = "sequence_name")
s16_boxes <- s16_site_ranges %>% dplyr::group_by(clade, box_run) %>%
  dplyr::summarise(xmin = min(start) - 12, xmax = max(end) + 12,
                   ymin = min(ypos) - 0.62, ymax = max(ypos) + 0.62,
                   .groups = "drop") %>%
  dplyr::mutate(ymid = (ymin + ymax) / 2)
s16_labels <- s16_boxes %>% dplyr::group_by(clade) %>%
  dplyr::summarise(y = mean(ymid), .groups = "drop") %>%
  dplyr::mutate(x = 600, label = paste0("Clade ", clade, "-like (n = ", as.integer(table(s16_targets$clade)[clade]), ")"))
s16_leaders <- s16_boxes %>%
  dplyr::transmute(clade, x = xmax, xend = 565, y = ymid, yend = ymid)
s16_leaders <- dplyr::bind_rows(s16_leaders,
                                s16_boxes %>% dplyr::group_by(clade) %>%
                                  dplyr::summarise(x = 565, xend = 565, y = min(ymid), yend = max(ymid), .groups = "drop"),
                                s16_labels %>% dplyr::transmute(clade, x = 565, xend = 590, y, yend = y))

# Figure S16A: boxes around the four remaining other elements with an HSF1 site upstream of a TATA box.
s16_other_pair_ranges <- point_other %>%
  dplyr::transmute(sequence_name, hsf_start = start, hsf_stop = stop) %>%
  dplyr::inner_join(tata_other %>% dplyr::transmute(sequence_name, tata_start = start, tata_end = end),
                    by = "sequence_name", relationship = "many-to-many") %>%
  dplyr::filter(tata_start > hsf_stop) %>%
  dplyr::mutate(gap_bp = tata_start - hsf_stop - 1L)
s16_remaining_four_rows <- s16_other_pair_ranges %>%
  dplyr::anti_join(s16_targets %>% dplyr::select(sequence_name), by = "sequence_name") %>%
  dplyr::distinct(sequence_name, gap_bp) %>%
  dplyr::inner_join(line_other, by = "sequence_name") %>% dplyr::arrange(ypos) %>%
  dplyr::mutate(box_run = cumsum(ypos != dplyr::lag(ypos, default = dplyr::first(ypos) - 2L) + 1L))
s16_remaining_four_ranges <- s16_other_pair_ranges %>%
  dplyr::inner_join(s16_remaining_four_rows %>% dplyr::select(sequence_name, ypos, box_run), by = "sequence_name")
s16_remaining_four_boxes <- s16_remaining_four_ranges %>% dplyr::group_by(box_run) %>%
  dplyr::summarise(xmin = min(hsf_start) - 12, xmax = max(tata_end) + 12,
                   ymin = min(ypos) - 0.62, ymax = max(ypos) + 0.62, .groups = "drop") %>%
  dplyr::mutate(ymid = (ymin + ymax) / 2)
s16_remaining_four_labels <- data.frame(x = 900, y = round(0.60 * nrow(line_other)),
                                        label = "Remaining elements (n = 4)\nHSF1-TATA gap: 254 bp")
s16_remaining_four_leaders <- dplyr::bind_rows(
  s16_remaining_four_boxes %>% dplyr::transmute(x = xmax, xend = 930, y = ymid, yend = ymid),
  data.frame(x = 930, xend = 930, y = mean(s16_remaining_four_boxes$ymid), yend = s16_remaining_four_labels$y),
  data.frame(x = 930, xend = 910, y = s16_remaining_four_labels$y, yend = s16_remaining_four_labels$y))


pS16A <- ggplot() + geom_segment(data = line_other, aes(x = 1, xend = LTR_length, y = ypos, yend = ypos), colour = "grey80", linewidth = 0.2) +
  geom_segment(data = tata_other, aes(x = start, xend = end, y = ypos, yend = ypos), colour = "#315B8A", alpha = 0.65, linewidth = 0.85, lineend = "butt") +
  geom_point(data = plot_other, aes(x = mid, y = ypos), shape = 17, size = 1.0, colour = "firebrick") +
  geom_rect(data = s16_boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, colour = clade),
            fill = NA, linewidth = 0.6, show.legend = FALSE) +
  geom_segment(data = s16_leaders, aes(x = x, xend = xend, y = y, yend = yend, colour = clade),
               linewidth = 0.4, show.legend = FALSE) +
  geom_text(data = s16_labels, aes(x = x, y = y, label = label, colour = clade),
            hjust = 0, size = 12, size.unit = "pt", lineheight = 1.05, fontface = "bold", show.legend = FALSE) +
  geom_rect(data = s16_remaining_four_boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            colour = "#505050", fill = NA, linewidth = 0.6, linetype = "dashed") +
  geom_segment(data = s16_remaining_four_leaders, aes(x = x, xend = xend, y = y, yend = yend),
               colour = "#505050", linewidth = 0.4) +
  geom_text(data = s16_remaining_four_labels, aes(x = x, y = y, label = label),
            colour = "#505050", hjust = 1, size = 12, size.unit = "pt", lineheight = 1.05, fontface = "bold") +
  scale_colour_manual(values = c(A = "#D55E00", B = "#0072B2"), guide = "none") +
  scale_x_continuous(limits = c(0, length_other), breaks = seq(0, length_other, 200), expand = expansion(add = c(0, 5))) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  labs(x = "Position within 5′ LTR (bp)", y = NULL, title = paste0("Other intact Ty3/gypsy (n = ", nrow(line_other), ")"), tag = "A") + theme_fimo
pS16B <- ggplot() + geom_segment(data = line_HSP70, aes(x = -upstream_length, xend = -1, y = ypos, yend = ypos), colour = "grey80", linewidth = 0.5) +
  geom_segment(data = tata_HSP70, aes(x = start_relative, xend = end_relative, y = ypos, yend = ypos), colour = "#315B8A", alpha = 0.65, linewidth = 1.3, lineend = "butt") +
  geom_point(data = point_HSP70, aes(x = mid_relative, y = ypos), shape = 17, size = 1.65, colour = "firebrick") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  scale_x_continuous(limits = c(-length_HSP70, 0), breaks = seq(-length_HSP70, 0, 100), expand = expansion(add = c(0, 5))) +
  scale_y_continuous(breaks = line_HSP70$ypos, labels = line_HSP70$display_id, expand = expansion(add = 0.5)) +
  labs(x = "Position relative to start codon (bp)", y = NULL, title = paste0("Hsp70 (n = ", nrow(line_HSP70), ")"), tag = "B") + theme_fimo +
  theme(axis.text.y = element_text(size = 12, colour = "black", hjust = 0, margin = margin(r = 4)))
pS16C <- ggplot() + geom_segment(data = line_A, aes(x = 1, xend = LTR_length, y = ypos, yend = ypos), colour = "grey80", linewidth = 0.35) +
  geom_segment(data = tata_copy_A, aes(x = start, xend = end, y = ypos, yend = ypos), colour = "#315B8A", alpha = 0.65, linewidth = 1.3, lineend = "butt") +
  geom_point(data = point_A, aes(x = mid, y = ypos), shape = 17, size = 1.65, colour = "firebrick") +
  scale_x_continuous(limits = c(0, length_A), breaks = seq(0, length_A, 100), expand = expansion(add = c(0, 5))) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  labs(x = "Position within 5′ LTR (bp)", y = NULL, title = paste0("Clade A (n = ", nrow(line_A), ")"), tag = "C") + theme_fimo +
  theme(plot.title = element_text(colour = clade_colors[["cladeA_cons"]]))
pS16D <- ggplot() + geom_segment(data = line_B, aes(x = 1, xend = LTR_length, y = ypos, yend = ypos), colour = "grey80", linewidth = 0.35) +
  geom_segment(data = tata_copy_B, aes(x = start, xend = end, y = ypos, yend = ypos), colour = "#315B8A", alpha = 0.65, linewidth = 1.3, lineend = "butt") +
  geom_point(data = point_B, aes(x = mid, y = ypos), shape = 17, size = 1.65, colour = "firebrick") +
  scale_x_continuous(limits = c(0, length_B), breaks = seq(0, length_B, 100), expand = expansion(add = c(0, 5))) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  labs(x = "Position within 5′ LTR (bp)", y = NULL, title = paste0("Clade B (n = ", nrow(line_B), ")"), tag = "D") + theme_fimo +
  theme(plot.title = element_text(colour = clade_colors[["cladeB_cons"]]))

# Combine A-D with the same width per bp in all four panels: A on the left, B on the top right, C and D below B.
assemble_s16_shared_bp <- function(plots, width_in = 12, height_in = 17) {
  metric_file <- tempfile("s16_layout_metrics_", fileext = ".png")
  grDevices::png(metric_file, width = width_in, height = height_in,
                 units = "in", res = 96, type = "cairo")
  metric_device <- grDevices::dev.cur()
  on.exit({grDevices::dev.off(metric_device); unlink(metric_file)}, add = TRUE)
  grobs <- lapply(plots, ggplot2::ggplotGrob)
  spans <- vapply(plots, function(p) diff(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range), numeric(1))
  panel <- lapply(grobs, function(g) g$layout[g$layout$name == "panel", ])
  extra_width <- vapply(seq_along(grobs), function(i) {
    j <- seq.int(panel[[i]]$l, panel[[i]]$r)
    grid::convertWidth(sum(grobs[[i]]$widths[-j]), "mm", valueOnly = TRUE)
  }, numeric(1))
  extra_height <- vapply(seq_along(grobs), function(i) {
    j <- seq.int(panel[[i]]$t, panel[[i]]$b)
    grid::convertHeight(sum(grobs[[i]]$heights[-j]), "mm", valueOnly = TRUE)
  }, numeric(1))
  width_mm <- width_in * 25.4
  outer_top_mm <- 2; outer_bottom_mm <- 4
  height_mm <- height_in * 25.4 - outer_top_mm - outer_bottom_mm
  mm_per_bp <- min(
    (width_mm - sum(extra_width[c(1, 3, 4)])) / sum(spans[c(1, 3, 4)]),
    (width_mm - sum(extra_width[c(1, 2)])) / sum(spans[c(1, 2)]))
  row_heights <- height_mm * c(1.25, 7) / 8.25
  plot_heights <- c(height_mm, row_heights[1], row_heights[2], row_heights[2])
  plot_widths <- extra_width + spans * mm_per_bp
  for (i in seq_along(grobs)) {
    grobs[[i]]$widths[seq.int(panel[[i]]$l, panel[[i]]$r)] <- grid::unit(spans[i] * mm_per_bp, "mm")
    grobs[[i]]$heights[seq.int(panel[[i]]$t, panel[[i]]$b)] <- grid::unit(plot_heights[i] - extra_height[i], "mm")
  }
  # Leave the unused part of B's row blank rather than stretching its 500 bp axis.
  top_right <- gtable::gtable(widths = grid::unit(c(plot_widths[2],
                                                    max(0, sum(plot_widths[c(3, 4)]) - plot_widths[2])), "mm"),
                              heights = grid::unit(row_heights[1], "mm"))
  top_right <- gtable::gtable_add_grob(top_right, grobs[[2]], t = 1, l = 1, clip = "off")
  layout <- gtable::gtable(widths = grid::unit(plot_widths[c(1, 3, 4)], "mm"),
                           heights = grid::unit(row_heights, "mm"))
  layout <- gtable::gtable_add_grob(layout, grobs[[1]], t = 1, b = 2, l = 1, clip = "off")
  layout <- gtable::gtable_add_grob(layout, top_right, t = 1, l = 2, r = 3, clip = "off")
  layout <- gtable::gtable_add_grob(layout, grobs[[3]], t = 2, l = 2, clip = "off")
  layout <- gtable::gtable_add_grob(layout, grobs[[4]], t = 2, l = 3, clip = "off")
  # Small margins at the top and bottom of the page.
  frame <- gtable::gtable(widths = grid::unit(width_mm, "mm"),
                          heights = grid::unit(c(outer_top_mm, height_mm, outer_bottom_mm), "mm"))
  frame <- gtable::gtable_add_grob(frame, layout, t = 2, l = 1, clip = "off")
  list(figure = patchwork::wrap_elements(full = frame) +
         ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0)),
       mm_per_bp = mm_per_bp, x_spans_bp = spans, panel_widths_mm = spans * mm_per_bp)
}
s16_shared_layout <- assemble_s16_shared_bp(list(pS16A, pS16B, pS16C, pS16D))
figS16 <- s16_shared_layout$figure
figS16
ggsave("FigureS16.pdf", figS16, width = 12, height = 17, bg = "white", device = cairo_pdf)


# Percentage of 5′ LTRs with at least one HSF1 site in clade A and clade B, compared with the other Ty3/gypsy elements.
data.frame(group = c("Clade A", "Clade B"), percent = hsf1_summary$n_with_site[2:3] / hsf1_summary$n_sequences[2:3] * 100, percent_other = hsf1_summary$n_with_site[4] / hsf1_summary$n_sequences[4] * 100,
           fold_enrichment = (hsf1_summary$n_with_site[2:3] / hsf1_summary$n_sequences[2:3]) / (hsf1_summary$n_with_site[4] / hsf1_summary$n_sequences[4]),
           odds_ratio = hsf1_summary$odds_ratio_vs_other[2:3], p = hsf1_summary$p_vs_other[2:3], p_adj = hsf1_summary$q_vs_other[2:3])


############################################################
#Part 6: coverage of the consensus by the copies (Figure 5A) and alignment of each copy to its consensus (Figure S14)

# BLASTN hits of each copy on the consensus of its own clade.
blast2 <- read.table("intactLTR_vs_consensus.blastn.tsv", header = FALSE, col.names = c("qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"))
blast2 <- blast2 %>% dplyr::mutate(qseqid = sub("\\([+-]\\)$", "", qseqid), sseqid = sub("[#[:space:]].*$", "", sseqid), cs = pmin(sstart, send), ce = pmax(sstart, send)) %>%
  dplyr::semi_join(members, by = c("qseqid", "sseqid"))

#Part 6A: nucleotide identity of each copy to the consensus of its clade (MAFFT alignments of Figure S14)
# Only positions with A/C/G/T in both sequences are compared; gaps and N are excluded.

calculate_tf2_copy_identity <- function(
    alignment_files = c(cladeA_cons = "cladeA_to_cons.aln.fa", cladeB_cons = "cladeB_to_cons.aln.fa")) {
  clades <- c("cladeA_cons", "cladeB_cons")
  result <- lapply(clades, function(clade) {
    aligned <- Biostrings::readDNAStringSet(alignment_files[[clade]])
    names(aligned) <- sub("\\(.*\\)$", "", sub("[[:space:]].*$", "", names(aligned)))
    mat <- as.matrix(aligned)
    reference <- mat[clade, ]
    ids <- setdiff(rownames(mat), clade)
    reference_valid <- reference %in% c("A", "C", "G", "T")
    rows <- lapply(ids, function(id) {
      bases <- mat[id, ]
      comparable <- reference_valid & bases %in% c("A", "C", "G", "T")
      n_compared <- sum(comparable)
      n_identical <- sum(reference[comparable] == bases[comparable])
      data.frame(clade = clade, copy_id = id, compared_bp = n_compared,
                 identical_bp = n_identical, mismatched_bp = n_compared - n_identical,
                 identity_pct = 100 * n_identical / n_compared,
                 compared_consensus_coverage_pct = 100 * n_compared / sum(reference_valid),
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  })
  result <- do.call(rbind, result)
  result <- result[order(match(result$clade, clades), result$copy_id), ]
  rownames(result) <- NULL
  result
}

tf2_copy_identity <- calculate_tf2_copy_identity()
tf2_clade_identity_summary <- tf2_copy_identity %>% dplyr::group_by(clade) %>%
  dplyr::summarise(n_copies = dplyr::n(), mean_identity_pct = mean(identity_pct),
                   min_identity_pct = min(identity_pct), max_identity_pct = max(identity_pct),
                   median_identity_pct = median(identity_pct), .groups = "drop")
print(tf2_copy_identity)             # One row and one identity for every copy.
print(tf2_clade_identity_summary)    # Per-clade mean, minimum, maximum and median.


# Proportion of copies covering each consensus position (overlapping hits of one copy merged first; 5 bp at both ends not shown).
blast2_A <- blast2 %>% dplyr::filter(sseqid == "cladeA_cons")
blast2_B <- blast2 %>% dplyr::filter(sseqid == "cladeB_cons")
covered_A <- unlist(IRanges::reduce(split(IRanges(blast2_A$cs, blast2_A$ce), blast2_A$qseqid)))
covered_B <- unlist(IRanges::reduce(split(IRanges(blast2_B$cs, blast2_B$ce), blast2_B$qseqid)))
coverage_data <- rbind(data.frame(sseqid = "cladeA_cons", pos = 1:cons_len[["cladeA_cons"]], prop = as.integer(coverage(covered_A, width = cons_len[["cladeA_cons"]])) / nrow(member_A)),
                       data.frame(sseqid = "cladeB_cons", pos = 1:cons_len[["cladeB_cons"]], prop = as.integer(coverage(covered_B, width = cons_len[["cladeB_cons"]])) / nrow(member_B))) %>%
  dplyr::filter(pos > 5, pos <= cons_len[sseqid] - 5) %>% dplyr::mutate(sseqid = factor(sseqid, levels = clade_order))

# Common layout of the ggplot figures.
theme_TF2 <- theme_classic(base_size = 13, base_family = "Arial") +
  theme(axis.title = element_text(size = 14), axis.title.y = element_text(margin = margin(r = 7)), axis.text = element_text(size = 12, colour = "black"),
        axis.line = element_blank(), axis.ticks = element_line(colour = "black", linewidth = 0.4), panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "grey96", linewidth = 0.3), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), strip.background = element_rect(fill = "#D3D3D3", colour = "black", linewidth = 0.5),
        strip.text = element_text(size = 13, face = "bold", margin = margin(6, 0, 6, 0)), panel.spacing = grid::unit(1, "lines"),
        plot.title = element_text(size = 14, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12),
        plot.tag = element_text(size = 18, face = "bold"), plot.tag.position = "topleft", plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(6, 8, 6, 8))
# Number of 5′ LTRs with an HSF1 site upstream of a candidate TATA box (TATA boxes as in Figure S16; each sequence counted once).
hsf_summary_group_order <- c("Other intact Ty3/gypsy", "Clade A", "Clade B")
hsf_summary_sequences <- list("Other intact Ty3/gypsy" = seq_other, "Clade A" = seq5_A, "Clade B" = seq5_B)
hsf_summary_hsf_ids <- list("Other intact Ty3/gypsy" = unique(point_other$sequence_name),
                            "Clade A" = unique(point_A$sequence_name), "Clade B" = unique(point_B$sequence_name))
hsf_summary_tata_ids <- list("Other intact Ty3/gypsy" = unique(tata_other$sequence_name),
                             "Clade A" = unique(tata_copy_A$sequence_name), "Clade B" = unique(tata_copy_B$sequence_name))
# HSF1-TATA gap = number of bases between the end of the HSF1 site and the start of the TATA box (also used in Figure 5B).
hsf1_tata_gap_bp <- function(hsf_end, tata_start) tata_start - hsf_end - 1L
hsf_summary_downstream_pairs <- function(hsf_hits, tata_hits, group_name) {
  h <- hsf_hits %>% dplyr::transmute(sequence_name, hsf_start = start, hsf_stop = stop,
                                     hsf_strand = strand, hsf_sequence = matched_sequence, hsf_pvalue = p.value)
  t <- tata_hits %>% dplyr::transmute(sequence_name, tata_start = start, tata_end = end)
  dplyr::inner_join(h, t, by = "sequence_name", relationship = "many-to-many") %>%
    dplyr::filter(tata_start > hsf_stop) %>%
    dplyr::mutate(group = group_name, gap_bp = hsf1_tata_gap_bp(hsf_stop, tata_start)) %>%
    dplyr::select(group, dplyr::everything())
}
hsf_summary_pairs <- dplyr::bind_rows(
  hsf_summary_downstream_pairs(point_other, tata_other, "Other intact Ty3/gypsy"),
  hsf_summary_downstream_pairs(point_A, tata_copy_A, "Clade A"),
  hsf_summary_downstream_pairs(point_B, tata_copy_B, "Clade B"))
hsf_summary_presence <- dplyr::bind_rows(lapply(hsf_summary_group_order, function(g) {
  seq_ids <- names(hsf_summary_sequences[[g]])
  data.frame(group = g, sequence_name = seq_ids,
             has_HSF1 = seq_ids %in% hsf_summary_hsf_ids[[g]], has_TATA = seq_ids %in% hsf_summary_tata_ids[[g]],
             has_HSF1_then_TATA = seq_ids %in% unique(hsf_summary_pairs$sequence_name[hsf_summary_pairs$group == g]))
}))
hsf_summary_summary <- hsf_summary_presence %>% dplyr::group_by(group) %>%
  dplyr::summarise(n_sequences = dplyr::n(), n_HSF1 = sum(has_HSF1), n_TATA = sum(has_TATA),
                   n_both = sum(has_HSF1 & has_TATA),
                   n_HSF1_then_TATA = sum(has_HSF1_then_TATA), .groups = "drop") %>%
  dplyr::mutate(group = factor(group, levels = hsf_summary_group_order),
                percent_HSF1 = n_HSF1 / n_sequences * 100,
                percent_TATA = n_TATA / n_sequences * 100,
                percent_both = n_both / n_sequences * 100,
                percent_HSF1_then_TATA = n_HSF1_then_TATA / n_sequences * 100) %>% dplyr::arrange(group)
# Fisher's exact tests: 5′ LTRs with an HSF1 site upstream of a TATA box in clade A or clade B vs the other Ty3/gypsy elements
# (BH over the two tests). The Hsp70 genes are the positive control.
hsf1_tata_background <- hsf_summary_summary[as.character(hsf_summary_summary$group) == "Other intact Ty3/gypsy", ]
hsf1_tata_fisher_tests <- do.call(rbind, lapply(c("Clade A", "Clade B"), function(group_name) {
  dat <- hsf_summary_summary[as.character(hsf_summary_summary$group) == group_name, ]
  present <- dat$n_HSF1_then_TATA
  absent <- dat$n_sequences - present
  background_present <- hsf1_tata_background$n_HSF1_then_TATA
  background_absent <- hsf1_tata_background$n_sequences - background_present
  contingency <- matrix(c(present, absent, background_present, background_absent),
                        nrow = 2, byrow = TRUE, dimnames = list(c(group_name, "Other intact Ty3/gypsy"),
                                                                c("Arrangement present", "Arrangement absent")))
  test <- stats::fisher.test(contingency, alternative = "two.sided")
  data.frame(comparison = paste(group_name, "vs other intact Ty3/gypsy elements"),
             n_present = present, n_absent = absent, n_total = dat$n_sequences,
             other_present = background_present, other_absent = background_absent,
             other_total = hsf1_tata_background$n_sequences,
             odds_ratio = unname(test$estimate), odds_ratio_CI_low = test$conf.int[1],
             odds_ratio_CI_high = test$conf.int[2], p_value = test$p.value,
             stringsAsFactors = FALSE)
}))
hsf1_tata_fisher_tests$FDR_BH <- stats::p.adjust(hsf1_tata_fisher_tests$p_value, method = "BH")
hsf1_tata_HSP70_pairs <- hsf_summary_downstream_pairs(point_HSP70, tata_HSP70, "Hsp70 upstream")
hsf1_tata_positive_control <- data.frame(group = "Hsp70 upstream (500 bp)",
                                         n_sequences = length(seq_HSP70), n_HSF1_then_TATA = length(unique(hsf1_tata_HSP70_pairs$sequence_name)))
print(hsf1_tata_fisher_tests)
print(hsf1_tata_positive_control)


# The seven other elements similar to clade A or B, boxed in Figure S16A.
hsf_summary_unassigned <- s16_targets
hsf_summary_unassigned_check <- dplyr::left_join(hsf_summary_unassigned,
                                                 hsf_summary_presence %>% dplyr::filter(group == "Other intact Ty3/gypsy"), by = "sequence_name")
# HSF1-TATA gap of each element (the shortest one if an element has several HSF1 sites upstream of a TATA box) in clade A,
# clade B, the seven other elements similar to clade A or B, and the four remaining other elements.
hsf1_tata_distance_groups <- c("Clade A", "Clade B", "Other: unassigned A/B", "Other: remaining")
hsf1_tata_gap_pairs <- hsf_summary_pairs %>%
  dplyr::mutate(distance_group = dplyr::case_when(
    group == "Clade A" ~ "Clade A", group == "Clade B" ~ "Clade B",
    sequence_name %in% hsf_summary_unassigned$sequence_name ~ "Other: unassigned A/B",
    TRUE ~ "Other: remaining"),
    distance_group = factor(distance_group, levels = hsf1_tata_distance_groups)) %>%
  dplyr::left_join(hsf_summary_unassigned %>% dplyr::rename(unassigned_clade = clade), by = "sequence_name")
hsf1_tata_gap_per_element <- hsf1_tata_gap_pairs %>%
  dplyr::group_by(distance_group, sequence_name) %>%
  dplyr::mutate(n_qualifying_pairs = dplyr::n()) %>%
  dplyr::arrange(gap_bp, hsf_pvalue, tata_start, hsf_start, .by_group = TRUE) %>%
  dplyr::slice_head(n = 1) %>% dplyr::ungroup()
hsf1_tata_gap_summary <- hsf1_tata_gap_per_element %>%
  dplyr::group_by(distance_group) %>%
  dplyr::summarise(n_elements = dplyr::n(), mean_gap_bp = mean(gap_bp),
                   median_gap_bp = median(gap_bp), min_gap_bp = min(gap_bp), max_gap_bp = max(gap_bp),
                   .groups = "drop")
hsf1_tata_gap_frequencies <- hsf1_tata_gap_per_element %>%
  dplyr::count(distance_group, gap_bp, name = "n_elements")
hsf1_tata_gap_all_pairs_summary <- hsf1_tata_gap_pairs %>%
  dplyr::group_by(distance_group) %>%
  dplyr::summarise(n_elements = dplyr::n_distinct(sequence_name), n_pairs = dplyr::n(),
                   mean_gap_bp = mean(gap_bp), median_gap_bp = median(gap_bp),
                   min_gap_bp = min(gap_bp), max_gap_bp = max(gap_bp), .groups = "drop")
hsf1_tata_gap_other_elements <- hsf1_tata_gap_per_element %>%
  dplyr::filter(group == "Other intact Ty3/gypsy") %>%
  dplyr::select(distance_group, sequence_name, unassigned_clade, hsf_start, hsf_stop,
                hsf_strand, tata_start, tata_end, gap_bp)
print(hsf1_tata_gap_summary)
print(hsf1_tata_gap_frequencies)
print(hsf1_tata_gap_all_pairs_summary)
print(hsf1_tata_gap_other_elements)


clade_labels <- c(cladeA_cons = "Clade A", cladeB_cons = "Clade B")
feature_colors <- c(LTR = "#4D4D4D", GAG = "#4A9E6B", INT = "#C9484F", RT = "#4A72B0", RNaseH = "#8A74B5")
p5A <- ggplot() + geom_line(data = coverage_data, aes(x = pos, y = prop, colour = sseqid), linewidth = 0.8) +
  geom_segment(data = data.frame(sseqid = factor(clade_order, levels = clade_order), end = as.numeric(cons_len)), aes(x = 1, xend = end, y = 0.25, yend = 0.25), colour = "grey40", linewidth = 0.6) +
  geom_rect(data = features, aes(xmin = start, xmax = end, ymin = 0.16, ymax = 0.34, fill = feature), colour = "black", linewidth = 0.35) +
  geom_text(data = dplyr::filter(features, feature != "RNaseH"), aes(x = mid, y = 0.25, label = label), colour = "white", fontface = "bold", size = 10, size.unit = "pt") +
  geom_text(data = dplyr::filter(features, feature == "RNaseH"), aes(x = mid, y = 0.25, label = label), colour = "white", fontface = "bold", size = 8, size.unit = "pt") +
  geom_segment(data = hsf_sites, aes(x = mid, xend = mid, y = 0.14, yend = 0.21), colour = "firebrick", linewidth = 0.9) +
  geom_segment(data = hsf_sites, aes(x = mid, xend = mid, y = 0.29, yend = 0.36), colour = "firebrick", linewidth = 0.9) +
  geom_segment(data = hsf_sites, aes(x = mid + 350, xend = mid, y = 0.54, yend = 0.37), colour = "firebrick", linewidth = 0.5, arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")) +
  geom_text(data = hsf_sites, aes(x = mid + 350, y = 0.61), label = "Predicted HSF1 site", colour = "firebrick", hjust = 0, size = 10, size.unit = "pt", fontface = "bold") +
  scale_colour_manual(values = clade_colors, guide = "none") + scale_fill_manual(values = feature_colors, guide = "none") + scale_size_identity() +
  scale_x_continuous(breaks = seq(0, max(cons_len), 1000), expand = expansion(mult = c(0.005, 0.005))) +
  scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1), labels = scales::label_percent(), expand = expansion(mult = c(0, 0.02))) +
  facet_wrap(~sseqid, ncol = 1, scales = "free_x", labeller = as_labeller(clade_labels)) +
  labs(x = "Position on consensus (bp)", y = "Intact copies covering position (%)", tag = "A") + theme_TF2 + theme(panel.grid.major = element_blank())

# Figure S14: MAFFT alignment; grey = match, red = mismatch, white = gap or N.
# Copies are ordered by the similarity of their gap pattern (5 bp at both ends not shown).
aln_A <- as.matrix(readDNAStringSet("cladeA_to_cons.aln.fa"))
aln_B <- as.matrix(readDNAStringSet("cladeB_to_cons.aln.fa"))
keep_A <- 6:(ncol(aln_A) - 5)
keep_B <- 6:(ncol(aln_B) - 5)
mat_A <- aln_A[rownames(aln_A) != "cladeA_cons", keep_A, drop = FALSE]
mat_B <- aln_B[rownames(aln_B) != "cladeB_cons", keep_B, drop = FALSE]
dend_A <- reorder(as.dendrogram(hclust(dist(mat_A == "-", method = "manhattan"), method = "average")), wts = rowSums(mat_A == "-"), agglo.FUN = mean)
dend_B <- reorder(as.dendrogram(hclust(dist(mat_B == "-", method = "manhattan"), method = "average")), wts = rowSums(mat_B == "-"), agglo.FUN = mean)
mat_A <- mat_A[order.dendrogram(dend_A), , drop = FALSE]
mat_B <- mat_B[order.dendrogram(dend_B), , drop = FALSE]
alignment_A <- expand.grid(copy = seq_len(nrow(mat_A)), pos = keep_A)
alignment_B <- expand.grid(copy = seq_len(nrow(mat_B)), pos = keep_B)
# Gaps and N are white, as they are excluded from the identity calculation.
classify_tf2_alignment_bases <- function(copy_bases, consensus_bases) {
  reference_bases <- matrix(consensus_bases, nrow(copy_bases), ncol(copy_bases), byrow = TRUE)
  comparable <- copy_bases %in% c("A", "C", "G", "T") &
    reference_bases %in% c("A", "C", "G", "T")
  ifelse(!comparable, "gap", ifelse(copy_bases == reference_bases, "match", "mismatch"))
}
alignment_A$type <- as.vector(classify_tf2_alignment_bases(mat_A, aln_A["cladeA_cons", keep_A]))
alignment_B$type <- as.vector(classify_tf2_alignment_bases(mat_B, aln_B["cladeB_cons", keep_B]))

pS14_aln_A <- ggplot(alignment_A, aes(x = pos, y = copy, fill = type)) + geom_raster() +
  scale_fill_manual(values = c(match = "grey70", mismatch = "#E41A1C", gap = "white"), breaks = c("match", "mismatch"), labels = c("Match", "Mismatch")) +
  scale_x_continuous(limits = c(0, cons_len[["cladeA_cons"]]), expand = c(0, 0)) + scale_y_reverse(expand = c(0, 0)) +
  labs(x = NULL, y = "Intact copies", fill = NULL, title = paste0("Clade A (n = ", nrow(mat_A), ")"), tag = "A") + theme_TF2 +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.title = element_text(colour = clade_colors[["cladeA_cons"]]))
pS14_aln_B <- ggplot(alignment_B, aes(x = pos, y = copy, fill = type)) + geom_raster() +
  scale_fill_manual(values = c(match = "grey70", mismatch = "#E41A1C", gap = "white"), breaks = c("match", "mismatch"), labels = c("Match", "Mismatch")) +
  scale_x_continuous(limits = c(0, cons_len[["cladeB_cons"]]), expand = c(0, 0)) + scale_y_reverse(expand = c(0, 0)) +
  labs(x = NULL, y = "Intact copies", fill = NULL, title = paste0("Clade B (n = ", nrow(mat_B), ")"), tag = "B") + theme_TF2 +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.title = element_text(colour = clade_colors[["cladeB_cons"]]))
pS14_track_A <- ggplot(features %>% dplyr::filter(sseqid == "cladeA_cons")) +
  annotate("segment", x = 1, xend = cons_len[["cladeA_cons"]], y = 0.5, yend = 0.5, colour = "grey40", linewidth = 0.4) +
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = feature), colour = "black", linewidth = 0.3) +
  geom_text(data = dplyr::filter(features, sseqid == "cladeA_cons", feature != "RNaseH"), aes(x = mid, y = 0.5, label = label), colour = "white", fontface = "bold", size = 12, size.unit = "pt") +
  geom_segment(data = dplyr::filter(features, sseqid == "cladeA_cons", feature == "RNaseH"), aes(x = mid, xend = mid, y = 1, yend = 1.24), colour = "#6F5890", linewidth = 0.35) +
  geom_text(data = dplyr::filter(features, sseqid == "cladeA_cons", feature == "RNaseH"), aes(x = mid, y = 1.5, label = label), colour = "#6F5890", fontface = "bold", size = 11, size.unit = "pt") + scale_size_identity() + scale_fill_manual(values = feature_colors, guide = "none") +
  scale_x_continuous(limits = c(0, cons_len[["cladeA_cons"]]), breaks = seq(0, cons_len[["cladeA_cons"]], 1000), expand = c(0, 0)) +
  coord_cartesian(ylim = c(-0.1, 1.85)) + labs(x = NULL, y = NULL) + theme_TF2 +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank())
pS14_track_B <- ggplot(features %>% dplyr::filter(sseqid == "cladeB_cons")) +
  annotate("segment", x = 1, xend = cons_len[["cladeB_cons"]], y = 0.5, yend = 0.5, colour = "grey40", linewidth = 0.4) +
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = feature), colour = "black", linewidth = 0.3) +
  geom_text(data = dplyr::filter(features, sseqid == "cladeB_cons", feature != "RNaseH"), aes(x = mid, y = 0.5, label = label), colour = "white", fontface = "bold", size = 12, size.unit = "pt") +
  geom_segment(data = dplyr::filter(features, sseqid == "cladeB_cons", feature == "RNaseH"), aes(x = mid, xend = mid, y = 1, yend = 1.24), colour = "#6F5890", linewidth = 0.35) +
  geom_text(data = dplyr::filter(features, sseqid == "cladeB_cons", feature == "RNaseH"), aes(x = mid, y = 1.5, label = label), colour = "#6F5890", fontface = "bold", size = 11, size.unit = "pt") + scale_size_identity() + scale_fill_manual(values = feature_colors, guide = "none") +
  scale_x_continuous(limits = c(0, cons_len[["cladeB_cons"]]), breaks = seq(0, cons_len[["cladeB_cons"]], 1000), expand = c(0, 0)) +
  coord_cartesian(ylim = c(-0.1, 1.85)) + labs(x = "Position on consensus (bp)", y = NULL) + theme_TF2 +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank())

pS14_aln_A <- pS14_aln_A + theme(axis.title.y = element_text(size = 14), plot.title = element_text(size = 16), plot.tag = element_text(size = 20))
pS14_aln_B <- pS14_aln_B + theme(axis.title.y = element_text(size = 14), plot.title = element_text(size = 16), plot.tag = element_text(size = 20))
pS14_track_A <- pS14_track_A + theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14))
pS14_track_B <- pS14_track_B + theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14))

# Legend in its own row on top; the alignment panels are scaled to the number of copies.
figS14 <- wrap_plots(guide_area(), pS14_aln_A, pS14_track_A, pS14_aln_B, pS14_track_B, ncol = 1, heights = c(0.45, 6, 1, 6 * nrow(mat_B) / nrow(mat_A), 1), guides = "collect") &
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0), legend.text = element_text(size = 13))
figS14
ggsave("FigureS14.pdf", figS14, width = 9, height = 9, bg = "white", device = cairo_pdf)

############################################################

#Part 6B: dot plot of the 5′ LTRs of the two consensus sequences (Figure 5B)

# Dot plot: a segment is drawn where a 15 bp window of the clade A 5′ LTR and a 15 bp window of the clade B 5′ LTR share at
# least 12 identical bases (same strand, and clade A against the reverse complement of clade B).
dot_A <- strsplit(as.character(ltr5_consensus[["cladeA_5LTR"]]), "")[[1]]
dot_B <- strsplit(as.character(ltr5_consensus[["cladeB_5LTR"]]), "")[[1]]
dot_B_rc <- strsplit(as.character(reverseComplement(ltr5_consensus[["cladeB_5LTR"]])), "")[[1]]
same_fw <- outer(dot_A, dot_B, "==") & outer(dot_A %in% c("A", "C", "G", "T"), dot_B %in% c("A", "C", "G", "T"), "&")
same_rc <- outer(dot_A, dot_B_rc, "==") & outer(dot_A %in% c("A", "C", "G", "T"), dot_B_rc %in% c("A", "C", "G", "T"), "&")
nx <- length(dot_A) - 14
ny <- length(dot_B) - 14
window_fw <- 0
window_rc <- 0
for (k in 0:14) {
  window_fw <- window_fw + same_fw[k + 1:nx, k + 1:ny]
  window_rc <- window_rc + same_rc[k + 1:nx, k + 1:ny]
}
dot_fw <- as.data.frame(which(window_fw >= 12, arr.ind = TRUE))
dot_rc <- as.data.frame(which(window_rc >= 12, arr.ind = TRUE))
# Clade A on x, clade B on y; reverse-complement windows keep the forward coordinates of clade B.
dot_hits <- rbind(data.frame(x = dot_fw$row, xend = dot_fw$row + 14, y = dot_fw$col, yend = dot_fw$col + 14, type = "Same strand"),
                  data.frame(x = dot_rc$row, xend = dot_rc$row + 14, y = length(dot_B) - dot_rc$col + 1, yend = length(dot_B) - dot_rc$col - 13, type = "Reverse complement"))

# HSF1 site (the most frequent motif sequence of each clade, Part 5) and candidate TATA box (TATAWAW) in the consensus 5′ LTRs,
# and the gap between them.
hsf5_A <- matchPattern(dominant_A$cons_sequence, ltr5_consensus[["cladeA_5LTR"]])
hsf5_B <- matchPattern(dominant_B$cons_sequence, ltr5_consensus[["cladeB_5LTR"]])
tata_A <- matchPattern("TATAWAW", ltr5_consensus[["cladeA_5LTR"]], fixed = "subject")
tata_B <- matchPattern("TATAWAW", ltr5_consensus[["cladeB_5LTR"]], fixed = "subject")
motif_positions <- data.frame(clade = c("Clade A", "Clade B"), LTR_length = c(length(dot_A), length(dot_B)), hsf_start = c(start(hsf5_A), start(hsf5_B)),
                              hsf_end = c(end(hsf5_A), end(hsf5_B)), tata_start = c(start(tata_A), start(tata_B)), tata_end = c(end(tata_A), end(tata_B)))
motif_positions$gap_bp <- hsf1_tata_gap_bp(motif_positions$hsf_end, motif_positions$tata_start)
motif_positions$motif_id <- hsf1_motif_id
motif_positions$hsf_sequence <- c(dominant_A$cons_sequence, dominant_B$cons_sequence)
motif_positions$consensus_hsf_start <- c(start(hsf_A), start(hsf_B))
motif_positions$consensus_hsf_end <- c(end(hsf_A), end(hsf_B))
motif_positions$consensus_hsf_mid <- c((start(hsf_A) + end(hsf_A)) / 2, (start(hsf_B) + end(hsf_B)) / 2)
write.table(motif_positions, "Figure5_HSF1_MA04862_positions.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
motif_positions
gap_A <- motif_positions$gap_bp[1]
gap_B <- motif_positions$gap_bp[2]

# Arcs marking the HSF1-TATA gaps.
arc_A <- data.frame(x = end(hsf5_A) + (start(tata_A) - end(hsf5_A)) * seq(0, 1, length.out = 100), y = 18 - 5 * sin(pi * seq(0, 1, length.out = 100)))
arc_B <- data.frame(x = 17 - 5 * sin(pi * seq(0, 1, length.out = 100)), y = end(hsf5_B) + (start(tata_B) - end(hsf5_B)) * seq(0, 1, length.out = 100))

# Arrows pointing to the TATA box and the HSF1 site.
tata_centre_x <- (start(tata_A) + end(tata_A)) / 2
tata_centre_y <- (start(tata_B) + end(tata_B)) / 2
tata_arrow_x <- tata_centre_x + (235 - tata_centre_x) * 12 / sqrt((235 - tata_centre_x)^2 + (255 - tata_centre_y)^2)
tata_arrow_y <- tata_centre_y + (255 - tata_centre_y) * 12 / sqrt((235 - tata_centre_x)^2 + (255 - tata_centre_y)^2)
hsf_centre_x <- (start(hsf5_A) + end(hsf5_A)) / 2
hsf_centre_y <- (start(hsf5_B) + end(hsf5_B)) / 2
hsf_arrow_x <- hsf_centre_x + (395 - hsf_centre_x) * 12 / sqrt((395 - hsf_centre_x)^2 + (125 - hsf_centre_y)^2)
hsf_arrow_y <- hsf_centre_y + (125 - hsf_centre_y) * 12 / sqrt((395 - hsf_centre_x)^2 + (125 - hsf_centre_y)^2)

p5B <- ggplot() +
  geom_vline(xintercept = start(tata_A), linetype = "dashed", colour = "#315B8A", linewidth = 0.35) +
  geom_hline(yintercept = start(tata_B), linetype = "dashed", colour = "#315B8A", linewidth = 0.35) +
  geom_vline(xintercept = end(hsf5_A), linetype = "dashed", colour = "firebrick", linewidth = 0.35) +
  geom_hline(yintercept = end(hsf5_B), linetype = "dashed", colour = "firebrick", linewidth = 0.35) +
  geom_segment(data = dot_hits, aes(x = x, y = y, xend = xend, yend = yend, colour = type), linewidth = 0.6) +
  annotate("rect", xmin = start(hsf5_A), xmax = end(hsf5_A), ymin = start(hsf5_B), ymax = end(hsf5_B), colour = "firebrick", fill = "firebrick", alpha = 0.08, linewidth = 0.3) +
  annotate("rect", xmin = start(tata_A), xmax = end(tata_A), ymin = start(tata_B), ymax = end(tata_B), colour = "#315B8A", fill = "#315B8A", alpha = 0.10, linewidth = 0.3) +
  annotate("segment", x = 235, y = 255, xend = tata_arrow_x, yend = tata_arrow_y, colour = "#315B8A", linewidth = 0.35, arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")) +
  annotate("text", x = 100, y = 270, label = "Predicted TATA box", hjust = 0, size = 10, size.unit = "pt", fontface = "bold", colour = "#315B8A") +
  annotate("segment", x = 395, y = 125, xend = hsf_arrow_x, yend = hsf_arrow_y, colour = "firebrick", linewidth = 0.35, arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")) +
  annotate("text", x = 520, y = 95, label = "Predicted HSF1 site", hjust = 1, size = 10, size.unit = "pt", fontface = "bold", colour = "firebrick") +
  geom_path(data = arc_A, aes(x, y), colour = "black", linewidth = 0.3) + geom_path(data = arc_B, aes(x, y), colour = "black", linewidth = 0.3) +
  annotate("segment", x = c(end(hsf5_A), start(tata_A)), xend = c(end(hsf5_A), start(tata_A)), y = 16, yend = 20, colour = "black", linewidth = 0.3) +
  annotate("segment", x = 15, xend = 19, y = c(end(hsf5_B), start(tata_B)), yend = c(end(hsf5_B), start(tata_B)), colour = "black", linewidth = 0.3) +
  annotate("text", x = end(hsf5_A) - 9, y = 18, label = paste0("HSF1-TATA gap: ", gap_A, " bp"), hjust = 1, vjust = 0.5, size = 9, size.unit = "pt") +
  annotate("text", x = 15, y = end(hsf5_B) - 19, label = paste0("HSF1-TATA gap: ", gap_B, " bp"), hjust = 0, vjust = 0.5, size = 9, size.unit = "pt") +
  scale_colour_manual(values = c("Same strand" = "black", "Reverse complement" = "grey65"), guide = "none") +
  annotate("rect", xmin = 0.03 * length(dot_A), xmax = 0.46 * length(dot_A), ymin = 0.85 * length(dot_B), ymax = 0.98 * length(dot_B), fill = "white", colour = NA) +
  annotate("segment", x = 0.045 * length(dot_A), xend = 0.11 * length(dot_A), y = c(0.95, 0.89) * length(dot_B), yend = c(0.95, 0.89) * length(dot_B), colour = c("black", "grey65"), linewidth = 0.6) +
  annotate("text", x = 0.13 * length(dot_A), y = c(0.95, 0.89) * length(dot_B), label = c("Same strand", "Reverse complement"), hjust = 0, size = 9, size.unit = "pt") +
  scale_x_continuous(breaks = c(seq(0, length(dot_A), 100), length(dot_A))) + scale_y_continuous(breaks = c(seq(0, length(dot_B), 100), length(dot_B))) +
  coord_fixed(ratio = 1, xlim = c(0, length(dot_A)), ylim = c(0, length(dot_B)), expand = FALSE) +
  labs(x = "Clade A 5′ LTR (bp)", y = "Clade B 5′ LTR (bp)", tag = "B") + theme_TF2 +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), panel.grid.major = element_blank(), legend.position = "none")


############################################################
#Part 7: insertion time (Figure 5C, Figure S17)

# Insertion time of each intact element estimated by LTR_retriever from the divergence between its two LTRs.
pass_list <- read.delim("GCA_037325665.1_US_GU_Lsax_2.0_genomic.fna.mod.pass.list", header = FALSE, comment.char = "#", quote = "")
colnames(pass_list) <- c("LTR_loc", "Category", "Motif", "TSD", "TSD_5", "TSD_3", "Internal", "Identity", "Strand", "SuperFamily", "TE_type", "Insertion_Time")

# Element ID as in the GFF3 (chr:start-end of the whole element including TSDs).
pass_list <- pass_list %>%
  dplyr::mutate(chr = sub(":.*", "", LTR_loc), pos1 = as.numeric(sub(".*:(\\d+)\\.\\..*", "\\1", LTR_loc)), pos2 = as.numeric(sub(".*\\.\\.", "", LTR_loc)),
                rr_start = dplyr::if_else(TSD_5 == "..", pmin(pos1, pos2), suppressWarnings(as.numeric(sub("\\.\\..*", "", TSD_5)))),
                rr_end = dplyr::if_else(TSD_3 == "..", pmax(pos1, pos2), suppressWarnings(as.numeric(sub(".*\\.\\.", "", TSD_3)))),
                qseqid = paste0(chr, ":", rr_start, "-", rr_end), insertion_time_Mya = as.numeric(Insertion_Time) / 1e6)

comparison_colors <- c("Other intact Gypsy LTRs" = "grey75", "Clade A" = "#D55E00", "Clade B" = "#0072B2")
# Clade A, clade B and all other intact Gypsy elements.
insertion_data <- pass_list %>%
  dplyr::mutate(comparison_group = dplyr::case_when(qseqid %in% member_A$qseqid ~ "Clade A", qseqid %in% member_B$qseqid ~ "Clade B", SuperFamily == "Gypsy" ~ "Other intact Gypsy LTRs")) %>%
  dplyr::filter(!is.na(comparison_group)) %>% dplyr::mutate(comparison_group = factor(comparison_group, levels = names(comparison_colors)))
table(insertion_data$comparison_group)

# Kruskal-Wallis test, then pairwise Wilcoxon tests (BH correction over the three comparisons), on all insertion times.
kruskal.test(insertion_time_Mya ~ comparison_group, data = insertion_data)
insertion_other <- insertion_data$insertion_time_Mya[insertion_data$comparison_group == "Other intact Gypsy LTRs"]
insertion_A <- insertion_data$insertion_time_Mya[insertion_data$comparison_group == "Clade A"]
insertion_B <- insertion_data$insertion_time_Mya[insertion_data$comparison_group == "Clade B"]
insertion_test <- data.frame(group1 = c("Other intact Gypsy LTRs", "Other intact Gypsy LTRs", "Clade A"), group2 = c("Clade A", "Clade B", "Clade B"), xmin = c(1, 1, 2), xmax = c(2, 3, 3),
                             p = c(wilcox.test(insertion_other, insertion_A, exact = FALSE)$p.value, wilcox.test(insertion_other, insertion_B, exact = FALSE)$p.value,
                                   wilcox.test(insertion_A, insertion_B, exact = FALSE)$p.value))
insertion_test$p_adj <- p.adjust(insertion_test$p, "BH")
insertion_test$label <- ifelse(insertion_test$p_adj < 0.001, "***", ifelse(insertion_test$p_adj < 0.01, "**", ifelse(insertion_test$p_adj < 0.05, "*", "ns")))
insertion_test

# Significant comparisons only; shorter brackets below the longer ones.
insertion_sig <- insertion_test %>% dplyr::filter(p_adj < 0.05) %>% dplyr::arrange(xmax - xmin) %>% dplyr::mutate(y = 1.55 + (dplyr::row_number() - 1) * 0.17)

p5C <- ggplot(insertion_data, aes(x = comparison_group, y = insertion_time_Mya, fill = comparison_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85, linewidth = 0.4) +
  geom_point(data = insertion_data %>% dplyr::filter(comparison_group != "Other intact Gypsy LTRs"), aes(colour = comparison_group), position = position_jitter(width = 0.12, height = 0, seed = 123), size = 1.4, alpha = 0.4) +
  scale_fill_manual(values = comparison_colors, guide = "none") + scale_colour_manual(values = comparison_colors, guide = "none") +
  scale_x_discrete(labels = c("Other intact Gypsy LTRs" = "Other intact\nTy3/gypsy LTRs", "Clade A" = "Clade A", "Clade B" = "Clade B")) +
  geom_segment(data = insertion_sig, aes(x = xmin, xend = xmax, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = insertion_sig, aes(x = xmin, xend = xmin, y = y, yend = y - 0.035), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = insertion_sig, aes(x = xmax, xend = xmax, y = y, yend = y - 0.035), inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = insertion_sig, aes(x = (xmin + xmax) / 2, y = y + 0.035, label = label), inherit.aes = FALSE, size = 14, size.unit = "pt", fontface = "bold") +
  coord_cartesian(ylim = c(0, 1.85)) + labs(x = NULL, y = "Estimated insertion time (Mya)", tag = "C") + theme_TF2 +
  theme(axis.text.x = element_text(size = 12, face = "bold"), panel.grid.major = element_blank())

# Figure S17: distribution of insertion times (bins of 0.05 Mya, relative frequency within each group).
hist_data <- insertion_data %>% dplyr::mutate(bin_mid = floor(insertion_time_Mya / 0.05) * 0.05 + 0.05 / 2) %>% dplyr::count(comparison_group, bin_mid) %>%
  dplyr::group_by(comparison_group) %>% dplyr::mutate(relative_frequency = n / sum(n)) %>% dplyr::ungroup()
windowsFonts(Arial = windowsFont("Arial"))
s17_group_counts <- insertion_data %>% dplyr::count(comparison_group, name = "n_elements")
s17_group_labels <- stats::setNames(
  paste0(ifelse(as.character(s17_group_counts$comparison_group) == "Other intact Gypsy LTRs",
                "Other intact Ty3/gypsy", as.character(s17_group_counts$comparison_group)),
         " (n = ", s17_group_counts$n_elements, ")"),
  as.character(s17_group_counts$comparison_group))
figS17 <- ggplot(hist_data, aes(x = bin_mid, y = relative_frequency, fill = comparison_group)) + geom_col(width = 0.05, colour = "black", linewidth = 0.3) +
  facet_wrap(~comparison_group, ncol = 1, labeller = as_labeller(s17_group_labels)) + scale_fill_manual(values = comparison_colors, guide = "none") +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), expand = expansion(mult = c(0, 0.02))) + coord_cartesian(ylim = c(0, 0.8)) +
  labs(x = "Estimated insertion time (Mya)", y = "Relative frequency") + theme_TF2
figS17
ggsave("FigureS17.pdf", figS17, width = 8, height = 9.5, device = cairo_pdf)

############################################################
#Part 8: relative copy number from pool-seq depth (Figure S18)

# Per-base pool-seq depth of the six pools (population x shore level) on the two consensus sequences and the two housekeeping genes.
depth <- read.delim("all_depth.txt", check.names = FALSE)
pools <- c("Northhigh", "Northlow", "Centralhigh", "Centrallow", "Southhigh", "Southlow")

# Housekeeping genes (single copy): exon positions converted to coordinates on the gene sequences used as reference.
hk <- read.delim("housekeeper_original.gtf", header = FALSE, quote = "", comment.char = "#")
colnames(hk) <- c("Chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
hk_pos <- hk %>% dplyr::filter(feature == "exon") %>%
  dplyr::mutate(origin = dplyr::case_when(Chr == "LG3_SUPER_4" ~ 48883704, Chr == "LG12_SUPER_3" ~ 54491592), Chr = dplyr::recode(Chr, "LG3_SUPER_4" = "g15582.t1", "LG12_SUPER_3" = "g13887.t1"),
                start = start - origin + 1, end = end - origin + 1) %>%
  dplyr::rowwise() %>% dplyr::mutate(pos = list(seq.int(start, end))) %>% tidyr::unnest(pos) %>% dplyr::ungroup() %>% dplyr::distinct(Chr, pos)

# Single-copy depth of each pool = mean exon depth of each housekeeping gene, averaged over the two genes.
hk_depth <- depth %>% dplyr::semi_join(hk_pos, by = c("Chr", "pos")) %>% tidyr::pivot_longer(dplyr::all_of(pools), names_to = "pool", values_to = "depth") %>%
  dplyr::group_by(Chr, pool) %>% dplyr::summarise(mean_depth = mean(depth), .groups = "drop") %>% dplyr::group_by(pool) %>% dplyr::summarise(housekeeper_depth = mean(mean_depth), .groups = "drop")
hk_depth

# RT domain of each consensus.
rt <- read.delim("tf2_RT.gtf", header = FALSE, quote = "", comment.char = "#")
colnames(rt) <- c("Chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
rt_pos <- rt %>% dplyr::filter(feature == "exon", Chr %in% clade_order) %>% dplyr::rowwise() %>% dplyr::mutate(pos = list(seq.int(start, end))) %>% tidyr::unnest(pos) %>% dplyr::ungroup() %>% dplyr::distinct(Chr, pos)

# Relative copy number = mean depth on the RT domain / single-copy depth, per pool.
copy_pool <- depth %>% dplyr::semi_join(rt_pos, by = c("Chr", "pos")) %>% tidyr::pivot_longer(dplyr::all_of(pools), names_to = "pool", values_to = "depth") %>%
  dplyr::group_by(Chr, pool) %>% dplyr::summarise(RT_depth = mean(depth), .groups = "drop") %>% dplyr::left_join(hk_depth, by = "pool") %>%
  dplyr::mutate(relative_copy_number = RT_depth / housekeeper_depth)
copy_pool

population_order <- c("North", "Central", "South")
copy_shore <- copy_pool %>% dplyr::mutate(population = factor(sub("(high|low)$", "", pool), levels = population_order),
                                          shore_level = factor(sub("^(North|Central|South)", "", pool), levels = c("high", "low"), labels = c("High shore", "Low shore")), Chr = factor(Chr, levels = clade_order))
copy_shore
windowsFonts(Arial = windowsFont("Arial"))
figS18 <- ggplot(copy_shore, aes(x = population, y = relative_copy_number, fill = Chr)) + geom_col(width = 0.65, colour = "black", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.1f", relative_copy_number)), vjust = -0.4, size = 11, size.unit = "pt") + scale_fill_manual(values = clade_colors, guide = "none") +
  facet_grid(shore_level ~ Chr, labeller = labeller(Chr = clade_labels)) + scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Relative copy number") + theme_TF2 + theme(strip.text.y = element_text(angle = 270, margin = margin(6, 6, 6, 6)))
figS18
ggsave("FigureS18.pdf", figS18, width = 9, height = 7, device = cairo_pdf)


############################################################
#Part 9: expression of the RT domain in each individual (Figure 5D, E)

# RNA-seq read counts on the RT domain of each consensus and number of mapped fragments of each individual.
rt_counts <- read.delim("tf2_RT_counts.txt", check.names = FALSE)
library_size <- read.delim("library_size.txt") %>% dplyr::filter(sample != "SPL47")
group <- read.table("../01counts_ind_filter/group.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_ids <- names(rt_counts)[-(1:6)]

stress_order <- c("No/short heat stress", "Long heat stress")
# 31°C for 3 hours = long heat stress; the other three treatments = no or short heat stress.
group <- group %>% dplyr::mutate(exposure_group = factor(dplyr::if_else(temperature == "31degree" & time == "3h", "Long heat stress", "No/short heat stress"), levels = stress_order),
                                 pool = paste0(population, shore_level), population = factor(population, levels = population_order))

# FPKM = fragments per kb of RT per million mapped fragments; expression per copy = FPKM / relative copy number of the pool
# (population x shore level) of the individual.
expression_data <- rt_counts %>% dplyr::select(Chr, Length, dplyr::all_of(sample_ids)) %>% tidyr::pivot_longer(dplyr::all_of(sample_ids), names_to = "id", values_to = "count") %>%
  dplyr::left_join(library_size, by = c("id" = "sample")) %>% dplyr::left_join(group, by = "id") %>% dplyr::left_join(copy_pool %>% dplyr::select(Chr, pool, relative_copy_number), by = c("Chr", "pool")) %>%
  dplyr::mutate(Chr = factor(Chr, levels = clade_order), FPKM = count / mapped_fragments * 1e6 / (Length / 1000), FPKM_per_copy = FPKM / relative_copy_number)

# Wilcoxon test between two populations within one clade and one exposure group.
test_populations <- function(clade, stress, pop1, pop2) {
  a <- expression_data[expression_data$Chr == clade & expression_data$exposure_group == stress & expression_data$population == pop1, ]
  b <- expression_data[expression_data$Chr == clade & expression_data$exposure_group == stress & expression_data$population == pop2, ]
  data.frame(Chr = clade, exposure_group = stress, population1 = pop1, population2 = pop2, n1 = nrow(a), n2 = nrow(b),
             p_FPKM = wilcox.test(a$FPKM, b$FPKM, exact = FALSE)$p.value, p_copy = wilcox.test(a$FPKM_per_copy, b$FPKM_per_copy, exact = FALSE)$p.value)
}

# Three pairwise comparisons per clade and exposure group; BH correction over the 12 tests of each measure.
population_tests <- rbind(test_populations("cladeA_cons", "No/short heat stress", "North", "Central"), test_populations("cladeA_cons", "No/short heat stress", "North", "South"),
                          test_populations("cladeA_cons", "No/short heat stress", "Central", "South"), test_populations("cladeA_cons", "Long heat stress", "North", "Central"),
                          test_populations("cladeA_cons", "Long heat stress", "North", "South"), test_populations("cladeA_cons", "Long heat stress", "Central", "South"),
                          test_populations("cladeB_cons", "No/short heat stress", "North", "Central"), test_populations("cladeB_cons", "No/short heat stress", "North", "South"),
                          test_populations("cladeB_cons", "No/short heat stress", "Central", "South"), test_populations("cladeB_cons", "Long heat stress", "North", "Central"),
                          test_populations("cladeB_cons", "Long heat stress", "North", "South"), test_populations("cladeB_cons", "Long heat stress", "Central", "South"))
population_tests$FDR_FPKM <- p.adjust(population_tests$p_FPKM, "BH")
population_tests$FDR_copy <- p.adjust(population_tests$p_copy, "BH")
population_tests$label_FPKM <- as.character(cut(population_tests$FDR_FPKM, c(0, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns"), right = FALSE))
population_tests$label_copy <- as.character(cut(population_tests$FDR_copy, c(0, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns"), right = FALSE))
population_tests

# Significance brackets above the boxes of the two populations compared.
test_steps <- expression_data %>% dplyr::group_by(Chr) %>%
  dplyr::summarise(step_FPKM = max(diff(range(log2(FPKM))), 0.1) * 0.1, step_copy = max(diff(range(log2(FPKM_per_copy))), 0.1) * 0.1, .groups = "drop")
test_tops <- expression_data %>% dplyr::group_by(Chr, exposure_group) %>% dplyr::summarise(top_FPKM = max(log2(FPKM)), top_copy = max(log2(FPKM_per_copy)), .groups = "drop")
population_brackets <- population_tests %>% dplyr::mutate(Chr = factor(Chr, levels = clade_order), exposure_group = factor(exposure_group, levels = stress_order)) %>%
  dplyr::left_join(test_steps, by = "Chr") %>% dplyr::left_join(test_tops, by = c("Chr", "exposure_group")) %>%
  dplyr::mutate(x = match(exposure_group, stress_order), x_left = x + (match(population1, population_order) - 2) * 0.8 / 3 + 0.02,
                x_right = x + (match(population2, population_order) - 2) * 0.8 / 3 - 0.02, y_FPKM = top_FPKM + step_FPKM, y_copy = top_copy + step_copy)
sig_FPKM <- population_brackets %>% dplyr::filter(FDR_FPKM < 0.05) %>% dplyr::group_by(Chr, exposure_group) %>%
  dplyr::arrange(x_right - x_left, .by_group = TRUE) %>% dplyr::mutate(y_FPKM = top_FPKM + step_FPKM * dplyr::row_number() * 0.8) %>% dplyr::ungroup()
sig_copy <- population_brackets %>% dplyr::filter(FDR_copy < 0.05) %>% dplyr::group_by(Chr, exposure_group) %>%
  dplyr::arrange(x_right - x_left, .by_group = TRUE) %>% dplyr::mutate(y_copy = top_copy + step_copy * dplyr::row_number() * 0.8) %>% dplyr::ungroup()

population_colors <- c(North = "#87aade", Central = "#ffe57f", South = "#ff9999")
p5D <- ggplot(expression_data, aes(x = exposure_group, y = log2(FPKM), fill = population)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.65, colour = "black", linewidth = 0.4, outlier.shape = NA) +
  geom_point(aes(colour = population), size = 1.6, alpha = 0.65, position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8, seed = 123)) +
  geom_segment(data = sig_FPKM, aes(x = x_left, xend = x_right, y = y_FPKM, yend = y_FPKM), inherit.aes = FALSE, linewidth = 0.35) +
  geom_segment(data = sig_FPKM, aes(x = x_left, xend = x_left, y = y_FPKM, yend = y_FPKM - step_FPKM * 0.15), inherit.aes = FALSE, linewidth = 0.35) +
  geom_segment(data = sig_FPKM, aes(x = x_right, xend = x_right, y = y_FPKM, yend = y_FPKM - step_FPKM * 0.15), inherit.aes = FALSE, linewidth = 0.35) +
  geom_text(data = sig_FPKM, aes(x = (x_left + x_right) / 2, y = y_FPKM + step_FPKM * 0.22, label = label_FPKM), inherit.aes = FALSE, size = 14, size.unit = "pt", fontface = "bold", vjust = 0.35) +
  scale_fill_manual(values = population_colors) + scale_colour_manual(values = population_colors, guide = "none") + scale_x_discrete(labels = c("No/short heat stress" = "No/short\nheat stress", "Long heat stress" = "Long\nheat stress"), expand = expansion(add = 0.5)) +
  facet_wrap(~Chr, nrow = 1, labeller = as_labeller(clade_labels)) + scale_y_continuous(expand = expansion(mult = c(0.04, 0.15))) +
  labs(x = NULL, y = expression(log[2]*"(FPKM)"), tag = "D") + theme_TF2 + theme(axis.text.x = element_text(face = "bold"))

p5E <- ggplot(expression_data, aes(x = exposure_group, y = log2(FPKM_per_copy), fill = population)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.65, colour = "black", linewidth = 0.4, outlier.shape = NA) +
  geom_point(aes(colour = population), size = 1.6, alpha = 0.65, position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8, seed = 123)) +
  geom_segment(data = sig_copy, aes(x = x_left, xend = x_right, y = y_copy, yend = y_copy), inherit.aes = FALSE, linewidth = 0.35) +
  geom_segment(data = sig_copy, aes(x = x_left, xend = x_left, y = y_copy, yend = y_copy - step_copy * 0.15), inherit.aes = FALSE, linewidth = 0.35) +
  geom_segment(data = sig_copy, aes(x = x_right, xend = x_right, y = y_copy, yend = y_copy - step_copy * 0.15), inherit.aes = FALSE, linewidth = 0.35) +
  geom_text(data = sig_copy, aes(x = (x_left + x_right) / 2, y = y_copy + step_copy * 0.22, label = label_copy), inherit.aes = FALSE, size = 14, size.unit = "pt", fontface = "bold", vjust = 0.35) +
  scale_fill_manual(values = population_colors) + scale_colour_manual(values = population_colors, guide = "none") + scale_x_discrete(labels = c("No/short heat stress" = "No/short\nheat stress", "Long heat stress" = "Long\nheat stress"), expand = expansion(add = 0.5)) +
  facet_wrap(~Chr, nrow = 1, labeller = as_labeller(clade_labels)) + scale_y_continuous(expand = expansion(mult = c(0.04, 0.15))) +
  labs(x = NULL, y = expression(log[2]*"(FPKM / relative copy number)"), tag = "E") + theme_TF2 + theme(axis.text.x = element_text(face = "bold"))


############################################################
#Figure 5: A on top, B and C in the middle, D and E at the bottom
theme_figure5 <- theme_TF2 + theme(axis.text = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11.5), axis.title.y = element_text(margin = margin(r = 2)), axis.text.y = element_text(margin = margin(r = 2)), strip.text = element_text(size = 12, face = "bold"),
                                   plot.tag = element_text(size = 16, face = "bold", hjust = 0, vjust = 1), plot.tag.position = c(0, 1), legend.text = element_text(size = 10.5), panel.spacing = grid::unit(0.8, "lines"), plot.margin = margin(6, 5, 5, 5))
figure5_A <- p5A + theme_figure5 + theme(legend.position = "none", panel.grid.major = element_blank())
figure5_B <- p5B + theme_figure5 + theme(legend.position = "none", panel.grid.major = element_blank())
figure5_C <- p5C + theme_figure5 + theme(axis.text.x = element_text(size = 9.5, face = "bold"), aspect.ratio = 1.08 * length(dot_B) / length(dot_A), legend.position = "none", panel.grid.major = element_blank())
figure5_D <- p5D + theme_figure5 + theme(axis.text.x = element_text(size = 9.5, face = "bold"), axis.title.y = element_text(margin = margin(r = 2)), legend.position = "bottom")
figure5_E <- p5E + theme_figure5 + theme(axis.text.x = element_text(size = 9.5, face = "bold"), axis.title.y = element_text(margin = margin(r = 2)), legend.position = "bottom")
figure5_B <- patchwork::free(figure5_B, type = "label", side = "l")
figure5_D <- patchwork::free(figure5_D, type = "label", side = "l")
figure5_BC <- (figure5_B | figure5_C) + plot_layout(widths = c(1.08, 1))
figure5_DE <- (figure5_D | figure5_E) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
figure5 <- (figure5_A / figure5_BC / figure5_DE) + plot_layout(heights = c(4.3, 4.0, 3.9))


############################################################
#Part 10: other intact Ty3/gypsy elements with both an HSF1 site and a TATA box (Figure S16A)

# Other intact Ty3/gypsy elements with both an HSF1 site (FIMO, p < 1e-4) and a candidate TATA box (TATAWAW).
hsf_tata_names <- intersect(unique(plot_other$sequence_name), unique(tata_other$sequence_name))
length(hsf_tata_names)

# One row per element: 5′ LTR length, and the positions of its HSF1 sites and TATA boxes (several sites separated by commas).
hsf_tata_other <- line_other %>% dplyr::filter(sequence_name %in% hsf_tata_names) %>% dplyr::select(sequence_name, LTR_length) %>%
  dplyr::left_join(plot_other %>% dplyr::group_by(sequence_name) %>%
                     dplyr::summarise(hsf_start = paste(start, collapse = ","), hsf_sequence = paste(matched_sequence, collapse = ","), .groups = "drop"), by = "sequence_name") %>%
  dplyr::left_join(tata_other %>% dplyr::group_by(sequence_name) %>% dplyr::summarise(tata_start = paste(start, collapse = ","), .groups = "drop"), by = "sequence_name")

# Intact element of each 5′ LTR (from its BED coordinates and the GFF3); identical LTR sequences can occur in several elements.
seq_hsf_tata_other <- seq_other[hsf_tata_names]
other_ltr_bed <- read.delim("otherLTR_5LTR.bed", header = FALSE, quote = "", col.names = c("chr", "start0", "end", "fasta_id", "score", "strand"))
other_ltr_source <- gff %>% dplyr::filter(type == "long_terminal_repeat") %>% dplyr::transmute(chr, ltr_start = start, ltr_end = end, parent = sub(".*Parent=([^;]+).*", "\\1", attribute))
ltr5_to_intact <- other_ltr_bed %>% dplyr::mutate(sequence_name = paste0(fasta_id, "(", strand, ")"), ltr_start = start0 + 1L, ltr_end = end) %>% dplyr::filter(sequence_name %in% hsf_tata_names) %>%
  dplyr::left_join(other_ltr_source, by = c("chr", "ltr_start", "ltr_end")) %>% dplyr::left_join(repeat_region %>% dplyr::select(parent, intact_id), by = "parent") %>% dplyr::select(intact_id, sequence_name) %>% dplyr::distinct()
table(table(ltr5_to_intact$sequence_name))

# Element position and strand, LTR_retriever classification and insertion time (Part 7), and the best BLASTN hit to the two Tf2
# consensus sequences (Part 3; NA = no hit to either consensus).
hsf_tata_other <- hsf_tata_other %>% dplyr::left_join(ltr5_to_intact, by = "sequence_name") %>%
  dplyr::left_join(repeat_region %>% dplyr::select(intact_id, chr, start, end, strand, parent), by = "intact_id") %>%
  dplyr::left_join(pass_list %>% dplyr::select(intact_id = qseqid, SuperFamily, TE_type, insertion_time_Mya), by = "intact_id") %>%
  dplyr::left_join(blast1_summary %>% dplyr::group_by(qseqid) %>% dplyr::slice_max(bitscore, n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>%
                     dplyr::select(intact_id = qseqid, best_consensus = sseqid, identity, aln_len, cov_intact), by = "intact_id") %>%
  dplyr::arrange(dplyr::desc(LTR_length))
hsf_tata_other

gff_hsf_tata <- gff[sub(".*(ID|Parent)=([^;]+).*", "\\2", gff$attribute) %in% hsf_tata_other$parent, ]

############################################################
#Part 11: general temperature response genes within 50 kb of clade A and clade B copies (Table S9)

analyse_tf2_flanks <- function(copy_members, full_gene_ranges, de_results, gene_annotation) {
  # General temperature response genes in any population (up or down). Genes annotated as Tf2 or overlapping a copy are
  # excluded from both the nearby and the background genes.
  de <- dplyr::bind_rows(de_results)
  background_all <- unique(de$gene[is.finite(de$padj) & is.finite(de$log2FoldChange)])
  DEG_union_all <- unique(de$gene[which(is.finite(de$padj) & is.finite(de$log2FoldChange) &
                                          de$padj < 0.05 & abs(de$log2FoldChange) >= 0.5)])
  ltr <- copy_members %>% dplyr::transmute(
    clade = dplyr::recode(sseqid, cladeA_cons = "Clade A", cladeB_cons = "Clade B"),
    intact_id = qseqid, chr, LTR_start = start, LTR_end = end)
  # Gene spans as in Part 1.
  all_gene_gr <- GenomicRanges::GRanges(
    as.character(GenomicRanges::seqnames(full_gene_ranges)),
    IRanges::IRanges(start(full_gene_ranges), end(full_gene_ranges)))
  all_gene_ids <- full_gene_ranges$gene_id
  body_gr <- GenomicRanges::GRanges(ltr$chr, IRanges::IRanges(ltr$LTR_start, ltr$LTR_end))
  common_seqlevels <- union(GenomeInfoDb::seqlevels(all_gene_gr), GenomeInfoDb::seqlevels(body_gr))
  GenomeInfoDb::seqlevels(all_gene_gr) <- common_seqlevels
  GenomeInfoDb::seqlevels(body_gr) <- common_seqlevels
  body_hits <- GenomicRanges::findOverlaps(all_gene_gr, body_gr, ignore.strand = TRUE)
  body_overlap_genes <- unique(all_gene_ids[S4Vectors::queryHits(body_hits)])
  tf2_annotated_genes <- gene_annotation$gene[grepl("\\bTf2(?:\\b|[-_])",
                                                    gene_annotation$description, ignore.case = TRUE, perl = TRUE)]
  excluded_genes <- union(body_overlap_genes, tf2_annotated_genes)
  background <- setdiff(background_all, excluded_genes)
  DEG_union <- intersect(DEG_union_all, background)
  eligible <- which(all_gene_ids %in% background)
  gene_gr <- all_gene_gr[eligible]
  gene_ids <- all_gene_ids[eligible]
  # 50 kb on each side of each copy.
  windows <- tidyr::crossing(ltr, side = c("left", "right")) %>%
    dplyr::mutate(
      window_start = dplyr::if_else(side == "left", pmax(1L, LTR_start - 50000L), LTR_end + 1L),
      window_end = dplyr::if_else(side == "left", LTR_start - 1L, LTR_end + 50000L)) %>%
    dplyr::filter(window_start <= window_end)
  window_gr <- GenomicRanges::GRanges(windows$chr,
                                      IRanges::IRanges(windows$window_start, windows$window_end))
  GenomeInfoDb::seqlevels(window_gr) <- common_seqlevels
  flank_hits <- GenomicRanges::findOverlaps(window_gr, gene_gr, ignore.strand = TRUE)
  nearby_genes <- data.frame(
    clade = windows$clade[S4Vectors::queryHits(flank_hits)],
    gene = gene_ids[S4Vectors::subjectHits(flank_hits)], stringsAsFactors = FALSE) %>%
    dplyr::distinct(clade, gene)
  nearby_DEGs <- nearby_genes %>% dplyr::filter(gene %in% DEG_union)
  
  # Two-sided Fisher's exact test for each clade; BH over the two clades.
  flank_enrichment <- dplyr::bind_rows(lapply(c("Clade A", "Clade B"), function(clade_name) {
    near <- nearby_genes$gene[nearby_genes$clade == clade_name]
    elsewhere <- setdiff(background, near)
    a <- sum(near %in% DEG_union); b <- length(near) - a
    c <- sum(elsewhere %in% DEG_union); d <- length(elsewhere) - c
    test <- stats::fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE), alternative = "two.sided")
    data.frame(clade = clade_name, window_kb = 50L,
               background_genes = length(background), background_DEGs = length(DEG_union),
               nearby_genes = length(near), nearby_DEGs = a,
               nearby_DEG_percent = 100 * a / length(near),
               elsewhere_genes = length(elsewhere), elsewhere_DEGs = c,
               elsewhere_DEG_percent = 100 * c / length(elsewhere),
               odds_ratio = unname(test$estimate), P_value = test$p.value)
  })) %>% dplyr::mutate(FDR_BH = stats::p.adjust(P_value, method = "BH"))
  
  list(enrichment = flank_enrichment, nearby_genes = nearby_genes, nearby_DEGs = nearby_DEGs)
}

tf2_flanks <- analyse_tf2_flanks(members, gene_ranges, list(deg_North, deg_Central, deg_South), annotation)
tf2_flanks$enrichment

figure5
ggsave("Figure5.pdf", figure5, width = 9, height = 12.2, device = cairo_pdf, bg = "white")

############################################################
#Part 12: expression of the HSF1-like gene (Figure S15)

# Counts and sample information.
counts <- read.delim("../01counts_ind_filter/counts_filter.txt", row.names = 1, check.names = FALSE)
group <- read.delim("../01counts_ind_filter/group.txt", stringsAsFactors = FALSE)
group <- group[match(colnames(counts), group$id), ]

# CPM of the HSF1-like gene MSTRG.5684.
group$CPM <- as.numeric(counts["MSTRG.5684", ]) / colSums(counts) * 1e6
group$gene <- "MSTRG.5684"
group$population <- factor(group$population, levels = c("North", "Central", "South"))
group$shore <- factor(group$shore_level, levels = c("high", "low"), labels = c("High shore", "Low shore"))
group$exposure_group <- factor(ifelse(group$temperature == "31degree" & group$time == "3h", "Long", "NoShort"), levels = c("NoShort", "Long"))

# Population colours as in the other figures.
population_colors <- c(North = "#87aade", Central = "#ffe57f", South = "#ff9999")

p <- ggplot(group, aes(exposure_group, CPM, color = population)) +
  geom_boxplot(position = position_dodge(width = 0.78), width = 0.66, outlier.shape = NA, linewidth = 0.6, fill = "white") +
  geom_point(aes(shape = shore, group = population), size = 1.8, position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.78, seed = 5684)) +
  facet_wrap(~gene) +
  scale_color_manual(values = population_colors) +
  scale_shape_manual(values = c("High shore" = 17, "Low shore" = 16)) +
  scale_x_discrete(labels = c(NoShort = "No/short heat stress", Long = "Long heat stress")) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.06))) +
  labs(x = NULL, y = "HSF1-like expression (CPM)", color = "Sampling location", shape = "Shore level") +
  guides(color = guide_legend(order = 1, override.aes = list(shape = 16)), shape = guide_legend(order = 2, override.aes = list(color = "grey30"))) +
  theme_bw(base_size = 11) +
  theme(axis.text = element_text(color = "black"), strip.background = element_rect(fill = "grey85", color = "grey35"), strip.text = element_text(size = 11), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(color = "grey92"), panel.grid.minor.y = element_line(color = "grey96"), legend.position = "right")

p

ggsave("FigureS15.pdf", p, width = 6.4, height = 4.6, bg = "white")
