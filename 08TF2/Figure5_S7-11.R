# =============================================================================
# Tf2 retrotransposons among the heat-induced genes (Figure 5, Figures S7-S11)
#
# Part 1: Tf2-annotated genes up-regulated at 31°C in at least one population (125 genes), and the intact LTR
#         retrotransposons (LTR_retriever) that they overlap (47 intact elements).
# Part 2: Neighbour-joining tree of the 47 elements, two clades (A and B) and their first consensus sequences (Figure S7).
# Part 3: All intact LTR elements of the genome are searched against the two consensus sequences (BLASTN, 80-80-80 rule)
#         to collect all members of clade A (138 copies) and clade B (64 copies); second consensus from all members.
# Part 4: Structure of the two consensus sequences (5′/3′ LTR and protein domains).
# Part 5: HSF1 binding motifs in the 5′ LTR of each copy and upstream of the Hsp70 genes (Figure S11).
# Part 6: Coverage of the consensus by the intact copies (Figure 5A) and copy-by-copy alignment (Figure S8).
# Part 7: Insertion time of clade A and clade B compared with the other intact Gypsy elements (Figure 5B, Figure S9).
# Part 8: Relative copy number of clade A and clade B in each population from pool-seq depth (Figure S10).
# Part 9: Expression of the reverse transcriptase (RT) of clade A and clade B in each individual (Figure 5C, 5D).
#
# Steps run outside R are given as comments where they occur (LTR_retriever, bedtools, BLASTN, TEsorter, MAFFT,
# read mapping, FIMO). Lines that write intermediate files are commented out: these files were produced once and are
# read back in the later steps.
#
# Input:  ../04DEseq2/*_generaltemp_allgenes.txt, ../05GO/gene_annotation.txt, ../05GO/Lsax_stringtie_merged.gtf,
#         ../01counts_ind_filter/group.txt, LTR_retriever output (*.pass.list, *.pass.list.gff3) and the files listed in the comments
# Output: Figure5.pdf, FigureS7.pdf, FigureS8.pdf, FigureS9.pdf, FigureS10.pdf, FigureS11.pdf
# =============================================================================

setwd("/data/08TF2/")
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
#Part 1: Tf2 genes up-regulated by heat and the intact LTR elements they belong to

# Genes up-regulated at 31°C (padj < 0.05, log2FC >= 1) in each population; union of the three populations.
deg_North <- read.delim("../04DEseq2/North_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
deg_Central <- read.delim("../04DEseq2/Central_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
deg_South <- read.delim("../04DEseq2/South_generaltemp_allgenes.txt", header = TRUE, stringsAsFactors = FALSE)
up_North <- deg_North$gene[which(deg_North$padj < 0.05 & deg_North$log2FoldChange >= 1)]
up_Central <- deg_Central$gene[which(deg_Central$padj < 0.05 & deg_Central$log2FoldChange >= 1)]
up_South <- deg_South$gene[which(deg_South$padj < 0.05 & deg_South$log2FoldChange >= 1)]
up_union <- unique(c(up_North, up_Central, up_South))

# Up-regulated genes annotated as Tf2 retrotransposon: 125 genes.
annotation <- read.delim("../05GO/gene_annotation.txt", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
tf2_genes <- unique(annotation$gene[annotation$gene %in% up_union & grepl("Tf2", annotation$description, ignore.case = TRUE)])
length(tf2_genes)

# Genomic position of each gene = span of its exons in the StringTie annotation.
gtf <- import("../05GO/Lsax_stringtie_merged.gtf")
gtf_exon <- gtf[gtf$type == "exon" & !is.na(gtf$gene_id) & gtf$gene_id != ""]
gene_ranges <- unlist(range(split(gtf_exon, gtf_exon$gene_id)), use.names = FALSE)
gene_ranges$gene_id <- names(split(gtf_exon, gtf_exon$gene_id))
tf2_gene_ranges <- gene_ranges[gene_ranges$gene_id %in% tf2_genes]
length(tf2_gene_ranges)

# Intact LTR retrotransposons identified in the genome with LTR_retriever (run beforehand on GCA_037325665.1, default settings).
# repeat_region = whole intact element (from TSD to TSD); long_terminal_repeat = its two LTRs.
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
# Sequences extracted from the genome with bedtools getfasta -> TF2_intact_LTR_80percent_overlap.fa


############################################################
#Part 2: Tree of the 47 elements, clade A and clade B, first consensus (Figure S7)

# Elements on the minus strand are reverse-complemented so that all sequences have the same orientation.
seed_seq <- readDNAStringSet("TF2_intact_LTR_80percent_overlap.fa")
names(seed_seq) <- sub("::.*", "", names(seed_seq))
seed_strand <- repeat_region$strand[match(names(seed_seq), repeat_region$intact_id)]
table(seed_strand)
seed_seq[seed_strand == "-"] <- reverseComplement(seed_seq[seed_strand == "-"])

# Neighbour-joining tree from the multiple alignment, rooted at the midpoint. The nodes are only rotated (same topology and
# branch lengths) so that the long elements (> 7 kb) are drawn at the top.
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

# Figure S7
clade_colors <- c(cladeA_cons = "#D55E00", cladeB_cons = "#0072B2")
tip_colors <- ifelse(tree_rot$tip.label %in% seed_clade_A, clade_colors[["cladeA_cons"]], ifelse(tree_rot$tip.label %in% seed_clade_B, clade_colors[["cladeB_cons"]], "grey40"))
cairo_pdf("FigureS7.pdf", width = 11, height = 13)
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
#Part 3: All intact elements of clade A and clade B in the genome

# All intact elements of the genome.
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

# One element usually has several hits on the same consensus (e.g. one for each LTR); overlapping hits are merged before the
# coverage is calculated. strand = orientation of the element relative to the consensus.
merged_length <- function(s, e) sum(width(IRanges::reduce(IRanges(pmin(s, e), pmax(s, e)))))
blast1_summary <- blast1 %>% dplyr::group_by(qseqid, sseqid) %>%
  dplyr::summarise(aln_len = sum(length), identity = sum(pident * length) / sum(length), cov_intact = merged_length(qstart, qend) / dplyr::first(qlen) * 100,
                   strand = ifelse(sum(length[send > sstart]) >= sum(length) / 2, "+", "-"), bitscore = sum(bitscore), .groups = "drop")

# 80-80-80 rule (identity >= 80%, alignment >= 80 bp, >= 80% of the element covered); an element matching both consensus
# sequences is assigned to the one with the highest score. Clade A: 138 copies, clade B: 64 copies.
hits <- blast1_summary %>% dplyr::filter(identity >= 80, aln_len >= 80, cov_intact >= 80) %>% dplyr::group_by(qseqid) %>% dplyr::slice_max(bitscore, n = 1, with_ties = FALSE) %>% dplyr::ungroup()
table(hits$sseqid)

# Genomic strand on which each copy has the same orientation as its consensus (the sequences were extracted with "?" read as +).
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

# Second consensus from all copies: columns covered by fewer than 30% of the copies are dropped (gaps are not counted as a base),
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

# Run on Linux with the second consensus sequences:
# - TEsorter consensus_cladeA_cladeB.fa -db gydb (and -db rexdb-metazoa); protein domain positions (GYDB) -> domain.txt
# - bedtools getfasta -s -nameOnly with cladeA_intact.bed / cladeB_intact.bed -> cladeA_intact.fa / cladeB_intact.fa
# - BLASTN of the copies of each clade against the consensus -> intactLTR_vs_consensus.blastn.tsv
# - mafft --add cladeA_intact.fa --keeplength cladeA_cons.fa > cladeA_to_cons.aln.fa (same for clade B)
# - RNA-seq reads and pool-seq reads mapped to the two consensus sequences plus two housekeeping genes (g15582.t1, g13887.t1):
#   per-base pool-seq depth -> all_depth.txt; RNA-seq read counts in the RT domain -> tf2_RT_counts.txt; mapped fragments -> library_size.txt


############################################################
#Part 4: Structure of the two consensus sequences

clade_order <- c("cladeA_cons", "cladeB_cons")
cons <- readDNAStringSet("consensus_cladeA_cladeB.fa")
names(cons) <- sub("[#[:space:]].*$", "", names(cons))
cons <- cons[clade_order]
cons_len <- setNames(width(cons), names(cons))

# 5′ and 3′ LTR of every copy (LTR_retriever): on the plus strand the left LTR is the 5′ LTR, on the minus strand the right LTR.
ltr_lengths <- members %>% dplyr::select(qseqid, sseqid, strand_cons) %>% dplyr::left_join(repeat_region %>% dplyr::select(qseqid = intact_id, parent), by = "qseqid") %>%
  dplyr::left_join(ltr_region, by = "parent") %>%
  dplyr::mutate(label = dplyr::if_else((strand_cons == "+" & grepl("^lLTR_", ltr_id)) | (strand_cons == "-" & grepl("^rLTR_", ltr_id)), "5′ LTR", "3′ LTR"),
                LTR_length = ltr_end - ltr_start + 1)

# LTRs drawn on the consensus with the median length of each clade; protein domains from TEsorter.
ltr_cons <- ltr_lengths %>% dplyr::group_by(sseqid, label) %>% dplyr::summarise(LTR_length = round(median(LTR_length)), .groups = "drop") %>%
  dplyr::mutate(start = dplyr::if_else(label == "5′ LTR", 1, as.numeric(cons_len[sseqid]) - LTR_length + 1), end = dplyr::if_else(label == "5′ LTR", LTR_length, as.numeric(cons_len[sseqid])), feature = "LTR")
ltr_cons
domain <- read.table("domain.txt", header = FALSE, col.names = c("sseqid", "start", "end", "feature"))
features <- dplyr::bind_rows(ltr_cons %>% dplyr::select(sseqid, start, end, feature, label), domain %>% dplyr::mutate(label = feature)) %>%
  dplyr::mutate(mid = (start + end) / 2, sseqid = factor(sseqid, levels = clade_order))


############################################################
#Part 5: HSF1 binding motifs (Figure S11)

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

# Run on Linux: FIMO with the HSF1 motif (JASPAR MA0486.1), p < 1e-4.
# fimo --o fimo_cladeA_5LTR --thresh 1e-4 MA0486.1.meme cladeA_5LTR.fa   -> cladeA_5LTR_fimo.tsv
# fimo --o fimo_cladeB_5LTR --thresh 1e-4 MA0486.1.meme cladeB_5LTR.fa   -> cladeB_5LTR_fimo.tsv
# Positive control: 2 kb upstream of the start codon of the eight Hsp70 genes (8HSP70_up2kb.fa) -> 8HSP70_up2kb_fimo.tsv
fimo_A <- read.delim("cladeA_5LTR_fimo.tsv", comment.char = "#", quote = "")
fimo_B <- read.delim("cladeB_5LTR_fimo.tsv", comment.char = "#", quote = "")
fimo_HSP70 <- read.delim("8HSP70_up2kb_fimo.tsv", comment.char = "#", quote = "")
seq_HSP70 <- readDNAStringSet("8HSP70_up2kb.fa")

# One line per sequence (copies labelled with their genomic position).
line_A <- data.frame(sequence_name = ltr5_A$fasta_id, LTR_length = ltr5_A$ltr_length, display_id = ltr5_A$qseqid) %>% dplyr::arrange(sequence_name) %>% dplyr::mutate(ypos = rev(seq_len(dplyr::n())))
line_B <- data.frame(sequence_name = ltr5_B$fasta_id, LTR_length = ltr5_B$ltr_length, display_id = ltr5_B$qseqid) %>% dplyr::arrange(sequence_name) %>% dplyr::mutate(ypos = rev(seq_len(dplyr::n())))
line_HSP70 <- data.frame(sequence_name = names(seq_HSP70), upstream_length = width(seq_HSP70)) %>% dplyr::arrange(sequence_name) %>%
  dplyr::mutate(ypos = rev(seq_len(dplyr::n())), display_id = sequence_name)

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

# Number of copies with at least one HSF1 site.
c(clade_A = dplyr::n_distinct(point_A$sequence_name), clade_B = dplyr::n_distinct(point_B$sequence_name))

# The most frequent motif sequence in each clade (in consensus orientation) is located in the 5′ LTR of the consensus;
# this position is marked in Figure 5A.
dominant_A <- point_A %>% dplyr::mutate(cons_sequence = dplyr::if_else(strand == "+", matched_sequence, as.character(reverseComplement(DNAStringSet(matched_sequence))))) %>%
  dplyr::count(cons_sequence, sort = TRUE) %>% dplyr::slice_head(n = 1)
dominant_B <- point_B %>% dplyr::mutate(cons_sequence = dplyr::if_else(strand == "+", matched_sequence, as.character(reverseComplement(DNAStringSet(matched_sequence))))) %>%
  dplyr::count(cons_sequence, sort = TRUE) %>% dplyr::slice_head(n = 1)
hsf_A <- matchPattern(dominant_A$cons_sequence, cons[["cladeA_cons"]])
hsf_B <- matchPattern(dominant_B$cons_sequence, cons[["cladeB_cons"]])
hsf_A <- hsf_A[end(hsf_A) <= ltr_cons$end[ltr_cons$sseqid == "cladeA_cons" & ltr_cons$label == "5′ LTR"]]
hsf_B <- hsf_B[end(hsf_B) <= ltr_cons$end[ltr_cons$sseqid == "cladeB_cons" & ltr_cons$label == "5′ LTR"]]
hsf_sites <- data.frame(sseqid = factor(clade_order, levels = clade_order), mid = c((start(hsf_A) + end(hsf_A)) / 2, (start(hsf_B) + end(hsf_B)) / 2))
hsf_sites

# Figure S11
theme_fimo <- theme_minimal(base_size = 13, base_family = "Arial") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 14, margin = margin(t = 7)),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), plot.tag = element_text(size = 18, face = "bold"), plot.tag.position = "topleft",
        plot.margin = margin(8, 10, 8, 8))

pS11A <- ggplot() + geom_segment(data = line_HSP70, aes(x = -upstream_length, xend = -1, y = ypos, yend = ypos), colour = "grey80", linewidth = 0.5) +
  geom_point(data = point_HSP70, aes(x = mid_relative, y = ypos), shape = 17, size = 2.1, colour = "firebrick") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  scale_x_continuous(limits = c(-max(line_HSP70$upstream_length), 0), breaks = seq(-max(line_HSP70$upstream_length), 0, 500), expand = expansion(add = c(0, 30))) +
  scale_y_continuous(breaks = line_HSP70$ypos, labels = line_HSP70$display_id, expand = expansion(add = 0.5)) +
  labs(x = "HSF1 motif position relative to start codon (bp)", y = NULL, title = "Hsp70", tag = "A") + theme_fimo +
  theme(axis.text.y = element_text(size = 12, colour = "black", margin = margin(r = 4)))
pS11B <- ggplot() + geom_segment(data = line_A, aes(x = 1, xend = LTR_length, y = ypos, yend = ypos), colour = "grey80", linewidth = 0.35) +
  geom_point(data = point_A, aes(x = mid, y = ypos), shape = 17, size = 2.1, colour = "firebrick") +
  scale_x_continuous(limits = c(0, max(line_A$LTR_length)), breaks = seq(0, max(line_A$LTR_length), 100), expand = expansion(add = c(0, 5))) +
  scale_y_continuous(breaks = line_A$ypos, labels = line_A$display_id, expand = expansion(add = 0.5)) +
  labs(x = "Position within 5′ LTR (bp)", y = NULL, title = "Clade A", tag = "B") + theme_fimo +
  theme(axis.text.y = element_text(size = 8, colour = "black", margin = margin(r = 4)), plot.title = element_text(colour = clade_colors[["cladeA_cons"]]))
pS11C <- ggplot() + geom_segment(data = line_B, aes(x = 1, xend = LTR_length, y = ypos, yend = ypos), colour = "grey80", linewidth = 0.35) +
  geom_point(data = point_B, aes(x = mid, y = ypos), shape = 17, size = 2.1, colour = "firebrick") +
  scale_x_continuous(limits = c(0, max(line_B$LTR_length)), breaks = seq(0, max(line_B$LTR_length), 100), expand = expansion(add = c(0, 5))) +
  scale_y_continuous(breaks = line_B$ypos, labels = line_B$display_id, expand = expansion(add = 0.5)) +
  labs(x = "Position within 5′ LTR (bp)", y = NULL, title = "Clade B", tag = "C") + theme_fimo +
  theme(axis.text.y = element_text(size = 8, colour = "black", margin = margin(r = 4)), plot.title = element_text(colour = clade_colors[["cladeB_cons"]]))

figS11 <- pS11A / (pS11B | pS11C) + plot_layout(heights = c(2, 18))
figS11
ggsave("FigureS11.pdf", figS11, width = 16, height = 24, device = cairo_pdf)


############################################################
#Part 6: Coverage of the consensus by the intact copies (Figure 5A) and alignment of each copy (Figure S8)

# BLASTN hits of each copy on the consensus of its own clade.
blast2 <- read.table("intactLTR_vs_consensus.blastn.tsv", header = FALSE, col.names = c("qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"))
blast2 <- blast2 %>% dplyr::mutate(qseqid = sub("\\([+-]\\)$", "", qseqid), sseqid = sub("[#[:space:]].*$", "", sseqid), cs = pmin(sstart, send), ce = pmax(sstart, send)) %>%
  dplyr::semi_join(members, by = c("qseqid", "sseqid"))

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
clade_labels <- c(cladeA_cons = "Clade A", cladeB_cons = "Clade B")
feature_colors <- c(LTR = "#4D4D4D", GAG = "#4A9E6B", PR = "#E9A23B", INT = "#C9484F", RT = "#4A72B0", RNaseH = "#8A74B5")
p5A <- ggplot() + geom_line(data = coverage_data, aes(x = pos, y = prop, colour = sseqid), linewidth = 0.8) +
  geom_segment(data = data.frame(sseqid = factor(clade_order, levels = clade_order), end = as.numeric(cons_len)), aes(x = 1, xend = end, y = 0.25, yend = 0.25), colour = "grey40", linewidth = 0.6) +
  geom_rect(data = features, aes(xmin = start, xmax = end, ymin = 0.16, ymax = 0.34, fill = feature), colour = "black", linewidth = 0.35) +
  geom_text(data = features, aes(x = mid, y = 0.25, label = label, size = dplyr::if_else(feature == "RNaseH", 3, 3.4)), colour = "white", fontface = "bold") +
  geom_segment(data = hsf_sites, aes(x = mid, xend = mid, y = 0.14, yend = 0.21), colour = "firebrick", linewidth = 0.9) +
  geom_segment(data = hsf_sites, aes(x = mid, xend = mid, y = 0.29, yend = 0.36), colour = "firebrick", linewidth = 0.9) +
  geom_segment(data = hsf_sites, aes(x = mid + 350, xend = mid, y = 0.54, yend = 0.37), colour = "black", linewidth = 0.5, arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")) +
  geom_text(data = hsf_sites, aes(x = mid + 350, y = 0.61), label = "Predicted HSF1 site", hjust = 0, size = 3.6) +
  scale_colour_manual(values = clade_colors, guide = "none") + scale_fill_manual(values = feature_colors, guide = "none") + scale_size_identity() +
  scale_x_continuous(breaks = seq(0, max(cons_len), 1000), expand = expansion(mult = c(0.005, 0.005))) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2), labels = scales::label_percent(), expand = expansion(mult = c(0, 0.02))) +
  facet_wrap(~sseqid, ncol = 1, scales = "free_x", labeller = as_labeller(clade_labels)) +
  labs(x = "Position on consensus (bp)", y = "Intact copies covering position (%)", tag = "A") + theme_TF2 + theme(panel.grid.major = element_blank())

# Figure S8: MAFFT alignment of every copy to its consensus; each position is a match, a mismatch (SNP) or a gap.
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
alignment_A$type <- as.vector(ifelse(mat_A == "-", "gap", ifelse(mat_A == matrix(aln_A["cladeA_cons", keep_A], nrow(mat_A), length(keep_A), byrow = TRUE), "match", "mismatch")))
alignment_B$type <- as.vector(ifelse(mat_B == "-", "gap", ifelse(mat_B == matrix(aln_B["cladeB_cons", keep_B], nrow(mat_B), length(keep_B), byrow = TRUE), "match", "mismatch")))

pS8_aln_A <- ggplot(alignment_A, aes(x = pos, y = copy, fill = type)) + geom_raster() +
  scale_fill_manual(values = c(match = "grey70", mismatch = "#E41A1C", gap = "white"), breaks = c("match", "mismatch"), labels = c("Match", "Mismatch")) +
  scale_x_continuous(limits = c(0, cons_len[["cladeA_cons"]]), expand = c(0, 0)) + scale_y_reverse(expand = c(0, 0)) +
  labs(x = NULL, y = "Intact copies", fill = NULL, title = paste0("Clade A (n = ", nrow(mat_A), ")"), tag = "A") + theme_TF2 +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.title = element_text(colour = clade_colors[["cladeA_cons"]]))
pS8_aln_B <- ggplot(alignment_B, aes(x = pos, y = copy, fill = type)) + geom_raster() +
  scale_fill_manual(values = c(match = "grey70", mismatch = "#E41A1C", gap = "white"), breaks = c("match", "mismatch"), labels = c("Match", "Mismatch")) +
  scale_x_continuous(limits = c(0, cons_len[["cladeB_cons"]]), expand = c(0, 0)) + scale_y_reverse(expand = c(0, 0)) +
  labs(x = NULL, y = "Intact copies", fill = NULL, title = paste0("Clade B (n = ", nrow(mat_B), ")"), tag = "B") + theme_TF2 +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.title = element_text(colour = clade_colors[["cladeB_cons"]]))
pS8_track_A <- ggplot(features %>% dplyr::filter(sseqid == "cladeA_cons")) +
  annotate("segment", x = 1, xend = cons_len[["cladeA_cons"]], y = 0.5, yend = 0.5, colour = "grey40", linewidth = 0.4) +
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = feature), colour = "black", linewidth = 0.3) +
  geom_text(aes(x = mid, y = 0.5, label = label), colour = "white", fontface = "bold", size = 3.8) + scale_fill_manual(values = feature_colors, guide = "none") +
  scale_x_continuous(limits = c(0, cons_len[["cladeA_cons"]]), breaks = seq(0, cons_len[["cladeA_cons"]], 1000), expand = c(0, 0)) +
  coord_cartesian(ylim = c(-0.1, 1.1)) + labs(x = NULL, y = NULL) + theme_TF2 +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank())
pS8_track_B <- ggplot(features %>% dplyr::filter(sseqid == "cladeB_cons")) +
  annotate("segment", x = 1, xend = cons_len[["cladeB_cons"]], y = 0.5, yend = 0.5, colour = "grey40", linewidth = 0.4) +
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = feature), colour = "black", linewidth = 0.3) +
  geom_text(aes(x = mid, y = 0.5, label = label), colour = "white", fontface = "bold", size = 3.8) + scale_fill_manual(values = feature_colors, guide = "none") +
  scale_x_continuous(limits = c(0, cons_len[["cladeB_cons"]]), breaks = seq(0, cons_len[["cladeB_cons"]], 1000), expand = c(0, 0)) +
  coord_cartesian(ylim = c(-0.1, 1.1)) + labs(x = "Position on consensus (bp)", y = NULL) + theme_TF2 +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank())

# Legend in its own row on top; the alignment panels are scaled to the number of copies.
figS8 <- wrap_plots(guide_area(), pS8_aln_A, pS8_track_A, pS8_aln_B, pS8_track_B, ncol = 1, heights = c(0.45, 6, 1, 6 * nrow(mat_B) / nrow(mat_A), 1), guides = "collect") &
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))
figS8
ggsave("FigureS8.pdf", figS8, width = 13, height = 11, device = cairo_pdf)


############################################################
#Part 7: Insertion time (Figure 5B, Figure S9)

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
insertion_test$label <- ifelse(insertion_test$p_adj < 0.0001, "****", ifelse(insertion_test$p_adj < 0.001, "***", ifelse(insertion_test$p_adj < 0.01, "**", ifelse(insertion_test$p_adj < 0.05, "*", "ns"))))
insertion_test

# Significant comparisons only; shorter brackets below the longer ones.
insertion_sig <- insertion_test %>% dplyr::filter(p_adj < 0.05) %>% dplyr::arrange(xmax - xmin) %>% dplyr::mutate(y = 1.55 + (dplyr::row_number() - 1) * 0.17)

p5B <- ggplot(insertion_data, aes(x = comparison_group, y = insertion_time_Mya, fill = comparison_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85, linewidth = 0.4) +
  geom_point(data = insertion_data %>% dplyr::filter(comparison_group != "Other intact Gypsy LTRs"), aes(colour = comparison_group), position = position_jitter(width = 0.12, height = 0, seed = 123), size = 1.4, alpha = 0.4) +
  scale_fill_manual(values = comparison_colors, guide = "none") + scale_colour_manual(values = comparison_colors, guide = "none") +
  scale_x_discrete(labels = c("Other intact Gypsy LTRs" = "Other intact\nGypsy LTRs", "Clade A" = "Clade A", "Clade B" = "Clade B")) +
  geom_segment(data = insertion_sig, aes(x = xmin, xend = xmax, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = insertion_sig, aes(x = xmin, xend = xmin, y = y, yend = y - 0.035), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = insertion_sig, aes(x = xmax, xend = xmax, y = y, yend = y - 0.035), inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = insertion_sig, aes(x = (xmin + xmax) / 2, y = y + 0.035, label = label), inherit.aes = FALSE, size = 4) +
  coord_cartesian(ylim = c(0, 1.85)) + labs(x = NULL, y = "Estimated insertion time (Mya)", tag = "B") + theme_TF2 +
  theme(axis.text.x = element_text(size = 12, face = "bold"), panel.grid.major = element_blank())

# Figure S9: distribution of insertion times (bins of 0.05 Mya, relative frequency within each group).
hist_data <- insertion_data %>% dplyr::mutate(bin_mid = floor(insertion_time_Mya / 0.05) * 0.05 + 0.05 / 2) %>% dplyr::count(comparison_group, bin_mid) %>%
  dplyr::group_by(comparison_group) %>% dplyr::mutate(relative_frequency = n / sum(n)) %>% dplyr::ungroup()
windowsFonts(Arial = windowsFont("Arial"))
figS9 <- ggplot(hist_data, aes(x = bin_mid, y = relative_frequency, fill = comparison_group)) + geom_col(width = 0.05, colour = "black", linewidth = 0.3) +
  facet_wrap(~comparison_group, ncol = 1) + scale_fill_manual(values = comparison_colors, guide = "none") +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), expand = expansion(mult = c(0, 0.02))) + coord_cartesian(ylim = c(0, 0.8)) +
  labs(x = "Estimated insertion time (Mya)", y = "Relative frequency") + theme_TF2
figS9
ggsave("FigureS9.pdf", figS9, width = 8, height = 10, device = cairo_pdf)


############################################################
#Part 8: Relative copy number from pool-seq depth (Figure S10)

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
# Per population: mean of the high- and low-shore pools.
copy_population <- copy_pool %>% dplyr::mutate(population = factor(sub("(high|low)$", "", pool), levels = population_order)) %>%
  dplyr::group_by(Chr, population) %>% dplyr::summarise(relative_copy_number = mean(relative_copy_number), .groups = "drop") %>% dplyr::mutate(Chr = factor(Chr, levels = clade_order))
copy_population
windowsFonts(Arial = windowsFont("Arial"))
figS10 <- ggplot(copy_population, aes(x = population, y = relative_copy_number, fill = Chr)) + geom_col(width = 0.65, colour = "black", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.1f", relative_copy_number)), vjust = -0.4, size = 4) + scale_fill_manual(values = clade_colors, guide = "none") +
  facet_wrap(~Chr, nrow = 1, labeller = as_labeller(clade_labels)) + scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Relative copy number") + theme_TF2
figS10
ggsave("FigureS10.pdf", figS10, width = 9, height = 5.5, device = cairo_pdf)


############################################################
#Part 9: RT expression in each individual (Figure 5C, 5D)

# RNA-seq read counts on the RT domain of each consensus and number of mapped fragments of each individual.
rt_counts <- read.delim("tf2_RT_counts.txt", check.names = FALSE)
library_size <- read.delim("library_size.txt") %>% dplyr::filter(sample != "SPL47")
group <- read.table("../01counts_ind_filter/group.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_ids <- names(rt_counts)[-(1:6)]

stress_order <- c("No/short heat stress", "Long heat stress")
# 31°C for 3 hours = long heat stress; the other three treatments = no or short heat stress.
group <- group %>% dplyr::mutate(exposure_group = factor(dplyr::if_else(temperature == "31degree" & time == "3h", "Long heat stress", "No/short heat stress"), levels = stress_order),
                                 pool = paste0(population, shore_level), population = factor(population, levels = population_order))

# FPKM = reads per million mapped fragments per kb of RT; FPKM / relative copy number = expression per copy, using the copy
# number of the pool (population x shore level) the individual comes from.
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
population_tests$label_FPKM <- as.character(cut(population_tests$FDR_FPKM, c(0, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("****", "***", "**", "*", "ns"), right = FALSE))
population_tests$label_copy <- as.character(cut(population_tests$FDR_copy, c(0, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("****", "***", "**", "*", "ns"), right = FALSE))
population_tests

# Significance brackets above the boxes of the two populations compared.
test_steps <- expression_data %>% dplyr::group_by(Chr) %>%
  dplyr::summarise(step_FPKM = max(diff(range(log2(FPKM))), 0.1) * 0.1, step_copy = max(diff(range(log2(FPKM_per_copy))), 0.1) * 0.1, .groups = "drop")
test_tops <- expression_data %>% dplyr::group_by(Chr, exposure_group) %>% dplyr::summarise(top_FPKM = max(log2(FPKM)), top_copy = max(log2(FPKM_per_copy)), .groups = "drop")
population_brackets <- population_tests %>% dplyr::mutate(Chr = factor(Chr, levels = clade_order), exposure_group = factor(exposure_group, levels = stress_order)) %>%
  dplyr::left_join(test_steps, by = "Chr") %>% dplyr::left_join(test_tops, by = c("Chr", "exposure_group")) %>%
  dplyr::mutate(x = match(exposure_group, stress_order), x_left = x + (match(population1, population_order) - 2) * 0.8 / 3 + 0.02,
                x_right = x + (match(population2, population_order) - 2) * 0.8 / 3 - 0.02, y_FPKM = top_FPKM + step_FPKM, y_copy = top_copy + step_copy)
sig_FPKM <- population_brackets %>% dplyr::filter(FDR_FPKM < 0.05)
sig_copy <- population_brackets %>% dplyr::filter(FDR_copy < 0.05)

population_colors <- c(North = "#87aade", Central = "#ffe57f", South = "#ff9999")
p5C <- ggplot(expression_data, aes(x = exposure_group, y = log2(FPKM), fill = population)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.65, colour = "black", linewidth = 0.4, outlier.shape = NA) +
  geom_point(aes(colour = population), size = 1.6, alpha = 0.65, position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8, seed = 123)) +
  geom_segment(data = sig_FPKM, aes(x = x_left, xend = x_right, y = y_FPKM, yend = y_FPKM), inherit.aes = FALSE, linewidth = 0.35) +
  geom_segment(data = sig_FPKM, aes(x = x_left, xend = x_left, y = y_FPKM, yend = y_FPKM - step_FPKM * 0.15), inherit.aes = FALSE, linewidth = 0.35) +
  geom_segment(data = sig_FPKM, aes(x = x_right, xend = x_right, y = y_FPKM, yend = y_FPKM - step_FPKM * 0.15), inherit.aes = FALSE, linewidth = 0.35) +
  geom_text(data = sig_FPKM, aes(x = (x_left + x_right) / 2, y = y_FPKM + step_FPKM * 0.22, label = label_FPKM), inherit.aes = FALSE, size = 4) +
  scale_fill_manual(values = population_colors) + scale_colour_manual(values = population_colors, guide = "none") + scale_x_discrete(expand = expansion(add = 0.5)) +
  facet_wrap(~Chr, nrow = 1, labeller = as_labeller(clade_labels)) + scale_y_continuous(expand = expansion(mult = c(0.04, 0.15))) +
  labs(x = NULL, y = expression(log[2]*"(FPKM)"), tag = "C") + theme_TF2 + theme(axis.text.x = element_text(face = "bold"))

p5D <- ggplot(expression_data, aes(x = exposure_group, y = log2(FPKM_per_copy), fill = population)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.65, colour = "black", linewidth = 0.4, outlier.shape = NA) +
  geom_point(aes(colour = population), size = 1.6, alpha = 0.65, position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8, seed = 123)) +
  geom_segment(data = sig_copy, aes(x = x_left, xend = x_right, y = y_copy, yend = y_copy), inherit.aes = FALSE, linewidth = 0.35) +
  geom_segment(data = sig_copy, aes(x = x_left, xend = x_left, y = y_copy, yend = y_copy - step_copy * 0.15), inherit.aes = FALSE, linewidth = 0.35) +
  geom_segment(data = sig_copy, aes(x = x_right, xend = x_right, y = y_copy, yend = y_copy - step_copy * 0.15), inherit.aes = FALSE, linewidth = 0.35) +
  geom_text(data = sig_copy, aes(x = (x_left + x_right) / 2, y = y_copy + step_copy * 0.22, label = label_copy), inherit.aes = FALSE, size = 4) +
  scale_fill_manual(values = population_colors) + scale_colour_manual(values = population_colors, guide = "none") + scale_x_discrete(expand = expansion(add = 0.5)) +
  facet_wrap(~Chr, nrow = 1, labeller = as_labeller(clade_labels)) + scale_y_continuous(expand = expansion(mult = c(0.04, 0.15))) +
  labs(x = NULL, y = expression(log[2]*"(FPKM / relative copy number)"), tag = "D") + theme_TF2 + theme(axis.text.x = element_text(face = "bold"))


############################################################
#Figure 5

figure5 <- ((p5A | p5B) + plot_layout(widths = c(2, 1))) / ((free(p5C, type = "label", side = "l") | p5D) + plot_layout(widths = c(1, 1))) +
  plot_layout(heights = c(1.15, 1), guides = "collect") & theme(legend.position = "bottom")
figure5
ggsave("Figure5.pdf", figure5, width = 16, height = 11, device = cairo_pdf)

