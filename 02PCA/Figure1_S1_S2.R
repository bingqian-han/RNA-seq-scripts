setwd("/data/02PCA/")

library(DESeq2)
library(scatterplot3d)
library(ggplot2)
library(patchwork)
library(vegan)

# Load the count matrix and sample information.
counts <- read.table("../01counts_ind_filter/counts_filter.txt", header = TRUE)
group <- read.table("../01counts_ind_filter/group.txt", header = TRUE, sep = "\t")
rownames(group) <- group$id
group <- group[colnames(counts), ]

# PCA of all retained genes.
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = group, design = ~1)
vsd <- vst(dds, blind = TRUE)
head(assay(vsd), 3)
pca <- prcomp(t(assay(vsd)), center = TRUE, scale. = FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)

# Combine PC1-PC3 scores with population and assign colours.
pcaData <- data.frame(pca$x[, 1:3], population = group$population)
pop_colors <- c("North" = "#87aade", "Central" = "#ffe57f", "South" = "#ff9999")  # blue, yellow, pink
pcaData$color <- pop_colors[pcaData$population]

# 3D PCA plot.
pdf("Figure1B.pdf", width = 5, height = 5.3)
par(cex.lab = 0.92, cex.axis = 0.95)
scatterplot3d(pcaData$PC1, pcaData$PC2, pcaData$PC3, color = pcaData$color, pch = 16, main = "", angle = 135,
              xlab = paste0("PC1 (", percentVar[1], "%)"), ylab = paste0("PC2 (", percentVar[2], "%)"), 
              zlab = paste0("PC3 (", percentVar[3], "%)"))
dev.off()


# Combine PC1-PC3 scores with sample information and assign population colours.
pcaData <- data.frame(pca$x[, 1:3], group)

# Panel S1A: PC1 vs PC2.
pS1A <- ggplot(pcaData, aes(PC1, PC2, colour = population, fill = temperature, shape = shore_level, size = time)) + geom_point(stroke = 1) +
  labs(title = "All samples", x = paste0("PC1 (", percentVar[1], "%)"), y = paste0("PC2 (", percentVar[2], "%)"))

# Panel S1B: PC1 vs PC3.
pS1B <- ggplot(pcaData, aes(PC1, PC3, colour = population, fill = temperature, shape = shore_level, size = time)) + geom_point(stroke = 1) +
  labs(title = "All samples", x = paste0("PC1 (", percentVar[1], "%)"), y = paste0("PC3 (", percentVar[3], "%)"))

# Combine panel S1A with S1B.
pS1 <- pS1A + pS1B + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
  scale_colour_manual(values = pop_colors, breaks = c("North", "Central", "South"), name = "Sampling location") &
  scale_fill_manual(values = c("17degree" = "white", "31degree" = "black"), labels = c("17\u00B0C", "31\u00B0C"), name = "Temperature") &
  scale_shape_manual(values = c("high" = 24, "low" = 21), labels = c("High shore", "Low shore"), name = "Shore level") &
  scale_size_manual(values = c("15min" = 2.8, "3h" = 4.5), labels = c("15 mins", "3 hours"), name = "Time") &
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 21, fill = "grey50", size = 4)),
         shape = guide_legend(order = 2, override.aes = list(fill = "grey50", colour = "black", size = 4)),
         size = guide_legend(order = 3, override.aes = list(shape = 21, fill = "grey50", colour = "black")),
         fill = guide_legend(order = 4, override.aes = list(shape = 21, colour = "black", size = 4, stroke = 1.2))) &
  theme_bw() & theme(panel.grid = element_blank(), plot.title = element_text(size = 16, face = "bold", hjust = 0.5), plot.tag = element_text(size = 14, face = "bold"),
                     axis.title = element_text(size = 11), axis.text = element_text(size = 9), legend.title = element_text(size = 14), legend.text = element_text(size = 11))
pS1

ggsave("FigureS1.pdf", pS1, width = 14, height = 6)


# PCA within each population.
vsd_North <- vst(dds[, group$population == "North"], blind = TRUE)
vsd_Central <- vst(dds[, group$population == "Central"], blind = TRUE)
vsd_South <- vst(dds[, group$population == "South"], blind = TRUE)
pca_North <- prcomp(t(assay(vsd_North)), center = TRUE, scale. = FALSE)
pca_Central <- prcomp(t(assay(vsd_Central)), center = TRUE, scale. = FALSE)
pca_South <- prcomp(t(assay(vsd_South)), center = TRUE, scale. = FALSE)
pv_North <- sprintf("%.2f", 100 * pca_North$sdev^2 / sum(pca_North$sdev^2))
pv_Central <- sprintf("%.2f", 100 * pca_Central$sdev^2 / sum(pca_Central$sdev^2))
pv_South <- sprintf("%.2f", 100 * pca_South$sdev^2 / sum(pca_South$sdev^2))

# Combine PC1-PC3 scores with sample information.
pcaData_North <- data.frame(pca_North$x[, 1:3], group[group$population == "North", ])
pcaData_Central <- data.frame(pca_Central$x[, 1:3], group[group$population == "Central", ])
pcaData_South <- data.frame(pca_South$x[, 1:3], group[group$population == "South", ])

# Panels A-C: PC1 vs PC2. Colour = temperature, shape = shore level, size = time.
pS2A <- ggplot(pcaData_North, aes(PC1, PC2, colour = temperature, shape = shore_level, size = time)) + geom_point() +
  labs(title = "North", x = paste0("PC1 (", pv_North[1], "%)"), y = paste0("PC2 (", pv_North[2], "%)"))
pS2B <- ggplot(pcaData_Central, aes(PC1, PC2, colour = temperature, shape = shore_level, size = time)) + geom_point() +
  labs(title = "Central", x = paste0("PC1 (", pv_Central[1], "%)"), y = paste0("PC2 (", pv_Central[2], "%)"))
pS2C <- ggplot(pcaData_South, aes(PC1, PC2, colour = temperature, shape = shore_level, size = time)) + geom_point() +
  labs(title = "South", x = paste0("PC1 (", pv_South[1], "%)"), y = paste0("PC2 (", pv_South[2], "%)"))

# Panels D-F: PC1 vs PC3.
pS2D <- ggplot(pcaData_North, aes(PC1, PC3, colour = temperature, shape = shore_level, size = time)) + geom_point() +
  labs(title = "North", x = paste0("PC1 (", pv_North[1], "%)"), y = paste0("PC3 (", pv_North[3], "%)"))
pS2E <- ggplot(pcaData_Central, aes(PC1, PC3, colour = temperature, shape = shore_level, size = time)) + geom_point() +
  labs(title = "Central", x = paste0("PC1 (", pv_Central[1], "%)"), y = paste0("PC3 (", pv_Central[3], "%)"))
pS2F <- ggplot(pcaData_South, aes(PC1, PC3, colour = temperature, shape = shore_level, size = time)) + geom_point() +
  labs(title = "South", x = paste0("PC1 (", pv_South[1], "%)"), y = paste0("PC3 (", pv_South[3], "%)"))

# Combine the six panels with a shared legend; everything after "&" applies to all panels.
p_figureS2 <- pS2A + pS2B + pS2C + pS2D + pS2E + pS2F + plot_layout(ncol = 3, guides = "collect") + plot_annotation(tag_levels = "A") &
  scale_colour_manual(values = c("17degree" = "#00A1D5", "31degree" = "#DC5E83"), labels = c("17\u00B0C", "31\u00B0C"), name = "Temperature") &
  scale_shape_manual(values = c("high" = 17, "low" = 19), labels = c("High shore", "Low shore"), name = "Shore level") &
  scale_size_manual(values = c("15min" = 2, "3h" = 4), labels = c("15 mins", "3 hours"), name = "Time") &
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2), size = guide_legend(order = 3)) &
  theme(plot.title = element_text(hjust = 0.5, vjust = 2, face = "bold", size = 16), panel.grid = element_blank(),
        panel.border = element_rect(colour = "#969696", fill = NA, linewidth = 0.8), panel.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(size = 14, face = "bold"), legend.text = element_text(size = 11), legend.title = element_text(size = 14))
p_figureS2

# Save the figure.
ggsave("FigureS2.pdf", p_figureS2, width = 14, height = 8.5)


# Define exposure groups: long exposure = 31°C for 3 h; all other treatments = no/short exposure.
pcaData_North$exposure <- ifelse(pcaData_North$temperature == "31degree" & pcaData_North$time == "3h", "Long exposure", "No/short exposure")
pcaData_Central$exposure <- ifelse(pcaData_Central$temperature == "31degree" & pcaData_Central$time == "3h", "Long exposure", "No/short exposure")
pcaData_South$exposure <- ifelse(pcaData_South$temperature == "31degree" & pcaData_South$time == "3h", "Long exposure", "No/short exposure")

# Panels C-E: PC1 vs PC2 within each population. Colour = exposure group, shape = shore level.
pMainC <- ggplot(pcaData_North, aes(PC1, PC2, colour = exposure, shape = shore_level)) + geom_point(size = 2) +
  labs(title = "North", x = paste0("PC1 (", pv_North[1], "%)"), y = paste0("PC2 (", pv_North[2], "%)"))
pMainD <- ggplot(pcaData_Central, aes(PC1, PC2, colour = exposure, shape = shore_level)) + geom_point(size = 2) +
  labs(title = "Central", x = paste0("PC1 (", pv_Central[1], "%)"), y = paste0("PC2 (", pv_Central[2], "%)"))
pMainE <- ggplot(pcaData_South, aes(PC1, PC2, colour = exposure, shape = shore_level)) + geom_point(size = 2) +
  labs(title = "South", x = paste0("PC1 (", pv_South[1], "%)"), y = paste0("PC2 (", pv_South[2], "%)"))

# Combine the three panels (labelled C-E) with a shared legend; everything after "&" applies to all panels.
p_figure1CDE <- pMainC + pMainD + pMainE + plot_layout(ncol = 3, guides = "collect") + plot_annotation(tag_levels = list(c("C", "D", "E"))) &
  scale_colour_manual(values = c("No/short exposure" = "#00A1D5", "Long exposure" = "#DC5E83"), breaks = c("No/short exposure", "Long exposure"), name = "Exposure group") &
  scale_shape_manual(values = c("high" = 17, "low" = 19), labels = c("High shore", "Low shore"), name = "Shore level") &
  guides(shape = guide_legend(order = 1), colour = guide_legend(order = 2)) &
  theme(plot.title = element_text(hjust = 0.5, vjust = 2, face = "bold", size = 16), panel.grid = element_blank(),
        panel.border = element_rect(colour = "#969696", fill = NA, linewidth = 0.8), panel.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(size = 14, face = "bold"), legend.text = element_text(size = 11), legend.title = element_text(size = 14))
p_figure1CDE

# Save the figure.
ggsave("Figure1CDE.pdf", p_figure1CDE, width = 13, height = 4)

# PERMANOVA within each population: Euclidean distances of VST expression values,
# testing exposure group and shore level (99999 permutations).
dist_North <- dist(t(assay(vsd_North)), method = "euclidean")
dist_Central <- dist(t(assay(vsd_Central)), method = "euclidean")
dist_South <- dist(t(assay(vsd_South)), method = "euclidean")

set.seed(123)
adonis_North <- adonis2(dist_North ~ exposure + shore_level, data = pcaData_North, permutations = 99999, by = "margin")
set.seed(123)
adonis_Central <- adonis2(dist_Central ~ exposure + shore_level, data = pcaData_Central, permutations = 99999, by = "margin")
set.seed(123)
adonis_South <- adonis2(dist_South ~ exposure + shore_level, data = pcaData_South, permutations = 99999, by = "margin")
adonis_North
adonis_Central
adonis_South

