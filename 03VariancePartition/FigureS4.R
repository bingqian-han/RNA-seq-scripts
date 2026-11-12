setwd("C:/path/to/RNAgithub/03VariancePartition")

# =============================================================================
# Figure S4: variance partitioning of gene expression across all samples (Supplementary Text 4).
# One linear mixed model per gene, with population, exposure group, shore level, population x exposure group and
# exposure group x shore level as random effects.
# Fits the models and saves variancePartition.rds (run on the server)
# Input: counts_filter.txt, group.txt (from 01counts_ind_filter)
# Output: variancePartition.rds, FigureS4.pdf
# =============================================================================

#Fit the models
library(variancePartition)
library(edgeR)
library(BiocParallel)
library(ggplot2)

# Counts and sample information.
counts_filter <- as.matrix(read.delim("../01counts_ind_filter/counts_filter.txt", header = TRUE, row.names = 1, check.names = FALSE))
storage.mode(counts_filter) <- "numeric"
group <- read.delim("../01counts_ind_filter/group.txt", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(group) <- group$id
group <- group[colnames(counts_filter), ]

# Population (named by country in the model: North = Norway, Central = Sweden, South = Spain), shore level and exposure group.
group$country <- factor(c(North = "Norway", Central = "Sweden", South = "Spain")[group$population])
group$shore_level <- factor(group$shore_level)
group$exposure_group <- factor(ifelse(group$temperature == "31degree" & group$time == "3h", "Long", "No/short"), levels = c("No/short", "Long"))

# Model terms.
form <- ~ (1 | country) + (1 | shore_level) + (1 | exposure_group) + (1 | exposure_group:shore_level) + (1 | country:exposure_group)
# TMM normalisation and voom weights, then the proportion of variance explained by each term for every gene.
param <- MulticoreParam(workers = 4, progressbar = TRUE)
dge <- DGEList(counts = counts_filter)
dge <- calcNormFactors(dge, method = "TMM")
vobj <- voomWithDreamWeights(dge, form, group, plot = FALSE, save.plot = FALSE, BPPARAM = param)
varPart <- fitExtractVarPartModel(vobj, form, group, BPPARAM = param)
saveRDS(varPart, "variancePartition.rds")


#Figure S4
library(variancePartition)
library(ggplot2)
varPart <- readRDS("variancePartition.rds")

# One violin per term across all genes; residuals in light grey.
p <- plotVarPart(varPart, col = c(rep("#8C8C8C", 5), "#D9D9D9")) +
  scale_x_discrete(limits = c("country", "exposure_group", "shore_level", "country:exposure_group", "exposure_group:shore_level", "Residuals"),
                   labels = c("country" = "Population", "exposure_group" = "Exposure group", "shore_level" = "Shore level", "country:exposure_group" = "Population \u00D7\nexposure group", "exposure_group:shore_level" = "Exposure group \u00D7\nshore level", "Residuals" = "Residuals")) +
  labs(x = NULL, y = "Variance explained (%)") +
  theme_bw(base_size = 10, base_family = "Arial") +
  theme(axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5, lineheight = 0.95, colour = "black", margin = margin(t = 4)),
        axis.text.y = element_text(size = 9, colour = "black"), axis.title.y = element_text(size = 10, margin = margin(r = 5)),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#EBEBEB", linewidth = 0.25),
        panel.border = element_rect(colour = "#666666", fill = NA, linewidth = 0.4), axis.ticks = element_line(linewidth = 0.3),
        legend.position = "none", plot.margin = margin(4, 4, 4, 4))

print(p)
ggsave("FigureS4.pdf", plot = p, width = 17, height = 8, units = "cm", device = cairo_pdf, bg = "white")