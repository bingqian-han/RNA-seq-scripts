###################################country:three countries
#library(RColorBrewer)
#display.brewer.all()
#brewer.pal(9, "Set1")[c(2:6, 9, 7, 8, 1)]

library(DESeq2)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(rgl)

setwd("C:\\Users\\hanbq\\Desktop\\RNA_R_202409/PCA/")
counts <- read.table("counts_filter.txt", header = TRUE)
group <- read.table("group_order.txt", header = TRUE,sep="\t")
rownames(group) <- group$id
group <- group[colnames(counts), ]
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = group,
                              design = ~group)
vsd <- vst(dds, blind = TRUE)
head(assay(vsd), 3)

# PCA for PC1 vs PC2
pcaData_1_2 <- plotPCA(vsd, intgroup = "group", ntop = 20716, pcsToUse = c(1, 2), returnData = TRUE)
percentVar_1_2 <- round(100 * attr(pcaData_1_2, "percentVar"),2)

# PCA for PC2 vs PC3
pcaData_2_3 <- plotPCA(vsd, intgroup = "group", ntop = 20716, pcsToUse = c(2, 3), returnData = TRUE)
percentVar_2_3 <- round(100 * attr(pcaData_2_3, "percentVar"),2)

# Combining variance
percentVar <- c(percentVar_1_2[1], percentVar_1_2[2], percentVar_2_3[2])

# Combine PCA results
pcaData <- pcaData_1_2[, c("PC1", "PC2", "group")]
pcaData$PC3 <- pcaData_2_3$PC3

# Add country information
pcaData$country <- group$country
country_colors <- c("Ersvika" = "#87aade",   # 深绿色
                    "Ursholmen" = "#ffe57f",   # 柔和的粉紫色
                    "Canabal" = "#ff9999")    # 深橙色


library(scatterplot3d)
pcaData$color <- country_colors[pcaData$country]

pdf("PCA_3D_plot111.pdf", width = 5, height = 5)

par(cex.lab = 0.92,   # 11pt
    cex.axis = 0.95)  # 9pt

scatterplot3d(
  pcaData$PC1, pcaData$PC2, pcaData$PC3,
  color = pcaData$color,
  pch = 16,
  xlab = paste0("PC1 (", percentVar[1], "% variance)"),
  ylab = paste0("PC2 (", percentVar[2], "% variance)"),
  zlab = paste0("PC3 (", percentVar[3], "% variance)"),
  main = "",
  angle = 135
)

dev.off()


pdf("PCA_3D_plot222.pdf", width = 5, height = 5)

par(cex.lab = 0.92, cex.axis = 0.95)

scatterplot3d(
  pcaData$PC1, pcaData$PC2, pcaData$PC3,
  color = pcaData$color,
  pch = 16,
  xlab = "",      # 不显示 PC1
  ylab = "",      # 不显示 PC2
  zlab = "",      # 不显示 PC3
  main = "",
  angle = 135
)

dev.off()



###################################country:norway
library(DESeq2)
library(dplyr)
library(ggplot2)

setwd("C:\\Users\\hanbq\\Desktop\\RNA_R_202409/PCA/")
counts <- read.table("counts_filter.txt",header = T)
group <- read.table("group_order.txt",h=T,sep = "\t")
rownames(group) <- group$id
group <- group[colnames(counts),]
counts <- counts[, 1:34]
group <-group[1:34, ]
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = group,
                              design = ~group)
vsd <- vst(dds, blind=TRUE)
head(assay(vsd), 3)

pcaData <- plotPCA(vsd, intgroup = "group", ntop = 20716, pcsToUse = c(1, 2), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"),2)  # 提取PC1和PC2的方差百分比

# Add descriptive columns for Shore Level, Temperature, and Time
pcaData <- pcaData %>%
  mutate(
    shore_level = factor(case_when(
      grepl("^NH", group) ~ "High shore",
      grepl("^NL", group) ~ "Low shore"
    )),
    exposure_group = factor(case_when(
      grepl("31°C", group) & grepl("3hour$", group) ~ "Long exposure",
      TRUE ~ "No/short exposure"
    ),
    levels = c("No/short exposure", "Long exposure"))
  )

# Define colors, shapes, and sizes for each category
exp_colors <- c(
  "No/short exposure" = "#00A1D5",  # 淡蓝
  "Long exposure" = "#DC5E83"       # 淡橙
)

shore_level_shapes <- c("High shore" = 17, "Low shore" = 19)

# Plot PCA with custom legends

#pdf(file="PCA_Norway_PC1_PC2.pdf",width=6,height=5)

p1 <- ggplot(pcaData, aes(x = PC1, y = PC2, 
                          color = exposure_group,
                          shape = shore_level)) +
  geom_point(size = 2) +
  ggtitle("North") +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2, face = "bold", size = 16),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "#969696", fill = NA, linewidth = 0.8),
    panel.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 14) 
  ) +
  scale_color_manual(name = "Exposure group", values = exp_colors) +
  scale_shape_manual(name = "Shore level", values = shore_level_shapes) +
  guides(
    shape = guide_legend(order = 1),   # ← shore level 在上
    color = guide_legend(order = 2)    # ← exposure group 在下
  ) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)"))

p1


#dev.off()



# PERMANOVA 分析
library(vegan)
expr_matrix <- t(assay(vsd))  # 样本在行
dist_matrix <- dist(expr_matrix, method = "euclidean")
group_info <- group[1:34, ]
group_info$shore_level <- ifelse(grepl("^NH", group_info$group), 
                                 "High shore", "Low shore")
group_info$temperature <- ifelse(grepl("17", group_info$group), "17°C", "31°C")
group_info$time <- ifelse(grepl("15min$", group_info$group), "15 mins", "3 hours")
group_info$exposure_group <- ifelse(
  group_info$temperature == "31°C" & group_info$time == "3 hours",
  "Long exposure",
  "No/short exposure"
)
group_info$shore_level <- factor(group_info$shore_level)
group_info$exposure_group <- factor(group_info$exposure_group)

adonis_result <- adonis2(dist_matrix ~ exposure_group + shore_level,
                         data = group_info,
                         permutations = 999,
                         method = "euclidean")

print(adonis_result)


###################################country:sweden
library(DESeq2)
library(dplyr)
library(ggplot2)

setwd("C:\\Users\\hanbq\\Desktop\\RNA_R_202409/PCA/")
counts <- read.table("counts_filter.txt",header = T)
group <- read.table("group_order.txt",h=T,sep = "\t")
rownames(group) <- group$id
group <- group[colnames(counts),]
counts <- counts[, 35:68]
group <-group[35:68, ]

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = group,
                              design = ~group)
vsd <- vst(dds, blind=TRUE)

pcaData <- plotPCA(vsd, intgroup = "group", ntop = 20716, pcsToUse = c(1, 2), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"),2)

# --- 修改为和 Norway 一样的分组方式 ---
pcaData <- pcaData %>%
  mutate(
    shore_level = factor(case_when(
      grepl("^WH", group) ~ "High shore",
      grepl("^WL", group) ~ "Low shore"
    )),
    exposure_group = factor(case_when(
      grepl("31°C", group) & grepl("3hour$", group) ~ "Long exposure",
      TRUE ~ "No/short exposure"
    ),
    levels = c("No/short exposure", "Long exposure"))
  )

# 颜色（与挪威一致）
exp_colors <- c(
  "No/short exposure" = "#00A1D5",
  "Long exposure" = "#DC5E83"
)

# 形状（与挪威一致）
shore_level_shapes <- c("High shore" = 17, "Low shore" = 19)

# --- PCA 绘图 ---
p2 <- ggplot(pcaData, aes(x = PC1, y = PC2,
                          color = exposure_group,
                          shape = shore_level)) +
  geom_point(size = 2) +
  ggtitle("Central") +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2, face = "bold", size = 16),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "#969696", fill = NA, linewidth = 0.8),
    panel.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 14) 
  ) +
  scale_color_manual(name = "Exposure group", values = exp_colors) +
  scale_shape_manual(name = "Shore level", values = shore_level_shapes) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)"))

p2


######################PERMANOVA 分析 Sweden
library(vegan)
expr_matrix <- t(assay(vsd))  # 样本在行
dist_matrix <- dist(expr_matrix, method = "euclidean")
group_info <- group[1:34, ]
group_info$shore_level <- ifelse(grepl("^WH", group_info$group), 
                                 "High shore", "Low shore")
group_info$temperature <- ifelse(grepl("17", group_info$group), "17°C", "31°C")
group_info$time <- ifelse(grepl("15min$", group_info$group), "15 mins", "3 hours")
group_info$exposure_group <- ifelse(
  group_info$temperature == "31°C" & group_info$time == "3 hours",
  "Long exposure",
  "No/short exposure"
)
group_info$shore_level <- factor(group_info$shore_level)
group_info$exposure_group <- factor(group_info$exposure_group)

adonis_result <- adonis2(dist_matrix ~ exposure_group + shore_level,
                         data = group_info,
                         permutations = 999,
                         method = "euclidean")

print(adonis_result)







###################################country:spain
setwd("C:\\Users\\hanbq\\Desktop\\RNA_R_202409/PCA/")
counts <- read.table("counts_filter.txt",header = T)
group <- read.table("group_order.txt",h=T,sep = "\t")

rownames(group) <- group$id
group <- group[colnames(counts),]
counts <- counts[, 69:101]
group <- group[69:101, ]

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = group,
                              design = ~group)
vsd <- vst(dds, blind=TRUE)

pcaData <- plotPCA(vsd, intgroup = "group", ntop = 20716, 
                   pcsToUse = c(1, 2), returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"),2)




pcaData <- pcaData %>%
  mutate(
    shore_level = factor(case_when(
      grepl("^PH", group) ~ "High shore",
      grepl("^PL", group) ~ "Low shore"
    )),
    exposure_group = factor(case_when(
      grepl("31°C", group) & grepl("3hour$", group) ~ "Long exposure",
      TRUE ~ "No/short exposure"
    ),
    levels = c("No/short exposure", "Long exposure"))
  )

exp_colors <- c(
  "No/short exposure" = "#00A1D5",
  "Long exposure" = "#DC5E83"
)

shore_level_shapes <- c("High shore" = 17, "Low shore" = 19)

# --- PC1 vs PC2 ---
p3 <- ggplot(pcaData, aes(x = PC1, y = PC2,
                               color = exposure_group,
                               shape = shore_level)) +
  geom_point(size = 2) +
  ggtitle("South") +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2, face = "bold", size = 16),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "#969696", fill = NA, linewidth = 0.8),
    panel.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 14) 
  ) +
  scale_color_manual(name = "Exposure group", values = exp_colors) +
  scale_shape_manual(name = "Shore level", values = shore_level_shapes) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)"))

p3

###########################PERMANOVA Spain
library(vegan)
expr_matrix <- t(assay(vsd))  # 样本在行
dist_matrix <- dist(expr_matrix, method = "euclidean")
group_info <- group[1:33, ]
group_info$shore_level <- ifelse(grepl("^PH", group_info$group), 
                                 "High shore", "Low shore")
group_info$temperature <- ifelse(grepl("17", group_info$group), "17°C", "31°C")
group_info$time <- ifelse(grepl("15min$", group_info$group), "15 mins", "3 hours")
group_info$exposure_group <- ifelse(
  group_info$temperature == "31°C" & group_info$time == "3 hours",
  "Long exposure",
  "No/short exposure"
)
group_info$shore_level <- factor(group_info$shore_level)
group_info$exposure_group <- factor(group_info$exposure_group)

adonis_result <- adonis2(dist_matrix ~ exposure_group + shore_level,
                         data = group_info,
                         permutations = 999,
                         method = "euclidean")

print(adonis_result)



library(ggplot2)
library(ggpubr)

combined_plot <- ggarrange(p1, p2, p3,
                           labels = c("C", "D", "E"),
                           ncol = 3, nrow = 1,
                           common.legend = TRUE,
                           legend = "right",
                           font.label = list(size = 14, face = "bold"))

# 保存为 PDF
ggsave("combined_plot0407.pdf", combined_plot, width = 13, height = 4)

combined_plot <- ggarrange(p1, p1, 
                           p1, p1, 
                           p1, p1,
                           labels = c("A", "B", "C", "D", "E", "F"),
                           ncol = 3, nrow = 2,
                           common.legend = TRUE,
                           legend = "right",
                           font.label = list(size = 14, face = "bold"))

# 保存为 PDF
ggsave("blank22.pdf", combined_plot, width = 14, height = 9)

