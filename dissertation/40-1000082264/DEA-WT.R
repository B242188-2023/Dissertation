library(DESeq2)
setwd("/Users/chenziyin/Downloads/Dissertation")

# 样本名称
WT_sampleNames <- c("WT-LD1", "WT-LD2", "WT-LL1", "WT-LL2", "WT-SD1", "WT-SD2", "WT-SL1", "WT-SL2")
for (name in WT_sampleNames) {
  file_path <- paste0(name, "_counts.txt")
  data <- read.table(file_path, header = TRUE, row.names = 1)
  assign(paste0(name, "_counts"), data)
}
WT_DataMatrix <- do.call(cbind, lapply(WT_sampleNames, function(name) {
  get(paste0(name, "_counts"))[6]
}))
'WT_DataMatrix <- cbind(WT-LD1_counts[6],WT-LD2_counts[6],WT-LL1_counts[6],
                       WT-LL2_counts[6], WT-SD1_counts[6],WT-SD2_counts[6],
                       WT-SL1_counts[6],WT-SL2_counts[6])'

sampleNames <- c("LD1", "LD2", "LL1", "LL2", "SD1", "SD2", "SL1", "SL2", "XD1", "XD2")
for (name in sampleNames) {
  file_path <- paste0(name, "_counts.txt")
  data <- read.table(file_path, header = TRUE, row.names = 1)
  assign(paste0(name, "_counts"), data)
}
DataMatrix <- do.call(cbind, lapply(sampleNames, function(name) {
  get(paste0(name, "_counts"))[6]
}))
#LD1 & LD2 <-> LL1 & LL2 
#LD1 & LD2 <-> SD1 & SD2
#LL1 & LL2 <-> SL1 & SL2
#SD1 & SD2 <-> SL1 & SL2
colnames(WT_DataMatrix) <- WT_sampleNames
colnames(DataMatrix) <- sampleNames

# 定义组信息
comparisons <- list(
  LD_LL = list(samples = c("LD1", "LD2", "LL1", "LL2"), group = c("LD", "LD", "LL", "LL")),
  LD_SD = list(samples = c("LD1", "LD2", "SD1", "SD2"), group = c("LD", "LD", "SD", "SD")),
  LL_SL = list(samples = c("LL1", "LL2", "SL1", "SL2"), group = c("LL", "LL", "SL", "SL")),
  SD_SL = list(samples = c("SD1", "SD2", "SL1", "SL2"), group = c("SD", "SD", "SL", "SL")),
  XD_SD  = list(samples = c("XD1", "XD2", "SD1", "SD2"), group = c("XD", "XD", "SD", "SD")),
  XD_LD  = list(samples = c("XD1", "XD2", "LD1", "LD2"), group = c("XD", "XD", "LD", "LD"))
  )

WT_comparisons <- list(
  WT_LD_LL = list(samples = c("WT_LD1", "WT_LD2", "WT_LL1", "WT_LL2"), group = c("WT_LD", "WT_LD", "WT_LL", "WT_LL")),
  WT_LD_SD = list(samples = c("WT_LD1", "WT_LD2", "WT_SD1", "WT_SD2"), group = c("WT_LD", "WT_LD", "WT_SD", "WT_SD")),
  WT_LL_SL = list(samples = c("WT_LL1", "WT_LL2", "WT_SL1", "WT_SL2"), group = c("WT_LL", "WT_LL", "WT_SL", "WT_SL")),
  WT_SD_SL = list(samples = c("WT_SD1", "WT_SD2", "WT_SL1", "WT_SL2"), group = c("WT_SD", "WT_SD", "WT_SL", "WT_SL"))
)



library(rtracklayer)

# 读取GFF3文件
gff3_file <- "AeUmbellulata_TA1851_v1.gff3"
gff3_data <- import(gff3_file, format = "gff3")

library(stringr)
# 提取基因ID和注释信息
gene_annotations <- gff3_data[gff3_data$type == "gene", ]
gene_annotations <- as.data.frame(gene_annotations)
# 计算基因长度
gene_annotations$gene_length <- gene_annotations$end - gene_annotations$start + 1

# 定义提取基因名称的函数
extract_gene_name <- function(note) {
  match <- str_match(note, "Similar to ([^:]+):")
  return(match[,2])
}
gene_annotations$gene_name <- sapply(gene_annotations$Note, extract_gene_name)
gene_annotations <- gene_annotations[, c("ID", "seqnames", "start", "end", "gene_length", "gene_name", "Ontology_term")]

gene_length_vector <- gene_annotations$gene_length
names(gene_length_vector) <- gene_annotations$ID


all_samples <- list()
all_groups <- list()

library(dplyr)
for (comparison in names(comparisons)) {
  samples <- comparisons[[comparison]]$samples
  group <- comparisons[[comparison]]$group
  
  DataMatrix_sub <- DataMatrix[, samples]
  colData <- data.frame(row.names = samples, condition = factor(group))
  
  dds <- DESeqDataSetFromMatrix(countData = DataMatrix_sub,
                              colData = colData,
                              design= ~condition)
  dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = TRUE)
  #注意，需将 Long 在前，short 在后，意为 long 相较于 short 中哪些基因上调/下调
  res <- results(dds1, contrast = c('condition', unique(group)))
  print(unique(group))
  #res
  res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  #write.table(res1, 'DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)

  ##筛选差异表达基因
  #首先对表格排个序，按 padj 值升序排序，相同 padj 值下继续按 log2FC 降序排序
  res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

  #log2FC≥1 & padj<0.01 标识 up，代表显著上调的基因
  #log2FC≤-1 & padj<0.01 标识 down，代表显著下调的基因
  #其余标识 none，代表非差异的基因
  #5e-324范围内，超出这个范围就显示为0.000000e+00
  res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.05 & res1$baseMean >= 100), 'sig'] <- 'up'
  res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.05 & res1$baseMean >= 100), 'sig'] <- 'down'
  res1[which(abs(res1$log2FoldChange) < 1 | res1$padj >= 0.05 | res1$baseMean < 100), 'sig'] <- 'none'
  #LL& SL
  #LD&SD top25
  
  res1$ensembl_gene_id <- rownames(res1)
  merged_results <- merge(gene_annotations, res1, by.x = "ID",
                          by.y = "ensembl_gene_id", all.x = F, all.y = T)
  res1$ensembl_gene_id <- NULL
  write.table(merged_results, file = paste0(comparison, '_gene_info.txt'), sep = '\t',
              col.names = NA, quote = F)
  
  #根据 up 和 down 分开输出
  res1_up <- subset(merged_results, sig == 'up')
  res1_down <- subset(merged_results, sig == 'down')
  
  write.table(merged_results, file = paste0(comparison, '.DESeq2.txt'),
              sep = '\t', col.names = NA, quote = FALSE)
  write.table(res1_up, file = paste0(comparison, '.DESeq2.up.txt'), sep = '\t',
              col.names = NA, quote = FALSE)
  write.table(res1_down, file = paste0(comparison, '.DESeq2.down.txt'),
              sep = '\t', col.names = NA, quote = FALSE)
  ##ggplot2 差异火山图
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)

  # 筛选前25个差异显著的基因
  top25 <- merged_results[order(-abs(merged_results$log2FoldChange), merged_results$padj), ][1:25, ]
  top25_up <- res1_up[order(-abs(res1_up$log2FoldChange), res1_up$padj), ][1:25, ]
  top25_down <- res1_down[order(-abs(res1_down$log2FoldChange), res1_down$padj), ][1:25, ]
  highlight_genes <- rbind(top25, top25_up, top25_down)
  write.table(highlight_genes,file = paste0(comparison,'.highlight.txt'),sep = '\t', col.name = NA, quote = F)

  p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.4, size = 2) +
    scale_color_manual(values = c("#ff4757", "#d2dae2", "#546de5"), limits = c('up', 'none', 'down')) +
    labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = paste(comparison, 'Comparison'), color = '') +
    theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          legend.key = element_rect(fill = 'transparent')) +
    geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +
    geom_hline(yintercept = 2, lty = 3, color = 'black') +
    xlim(-15, 15) + ylim(0, 250)
  p <- p + geom_text(data = highlight_genes, aes(label = gene_name),
                     vjust = 2, hjust = 0.5, check_overlap = TRUE, size = 3)
  p
  ggsave(filename = paste0(comparison, '_volcano_plot.png'), plot = p)
  }
  

  
  
for (WT_comparison in names(WT_comparisons)) {
  samples <- WT_comparisons[[WT_comparison]]$samples
  samples <- str_replace_all(samples, "_", "-")
  group <- WT_comparisons[[WT_comparison]]$group
  
  DataMatrix_sub <- WT_DataMatrix[, samples]
  colData <- data.frame(row.names = samples, condition = factor(group))
  
  dds <- DESeqDataSetFromMatrix(countData = DataMatrix_sub,
                                colData = colData,
                                design= ~condition)
  dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = TRUE)
  #注意，需将 Long 在前，short 在后，意为 long 相较于 short 中哪些基因上调/下调
  res <- results(dds1, contrast = c('condition', unique(group)))
  print(unique(group))
  #res
  res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  #write.table(res1, 'DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)
  
  ##筛选差异表达基因
  #首先对表格排个序，按 padj 值升序排序，相同 padj 值下继续按 log2FC 降序排序
  res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
  
  #log2FC≥1 & padj<0.01 标识 up，代表显著上调的基因
  #log2FC≤-1 & padj<0.01 标识 down，代表显著下调的基因
  #其余标识 none，代表非差异的基因
  #5e-324范围内，超出这个范围就显示为0.000000e+00
  res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.05 & res1$baseMean >= 100), 'sig'] <- 'up'
  res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.05 & res1$baseMean >= 100), 'sig'] <- 'down'
  res1[which(abs(res1$log2FoldChange) < 1 | res1$padj >= 0.05 | res1$baseMean < 100), 'sig'] <- 'none'
  #LL& SL
  #LD&SD top25
  
  res1$ensembl_gene_id <- rownames(res1)
  merged_results <- merge(gene_annotations, res1, by.x = "ID",
                          by.y = "ensembl_gene_id", all.x = F, all.y = T)
  res1$ensembl_gene_id <- NULL
  write.table(merged_results, file = paste0(WT_comparison, '_gene_info.txt'), sep = '\t',
              col.names = NA, quote = F)
  
  #根据 up 和 down 分开输出
  res1_up <- subset(merged_results, sig == 'up')
  res1_down <- subset(merged_results, sig == 'down')
  
  write.table(merged_results, file = paste0(WT_comparison, '.DESeq2.txt'),
              sep = '\t', col.names = NA, quote = FALSE)
  write.table(res1_up, file = paste0(WT_comparison, '.DESeq2.up.txt'), sep = '\t',
              col.names = NA, quote = FALSE)
  write.table(res1_down, file = paste0(WT_comparison, '.DESeq2.down.txt'),
              sep = '\t', col.names = NA, quote = FALSE)
  ##ggplot2 差异火山图
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
  
  # 筛选前25个差异显著的基因
  top25 <- merged_results[order(-abs(merged_results$log2FoldChange), merged_results$padj), ][1:25, ]
  top25_up <- res1_up[order(-abs(res1_up$log2FoldChange), res1_up$padj), ][1:25, ]
  top25_down <- res1_down[order(-abs(res1_down$log2FoldChange), res1_down$padj), ][1:25, ]
  highlight_genes <- rbind(top25, top25_up, top25_down)
  write.table(highlight_genes,file = paste0(WT_comparison,'.highlight.txt'),sep = '\t', col.name = NA, quote = F)
  
  #highlight_genes <- rbind(res1_up,res1_down)
  #target <- c('AeUmb.TA1851.r1.4UG0030010',
  #'AeUmb.TA1851.r1.5UG0039190',
  #'AeUmb.TA1851.r1.2UG0045860',
  #'AeUmb.TA1851.r1.6UG0050370')
  #target_genes <- merged_results %>% filter(merged_results$ensembl_gene_id %in% target)
  
  p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.4, size = 2) +
    scale_color_manual(values = c("#ff4757", "#d2dae2", "#546de5"), limits = c('up', 'none', 'down')) +
    labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = paste(WT_comparison, 'Comparison'), color = '') +
    theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          legend.key = element_rect(fill = 'transparent')) +
    geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +
    geom_hline(yintercept = 2, lty = 3, color = 'black') +
    xlim(-15, 15) + ylim(0, 250)
  p <- p + geom_text(data = highlight_genes, aes(label = gene_name),
                     vjust = 2, hjust = 0.5, check_overlap = TRUE, size = 3)
  ggsave(filename = paste0(WT_comparison, '_volcano_plot.png'), plot = p)
}



# 合并所有组的数据和组信息
combined_data <- do.call(cbind, all_samples)
combined_colData <- do.call(rbind, all_groups)

# 重新创建DESeqDataSet对象并进行PCA分析
combined_dds <- DESeqDataSetFromMatrix(countData = combined_data,
                                       colData = combined_colData,
                                       design = ~ condition)
vsd_combined <- vst(combined_dds, blind=FALSE)
pcaData_combined <- plotPCA(vsd_combined, intgroup=c("condition"), returnData=TRUE)
percentVar_combined <- round(100 * attr(pcaData_combined, "percentVar"))



# 提取不同光周期组的PC1和PC2得分以及组信息
pc1_scores <- pcaData_combined$PC1
pc2_scores <- pcaData_combined$PC2
conditions <- pcaData_combined$condition

# 将因子变量转换为数值
pcaData_combined$condition_numeric <- as.numeric(as.factor(pcaData_combined$condition))
condition_numeric <- pcaData_combined$condition_numeric

# 进行PC1的回归分析和相关性分析
lm_model_PC1 <- lm(pc1_scores ~ condition_numeric)
correlation_PC1 <- cor(pc1_scores, condition_numeric)

# 进行PC2的回归分析和相关性分析
lm_model_PC2 <- lm(pc2_scores ~ condition_numeric)
correlation_PC2 <- cor(pc2_scores, condition_numeric)

# 打印回归模型摘要
summary(lm_model_PC1)
summary(lm_model_PC2)

# 打印相关性
print(correlation_PC1)
print(correlation_PC2)

library(gridExtra)
# 可视化PC1的回归和相关性结果
pca_plot_PC1 <- ggplot(pcaData_combined, aes(x=condition_numeric, y=PC1, color=conditions)) +
  geom_point(size=3) +
  geom_smooth(method="lm", se=FALSE, color="black") +
  xlab("Conditions (Numeric)") +
  ylab(paste0("PC1: ", percentVar_combined[1], "% variance")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size=16), # 调整 x 轴标签的字体大小
    axis.title.y = element_text(size=16), # 调整 y 轴标签的字体大小
    axis.text = element_text(size=20)
  ) +
  ggtitle("PC1 of RNA-seq data with Regression") +
  annotate("text", x=mean(condition_numeric), y=max(pcaData_combined$PC1), 
           label=paste0("r = ", round(cor(pcaData_combined$PC1, condition_numeric), 2)), hjust=0, vjust=1) +
  annotation_custom(grid::textGrob("(B)", x=unit(0, "npc"), y=unit(1, "npc"), 
                                   just=c("left", "top"), gp=grid::gpar(col="black", fontsize=14, fontface="bold")))

# 可视化PC2的回归和相关性结果
pca_plot_PC2 <- ggplot(pcaData_combined, aes(x=condition_numeric, y=PC2, color=conditions)) +
  geom_point(size=3) +
  geom_smooth(method="lm", se=FALSE, color="black") +
  xlab("Conditions (Numeric)") +
  ylab(paste0("PC2: ", percentVar_combined[2], "% variance")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size=16), # 调整 x 轴标签的字体大小
    axis.title.y = element_text(size=16), # 调整 y 轴标签的字体大小
    axis.text = element_text(size=20)
  ) +
  ggtitle("PC2 of RNA-seq data with Regression") +
  annotate("text", x=mean(condition_numeric), y=max(pcaData_combined$PC2), 
           label=paste0("r = ", round(cor(pcaData_combined$PC2, condition_numeric), 2)), hjust=0, vjust=1) +
  annotation_custom(grid::textGrob("(C)", x=unit(0, "npc"), y=unit(1, "npc"), 
                                   just=c("left", "top"), gp=grid::gpar(col="black", fontsize=14, fontface="bold")))

# 绘制合并数据的PCA图
pca_plot_combined <- ggplot(pcaData_combined, aes(x=PC1, y=PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar_combined[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_combined[2], "% variance")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size=16), # 调整 x 轴标签的字体大小
    axis.title.y = element_text(size=16), # 调整 y 轴标签的字体大小
    axis.text = element_text(size=20)
  ) +
  ggtitle("PCA of RNA-seq data") +
  annotation_custom(grid::textGrob("(A)", x=unit(0, "npc"), y=unit(1, "npc"), 
                                   just=c("left", "top"), gp=grid::gpar(col="black", fontsize=14, fontface="bold")))

combined_plot <- grid.arrange(
  pca_plot_combined, 
  arrangeGrob(pca_plot_PC1, pca_plot_PC2, ncol = 1),
  ncol = 2
)
# 保存合并数据的PCA图
png("combined_PCA_plots.png", width = 15, height = 8, units = "in", res = 300)
grid.arrange(pca_plot_combined, arrangeGrob(pca_plot_PC1, pca_plot_PC2, ncol = 1), ncol = 2)
dev.off()




library(VennDiagram)
# 读取保存的文件
LD_LL_up <- read.table("LD_LL.DESeq2.up.txt", sep = '\t', header = TRUE, quote = "")
LD_LL_down <- read.table("LD_LL.DESeq2.down.txt", sep = '\t', header = TRUE, quote = "")
LD_SD_up <- read.table("LD_SD.DESeq2.up.txt", sep = '\t', header = TRUE, quote = "")
LD_SD_down <- read.table("LD_SD.DESeq2.down.txt", sep = '\t', header = TRUE, quote = "")
LL_SL_up <- read.table("LL_SL.DESeq2.up.txt", sep = '\t', header = TRUE, quote = "")
LL_SL_down <- read.table("LL_SL.DESeq2.down.txt", sep = '\t', header = TRUE, quote = "")
SD_SL_up <- read.table("SD_SL.DESeq2.up.txt", sep = '\t', header = TRUE, quote = "")
SD_SL_down <- read.table("SD_SL.DESeq2.down.txt", sep = '\t', header = TRUE, quote = "")

LD_LL_up <- LD_LL_up[!is.na(LD_LL_up$gene_name), ]
LD_LL_down <- LD_LL_down[!is.na(LD_LL_down$gene_name), ]
LD_SD_up <- LD_SD_up[!is.na(LD_SD_up$gene_name), ]
LD_SD_down <- LD_SD_down[!is.na(LD_SD_down$gene_name), ]
LL_SL_up <- LL_SL_up[!is.na(LL_SL_up$gene_name), ]
LL_SL_down <- LL_SL_down[!is.na(LL_SL_down$gene_name), ]
SD_SL_up <- SD_SL_up[!is.na(SD_SL_up$gene_name), ]
SD_SL_down <- SD_SL_down[!is.na(SD_SL_down$gene_name), ]


WT_LD_LL_up <- read.table("WT_LD_LL.DESeq2.up.txt", sep = '\t', header = TRUE, quote = "")
WT_LD_LL_down <- read.table("WT_LD_LL.DESeq2.down.txt", sep = '\t', header = TRUE, quote = "")
WT_LD_SD_up <- read.table("WT_LD_SD.DESeq2.up.txt", sep = '\t', header = TRUE, quote = "")
WT_LD_SD_down <- read.table("WT_LD_SD.DESeq2.down.txt", sep = '\t', header = TRUE, quote = "")
WT_LL_SL_up <- read.table("WT_LL_SL.DESeq2.up.txt", sep = '\t', header = TRUE, quote = "")
WT_LL_SL_down <- read.table("WT_LL_SL.DESeq2.down.txt", sep = '\t', header = TRUE, quote = "")
WT_SD_SL_up <- read.table("WT_SD_SL.DESeq2.up.txt", sep = '\t', header = TRUE, quote = "")
WT_SD_SL_down <- read.table("WT_SD_SL.DESeq2.down.txt", sep = '\t', header = TRUE, quote = "")
XD_LD_up <- read.table("XD_LD.DESeq2.up.txt", sep = '\t', header = TRUE, quote = "")
XD_LD_down <- read.table("XD_LD.DESeq2.down.txt", sep = '\t', header = TRUE, quote = "")
XD_SD_up <- read.table("XD_SD.DESeq2.up.txt", sep = '\t', header = TRUE, quote = "")
XD_SD_down <- read.table("XD_SD.DESeq2.down.txt", sep = '\t', header = TRUE, quote = "")



WT_LD_LL_up <- WT_LD_LL_up[!is.na(WT_LD_LL_up$gene_name), ]
WT_LD_LL_down <- WT_LD_LL_down[!is.na(WT_LD_LL_down$gene_name), ]
WT_LD_SD_up <- WT_LD_SD_up[!is.na(WT_LD_SD_up$gene_name), ]
WT_LD_SD_down <- WT_LD_SD_down[!is.na(WT_LD_SD_down$gene_name), ]
WT_LL_SL_up <- WT_LL_SL_up[!is.na(WT_LL_SL_up$gene_name), ]
WT_LL_SL_down <- WT_LL_SL_down[!is.na(WT_LL_SL_down$gene_name), ]
WT_SD_SL_up <- WT_SD_SL_up[!is.na(WT_SD_SL_up$gene_name), ]
WT_SD_SL_down <- WT_SD_SL_down[!is.na(WT_SD_SL_down$gene_name), ]
XD_LD_up <- XD_LD_up[!is.na(XD_LD_up$gene_name), ]
XD_LD_down <- XD_LD_down[!is.na(XD_LD_down$gene_name), ]
XD_SD_up <- XD_SD_up[!is.na(XD_SD_up$gene_name), ]
XD_SD_down <- XD_SD_down[!is.na(XD_SD_down$gene_name), ]


LD_LL_genes <- unique(c(LD_LL_up$gene_name, LD_LL_down$gene_name))
LD_SD_genes <- unique(c(LD_SD_up$gene_name, LD_SD_down$gene_name))
LL_SL_genes <- unique(c(LL_SL_up$gene_name, LL_SL_down$gene_name))
SD_SL_genes <- unique(c(SD_SL_up$gene_name, SD_SL_down$gene_name))

LD_SD_up_genes <- unique(LD_SD_up$gene_name)
LD_SD_down_genes <- unique(LD_SD_down$gene_name)
LL_SL_up_genes <- unique(LL_SL_up$gene_name)
LL_SL_down_genes <- unique(LL_SL_down$gene_name)

# 只在LD_SD和LL_SL中表达的基因
LD_SD_up_unique_genes <- setdiff(LD_SD_up_genes, c(LD_LL_genes, LL_SL_genes, SD_SL_genes, LD_SD_down_unique_genes))
LD_SD_down_unique_genes <- setdiff(LD_SD_down_genes, c(LD_LL_genes, LL_SL_genes, SD_SL_genes, LD_SD_up_unique_genes))
LL_SL_up_unique_genes <- setdiff(LL_SL_up_genes, c(LD_LL_genes, LD_SD_genes, SD_SL_genes, LL_SL_down_unique_genes))
LL_SL_down_unique_genes <- setdiff(LL_SL_down_genes, c(LD_LL_genes, LD_SD_genes, SD_SL_genes, LL_SL_up_unique_genes))

write.table(LD_SD_up_unique_genes, file = 'Venn_LD_SD_up_unique_genes.txt',
            sep = '\t', col.names = NA, quote = FALSE)
write.table(LD_SD_down_unique_genes, file = 'Venn_LD_SD_down_unique_genes.txt',
            sep = '\t', col.names = NA, quote = FALSE)
write.table(LL_SL_up_unique_genes, file = 'Venn_LL_SL_up_unique_genes.txt',
            sep = '\t', col.names = NA, quote = FALSE)
write.table(LL_SL_down_unique_genes, file = 'Venn_LL_SL_down_unique_genes.txt',
            sep = '\t', col.names = NA, quote = FALSE)

gene_list <- list(
  "LD vs LL" = LD_LL_genes,
  "LD vs SD" = LD_SD_genes,
  "LL vs SL" = LL_SL_genes,
  "SD vs SL" = SD_SL_genes
)

gene_list_1 <- list(
  "LD_SD Up" = LD_SD_up_genes,
  "LD_SD Down" = LD_SD_down_genes,
  "LL_SL Up" = LL_SL_up_genes,
  "LL_SL Down" = LL_SL_down_genes
)

gene_list_2 <- list(
  "LD_SD" = LD_SD_genes,
  "LL_SL" = LL_SL_genes
)

venn.diagram(
  x = gene_list,
  filename = "venn_diagram.png",
  fill = c("#D72E25", "#F46D43", "#1A9850", "#A6D96A"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("#D72E25", "#F46D43", "#1A9850", "#A6D96A")
)

venn.diagram(
  x = gene_list_1,
  filename = "LD_SD & LL_SL_venn_diagram.png",
  fill = c("#D72E25", "#F46D43", "#1A9850", "#A6D96A"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.0,
  cat.col = c("#D72E25", "#F46D43", "#1A9850", "#A6D96A")
)

venn.diagram(
  x = gene_list_2,
  filename = "Daytime Genes (LD_SD) vs. Nighttime Genes (LL_SL) venn diagram.png",
  fill = c("#D72E25", "#1A9850"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.0,
  cat.col = c("#D72E25", "#1A9850")
)


LD_SD_highlight_75 <- read.table("LD_SD.highlight.txt", sep = '\t', header = TRUE)
LL_SL_highlight_75 <- read.table("LL_SL.highlight.txt", sep = '\t', header = TRUE)

LD_SD_genes <- unique(LD_SD_highlight_75$ID)
LL_SL_genes <- unique(LL_SL_highlight_75$ID)

# 保存为 CSV 文件
write.csv(LD_SD_genes, "LD_SD_genes.csv", row.names = FALSE)
write.csv(LL_SL_genes, "LL_SL_genes.csv", row.names = FALSE)

write.csv(unique_groups, "unique_groups.csv", row.names = FALSE)


library(clusterProfiler)
library(dplyr)
library(tidyr)


# 提取基因ID（字符向量）
up_genes <- as.character(LL_SL_up$gene_name)  
down_genes <- as.character(LL_SL_down$gene_name)  
up_genes <- up_genes[!is.na(up_genes)]
down_genes <- down_genes[!is.na(down_genes)]

#up_genes <- as.character(SD_SL_up$gene_name)   
#down_genes <- as.character(SD_SL_down$gene_name)   
#up_genes <- up_genes[!is.na(up_genes)]
#down_genes <- down_genes[!is.na(down_genes)]

#total_gene <- c(up_genes,down_genes)
library(org.At.tair.db)
go_up <- enrichGO(gene = up_genes,
           OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
           keyType = "SYMBOL",  # 根据基因名称进行富集分析
           ont = "BP",  # 生物过程
           pvalueCutoff = 0.05,
           pAdjustMethod = "none",
           qvalueCutoff = 1)
dim(go_up)
LL_SL_go_up <- go_up@result
write.csv(LL_SL_go_up, "LL_SL_go_up.csv", row.names = FALSE)

go_down <- enrichGO(gene = down_genes,
                  OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                  keyType = "SYMBOL",  # 根据基因名称进行富集分析
                  ont = "BP",  # 生物过程
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "none",
                  qvalueCutoff = 1)
dim(go_down)
LL_SL_go_down <- go_down@result
write.csv(LL_SL_go_down, "LL_SL_go_down.csv", row.names = FALSE)


#结果可视化 
LL_SL_up_bar <- barplot(go_up,showCategory=15,drop=T) 
LL_SL_up_bar <- LL_SL_up_bar + theme(
  axis.text.y = element_text(size = 15)) + # 调整 y 轴标签的字体大小
  ggtitle("LL_SL_Upregulated_Bar")

LL_SL_up_bar
ggsave(filename = "LL_SL_GO_Enrichment_Upregulated_Bar_Font.png", plot = LL_SL_up_bar)
LL_SL_up_dot <- dotplot(go_up,showCategory=15)
LL_SL_up_dot <- LL_SL_up_dot + theme(
  axis.text.y = element_text(size = 15)) +  # 调整 y 轴标签的字体大小
  ggtitle("LL_SL_Upregulated_Dot")

LL_SL_up_dot
ggsave(filename = "LL_SL_GO_Enrichment_Upregulated_Dot_Font.png", plot = LL_SL_up_dot)


LL_SL_down_bar <- barplot(go_down,showCategory=15,drop=T) 
LL_SL_down_bar <- LL_SL_down_bar + theme(
  axis.text.y = element_text(size = 15)) +  # 调整 y 轴标签的字体大小
  ggtitle("LL_SL_Downregulated_Bar")

LL_SL_down_bar
ggsave(filename = "LL_SL_GO_Enrichment_Downregulated_Bar_Font.png", plot = LL_SL_down_bar)
LL_SL_down_dot <- dotplot(go_down,showCategory=15)
LL_SL_down_dot <- LL_SL_down_dot + theme(
  axis.text.y = element_text(size = 15)) +  # 调整 y 轴标签的字体大小
  ggtitle("LL_SL_Downregulated_Dot")

LL_SL_down_dot
ggsave(filename = "LL_SL_GO_Enrichment_Downregulated_Dot_Font.png", plot = LL_SL_down_dot)





# 提取基因ID（确保是字符向量）
up_genes <- as.character(LD_SD_up$gene_name)   
down_genes <- as.character(LD_SD_down$gene_name)   
up_genes <- up_genes[!is.na(up_genes)]
down_genes <- down_genes[!is.na(down_genes)]
library(org.At.tair.db)
go_up <- enrichGO(gene = up_genes,
                  OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                  keyType = "SYMBOL",  # 根据基因名称进行富集分析
                  ont = "BP",  # 生物过程
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "none",
                  qvalueCutoff = 1)
dim(go_up)
LD_SD_go_up <- go_up@result
write.csv(LD_SD_go_up, "LD_SD_go_up.csv", row.names = FALSE)


go_down <- enrichGO(gene = down_genes,
                    OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                    keyType = "SYMBOL",  # 根据基因名称进行富集分析
                    ont = "BP",  # 生物过程
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "none",
                    qvalueCutoff = 1)

dim(go_down)
LD_SD_go_down <- go_down@result
write.csv(LD_SD_go_down, "LD_SD_go_down.csv", row.names = FALSE)


#结果可视化 
LD_SD_up_bar <- barplot(go_up,showCategory=15,drop=T) 
LD_SD_up_bar <- LD_SD_up_bar + theme(
  axis.text.y = element_text(size = 15)) +  # 调整 y 轴标签的字体大小
  ggtitle("LD_SD_Upregulated_Bar")

LD_SD_up_bar
ggsave(filename = "LD_SD_GO_Enrichment_Upregulated_Bar_Font.png", plot = LD_SD_up_bar)
LD_SD_up_dot <- dotplot(go_up,showCategory=15)
LD_SD_up_dot <- LD_SD_up_dot + theme(
  axis.text.y = element_text(size = 15)) +  # 调整 y 轴标签的字体大小
  ggtitle("LD_SD_Upregulated_Dot")

LD_SD_up_dot
ggsave(filename = "LD_SD_GO_Enrichment_Upregulated_Dot_Font.png", plot = LD_SD_up_dot)


LD_SD_down_bar <- barplot(go_down,showCategory=15,drop=T) 
LD_SD_down_bar <- LD_SD_down_bar + theme(
  axis.text.y = element_text(size = 15)) +  # 调整 y 轴标签的字体大小
  ggtitle("LD_SD_Downregulated_Bar")

LD_SD_down_bar
ggsave(filename = "LD_SD_GO_Enrichment_Downregulated_Bar_Font.png", plot = LD_SD_down_bar)
LD_SD_down_dot <- dotplot(go_down,showCategory=15)
LD_SD_down_dot <- LD_SD_down_dot + theme(
  axis.text.y = element_text(size = 15)) +  # 调整 y 轴标签的字体大小
  ggtitle("LD_SD_Downregulated_Dot")

LD_SD_down_dot
ggsave(filename = "LD_SD_GO_Enrichment_Downregulated_Dot_Font.png", plot = LD_SD_down_dot)







# 创建一个空的向量来存储每个样本的组信息
group_info <- c()

# 遍历每个比较，并按样本的顺序填充 group_info
for (comparison in names(comparisons)) {
  # 获取当前比较的样本和组信息
  samples <- comparisons[[comparison]]$samples
  groups <- comparisons[[comparison]]$group
  
  # 确保样本和组信息对应，并将其填充到 group_info 中
  for (i in 1:length(samples)) {
    sample_name <- samples[i]
    group_info[sample_name] <- groups[i]
  }
}


group_info <- group_info[sampleNames]

# 创建所有样本和对应的组信息
colData <- data.frame(
  row.names = sampleNames,   # 样本名
  condition = factor(group_info)   # 组信息
)

# 创建 DESeqDataSet
dds_full <- DESeqDataSetFromMatrix(
  countData = DataMatrix,   # 原始表达矩阵
  colData = colData,        # 样本信息
  design = ~condition       
)

# 运行 DESeq2 分析
dds_full <- DESeq(dds_full)

# 获取归一化后的表达矩阵
normalized_counts <- counts(dds_full, normalized = TRUE)

# 保存总归一化矩阵
write.table(normalized_counts, "Normalized_Expression_Total.txt", sep = "\t", quote = FALSE, col.names = NA)




for (WT_comparison in names(WT_comparisons)) {
  samples <- WT_comparisons[[WT_comparison]]$samples
  samples <- str_replace_all(samples, "_", "-")
  group <- WT_comparisons[[WT_comparison]]$group
  
  # 确保样本和组信息对应，并将其填充到 group_info 中
  for (i in 1:length(samples)) {
    sample_name <- samples[i]
    group_info[sample_name] <- groups[i]
  }
}

# 按样本顺序排列 group_info
group_info <- group_info[WT_sampleNames]

# 创建 colData 数据框，包含所有样本和对应的组信息
colData <- data.frame(
  row.names = WT_sampleNames,   # 样本名
  condition = factor(group_info)   # 组信息，转化为因子（factor）
)

# 创建 DESeqDataSet
dds_full <- DESeqDataSetFromMatrix(
  countData = WT_DataMatrix,   # 你的原始表达矩阵
  colData = colData,        # 你的样本信息（包括组信息）
  design = ~condition       # 设计公式，基于组信息
)

# 运行 DESeq2 分析
dds_full <- DESeq(dds_full)

# 获取归一化后的表达矩阵
WT_normalized_counts <- counts(dds_full, normalized = TRUE)

# 保存总归一化矩阵
write.table(WT_normalized_counts, "WT_Normalized_Expression_Total.txt", sep = "\t", quote = FALSE, col.names = NA)

countData <- cbind(DataMatrix, WT_DataMatrix)
all_sampleNames <- colnames(countData)

k2 <- read.csv("photoperiod_response_expression_matrix.csv")
countData_2685 <- countData[rownames(countData) %in% k2$X, ]



# 初始化 group_info
group_info <- c()

# DataMatrix 部分
for (comparison in names(comparisons)) {
  samples <- comparisons[[comparison]]$samples
  groups <- comparisons[[comparison]]$group
  
  for (i in 1:length(samples)) {
    sample_name <- samples[i]
    group_info[sample_name] <- groups[i]
  }
}

# WT_DataMatrix 部分
for (WT_comparison in names(WT_comparisons)) {
  samples <- WT_comparisons[[WT_comparison]]$samples
  # 把样本名里的 _ 统一成 -
  samples <- str_replace_all(samples, "_", "-")
  groups <- WT_comparisons[[WT_comparison]]$group
  
  for (i in 1:length(samples)) {
    sample_name <- samples[i]
    group_info[sample_name] <- groups[i]
  }
}

# 按照合并后的样本顺序排列 group_info
group_info <- group_info[all_sampleNames]

# 创建 colData
colData <- data.frame(
  row.names = all_sampleNames,
  condition = factor(group_info)
)

dds <- DESeqDataSetFromMatrix(
  countData = countData_2685,
  colData = colData,
  design = ~condition
)

# 统一归一化
dds <- DESeq(dds)
sample_normalized_counts_2685 <- counts(dds, normalized=TRUE)

library(tibble)
countData_2685_GeneName <- as.data.frame(sample_normalized_counts_2685) %>%
  rownames_to_column("gene_id") %>%
  mutate(gene_name = gene_annotations$gene_name[match(gene_id, gene_annotations$ID)]) %>%
  select(gene_name, everything()) %>%   # 把 gene_name 放到第一列
  column_to_rownames("gene_id")         # 把 gene_id 恢复成行名
write.table(countData_2685_GeneName, "Normalized_Expression_2685_with_GeneName.txt", sep = "\t", quote = FALSE, col.names = NA)


write.table(sample_normalized_counts_2685, "Sample_Normalized_Expression_2685.txt", sep = "\t", quote = FALSE, col.names = NA)

test <- read.table("HC_TPM_Expression_Matrix_2685.txt", sep = '\t', header = TRUE, quote = "")




#BiocManager::install("GENIE3")
library(GENIE3)

# 读取原始 count 矩阵
#countData <- cbind(DataMatrix, WT_DataMatrix)
  
# 确保基因ID在基因长度表里存在
common_genes <- intersect(rownames(countData_2685), gene_annotations$ID)

# 保留两者共有的基因
countData_2685 <- countData_2685[common_genes, ]
gene_annotations <- gene_annotations[gene_annotations$ID %in% common_genes, ]

# 确保顺序对应
gene_annotations <- gene_annotations[match(rownames(countData_2685), gene_annotations$ID), ]


# 转换 gene length (kb)
gene_length_kb <- gene_annotations$gene_length / 1000  

# 计算 RPK
rpk <- countData_2685 / gene_length_kb

# 计算每个样本的 scaling factor（Per Million Scaling Factor）
scaling_factors <- colSums(rpk) / 1e6  

# 计算 TPM
tpm <- sweep(rpk, 2, scaling_factors, FUN = "/")

head(tpm)

write.table(tpm, "TPM_Expression_Matrix_2685.txt", sep = "\t", quote = FALSE, col.names = NA)

hc_genes <- rownames(tpm)[apply(tpm, 1, function(x) any(x > 5))]

# 查看筛出来的高表达基因数量
length(hc_genes)

# 提取这部分基因对应的 TPM 矩阵
tpm_hc <- tpm[hc_genes, ]

write.table(tpm_hc, "HC_TPM_Expression_Matrix_2685.txt", sep = "\t", quote = FALSE, col.names = NA)

tpm_hc_GeneName <- as.data.frame(tpm_hc) %>%
  rownames_to_column("gene_id") %>%
  mutate(gene_name = gene_annotations$gene_name[match(gene_id, gene_annotations$ID)]) %>%
  select(gene_name, everything()) %>%   # 把 gene_name 放到第一列
  column_to_rownames("gene_id")         # 把 gene_id 恢复成行名
write.table(tpm_hc_GeneName, "HC_TPM_Expression_Matrix_2685_with_GeneName.txt", sep = "\t", quote = FALSE, col.names = NA)

#highly_expressed_total_sample <- countData_2685[rownames(countData_2685) %in% all_dark_condition_gene, ]
#write.csv(highly_expressed_total_sample, "highly_expressed_total_samples_expression_matrix.csv", row.names = TRUE, quote = F)


#all_dark_condition <- rbind(LD_SD_down,LD_SD_up,LD_LL_up,SD_SL_up)
# updated_dataset <- rbind(LD_SD_up,LD_SD_down,LL_SL_up,LL_SL_down,WT_LD_SD_up,
#                          WT_LD_SD_down,WT_LL_SL_up,WT_LL_SL_down,XD_LD_up,
#                          XD_LD_down,XD_SD_up,XD_SD_down)
#day_time_condition <- rbind(LL_SL_up, LL_SL_down)

#SD_SL_genes <- unique(SD_SL_highlight_75$ID)
#LD_LL_genes <- unique(LD_LL_highlight_75$ID)

#all_dark_condition_gene <- unique(all_dark_condition$ID)
hc_genes <- unique(hc_genes)
write.table(hc_genes, "hc_genes.txt", row.names = FALSE, quote = F)

#day_time_condition <- unique(day_time_condition$ID)

#combined_genes <- c(SD_SL_genes, LD_LL_genes, LL_SL_genes, LD_SD_genes)
#combined_genes <- unique(combined_genes)

#LL_SL_expression_data <- DataMatrix[rownames(DataMatrix) %in% LL_SL_genes, ]
#LD_SD_expression_data <- DataMatrix[rownames(DataMatrix) %in% LD_SD_genes, ]
#expression_data <- normalized_counts[rownames(normalized_counts) %in% combined_genes, ]

#vu <- normalized_counts[rownames(normalized_counts) %in% all_dark_condition_gene, ]
#day_time_expression_data <- normalized_counts[rownames(normalized_counts) %in% day_time_condition, ]
#updated_dataset_expression_data1 <- normalized_counts[rownames(normalized_counts) %in% updated_dataset_gene, ]
#updated_dataset_expression_data2 <- WT_normalized_counts[rownames(WT_normalized_counts) %in% updated_dataset_gene, ]
#updated_dataset_expression_data <- cbind(updated_dataset_expression_data1,updated_dataset_expression_data2)

#dark_TF_list <- read.table("Dark_TF.list.txt", header = FALSE, stringsAsFactors = FALSE)
#hc_TF_list <- read.table("hc_TF.list.txt", header = FALSE, stringsAsFactors = FALSE)
photoperiod_response_TF <- read.table("TF_ids_167.list.txt", header = FALSE, stringsAsFactors = FALSE)


#hc_TF_list$V1 <- sub("\\.\\d+$", "", hc_TF_list$V1)
photoperiod_response_TF$V1 <- sub("\\.\\d+$", "", photoperiod_response_TF$V1)
#Dark_TF_list$V1 <- sub("\\.\\d+$", "", Dark_TF_list$V1)


#LL_SL_TF_ids <- as.character(LL_SL_TF_list$V1)
#LD_SD_TF_ids <- as.character(LD_SD_TF_list$V1)
#SD_SL_TF_ids <- as.character(SD_SL_TF_list$V1)
#LD_LL_TF_ids <- as.character(LD_LL_TF_list$V1)
#TF_ids <- c(LL_SL_TF_ids,LD_SD_TF_ids,SD_SL_TF_ids,LD_LL_TF_ids)

#TF_list <- rbind(LL_SL_TF_list, LD_SD_TF_list, SD_SL_TF_list, LD_LL_TF_list)
#TF_list <- TF_list[!duplicated(TF_list$V1), ]
#Dark_TF_list <- Dark_TF_list[!duplicated(Dark_TF_list$V1), ]
#hc_TF_list <- hc_TF_list[!duplicated(hc_TF_list$V1), ]
photoperiod_response_TF <- photoperiod_response_TF[!duplicated(photoperiod_response_TF$V1), ]

#TF_ids <- as.character(TF_list$V1)
#Dark_TF_ids <- as.character(Dark_TF_list$V1)
#hc_TF_ids <- as.character(hc_TF_list$V1)
photoperiod_response_TF_ids <- as.character(photoperiod_response_TF$V1)

valid_regulators <- intersect(photoperiod_response_TF_ids, rownames(tpm_hc))
write.table(valid_regulators, file="valid_regulators.txt", 
            quote=FALSE, row.names=FALSE, col.names=FALSE)

missing_regulators <- as.data.frame(setdiff(photoperiod_response_TF_ids, rownames(tpm_hc)))
colnames(missing_regulators) <- c("missing_reguulators_id")
missing_regulators$missing_reguulators_name <- gene_annotations$gene_name[match(missing_regulators$missing_reguulators_id, gene_annotations$ID)]
write.table(missing_regulators, file="missing_regulators.txt", 
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# 构建基因调控网络
#LL_SL_weight_matrix <- GENIE3(as.matrix(LL_SL_expression_data), regulators = LL_SL_TF_ids)
#LD_SD_weight_matrix <- GENIE3(as.matrix(LD_SD_expression_data), regulators = LD_SD_TF_ids)

# 获取有效调控因子
#valid_regulators <- intersect(updated_dataset_TF_ids, rownames(updated_dataset_expression_data))
set.seed(2)
#weight_matrix <- GENIE3(as.matrix(updated_dataset_expression_data), regulators = valid_regulators)

weight_matrix <- GENIE3(as.matrix(tpm_hc),
                        regulators = valid_regulators,
                        K = "sqrt",
                        nTrees = 1000,
                        nCores = 4,    # 根据电脑核心数来调
                        verbose = TRUE)


n_bootstrap <- 100
all_weight_matrices <- list()

for (i in 1:n_bootstrap) {
  cat("Bootstrap run:", i, "\n")
  set.seed(1000 + i)  # 每次不同种子保证多样性
  
  w_mat <- GENIE3(as.matrix(tpm_hc),
                  regulators = valid_regulators,
                  K = "sqrt",
                  nTrees = 1000,
                  nCores = 4,
                  verbose = FALSE)
  
  all_weight_matrices[[i]] <- w_mat
}

# 计算所有权重矩阵的平均矩阵
# 权重矩阵是稀疏矩阵，先转换成普通矩阵
weight_array <- simplify2array(lapply(all_weight_matrices, as.matrix))
mean_weight_matrix <- apply(weight_array, c(1, 2), mean)

library(Matrix)
regulator_names <- rownames(all_weight_matrices[[1]])
target_names <- colnames(all_weight_matrices[[1]])

dimnames(mean_weight_matrix) <- list(regulator_names, target_names)
mean_weight_matrix_sparse <- Matrix(mean_weight_matrix, sparse = TRUE)

# 保存结果
saveRDS(mean_weight_matrix_sparse, "GENIE3_bootstrap_mean_weight_matrix.rds")

# 导出为链接列表
weight_table <- getLinkList(as.matrix(mean_weight_matrix_sparse))
weight_table_filtered <- weight_table[weight_table$weight > 0.005, ]

#write.table(weight_table_filtered, "GENIE3_bootstrap_links.txt", sep="\t", quote=FALSE, row.names=FALSE)

# 提取边列表
#LL_SL_link_list <- getLinkList(LL_SL_weight_matrix, threshold = 0.001)
#LD_SD_link_list <- getLinkList(LD_SD_weight_matrix, threshold = 0.001)
#link_list <- getLinkList(weight_matrix, threshold = 0.001)

#gene_annotations$gene_name[is.na(gene_annotations$gene_name)] <- "unknown"
#na_indices <- which(is.na(gene_annotations$gene_name))
#unique_na <- paste0("unknown_", seq_along(na_indices))
#gene_annotations$gene_name[na_indices] <- unique_na

#LL_SL_link_list$regulatoryGene <- gene_annotations$gene_name[match(LL_SL_link_list$regulatoryGene, gene_annotations$ID)]
#LL_SL_link_list$targetGene <- gene_annotations$gene_name[match(LL_SL_link_list$targetGene, gene_annotations$ID)]

#LD_SD_link_list$regulatoryGene <- gene_annotations$gene_name[match(LD_SD_link_list$regulatoryGene, gene_annotations$ID)]
#LD_SD_link_list$targetGene <- gene_annotations$gene_name[match(LD_SD_link_list$targetGene, gene_annotations$ID)]


# 添加原始的 seq ID 列
weight_table_filtered$regulatorySeqID <- weight_table_filtered$regulatoryGene
weight_table_filtered$targetSeqID <- weight_table_filtered$targetGene

# 用 gene ID 替换 seq ID，同时保留原始 seq ID
weight_table_filtered$regulatoryGene <- gene_annotations$gene_name[match(weight_table_filtered$regulatorySeqID, gene_annotations$ID)]
weight_table_filtered$targetGene <- gene_annotations$gene_name[match(weight_table_filtered$targetSeqID, gene_annotations$ID)]
# 保存结果
#write.csv(LL_SL_link_list, "LL_SL_gene_regulatory_network.csv", row.names = FALSE)
#write.csv(LD_SD_link_list, "LD_SD_gene_regulatory_network.csv", row.names = FALSE)
#write.table(weight_table_filtered, "GENIE3_bootstrap_links_with_seqID.csv", sep = "\t", row.names = FALSE, quote = FALSE)

weight_table_filtered_clean <- weight_table_filtered[!is.na(weight_table_filtered$regulatoryGene) & !is.na(weight_table_filtered$targetGene), ]
write.table(weight_table_filtered_clean, "GENIE3_bootstrap_links_cleaned_with_seqID_2685.csv", sep = "\t", row.names = FALSE, quote = FALSE)

weight_table_filtered_clean <- read.table("GENIE3_bootstrap_links_cleaned_with_seqID_2685_1.csv", sep = '\t', header = TRUE, quote = "")



# 计算 regulator 的 degree
#regulator_name_degree <- table(weight_table_filtered_clean$regulatoryGene)
#regulator_name_degree <- sort(regulator_name_degree, decreasing = TRUE)

regulator_id_degree <- table(weight_table_filtered_clean$regulatorySeqID)
regulator_id_degree <- sort(regulator_id_degree, decreasing = TRUE)

# top 50 regulator
#top50_regulators_names <- as.data.frame(head(regulator_name_degree, 50))
top50_regulators_ids <- as.data.frame(head(regulator_id_degree, 50))
top50_regulators_ids$GeneName <- gene_annotations$gene_name[match(top50_regulators_ids$Var1, gene_annotations$ID)]

# 计算 target gene 的 degree
#target_degree_name <- table(weight_table_filtered_clean$targetGene)
#target_degree_name <- sort(target_degree_name, decreasing = TRUE)

#target_degree_id <- table(weight_table_filtered_clean$targetSeqID)
#target_degree_id <- sort(target_degree_id, decreasing = TRUE)

# top 50 target
#top50_targets_names <- as.data.frame(head(target_degree_name, 50))
#top50_targets_ids <- as.data.frame(head(target_degree_id, 50))
#top50_targets_ids$GeneName <- gene_annotations$gene_name[match(top50_targets_ids$Var1, gene_annotations$ID)]



# 保存成表格
#write.table(as.data.frame(top50_regulators_names), "top50_regulators_name_degree_2685.txt", sep="\t", quote=FALSE, col.names=NA)
write.table(as.data.frame(top50_regulators_ids), "top50_regulators_id_degree_2685_add3.txt", sep="\t", quote=FALSE, col.names=NA)
#write.table(as.data.frame(top50_targets_names), "top50_targets_name_degree_2685.txt", sep="\t", quote=FALSE, col.names=NA)
#write.table(as.data.frame(top50_targets_ids), "top50_targets_is_degree_2685.txt", sep="\t", quote=FALSE, col.names=NA)

#gene_names <- c("BBX19", "BBX24", "PIF1", "BT2", "PRR73", "LHY", "RVE1",
#               "RVE6", "RVE8", "LNK1")
#target <- link_list[link_list$regulatoryGene %in% gene_names | link_list$targetGene %in% gene_names, ]
#write.csv(target, "Dark_targeted_gene_regulatory_network.csv", row.names = FALSE)
#hc_TF_list$gene_name <- gene_annotations$gene_name[match(hc_TF_list$V1, gene_annotations$ID)]
#colnames(hc_TF_list) <- c("SeqID", "Family", "GeneName")
#write.table(hc_TF_list, "TF_list_with_seqID.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#hc_top50 <- read.csv("hc_top50.csv", skip = 1, header = TRUE)
#hc_top50$GeneName <- gene_annotations$gene_name[match(hc_top50$Name, gene_annotations$ID)]

#BiocManager::install("topGO")
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(org.At.tair.db)

# top15 target id
top15_regulatory_ids <- top50_regulators_ids$GeneName[1:15]

for (gene in top15_regulatory_ids) {
  
  # 找所有调控这个regulator的 target genes
  target_ids <- weight_table_filtered_clean$targetGene[weight_table_filtered_clean$regulatoryGene == gene]
  
  # 去重
  target_ids <- as.character(unique(target_ids))
  
  if (length(target_ids) < 10) {
    cat(paste0("Regulatory ", gene, " has less than 10 targets, skipping.\n"))
    next
  }
  
  # GO 富集
  regulatoryID_GENIE3_go <- enrichGO(gene = target_ids,
                                       OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                                       keyType = "SYMBOL",  # 根据基因名称进行富集分析
                                       ont = "BP",  # 生物过程
                                       pvalueCutoff = 0.05,
                                       pAdjustMethod = "BH",
                                       qvalueCutoff = 0.05)
  
  regulatoryID_GENIE3_go_result <- regulatoryID_GENIE3_go@result
  
  # 保存富集结果
  write.table(regulatoryID_GENIE3_go_result,
              file = paste0("GO_results_regulatoryID_", gene, ".txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  df <- regulatoryID_GENIE3_go@result[1:10, ]
  
  # 用ggplot画点图
  g <- ggplot(df, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
    geom_point(aes(size=Count, color=pvalue)) +
    theme_minimal() +
    labs(title = "Top 10 GO terms", x = "Gene Ratio", y = "GO Term", color = "p value",
         size = "Gene Count") +
    scale_color_gradient(low = "blue", high = "red")
  
  ggsave(filename = paste0("GO_results_regulatoryID_", gene, "_dotplot.png"),
         plot = g)
  
}



# 统计up、down和none的个数
ldsd <- read.table("LD_SD.DESeq2.txt", sep = '\t', header = TRUE, quote = "")
llsl <- read.table("LL_SL.DESeq2.txt", sep = '\t', header = TRUE, quote = "")

result_summary <- data.frame(
  sig_type = names(table(ldsd$sig)),
  LD_SD_count = as.vector(table(ldsd$sig)),
  LL_SL_count = as.vector(table(llsl$sig))
)

write.table(result_summary, "DEG_sig_summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)














