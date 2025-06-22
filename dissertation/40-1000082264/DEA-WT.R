library(DESeq2)
setwd("/Users/chenziyin/Downloads/Dissertation")

# 样本名称
sampleNames <- c("WT-LD1", "WT-LD2", "WT-LL1", "WT-LL2", "WT-SD1", "WT-SD2",
                 "WT-SL1", "WT-SL2", "XD1", "XD2", "LD1", "LD2", "SD1", "SD2")
for (name in sampleNames) {
  file_path <- paste0(name, "_counts.txt")
  data <- read.table(file_path, header = TRUE, row.names = 1)
  assign(paste0(name, "_counts"), data)
}
DataMatrix <- cbind(
  get("WT-LD1_counts")[6], get("WT-LD2_counts")[6],
  get("WT-LL1_counts")[6], get("WT-LL2_counts")[6],
  get("WT-SD1_counts")[6], get("WT-SD2_counts")[6],
  get("WT-SL1_counts")[6], get("WT-SL2_counts")[6],
  get("XD1_counts")[6],    get("XD2_counts")[6],
  get("LD1_counts")[6], get("LD2_counts")[6],
  get("SD1_counts")[6], get("SD2_counts")[6]
)
# DataMatrix <- DataMatrix[rowSums(DataMatrix) > 0,]
#LD1 & LD2 <-> LL1 & LL2 
#LD1 & LD2 <-> SD1 & SD2
#LL1 & LL2 <-> SL1 & SL2
#SD1 & SD2 <-> SL1 & SL2
colnames(DataMatrix) <- c("WT_LD1", "WT_LD2", "WT_LL1", "WT_LL2", "WT_SD1",
                          "WT_SD2", "WT_SL1", "WT_SL2", "XD1", "XD2",
                          "LD1", "LD2", "SD1", "SD2")



# 假设非CDS基因的ID存储在non_cds_genes向量中（需提前定义）
non_cds_counts <- counts[rownames(counts) %in% non_cds_genes, ]
cat("非CDS基因数量:", nrow(non_cds_counts), "\n")
cat("非CDS基因计数分布:\n")
print(summary(rowSums(non_cds_counts)))


# 定义组信息
comparisons <- list(
  WT_LD_LL = list(samples = c("WT_LD1", "WT_LD2", "WT_LL1", "WT_LL2"), group = c("WT_LD", "WT_LD", "WT_LL", "WT_LL")),
  WT_LD_SD = list(samples = c("WT_LD1", "WT_LD2", "WT_SD1", "WT_SD2"), group = c("WT_LD", "WT_LD", "WT_SD", "WT_SD")),
  WT_LL_SL = list(samples = c("WT_LL1", "WT_LL2", "WT_SL1", "WT_SL2"), group = c("WT_LL", "WT_LL", "WT_SL", "WT_SL")),
  WT_SD_SL = list(samples = c("WT_SD1", "WT_SD2", "WT_SL1", "WT_SL2"), group = c("WT_SD", "WT_SD", "WT_SL", "WT_SL")),
  XD_LD = list(samples = c("XD1", "XD2", "LD1", "LD2"), group = c("XD", "XD", "LD", "LD")),
  XD_SD = list(samples = c("XD1", "XD2", "SD1", "SD2"), group = c("XD", "XD", "SD", "SD"))
  )

'comparisons_day_dark <- list(
  daytime = list(samples = c("WT-LL1", "WT-LL2", "WT-SL1", "WT-SL2"), group = c("LL", "LL", "SL", "SL")),
  nighttime = list(samples = c("WT-LD1", "WT-LD2", "WT-SD1", "WT-SD2"), group = c("LD", "LD", "SD", "SD"))
)'



library(rtracklayer)

# 读取GFF3文件
gff3_file <- "AeUmbellulata_TA1851_v1.gff3"
gff3_data <- import(gff3_file, format = "gff3")

library(stringr)
# 提取基因ID和注释信息
gene_annotations <- gff3_data[gff3_data$type == "gene", ]
gene_annotations <- as.data.frame(gene_annotations)
# 定义提取基因名称的函数
extract_gene_name <- function(note) {
  match <- str_match(note, "Similar to ([^:]+):")
  return(match[,2])
}
gene_annotations$gene_name <- sapply(gene_annotations$Note, extract_gene_name)
gene_annotations <- gene_annotations[, c("ID", "seqnames", "gene_name", "Ontology_term")]

all_samples <- list()
all_groups <- list()


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
  write.table(res1, 'DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)

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

  #highlight_genes <- rbind(res1_up,res1_down)
  #target <- c('AeUmb.TA1851.r1.4UG0030010',
              #'AeUmb.TA1851.r1.5UG0039190',
              #'AeUmb.TA1851.r1.2UG0045860',
              #'AeUmb.TA1851.r1.6UG0050370')
  #target_genes <- merged_results %>% filter(merged_results$ensembl_gene_id %in% target)

  #默认情况下，横轴展示 log2FoldChange，纵轴展示 -log10 转化后的 padj
  p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.4, size = 2) +
    scale_color_manual(values = c("#ff4757", "#d2dae2", "#546de5"), limits = c('up', 'none', 'down'),
                       labels = c("up", "none", "down")) +
    labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = paste(comparison, 'Comparison'), color = '') +
    theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          legend.key = element_rect(fill = 'transparent')) +
    geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +
    geom_hline(yintercept = 2, lty = 3, color = 'black') +
    xlim(-16.5, 16.5) + ylim(0, 250)
  
  
  # 
  # # 检查当前比较是否是 LD_SD 或 LL_SL
  # if (comparison %in% c("LD_SD", "LL_SL")) {
  #   
  #   # 提取所有的RVE1基因数据，包括seqname信息
  #   RVE1_data <- merged_results[merged_results$gene_name == "RVE1", ]
  #   
  #   if (nrow(RVE1_data) > 0) {
  #     # 为 RVE1 基因添加包含 seqname 的 label
  #     RVE1_data$label <- paste(RVE1_data$gene_name, "(", RVE1_data$seqname, ")")
  #   }
  #   
  #   
  #   # 提取 MADS15 的数据
  #   MADS15_data <- merged_results[tolower(merged_results$gene_name) == "mads15", ]
  #   
  #   highlight_genes <- rbind(highlight_genes, MADS15_data)
  #   
  #   #p <- p + geom_text(data = RVE1_data, aes(label = label), 
  #                      #vjust = 2, hjust = 0.5, check_overlap = TRUE, size = 3, color = "black")
  #   p <- p + geom_point(data = RVE1_data, aes(x = log2FoldChange, y = -log10(padj)), color = "black", size = 2) + 
  #     geom_text_repel(data = RVE1_data, aes(label = label), box.padding = 0.8, point.padding = 0.8, 
  #                     segment.color = 'black', color = "black", size = 3, max.overlaps = 3)
  # }
  # p <- p + geom_text(data = highlight_genes, aes(label = gene_name),
  #                    vjust = 2, hjust = 0.5, check_overlap = TRUE, size = 3)
  # 
  
  p <- p + geom_text(data = highlight_genes, aes(label = highlight_genes$gene_name),
                     vjust = 2, hjust = 0.5, check_overlap = TRUE, size = 3)
  
  
  
  # 显式打印图像并保存
  print(p)
  ggsave(filename = paste0(comparison, '_volcano_plot.png'), plot = p)
  
  # DataMatrix_up <- DataMatrix_sub[rownames(DataMatrix_sub) %in% res1_up$ID, ]
  # DataMatrix_down <- DataMatrix_sub[rownames(DataMatrix_sub) %in% res1_down$ID, ]
  # DataMatrix <- rbind(DataMatrix_up, DataMatrix_down)
  # Names_up <- res1_up$gene_name[match(rownames(DataMatrix_up), res1_up$ID)]
  # Names_down <- res1_down$gene_name[match(rownames(DataMatrix_down), res1_down$ID)]
  # Names <- c(Names_up, Names_down)
  # Names[is.na(Names)] <- "unknown"
  # 
  # unique_groups <- as.data.frame(unique(Names))
  # colnames(unique_groups) <- 'Groups'
  # Names <- make.unique(Names)
  # rownames(DataMatrix) <- Names
  # 
  # group_annotation <- sapply(rownames(DataMatrix), function(x) {
  #   # 检查是否包含唯一名称
  #   match_found <- sapply(unique_groups$Groups, function(y) grepl(paste0("\\b",y,"\\b"),x))
  #   if (any(match_found)) {
  #     return(unique_groups$Groups[which(match_found)[1]])
  #   } else if (grepl("unknown", x)) {
  #     return("Unknown")
  #   } else {
  #     return(x)
  #   }
  }
#)
  
  row_annotation <- data.frame(Group = group_annotation)
  
  
  heatmap <- pheatmap(DataMatrix,
           cluster_rows = TRUE, # Clustering of rows
           cluster_cols = TRUE, # Clustering of cols
           scale="row", # Standardisation with behavioural benchmarks
           main = paste(comparison, 'Comparison'),
           show_rownames = F,
           clustering_method = "ward.D2", # Using the square of the Euclidean distance is less subject to outliers.
           border = F,
           annotation_row = row_annotation, # 添加行注释
           color = colorRampPalette(c("#F9E07F","#FFDCA2", "#F9AD6A", "#FD8F52","#e23e35"))(50),
           cutree_rows = 4,
           fontsize_col = 10,
           angle_col = 45,
           annotation_legend = FALSE,
           treeheight_row = 50)
  
  # 保存热图
  ggsave(filename = paste0(comparison, '_heatmap.png'), plot = heatmap)
  # 将样本数据和组信息添加到列表中
  all_samples[[comparison]] <- DataMatrix_sub
  all_groups[[comparison]] <- colData
}

DataMatrix_subset <- DataMatrix[, 1:8]
colData <- data.frame(
  row.names = colnames(DataMatrix_subset),
  condition = factor(c(rep("WT_LD", 2), rep("WT_LL", 2),
                       rep("WT_SD", 2), rep("WT_SL", 2)))
)
# 使用所有样本数据（DataMatrix 已定义）
dds_all <- DESeqDataSetFromMatrix(
  countData = DataMatrix_subset,
  colData = colData,
  design = ~ condition  # 设计公式单独作为参数
)
vsd_combined <- vst(dds_all, blind=FALSE)
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
  xlab("Conditions") +
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
  xlab("Conditions") +
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
png("WT_PCA_plots.png", width = 15, height = 8, units = "in", res = 300)
grid.arrange(pca_plot_combined, arrangeGrob(pca_plot_PC1, pca_plot_PC2, ncol = 1), ncol = 2)
dev.off()




library(VennDiagram)
# 读取保存的文件
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

WT_LD_LL_genes <- unique(c(WT_LD_LL_up$gene_name, WT_LD_LL_down$gene_name))
WT_LD_SD_genes <- unique(c(WT_LD_SD_up$gene_name, WT_LD_SD_down$gene_name))
WT_LL_SL_genes <- unique(c(WT_LL_SL_up$gene_name, WT_LL_SL_down$gene_name))
WT_SD_SL_genes <- unique(c(WT_SD_SL_up$gene_name, WT_SD_SL_down$gene_name))
XD_LD_genes <- unique(c(XD_LD_up$gene_name, XD_LD_down$gene_name))
XD_SD_genes <- unique(c(XD_SD_up$gene_name, XD_SD_down$gene_name))

WT_LD_SD_up_genes <- unique(WT_LD_SD_up$gene_name)
WT_LD_SD_down_genes <- unique(WT_LD_SD_down$gene_name)
WT_LL_SL_up_genes <- unique(WT_LL_SL_up$gene_name)
WT_LL_SL_down_genes <- unique(WT_LL_SL_down$gene_name)
XD_LD_up_genes <- unique(XD_LD_up$gene_name)
XD_LD_down_genes <- unique(XD_LD_down$gene_name)
XD_SD_up_genes <- unique(XD_SD_up$gene_name)
XD_SD_down_genes <- unique(XD_SD_down$gene_name)

# 只在LD_SD和LL_SL中表达的基因
# WT_LD_SD_up_unique_genes <- setdiff(WT_LD_SD_up_genes, c(WT_LD_LL_genes, WT_LL_SL_genes, WT_SD_SL_genes, WT_LD_SD_down_unique_genes))
# WT_LD_SD_down_unique_genes <- setdiff(WT_LD_SD_down_genes, c(WT_LD_LL_genes, WT_LL_SL_genes, WT_SD_SL_genes, WT_LD_SD_up_unique_genes))
# WT_LL_SL_up_unique_genes <- setdiff(WT_LL_SL_up_genes, c(WT_LD_LL_genes, WT_LD_SD_genes, WT_SD_SL_genes, WT_LL_SL_down_unique_genes))
# WT_LL_SL_down_unique_genes <- setdiff(WT_LL_SL_down_genes, c(WT_LD_LL_genes, WT_LD_SD_genes, WT_SD_SL_genes, WT_LL_SL_up_unique_genes))
# 
# write.table(WT_LD_SD_up_unique_genes, file = 'Venn_WT_LD_SD_up_unique_genes.txt',
#             sep = '\t', col.names = NA, quote = FALSE)
# write.table(WT_LD_SD_down_unique_genes, file = 'Venn_WT_LD_SD_down_unique_genes.txt',
#             sep = '\t', col.names = NA, quote = FALSE)
# write.table(WT_LL_SL_up_unique_genes, file = 'Venn_WT_LL_SL_up_unique_genes.txt',
#             sep = '\t', col.names = NA, quote = FALSE)
# write.table(WT_LL_SL_down_unique_genes, file = 'Venn_WT_LL_SL_down_unique_genes.txt',
#             sep = '\t', col.names = NA, quote = FALSE)

gene_list <- list(
  "LD vs LL" = WT_LD_LL_genes,
  "LD vs SD" = WT_LD_SD_genes,
  "LL vs SL" = WT_LL_SL_genes,
  "SD vs SL" = WT_SD_SL_genes
)

gene_list_1 <- list(
  "WT_LD_SD Up" = WT_LD_SD_up_genes,
  "WT_LD_SD Down" = WT_LD_SD_down_genes,
  "WT_LL_SL Up" = WT_LL_SL_up_genes,
  "WT_LL_SL Down" = WT_LL_SL_down_genes
)

gene_list_2 <- list(
  "WT_LD_SD" = WT_LD_SD_genes,
  "WT_LL_SL" = WT_LL_SL_genes
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

WT_LD_SD_genes <- unique(LD_SD_highlight_75$ID)
WT_LL_SL_genes <- unique(LL_SL_highlight_75$ID)

# 保存为 CSV 文件
write.csv(WT_LD_SD_genes, "WT_LD_SD_genes.csv", row.names = FALSE)
write.csv(WT_LL_SL_genes, "WT_LL_SL_genes.csv", row.names = FALSE)

write.csv(unique_groups, "unique_groups.csv", row.names = FALSE)


library(clusterProfiler)
library(dplyr)
library(tidyr)


# 提取基因ID（确保是字符向量）
up_genes <- as.character(WT_LL_SL_up$gene_name)  # 使用基因名称而不是基因ID
down_genes <- as.character(WT_LL_SL_down$gene_name)  # 使用基因名称而不是基因ID
up_genes <- up_genes[!is.na(up_genes)]
down_genes <- down_genes[!is.na(down_genes)]

# up_genes <- as.character(WT_SD_SL_up$gene_name)  # 使用基因名称而不是基因ID
# down_genes <- as.character(WT_SD_SL_down$gene_name)  # 使用基因名称而不是基因ID
# up_genes <- up_genes[!is.na(up_genes)]
# down_genes <- down_genes[!is.na(down_genes)]

#total_gene <- c(up_genes,down_genes)
library(org.At.tair.db)
go_up <- enrichGO(gene = up_genes,
           OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
           keyType = "SYMBOL",  # 根据基因名称进行富集分析
           ont = "BP",  # 生物过程
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           qvalueCutoff = 0.05)
dim(go_up)
WT_LL_SL_go_up <- go_up@result
write.csv(WT_LL_SL_go_up, "WT_LL_SL_go_up.csv", row.names = FALSE)
# WT_SD_SL_go_up <- go_up@result
# write.csv(WT_SD_SL_go_up, "WT_SD_SL_go_up.csv", row.names = FALSE)


go_down <- enrichGO(gene = down_genes,
                  OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                  keyType = "SYMBOL",  # 根据基因名称进行富集分析
                  ont = "BP",  # 生物过程
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "none",
                  qvalueCutoff = 1)
dim(go_down)
WT_LL_SL_go_down <- go_down@result
write.csv(WT_LL_SL_go_down, "WT_LL_SL_go_down.csv", row.names = FALSE)

# WT_SD_SL_go_down <- go_down@result
# write.csv(WT_SD_SL_go_down, "WT_SD_SL_go_down.csv", row.names = FALSE)


#结果可视化 
WT_LL_SL_up_bar <- barplot(go_up,showCategory=15,drop=T) 
WT_LL_SL_up_bar <- WT_LL_SL_up_bar + theme(
  axis.text.y = element_text(size = 10)) + # 调整 y 轴标签的字体大小
  ggtitle("WT_LL_SL_upregulated_Bar")


WT_LL_SL_up_bar
ggsave(filename = "WT_LL_SL_GO_Enrichment_Upregulated_Bar_Font.png", plot = WT_LL_SL_up_bar)
WT_LL_SL_up_dot <- dotplot(go_up,showCategory=15)
WT_LL_SL_up_dot <- WT_LL_SL_up_dot + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("WT_LL_SL_upregulated_Dot")

WT_LL_SL_up_dot
ggsave(filename = "WT_LL_SL_GO_Enrichment_Upregulated_Dot_Font.png", plot = WT_LL_SL_up_dot)



WT_LL_SL_down_bar <- barplot(go_down,showCategory=15,drop=T) 
WT_LL_SL_down_bar <- WT_LL_SL_down_bar + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("WT_LL_SL_downregulated_Bar")

WT_LL_SL_down_bar
ggsave(filename = "WT_LL_SL_GO_Enrichment_Downregulated_Bar_Font.png", plot = WT_LL_SL_down_bar)
WT_LL_SL_down_dot <- dotplot(go_down,showCategory=15)
WT_LL_SL_down_dot <- WT_LL_SL_down_dot + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("WT_LL_SL_downregulated_Dot")

WT_LL_SL_down_dot
ggsave(filename = "WT_LL_SL_GO_Enrichment_Downregulated_Dot_Font.png", plot = WT_LL_SL_down_dot)





# 提取基因ID（确保是字符向量）
up_genes <- as.character(WT_LD_SD_up$gene_name)  # 使用基因名称而不是基因ID
down_genes <- as.character(WT_LD_SD_down$gene_name)  # 使用基因名称而不是基因ID
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
WT_LD_SD_go_up <- go_up@result
write.csv(WT_LD_SD_go_up, "WT_LD_SD_go_up.csv", row.names = FALSE)


go_down <- enrichGO(gene = down_genes,
                    OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                    keyType = "SYMBOL",  # 根据基因名称进行富集分析
                    ont = "BP",  # 生物过程
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "none",
                    qvalueCutoff = 1)

dim(go_down)
WT_LD_SD_go_down <- go_down@result
write.csv(WT_LD_SD_go_down, "WT_LD_SD_go_down.csv", row.names = FALSE)


#结果可视化 
WT_LD_SD_up_bar <- barplot(go_up,showCategory=15,drop=T) 
WT_LD_SD_up_bar <- WT_LD_SD_up_bar + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("WT_LD_SD_upregulated_Bar")

WT_LD_SD_up_bar
ggsave(filename = "WT_LD_SD_GO_Enrichment_Upregulated_Bar_Font.png", plot = WT_LD_SD_up_bar)
WT_LD_SD_up_dot <- dotplot(go_up,showCategory=15)
WT_LD_SD_up_dot <- WT_LD_SD_up_dot + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("WT_LD_SD_upregulated_Dot")

WT_LD_SD_up_dot
ggsave(filename = "WT_LD_SD_GO_Enrichment_Upregulated_Dot_Font.png", plot = WT_LD_SD_up_dot)


WT_LD_SD_down_bar <- barplot(go_down,showCategory=15,drop=T) 
WT_LD_SD_down_bar <- WT_LD_SD_down_bar + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("WT_LD_SD_downregulated_Bar")

WT_LD_SD_down_bar
ggsave(filename = "WT_LD_SD_GO_Enrichment_Downregulated_Bar_Font.png", plot = WT_LD_SD_down_bar)
WT_LD_SD_down_dot <- dotplot(go_down,showCategory=15)
WT_LD_SD_down_dot <- WT_LD_SD_down_dot + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("WT_LD_SD_downregulated_Dot")

WT_LD_SD_down_dot
ggsave(filename = "WT_LD_SD_GO_Enrichment_Downregulated_Dot_Font.png", plot = WT_LD_SD_down_dot)






# 提取基因ID（确保是字符向量）
up_genes <- as.character(XD_SD_up$gene_name)  # 使用基因名称而不是基因ID
down_genes <- as.character(XD_SD_down$gene_name)  # 使用基因名称而不是基因ID
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
XD_SD_go_up <- go_up@result
write.csv(XD_SD_go_up, "XD_SD_go_up.csv", row.names = FALSE)


go_down <- enrichGO(gene = down_genes,
                    OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                    keyType = "SYMBOL",  # 根据基因名称进行富集分析
                    ont = "BP",  # 生物过程
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "none",
                    qvalueCutoff = 1)

dim(go_down)
XD_SD_go_down <- go_down@result
write.csv(XD_SD_go_down, "XD_SD_go_down.csv", row.names = FALSE)


#结果可视化 
XD_SD_go_up_bar <- barplot(go_up,showCategory=15,drop=T) 
XD_SD_go_up_bar <- XD_SD_go_up_bar + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("XD_SD_upregulated_Bar")

XD_SD_go_up_bar
ggsave(filename = "XD_SD_GO_Enrichment_Upregulated_Bar_Font.png", plot = XD_SD_go_up_bar)
XD_SD_up_dot <- dotplot(go_up,showCategory=15)
XD_SD_up_dot <- XD_SD_up_dot + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("XD_SD_upregulated_Dot")

XD_SD_up_dot
ggsave(filename = "XD_SD_Enrichment_Upregulated_Dot_Font.png", plot = XD_SD_up_dot)


XD_SD_down_bar <- barplot(go_down,showCategory=15,drop=T) 
XD_SD_down_bar <- XD_SD_down_bar + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("XD_SD_downregulated_Bar")

XD_SD_down_bar
ggsave(filename = "XD_SD_GO_Enrichment_Downregulated_Bar_Font.png", plot = XD_SD_down_bar)
XD_SD_down_dot <- dotplot(go_down,showCategory=15)
XD_SD_down_dot <- XD_SD_down_dot + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("XD_SD_downregulated_Dot")

XD_SD_down_dot
ggsave(filename = "XD_SD_GO_Enrichment_Downregulated_Dot_Font.png", plot = XD_SD_down_dot)



# 提取基因ID（确保是字符向量）
up_genes <- as.character(XD_LD_up$gene_name)  # 使用基因名称而不是基因ID
down_genes <- as.character(XD_LD_down$gene_name)  # 使用基因名称而不是基因ID
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
XD_LD_go_up <- go_up@result
write.csv(XD_SD_go_up, "XD_LD_go_up.csv", row.names = FALSE)


go_down <- enrichGO(gene = down_genes,
                    OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                    keyType = "SYMBOL",  # 根据基因名称进行富集分析
                    ont = "BP",  # 生物过程
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "none",
                    qvalueCutoff = 1)

dim(go_down)
XD_LD_go_down <- go_down@result
write.csv(XD_LD_go_down, "XD_LD_go_down.csv", row.names = FALSE)


#结果可视化 
XD_LD_go_up_bar <- barplot(go_up,showCategory=15,drop=T) 
XD_LD_go_up_bar <- XD_LD_go_up_bar + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("XD_LD_upregulated_Bar")

XD_LD_go_up_bar
ggsave(filename = "XD_LD_GO_Enrichment_Upregulated_Bar_Font.png", plot = XD_LD_go_up_bar)
XD_LD_up_dot <- dotplot(go_up,showCategory=15)
XD_LD_up_dot <- XD_LD_up_dot + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("XD_LD_upregulated_Dot")

XD_LD_up_dot
ggsave(filename = "XD_LD_Enrichment_Upregulated_Dot_Font.png", plot = XD_LD_up_dot)


XD_LD_down_bar <- barplot(go_down,showCategory=15,drop=T) 
XD_LD_down_bar <- XD_LD_down_bar + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("XD_LD_downregulated_Bar")

XD_LD_down_bar
ggsave(filename = "XD_LD_GO_Enrichment_Downregulated_Bar_Font.png", plot = XD_LD_down_bar)
XD_LD_down_dot <- dotplot(go_down,showCategory=15)
XD_LD_down_dot <- XD_LD_down_dot + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("XD_LD_downregulated_Dot")

XD_LD_down_dot
ggsave(filename = "XD_LD_GO_Enrichment_Downregulated_Dot_Font.png", plot = XD_LD_down_dot)




# 创建一个空的向量来存储每个样本的组信息
group_info <- c()
sampleNames_standard <- gsub("-", "_", sampleNames)

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

# 按样本顺序排列 group_info
group_info <- group_info[sampleNames_standard]

# 创建 colData 数据框，包含所有样本和对应的组信息
colData <- data.frame(
  row.names = sampleNames_standard,   # 样本名
  condition = factor(group_info)   # 组信息，转化为因子（factor）
)
rownames(colData) <- gsub("-", "_", rownames(colData))

# 创建 DESeqDataSet
dds_full <- DESeqDataSetFromMatrix(
  countData = DataMatrix,   # 你的原始表达矩阵
  colData = colData,        # 你的样本信息（包括组信息）
  design = ~condition       # 设计公式，基于组信息
)

# 运行 DESeq2 分析
dds_full <- DESeq(dds_full)

# 获取归一化后的表达矩阵
normalized_counts <- counts(dds_full, normalized = TRUE)

# 保存总归一化矩阵
write.table(normalized_counts, "WT_Normalized_Expression_Total.txt", sep = "\t", quote = FALSE, col.names = NA)










#BiocManager::install("GENIE3")
library(GENIE3)
WT_SD_SL_highlight_75 <- read.table("WT_SD_SL.highlight.txt", sep = '\t', header = TRUE)
WT_LD_LL_highlight_75 <- read.table("WT_LD_LL.highlight.txt", sep = '\t', header = TRUE)

all_dark_condition <- rbind(WT_LD_SD_down,WT_LD_SD_up,WT_LD_LL_up,WT_SD_SL_up)
day_time_condition <- rbind(WT_LL_SL_up, WT_LL_SL_down)

WT_SD_SL_genes <- unique(WT_SD_SL_highlight_75$ID)
WT_LD_LL_genes <- unique(WT_LD_LL_highlight_75$ID)

all_dark_condition_gene <- unique(all_dark_condition$ID)
write.table(all_dark_condition_gene, "WT_all_dark_condition_gene.txt", row.names = FALSE, quote = F)

day_time_condition <- unique(day_time_condition$ID)

combined_genes <- c(WT_SD_SL_genes, WT_LD_LL_genes, WT_LL_SL_genes, WT_LD_SD_genes)
combined_genes <- unique(combined_genes)

#LL_SL_expression_data <- DataMatrix[rownames(DataMatrix) %in% WT_LL_SL_genes, ]
#LD_SD_expression_data <- DataMatrix[rownames(DataMatrix) %in% WT_LD_SD_genes, ]
expression_data <- normalized_counts[rownames(normalized_counts) %in% combined_genes, ]

all_dark_expression_data <- normalized_counts[rownames(normalized_counts) %in% all_dark_condition_gene, ]
# day_time_expression_data <- normalized_counts[rownames(normalized_counts) %in% day_time_condition, ]

# LL_SL_TF_list <- read.table("LL_SL_TF.list.txt", header = FALSE, stringsAsFactors = FALSE)
# LD_SD_TF_list <- read.table("LD_SD_TF.list.txt", header = FALSE, stringsAsFactors = FALSE)
# SD_SL_TF_list <- read.table("SD_SL_TF.list.txt", header = FALSE, stringsAsFactors = FALSE)
# LD_LL_TF_list <- read.table("LD_LL_TF.list.txt", header = FALSE, stringsAsFactors = FALSE)
Dark_TF_list <- read.table("WT_TF.list.txt", header = FALSE, stringsAsFactors = FALSE)

# LL_SL_TF_list$V1 <- sub("\\.\\d+$", "", LL_SL_TF_list$V1)
# LD_SD_TF_list$V1 <- sub("\\.\\d+$", "", LD_SD_TF_list$V1)
# SD_SL_TF_list$V1 <- sub("\\.\\d+$", "", SD_SL_TF_list$V1)
# LD_LL_TF_list$V1 <- sub("\\.\\d+$", "", LD_LL_TF_list$V1)
Dark_TF_list$V1 <- sub("\\.\\d+$", "", Dark_TF_list$V1)

# daytime_TF_list <- LL_SL_TF_list

#LL_SL_TF_ids <- as.character(LL_SL_TF_list$V1)
#LD_SD_TF_ids <- as.character(LD_SD_TF_list$V1)
#SD_SL_TF_ids <- as.character(SD_SL_TF_list$V1)
#LD_LL_TF_ids <- as.character(LD_LL_TF_list$V1)
#TF_ids <- c(LL_SL_TF_ids,LD_SD_TF_ids,SD_SL_TF_ids,LD_LL_TF_ids)

# TF_list <- rbind(LL_SL_TF_list, LD_SD_TF_list, SD_SL_TF_list, LD_LL_TF_list)
# TF_list <- TF_list[!duplicated(TF_list$V1), ]
Dark_TF_list <- Dark_TF_list[!duplicated(Dark_TF_list$V1), ]
# daytime_TF_list <- daytime_TF_list[!duplicated(daytime_TF_list$V1),]

# TF_ids <- as.character(TF_list$V1)
Dark_TF_ids <- as.character(Dark_TF_list$V1)
# daytime_TF_ids <- as.character(daytime_TF_list$V1)

# 构建基因调控网络
#LL_SL_weight_matrix <- GENIE3(as.matrix(LL_SL_expression_data), regulators = LL_SL_TF_ids)
#LD_SD_weight_matrix <- GENIE3(as.matrix(LD_SD_expression_data), regulators = LD_SD_TF_ids)

# 获取有效调控因子
valid_regulators <- intersect(Dark_TF_ids, rownames(all_dark_expression_data))
set.seed(2)
weight_matrix <- GENIE3(as.matrix(all_dark_expression_data), regulators = valid_regulators)

# 提取边列表
#LL_SL_link_list <- getLinkList(LL_SL_weight_matrix, threshold = 0.001)
#LD_SD_link_list <- getLinkList(LD_SD_weight_matrix, threshold = 0.001)
link_list <- getLinkList(weight_matrix, threshold = 0.001)

#gene_annotations$gene_name[is.na(gene_annotations$gene_name)] <- "unknown"
na_indices <- which(is.na(gene_annotations$gene_name))
unique_na <- paste0("unknown_", seq_along(na_indices))
gene_annotations$gene_name[na_indices] <- unique_na

#LL_SL_link_list$regulatoryGene <- gene_annotations$gene_name[match(LL_SL_link_list$regulatoryGene, gene_annotations$ID)]
#LL_SL_link_list$targetGene <- gene_annotations$gene_name[match(LL_SL_link_list$targetGene, gene_annotations$ID)]

#LD_SD_link_list$regulatoryGene <- gene_annotations$gene_name[match(LD_SD_link_list$regulatoryGene, gene_annotations$ID)]
#LD_SD_link_list$targetGene <- gene_annotations$gene_name[match(LD_SD_link_list$targetGene, gene_annotations$ID)]


# 添加原始的 seq ID 列
link_list$regulatorySeqID <- link_list$regulatoryGene
link_list$targetSeqID <- link_list$targetGene

# 用 gene ID 替换 seq ID，同时保留原始 seq ID
link_list$regulatoryGene <- gene_annotations$gene_name[match(link_list$regulatorySeqID, gene_annotations$ID)]
link_list$targetGene <- gene_annotations$gene_name[match(link_list$targetSeqID, gene_annotations$ID)]
# 保存结果
#write.csv(LL_SL_link_list, "LL_SL_gene_regulatory_network.csv", row.names = FALSE)
#write.csv(LD_SD_link_list, "LD_SD_gene_regulatory_network.csv", row.names = FALSE)
write.table(link_list, "WT_Dark_gene_regulatory_network_with_seqID.csv", sep = "\t", row.names = FALSE, quote = FALSE)
# gene_names <- c("BBX19", "BBX24", "PIF1", "BT2", "PRR73", "LHY", "RVE1",
#                "RVE6", "RVE8", "LNK1")
# target <- link_list[link_list$regulatoryGene %in% gene_names | link_list$targetGene %in% gene_names, ]
# write.csv(target, "Dark_targeted_gene_regulatory_network.csv", row.names = FALSE)


library(magick)
# 读取图像
img1 <- image_read("SD_SL_volcano_plot.png")
img2 <- image_read("LL_SL_volcano_plot.png")
img3 <- image_read("LD_SD_volcano_plot.png")
img4 <- image_read("LD_LL_volcano_plot.png")

img5 <- image_read("SD_SL_heatmap.png")
img6 <- image_read("LL_SL_heatmap.png")
img7 <- image_read("LD_SD_heatmap.png")
img8 <- image_read("LD_LL_heatmap.png")

img9 <- image_read("venn_diagram.png")

# 拼接图像，4行2列
combined_image <- image_append(c(image_append(c(img1, img2), stack = TRUE), 
                                 image_append(c(img3, img4), stack = TRUE)))

# 显示组合图像
print(combined_image)
image_write(combined_image, path = "combined_volcano_plots.png")

combined_image <- image_append(c(image_append(c(img5, img6), stack = TRUE), 
                                 image_append(c(img7, img8), stack = TRUE)))

# 显示组合图像
print(combined_image)
image_write(combined_image, path = "combined_heatmap_plots.png")

library(cowplot)
# 使用ggplot添加标签
final_plot <- ggdraw() + 
  draw_image("combined_volcano_plots.png", x = 0, y = 0, width = 1, height = 1) +
  draw_label("(A)", x = 0.05, y = 1, hjust = 1, vjust = 1, size = 15, fontface = 'bold') +
  draw_label("(B)", x = 0.55, y = 1, hjust = 1, vjust = 1, size = 15, fontface = 'bold') +
  draw_label("(C)", x = 0.05, y = 0.75, hjust = 1, vjust = 1, size = 15, fontface = 'bold') +
  draw_label("(D)", x = 0.55, y = 0.75, hjust = 1, vjust = 1, size = 15, fontface = 'bold')

# 显示并保存最终图像
print(final_plot)
ggsave(filename = "final_combined_figure.png", plot = final_plot, width = 10, height = 15)



img11 <- image_read("LD_SD_GO_Enrichment_Downregulated_Dot_Font.png")
img12 <- image_read("LD_SD_GO_Enrichment_Downregulated_Bar_Font.png")
img13 <- image_read("LD_SD_GO_Enrichment_Upregulated_Dot_Font.png")
img14 <- image_read("LD_SD_GO_Enrichment_Upregulated_Bar_Font.png")

img15 <- image_read("LL_SL_GO_Enrichment_Downregulated_Dot_Font.png")
img16 <- image_read("LL_SL_GO_Enrichment_Downregulated_Bar_Font.png")
img17 <- image_read("LL_SL_GO_Enrichment_Upregulated_Dot_Font.png")
img18 <- image_read("LL_SL_GO_Enrichment_Upregulated_Bar_Font.png")

combined_image1 <- image_append(c(image_append(c(img11, img13, img15, img17), stack = TRUE), 
                                image_append(c(img12, img14, img16, img18), stack = TRUE)))
image_write(combined_image1, path = "combined_GO_plots.png")

final_plot <- ggdraw() + 
  draw_image("combined_GO_plots.png", x = 0, y = 0, width = 1, height = 1) +
  draw_label("(A)", x = 0.03, y = 0.85, hjust = 1, vjust = 1, size = 15, fontface = 'bold') +
  draw_label("(B)", x = 0.35, y = 0.85, hjust = 1, vjust = 1, size = 15, fontface = 'bold') +
  draw_label("(C)", x = 0.03, y = 0.7, hjust = 1, vjust = 1, size = 15, fontface = 'bold') +
  draw_label("(D)", x = 0.35, y = 0.7, hjust = 1, vjust = 1, size = 15, fontface = 'bold') +
  draw_label("(E)", x = 0.03, y = 0.65, hjust = 1, vjust = 1, size = 15, fontface = 'bold') +
  draw_label("(F)", x = 0.35, y = 0.65, hjust = 1, vjust = 1, size = 15, fontface = 'bold') +
  draw_label("(G)", x = 0.03, y = 0.5, hjust = 1, vjust = 1, size = 15, fontface = 'bold') +
  draw_label("(H)", x = 0.35, y = 0.5, hjust = 1, vjust = 1, size = 15, fontface = 'bold')
ggsave(filename = "final_GO_combined_figure.png", plot = final_plot, width = 10, height = 15)





# Venn daytime gene vs night time gene
# 初始化列表来存储显著基因
significant_gene_ids <- list()
significant_gene_ids_up_down <- list()

for (comparison in names(comparisons_day_dark)) {
  samples <- comparisons_day_dark[[comparison]]$samples
  group <- comparisons_day_dark[[comparison]]$group
  
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
  
  # 筛选显著基因（上调或下调）
  significant_gene_ids[[comparison]] <- rownames(subset(res1, sig %in% c('up', 'down')))
  
  # 分别存储上调和下调基因
  significant_gene_ids_up_down[[paste0(comparison, "_up")]] <- rownames(subset(res1, sig == 'up'))
  significant_gene_ids_up_down[[paste0(comparison, "_down")]] <- rownames(subset(res1, sig == 'down'))
}

# 重命名组别为 Daytime Genes 和 Nighttime Genes
names(significant_gene_ids) <- c("Daytime Genes", "Nighttime Genes")

# 重命名上下调组别
names(significant_gene_ids_up_down) <- c(
  "Daytime Upregulated", "Daytime Downregulated",
  "Nighttime Upregulated", "Nighttime Downregulated"
)

# 绘制韦恩图1：两个圆，表示 Daytime 和 Nighttime 显著基因

venn.diagram(
  x = gene_list_2,
  category.names = c("Nighttime Genes", "Daytime Genes"),
  filename = "venn_diagram_day_night.png",
  main = "Significant Genes (Daytime vs. Nighttime)",
  fill = c("#D72E25", "#A6D96A"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("#D72E25", "#A6D96A"),
  cat.pos = c(180, 180),      # 将标签移到圆圈的正上方
  cat.dist = c(-0.42, -0.35)
)


# 绘制韦恩图2：四个圆，表示 Daytime 和 Nighttime 上下调基因
venn_plot <- venn.diagram(
  x = gene_list_1,
  category.names = c("Nighttime Up", "Nighttime Down",
                     "Daytime Up","Daytime Down"),
  #filename = "venn_diagram_up_down(day_night).png",
  filename = NULL,
  main = "Upregulated and Downregulated Genes (Daytime vs. Nighttime)",
  fill = c("#D72E25", "#F46D43", "#1A9850", "#A6D96A"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("#D72E25", "#F46D43", "#1A9850", "#A6D96A"),
  cat.pos = c(-25, 23, 0, 160),  # 调整标签位置以避免重叠
  cat.dist = c(0.23, 0.23, 0.1, -0.1), # 控制标签距离圆圈的距离
  #margin = 0.15
)
# 增加横向边距
grid.newpage()  # 新建页面
pushViewport(viewport(width = 0.8, height = 1))  # 调整绘图区域
grid.draw(venn_plot)  # 绘制 Venn 图

# 保存为文件
png("WT_venn_diagram_up_down(day_night).png", width = 2000, height = 1500, res = 300)
grid.newpage()
pushViewport(viewport(width = 0.8, height = 1))  # 再次设置横向边距
grid.draw(venn_plot)
dev.off()



# 统计每个比较中上调和下调基因的数量
gene_counts <- data.frame(
  Comparison = c("LD vs SD", "LD vs SD", "LL vs SL", "LL vs SL"),
  Regulation = c("Upregulated", "Downregulated", "Upregulated", "Downregulated"),
  Count = c(length(WT_LD_SD_up_genes), 
            length(WT_LD_SD_down_genes), 
            length(WT_LL_SL_up_genes), 
            length(WT_LL_SL_down_genes))
)


nighttime_17 <- setdiff(
  intersect(WT_LD_SD_up_genes, WT_LD_SD_down_genes),
  union(WT_LL_SL_up_genes, WT_LL_SL_down_genes)
)

day_night_13 <- setdiff(
  Reduce(intersect, list(WT_LL_SL_up_genes, WT_LD_SD_up_genes, WT_LD_SD_down_genes)),
  union(WT_LL_SL_down_genes, setdiff(WT_LD_SD_up_genes, WT_LL_SL_up_genes))
)

#57
day_night_12 <- setdiff(
  intersect(WT_LL_SL_up_genes, WT_LD_SD_down_genes),
  union(WT_LD_SD_up_genes, WT_LL_SL_down_genes)
)

#754
WT_nighttime_up_only_227 <- setdiff(WT_LD_SD_up_genes,
                                 union(WT_LD_SD_down_genes, union(WT_LL_SL_up_genes, WT_LL_SL_down_genes)))
#851
WT_nighttime_down_only_485 <- setdiff(WT_LD_SD_down_genes,
                                   union(WT_LD_SD_up_genes, union(WT_LL_SL_up_genes, WT_LL_SL_down_genes)))
#297
WT_daytime_up_only_475 <- setdiff(WT_LL_SL_up_genes,
                           union(WT_LL_SL_down_genes, union(WT_LD_SD_up_genes, WT_LD_SD_down_genes)))
#532
WT_daytime_down_only_295 <- setdiff(WT_LL_SL_down_genes, union(WT_LL_SL_up_genes,
                                                     union(WT_LD_SD_up_genes, WT_LD_SD_down_genes)))
daytime_up_down_only_7 <- setdiff(
  intersect(WT_LL_SL_up_genes, WT_LL_SL_down_genes),
  union(WT_LD_SD_up_genes, WT_LD_SD_down_genes)
)

#189
nighttime_up_daytime_up_142 <- setdiff(
  intersect(WT_LD_SD_up_genes, WT_LL_SL_up_genes),
  union(WT_LD_SD_down_genes, WT_LL_SL_down_genes)
)

#149
nighttime_down_daytime_down_108 <- setdiff(
  intersect(WT_LD_SD_down_genes, WT_LL_SL_down_genes),
  union(WT_LD_SD_up_genes, WT_LL_SL_up_genes)
)

#121
nighttime_up_daytime_down_109 <- setdiff(
  intersect(WT_LD_SD_up_genes, WT_LL_SL_down_genes),
  union(WT_LD_SD_down_genes, WT_LL_SL_up_genes)
)

nighttime_up_down_daytime_down_8 <- setdiff(
  Reduce(intersect, list(WT_LD_SD_up_genes, WT_LD_SD_down_genes, WT_LL_SL_down_genes)),
  WT_LL_SL_up_genes
)
daytime_up_down_nighttime_up_5 <- setdiff(
  Reduce(intersect, list(WT_LL_SL_up_genes, WT_LL_SL_down_genes, WT_LD_SD_up_genes)),
  WT_LD_SD_down_genes
)
daytime_up_down_nighttime_down_2 <- setdiff(
  Reduce(intersect, list(WT_LL_SL_up_genes, WT_LL_SL_down_genes, WT_LD_SD_down_genes)),
  WT_LD_SD_up_genes
)
all_groups_8 <- Reduce(intersect, list(WT_LD_SD_up_genes, WT_LD_SD_down_genes,
                                     WT_LL_SL_up_genes, WT_LL_SL_down_genes))

WT_LD_SD_up$condition <- "nighttime"
WT_LD_SD_down$condition <- "nighttime"
WT_LL_SL_up$condition <- "daytime"
WT_LL_SL_down$condition <- "daytime"

all_data <- rbind(WT_LD_SD_up, WT_LD_SD_down, WT_LL_SL_up, WT_LL_SL_down)
matched_genes_17 <- all_data[all_data$gene_name %in% nighttime_17, ]
matched_genes_13 <- all_data[all_data$gene_name %in% day_night_13, ]
matched_genes_12 <- all_data[all_data$gene_name %in% day_night_12, ]

WT_matched_gene_227 <- all_data[all_data$gene_name %in% WT_nighttime_up_only_227,]
WT_matched_gene_485 <- all_data[all_data$gene_name %in% WT_nighttime_down_only_485,]
WT_matched_gene_475 <- all_data[all_data$gene_name %in% WT_daytime_up_only_475,]
WT_matched_gene_295 <- all_data[all_data$gene_name %in% WT_daytime_down_only_295,]
matched_gene_7 <- all_data[all_data$gene_name %in% daytime_up_down_only_7,]
matched_gene_142 <- all_data[all_data$gene_name %in% nighttime_up_daytime_up_142,]
matched_gene_108 <- all_data[all_data$gene_name %in% nighttime_down_daytime_down_108,]
matched_gene_109 <- all_data[all_data$gene_name %in% nighttime_up_daytime_down_109,]
matched_gene_8 <- all_data[all_data$gene_name %in% nighttime_up_down_daytime_down_8,]
matched_gene_5 <- all_data[all_data$gene_name %in% daytime_up_down_nighttime_up_5,]
matched_gene_2 <- all_data[all_data$gene_name %in% daytime_up_down_nighttime_down_2,]
matched_gene_8.1 <- all_data[all_data$gene_name %in% all_groups_8,]

write.csv(matched_genes_12, file = "matched_genes_WT_daytime_up_nighttime_down_12.csv", row.names = FALSE)
write.csv(WT_matched_gene_485, file = "matched_genes_WT_nighttime_down_only_485.csv", row.names = FALSE)
write.csv(WT_matched_gene_475, file = "matched_genes_WT_daytime_up_only_475.csv", row.names = FALSE)
write.csv(WT_matched_gene_295, file = "matched_genes_WT_daytime_down_only_295.csv", row.names = FALSE)



write.csv(matched_gene_596, "nighttime_up_only_596.csv", row.names = FALSE)
write.csv(matched_gene_554, "nighttime_down_only_554.csv", row.names = FALSE)
write.csv(matched_gene_219, "daytime_up_only_219.csv", row.names = FALSE)
write.csv(matched_gene_381, "daytime_down_only_381.csv", row.names = FALSE)
write.csv(matched_gene_7, "daytime_up_down_only_7.csv", row.names = FALSE)
write.csv(matched_gene_142, "nighttime_up_daytime_up_142.csv", row.names = FALSE)
write.csv(matched_gene_108, "nighttime_down_daytime_down_108.csv", row.names = FALSE)
write.csv(matched_gene_109, "nighttime_up_daytime_down_109.csv", row.names = FALSE)
write.csv(matched_gene_8, "nighttime_up_down_daytime_down_8.csv", row.names = FALSE)
write.csv(matched_gene_5, "daytime_up_down_nighttime_up_5.csv", row.names = FALSE)
write.csv(matched_gene_2, "daytime_up_down_nighttime_down_2.csv", row.names = FALSE)
write.csv(matched_gene_8.1, "all_groups_8.csv", row.names = FALSE)




# 绘制 bar plot
bar_plot <- ggplot(gene_counts, aes(x = Comparison, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.9), # 对齐文本
            vjust = -0.5, # 调整标签位置
            size = 4) +   # 设置标签字体大小
  labs(title = "Differentially Expressed Genes", 
       x = "Comparison", 
       y = "Number of Genes") +
  scale_fill_manual(values = c("Upregulated" = "#F46D43", "Downregulated" = "#A6D96A")) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.title = element_blank(),
    panel.grid.major = element_line(color = "grey80"),  # 添加主要网格线
    panel.grid.minor = element_line(color = "grey90")   # 添加次要网格线
  )

ggsave(filename = 'bar_plot_up_down(day_night).png', plot = bar_plot)


#LL_SL_DESeq <- read.table("LL_SL.DESeq2.txt", header = T, sep = '\t', quote = "")
#LD_SD_DESeq <- read.table("LD_SD.DESeq2.txt", header = T, sep = '\t', quote = "")

BHLH062_ids <- c("AeUmb.TA1851.r1.2UG0002210", "AeUmb.TA1851.r1.2UG0045890")

#LL_SL_filtered <- LL_SL_DESeq[LL_SL_DESeq$ID %in% BHLH062_ids, ]
#LD_SD_filtered <- LD_SD_DESeq[LD_SD_DESeq$ID %in% BHLH062_ids, ]


BHLH062_ids
expression_levels_47 <- normalized_counts[rownames(normalized_counts) %in% matched_genes_47$ID, ]

BHLH062_expression_levels <- normalized_counts[rownames(normalized_counts) %in% BHLH062_ids, ]

expression_levels_47_BHLH062 <- rbind(expression_levels_47, BHLH062_expression_levels)
write.csv(expression_levels_47_BHLH062, "expression_levels_47_BHLH062.csv", row.names = TRUE)



expression_levels_109 <- normalized_counts[rownames(normalized_counts) %in% matched_gene_109$ID, ]
write.csv(expression_levels_109, "expression_levels_109.csv", row.names = TRUE)


expression_levels_47_BHLH062_with_id <- as.data.frame(expression_levels_47_BHLH062)
# 将 expression_levels_47_BHLH062 的行名保存为一个新列
expression_levels_47_BHLH062_with_id <- expression_levels_47_BHLH062_with_id %>%
  mutate(gene_id = rownames(expression_levels_47_BHLH062_with_id))

# 计算长光周期（LD, LL）和短光周期（SD, SL）的表达变化
data_opposite <- expression_levels_47_BHLH062_with_id %>%
  rowwise() %>%
  mutate(
    # 计算长光周期变化：从LD到LL
    change_LD_LL = mean(c(LD2, LD1)) - mean(c(LL1, LL2)),
    # 计算短光周期变化：从SD到SL
    change_SD_SL = mean(c(SD2, SD1)) - mean(c(SL1, SL2))
  ) %>%
  ungroup()  # 解除 rowwise 模式

# 将行名保存到一个新的列中
data_opposite <- data_opposite %>%
  mutate(gene_id = rownames(expression_levels_47_BHLH062))

# 筛选符合条件的基因：长光周期下表达上升，短光周期下表达下降
data_opposite_expression <- data_opposite %>%
  filter(
    (change_LD_LL > 0 & change_SD_SL < 0) | (change_LD_LL < 0 & change_SD_SL > 0)
  )

data_opposite_expression <- as.data.frame(data_opposite_expression)

# 将 gene_id 列重新设置为行名
rownames(data_opposite_expression) <- data_opposite_expression$gene_id

# 移除 gene_id 列
data_opposite_expression <- data_opposite_expression %>% select(-gene_id)

write.csv(data_opposite_expression, "data_opposite_expression.csv", row.names = TRUE)




expression_levels_596 <- normalized_counts[rownames(normalized_counts) %in% matched_gene_596$ID, ]
write.csv(expression_levels_596, "expression_levels_596.csv", row.names = TRUE)

expression_levels_219 <- normalized_counts[rownames(normalized_counts) %in% matched_gene_219$ID, ]
write.csv(expression_levels_219, "expression_levels_219.csv", row.names = TRUE)

expression_levels_381 <- normalized_counts[rownames(normalized_counts) %in% matched_gene_381$ID, ]
write.csv(expression_levels_381, "expression_levels_381.csv", row.names = TRUE)

expression_levels_554 <- normalized_counts[rownames(normalized_counts) %in% matched_gene_554$ID, ]
write.csv(expression_levels_554, "expression_levels_554.csv", row.names = TRUE)






library(clusterProfiler)
library(dplyr)
library(tidyr)

daytime_up_219_47 <- c(daytime_up_only_219, day_night_47)
nighttime_down_554_47 <- c(nighttime_down_only_554, day_night_47)

library(org.At.tair.db)
daytime_up <- enrichGO(gene = daytime_up_219_47,
                  OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                  keyType = "SYMBOL",  # 根据基因名称进行富集分析
                  ont = "BP",  # 生物过程
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "none",
                  qvalueCutoff = 1)
dim(daytime_up)
daytime_up_219_47_go <- daytime_up@result
write.csv(daytime_up_219_47_go, "daytime_up_219_47_go.csv", row.names = FALSE)


nighttime_down <- enrichGO(gene = nighttime_down_554_47,
                    OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                    keyType = "SYMBOL",  # 根据基因名称进行富集分析
                    ont = "BP",  # 生物过程
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "none",
                    qvalueCutoff = 1)

dim(nighttime_down)
nighttime_down_554_47_go <- nighttime_down@result
write.csv(nighttime_down_554_47_go, "nighttime_down_554_47_go.csv", row.names = FALSE)


#结果可视化 
daytime_up_219_47_go_bar <- barplot(daytime_up,showCategory=15,drop=T) 
daytime_up_219_47_go_bar <- daytime_up_219_47_go_bar + theme(
  axis.text.y = element_text(size = 8)) +  # 调整 y 轴标签的字体大小
  ggtitle("Daytime_up_Bar")

daytime_up_219_47_go_bar
ggsave(filename = "daytime_up_219_47_go_bar.png", plot = daytime_up_219_47_go_bar)

daytime_up_219_47_go_dot <- dotplot(daytime_up,showCategory=15)
daytime_up_219_47_go_dot <- daytime_up_219_47_go_dot + theme(
  axis.text.y = element_text(size = 8)) +  # 调整 y 轴标签的字体大小
  ggtitle("Daytime_up_Dot")

daytime_up_219_47_go_dot
ggsave(filename = "daytime_up_219_47_go_dot.png", plot = daytime_up_219_47_go_dot)


nighttime_down_554_47_bar <- barplot(nighttime_down,showCategory=15,drop=T) 
nighttime_down_554_47_bar <- nighttime_down_554_47_bar + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("Nighttime_down_554_47_bar")

nighttime_down_554_47_bar
ggsave(filename = "nighttime_down_554_47_bar.png", plot = nighttime_down_554_47_bar)

nighttime_down_554_47_dot <- dotplot(nighttime_down,showCategory=15)
nighttime_down_554_47_dot <- nighttime_down_554_47_dot + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("Nighttime_down_554_47_dot")

nighttime_down_554_47_dot
ggsave(filename = "Nighttime_down_554_47_dot.png", plot = nighttime_down_554_47_dot)







nighttime_up_daytime_down_596_109 <- c(nighttime_up_daytime_down_109,
                                       nighttime_up_only_596)
daytime_down_381_109 <- c(nighttime_up_daytime_down_109, daytime_down_only_381)

library(org.At.tair.db)
nighttime_up_daytime_down <- enrichGO(gene = nighttime_up_daytime_down_596_109,
                       OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                       keyType = "SYMBOL",  # 根据基因名称进行富集分析
                       ont = "BP",  # 生物过程
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "none",
                       qvalueCutoff = 1)
dim(nighttime_up_daytime_down)
nighttime_up_daytime_down_go <- nighttime_up_daytime_down@result
write.csv(nighttime_up_daytime_down_go, "nighttime_up_daytime_down_go.csv",
          row.names = FALSE)


daytime_down <- enrichGO(gene = daytime_down_381_109,
                           OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                           keyType = "SYMBOL",  # 根据基因名称进行富集分析
                           ont = "BP",  # 生物过程
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "none",
                           qvalueCutoff = 1)
dim(daytime_down)
daytime_down_go <- daytime_down@result
write.csv(daytime_down_go, "daytime_down_go.csv",
          row.names = FALSE)

#结果可视化 
nighttime_up_daytime_down_go_bar <- barplot(nighttime_up_daytime_down,showCategory=15,drop=T) 
nighttime_up_daytime_down_go_bar <- nighttime_up_daytime_down_go_bar + theme(
  axis.text.y = element_text(size = 8)) +  # 调整 y 轴标签的字体大小
  ggtitle("Nighttime_Up_Daytime_Down_Bar")

nighttime_up_daytime_down_go_bar
ggsave(filename = "nighttime_up_daytime_down_go_bar.png", plot = nighttime_up_daytime_down_go_bar)

nighttime_up_daytime_down_go_dot <- dotplot(nighttime_up_daytime_down,showCategory=15)
nighttime_up_daytime_down_go_dot <- nighttime_up_daytime_down_go_dot + theme(
  axis.text.y = element_text(size = 8)) +  # 调整 y 轴标签的字体大小
  ggtitle("Nighttime_Up_Daytime_Down_Dot")

nighttime_up_daytime_down_go_dot
ggsave(filename = "nighttime_up_daytime_down_go_dot.png", plot = nighttime_up_daytime_down_go_dot)


daytime_down_bar <- barplot(daytime_down,showCategory=15,drop=T) 
daytime_down_bar <- daytime_down_bar + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("Daytime_Down_bar")

daytime_down_bar
ggsave(filename = "daytime_down_bar.png", plot = daytime_down_bar)

daytime_down_bar_dot <- dotplot(daytime_down,showCategory=15)
daytime_down_bar_dot <- daytime_down_bar_dot + theme(
  axis.text.y = element_text(size = 10)) +  # 调整 y 轴标签的字体大小
  ggtitle("Daytime_Down_dot")

daytime_down_bar_dot
ggsave(filename = "daytime_down_bar_dot.png", plot = daytime_down_bar_dot)



compare_pathway_genes <- function(df1, df2, row_num) {
  # 提取基因
  genes_df1 <- strsplit(df1$geneID[row_num], "/")[[1]]
  genes_df2 <- strsplit(df2$geneID[row_num], "/")[[1]]
  
  # 计算差异
  list(
    Term = df1$Description[row_num],
    Unique_in_df1 = setdiff(genes_df1, genes_df2),
    Unique_in_df2 = setdiff(genes_df2, genes_df1),
    Common = intersect(genes_df1, genes_df2)
  )
}
# 对比第一行
result_row1 <- compare_pathway_genes(LD_SD_go_up, WT_LL_SL_go_up, 1)
# 转换为长格式 data.frame
result_df <- data.frame(
  Gene = unlist(result_row1[c("Unique_in_df1", "Unique_in_df2", "Common")]),
  Category = rep(c("Unique_in_df1", "Unique_in_df2", "Common"),
                 times = sapply(result_row1[c("Unique_in_df1", "Unique_in_df2", "Common")], length))
)

# 可选：添加描述（Term）
result_df$Term <- result_row1$Term
result_df$Category <- recode(result_df$Category,
                             Unique_in_df1 = "LD_SD_only",
                             Unique_in_df2 = "WT_LL_SL_only",
                             Common = "Both")
write.csv(result_df, "ribosome_biogenesis_genes.csv", row.names = FALSE)

# 对比第二行
result_row2 <- compare_pathway_genes(LD_SD_go_up, WT_LL_SL_go_up, 2)
result_df2 <- data.frame(
  Gene = unlist(result_row2[c("Unique_in_df1", "Unique_in_df2", "Common")]),
  Category = rep(c("Unique_in_df1", "Unique_in_df2", "Common"),
                 times = sapply(result_row1[c("Unique_in_df1", "Unique_in_df2", "Common")], length))
)

# 可选：添加描述（Term）
result_df2$Term <- result_row2$Term
result_df2$Category <- recode(result_df$Category,
                             Unique_in_df1 = "LD_SD_only",
                             Unique_in_df2 = "WT_LL_SL_only",
                             Common = "Both")
write.csv(result_df2, "ribonucleoprotein_complex_biogenesis_genes.csv", row.names = FALSE)


# log2foldchange SD vs XD
# log2foldchange LD vs XD

library(tibble)
matched_47 <- read.csv("matched_genes_day_UP_night_DOWN_47.csv")

XD_SD <- read.table("XD_SD.DESeq2.txt", header = TRUE, sep = "\t", quote = "")
XD_LD <- read.table("XD_LD.DESeq2.txt", header = TRUE, sep = "\t", quote = "")

# 找到 log2FoldChange 列的位置
log2fc_pos <- which(colnames(matched_47) == "log2FoldChange")

# 在 log2FoldChange 后插入两列
matched_47 <- matched_47 %>%
  add_column(
    XD_SD_log2fc = XD_SD$log2FoldChange[match(matched_47$ID, XD_SD$ID)],
    XD_LD_log2fc = XD_LD$log2FoldChange[match(matched_47$ID, XD_LD$ID)],
    .after = log2fc_pos
  )
write.csv(matched_47, "LD_SD_XD.csv")


library(GENIE3)
LD_SD_up <- read.table("LD_SD.DESeq2.up.txt", header = TRUE, sep = "\t", quote = "")
LD_SD_down <- read.table("LD_SD.DESeq2.down.txt", header = TRUE, sep = "\t", quote = "")
LD_SD <- rbind(LD_SD_up,LD_SD_down)

XD_SD_up <-read.table("XD_SD.DESeq2.up.txt", header = TRUE, sep = "\t", quote = "")
XD_SD_down <-read.table("XD_SD.DESeq2.down.txt", header = TRUE, sep = "\t", quote = "")
XD_SD <- rbind(XD_SD_down, XD_SD_up)

XD_LD_up <- read.table("XD_LD.DESeq2.up.txt", header = TRUE, sep = "\t", quote = "")
XD_LD_down <-read.table("XD_LD.DESeq2.down.txt", header = TRUE, sep = "\t", quote = "")
XD_LD <- rbind(XD_LD_down, XD_LD_up)

LD_SD_XD_condition <- rbind(LD_SD,XD_SD,XD_LD)


LD_SD_XD_condition_gene <- unique(LD_SD_XD_condition$ID)
write.table(LD_SD_XD_condition_gene, "LD_SD_XD_condition_gene.txt", row.names = FALSE, quote = F)


LD_SD_XD_expression_data <- normalized_counts[rownames(normalized_counts) %in% LD_SD_XD_condition_gene, ]

LD_SD_XD_TF_list <- read.table("LD_SD_XD_TF.list.txt", header = FALSE, stringsAsFactors = FALSE)


LD_SD_XD_TF_list$V1 <- sub("\\.\\d+$", "", LD_SD_XD_TF_list$V1)

LD_SD_XD_TF_list <- LD_SD_XD_TF_list[!duplicated(LD_SD_XD_TF_list$V1), ]

LD_SD_XD_TF_ids <- as.character(LD_SD_XD_TF_list$V1)

# 构建基因调控网络
# 获取有效调控因子
valid_regulators <- intersect(LD_SD_XD_TF_ids, rownames(LD_SD_XD_expression_data))
set.seed(6)
weight_matrix <- GENIE3(as.matrix(LD_SD_XD_expression_data), regulators = valid_regulators)

# 提取边列表

LD_SD_XD_link_list <- getLinkList(weight_matrix, threshold = 0.001)

#gene_annotations$gene_name[is.na(gene_annotations$gene_name)] <- "unknown"
na_indices <- which(is.na(gene_annotations$gene_name))
unique_na <- paste0("unknown_", seq_along(na_indices))
gene_annotations$gene_name[na_indices] <- unique_na


# 添加原始的 seq ID 列
LD_SD_XD_link_list$regulatorySeqID <- LD_SD_XD_link_list$regulatoryGene
LD_SD_XD_link_list$targetSeqID <- LD_SD_XD_link_list$targetGene

# 用 gene ID 替换 seq ID，同时保留原始 seq ID
LD_SD_XD_link_list$regulatoryGene <- gene_annotations$gene_name[match(LD_SD_XD_link_list$regulatorySeqID, gene_annotations$ID)]
LD_SD_XD_link_list$targetGene <- gene_annotations$gene_name[match(LD_SD_XD_link_list$targetSeqID, gene_annotations$ID)]
# 保存结果

write.table(link_list, "LD_SD_XD_regulatory_network_with_seqID.csv", sep = "\t", row.names = FALSE, quote = FALSE)










