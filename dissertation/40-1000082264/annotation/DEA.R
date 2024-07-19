library(DESeq2)
#setwd("/Users/chenziyin/Downloads/Dissertation")
setwd('/localdisk/home/s2558632/dissertation/40-1000082264/annotation')
# 样本名称
sampleNames <- c("LD1", "LD2", "LL1", "LL2", "SD1", "SD2", "SL1", "SL2")

for (name in sampleNames) {
  file_path <- paste0(name, "_counts.txt")
  data <- read.table(file_path, header = TRUE, row.names = 1)
  assign(paste0(name, "_counts"), data)
}
DataMatrix <- cbind(LD1_counts[6],LD2_counts[6],LL1_counts[6],LL2_counts[6],
                    SD1_counts[6],SD2_counts[6],SL1_counts[6],SL2_counts[6])
#LD1 & LD2 <-> LL1 & LL2 
#LD1 & LD2 <-> SD1 & SD2
#LL1 & LL2 <-> SL1 & SL2
#SD1 & SD2 <-> SL1 & SL2
colnames(DataMatrix) <- sampleNames

# 定义组信息
comparisons <- list(
  LD_LL = list(samples = c("LD1", "LD2", "LL1", "LL2"), group = c("LD", "LD", "LL", "LL")),
  LD_SD = list(samples = c("LD1", "LD2", "SD1", "SD2"), group = c("LD", "LD", "SD", "SD")),
  LL_SL = list(samples = c("LL1", "LL2", "SL1", "SL2"), group = c("LL", "LL", "SL", "SL")),
  SD_SL = list(samples = c("SD1", "SD2", "SL1", "SL2"), group = c("SD", "SD", "SL", "SL"))
)


library(rtracklayer)

# 读取GFF3文件
gff3_file <- "../alignment/AeUmbellulata_TA1851_v1.gff3"
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
    scale_color_manual(values = c("#ff4757", "#d2dae2", "#546de5"), limits = c('up', 'none', 'down')) +
    labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = paste(comparison, 'Comparison'), color = '') +
    theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          legend.key = element_rect(fill = 'transparent')) +
    geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +
    geom_hline(yintercept = 2, lty = 3, color = 'black') +
    xlim(-10, 10) + ylim(0, 250)
  p <- p + geom_text(data = highlight_genes, aes(label = gene_name),
                     vjust = 2, hjust = 0.5, check_overlap = TRUE, size = 3)
  
  #p <- p + geom_text(data = target_genes, aes(label = ensembl_gene_id),
                     #vjust = 2, hjust = 0.5, check_overlap = TRUE, size = 3)
  
  # 显式打印图像并保存
  print(p)
  ggsave(filename = paste0(comparison, '_volcano_plot.png'), plot = p)
  
  datamatrix_up <- DataMatrix_sub[rownames(DataMatrix_sub) %in% res1_up$ID, ]
  datamatrix_down <- DataMatrix_sub[rownames(DataMatrix_sub) %in% res1_down$ID, ]
  datamatrix <- rbind(datamatrix_up, datamatrix_down)
  Names_up <- res1_up$gene_name[match(rownames(datamatrix_up), res1_up$ID)]
  Names_down <- res1_down$gene_name[match(rownames(datamatrix_down), res1_down$ID)]
  Names <- c(Names_up, Names_down)
  Names[is.na(Names)] <- "unknown"
  
  unique_groups <- as.data.frame(unique(Names))
  colnames(unique_groups) <- 'Groups'
  Names <- make.unique(Names)
  rownames(datamatrix) <- Names
  
  group_annotation <- sapply(rownames(datamatrix), function(x) {
    # 检查是否包含唯一名称
    match_found <- sapply(unique_groups$Groups, function(y) grepl(paste0("\\b",y,"\\b"),x))
    if (any(match_found)) {
      return(unique_groups$Groups[which(match_found)[1]])
    } else if (grepl("unknown", x)) {
      return("Unknown")
    } else {
      return(x)
    }
  })
  
  row_annotation <- data.frame(Group = group_annotation)
  
  
  heatmap <- pheatmap(datamatrix,
           cluster_rows = TRUE, # Clustering of rows
           cluster_cols = TRUE, # Clustering of cols
           scale="row", # Standardisation with behavioural benchmarks
           main = paste(comparison, 'Comparison'),
           show_rownames = F,
           clustering_method = "ward.D2", # Using the square of the Euclidean distance is less subject to outliers.
           border = F,
           annotation_row = row_annotation, # 添加行注释
           color = colorRampPalette(c("#D72E25", "#F46D43", "#FEE08B","#FFFFBF","#A6D96A","#1A9850"))(50),
           cutree_rows = 4,
           fontsize_col = 10,
           angle_col = 45,
           annotation_legend = FALSE)
  
  # 保存热图
  ggsave(filename = paste0(comparison, '_heatmap.png'), plot = heatmap)
}
library(VennDiagram)
# 读取保存的文件
LD_LL_up <- read.table("LD_LL.DESeq2.up.txt", sep = '\t', header = TRUE)
LD_LL_down <- read.table("LD_LL.DESeq2.down.txt", sep = '\t', header = TRUE)
LD_SD_up <- read.table("LD_SD.DESeq2.up.txt", sep = '\t', header = TRUE)
LD_SD_down <- read.table("LD_SD.DESeq2.down.txt", sep = '\t', header = TRUE)
LL_SL_up <- read.table("LL_SL.DESeq2.up.txt", sep = '\t', header = TRUE)
LL_SL_down <- read.table("LL_SL.DESeq2.down.txt", sep = '\t', header = TRUE)
SD_SL_up <- read.table("SD_SL.DESeq2.up.txt", sep = '\t', header = TRUE)
SD_SL_down <- read.table("SD_SL.DESeq2.down.txt", sep = '\t', header = TRUE)

LD_LL_up <- LD_LL_up[!is.na(LD_LL_up$gene_name), ]
LD_LL_down <- LD_LL_down[!is.na(LD_LL_down$gene_name), ]
LD_SD_up <- LD_SD_up[!is.na(LD_SD_up$gene_name), ]
LD_SD_down <- LD_SD_down[!is.na(LD_SD_down$gene_name), ]
LL_SL_up <- LL_SL_up[!is.na(LL_SL_up$gene_name), ]
LL_SL_down <- LL_SL_down[!is.na(LL_SL_down$gene_name), ]
SD_SL_up <- SD_SL_up[!is.na(SD_SL_up$gene_name), ]
SD_SL_down <- SD_SL_down[!is.na(SD_SL_down$gene_name), ]

LD_LL_genes <- unique(c(LD_LL_up$gene_name, LD_LL_down$gene_name))
LD_SD_genes <- unique(c(LD_SD_up$gene_name, LD_SD_down$gene_name))
LL_SL_genes <- unique(c(LL_SL_up$gene_name, LL_SL_down$gene_name))
SD_SL_genes <- unique(c(SD_SL_up$gene_name, SD_SL_down$gene_name))

gene_list <- list(
  "LD_LL" = LD_LL_genes,
  "LD_SD" = LD_SD_genes,
  "LL_SL" = LL_SL_genes,
  "SD_SL" = SD_SL_genes
)

venn.diagram(
  x = gene_list,
  filename = "venn_diagram.png",
  fill = c("#D72E25", "#F46D43", "#1A9850", "#A6D96A"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("#D72E25", "#F46D43", "#1A9850", "#A6D96A"),
  main = "Venn Diagram of Differentially Expressed Genes"
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


# 提取基因ID（确保是字符向量）
up_genes <- as.character(LL_SL_up$gene_name)  # 使用基因名称而不是基因ID
down_genes <- as.character(LL_SL_down$gene_name)  # 使用基因名称而不是基因ID
up_genes <- up_genes[!is.na(up_genes)]
down_genes <- down_genes[!is.na(down_genes)]

#up_genes <- as.character(LD_SD_up$gene_name)  # 使用基因名称而不是基因ID
#down_genes <- as.character(LD_SD_down$gene_name)  # 使用基因名称而不是基因ID
#up_genes <- up_genes[!is.na(up_genes)]
#down_genes <- down_genes[!is.na(down_genes)]

total_gene <- c(up_genes,down_genes)
library(org.At.tair.db)
go_up <- enrichGO(gene = up_genes,
           OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
           keyType = "SYMBOL",  # 根据基因名称进行富集分析
           ont = "BP",  # 生物过程
           pvalueCutoff = 0.4,
           pAdjustMethod = "BH",
           qvalueCutoff = 0.4)
dim(go_up)

go_down <- enrichGO(gene = down_genes,
                  OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                  keyType = "SYMBOL",  # 根据基因名称进行富集分析
                  ont = "BP",  # 生物过程
                  pvalueCutoff = 0.4,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.4)
dim(go_down)


#结果可视化 
LL_SL_up_bar <- barplot(go_up,showCategory=10,drop=T) 
LL_SL_up_bar <- LL_SL_up_bar + theme(
  axis.text.y = element_text(size = 8)  # 调整 y 轴标签的字体大小
)
LL_SL_up_bar
ggsave(filename = "LL_SL_GO_Enrichment_Upregulated_Bar_Adjusted_Font.png", plot = LL_SL_up_bar)
LL_SL_up_dot <- dotplot(go_up,showCategory=10)
LL_SL_up_dot <- LL_SL_up_dot + theme(
  axis.text.y = element_text(size = 8)  # 调整 y 轴标签的字体大小
)
LL_SL_up_dot
ggsave(filename = "LL_SL_GO_Enrichment_Upregulated_Dot_Adjusted_Font.png", plot = LL_SL_up_dot)


LL_SL_down_bar <- barplot(go_down,showCategory=15,drop=T) 
LL_SL_down_bar <- LL_SL_down_bar + theme(
  axis.text.y = element_text(size = 8)  # 调整 y 轴标签的字体大小
)
LL_SL_down_bar
ggsave(filename = "LL_SL_GO_Enrichment_Downregulated_Bar_Adjusted_Font.png", plot = LL_SL_down_bar)
LL_SL_down_dot <- dotplot(go_down,showCategory=15)
LL_SL_down_dot <- LL_SL_down_dot + theme(
  axis.text.y = element_text(size = 8)  # 调整 y 轴标签的字体大小
)
LL_SL_down_dot
ggsave(filename = "LL_SL_GO_Enrichment_Downregulated_Dot_Adjusted_Font.png", plot = LL_SL_down_dot)





# 提取基因ID（确保是字符向量）
up_genes <- as.character(LD_SD_up$gene_name)  # 使用基因名称而不是基因ID
down_genes <- as.character(LD_SD_down$gene_name)  # 使用基因名称而不是基因ID
up_genes <- up_genes[!is.na(up_genes)]
down_genes <- down_genes[!is.na(down_genes)]
library(org.At.tair.db)
go_up <- enrichGO(gene = up_genes,
                  OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                  keyType = "SYMBOL",  # 根据基因名称进行富集分析
                  ont = "BP",  # 生物过程
                  pvalueCutoff = 0.4,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.4)
dim(go_up)

go_down <- enrichGO(gene = down_genes,
                    OrgDb = org.At.tair.db,  # 使用拟南芥的注释包
                    keyType = "SYMBOL",  # 根据基因名称进行富集分析
                    ont = "BP",  # 生物过程
                    pvalueCutoff = 0.4,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.4)
dim(go_down)


#结果可视化 
LD_SD_up_bar <- barplot(go_up,showCategory=10,drop=T) 
LD_SD_up_bar <- LL_SL_up_bar + theme(
  axis.text.y = element_text(size = 8)  # 调整 y 轴标签的字体大小
)
LD_SD_up_bar
ggsave(filename = "LD_SD_GO_Enrichment_Upregulated_Bar_Adjusted_Font.png", plot = LD_SD_up_bar)
LD_SD_up_dot <- dotplot(go_up,showCategory=10)
LD_SD_up_dot <- LD_SD_up_dot + theme(
  axis.text.y = element_text(size = 8)  # 调整 y 轴标签的字体大小
)
LD_SD_up_dot
ggsave(filename = "LD_SD_GO_Enrichment_Upregulated_Dot_Adjusted_Font.png", plot = LD_SD_up_dot)


LD_SD_down_bar <- barplot(go_down,showCategory=15,drop=T) 
LD_SD_down_bar <- LD_SD_down_bar + theme(
  axis.text.y = element_text(size = 8)  # 调整 y 轴标签的字体大小
)
LD_SD_down_bar
ggsave(filename = "LD_SD_GO_Enrichment_Downregulated_Bar_Adjusted_Font.png", plot = LD_SD_down_bar)
LD_SD_down_dot <- dotplot(go_down,showCategory=15)
LD_SD_down_dot <- LD_SD_down_dot + theme(
  axis.text.y = element_text(size = 8)  # 调整 y 轴标签的字体大小
)
LD_SD_down_dot
ggsave(filename = "LD_SD_GO_Enrichment_Downregulated_Dot_Adjusted_Font.png", plot = LD_SD_down_dot)


BiocManager::install("GENIE3")
library(GENIE3)
SD_SL_highlight_75 <- read.table("SD_SL.highlight.txt", sep = '\t', header = TRUE)
LD_LL_highlight_75 <- read.table("LD_LL.highlight.txt", sep = '\t', header = TRUE)

SD_SL_genes <- unique(SD_SL_highlight_75$ID)
LD_LL_genes <- unique(LD_LL_highlight_75$ID)

combined_genes <- c(SD_SL_genes, LD_LL_genes, LL_SL_genes, LD_SD_genes)
combined_genes <- unique(combined_genes)

#LL_SL_expression_data <- DataMatrix[rownames(DataMatrix) %in% LL_SL_genes, ]
#LD_SD_expression_data <- DataMatrix[rownames(DataMatrix) %in% LD_SD_genes, ]
expression_data <- DataMatrix[rownames(DataMatrix) %in% combined_genes, ]

LL_SL_TF_list <- read.table("LL_SL_TF.list.txt", header = FALSE, stringsAsFactors = FALSE)
LD_SD_TF_list <- read.table("LD_SD_TF.list.txt", header = FALSE, stringsAsFactors = FALSE)
SD_SL_TF_list <- read.table("SD_SL_TF.list.txt", header = FALSE, stringsAsFactors = FALSE)
LD_LL_TF_list <- read.table("LD_LL_TF.list.txt", header = FALSE, stringsAsFactors = FALSE)
LL_SL_TF_list$V1 <- sub("\\.\\d+$", "", LL_SL_TF_list$V1)
LD_SD_TF_list$V1 <- sub("\\.\\d+$", "", LD_SD_TF_list$V1)
SD_SL_TF_list$V1 <- sub("\\.\\d+$", "", SD_SL_TF_list$V1)
LD_LL_TF_list$V1 <- sub("\\.\\d+$", "", LD_LL_TF_list$V1)

#LL_SL_TF_ids <- as.character(LL_SL_TF_list$V1)
#LD_SD_TF_ids <- as.character(LD_SD_TF_list$V1)
#SD_SL_TF_ids <- as.character(SD_SL_TF_list$V1)
#LD_LL_TF_ids <- as.character(LD_LL_TF_list$V1)
#TF_ids <- c(LL_SL_TF_ids,LD_SD_TF_ids,SD_SL_TF_ids,LD_LL_TF_ids)

TF_list <- rbind(LL_SL_TF_list, LD_SD_TF_list, SD_SL_TF_list, LD_LL_TF_list)
TF_list <- TF_list[!duplicated(TF_list$V1), ]
TF_ids <- as.character(TF_list$V1)

# 构建基因调控网络
#LL_SL_weight_matrix <- GENIE3(as.matrix(LL_SL_expression_data), regulators = LL_SL_TF_ids)
#LD_SD_weight_matrix <- GENIE3(as.matrix(LD_SD_expression_data), regulators = LD_SD_TF_ids)

weight_matrix <- GENIE3(as.matrix(expression_data), regulators = TF_ids)

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
link_list$regulatoryGene <- gene_annotations$gene_name[match(link_list$regulatoryGene, gene_annotations$ID)]
link_list$targetGene <- gene_annotations$gene_name[match(link_list$targetGene, gene_annotations$ID)]

# 保存结果
#write.csv(LL_SL_link_list, "LL_SL_gene_regulatory_network.csv", row.names = FALSE)
#write.csv(LD_SD_link_list, "LD_SD_gene_regulatory_network.csv", row.names = FALSE)
write.csv(link_list, "gene_regulatory_network.csv", row.names = FALSE)



