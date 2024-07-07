library(DESeq2)
setwd("/Users/chenziyin/Downloads/Dissertation")

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
  res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.01 & res1$baseMean >= 100), 'sig'] <- 'up'
  res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.01 & res1$baseMean >= 100), 'sig'] <- 'down'
  res1[which(abs(res1$log2FoldChange) < 1 | res1$padj >= 0.01 | res1$baseMean < 100), 'sig'] <- 'none'
  
  
  res1$ensembl_gene_id <- rownames(res1)
  merged_results <- merge(gene_annotations, res1, by.x = "ID",
                          by.y = "ensembl_gene_id", all.x = TRUE, all.y = T)
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

  # 筛选前10个差异显著的基因
  top5 <- merged_results[order(merged_results$padj, -abs(merged_results$log2FoldChange)), ][1:5, ]
  top5_up <- res1_up[order(res1_up$padj, -abs(res1_up$log2FoldChange)), ][1:5, ]
  top5_down <- res1_down[order(res1_down$padj, -abs(res1_down$log2FoldChange)), ][1:5, ]
  highlight_genes <- rbind(top5, top5_up, top5_down)
  #write.table(highlight_genes,file = paste0(comparison,'.highlight.txt'),sep = '\t', col.name = NA, quote = F)
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
           fontsize_col = 5,
           angle_col = 45)
  
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

LD_LL_genes <- unique(c(LD_LL_up$ID, LD_LL_down$ID))
LD_SD_genes <- unique(c(LD_SD_up$ID, LD_SD_down$ID))
LL_SL_genes <- unique(c(LL_SL_up$ID, LL_SL_down$ID))
SD_SL_genes <- unique(c(SD_SL_up$ID, SD_SL_down$ID))

# 保存为 CSV 文件
write.csv(LD_LL_genes, "LD_LL_genes.csv", row.names = FALSE)
write.csv(LD_SD_genes, "LD_SD_genes.csv", row.names = FALSE)
write.csv(LL_SL_genes, "LL_SL_genes.csv", row.names = FALSE)
write.csv(SD_SL_genes, "SD_SL_genes.csv", row.names = FALSE)













