################################################################
### Characteristics of malignant breast cancer cell clusters ###
################################################################

#################
### Load data ###
#################
cancer.so <- readRDS("cancer_SO.rds")

###############
### Gene 70 ###
###############
### We used AddModuleScore function in Seurat
### Load Gene70
data(sig.gene70)
gene70 <- sig.gene70$NCBI.gene.symbol
gene70 <- na.omit(gene70)

cancer.so <- AddModuleScore(cancer.so, features = list(gene70), name = 'Gene70', assay = 'SCT')

#######################################################
### Identify DEGs in cluster 8 and 19 vs. cluster 1 ###
#######################################################
cancer.so <- SetIdent(cancer.so, value = 'seurat_clusters')
DEG_1_8 <- FindMarkers(cancer.so, ident.1 = c(8), ident.2 = c(1), logfc.threshold = 0, min.diff.pct = 0, min.pct = 0)
DEG_1_19 <- FindMarkers(cancer.so, ident.1 = c(19), ident.2 = c(1), logfc.threshold = 0, min.diff.pct = 0, min.pct = 0)

DEG_1_8$p_val <- ifelse(DEG_1_8$p_val == 0, min(DEG_1_8$p_val[DEG_1_8$p_val > 0]), DEG_1_8$p_val)
DEG_1_19$p_val <- ifelse(DEG_1_19$p_val == 0, min(DEG_1_19$p_val[DEG_1_19$p_val > 0]), DEG_1_19$p_val)

DEG_1_8$gene <- rownames(DEG_1_8)
DEG_1_19$gene <- rownames(DEG_1_19)

### Save files
write.table(DEG_1_8, 'DEG_1_8.txt')
write.table(DEG_1_19, 'DEG_1_19.txt')

#####################
### Visualization ###
#####################
png(
  filename  = 'paper_figure/Fig_3a.png',
  width     = 13,
  height    = 4,
  unit = 'in',
  res = 300
)
VlnPlot(cancer.so, features = 'Gene701', pt.size = 0, group.by = 'seurat_clusters', split.by = 'pam50', cols = c('#EB8677','#7DC0A6','#919FC7','#Da8EC0','#8CBA54')) + 
  ggtitle("") + ylab("Module score of Gene70") + xlab("Clusters") + labs(fill = "PAM50\nclassification") + 
  scale_y_continuous(limits = c(-0.12, 0.4), breaks = seq(-0.1, 0.3, by = 0.1), expand = expansion(mult = c(0, 0))) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15), axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18)) + 
  stat_compare_means(comparisons = list(c('1','8'), c('1','19')), label = "p.val", size = 5) 
dev.off()

### Fig 3b-1
### Cluster 8 vs. cluster 1
# add a column of NAs
DEG_1_8$DE <- "No-sig."
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
DEG_1_8$DE[DEG_1_8$avg_log2FC > 0.451 & DEG_1_8$p_val < 0.05] <- "Up-reg."
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
DEG_1_8$DE[DEG_1_8$avg_log2FC < -0.451 & DEG_1_8$p_val < 0.05] <- "Down-reg."
label <- c('EEF1A1','RPS4X','RPL19','TUBA1B','HIST1H4C','STMN1')### select top3 genes of each side
DEG_1_8$label <- ifelse(DEG_1_8$gene %in% label, DEG_1_8$gene, NA)
pos <- position_jitter(width = 0.3, seed = 2)

png(
  filename  = 'paper_figure/Fig_3b-1.png',
  width     = 6,
  height    = 4,
  unit = 'in',
  res = 300
)

ggplot(DEG_1_8, aes(x = avg_log2FC, y = -log10(p_val), col = DE)) + geom_point() + 
  geom_vline(xintercept = c(-0.451, 0.451), linetype = 'dashed') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  scale_color_manual(values = c('#919FC7','grey80','#EB8677')) + 
  theme_classic() + xlab(expression('log'[2]*'(Fold change)')) + ylab(expression('-log'[10]*'(p-value)')) + 
  scale_y_continuous(expand = expansion(add = c(0, 10))) +
  ggtitle("Cluster 8 vs. cluster 1") +
  geom_label_repel(aes(label = label),
                   point.padding = 0.5,
                   segment.color = 'grey80', 
                   nudge_x = 0.5,
                   show_guide  = FALSE) +
  theme(axis.text = element_text(size = 18, colour = 'black'), axis.title = element_text(size = 20), 
        title = element_text(size = 18), legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))

dev.off()

### Fig 3b-2 
### Cluster 19 vs. cluster 1
# add a column of NAs
DEG_1_19$DE <- "No-sig."
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
DEG_1_19$DE[DEG_1_19$avg_log2FC > 0.451 & DEG_1_19$p_val < 0.05] <- "Up-reg."
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
DEG_1_19$DE[DEG_1_19$avg_log2FC < -0.451 & DEG_1_19$p_val < 0.05] <- "Down-reg."
label <- c('PABPC1','RPL30','VIM','TFF3','PIP','TFF1')### select top3 genes of each side
DEG_1_19$label <- ifelse(DEG_1_19$gene %in% label, DEG_1_19$gene, NA)
pos <- position_jitter(width = 0.3, seed = 2)

png(
  filename  = 'paper_figure/Fig_3b-2.png',
  width     = 6,
  height    = 4,
  unit = 'in',
  res = 300
)

ggplot(DEG_1_19, aes(x = avg_log2FC, y = -log10(p_val), col = DE)) + geom_point() + 
  geom_vline(xintercept = c(-0.451, 0.451), linetype = 'dashed') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  scale_color_manual(values = c('#919FC7','grey80','#EB8677')) + 
  theme_classic() + xlab(expression('log'[2]*'(Fold change)')) + ylab(expression('-log'[10]*'(p-value)')) + 
  scale_y_continuous(expand = expansion(add = c(0, 10))) +
  ggtitle("Cluster 19 vs. cluster 1") +
  geom_label_repel(aes(label = label),
                   point.padding = 0.5,
                   segment.color = 'grey50', 
                   nudge_x = 0.5,
                   show_guide  = FALSE) +
  theme(axis.text = element_text(size = 18, colour = 'black'), axis.title = element_text(size = 20), 
        title = element_text(size = 18), legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))

dev.off()

png(
  filename  = 'paper_figure/Supp_Fig_2a.png',
  width     = 4.8,
  height    = 4,
  unit = 'in',
  res = 300
)
FeaturePlot(cancer.so, features = c("Gene701"), order = T, cols = c('grey90',"#D4524E"), raster = F, label = T) + 
  NoAxes() + ggtitle("Gene 70") + 
  theme(text = element_text(family = "Sans"), plot.title = element_text(face = 'plain', size = 18)) 
dev.off()

png(
  filename  = 'paper_figure/Supp_Fig_2b.png',
  width     = 7,
  height    = 4,
  unit = 'in',
  res = 300
)
VlnPlot(cancer.so, features = 'Gene701', pt.size = 0, group.by = 'pam50', cols = c('#EB8677','#7DC0A6','#919FC7','#Da8EC0','#8CBA54')) + 
  ggtitle("") + ylab("Module score of Gene70") + xlab("Clusters") + labs(fill = "PAM50\nclassification") + 
  scale_y_continuous(limits = c(-0.12, 0.4), breaks = seq(-0.1, 0.3, by = 0.1), expand = expansion(mult = c(0, 0))) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15), axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18)) + 
  stat_compare_means(comparisons = list(c('LumA','Normal'), c('Basal','Normal')), label = "p.val", size = 5) 
dev.off()