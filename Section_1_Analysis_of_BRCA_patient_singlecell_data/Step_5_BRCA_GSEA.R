####################################
### Gene set enrichment analysis ###
####################################

#################
### Load data ###
#################
cancer.so <- readRDS("cancer_SO.rds")

### We used AUCell for GSEA 
CellRank <- AUCell_buildRankings(as.matrix(cancer.so@assays$SCT@counts))

### Get human gene sets data
### Here, we used KEGG and Reactome gene sets
gs.go.kegg <- getGeneSets(species = 'Homo sapiens', library = 'C2', subcategory = 'KEGG')
gs.go.reactome <- getGeneSets(species = 'Homo sapiens', library = 'C2', subcategory = 'REACTOME')

### Calculate the gene set AUC score for KEGG gene set
cells_AUC_kegg <- AUCell_calcAUC(gs.go.kegg, CellRank, aucMaxRank=nrow(CellRank)*0.05, nCores = 10)
save(cells_AUC_kegg, file = 'cells_AUC_kegg.RData')

### Calculate the gene set AUC score for REACTOME gene set
cells_AUC_reactome <- AUCell_calcAUC(gs.go.reactome, CellRank, aucMaxRank = nrow(CellRank)*0.05, nCores = 10)
save(cells_AUC_reactome, file = 'cells_AUC_reactome.RData')

################################################
### Step 2, Exploring the enrichment results ###
################################################

### We first added the enrichment results to seurat object
cancer.so[['KEGG']] <- CreateAssayObject(cells_AUC_kegg@assays@data@listData[["AUC"]])
cancer.so[['Reactome']] <- CreateAssayObject(cells_AUC_reactome@assays@data@listData[["AUC"]])

cancer.so <- ScaleData(cancer.so, assay = "KEGG")
cancer.so <- ScaleData(cancer.so, assay = "Reactome")

### Next, we identified the activated pathways in malignant breast cancer cells (cluster 8 and 19) vs. normal-like cells (cluster 1)
### set seurat_clusters as active ident
cancer.so <- SetIdent(cancer.so, value = 'seurat_clusters') 

PE8_1_R <- FindMarkers(cancer.so, assay = 'Reactome', only.pos = T, logfc.threshold = 0.01, ident.1 = c(8), ident.2 = c(1))
PE8_1_K <- FindMarkers(cancer.so, assay = 'KEGG', only.pos = T, logfc.threshold = 0.01, ident.1 = c(8), ident.2 = c(1))

PE19_1_R <- FindMarkers(cancer.so, assay = 'Reactome', only.pos = T, logfc.threshold = 0.01, ident.1 = c(19), ident.2 = c(1))
PE19_1_K <- FindMarkers(cancer.so, assay = 'KEGG', only.pos = T, logfc.threshold = 0.01, ident.1 = c(19), ident.2 = c(1))

reorder_markers <- function(PE){
  PE <- PE[order(PE$p_val_adj),]
  PE$p_val_rank <- 1:nrow(PE)
  PE$name <- rownames(PE)
  
  return(PE)
}

PE8_1_R <- reorder_markers(PE8_1_R)
PE8_1_K <- reorder_markers(PE8_1_K)
PE19_1_R <- reorder_markers(PE19_1_R)
PE19_1_K <- reorder_markers(PE19_1_K)

over_reactome <- PE8_1_R$name[PE8_1_R$name %in% PE19_1_R$name]
PE_reactome <- cbind(PE8_1_R[PE8_1_R$name %in% over_reactome,1:5], PE19_1_R[PE19_1_R$name %in% over_reactome,1:5])
colnames(PE_reactome) <- c('p_val_c8','logFC_c8','pct_1_c8','pct_2_c8', 'p_val_adj_c8', 'p_val_c19','logFC_c19','pct_1_c19','pct_2_c19', 'p_val_adj_c19')

over_kegg <- PE8_1_K$name[PE8_1_K$name %in% PE19_1_K$name]
PE_kegg <- cbind(PE8_1_K[PE8_1_K$name %in% over_kegg,1:5], PE19_1_K[PE19_1_K$name %in% over_kegg,1:5])
colnames(PE_kegg) <- c('p_val_c8','logFC_c8','pct_1_c8','pct_2_c8', 'p_val_adj_c8', 'p_val_c19','logFC_c19','pct_1_c19','pct_2_c19', 'p_val_adj_c19')

### Save files
saveRDS(cancer.so, 'cancerSO_path.rds')
write.table(PE_reactome, 'PE_reactome.txt')
write.table(PE_kegg, 'PE_kegg.txt')

#####################
### Visualization ###
#####################
### Venn diagram (KEGG)
my_venndiag <- draw.pairwise.venn(area1 = 53, area2 = 50, cross.area = 17, category = c("C8", "C19"),
                                  fill = c("#EB8677", "#919FC7"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.05, 2),
                                  cat.cex = 4, cex = 4, fontfamily = 'Arial', cat.fontfamily = 'Arial')

png(
  filename  = 'paper_figure/Fig_2c-1.png',
  width     = 6,
  height    = 6,
  unit = 'in',
  res = 300
)

grid.arrange(gTree(children = my_venndiag), # Add title & subtitle
             top = textGrob("KEGG (vs. C1)", gp=gpar(fontsize=40)))

dev.off()

### Venn diagram (Reactome)
my_venndiag <- draw.pairwise.venn(area1 = 599, area2 = 428, cross.area = 137, category = c("C8", "C19"),
                                  fill = c("#EB8677", "#919FC7"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.05, 2),
                                  cat.cex = 4, cex = 4, fontfamily = 'Arial', cat.fontfamily = 'Arial')

png(
  filename  = 'paper_figure/Fig_2c-2.png',
  width     = 6,
  height    = 6,
  unit = 'in',
  res = 300
)

grid.arrange(gTree(children = my_venndiag), # Add title & subtitle
             top = textGrob("Reactome (vs. C1)", gp=gpar(fontsize=40)))

dev.off()

### VlnPlot of KEGG and Reactome 
png(
  filename  = 'paper_figure/Fig_2d-1.png',
  width     = 5.5,
  height    = 3.5,
  unit = 'in',
  res = 300
)

VlnPlot(cancer.so, features = c('KEGG-STEROID-BIOSYNTHESIS'),  ncol = 1, pt.size = 0, idents = c(1,8,19)) +
  ggtitle("KEGG (vs. C1)\nsteroid biosynthesis") + scale_fill_manual(values = c('#8CBA54','#EB8677','#919FC7')) +
  ylab("Expression\nlevel") + xlab("Clusters") + labs(fill = "Clusters") + 
  theme(text = element_text(family = "Sans"), axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18), axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20), legend.position = "none", plot.title = element_text(size = 22, face = 'plain'))

dev.off()

png(
  filename  = 'paper_figure/Fig_2d-2.png',
  width     = 5.5,
  height    = 3.5,
  unit = 'in',
  res = 300
)

VlnPlot(cancer.so, features = c('REACTOME-CHOLESTEROL-BIOSYNTHESIS'),  ncol = 1, pt.size = 0, idents = c(1,8,19)) + 
  ggtitle("Reactome (vs. C1)\ncholesterol biosynthesis") + scale_fill_manual(values = c('#8CBA54','#EB8677','#919FC7')) +
  ylab("Expression\nlevel") + xlab("Clusters") + labs(fill = "Clusters") + 
  theme(text = element_text(family = "Sans"), axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18), axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20), legend.position = "none", plot.title = element_text(size = 22, face = 'plain'))

dev.off()

### Dot plot of GSEA
### create dataset
bubble <- data.frame(matrix(nrow = 4, ncol = 2))
bubble$X1 <- as.numeric(c( 0.649, 0.837, 0.706,0.749))
bubble$X2 <- as.numeric(c( 9.718608e-20, 1.434073e-42, 1.315142e-25,1.668163e-03))

colnames(bubble) <- c('pct','pvaladj')
bubble$group <- as.factor(c(8,8,19,19))
bubble$db <- c('kegg', 'reactome','kegg','reactome')
bubble$pw <- c('KEGG\nsteroid biosynthesis', 'Reactome\ncholesterol biosynthesis','KEGG\nsteroid biosynthesis','Reactome\ncholesterol biosynthesis')
bubble$logpvaladj <- -log10(bubble$pvaladj)
bubble$order <- c()

png(
  filename  = 'paper_figure/Fig_2e.png',
  width     = 9.8,
  height    = 3.2,
  unit = 'in',
  res = 300
)

ggplot(bubble, aes(x=group, y=pw, size=pct, color=logpvaladj)) +
  geom_point(alpha=0.8) + theme_bw() + ylab('Pathways') + xlab('Clusters\n(vs. C1)') + 
  theme(text = element_text(family = 'Sans'),
        legend.direction = "horizontal", legend.box = "vertical", legend.justification = 'top',
        legend.text = element_text(size = 16), legend.title = element_text(size = 18),
        axis.text = element_text(size = 20, color = 'black'), axis.title = element_text(size = 24)) +
  scale_size(limits = c(0.4,1), breaks = seq(0.45, 0.8, by = 0.15), range = c(7, 21), name = "Percentage") +
  scale_colour_gradientn(name = expression("-log"[10]*'(p-val-adj)'),
                         colors = c('#919FC7', '#EB8677'), limits = c(0, 45)) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0, barwidth = 10),
         size = guide_legend(title.position="top", title.hjust = 0)) +
  scale_y_discrete(limits=rev)

dev.off()
