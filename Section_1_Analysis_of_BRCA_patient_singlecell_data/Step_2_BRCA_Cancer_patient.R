############################################################
### Analysis of cancer cells from breast cancer patients ###
############################################################

### Basic setting
set.seed(1234)
setwd('~/')

####################
### Load dataset ###
####################
brca.so <- readRDS('brca_SO.rds')

########################################
### Re-clustering breast cancer cell ###
########################################
### Here, we only focus on the cancer cell
Idents(brca.so) <- 'cell_type'
cancer.so <- subset(brca.so, idents = c('Cancer_cell'))

### store mitochondrial percentage in object meta data
cancer.so <- PercentageFeatureSet(cancer.so, pattern = "^MT-", col.name = "percent.mt")

cancer.so <- SCTransform(cancer.so, vars.to.regress = "percent.mt", verbose = FALSE, variable.features.n = 2000)
cancer.so <- RunPCA(cancer.so, npcs = 50, verbose = FALSE, assay.use = "SCT")
cancer.so <- RunHarmony(cancer.so, group.by.vars = 'patient', assay.use = "SCT")

### set npcs to 30
pcs <- 20
cancer.so <- FindNeighbors(cancer.so, reduction = "harmony", dims = 1:pcs)
cancer.so <- FindClusters(cancer.so, resolution = 0.7, method = 'igraph', random.seed = 1234, algorithm = 4)
cancer.so <- RunUMAP(cancer.so, reduction = "harmony", dims = 1:pcs)

############################
### PAM50 classification ###
############################
### We use genefu for PAM50 classification on breast cancer cell clusters
### Load PAM50 data
data(pam50.robust)

### Normalization
DefaultAssay(cancer.so) <- 'RNA'
cancer.so <- NormalizeData(cancer.so)
cancer.so <- ScaleData(cancer.so)

### Extract cluster average data
cancer.so.clavg <- AverageExpression(cancer.so, assays = 'RNA', group.by = 'seurat_clusters', layer = 'data')[['RNA']]

### Modify PAM50 gene list
pam50.gene <- pam50.robust[["centroids.map"]][["probe.centroids"]]
pam50.gene[11] <- 'NUF2'
pam50.gene[26] <- 'NDC80'
pam50.gene[40] <- 'ORC6'

# extract the 50 gene in cancer.so.clavg
sum(cancer.so@assays[["SCT"]]@data@Dimnames[[1]] %in% pam50.gene)
cancer.so.clavg <- as.data.frame(cancer.so.clavg[rownames(cancer.so.clavg) %in% pam50.gene,])

dannot <- as.data.frame(rownames(cancer.so.clavg))
colnames(dannot) <- 'Gene.Symbol'

### Gene symbol to EntrezGene ID
hs <- org.Hs.eg.db
my.symbols <- dannot$Gene.Symbol
for(i in 1:nrow(dannot)){
  dannot$EntrezGene.ID[i] <- AnnotationDbi::select(hs, my.symbols[i], keytype = 'SYMBOL', columns = "ENTREZID")$ENTREZID
}

annot <- as.matrix(dannot)
colnames(dannot) <- c('probe','EntrezGene.ID')
cancer.so.clavg.mtx <- as.matrix(t(cancer.so.clavg))

### PAM50 subtype prediction using genefu
SubtypePredictions <- molecular.subtyping(sbt.model = "pam50", data = cancer.so.clavg.mtx, annot = dannot, do.mapping = TRUE)

# sort the pam50 matrix data
SP <- as.data.frame(SubtypePredictions[['subtype']])
SP$cluster <- as.factor(rownames(SP))
SP <- as.data.frame(SP[order(SP$`SubtypePredictions[["subtype"]]`),])
colnames(SP) <- c('subtype', 'clusters')

### Put prediction results back to Seurat object
SP <- SP[order(SP$cluster),]
new.cluster.ids <- as.character(SP$subtype)
names(new.cluster.ids) <- rownames(SP)

### set identity as patient
cancer.so <- SetIdent(cancer.so, value = 'seurat_clusters')
cancer.so <- RenameIdents(cancer.so, new.cluster.ids)
cancer.so@meta.data[['pam50']] <- cancer.so@active.ident
cancer.so <- SetIdent(cancer.so, value = 'pam50')

### Save PAM50 annotated seurat object
saveRDS(cancer.so, file = 'cancer_SO.rds')

#####################
### Visualization ### 
#####################
png(
  filename  = 'paper_figure/Fig_2a.png',
  width     = 5,
  height    = 4,
  unit = 'in',
  res = 300
)

DimPlot(cancer.so, reduction = 'umap', raster = F, label = F) + labs(color = "Clusters of\ncancer cells") & NoAxes()

dev.off()


### plot heatmap of PAM50 genes expression level in each cancer cell cluster
cancer.so.clavg.mtx <- cancer.so.clavg.mtx[match(rownames(SP), rownames(cancer.so.clavg.mtx)),]
annot_SP <- as.data.frame(SP$subtype)
rownames(annot_SP) <- rownames(cancer.so.clavg.mtx)

mat <- scale(cancer.so.clavg.mtx)
row_ha = rowAnnotation(Subtype = SP$subtype,
                       col = list(Subtype = c("Basal" = "#EB8677", "Her2" = "#7DC0A6", "LumA" = "#919FC7",
                                              "LumB" = "#Da8EC0", "Normal" = "#8CBA54")),
                       gp = gpar(col = "white", lwd = 2),
                       annotation_name_gp = gpar(fontsize = 0))

mat2 <- Heatmap(mat, name = "z-score", rect_gp = gpar(col = "white", lwd = 2), cluster_rows = FALSE, 
                left_annotation = row_ha, row_title = 'Clusters', row_title_side = 'right', 
                column_title = "PAM50 genes", column_names_rot = 45, column_title_side = 'bottom',
                column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 12),
                column_title_gp = gpar(fontsize = 14), row_title_gp = gpar(fontsize = 14))

png(
  filename  = 'paper_figure/Fig_2b.png',
  width     = 10,
  height    = 6,
  unit = 'in',
  res = 300
)

mat2

dev.off()

png(
  filename  = 'paper_figure/Fig_2c.png',
  width     = 5,
  height    = 4,
  unit = 'in',
  res = 300
)

DimPlot(cancer.so, label = F, cols = c('#EB8677','#7DC0A6','#919FC7','#Da8EC0','#8CBA54')) + 
  labs(color = "PAM50\nclassification") & NoAxes()

dev.off()
