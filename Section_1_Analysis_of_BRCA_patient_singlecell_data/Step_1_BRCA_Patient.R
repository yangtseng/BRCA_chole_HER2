###############################################################
### Analysis of 31 breast cancer patients' single-cell data ###
###############################################################

### Basic setting
set.seed(1234)
setwd('~/')

####################
### Load dataset ###
####################

brca <- readRDS("patient_dataset/1863-counts_cells_cohort1.rds")
### [25288 * 175942]
brca_biokey <- read.csv('patient_dataset/pre_treatment.csv')
### [9 * 175942]

#######################
### Quality control ###
#######################

### QC using scater pre-processing pipeline
brca.qc <- perCellQCMetrics(brca, subsets = list(Mito=grep("MT-", rownames(brca))))
### Remove outlier 3 nmad
ol <- isOutlier(metric = brca.qc$detected, nmad = 3, log = T)
### no outlier
### Remove mito_percent > 10%
mito_filter <- brca.qc$subsets_Mito_percent >= 10

# overlap between ol and mito_filter
filter <- ol | mito_filter
### 35512

brca <- brca[,!filter]
brca_biokey <- brca_biokey[!filter,]

### Remove genes expressed low
gene_filter <- rowSums(brca > 0) >= 10 ### 2206
cell_filter <- colSums(brca > 0) >= 100 ### 0

brca <- brca[gene_filter, ]
brca <- brca[ ,cell_filter]
brca_biokey <- brca_biokey[cell_filter,]

########################################
### Single-cell analysis with Seurat ###
########################################

### Create seurat object
brca.so <- CreateSeuratObject(brca, project = 'brca')

### Add clinical information
brca.so <- AddMetaData(brca.so, as.factor(brca_biokey$cellType), 'cell_type')
brca.so <- AddMetaData(brca.so, as.factor(brca_biokey$expansion), 'expansion')
brca.so <- AddMetaData(brca.so, as.factor(brca_biokey$patient_id), 'patient')
brca.so <- AddMetaData(brca.so, as.factor(brca_biokey$BC_type), 'BC_type')
brca.so <- AddMetaData(brca.so, as.factor(brca_biokey$timepoint), 'timepoint')

### store mitochondrial percentage in object meta data
brca.so <- PercentageFeatureSet(brca.so, pattern = "^MT-", col.name = "percent.mt")

brca.so <- SCTransform(brca.so, vars.to.regress = "percent.mt", verbose = FALSE, variable.features.n = 2000)
brca.so <- RunPCA(brca.so, npcs = 50, verbose = FALSE, assay.use = "SCT")
brca.so <- RunHarmony(brca.so, group.by.vars = 'patient', assay.use = "SCT")

### set npcs to 30
pcs <- 30
brca.so <- FindNeighbors(brca.so, reduction = "harmony", dims = 1:pcs)
brca.so <- FindClusters(brca.so, resolution = 0.5, method = 'igraph', random.seed = 1234, algorithm = 4)
brca.so <- RunUMAP(brca.so, reduction = "harmony", dims = 1:pcs)

brca.so <- SetIdent(brca.so, value = 'seurat_clusters') 
markers <- FindAllMarkers(brca.so)

### Save RDS file
saveRDS(brca.so, file = "brca_SO.rds")

#####################
### Visualization ###
#####################

png(
  filename  = 'paper_figure/Fig_1a.png',
  width     = 5,
  height    = 4,
  unit = 'in',
  res = 300
)

DimPlot(brca.so, raster = F) + labs(color = "Clusters") & NoAxes()

dev.off()

png(
  filename  = 'paper_figure/Fig_1b.png',
  width     = 10,
  height    = 4,
  unit = 'in',
  res = 300
)

FeaturePlot(brca.so, features = c('EPCAM','KRT18','KRT19'), ncol = 3, order = F, cols = c('grey80',"#D4524E"), raster = F) & 
  NoAxes() & NoLegend() & theme(text = element_text(family = "Sans"),  plot.title = element_text(size = 20, face = 'italic'))

dev.off()

png(
  filename  = 'paper_figure/Fig_1c.png',
  width     = 10,
  height    = 24,
  unit = 'in',
  res = 300
)

FeaturePlot(brca.so, features = c('CD19','IGHG1','IGKC','CD3D','CD3E','CD3G', 'CD33','CD68','ITGAX',
                                  'COL1A2','DCN','FAP','CD34','CLDN5','PECAM1','CLEC4C','IRF7','IRF8',
                                  'TPSAB1','TPSB2','MS4A2'), ncol = 3, cols = c('grey80',"#D4524E"), raster = F) & NoAxes() & NoLegend() & 
  theme(text = element_text(family = 'Sans'), title = element_text(size = 24, face = 'bold.italic'))

dev.off()

brca.so <- SetIdent(brca.so, value = 'cell_type') 

png(
  filename  = 'paper_figure/Fig_1d.png',
  width     = 5.5,
  height    = 4,
  unit = 'in',
  res = 300
)
DimPlot(brca.so, raster = F) + scale_color_brewer(palette='Set2') + labs(color = "Cell types") & NoAxes()
dev.off()

png(
  filename  = 'paper_figure/Fig_1e.png',
  width     = 8,
  height    = 12,
  unit = 'in',
  res = 300
)

dittoBarPlot(
  object = brca.so,
  var = brca.so@meta.data$cell_type,
  group.by = "patient") + coord_flip() + xlab("Patients") +
  theme(text = element_text(family = 'Sans'), legend.position = "bottom", legend.justification = 'right',
        axis.title = element_text(size = 26), axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 20), legend.text = element_text(size = 20),
        legend.title = element_text(size = 24)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(fill=guide_legend(ncol=3, title.position="top")) +
  scale_fill_brewer(palette="Set2", name = 'Cell types')

dev.off()
