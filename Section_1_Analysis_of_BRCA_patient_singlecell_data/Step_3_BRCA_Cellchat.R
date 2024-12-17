#######################################################
### Cell-cell communication using CellChat analysis ###
#######################################################

### Cell-cell communication between breast cancer subtypes and other cell types

####################
### Load dataset ###
####################
brca.so <- readRDS('brca_SO.rds')
cancer.so <- readRDS("cancer_SO.rds")

##################################################################
### Adjust the metadata to breast cancer subtypes + cell types ###
##################################################################
subtypes <- as.character(brca.so@meta.data[['cell_type']])
names(subtypes) <- brca.so@assays[["SCT"]]@data@Dimnames[[2]]

pam50 <- as.character(cancer.so@meta.data[["pam50"]])
pam50 <- paste0('Cancer_cell_',pam50)
names(pam50) <- cancer.so@assays[['SCT']]@data@Dimnames[[2]]

### Add cancer cell subtypes information onto original cell types annotation
for(i in 1:length(subtypes)){
  if(names(subtypes[i]) %in% names(pam50)){
    subtypes[i] <- pam50[names(subtypes[i])]
  }
}

subtypes <- as.factor(subtypes)
brca.so@meta.data[['breast_subtypes']] <- subtypes
brca.so <- SetIdent(brca.so, value = 'breast_subtypes') 

### We then work on the breast_subtyeps and conduct the CellChat analysis
cellchat <- createCellChat(object = brca.so, group.by = 'breast_subtypes', assay = 'RNA')

### interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

future::plan("multicore", workers = 12) # do parallel
options(future.globals.maxSize= 1610612736)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# The number of highly variable ligand-receptor pairs used for signaling inference is 2478 

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

### Compute and visualize the network centrality scores
### Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

### Identify and visualize outgoing communication pattern of secreting cells
selectK(cellchat, pattern = "outgoing") 
nPatterns = 3 # According to cophenetic and silhouette score
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

selectK(cellchat, pattern = "incoming")
nPatterns = 7 # According to cophenetic and silhouette score
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
# Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
# Classification learning of the signaling networks for a single dataset

############################
### Save Cellchat object ###
############################
saveRDS(cellchat, file = "cellchat.rds")

#####################
### Visualization ###
#####################
png("paper_figure/Fig_2d.png",
    width = 21,
    height = 16, 
    units = 'cm',
    res = 600)

par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

dev.off()

png("paper_figure/Fig_2d-2.png",
    width = 21,
    height = 16, 
    units = 'cm',
    res = 600)

par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

dev.off()

# river plot
p <- netAnalysis_river(cellchat, pattern = "outgoing", font.size = 4.5, font.size.title = 26)
png(
  'paper_figure/Fig_2e.png',
  width = 42,
  height = 42,
  units = 'cm',
  res = 600
)
p
dev.off()

# river plot
p <- netAnalysis_river(cellchat, pattern = "incoming", font.size = 4, cutoff = 0.5)
#> Please make sure you have load `library(ggalluvial)` when running this function
png(
  'paper_figure/Fig_2f.png',
  width = 35,
  height = 35,
  units = 'cm',
  res = 600
)
p
dev.off()

mat <- cellchat@net$weight
png("paper_figure/Supp_Fig_1a.png",
    width = 31,
    height = 26, 
    units = 'cm',
    res = 600)
par(mfrow = c(3,4), xpd=TRUE, mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# Visualization in 2D-space
png(
  'paper_figure/Supp_Fig_1b.png',
  width = 15,
  height = 15,
  units = 'cm',
  res = 600
)
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
dev.off()

