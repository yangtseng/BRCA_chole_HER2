###################################################################################
### Cholesterol biosynthesis ranking in breast cancer cell lines data from CCLE ###
###################################################################################

#################
### Load data ###
#################
CCLE <- read.csv('~/R/BRCA/brand_new/CCLE_new/Expression_Public_23Q4_subsetted.csv', row.names = 1)
CCLE_meta <- read.table('CCLE_metadata.txt')
rownames(CCLE) <- CCLE$cell_line_display_name
CCLE <- CCLE[,-c(1:8)]

############
### GSVA ###
############
React_Chole_biosyn <- c('ACAT2',	'ARV1','CYP51A1',	'DHCR24',	'DHCR7',	'EBP',	'FDFT1',
                        'FDPS',	'GGPS1',	'HMGCR',	'HMGCS1',	'HSD17B7',	'IDI1',	'IDI2',
                        'LBR',	'LSS',	'MSMO1',	'MVD',	'MVK',	'NSDHL',	'PLPP6',	'PMVK',
                        'SC5D',	'SQLE',	'SREBF1',	'SREBF2',	'TM7SF2')
KEGG_steroid <- c('CEL',	'CYP27B1',	'CYP51A1',	'DHCR24',	'DHCR7',	'EBP',	'FDFT1',	
                  'HSD17B7',	'LIPA',	'LSS',	'MSMO1',	'NSDHL',	'SC5D',	'SOAT1',
                  'SOAT2',	'SQLE',	'TM7SF2')

### Separate cell lines into basal and lumA
basal_P <- CCLE_meta$name[grepl("Basal", CCLE_meta$lineage_6) & grepl("ERneg HER2neg", CCLE_meta$lineage_5)]
lum_P <- CCLE_meta$name[grepl("Luminal", CCLE_meta$lineage_6)]
non_HER2 <- c(basal_P, lum_P)

CCLE_nonHER2 <- CCLE[rownames(CCLE) %in% non_HER2,] #24

### GSVA
gsva.chole_react <- gsva(as.matrix(t(CCLE_nonHER2)), list(React_Chole_biosyn), kcdf = 'Poisson', verbose=FALSE)

gsva.chole_kegg <- gsva(as.matrix(t(CCLE_nonHER2)), list(KEGG_steroid), kcdf = 'Poisson', verbose=FALSE)
### Combine two column
gsva.chole <- cbind(t(gsva.chole_react), t(gsva.chole_kegg))
gsva.chole <- as.data.frame(gsva.chole)
colnames(gsva.chole) <- c('Reactome','KEGG')
gsva.chole <- gsva.chole[order(gsva.chole$Reactome, decreasing = T),]

#####################
### Visualization ###
#####################
### Draw heatmap by rank
p <- Heatmap(gsva.chole, rect_gp = gpar(col = 'white',lwd = 2), 
             # column_split = factor(avg$markertype, levels=unique(avg$markertype)), 
             # top_annotation = column_ha, 
             column_title = 'Cholesterol\nbiosynthesis',
             cluster_rows = F, cluster_columns = F, column_names_rot = 45, column_names_centered = F,
             heatmap_height = unit(0.8, "cm")*24, 
             heatmap_width = unit(1.5, "cm")*4, 
             row_names_gp = gpar(fontsize = 16),
             column_names_gp = gpar(fontsize = 16),
             column_title_gp = gpar(fontsize = 16),
             name = "GSVA\nscore",    
             col = colorRamp2(c(-1, 0, 1), c("#919FB7", "white", "#E58C82")),
             heatmap_legend_param = list(legend_direction = "vertical", 
                                         legend_height = unit(4, "cm")), 
             show_heatmap_legend = T)
png(
  filename  = 'paper_figure/Fig_6a.png',
  width     = 4,
  height    = 9,
  unit = 'in',
  res = 300
)
p
dev.off()