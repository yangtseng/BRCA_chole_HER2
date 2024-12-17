###################################################################
### Analysis of CCLE Proteomics data (breast cancer cell lines) ###
###################################################################

#################
### Load data ###
#################
CCLE_P <- read.csv('~/R/BRCA/brand_new/CCLE_proteomics_scaled.csv')
CCLE <- read.csv('~/R/BRCA/brand_new/CCLE_new/Expression_Public_23Q4_subsetted.csv')
rownames(CCLE) <- CCLE$cell_line_display_name
### Extract metadata from CCLE transcriptomics dataset
CCLE_meta <- CCLE[,c(1:8)]

###########################
### Data pre-processing ###
###########################
### Remove NA 
CCLE_P <- na.omit(CCLE_P) ### 7365*87

rownames(CCLE_P) <- CCLE_P$ProteinId
CCLE_Protein_meta <- CCLE_P[,c(1:3)] ### ProteinID and Gene information
CCLE_P <- CCLE_P[,c(4:87)]

### Remove bridges from the dataset
### Bridge: Mixed sample of six cell lines (HCC1806, Hs578T, MCF7, MCF 10A, MDAMB231, SKBR3)
### Thus, it was possible to compare sets to each other
CCLE_P <- CCLE_P[,!grepl("Bridge",colnames(CCLE_P))]

### manually assign BRCA subtypes of each cell line
CCLE_P_cellline <- data.frame(Cellline = colnames(CCLE_P))
CCLE_P_cellline$Cell <- gsub("[.]", "-", CCLE_P_cellline$Cellline)
CCLE_P_cellline$Cell <- gsub("--", " (", CCLE_P_cellline$Cell)
CCLE_P_cellline$Cell <- ifelse(grepl("repeat-", CCLE_P_cellline$Cell), sub("-([^-]*)$",")\\1", CCLE_P_cellline$Cell, 2), CCLE_P_cellline$Cell)
CCLE_P_cellline$name <- sub(" .*", "", CCLE_P_cellline$Cell)

### Annotating the CCLE cell line types
### manually assign BRCA status
CCLE_meta$lineage_5[1] <- "ERpos HER2neg"
CCLE_meta$lineage_5[4] <- "ERneg HER2neg"
CCLE_meta$lineage_5[8] <- "Non Cancerous"
CCLE_meta$lineage_5[23] <- "ERneg HER2pos"
CCLE_meta$lineage_5[24] <- "ERneg HER2pos"

### Simple modification
colnames(CCLE_meta)[2] <- "name"
CCLE_meta$name[34] <- "ZR-75-1"
CCLE_meta$name[25] <- "ZR-75-30"
CCLE_meta$name[49] <- "AU-565"
CCLE_meta$name[2] <- "BT-483"
CCLE_meta$name[3] <- "CAL85-1"
CCLE_meta$name[57] <- "BT-474"
CCLE_meta$name[28] <- "EFM-19"
CCLE_meta$name[15] <- "SUM1315"
CCLE_meta$name[16] <- "SUM149"
CCLE_meta$name[17] <- "SUM159"

CCLE_P_cellline[3,] <- c('Hs578T', "HS578T","HS578T")
CCLE_P_cellline[43,] <- c('MDA.MB.175VII','MDAMB175VII','MDAMB175VII')
CCLE_P_cellline[44,] <- c('MDA.MB.415','MDAMB415','MDAMB415')
CCLE_P_cellline[50,] <- c('MDA.MB.330','MDAMB330','MDAMB330')
CCLE_P_cellline[53,] <- c('UACC.812','UACC812','UACC812')

CCLE_P_cellline <- left_join(CCLE_P_cellline, CCLE_meta)

### Fill in the information manually
CCLE_P_cellline[4,4] <- "ACH-001357"
CCLE_P_cellline[32,4] <- "ACH-001357"
CCLE_P_cellline[42,4] <- "ACH-001357"
CCLE_P_cellline[60,4] <- "ACH-001357"
CCLE_P_cellline[61,4] <- "ACH-001357"
CCLE_P_cellline[4,8] <- "Non cancerous"
CCLE_P_cellline[32,8] <- "Non cancerous"
CCLE_P_cellline[42,8] <- "Non cancerous"
CCLE_P_cellline[60,8] <- "Non cancerous"
CCLE_P_cellline[61,8] <- "Non cancerous"

CCLE_P_cellline[31,8] <- "Non cancerous"
CCLE_P_cellline[56,8] <- "Non cancerous"
CCLE_P_cellline[57,8] <- "Non cancerous"
CCLE_P_cellline[62,8] <- "Non cancerous"

### Save the cell lines with correct annotation
CCLE_P_cellline_final <- CCLE_P_cellline[,c(1,2,3,8)]
CCLE_P_cellline_final <- na.omit(CCLE_P_cellline_final)
### 59

### re-construct the normalized, scaled proteomics matrix with annotated celllines
CCLE_P <- CCLE_P[,colnames(CCLE_P) %in% CCLE_P_cellline_final$Cellline]

### Save files
write.table(CCLE_P_cellline, 'CCLE_metadata.txt')

###########
### PCA ###
########### 
mat <- as.matrix(CCLE_P)
# mat <- as.matrix(CCLE[,colnames(CCLE) %in% gene70])
class(mat) <- "numeric"
mat <- mat[!duplicated(rownames(mat)), !duplicated(colnames(mat))]
p <- pca(mat)
mat.pca <- cbind(p$rotated$PC1, p$rotated$PC2)
mat.pca <- as.data.frame(mat.pca)
rownames(mat.pca) <- p$yvars
celllines <- c("BT20",'HCC38','BT.474', 'MDAMB231','T47D', "MCF7")
mat.pca$gene <- ifelse(rownames(mat.pca) %in% celllines, rownames(mat.pca), NA)
mat.pca$type<- CCLE_P_cellline_final$lineage_5
for(i in 1:59){
  if(mat.pca$type[i] == 'ERpos HER2pos'){
    mat.pca$Types[i] <- "ER+ HER2+"
  }else if(mat.pca$type[i] == 'ERpos HER2neg'){
    mat.pca$Types[i] <- "ER+ HER2-"
  }else if(mat.pca$type[i] == 'ERneg HER2neg'){
    mat.pca$Types[i] <- "ER- HER2-"
  }else if(mat.pca$type[i] == 'ERneg HER2pos'){
    mat.pca$Types[i] <- "ER- HER2+"
  }else{mat.pca$Types[i] <- "Non cancerous"}
}

#########################################
### Corr between HER2 and Cholesterol ###
#########################################

### Separate cell lines into basal and lumA
basal_P <- CCLE_P_cellline$Cellline[grepl("Basal", CCLE_P_cellline$lineage_6) & grepl("ERneg HER2neg", CCLE_P_cellline$lineage_5)]
lum_P <- CCLE_P_cellline$Cellline[grepl("Luminal", CCLE_P_cellline$lineage_6)]

basal_P <- CCLE_P[,colnames(CCLE_P) %in% basal_P]
lum_P <- CCLE_P[,colnames(CCLE_P) %in% lum_P]

### Take Cholesterol biosynthesis genes and ERBB family members
genes <- c('sp|P04035-3|HMDH_HUMAN', #HMGCR
           'sp|Q03426|KIME_HUMAN', #MVK
           'sp|P53602|MVD1_HUMAN', #MVD
           'sp|O95749|GGPPS_HUMAN', #GGPS1
           'sp|P37268|FDFT_HUMAN', #FDFT1
           'sp|Q14534|ERG1_HUMAN', #SQLE
           'sp|P00533|EGFR_HUMAN',
           'sp|P04626|ERBB2_HUMAN', 
           'sp|P21860|ERBB3_HUMAN',
           'sp|Q15303|ERBB4_HUMAN')

### Basal (n=24)
mat_basal <- as.matrix(basal_P[rownames(basal_P) %in% genes,])

### change rownames to leave only protein names
rownames(mat_basal) <- c('GGPPS','EGFR','HMDH','HER2','ERBB3','FDFT','MVD1','KIME','ERG1','ERBB4')
protein_level <- c("HMDH",'KIME','MVD1','GGPPS','FDFT','ERG1','EGFR','HER2','ERBB3','ERBB4')

cormat_basal <- round(cor(t(mat_basal), method = c("spearman")), 2)
cormat_basal <- cormat_basal[order(match(rownames(cormat_basal), protein_level)), order(match(colnames(cormat_basal), protein_level))]

### Luminal (n=14)
mat_lum <- as.matrix(lum_P[rownames(lum_P) %in% genes,])

### change rownames to leave only protein names
rownames(mat_lum) <- c('GGPPS','EGFR','HMDH','HER2','ERBB3','FDFT','MVD1','KIME','ERG1','ERBB4')
cormat_lum <- round(cor(t(mat_lum), method = c("spearman")), 2)
cormat_lum <- cormat_lum[order(match(rownames(cormat_lum), protein_level)), order(match(colnames(cormat_lum), protein_level))]

#####################
### Visualization ###
#####################
### PCA
p <- ggplot(mat.pca, aes(x = V1, y = V2, label = gene, colour = Types)) + geom_point() + 
  xlab("PC 1") + ylab("PC 2") + theme_bw() + 
  theme(text = element_text(family = "Sans", size = 14),
        axis.text = element_text(size = 20, colour = 'black'),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = 'top',
        legend.justification = 'right',
        legend.margin=margin(0, 0, -10, 0),
        legend.spacing.x = unit(0.3, 'mm'),
        legend.spacing.y = unit(0.3, 'mm'),
        panel.border = element_rect(size = 0.8, colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  geom_label_repel(aes(label = gene),
                   point.padding = 0.5,
                   min.segment.length = 0, 
                   box.padding = 1,
                   segment.color = 'grey80', 
                   show_guide  = F,
                   nudge_x = 0.5,
                   xlim = c(-200,NA),
                   nudge_y = 0.5) +
  scale_colour_brewer(palette = "Set2") +
  guides(colour=guide_legend(ncol=2,byrow=TRUE))

png(
  filename  = 'paper_figure/Fig_5b.png',
  width     = 4.8,
  height    = 5,
  unit = 'in',
  res = 300
)
p
dev.off()

###############
### Heatmap ###
###############
### Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri_basal <- get_upper_tri(cormat_basal)
upper_tri_lum <- get_upper_tri(cormat_lum)


ht_opt(RESET = T)

draw_heatmap <- function(upper_tri){
  col_fun = colorRamp2(c(-1,0, 1), col = c("#A5B6C5","white", "#F0988C"))
  row_ha = rowAnnotation(Types = c(rep('Cholesterol biosynthesis\nrelated proteins', 6), rep('ERBB family', 4)),
                         col = list(Types = c("Cholesterol biosynthesis\nrelated proteins" = "#FEE699", 
                                              "ERBB family" = "#EA8576")),
                         show_annotation_name = F,
                         simple_anno_size = unit(2, "mm"),
                         annotation_legend_param = list(title = "Types", 
                                                        title_gp = gpar(fontsize = 10), 
                                                        labels_gp = gpar(fontsize = 8)))
  column_ha = HeatmapAnnotation(Types = c(rep('Cholesterol biosynthesis\nrelated proteins', 6), rep('ERBB family', 4)), 
                                col = list(Types = c("Cholesterol biosynthesis\nrelated proteins" = "#FEE699", 
                                                     "ERBB family" = "#EA8576")),
                                show_annotation_name = F,
                                simple_anno_size = unit(2, "mm"),
                                annotation_legend_param = list(title = "Types", 
                                                               title_gp = gpar(fontsize = 10), 
                                                               labels_gp = gpar(fontsize = 8)))
  
  mat2 <- Heatmap(upper_tri, na_col = 'white',
                  heatmap_legend_param = list(title = 'Spearman correlation',
                                              legend_width = unit(3, "cm"),
                                              direction = "horizontal",
                                              labels_gp = gpar(fontsize = 8),
                                              title_gp = gpar(fontsize = 10),
                                              title_position = "topleft",
                                              border = "black"),
                  col = col_fun, 
                  top_annotation = column_ha,
                  right_annotation = row_ha,
                  rect_gp = gpar(col = "white", lwd = 1), 
                  cluster_columns = F, cluster_rows = F,
                  column_names_rot = 45, 
                  column_names_side = 'top',
                  column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 10),
                  column_title_gp = gpar(fontsize = 14), row_title_gp = gpar(fontsize = 14), 
                  cell_fun = function(i, j, x, y, width, height, fill){
                    if(is.na(upper_tri[j,i]) == T){
                      grid.text(sprintf("%.1f", upper_tri[j, i]), x, y, gp = gpar(fontsize = 20, col = NA))
                    }else{
                      grid.text(sprintf("%.1f", upper_tri[j, i]), x, y, gp = gpar(fontsize = 8))}})
  return(mat2)
}

hmp_lum <- draw_heatmap(upper_tri_lum)
png(
  filename  = 'paper_figure/Fig_5c.png',
  width     = 6,
  height    = 3.5,
  unit = 'in',
  res = 300
)
draw(hmp_lum, merge_legends = T, padding = unit(c(2, 10, 2, 2), "mm"))
dev.off()

hmp_basal <- draw_heatmap(upper_tri_basal)
png(
  filename  = 'paper_figure/Fig_5e.png',
  width     = 6,
  height    = 3.5,
  unit = 'in',
  res = 300
)
draw(hmp_basal, merge_legends = T, padding = unit(c(2, 10, 2, 2), "mm"))
dev.off()

###################################
### Protein protein correlation ### 
###################################

### Try basal_P (n=24) and lum_P (n=14)

### 'sp|P04035-3|HMDH_HUMAN', #HMGCR
### 'sp|P37268|FDFT_HUMAN', #FDFT1
### 'sp|Q14534|ERG1_HUMAN', #SQLE
### 'sp|P04626|ERBB2_HUMAN' #HER2

### Luminal (n=14)
df <- as.data.frame(t(lum_P))
p1 <- ggplot(df, aes(x = `sp|P04035-3|HMDH_HUMAN`, y = `sp|P04626|ERBB2_HUMAN`)) + 
  geom_point() + 
  xlab(expression("HMDH")) + 
  ylab(expression("HER2")) + 
  theme_bw() + 
  theme(text = element_text(family = "Sans", size = 10),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = 'top',
        legend.justification = 'right',
        legend.direction = 'horizontal',
        legend.margin=margin(0, 0, -10, 0),
        legend.spacing.x = unit(0.3, 'mm'),
        legend.spacing.y = unit(0.3, 'mm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = c('#949FC4', '#DD8B7B')) +
  geom_smooth(method=lm, se=T, color = "black", linetype = "dashed") +
  stat_cor(method = "spearman", label.x = 5, label.y = -5, size = 4) 

p1
png(
  filename  = 'paper_figure/Fig_7d.png',
  width     = 3.5,
  height    = 3,
  unit = 'in',
  res = 300
)
p1
dev.off()

### Basal (n=24)
df <- as.data.frame(t(basal_P))
p1 <- ggplot(df, aes(x = `sp|P37268|FDFT_HUMAN`, y = `sp|P04626|ERBB2_HUMAN`)) + 
  geom_point() + 
  xlab(expression("FDFT")) + 
  ylab(expression("HER2")) + 
  theme_bw() + 
  theme(text = element_text(family = "Sans", size = 10),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = 'top',
        legend.justification = 'right',
        legend.direction = 'horizontal',
        legend.margin=margin(0, 0, -10, 0),
        legend.spacing.x = unit(0.3, 'mm'),
        legend.spacing.y = unit(0.3, 'mm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = c('#949FC4', '#DD8B7B')) +
  geom_smooth(method=lm, se=T, color = "black", linetype = "dashed") +
  stat_cor(method = "spearman", label.x = 13, label.y = 7, size = 4) 

png(
  filename  = 'paper_figure/Fig_7f.png',
  width     = 3.5,
  height    = 3,
  unit = 'in',
  res = 300
)
p1
dev.off()

p1 <- ggplot(df, aes(x = `sp|Q14534|ERG1_HUMAN`, y = `sp|P04626|ERBB2_HUMAN`)) + 
  geom_point() + 
  xlab(expression("ERG1")) + 
  ylab(expression("HER2")) + 
  theme_bw() + 
  theme(text = element_text(family = "Sans", size = 10),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = 'top',
        legend.justification = 'right',
        legend.direction = 'horizontal',
        legend.margin=margin(0, 0, -10, 0),
        legend.spacing.x = unit(0.3, 'mm'),
        legend.spacing.y = unit(0.3, 'mm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = c('#949FC4', '#DD8B7B')) +
  geom_smooth(method=lm, se=T, color = "black", linetype = "dashed") +
  stat_cor(method = "spearman", label.x = 13, label.y = 7, size = 4) 

png(
  filename  = 'paper_figure/Fig_7g.png',
  width     = 3.5,
  height    = 3,
  unit = 'in',
  res = 300
)
p1
dev.off()

### HER2 Surface w/ cholesterol biosynthesis
# MDAMB468:0.710063
# MCF7:0.6037967
# T47D:0.2987737
# MDAMB231:-1.009477

chole_her2 <- data.frame(chole=c(0.710063,0.6037967,0.2987737,-1.009477), her2=c(5.1,12.1,26.8,21.8), 
                         rank=c(1,2,3,4), cellline = c("MDAMB468", 'MCF7','T47D','MDAMB231'),
                         types = c('ER- HER2-', 'ER+ HER2-','ER+ HER2-','ER- HER2-'))
ggplot(chole_her2, aes(x = rank, y = her2, label = cellline)) + geom_point() + geom_text(hjust=0.5, vjust=-1)
library(ggrepel)
p <- ggplot(chole_her2, aes(x = rank, y = her2)) + geom_point(aes(colour = types)) + 
  xlab("Rank of\nCholesterol Biosynthesis") + ylab("Cell Surface HER2\nintensity (MFI)") + theme_bw() + 
  theme(text = element_text(family = "Sans", size = 14),
        axis.text = element_text(size = 20, colour = 'black'),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = 'top',
        legend.justification = 'right',
        legend.margin=margin(0, 0, -10, 0),
        legend.spacing.x = unit(0.3, 'mm'),
        legend.spacing.y = unit(0.3, 'mm'),
        panel.border = element_rect(size = 0.8, colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  geom_smooth(aes(rank, her2), method = "lm",se = FALSE,color='grey80') +
  geom_label_repel(aes(label = cellline, colour = types),
                   point.padding = 0.5,
                   min.segment.length = 0, 
                   box.padding = 1,
                   segment.color = 'grey80', 
                   show_guide  = F,
                   nudge_x = 0.5,
                   xlim = c(1,4),
                   nudge_y = 0.5) +
  scale_colour_manual(values = c('#7BBDA3','#909EC6'))
p + annotate("text", label = expression("R"^2), x= 1.2, y=25) + annotate("text",label = "= 0.741", x = 1.5, y = 25)
png(
  filename  = 'surfaceHER2.png',
  width     = 5.2,
  height    = 4,
  unit = 'in',
  res = 300
)
p+ annotate("text", label = expression("R"^2), x= 1.2, y=25) + annotate("text",label = "= 0.741", x = 1.5, y = 25)
dev.off()
fit <- lm(rank ~ her2, data = chole_her2) 
summary(fit)

### PRISM
