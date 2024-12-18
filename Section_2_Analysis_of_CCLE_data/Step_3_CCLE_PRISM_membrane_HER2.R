#############################################################
### CCLE with PRISM score and Cell surface HER2 Intensity ###
#############################################################

# We took the rank of cholesterol biosynthesis from GSVA results using CCLE transcriptomics data
# We first used PRISM score of Neratinib (HER2 inhibitor)
#################
### Load data ###
#################
dt.prism <- read.csv('primary-screen-replicate-collapsed-logfold-change.csv', header = T, row.names = 1)
drug.info <- read.csv('primary-screen-replicate-treatment-info.csv', header = T, row.names = 1)
# remove NA in drug.info
drug.info <- drug.info[complete.cases(drug.info$name),]

### BRD-K85606544-001-09-1 Neratinib
drug.info <- drug.info[drug.info$broad_id %in% 'BRD-K85606544-001-09-1',]

### we collect the nonHER2 breast cancer cell lines manually from PRISM database
dt.nonHER2 <- c('ACH-000849', #MDAMB468
                'ACH-000019', #MCF7, LumA
                'ACH-000288', #BT549
                'ACH-000147', #T47D, LumA
                'ACH-000223', #HCC1937
                'ACH-000856', #CAL51
                'ACH-000783', #CAMA1, LumB
                'ACH-000276', #HCC38
                'ACH-000573', #MDAMB436
                'ACH-000374') #HCC1143

dt.prism <- dt.prism[rownames(dt.prism) %in% dt.nonHER2, ]
dt.prism.drug <- colnames(dt.prism)
dt.prism.drug <- substr(dt.prism.drug, 1, 22)
dt.prism.drug <- gsub("\\.", '-', dt.prism.drug)
colnames(dt.prism) <- dt.prism.drug
dt.prism <- t(dt.prism)
dt.prism <- dt.prism[rownames(dt.prism) %in% 'BRD-K85606544-001-09-1',]
dt.prism <- as.data.frame(dt.prism)
dt.nonHER2_cellline_name <- c("MCF7",'T47D','HCC1937','HCC38','BT549','HCC1143','MDAMB436','CAMA1','MDAMB468','CAL51')
rownames(dt.prism) <- dt.nonHER2_cellline_name
dt.prism$cellline <- rownames(dt.prism)

### Ranked by cholesterol biosynthesis
dt.prism$rank <- c(3,4,5,7,1,9,10,8,2,6)
dt.prism$types <- c(rep("ER- HER2-", 10))
dt.prism$types[1] <- c("ER+ HER2-")
dt.prism$types[2] <- c("ER+ HER2-")
dt.prism$types[8] <- c("ER+ HER2-")
dt.prism.new <- dt.prism[c(1,2,3,4,5,6,7,8,9,10),]
fit <- lm(rank ~ dt.prism, data = dt.prism.new) 
summary(fit)

###################################
### Cell surface HER2 intensity ###
###################################
chole_her2 <- data.frame(chole=c(0.710063,0.6037967,0.2987737,-1.009477), her2=c(5.1,12.1,26.8,21.8), 
                         rank=c(1,2,3,4), cellline = c("MDAMB468", 'MCF7','T47D','MDAMB231'),
                         types = c('ER- HER2-', 'ER+ HER2-','ER+ HER2-','ER- HER2-'))
fit <- lm(rank ~ her2, data = chole_her2) 
summary(fit)

#####################
### Visualization ###
#####################
p <- ggplot(dt.prism.new, aes(x = rank, y = dt.prism)) + geom_point(aes(colour = types)) + 
  xlab("Rank of\nCholesterol Biosynthesis") + ylab("PRISM score\n Neratinib") + theme_bw() + 
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
  geom_smooth(aes(rank, dt.prism), method = "lm",se = FALSE,color='grey80') +
  geom_label_repel(aes(label = cellline, colour = types),
                   point.padding = 0.5,
                   min.segment.length = 0, 
                   box.padding = 1,
                   segment.color = 'grey80', 
                   show_guide  = F,
                   nudge_x = 0.5,
                   xlim = c(1,10),
                   nudge_y = 0.5) +
  scale_colour_manual(values = c('#7BBDA3','#909EC6')) +
  scale_x_continuous(breaks = c(1:10))

png(
  filename  = 'paper_figure/Fig_6b.png',
  width     = 6,
  height    = 4,
  unit = 'in',
  res = 300
)
p+ annotate("text",label = "R = 0.256", x = 9, y = -1.3)
dev.off()

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

png(
  filename  = 'paper_figure/Fig_6c.png',
  width     = 5.2,
  height    = 4,
  unit = 'in',
  res = 300
)
p + annotate("text",label = "R = 0.861", x = 1.4, y = 25)
dev.off()
