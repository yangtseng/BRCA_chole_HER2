#############################
### synergy finder MDA NB ###
#############################

#################
### Load data ###
#################

#############################
### synergy finder MDA NB ###
#############################
MDA_NB <- read.table('Colony/MDA_NB.txt', header = 1)
colnames(MDA_NB) <- c('BlockId','Row','Col','Response','Replicate','DrugRow','DrugCol','ConcRow','ConcCol','conc_r_unit','conc_c_unit')

### Reshape the matrix
MDA_NB_reshape <- ReshapeData(MDA_NB, data_type = "viability")

PlotDoseResponse(MDA_NB_reshape)
MDA_NB_score <- CalculateSynergy(data = MDA_NB_reshape, method = "ZIP", correct_baseline = 'all',
                                 Emin = 0, Emax = 100)

#################################
### synergy finder MDA statin ###
#################################
MDA_statin <- read.table('Colony/MDA_statin.txt', header = 1)
colnames(MDA_statin) <- c('BlockId','Row','Col','Response','Replicate','DrugRow','DrugCol','ConcRow','ConcCol','conc_r_unit','conc_c_unit')

MDA_statin_reshape <- ReshapeData(MDA_statin, data_type = "viability")

PlotDoseResponse(MDA_statin_reshape)
MDA_statin_score <- CalculateSynergy(data = MDA_statin_reshape, method = "ZIP", correct_baseline = 'all',
                                 Emin = 0, Emax = 100)

##############################
### synergy finder T47D NB ###
##############################
T47D_NB <- read.table('Colony/T47D_NB.txt', header = 1)
colnames(T47D_NB) <- c('BlockId','Row','Col','Response','Replicate','DrugRow','DrugCol','ConcRow','ConcCol','conc_r_unit','conc_c_unit')
T47D_NB_reshape <- ReshapeData(T47D_NB)

PlotDoseResponse(T47D_NB_reshape)
T47D_NB_response <- T47D_NB_reshape$response[T47D_NB_reshape$response$block_id == 1, c("conc1", "conc2", "response")]
T47D_NB_ZIP <- ZIP(T47D_NB_response)

T47D_NB_score <- CalculateSynergy(data = T47D_NB_reshape, method = "ZIP", correct_baseline = 'all',
                                  Emin = 0, Emax = 100)


##################################
### synergy finder T47D statin ###
##################################
T47D_statin <- read.table('Colony/T47D_statin.txt', header = 1)
colnames(T47D_statin) <- c('BlockId','Row','Col','Response','Replicate','DrugRow','DrugCol','ConcRow','ConcCol','conc_r_unit','conc_c_unit')
T47D_statin_reshape <- ReshapeData(T47D_statin, data_type = "viability")

PlotDoseResponse(T47D_statin_reshape)
T47D_statin_score <- CalculateSynergy(data = T47D_statin_reshape, method = "ZIP", correct_baseline = 'all',
                                      Emin = 0, Emax = 100)

#####################
### Visualization ###
#####################
### MDAMB231 NB598 + Neratinib
p <- PlotSynergy(MDA_NB_score, type = "2D", grid = F)

p1 <- p$`1` + 
  theme(text = element_text(family = 'Sans'), 
        axis.text = element_text(size = 20, color = 'black'), 
        axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 0, b = 0, l = 0)), 
        axis.title.x = element_text(size = 22, angle = 180, margin = margin(t = 20, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
        legend.title = element_text(size = 18, hjust = 0.5, angle = 90),
        legend.position = 'right',
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(0.4,'cm'),
        legend.key = element_rect(color="black"),
        plot.title = element_blank(), 
        plot.subtitle = element_blank(),
        plot.margin = margin(1,0,1,0.2, "cm"),
        panel.grid.major = element_line(size = 0.8, colour = 'black', linetype = 'dashed'),
        panel.grid.minor = element_blank(),
        panel.ontop = T) +
  scale_fill_gradientn(colours = c("#A5B6C5", 'white',"#F0988C"), 
                       values = rescale(c(-10,0,45)),
                       breaks = c(-10, 0, 10, 20, 30, 40),
                       labels = c("-10", "0", "10", "20", "30", "40"),
                       name = "Synergy score") + 
  ylab(expression(paste('NB-598 (',mu,"M)", sep = ""))) + 
  xlab("Neratinib (nM)\n") +
  guides(fill=guide_colorbar("Synergy score",
                             direction = "vertical", title.position = "left",
                             label.position="right", label.hjust = 0.5, label.vjust = 0.5,
                             label.theme = element_text(angle = 90, size = 16),
                             frame.colour = "black", ticks.colour = "black")) + 
  coord_flip() 
p1

### add geom_text
seg <- data.frame(x = 0.75, xend = 1.25, y = 1, yend = 1)
row1 <- ggplot(seg) + theme_void() + geom_richtext(aes(x = 1 ,y = 1, label = p[["1"]][["labels"]][["subtitle"]]), 
                                                   fill = "white", label.size = NA, size = 6, angle = 90, hjust = 0.35)

p2 <- ggarrange(p1, row1, ncol = 2, nrow = 1, widths = c(11, 1), heights = c(12, 12))
p2

png(
  filename  = 'paper_figure/Fig_6e-1.png',
  width     = 5,
  height    = 5.5,
  unit = 'in',
  res = 300
)
p2
dev.off()

### MDAMB231 Atorvastatin + Neratinib
p <- PlotSynergy(MDA_statin_score, type = "2D", grid = F)
p1 <- p$`1` + 
  theme(text = element_text(family = 'Sans'), 
        axis.text = element_text(size = 20, color = 'black'), 
        axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(size = 22, angle = 180, margin = margin(t = 20, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
        legend.title = element_text(size = 18, hjust = 0.5, angle = 90),
        legend.position = 'right',
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(0.4,'cm'),
        legend.key = element_rect(color="black"),
        plot.title = element_blank(), 
        plot.subtitle = element_blank(),
        plot.margin = margin(1,0,1,0.2, "cm"),
        panel.grid.major = element_line(size = 0.8, colour = 'black', linetype = 'dashed'),
        panel.grid.minor = element_blank(),
        panel.ontop = T) +
  scale_fill_gradient2(low = '#A5B6C5', mid = 'white', high = '#F0988C', 
                       limits = c(-20, 20),
                       breaks = c(-20, -10, 0, 10, 20),
                       labels = c("-20", "-10", "0", "10", "20"),
                       name = "Synergy\nScore") +
  ylab(expression(paste("Atorvastatin (",mu,"M)", sep = ""))) +
  guides(fill=guide_colorbar("Synergy score",
                             direction = "vertical", title.position = "left",
                             label.position="right", label.hjust = 0.5, label.vjust = 0.5,
                             label.theme = element_text(angle = 90, size = 16),
                             frame.colour = "black", ticks.colour = "black")) + 
  coord_flip() 
p1

### add geom_text
seg <- data.frame(x = 0.75, xend = 1.25, y = 1, yend = 1)
row1 <- ggplot(seg) + theme_void() + geom_richtext(aes(x = 1 ,y = 1, label = p[["1"]][["labels"]][["subtitle"]]), 
                                                   fill = "white", label.size = NA, size = 6, angle = 90, hjust = 0.35)

p2 <- ggarrange(p1, row1, ncol = 2, nrow = 1, widths = c(11, 1), heights = c(12, 12))
p2

png(
  filename  = 'paper_figure/Fig_6e-2.png',
  width     = 5,
  height    = 5.5,
  unit = 'in',
  res = 300
)
p2
dev.off()

### T47D NB598 + Neratinib
p <- PlotSynergy(T47D_NB_score, type = "2D", grid = F)
p
p1 <- p$`1` + 
  theme(text = element_text(family = 'Sans'), 
        axis.text = element_text(size = 20, color = 'black'), 
        axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(size = 22, angle = 180, margin = margin(t = 20, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
        legend.title = element_text(size = 18, hjust = 0.5, angle = 90),
        legend.position = 'right',
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(0.4,'cm'),
        legend.key = element_rect(color="black"),
        plot.title = element_blank(), 
        plot.subtitle = element_blank(),
        plot.margin = margin(1,0,1,0.2, "cm"),
        panel.grid.major = element_line(size = 0.8, colour = 'black', linetype = 'dashed'),
        panel.grid.minor = element_blank(),
        panel.ontop = T) +
  scale_fill_gradient2(low = '#A5B6C5', mid = 'white', high = '#F0988C', 
                       limits = c(-20, 20),
                       breaks = c(-20, -10, 0, 10, 20),
                       labels = c("-20", "-10", "0", "10", "20"),
                       name = "Synergy\nScore") +
  ylab(expression(paste("NB-598 (",mu,"M)", sep = ""))) +
  guides(fill=guide_colorbar("Synergy score",
                             direction = "vertical", title.position = "left",
                             label.position="right", label.hjust = 0.5, label.vjust = 0.5,
                             label.theme = element_text(angle = 90, size = 16),
                             frame.colour = "black", ticks.colour = "black")) + 
  coord_flip() 
p1

### add geom_text
seg <- data.frame(x = 0.75, xend = 1.25, y = 1, yend = 1)
row1 <- ggplot(seg) + theme_void() + geom_richtext(aes(x = 1 ,y = 1, label = p[["1"]][["labels"]][["subtitle"]]), 
                                                   fill = "white", label.size = NA, size = 6, angle = 90, hjust = 0.35)

p2 <- ggarrange(p1, row1, ncol = 2, nrow = 1, widths = c(11, 1), heights = c(12, 12))
p2

png(
  filename  = 'paper_figure/Fig_6d-2.png',
  width     = 5,
  height    = 5.5,
  unit = 'in',
  res = 300
)
p2
dev.off()

### T47D Atorvastatin + Neratinib
p <- PlotSynergy(T47D_statin_score, type = "2D", grid = F)
p1 <- p$`1` + 
  theme(text = element_text(family = 'Sans'), 
        axis.text = element_text(size = 20, color = 'black'), 
        axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(size = 22, angle = 180, margin = margin(t = 20, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
        legend.title = element_text(size = 18, hjust = 0.5, angle = 90),
        legend.position = 'right',
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(0.4,'cm'),
        legend.key = element_rect(color="black"),
        plot.title = element_blank(), 
        plot.subtitle = element_blank(),
        plot.margin = margin(1,0,1,0.2, "cm"),
        panel.grid.major = element_line(size = 0.8, colour = 'black', linetype = 'dashed'),
        panel.grid.minor = element_blank(),
        panel.ontop = T) +
  scale_fill_gradient2(low = '#A5B6C5', mid = 'white', high = '#F0988C', 
                       limits = c(-20, 20),
                       breaks = c(-20, -10, 0, 10, 20),
                       labels = c("-20", "-10", "0", "10", "20"),
                       name = "Synergy\nScore") +
  ylab(expression(paste("Atorvastatin (",mu,"M)", sep = ""))) +
  guides(fill=guide_colorbar("Synergy score",
                             direction = "vertical", title.position = "left",
                             label.position="right", label.hjust = 0.5, label.vjust = 0.5,
                             label.theme = element_text(angle = 90, size = 16),
                             frame.colour = "black", ticks.colour = "black")) + 
  coord_flip() 
p1

### add geom_text
seg <- data.frame(x = 0.75, xend = 1.25, y = 1, yend = 1)
row1 <- ggplot(seg) + theme_void() + geom_richtext(aes(x = 1 ,y = 1, label = p[["1"]][["labels"]][["subtitle"]]), 
                                                   fill = "white", label.size = NA, size = 6, angle = 90, hjust = 0.35)

p2 <- ggarrange(p1, row1, ncol = 2, nrow = 1, widths = c(11, 1), heights = c(12, 12))
p2

png(
  filename  = 'paper_figure/Fig_6d-1.png',
  width     = 5,
  height    = 5.5,
  unit = 'in',
  res = 300
)
p2
dev.off()
