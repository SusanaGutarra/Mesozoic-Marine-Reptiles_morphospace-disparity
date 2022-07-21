
#LOAD PACKAGES

library(ggplot2)
library(factoextra) # R package to help in the interpretation of PCA (ggplot2-based visualization)
library(ggExtra)
library(dplyr)


#################
#################



rm(list = ls())

pca_aq_reptiles<-read.csv("FossilAqreptiles_Output.csv", 
                          row.names = 1,
                          header=T)


pca_aq_reptiles


#PCA 1-2
#########

# Calculate the hulls for each group
hull_cyl <-pca_aq_reptiles %>%
  group_by(Group_1) %>%
  slice(chull(Dim.1,Dim.2))


#Make plot with dots and polygons
pca_plot <- ggplot(pca_aq_reptiles,aes(x = Dim.1, y = Dim.2, fill = Group_1, colour=Group_1, shape=Group_1
                                       )) + #I need to specify color here so that ggMarginal can take it
  #geom_point(size=2) +
  geom_polygon(data = hull_cyl, alpha = 0.2, linetype=0)+      #0 is the size of the polygon line (before I had size=0)
  scale_colour_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","blue"))+
  scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","blue"))+
  scale_shape_manual(values= c(24,23,21,25,21,24,21,25,22))+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #horizontal reference line in O
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #vertical reference line in O
  labs(x = "PC1 (32.05%)",
       y = "PC2 (24.34%)")+
  theme_light()

pca_plot


pca_plot <- pca_plot + geom_point(color="black",stroke=0.3, size=2.2) #add points in layers! otherwise it messes up the colors of histograms, stroke indicates thickness of line

pca_plot



#Format plot:
pca_plot <- pca_plot  + theme (text=element_text(),
                               axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                               axis.text.y = element_text(size=12,color="#303030"),
                               axis.title.x  = element_text(size=14,color="#303030", margin = margin(t = 10, r = 0, b = 0, l = 0)), #margin is added to distance title from axis
                               axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                               strip.text = element_text(size=12,color="#303030"), #Change size of facet title
                               panel.grid.major.x = element_blank(),
                               panel.grid.minor.x = element_blank(),
                               panel.grid.major.y = element_blank(),
                               panel.grid.minor.y = element_blank(),
                               axis.ticks = element_line(colour = "#303030"),
                               panel.border = element_rect(color="#303030", fill= NA, size=0.6),
                               aspect.ratio = 1,
                               legend.position = "none")

pca_plot

#I will reduce to 1/2 size, so I use font x2
ggsave("PCA 1-2.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 13, height = 13, units = "cm",
       dpi = 300, limitsize = TRUE) 



 #Add marginal plots with ggMarginal (ggplot2-based package)
#https://cran.r-project.org/web/packages/ggExtra/vignettes/ggExtra.html
#https://cran.r-project.org/web/packages/ggExtra/readme/README.html
#https://www.rdocumentation.org/packages/ggExtra/versions/0.9/topics/ggMarginal
#https://www.trifields.jp/how-to-plot-the-scatter-plot-and-marginal-distribution-using-ggplot2-in-r-2992


#This function uses the plot as an argument (not the dataframe)

ggMarginal(pca_plot, colour =  "black", groupFill =  TRUE, #The aesthetics feed from the ones in the plot, could also take color of line: groupColor = TRUE
           type = "violin",               #here we can also use boxplot, histogram, density...
           size = 2,                      #control size of marginal plots compared to the main plot   
           xparams = list(size = 0.1),    #adjust parameters of the x - marginal plot. e.g: size controls line thickness
           yparams = list(size = 0.1))    #adjust parameters of the y - marginal plot  e.g: size controls line thickness                               



#################
#################

#PCA 3-2 (FIGURE supplementary)
#-----------------------------# 

# Calculate the hulls for each group
hull_cyl <-pca_aq_reptiles %>%
  group_by(Group_1) %>%
  slice(chull(Dim.3,Dim.2))


#Make plot with dots and polygons
pca_plot <- ggplot(pca_aq_reptiles,aes(x = Dim.3, y = Dim.2, fill = Group_1, colour=Group_1, shape=Group_1)) + #I need to specify color here so that ggMarginal can take it
  geom_point(size=2) +
  geom_polygon(data = hull_cyl, alpha = 0.2, size=0)+      #0 is the size of the polygon line
  scale_colour_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","blue"))+
  scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","blue"))+
  scale_shape_manual(values= c(24,23,21,25,21,24,21,25,22))+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #horizontal reference line in O
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #vertical reference line in O
  labs(x = "PC3 (10.93%)",
       y = "PC2 (24.34%)")+
  theme_light()

pca_plot

pca_plot <- pca_plot + geom_point(color="black", stroke=0.3, size=2.2) #add the black rim of the points in layers! otherwise it messes up the colors of histograms

pca_plot


#Format plot:
pca_plot <- pca_plot  + theme (text=element_text(),
                               axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                               axis.text.y = element_text(size=12,color="#303030"),
                               axis.title.x  = element_text(size=14,color="#303030", margin = margin(t = 10, r = 0, b = 0, l = 0)), #margin is added to distance title from axis
                               axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                               strip.text = element_text(size=12,color="#303030"), #Change size of facet title
                               panel.grid.major.x = element_blank(),
                               panel.grid.minor.x = element_blank(),
                               panel.grid.major.y = element_blank(),
                               panel.grid.minor.y = element_blank(),
                               axis.ticks = element_line(colour = "#303030"),
                               panel.border = element_rect(color="#303030", fill= NA, size=0.6),
                               aspect.ratio = 1,
                               legend.position = "none")

pca_plot


#I will reduce to 1/2 size, so I use font x2
ggsave("PCA 3-2.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 13, height = 13, units = "cm",
       dpi = 300, limitsize = TRUE) 



#Add marginal plots with ggMarginal (ggplot2-based package)
#---------------------------------------------------------#

#https://cran.r-project.org/web/packages/ggExtra/vignettes/ggExtra.html
#https://cran.r-project.org/web/packages/ggExtra/readme/README.html


#This function uses the plot as an argument (not the dataframe)

ggMarginal(pca_plot, colour =  "black", groupFill =  TRUE, #The aesthetics feed from the ones in the plot, could also take color of line: groupColor = TRUE
           type = "violin",               #here we can also use boxplot, histogram, density...
           size = 2,                      #control size of marginal plots compared to the main plot   
           xparams = list(size = 0.1),    #adjust parameters of the x - marginal plot. e.g: size controls line thickness
           yparams = list(size = 0.1))    #adjust parameters of the y - marginal plot  e.g: size controls line thickness                               




####################
####################


##PCA1-2 WITH NUMBERED LABELS (FIGURE)
#------------------------------------#

#PCA 1-2
#########

dframe1<-read.csv("FossilAqreptiles_Output.csv", 
                          row.names = 1,
                          header=T)

dframe1


Group <- factor(dframe1$Group_1, 
                      levels=c("Ichthyosauromorpha", "Mosasauroidea", "Pantestudines", "Rhynchocephalia", "Sauropterygia","Saurosphargidae","Thalattosauria", "Thalattosuchia", "Tanystropheidae"
                                ),
                      labels = c("Ichthyosauromorpha", "Mosasauroidea", "Pantestudines", "Rhynchocephalia", "Sauropterygia","Saurosphargidae","Thalattosauria", "Thalattosuchia", "Tanystropheidae"
                                ))

Group


pca_plot <- ggplot(dframe1, aes(x = Dim.1, y = Dim.2, fill = Group,shape=Group,color=Group))+
  geom_point(size=2, colour="black", stroke=0.4) +
  scale_colour_manual(values= c("#FF2626","#FF8300","#0083ab","#c061ff","#19D27E","darkgrey","#FFD300","#3944BC", "#add8f0"))+
  scale_fill_manual(values= c("#FF2626","#FF8300","#0083ab","#c061ff","#19D27E","darkgrey","#FFD300","#3944BC", "#add8f0")) +
  scale_shape_manual(values= c(24,23,21,25,21,24,21,25,22))+
  scale_y_continuous (labels = scales::number_format(accuracy = 0.1))+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #horizontal reference line in O
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #vertical reference line in O
  labs(x = "PC1 (32.05%)",
       y = "PC2 (24.34%)")+
  geom_text(aes(label=Num), size=3.5, hjust=-0.4)+
  #ggrepel::geom_text_repel(aes(label=Num, color=Group), 
  #                         force=0.1,
  #                         segment.size=0.25,
  #                         segment.alpha=0.4,
  #                         size=2.8, hjust=0) +
  #geom_text(aes(label=Name, x=Inf,hjust=0)) +
  theme_light()

pca_plot


#Format plot:
pca_plot  + theme (text=element_text(),
                   axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                   axis.text.y = element_text(size=12,color="#303030"),
                   axis.title.x  = element_text(size=14,color="#303030", margin = margin(t = 10, r = 0, b = 0, l = 0)), #margin is added to distance title from axis
                   axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 5, b = 0, l = 0)), #margin is added to distance title from axis
                   strip.text = element_text(size=12,color="#303030"), #Change size of facet title
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.ticks = element_line(colour = "#303030"),
                   panel.border = element_rect(color="#303030", fill= NA, size=0.6),
                   aspect.ratio = 1,
                   legend.position = "none")


#will reduce to half size, so I use double size of font
ggsave("PCA 1-2-labels_2.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 16, height = 16, units = "cm",
       dpi = 300, limitsize = TRUE) 



#Zoom of the right-inferior quadrant
pca_plot <- ggplot(dframe1, aes(x = Dim.1, y = Dim.2,fill = Group,shape=Group,color=Group))+
  #geom_point(aes(fill = Locomotion1, colour=Locomotion1, shape=Locomotion1), size=1, stroke=0.4)+
  xlim(-1.9,2)+
  ylim(-3,-0.4)+
  geom_point(size=2.8, colour="black") +
  scale_colour_manual(values= c("#FF2626","#FF8300","#0083ab","#c061ff","#19D27E","darkgrey","#FFD300","#3944BC", "#add8f0"))+
  scale_shape_manual(values= c(24,23,21,25,21,24,21,25,22))+
  scale_fill_manual(values= c("#FF2626","#FF8300","#0083ab","#c061ff","#19D27E","darkgrey","#FFD300","#3944BC", "#add8f0")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #horizontal reference line in O
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #vertical reference line in O
  labs(x = "PCA 1 (32.05%)",
       y = "PCA 2 (24.34%)")+
  #geom_text(aes(label=Num), size=5, hjust=-0.45, vjust=0.45)+
  ggrepel::geom_text_repel(aes(label=Num,color=Group), 
                           force=0.1,
                           segment.size=0.4,
                           segment.alpha=0.4,
                           size=3.5) +
  #geom_text(aes(label=Name, x=Inf,hjust=0)) +
  theme_light()

pca_plot
  
  

#Format plot:
pca_plot  + theme (text=element_text(),
                   axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                   axis.text.y = element_text(size=12,color="#303030"),
                   axis.title.x  = element_blank(), #margin is added to distance title from axis
                   axis.title.y  = element_blank(), #margin is added to distance title from axis
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.ticks = element_line(colour = "#303030"),
                   panel.border = element_rect(color="#303030", fill= NA, size=0.6),
                   aspect.ratio = 1,
                   legend.position = "none")


#will reduce to half size, so I use double size of font
ggsave("PCA 1-2-labels_2_zoom.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, height = 13, width=13 , units = "cm",
       dpi = 300, limitsize = TRUE) 



#GEOM_DENSITY_RIDGES
#-------------------#


library(ggridges)
#https://stackoverflow.com/questions/18165578/subset-and-ggplot2
#https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html


## Make subsets first, so that we can exclude Tanystropheids and Saurosphargidae

#example:
#ggplot(subset(dat,ID %in% c("P1" , "P3"))) + 
#  geom_line(aes(Value1, Value2, group=ID, colour=ID))


#Create a subset of the data only with the group categories that I want to plot
pca_aq_reptiles2<- subset(pca_aq_reptiles,Group_1 %in% c("Thalattosuchia", "Thalattosauria", "Sauropterygia", "Rhynchocephalia","Pantestudines","Mosasauroidea","Ichthyosauromorpha"))
pca_aq_reptiles2

#Making a factor and indicating the right order (levels)
Group <- factor(pca_aq_reptiles2$Group_1, 
                  levels = c("Thalattosuchia", "Thalattosauria", "Sauropterygia", "Rhynchocephalia","Pantestudines","Mosasauroidea","Ichthyosauromorpha"))
Group


#PCA1-Ridges
############


#PCA1-Groups smooth histogram-1
plot <- ggplot(pca_aq_reptiles2, aes(x=Dim.1,fill = Group, y= Group)) + 
  geom_density_ridges(alpha=0.5,scale = 0.95, size=0.1,rel_min_height = 0.01) +  
  scale_fill_manual(values= c("#3944bc","#FFD300", "#19D27E","#c061ff","#009ccc","#FF8300","#FF2626"))+ 
  labs(x = "PCA 1 (32.05%)")+
  theme_light()

plot


#PCA1-Groups smooth histogram-1 + plot points
plot <- ggplot(pca_aq_reptiles2, aes(x=Dim.1,fill = Group, y = Group)) + 
  geom_density_ridges(alpha=0.5, scale = 0.95, size=0.1,rel_min_height = 0.01,
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1) +  
  scale_fill_manual(values= c("#3944bc","#FFD300", "#19D27E","#c061ff","#009ccc","#FF8300","#FF2626"))+
  labs(x = "PCA 1 (32.05%)")+
  theme_light()

plot



#PCA1-Groups smooth histogram-2 (stat density)
plot <- ggplot(pca_aq_reptiles2, aes(x=Dim.1,fill = Group, y = Group,height = stat(density))) + 
  geom_density_ridges(alpha=0.5,scale = 0.95, size=0.1,rel_min_height = 0.01,stat = "density") +  
  scale_fill_manual(values= c("#3944bc","#FFD300", "#19D27E","#c061ff","#009ccc","#FF8300","#FF2626"))+
  labs(x = "PCA 1 (32.05%)")+
  theme_light()

plot


#PCA1-Groups bar histogram
plot <- ggplot(pca_aq_reptiles2, aes(x=Dim.1,fill = Group, y = Group,height = stat(density))) + 
  geom_density_ridges(alpha=0.5,scale = 0.95, size=0.1,stat = "binline", bins=20, draw_baseline = FALSE) +
  scale_fill_manual(values= c("#3944bc","#FFD300", "#19D27E","#c061ff","#009ccc","#FF8300","#FF2626"))+
  labs(x = "PCA 1 (32.05%)")+
  theme_light()

plot

#PCA1-Groups bar histogram - grey (FIGURE)
plot <- ggplot(pca_aq_reptiles2, aes(x=Dim.1, y = Group, color=Group, height = stat(density))) + 
  geom_density_ridges(fill= "grey", 
                      alpha=0.5,scale = 0.9, 
                      size=0.2,stat = "binline", 
                      bins=20, 
                      draw_baseline = FALSE) +
  scale_color_manual(values= c("#3944bc","#FFD300", "#19D27E","#c061ff","#009ccc","#FF8300","#FF2626"))+
  labs(x = "PCA 1 (32.05%)")+
  theme_light()

plot



#Format WITH labels
plot +  theme (text=element_text(),
         axis.text.x = element_text(size=6, color="#303030", hjust=0.5),
         axis.text.y = element_text(size=6,color="#303030"),
         axis.title.x  = element_text(size=7,color="#303030", margin = margin(t=10, r = 0, b = 0, l = 0)), #margin is added to distance title from axis
         axis.title.y  = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.y = element_blank(),
         panel.grid.major.y = element_line(colour = "lightgrey"),
         axis.ticks = element_line(colour = "#303030"),
         panel.border = element_rect(color="#303030", fill= NA, size=0.3),
         aspect.ratio = 0.3,
         legend.position = "none")



#Format WITHOUT labels
plot +   theme (text=element_text(),
         axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y  = element_blank(),
         axis.title.x  = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.y = element_blank(),
         panel.grid.major.y = element_line(colour = "lightgrey"),
         axis.ticks =element_blank(),
         panel.border = element_blank(), #element_rect(color="#303030", fill= NA, size=0.3),
         aspect.ratio = 0.3,
         legend.position = "none")



ggsave("PCA 1-Groups histograms.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 6.5, height = 2.5, units = "cm",
       dpi = 300, limitsize = TRUE) 




#PCA2-Ridges
############


#PCA2-Groups bar histogram - grey
plot <- ggplot(pca_aq_reptiles2, aes(x=Dim.2, y = Group_1,height = stat(density))) + 
  geom_density_ridges(fill= "grey", alpha=0.5,scale = 0.95, size=0.1,stat = "binline", bins=20, draw_baseline = FALSE) +
  #scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","#FFD300","blue"))+
  labs(x = "PCA 2 (23.96%)")+
  theme_light()

plot


#PCA1-Groups bar histogram - grey (FIGURE)
plot <- ggplot(pca_aq_reptiles2, aes(x=Dim.2, y = Group, color=Group, height = stat(density))) + 
  geom_density_ridges(fill= "grey", 
                      alpha=0.5,scale = 0.9, 
                      size=0.2,stat = "binline", 
                      bins=20, 
                      draw_baseline = FALSE) +
  scale_color_manual(values= c("#3944bc","#FFD300", "#19D27E","#c061ff","#009ccc","#FF8300","#FF2626"))+
  labs(x = "PCA 2 (24.39%)")+
  theme_light()

plot


#Format WITHOUT labels
plot +coord_flip() + #It rotates the axis
     theme (text=element_text(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y  = element_blank(),
            axis.title.x  = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(colour = "lightgrey"),
            axis.ticks =element_blank(),
            panel.border = element_blank(),
            aspect.ratio = 3.33,
            legend.position = "none")


ggsave("PCA 2-Groups histograms.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 2.3, height = 6.5, units = "cm",
       dpi = 300, limitsize = TRUE) 



#PCA3-Ridges
############


#PCA3-Groups bar histogram - grey
plot <- ggplot(pca_aq_reptiles2, aes(x=Dim.3, y = Group_1,height = stat(density))) + 
  geom_density_ridges(fill= "grey", alpha=0.5,scale = 0.95, size=0.1,stat = "binline", bins=20, draw_baseline = FALSE) +
  #scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","#FFD300","blue"))+
  labs(x = "PCA 3 (10.93%)")+
  theme_light()

plot


#PCA3-Groups bar histogram - grey (FIGURE)
plot <- ggplot(pca_aq_reptiles2, aes(x=Dim.3, y = Group, color=Group, height = stat(density))) + 
  geom_density_ridges(fill= "grey", 
                      alpha=0.5,scale = 0.9, 
                      size=0.2,stat = "binline", 
                      bins=20, 
                      draw_baseline = FALSE) +
  scale_color_manual(values= c("#3944bc","#FFD300", "#19D27E","#c061ff","#009ccc","#FF8300","#FF2626"))+
  labs(x = "PCA 3 (10.93%)")+
  theme_light()

plot


#Format WITHOUT labels
plot +   theme (text=element_text(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.title.y  = element_blank(),
                axis.title.x  = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.major.y = element_line(colour = "lightgrey"),
                axis.ticks =element_blank(),
                panel.border = element_blank(), #element_rect(color="#303030", fill= NA, size=0.3),
                aspect.ratio = 0.3,
                legend.position = "none")

ggsave("PCA 3-Groups histograms.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 6.5, height = 2.5, units = "cm",
       dpi = 300, limitsize = TRUE) 



#################
#################


##MAKE NICE HEX-BIN and DENSITY PLOTS
#-----------------------------------#

#http://www.sthda.com/english/articles/32-r-graphics-essentials/131-plot-two-continuous-variables-scatter-graph-and-alternatives/

#Hexagonal bins PC1-2
pca_plot2 <-ggplot(pca_aq_reptiles, aes(x = Dim.1, y = Dim.2)) +
  geom_hex(bins = 10, color = "white")+    #the size of the bin is given by bin=n
  scale_fill_gradient(low =  "grey", high = "#FC4E07")+   # low was "#00AFBB"
  labs(x = "PC1 (32.05%)",
       y = "PC2 (24.34%)")+
  theme_light()

pca_plot2

pca_plot2  + theme (text=element_text(),
                               axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                               axis.text.y = element_text(size=12,color="#303030"),
                               axis.title.x  = element_text(size=14,color="#303030", margin = margin(t = 10, r = 0, b = 0, l = 0)), #margin is added to distance title from axis
                               axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                               strip.text = element_text(size=12,color="#303030"), #Change size of facet title
                               panel.grid.major.x = element_blank(),
                               panel.grid.minor.x = element_blank(),
                               panel.grid.major.y = element_blank(),
                               panel.grid.minor.y = element_blank(),
                               axis.ticks = element_line(colour = "#303030"),
                               panel.border = element_rect(color="#303030", fill= NA, size=0.5),
                               aspect.ratio = 1,
                               legend.position = "none")

#save 7x8 for presentation

# Will plot half the size so, saving with 2X size of font and dimensions
ggsave("PCA 1-2_hexbin-ggplot.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 13, height = 13, units = "cm",
       dpi = 300, limitsize = TRUE) 

#ggsave("LEGEND-hexbin.pdf", plot = last_plot(), device = "pdf", path = NULL,
#       scale = 1, width = 15, height = 15, units = "cm",
#       dpi = 300, limitsize = TRUE) 



#Hexagonal bins PC3-2
pca_plot2 <-ggplot(pca_aq_reptiles, aes(x = Dim.3, y = Dim.2)) +
  geom_hex(bins = 10, color = "white")+    #the size of the bin is given by bin=n
  scale_fill_gradient(low =  "grey", high = "#FC4E07")+   # low was "#00AFBB"
  labs(x = "PC3 (10.93%)",
       y = "PC2 (24.34%)")+
  theme_light()

pca_plot2

pca_plot2  + theme (text=element_text(),
                    axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                    axis.text.y = element_text(size=12,color="#303030"),
                    axis.title.x  = element_text(size=14,color="#303030", margin = margin(t = 10, r = 0, b = 0, l = 0)), #margin is added to distance title from axis
                    axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                    strip.text = element_text(size=12,color="#303030"), #Change size of facet title
                    panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    panel.grid.major.y = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    axis.ticks = element_line(colour = "#303030"),
                    panel.border = element_rect(color="#303030", fill= NA, size=0.5),
                    aspect.ratio = 1,
                    legend.position = "none")


# Will plot half the size so, saving with 2X size of font and dimensions
ggsave("PCA 3-2_hexbin-ggplot.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 13, height = 13, units = "cm",
       dpi = 300, limitsize = TRUE) 

#ggsave("LEGEND-hexbin3-2.pdf", plot = last_plot(), device = "pdf", path = NULL,
#       scale = 1, width = 15, height = 15, units = "cm",
#       dpi = 300, limitsize = TRUE) 


# 2d density estimation

pca_plot3 <- ggplot(pca_aq_reptiles, aes(x = Dim.1, y = Dim.2)) +
  scale_x_continuous(limits=c(-4.8,6.5))+
  scale_y_continuous(limits=c(-4.1,5))+
  geom_point(color = "lightgray")+
  geom_density_2d(color="blue",size=0.5)+
  labs(x = "PCA 1 (32.05%)",
       y = "PCA 2 (24.34%)")+
  theme_light()

pca_plot3 

pca_plot3  + theme (text=element_text(),
                   axis.text.x = element_text(size=10, color="#303030", hjust=0.5),
                   axis.text.y = element_text(size=10,color="#303030"),
                   axis.title.x  = element_text(size=12,color="#303030", margin = margin(t = 20, r = 0, b = 0, l = 0)), #margin is added to distance title from axis
                   axis.title.y  = element_text(size=12,color="#303030", margin = margin(t = 0, r = 20, b = 0, l = 0)), #margin is added to distance title from axis
                   strip.text = element_text(size=12,color="#303030"), #Change size of facet title
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.ticks = element_line(colour = "#303030"),
                   panel.border = element_rect(color="#303030", fill= NA, size=1.2),
                   aspect.ratio = 1)

#legend.position = "none")
#save 7x8

ggsave("PCA 1-2 Density lines.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 15, height = 15, units = "cm",
       dpi = 300, limitsize = TRUE) 


# Use different geometry and change the gradient color
pca_plot3 + stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_gradient(low =  "grey", high = "#FC4E07")
  

# Same but plotting it from scratch

pca_plot3 <- ggplot(pca_aq_reptiles, aes(x = Dim.1, y = Dim.2)) +
  scale_x_continuous(limits=c(-4.8,6.5))+
  scale_y_continuous(limits=c(-4.1,5))+
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_gradient(low =  "grey", high = "#FC4E07")+
  #geom_density_2d(color="blue",size=0.5)+
  labs(x = "PCA 1 (30.89%)",
       y = "PCA 2 (26.64%)")+
  theme_light()

pca_plot3


#####################
#####################


