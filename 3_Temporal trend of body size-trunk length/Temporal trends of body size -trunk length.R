#LOAD PACKAGES
library(FactoMineR) #R package for computing PCA
#https://cran.r-project.org/web/packages/FactoMineR/FactoMineR.pdf
library(ape)
library(permute)
library(lattice)
library(vegan)
library(dispRity)
library(ggplot2)
library(factoextra) # R package to help in the interpretation of PCA (ggplot2-based visualization)
library(dendextend)
library(maps)
library(phytools)
#install.packages("summarize", repos="http://R-Forge.R-project.org")
library(summarize)

#install.packages("RVAideMemoire") #Pairwise F test to compare variances
library(RVAideMemoire)



rm(list = ls())


# OPEN FILE

dframe1<-read.csv("MMreptiles_Trunk.csv",   #MMreptiles_Trunk.csv with 8 time bins
                  row.names = 1,
                  header=T)

dframe1



#1. TRUNK SIZE (ALL SIZES THROGH TIME)
#------------------------------------#
#Before plotting FILTER THE DATA THAT HAVE <NA> in a given column
dframe2  <- dframe1[!is.na(dframe1$Time), ]
dframe2


pca_plot <-ggplot(dframe1, aes(x=Age, y=TRUNK, fill = Group, colour=Group, shape=Group)) + 
  geom_point(size=2.8,stroke =0.1) +
  scale_x_reverse (limits=c(255,66))+
  #scale_y_log10()+
  scale_colour_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","#3944bc"))+
  scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","#3944bc"))+
  scale_shape_manual(values= c(24,23,21,25,21,24,21,25,22))+
  #geom_vline(xintercept = 247.2, linetype = "dashed", size = 0.3, col="gray80")+  #vertical reference line
  #geom_vline(xintercept = 201.3,  size = 0.5, col="gray80")+  #vertical reference line
  #geom_vline(xintercept = 237, linetype = "dashed", size = 0.3, col="gray80")+  #vertical reference line
  #geom_vline(xintercept = 174.1, linetype = "dashed", size = 0.3, col="gray80")+  #vertical reference line
  #geom_vline(xintercept = 163.5, linetype = "dashed", size = 0.3, col="gray80")+  #vertical reference line
  #geom_vline(xintercept = 145, size = 0.5, col="gray80")+  #vertical reference line
  #geom_vline(xintercept = 100.5, linetype = "dashed", size = 0.3, col="gray80")+  #vertical reference line
  labs(x = "Time (Ma)",
       y = "Trunk length (cm)")+
  theme_light()

pca_plot

pca_plot <- pca_plot + geom_point(color="black", size=2.8,stroke =0.4) #add the black rim of the points in layers

pca_plot


#Format plot:
pca_plot  + theme (text=element_text(),
                   axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                   axis.text.y = element_text(size=12,color="#303030"),
                   axis.title.x  = element_blank(), #margin is added to distance title from axis
                   axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.ticks = element_line(colour = "#303030"),
                   panel.border = element_rect(color="#303030", fill= NA, size=0.6),
                   aspect.ratio = 0.5, #was 0.5
                   legend.position = "none")


 ggsave("Z_Trunk size-1.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 16, height = 9, units = "cm",
       dpi = 300, limitsize = TRUE) 


#will reduce to half size, so I use double size of font
#ggsave("Z_Trunk size-1-log10.pdf", plot = last_plot(), device = "pdf", path = NULL,
#       scale = 1, width = 14, height = 7, units = "cm",
#       dpi = 300, limitsize = TRUE) 



#1.2_TRUNK SIZE ALL TAXA 11 TIME BINS - BOXPLOT 
#--------------------------------------------#

#Before plotting FILTER THE DATA THAT HAVE <NA> in a given column
dframe2  <- dframe1[!is.na(dframe1$Time_4), ]
dframe2


pca_plot <-ggplot(dframe2, aes(x=Time4, y=TRUNK, group=Time_4)) + 
  geom_boxplot(aes(group = cut_width(Time_4, 0.25)),
                   fill = "lightgrey", 
                   lwd=0.4,
                   outlier.alpha = 0.2,
                   alpha=0.4) + 
  geom_jitter(width = 0.2, alpha=0.25)+
  scale_x_reverse (limits=c(255, 66))+
  labs(x = "Time (Ma)",
       y = "Trunk length (cm)")+
  theme_light()

pca_plot


#Format plot:
pca_plot  + theme (text=element_text(),
                   axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                   axis.text.y = element_text(size=12,color="#303030"),
                   axis.title.x  = element_blank(),
                   axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                   strip.text = element_text(size=12,color="#303030"), #Change size of facet title
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.ticks = element_line(colour = "#303030"),
                   panel.border = element_rect(color="#303030", fill= NA, size=0.6),
                   aspect.ratio = 0.5,
                   legend.position = "none")


ggsave("Z_Trunk size 11 time bins.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 16, height = 9, units = "cm",
       dpi = 300, limitsize = TRUE) 

#ggsave("Z_Trunk size-log10.pdf", plot = last_plot(), device = "pdf", path = NULL,
#       scale = 1, width = 14, height = 7, units = "cm",
#       dpi = 300, limitsize = TRUE) 


#1.2_TRUNK SIZE ALL TAXA BROAD TIME BINS - BOXPLOT 
#--------------------------------------------#

#Before plotting FILTER THE DATA THAT HAVE <NA> in a given column
dframe2  <- dframe1[!is.na(dframe1$Time), ]
dframe2

library(forcats)
library(dplyr)

pca_plot <- dframe2 %>%
  mutate(Time = fct_relevel(Time, "Triassic", "Jurassic", "Cretaceous")) %>%  # Reorder following a precise order
  ggplot( aes(x=Time, y=TRUNK)) +
  geom_boxplot(fill = "lightgrey", 
               lwd=0.4,
               outlier.alpha = 0.2,
               alpha=0.4) + 
  geom_jitter(width = 0.2, alpha=0.25)+
  labs(y = "Trunk length (cm)")+
  theme_light()

pca_plot

#Format plot:
pca_plot  + theme (text=element_text(),
                   axis.text.x = element_text(size=12, color="#303030", hjust=1, angle=45),
                   axis.text.y = element_text(size=12,color="#303030"),
                   axis.title.x  = element_blank(),
                   axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                   strip.text = element_text(size=12,color="#303030"), #Change size of facet title
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.ticks = element_line(colour = "#303030"),
                   panel.border = element_rect(color="#303030", fill= NA, size=0.6),
                   aspect.ratio = 1.4,
                   legend.position = "none")


ggsave("Z_Trunk size BROAD time bins.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, height = 9, width = 8.8, units = "cm",
       dpi = 300, limitsize = TRUE) 

#ggsave("Z_Trunk size-log10.pdf", plot = last_plot(), device = "pdf", path = NULL,
#       scale = 1, width = 14, height = 7, units = "cm",
#       dpi = 300, limitsize = TRUE) 




#1.3_TRUNK SIZE ALL TAXA 13 TIME BINS - BOXPLOT 
#--------------------------------------------#

#Before plotting FILTER THE DATA THAT HAVE <NA> in a given column
dframe2  <- dframe1[!is.na(dframe1$Time_3), ]
dframe2


pca_plot <-ggplot(dframe2, aes(x=Time3, y=TRUNK, group=Time_3)) + 
  geom_boxplot(aes(group = cut_width(Time_3, 0.25)),
               fill = "lightgrey", 
               lwd=0.3,
               outlier.alpha = 0.2,
               alpha=0.4) + 
  geom_jitter(width = 0.2, alpha=0.25)+
  scale_x_reverse (limits=c(255, 66))+
  labs(x = "Time (Ma)",
       y = "Trunk length (cm)")+
  theme_light()

pca_plot


#Format plot:
pca_plot  + theme (text=element_text(),
                   axis.text.x = element_text(size=9.5, color="#303030", hjust=0.5),
                   axis.text.y = element_text(size=9.5,color="#303030"),
                   axis.title.x  = element_blank(),
                   axis.title.y  = element_text(size=11,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                   strip.text = element_text(size=9.5,color="#303030"), #Change size of facet title
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.ticks = element_line(colour = "#303030"),
                   panel.border = element_rect(color="#303030", fill= NA, size=0.5),
                   aspect.ratio = 0.5,
                   legend.position = "none")

#will reduce to half size, so I use double size of font
ggsave("Z_Trunk size 13 time bins.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 14, height = 7, units = "cm",
       dpi = 300, limitsize = TRUE) 

#ggsave("Z_Trunk size-log10.pdf", plot = last_plot(), device = "pdf", path = NULL,
#       scale = 1, width = 14, height = 7, units = "cm",
#       dpi = 300, limitsize = TRUE) 



##1.6_TRUNK SIZE MAIN GROUPS - PANNELS MIXED POINTS + BOXPLOT (WITH 11 TIME BINS - FIGURE-)
#-------------------------------------------------------------------------------------#

#Before plotting FILTER THE DATA THAT HAVE <NA> in a given column
dframe2  <- dframe1[!is.na(dframe1$Time_4), ]
dframe2


#Create a subset of the data only with the group categories that I want to plot
dframe2<- subset(dframe2,
                 Group %in% c("Thalattosuchia", "Thalattosauria", "Sauropterygia","Pantestudines","Mosasauroidea","Ichthyosauromorpha"))
dframe2

#Making a factor and indicating the right order (levels)
Group1 <- factor(dframe2$Group, 
                levels = c("Ichthyosauromorpha","Mosasauroidea", "Pantestudines", "Sauropterygia", "Thalattosauria","Thalattosuchia"))
Group1


pca_plot <- dframe2 %>%
  ggplot(aes(fill = Group1, shape=Group1)) +
  geom_point(aes(x=Age, y=TRUNK, fill = Group1, shape=Group1),size=3, stroke=0.4, alpha=0.35) +
  geom_boxplot(aes(x=Time4, y=TRUNK, group = cut_width(Time_4, 0.25)),lwd=0.4, outlier.shape = NA) + 
  scale_x_reverse (limits=c(255, 66))+
  #scale_y_log10()+
  scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#19D27E","#FFD300","#3944bc"))+
  scale_shape_manual(values= c(24,23,21,21,25,22))+
  facet_wrap(~ factor(Group1),
             ncol =  1) + #all the plots in one single column)+
  labs(x = "Time (Mya)",
       y = "Trunk (cm)")+
  theme_light()

pca_plot


#Format plot:
pca_plot  + theme (text=element_text(),
                   axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                   axis.text.y = element_text(size=12,color="#303030"),
                   axis.title.x  = element_text(size=14,color="#303030", margin = margin(t = 10, r = 0, b = 0, l = 0)), #margin is added to distance title from axis
                   axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                   strip.text = element_blank(), # pannel titles
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.ticks = element_line(colour = "#303030"),
                   panel.border = element_rect(color="#303030", fill= NA, size=0.6),
                   aspect.ratio = 0.6,
                   legend.position = "none")

#will reduce to half size, so I use double size of font
ggsave("Z_Trunk length_Pannels-2-boxplot and dots.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 12, height = 28, units = "cm",
       dpi = 300, limitsize = TRUE) 



##1.7_TRUNK SIZE MAIN GROUPS - PANNELS MIXED POINTS + BOXPLOT (WITH 13 TIME BINS)
#------------------------------------------------------------------------------#

#Before plotting FILTER THE DATA THAT HAVE <NA> in a given column
dframe2  <- dframe1[!is.na(dframe1$Time_3), ]
dframe2

#Create a subset of the data only with the group categories that I want to plot
dframe2<- subset(dframe2,
                 Group %in% c("Thalattosuchia", "Thalattosauria", "Sauropterygia","Pantestudines","Mosasauroidea","Ichthyosauromorpha"))
dframe2

#Making a factor and indicating the right order (levels)
Group1 <- factor(dframe2$Group, 
                 levels = c("Ichthyosauromorpha","Mosasauroidea", "Pantestudines", "Sauropterygia", "Thalattosauria","Thalattosuchia"))
Group1

pca_plot <- dframe2 %>%
  ggplot(aes(fill = Group1, shape=Group1)) +
  geom_point(aes(x=Age, y=TRUNK, fill = Group1, shape=Group1),size=2.5, stroke=0.3, alpha=0.35) +
  geom_boxplot(aes(x=Time3, y=TRUNK, group = cut_width(Time_3, 0.22)),lwd=0.4, outlier.shape = NA) + 
  scale_x_reverse (limits=c(255, 66))+
  #scale_y_log10()+
  scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#19D27E","#FFD300","#3944bc"))+
  scale_shape_manual(values= c(24,23,21,21,25,22))+
  facet_wrap(~ factor(Group1),
             ncol =  1) + #all the plots in one single column)+
  labs(x = "Time (Mya)",
       y = "Trunk (cm)")+
  theme_light()

pca_plot


#Format plot:
pca_plot  + theme (text=element_text(),
                   axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                   axis.text.y = element_text(size=12,color="#303030"),
                   axis.title.x  = element_text(size=14,color="#303030", margin = margin(t = 10, r = 0, b = 0, l = 0)), #margin is added to distance title from axis
                   axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                   strip.text = element_blank(), # pannel titles
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.ticks = element_line(colour = "#303030"),
                   panel.border = element_rect(color="#303030", fill= NA, size=0.6),
                   aspect.ratio = 0.6,
                   legend.position = "none")

#will reduce to half size, so I use double size of font
ggsave("Z_Trunk length_Pannels-2-boxplot and dots.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 10, height = 25, units = "cm",
       dpi = 300, limitsize = TRUE) 


############################
############################


#IQR, MEDIAN and VARIANCE PER GROUP
#----------------------------------#
#----------------------------------#


dframe1<-read.csv("MMreptiles_Trunk.csv", #This file contains repeated taxa
                  row.names = 1,
                  header=T)
dframe1



#FILTER THE DATA THAT HAVE <NA> in the Time column (there has to be 194 obs)
dframe2  <- dframe1[!is.na(dframe1$Time), ]
dframe2


medIQR(TRUNK ~ Group, data = dframe2)
write.csv((medIQR(TRUNK ~ Group, data = dframe2)), "TRUNK_IQR_per_group.csv")

meanSD(TRUNK ~ Group, data = dframe2)
write.csv((meanSD(TRUNK ~ Group, data = dframe2)), "TRUNK_Var_per_group.csv")




#F test to compare variances
#https://rdrr.io/cran/RVAideMemoire/man/pairwise.var.test.html
pairwise.var.test(dframe2$TRUNK, dframe2$Group,  #Pairwise comparisons using F tests to compare two variances 
                  p.method = "bonferroni",
                  alternative = c("two.sided"))

#This comparison doesn't make much sense


#Create a subset of the data for each group
SA<- subset(dframe2,Group %in% c("Sauropterygia"))
SA
I<- subset(dframe2,Group %in% c("Ichthyosauromorpha"))
I
PT<- subset(dframe2,Group %in% c("Pantestudines"))
PT
MO<- subset(dframe2,Group %in% c("Mosasauroidea"))
MO
TSAU<- subset(dframe2,Group %in% c("Thalattosauria"))
TSAU
TSU<- subset(dframe2,Group %in% c("Thalattosuchia"))
TSU


# Shapiro-Wilk normality test 
shapiro.test(SA$TRUNK)
shapiro.test(I$TRUNK)
shapiro.test(MO$TRUNK)
shapiro.test(PT$TRUNK)
shapiro.test(TSAU$TRUNK)
shapiro.test(TSU$TRUNK)


#Unpaired two sample (non-parametric) test, when normality can't be assumed
#http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
#Not include these comparisons!

res <- wilcox.test(SA$TRUNK, I$TRUNK)
res

res <- wilcox.test(SA$TRUNK, MO$TRUNK)
res

res <- wilcox.test(I$TRUNK, MO$TRUNK)
res

res <- wilcox.test(I$TRUNK, PT$TRUNK)
res

res <- wilcox.test(SA$TRUNK, PT$TRUNK)
res

res <- wilcox.test(MO$TRUNK, PT$TRUNK)
res



#IQR, MEDIAN and VARIANCE PER BROAD TIME BIN
#---------------------------------------------#
#---------------------------------------------#

dframe1<-read.csv("MMreptiles_Trunk.csv", #This file contains repeated taxa
                  row.names = 1,
                  header=T)
dframe1


#FILTER THE DATA THAT HAVE <NA> in the Time column (there has to be 194 obs)
dframe2  <- dframe1[!is.na(dframe1$Time), ]
dframe2


#Create factor and specify order for time bins
A<-factor(dframe2$Time,
          levels = c("Triassic", "Jurassic","Cretaceous"))

A

medIQR(TRUNK ~ A, data = dframe2)
write.csv((medIQR(TRUNK ~ A, data = dframe2)), "TRUNK_IQR_TimebinsBroad-All.csv")

meanSD(TRUNK ~ A, data = dframe2)
write.csv((meanSD(TRUNK ~ A, data = dframe2)), "TRUNK_Var_TimebinsBroad-All.csv")



#Creating subsets
T<- subset(dframe1,Time %in% c("Triassic"))
T
J<- subset(dframe1,Time %in% c("Jurassic"))
J
C<- subset(dframe1,Time %in% c("Cretaceous"))
C

# Shapiro-Wilk normality test 
shapiro.test(T$TRUNK)
shapiro.test(J$TRUNK)
shapiro.test(C$TRUNK)

#If p-values are greater than the significance level 0.05 
#then the distribution of the data is not significantly different from normal distribution. 
#Hence, we can assume the normality



#F test to compare variances
#https://rdrr.io/cran/RVAideMemoire/man/pairwise.var.test.html

pairwise.var.test(dframe2$TRUNK, dframe2$Time,  #Pairwise comparisons using F tests to compare two variances 
                  p.method = "bonferroni",
                  alternative = c("two.sided"))

#Need to do Ansari test instead because most samples are not normally distributed
res<-ansari.test(J$TRUNK,T$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


res<-ansari.test(C$TRUNK,J$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


res<-ansari.test(C$TRUNK,T$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


#NOT DOING T-TEST!
#Unpaired two sample (non-parametric) test, when normality can't be assumed
#http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r


#res <- wilcox.test(T$TRUNK, J$TRUNK)
#res


#ALL GROUPS TOGETHER per 13 time-bins
#-----------------------------------#
#-----------------------------------#

dframe1<-read.csv("MMreptiles_Trunk.csv", #This file contains repeated taxa
                  row.names = 1,
                  header=T)
dframe1


#FILTER THE DATA THAT HAVE <NA> in the Time column (there has to be 194 obs)
dframe2  <- dframe1[!is.na(dframe1$Time_3), ]
dframe2


#Create factor and specify order for time bins
A<-factor(dframe2$Time_3,
          levels = c("Ind-Ole", "Ani-Lad","Car","Nor-Rha", "Het-Sin", "Pli-Toa", "Aal-Call", "Oxf", "Kim-Tith",
                     "Ber-Bar", "Apt-Alb","Cen-Con", "San-Maas"))

A

medIQR(TRUNK ~ A, data = dframe2)
write.csv((medIQR(TRUNK ~ A, data = dframe2)), "TRUNK_IQR_13Timebins-All.csv")

meanSD(TRUNK ~ A, data = dframe2)
write.csv((meanSD(TRUNK ~ A, data = dframe2)), "TRUNK_Var_13Timebins-All.csv")


#Creating subsets
IO<- subset(dframe2,Time_3 %in% c("Ind-Ole"))
IO
AL<- subset(dframe2,Time_3 %in% c("Ani-Lad"))
AL
C<- subset(dframe2,Time_3 %in% c("Car"))
C
NR<- subset(dframe2,Time_3 %in% c("Nor-Rha"))
NR
HS<- subset(dframe2,Time_3 %in% c("Het-Sin"))
HS
PT<- subset(dframe2,Time_3 %in% c("Pli-Toa"))
PT
AAC<- subset(dframe2,Time_3 %in% c("Aal-Call"))
AAC
O<- subset(dframe2,Time_3 %in% c("Oxf"))
O
KT<- subset(dframe2,Time_3 %in% c("Kim-Tith"))
KT
BB<- subset(dframe2,Time_3 %in% c("Ber-Bar"))
BB
AA<- subset(dframe2,Time_3 %in% c("Apt-Alb"))
AA
CC<- subset(dframe2,Time_3 %in% c("Cen-Con"))
CC
SM<- subset(dframe2,Time_3 %in% c("San-Maas"))
SM



#ALL GROUPS TOGETHER per 11 time-bins
#-----------------------------------#
#-----------------------------------#

dframe1<-read.csv("MMreptiles_Trunk.csv", #This file contains repeated taxa
                  row.names = 1,
                  header=T)
dframe1


#FILTER THE DATA THAT HAVE <NA> in the Time column (there has to be 194 obs)
dframe2  <- dframe1[!is.na(dframe1$Time_4), ]
dframe2


#Create factor and specify order for time bins
A<-factor(dframe2$Time_4,
          levels = c("Ind-Ole", "Ani-Lad","Car","Nor-Rha", "Het-Sin", "Pli-Toa", "Aal-Call", "Oxf-Tith",
                     "Ber-Alb","Cen-Con", "San-Maas"))

A

medIQR(TRUNK ~ A, data = dframe2)
write.csv((medIQR(TRUNK ~ A, data = dframe2)), "TRUNK_IQR_11Timebins-All.csv")

meanSD(TRUNK ~ A, data = dframe2)
write.csv((meanSD(TRUNK ~ A, data = dframe2)), "TRUNK_Var_11Timebins-All.csv")


#F test to compare variances OF NORMALLY DISTRIBUTED SAMPLES
#https://rdrr.io/cran/RVAideMemoire/man/pairwise.var.test.html

#pairwise.var.test(dframe2$TRUNK, dframe2$Time_4,  #Pairwise comparisons using F tests to compare two variances 
#                  p.method = "bonferroni",
#                  alternative = c("two.sided"))


#Creating subsets
IO<- subset(dframe2,Time_4 %in% c("Ind-Ole"))
IO
AL<- subset(dframe2,Time_4 %in% c("Ani-Lad"))
AL
C<- subset(dframe2,Time_4 %in% c("Car"))
C
NR<- subset(dframe2,Time_4 %in% c("Nor-Rha"))
NR
HS<- subset(dframe2,Time_4 %in% c("Het-Sin"))
HS
PT<- subset(dframe2,Time_4 %in% c("Pli-Toa"))
PT
AAC<- subset(dframe2,Time_4 %in% c("Aal-Call"))
AAC
OT<- subset(dframe2,Time_4 %in% c("Oxf-Tith"))
OT
BA<- subset(dframe2,Time_4 %in% c("Ber-Alb"))
BA
CC<- subset(dframe2,Time_4 %in% c("Cen-Con"))
CC
SM<- subset(dframe2,Time_4 %in% c("San-Maas"))
SM


# Shapiro-Wilk normality test 
shapiro.test(IO$TRUNK)
shapiro.test(AL$TRUNK)
shapiro.test(C$TRUNK)
shapiro.test(NR$TRUNK)
shapiro.test(HS$TRUNK)
shapiro.test(PT$TRUNK)
shapiro.test(AAC$TRUNK)
shapiro.test(OT$TRUNK)
shapiro.test(BA$TRUNK)
shapiro.test(CC$TRUNK)
shapiro.test(SM$TRUNK)
#If p-values are greater than the significance level 0.05 
#then the distribution of the data is not significantly different from normal distribution. 
#Hence, we can assume the normality

#

#COMPARING THE SPREAD OF THE DISTRIBUTIONS
#F test for comparing variances: assumes normal distribution
#I need to use a non-parametric test because most samples are not normal

var.test(IO$TRUNK, AL$TRUNK) 
res<-ansari.test(AL$TRUNK,IO$TRUNK,
            alternative = c("greater"),
            exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(C$TRUNK, AL$TRUNK) 
res<-ansari.test(C$TRUNK,AL$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(NR$TRUNK, C$TRUNK) 
res<-ansari.test(NR$TRUNK,C$TRUNK,
                 alternative = c("less"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res

var.test(NR$TRUNK, HS$TRUNK) 
res<-ansari.test(HS$TRUNK,NR$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(NR$TRUNK, HS$TRUNK) 
res<-ansari.test(HS$TRUNK,NR$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(PT$TRUNK,HS$TRUNK) 
res<-ansari.test(PT$TRUNK,HS$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res



var.test(AAC$TRUNK,PT$TRUNK) 
res<-ansari.test(AAC$TRUNK,PT$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res



var.test(OT$TRUNK,AAC$TRUNK) 
res<-ansari.test(OT$TRUNK,AAC$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(BA$TRUNK,OT$TRUNK) 
res<-ansari.test(BA$TRUNK,OT$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(CC$TRUNK,BA$TRUNK) 
res<-ansari.test(CC$TRUNK,BA$TRUNK,
                 alternative = c("less"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(SM$TRUNK,CC$TRUNK) 
res<-ansari.test(SM$TRUNK,CC$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


#Wilcox rank test

res <- wilcox.test(AL$TRUNK, IO$TRUNK,exact = FALSE, alternative = "greater")
res

res <- wilcox.test(C$TRUNK, AL$TRUNK,exact = FALSE, alternative = "greater")
res

res <- wilcox.test(NR$TRUNK, C$TRUNK,exact = FALSE, alternative = "less")
res

res <- t.test(HS$TRUNK, NR$TRUNK,exact = FALSE, alternative = "greater")
res

res <- wilcox.test(PT$TRUNK, HS$TRUNK,exact = FALSE, alternative = "greater")
res
  
res <- wilcox.test(AAC$TRUNK, PT$TRUNK,exact = FALSE, alternative = "greater")
res

res <- wilcox.test(OT$TRUNK, AAC$TRUNK,exact = FALSE, alternative = "less")
res

res <- wilcox.test(BA$TRUNK, OT$TRUNK,exact = FALSE, alternative = "greater")
res

res <- wilcox.test(CC$TRUNK, BA$TRUNK,exact = FALSE, alternative = "greater")
res

res <- wilcox.test(SM$TRUNK, CC$TRUNK,exact = FALSE, alternative = "less")
res



#SAUROPTERYGIA per 11 time-bins
#------------------------------#
#------------------------------#

dframe1<-read.csv("MMreptiles_Trunk.csv", #This file contains repeated taxa
                  row.names = 1,
                  header=T)
dframe1

#FILTER THE DATA THAT HAVE <NA> in the Time column
dframe2  <- dframe1[!is.na(dframe1$Time_4), ]
dframe2


#Create a subset of the data
dframe2<- subset(dframe2,Group %in% c("Sauropterygia"))
dframe2

#Create factor and specify order for time bins
A<-factor(dframe2$Time_4,
          levels = c("Ind-Ole", "Ani-Lad","Car","Nor-Rha", "Het-Sin", "Pli-Toa", "Aal-Call", "Oxf-Tith",
                     "Ber-Alb","Cen-Con", "San-Maas"))

A


medIQR(TRUNK ~ A, data = dframe2)
write.csv((medIQR(TRUNK ~ A, data = dframe2)), "TRUNK_IQR_Timebins_Sauropterygia.csv")

meanSD(TRUNK ~ A, data = dframe2)
write.csv((meanSD(TRUNK ~ A, data = dframe2)), "TRUNK_Var_Timebins_Sauropterygia.csv")


#Creating subsets
IO<- subset(dframe2,Time_4 %in% c("Ind-Ole"))
IO
AL<- subset(dframe2,Time_4 %in% c("Ani-Lad"))
AL
C<- subset(dframe2,Time_4 %in% c("Car"))
C
NR<- subset(dframe2,Time_4 %in% c("Nor-Rha"))
NR
HS<- subset(dframe2,Time_4 %in% c("Het-Sin"))
HS
PT<- subset(dframe2,Time_4 %in% c("Pli-Toa"))
PT
AAC<- subset(dframe2,Time_4 %in% c("Aal-Call"))
AAC
OT<- subset(dframe2,Time_4 %in% c("Oxf-Tith"))
OT
BA<- subset(dframe2,Time_4 %in% c("Ber-Alb"))
BA
CC<- subset(dframe2,Time_4 %in% c("Cen-Con"))
CC
SM<- subset(dframe2,Time_4 %in% c("San-Maas"))
SM

# Shapiro-Wilk normality test 

shapiro.test(IO$TRUNK)
shapiro.test(AL$TRUNK)
shapiro.test(C$TRUNK)
shapiro.test(NR$TRUNK)
shapiro.test(HS$TRUNK)
shapiro.test(PT$TRUNK)
shapiro.test(AAC$TRUNK)
shapiro.test(OT$TRUNK)
shapiro.test(BA$TRUNK)
shapiro.test(CC$TRUNK)
shapiro.test(SM$TRUNK)



#If p-values are greater than the significance level 0.05 
#then the distribution of the data is not significantly different from normal distribution. 
#Hence, we can assume the normality


#COMPARING THE SPREAD OF THE DISTRIBUTIONS
#F test for comparing variances: assumes normal distribution
#I need to use a non-parametric test because most samples are not normal

var.test(C$TRUNK, AL$TRUNK) 
res<-ansari.test(C$TRUNK,AL$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res

var.test(HS$TRUNK, C$TRUNK) 
res<-ansari.test(HS$TRUNK,C$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res

var.test(PT$TRUNK, HS$TRUNK) 
res<-ansari.test(PT$TRUNK,HS$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res

var.test(AAC$TRUNK, PT$TRUNK) 
res<-ansari.test(AAC$TRUNK, PT$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(OT$TRUNK, AAC$TRUNK) 
res<-ansari.test(OT$TRUNK, AAC$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(BA$TRUNK, OT$TRUNK) 
res<-ansari.test(BA$TRUNK, OT$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(CC$TRUNK, BA$TRUNK) 
res<-ansari.test(CC$TRUNK, CC$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res

var.test(SM$TRUNK, CC$TRUNK) 
res<-ansari.test(SM$TRUNK, CC$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res




#8 TIME BINS
#------------#

dframe1<-read.csv("MMreptiles_Trunk.csv", #This file contains repeated taxa
                  row.names = 1,
                  header=T)
dframe1

#FILTER THE DATA THAT HAVE <NA> in the Time column
dframe2  <- dframe1[!is.na(dframe1$Time_2), ]
dframe2


#Create a subset of the data
dframe2<- subset(dframe2,Group %in% c("Sauropterygia"))
dframe2


A<-factor(dframe2$Time_2,
          levels = c("Early_Triassic", "Mid_Triassic", "Late_Triassic", "Early_Jurassic", "Mid_Jurassic", "Late_Jurassic", "Early_Cretaceous", "Late_Cretaceous"))

A



medIQR(TRUNK ~ A, data = dframe2)
write.csv((medIQR(TRUNK ~ A, data = dframe2)), "TRUNK_IQR_8Timebins_Sauropterygia.csv")

meanSD(TRUNK ~ A, data = dframe2)
write.csv((meanSD(TRUNK ~ A, data = dframe2)), "TRUNK_Var_8Timebins_Sauropterygia.csv")



#Creating subsets
ET<- subset(dframe2,Time_2 %in% c("Early_Triassic"))
ET
MT<- subset(dframe2,Time_2 %in% c("Mid_Triassic"))
MT
LT<- subset(dframe2,Time_2 %in% c("Late_Triassic"))
LT
EJ<- subset(dframe2,Time_2 %in% c("Early_Jurassic"))
EJ
MJ<- subset(dframe2,Time_2 %in% c("Mid_Jurassic"))
MJ
LJ<- subset(dframe2,Time_2 %in% c("Late_Jurassic"))
LJ
EC<- subset(dframe2,Time_2 %in% c("Early_Cretaceous"))
EC
LC<- subset(dframe2,Time_2 %in% c("Late_Cretaceous"))
LC


# Shapiro-Wilk normality test 
shapiro.test(MT$TRUNK)
shapiro.test(LT$TRUNK)
shapiro.test(EJ$TRUNK)
shapiro.test(MJ$TRUNK)
shapiro.test(LJ$TRUNK)
shapiro.test(EC$TRUNK)
shapiro.test(LC$TRUNK)


#If p-values are greater than the significance level 0.05 
#then the distribution of the data is not significantly different from normal distribution. 
#Hence, we can assume the normality


#COMPARING THE SPREAD OF THE DISTRIBUTIONS
#F test for comparing variances: assumes normal distribution
#I need to use a non-parametric test because most samples are not normal

var.test(LT$TRUNK, MT$TRUNK) 
res<-ansari.test(LT$TRUNK,MT$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res

var.test(EJ$TRUNK, LT$TRUNK) 
res<-ansari.test(EJ$TRUNK,LT$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(MJ$TRUNK, EJ$TRUNK) 
res<-ansari.test(MJ$TRUNK,EJ$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res



var.test(LJ$TRUNK, MJ$TRUNK) 
res<-ansari.test(LJ$TRUNK,MJ$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(EC$TRUNK, LJ$TRUNK) 
res<-ansari.test(EC$TRUNK,LJ$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res

var.test(LC$TRUNK, EC$TRUNK) 
res<-ansari.test(LC$TRUNK,EC$TRUNK,
                 alternative = c("less"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(EC$TRUNK, EJ$TRUNK) 
res<-ansari.test(EC$TRUNK,EJ$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res



#Testing homodasticity (equal variances)

library(tidyverse)

dframe2 <- subset(dframe2, A != "Early_Triassic") #Drop an unused Time-bin
res <- bartlett.test(TRUNK ~ Time_2, data = dframe2) #Bartlett’s test with one independent variable. It gives you a lump value, not by pairs
res
# if p< 0.05, variance of trunk is significantly different between the time bins

library(car)
res <-leveneTest(TRUNK ~ Time_2, data = dframe2)
res
# if p< 0.05, variance of trunk is significantly different between the time bins



#Unpaired two sample (non-parametric) test, when normality can't be assumed
#http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r


res <- wilcox.test(MT$TRUNK, LT$TRUNK,exact = FALSE, alternative = "less")
res

res <- wilcox.test(LT$TRUNK, EJ$TRUNK,exact = FALSE, alternative = "less")
res

res <- wilcox.test(EJ$TRUNK, MJ$TRUNK,exact = FALSE, alternative = "less")
res

res <- wilcox.test(MJ$TRUNK, LJ$TRUNK,exact = FALSE, alternative = "less")
res

res <- wilcox.test(EC$TRUNK, EJ$TRUNK,exact = FALSE, alternative = "less")
res

## t-test
#Brief description / how to use it:
#by default the Welch test is used (in which, unlike Student t-test, equal variance is not assumed)
#t.test(extra ~ group, sleep)                   #one column "extra" records the measurement and the other "group" specifies the grouping
# t.test(sleep_wide$group1, sleep_wide$group2)   # Same for wide data (two separate vectors)


t.test(LJ$TRUNK, EC$TRUNK,
       var.equal=F)
t.test(EC$TRUNK, LC$TRUNK,
       var.equal=F)

t.test(EC$TRUNK, EJ$TRUNK,
       var.equal=F)



#ICHTHYOSAUROMORPHA per 11 time bins
#-----------------------------------#
#-----------------------------------#


dframe1<-read.csv("MMreptiles_Trunk.csv", #This file contains repeated taxa
                  row.names = 1,
                  header=T)
dframe1

#FILTER THE DATA THAT HAVE <NA> in the Time column
dframe2  <- dframe1[!is.na(dframe1$Time_4), ]
dframe2


#Create a subset of the data
dframe2<- subset(dframe2,Group %in% c("Ichthyosauromorpha"))
dframe2

#Create factor and specify order for time bins
A<-factor(dframe2$Time_4,
          levels = c("Ind-Ole", "Ani-Lad","Car","Nor-Rha", "Het-Sin", "Pli-Toa", "Aal-Call", "Oxf-Tith",
                     "Ber-Alb","Cen-Con", "San-Maas"))

A

medIQR(TRUNK ~ A, data = dframe2)
write.csv((medIQR(TRUNK ~ A, data = dframe2)), "TRUNK_IQR_Timebins_Ichthyosauromorpha.csv")

meanSD(TRUNK ~ A, data = dframe2)
write.csv((meanSD(TRUNK ~ A, data = dframe2)), "TRUNK_Var_Timebins_Ichthyosauromorpha.csv")


#Creating subsets
IO<- subset(dframe2,Time_4 %in% c("Ind-Ole"))
IO
AL<- subset(dframe2,Time_4 %in% c("Ani-Lad"))
AL
C<- subset(dframe2,Time_4 %in% c("Car"))
C
NR<- subset(dframe2,Time_4 %in% c("Nor-Rha"))
NR
HS<- subset(dframe2,Time_4 %in% c("Het-Sin"))
HS
PT<- subset(dframe2,Time_4 %in% c("Pli-Toa"))
PT
AAC<- subset(dframe2,Time_4 %in% c("Aal-Call"))
AAC
OT<- subset(dframe2,Time_4 %in% c("Oxf-Tith"))
OT
BA<- subset(dframe2,Time_4 %in% c("Ber-Alb"))
BA


# Shapiro-Wilk normality test 

shapiro.test(IO$TRUNK)
shapiro.test(AL$TRUNK)
shapiro.test(C$TRUNK)
shapiro.test(NR$TRUNK)
shapiro.test(HS$TRUNK)
shapiro.test(PT$TRUNK)
shapiro.test(AAC$TRUNK)
shapiro.test(OT$TRUNK)
shapiro.test(BA$TRUNK)


#If p-values are greater than the significance level 0.05 
#then the distribution of the data is not significantly different from normal distribution. 
#Hence, we can assume the normality


#COMPARING THE SPREAD OF THE DISTRIBUTIONS
#F test for comparing variances: assumes normal distribution
#I need to use a non-parametric test because most samples are not normal

var.test(AL$TRUNK, IO$TRUNK) 
res<-ansari.test(AL$TRUNK, IO$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res

var.test(C$TRUNK, AL$TRUNK) 
res<-ansari.test(C$TRUNK, AL$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res

var.test(HS$TRUNK, C$TRUNK) 
res<-ansari.test(HS$TRUNK, C$TRUNK,
                 alternative = c("less"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(PT$TRUNK, HS$TRUNK) 
res<-ansari.test(PT$TRUNK, HS$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res

var.test(AAC$TRUNK, PT$TRUNK) 
res<-ansari.test(AAC$TRUNK, PT$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res

var.test(OT$TRUNK, AAC$TRUNK) 
res<-ansari.test(OT$TRUNK, AAC$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(BA$TRUNK, OT$TRUNK) 
res<-ansari.test(BA$TRUNK, OT$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res



#8 time bins
#----------#

dframe1<-read.csv("MMreptiles_Trunk.csv", #This file contains repeated taxa
                  row.names = 1,
                  header=T)
dframe1

#FILTER THE DATA THAT HAVE <NA> in the Time column
dframe2  <- dframe1[!is.na(dframe1$Time_2), ]
dframe2


#Create a subset of the data
dframe2<- subset(dframe2,Group %in% c("Ichthyosauromorpha"))
dframe2

A<-factor(dframe2$Time_2,
          levels = c("Early_Triassic", "Mid_Triassic", "Late_Triassic", "Early_Jurassic", "Mid_Jurassic", "Late_Jurassic", "Early_Cretaceous", "Late_Cretaceous"))

A


medIQR(TRUNK ~ A, data = dframe2)
write.csv((medIQR(TRUNK ~ A, data = dframe2)), "TRUNK_IQR_8Timebins_Ichthyosauromorpha.csv")

meanSD(TRUNK ~ A, data = dframe2)
write.csv((meanSD(TRUNK ~ A, data = dframe2)), "TRUNK_Var_8Timebins_Ichthyosauromorpha.csv")


#Creating subsets
ET<- subset(dframe2,Time_2 %in% c("Early_Triassic"))
ET
MT<- subset(dframe2,Time_2 %in% c("Mid_Triassic"))
MT
LT<- subset(dframe2,Time_2 %in% c("Late_Triassic"))
LT
EJ<- subset(dframe2,Time_2 %in% c("Early_Jurassic"))
EJ
MJ<- subset(dframe2,Time_2 %in% c("Mid_Jurassic"))
MJ
LJ<- subset(dframe2,Time_2 %in% c("Late_Jurassic"))
LJ
EC<- subset(dframe2,Time_2 %in% c("Early_Cretaceous"))
EC



# Shapiro-Wilk normality test 
shapiro.test(ET$TRUNK)
shapiro.test(MT$TRUNK)
shapiro.test(LT$TRUNK)
shapiro.test(EJ$TRUNK)
shapiro.test(MJ$TRUNK)
shapiro.test(LJ$TRUNK)
shapiro.test(EC$TRUNK)

#If p-values are greater than the significance level 0.05 
#then the distribution of the data is not significantly different from normal distribution. 
#Hence, we can assume the normality


var.test(MT$TRUNK, ET$TRUNK) 
res<-ansari.test(MT$TRUNK, ET$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res


var.test(LT$TRUNK, MT$TRUNK) 
res<-ansari.test(LT$TRUNK, MT$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res



var.test(LT$TRUNK, ET$TRUNK) 
res<-ansari.test(LT$TRUNK, ET$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res



var.test(EJ$TRUNK, LT$TRUNK) 
res<-ansari.test(EJ$TRUNK, LT$TRUNK,
                 alternative = c("less"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res

var.test(LJ$TRUNK, EJ$TRUNK) 
res<-ansari.test(LJ$TRUNK, EJ$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res

var.test(EC$TRUNK, LJ$TRUNK) 
res<-ansari.test(EC$TRUNK, LJ$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res



#Testing homodasticity (equal variances)

dframe2 <- subset(dframe2, A != "Mid_Jurassic") #Drop an unused Time-bin
dframe2 <- subset(dframe2, A != "Early-Late_Cretaceous") #Drop an unused Time-bin
dframe2 <- subset(dframe2, A != "Late-Late_Cretaceous") #Drop an unused Time-bin


res <- bartlett.test(TRUNK ~ Time_2, data = dframe2) #Bartlett’s test with one independent variable. It gives you a lump value, not by pairs
res
# if p< 0.05, variance of trunk is significantly different between the time bins

res <-leveneTest(TRUNK ~ Time_2, data = dframe2)
res
# if p< 0.05, variance of trunk is significantly different between the time bins


## t-test
#Brief description / how to use it:
#by default the Welch test is used (in which, unlike Student t-test, equal variance is not assumed)
#t.test(extra ~ group, sleep)                   #one column "extra" records the measurement and the other "group" specifies the grouping
# t.test(sleep_wide$group1, sleep_wide$group2)   # Same for wide data (two separate vectors)

t.test(ET$TRUNK, MT$TRUNK,
       var.equal=F)
t.test(MT$TRUNK, LT$TRUNK,
       var.equal=F)
t.test(ET$TRUNK, LT$TRUNK,
       var.equal=F)

t.test(LT$TRUNK, EJ$TRUNK,
       var.equal=F)
t.test(EJ$TRUNK, LJ$TRUNK,
       var.equal=F)

t.test(LJ$TRUNK, EC$TRUNK,
       var.equal=F)

#Unpaired two sample (non-parametric) test, when normality can't be assumed
#http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r

res <- wilcox.test(ET$TRUNK, MT$TRUNK, exact = FALSE, alternative = "less")
res

res <- wilcox.test(MT$TRUNK, LT$TRUNK, exact = FALSE, alternative = "less")
res

res <- wilcox.test(ET$TRUNK, LT$TRUNK, exact = FALSE, alternative = "less")
res

res <- wilcox.test(LT$TRUNK, EJ$TRUNK, exact = FALSE, alternative = "less")
res

res <- wilcox.test(EJ$TRUNK, LJ$TRUNK,exact = FALSE, alternative = "less")
res


res <- wilcox.test(LJ$TRUNK, EC$TRUNK,exact = FALSE, alternative = "less")
res



#MOSASAUROIDEA per 11 time bin
#---------------------------#

dframe1<-read.csv("MMreptiles_Trunk.csv", #This file contains repeated taxa
                 row.names = 1,
                 header=T)
dframe1

#FILTER THE DATA THAT HAVE <NA> in the Time column
dframe2  <- dframe1[!is.na(dframe1$Time_4), ]
dframe2


#Create a subset of the data
dframe2<- subset(dframe2,Group %in% c("Mosasauroidea"))
dframe2


#Create factor and specify order for time bins
A<-factor(dframe2$Time_4,
          levels = c("Ind-Ole", "Ani-Lad","Car","Nor-Rha", "Het-Sin", "Pli-Toa", "Aal-Call", "Oxf-Tith",
                     "Ber-Alb","Cen-Con", "San-Maas"))

A


medIQR(TRUNK ~ A, data = dframe2)
write.csv((medIQR(TRUNK ~ A, data = dframe2)), "TRUNK_IQR_Timebins_Mosasauroidea.csv")

meanSD(TRUNK ~ A, data = dframe2)
write.csv((meanSD(TRUNK ~ A, data = dframe2)), "TRUNK_Var_Timebins_Mosasauroidea.csv")


#Creating subsets
CC<- subset(dframe2,Time_4 %in% c("Cen-Con"))
CC
SM<- subset(dframe2,Time_4 %in% c("San-Maas"))
SM

# Shapiro-Wilk normality test 
shapiro.test(CC$TRUNK)
shapiro.test(SM$TRUNK)

#If p-values are greater than the significance level 0.05 
#then the distribution of the data is not significantly different from normal distribution. 
#Hence, we can assume the normality


var.test(SM$TRUNK, CC$TRUNK) 
res<-ansari.test(SM$TRUNK, CC$TRUNK,
                 alternative = c("less"),
                 exact = FALSE, conf.int = F, conf.level = 0.95) 
res

res<-mood.test(SM$TRUNK, CC$TRUNK,
                 alternative = c("greater"),
                 exact = NULL, conf.int = FALSE, conf.level = 0.95) 
res



#Unpaired two sample (non-parametric) test, when normality can't be assumed
#http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r

t.test(CC$TRUNK, SM$TRUNK,
       var.equal=F)

res <- wilcox.test(CC$TRUNK, SM$TRUNK,exact = FALSE, alternative = "less")
res



#PANTESTUDINES per 11 time bin
#---------------------------#

rm(list = ls())

dframe1<-read.csv("MMreptiles_Trunk.csv", #This file contains repeated taxa
                  row.names = 1,
                  header=T)
dframe1

#FILTER THE DATA THAT HAVE <NA> in the Time column
dframe2  <- dframe1[!is.na(dframe1$Time_4), ]
dframe2


#Create a subset of the data
dframe2<- subset(dframe2,Group %in% c("Pantestudines"))
dframe2


#Create factor and specify order for time bins
A<-factor(dframe2$Time_4,
          levels = c("Ind-Ole", "Ani-Lad","Car","Nor-Rha", "Het-Sin", "Pli-Toa", "Aal-Call", "Oxf-Tith",
                     "Ber-Alb","Cen-Con", "San-Maas"))

A

medIQR(TRUNK ~ A, data = dframe2)
write.csv((medIQR(TRUNK ~ A, data = dframe2)), "TRUNK_IQR_Timebins_Pantestudines.csv")

meanSD(TRUNK ~ A, data = dframe2)
write.csv((meanSD(TRUNK ~ A, data = dframe2)), "TRUNK_Var_Timebins_Pantestudines.csv")


#Creating subsets
C<- subset(dframe2,Time_4 %in% c("Car"))
C
NR<- subset(dframe2,Time_4 %in% c("Nor-Rha"))
NR
AAC<- subset(dframe2,Time_4 %in% c("Aal-Call"))
AAC
OT<- subset(dframe2,Time_4 %in% c("Oxf-Tith"))
OT
BA<- subset(dframe2,Time_4 %in% c("Ber-Alb"))
BA
CC<- subset(dframe2,Time_4 %in% c("Cen-Con"))
CC
SM<- subset(dframe2,Time_4 %in% c("San-Maas"))
SM


#8 time bin
#----------#

rm(list = ls())

dframe1<-read.csv("MMreptiles_Trunk.csv", #This file contains repeated taxa
                  row.names = 1,
                  header=T)
dframe1

#FILTER THE DATA THAT HAVE <NA> in the Time column
dframe2  <- dframe1[!is.na(dframe1$Time_2), ]
dframe2


#Create a subset of the data
dframe2<- subset(dframe2,Group %in% c("Pantestudines"))
dframe2

A<-factor(dframe2$Time_2,
          levels = c("Early_Triassic", "Mid_Triassic", "Late_Triassic", "Early_Jurassic", "Mid_Jurassic", "Late_Jurassic", "Early_Cretaceous", "Late_Cretaceous"))

A

medIQR(TRUNK ~ A, data = dframe2)
write.csv((medIQR(TRUNK ~ A, data = dframe2)), "TRUNK_IQR_8Timebins_Pantestudines.csv")

meanSD(TRUNK ~ A, data = dframe2)
write.csv((meanSD(TRUNK ~ A, data = dframe2)), "TRUNK_Var_8Timebins_Pantestudines.csv")



#Creating subsets
ET<- subset(dframe2,Time_2 %in% c("Early_Triassic"))
ET
MT<- subset(dframe2,Time_2 %in% c("Mid_Triassic"))
MT
LT<- subset(dframe2,Time_2 %in% c("Late_Triassic"))
LT
EJ<- subset(dframe2,Time_2 %in% c("Early_Jurassic"))
EJ
MJ<- subset(dframe2,Time_2 %in% c("Mid_Jurassic"))
MJ
LJ<- subset(dframe2,Time_2 %in% c("Late_Jurassic"))
LJ
EC<- subset(dframe2,Time_2 %in% c("Early_Cretaceous"))
EC
LC<- subset(dframe2,Time_2 %in% c("Late_Cretaceous"))
LC



# Shapiro-Wilk normality test 
shapiro.test(LT$TRUNK)
shapiro.test(EJ$TRUNK)
shapiro.test(MJ$TRUNK)
shapiro.test(LJ$TRUNK)
shapiro.test(EC$TRUNK)
shapiro.test(LC$TRUNK)

#If p-values are greater than the significance level 0.05 
#then the distribution of the data is not significantly different from normal distribution. 
#Hence, we can assume the normality

var.test(LT$TRUNK, LJ$TRUNK)
var.test(LJ$TRUNK, EC$TRUNK) 
var.test(EC$TRUNK, LC$TRUNK) 




#Testing homodasticity (equal variances)
dframe2 <- subset(dframe2, A != "Early_Triassic") #Drop an unused Time-bin
dframe2 <- subset(dframe2, A != "Mid_Triassic") #Drop an unused Time-bin
dframe2 <- subset(dframe2, A != "Early_Jurassic") #Drop an unused Time-bin

res <- bartlett.test(TRUNK ~ Time_2, data = dframe2) #Bartlett’s test with one independent variable. It gives you a lump value, not by pairs
res


# if p< 0.05, variance of trunk is significantly different between the time bins

res <-leveneTest(TRUNK ~ Time_2, data = dframe2)
res
# if p< 0.05, variance of trunk is significantly different between the time bins


## t-test
#Brief description / how to use it:
#by default the Welch test is used (in which, unlike Student t-test, equal variance is not assumed)
#t.test(extra ~ group, sleep)                   #one column "extra" records the measurement and the other "group" specifies the grouping
# t.test(sleep_wide$group1, sleep_wide$group2)   # Same for wide data (two separate vectors)


t.test(LT$TRUNK, LJ$TRUNK,
       var.equal=F)
t.test(LJ$TRUNK, EC$TRUNK,
       var.equal=F)
t.test(EC$TRUNK, LC$TRUNK,
       var.equal=F)





#THALATTOSUCHIA per time bin
#---------------------------#

rm(list = ls())

dframe1<-read.csv("MMreptiles_Trunk.csv", #This file contains repeated taxa
                  row.names = 1,
                  header=T)
dframe1

#FILTER THE DATA THAT HAVE <NA> in the Time column
dframe2  <- dframe1[!is.na(dframe1$Time_4), ]
dframe2


#Create a subset of the data
dframe2<- subset(dframe2,Group %in% c("Thalattosuchia"))
dframe2

A<-factor(dframe2$Time_2,
          levels = c("Early_Triassic", "Mid_Triassic", "Late_Triassic", "Early_Jurassic", "Mid_Jurassic", "Late_Jurassic", "Early_Cretaceous", "Late_Cretaceous"))

A

medIQR(TRUNK ~ A, data = dframe2)
write.csv((medIQR(TRUNK ~ A, data = dframe2)), "TRUNK_IQR_8Timebins_Thalattosuchia.csv")

meanSD(TRUNK ~ A, data = dframe2)
write.csv((meanSD(TRUNK ~ A, data = dframe2)), "TRUNK_Var_8Timebins_Thalattosuchia.csv")



#Creating subsets
EJ<- subset(dframe2,Time_2 %in% c("Early_Jurassic"))
EJ
MJ<- subset(dframe2,Time_2 %in% c("Mid_Jurassic"))
MJ
LJ<- subset(dframe2,Time_2 %in% c("Late_Jurassic"))
LJ
EC<- subset(dframe2,Time_2 %in% c("Early_Cretaceous"))
EC


# Shapiro-Wilk normality test 
shapiro.test(EJ$TRUNK)
shapiro.test(MJ$TRUNK)
shapiro.test(LJ$TRUNK)
shapiro.test(ELC$TRUNK)

var.test(LJ$TRUNK, EJ$TRUNK)

#If p-values are greater than the significance level 0.05 
#then the distribution of the data is not significantly different from normal distribution. 
#Hence, we can assume the normality


#Testing homodasticity (equal variances)

res <-leveneTest(TRUNK ~ Time_2, data = dframe2)
res
# if p< 0.05, variance of trunk is significantly different between the time bins


## t-test
#Brief description / how to use it:
#by default the Welch test is used (in which, unlike Student t-test, equal variance is not assumed)
#t.test(extra ~ group, sleep)                   #one column "extra" records the measurement and the other "group" specifies the grouping
# t.test(sleep_wide$group1, sleep_wide$group2)   # Same for wide data (two separate vectors)


t.test(EJ$TRUNK, LJ$TRUNK,
       var.equal=F)




#THALATTOSAURIA per 8 time bin
#---------------------------#

rm(list = ls())

dframe1<-read.csv("MMreptiles_Trunk.csv", #This file contains repeated taxa
                  row.names = 1,
                  header=T)
dframe1

#FILTER THE DATA THAT HAVE <NA> in the Time column
dframe2  <- dframe1[!is.na(dframe1$Time_2), ]
dframe2

#Create a subset of the data
dframe2<- subset(dframe2,Group %in% c("Thalattosauria"))
dframe2

A<-factor(dframe2$Time_2,
          levels = c("Early_Triassic", "Mid_Triassic", "Late_Triassic", "Early_Jurassic", "Mid_Jurassic", "Late_Jurassic", "Early_Cretaceous", "Late_Cretaceous"))

A

medIQR(TRUNK ~ A, data = dframe2)
write.csv((medIQR(TRUNK ~ A, data = dframe2)), "TRUNK_IQR_8Timebins_Thalattosauria.csv")

meanSD(TRUNK ~ A, data = dframe2)
write.csv((meanSD(TRUNK ~ A, data = dframe2)), "TRUNK_Var_8Timebins_Thalattosauria.csv")


#Creating subsets
MT<- subset(dframe2,Time_2 %in% c("Mid_Triassic"))
MT
LT<- subset(dframe2,Time_2 %in% c("Late_Triassic"))
LT

# Shapiro-Wilk normality test 
shapiro.test(MT$TRUNK)
shapiro.test(LT$TRUNK)

#If p-values are greater than the significance level 0.05 
#then the distribution of the data is not significantly different from normal distribution. 
#Hence, we can assume the normality

var.test(MT$TRUNK, LT$TRUNK)



#Testing homodasticity (equal variances)

res <- bartlett.test(TRUNK ~ Time_2, data = dframe2) #Bartlett’s test with one independent variable. It gives you a lump value, not by pairs
res
# if p< 0.05, variance of trunk is significantly different between the time bins


## t-test
#Brief description / how to use it:
#by default the Welch test is used (in which, unlike Student t-test, equal variance is not assumed)
#t.test(extra ~ group, sleep)                   #one column "extra" records the measurement and the other "group" specifies the grouping
# t.test(sleep_wide$group1, sleep_wide$group2)   # Same for wide data (two separate vectors)


t.test(MT$TRUNK, LT$TRUNK,
       var.equal=F)



##########################
###########################
