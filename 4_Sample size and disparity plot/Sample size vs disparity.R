
#LOAD PACKAGES
library(ggplot2)
#install.packages("summarize", repos="http://R-Forge.R-project.org")
library(summarize)

rm(list = ls())

#LOCOMOTORY DISPARITY: Disparity and sample size (11 TIME BINS)
#--------------------------------------------------------------#

#PCA_Morphospace_Disparity-PManova: calculate the disparity for all taxa divided into 11 time bins


#assigning numerical values to the ages
ages<- c(251, 242.1, 232, 214.1, 196, 182.4, 168.8, 154.2, 122.3, 93.4, 76.1)

#Creatig a column with the median (disparity obs)
#Disparity <- c(1.5,6.11,8.03,2.79,8.67,12.43,13.59,14.68,17.83,14.66,15.63)

Disparity <- c(1.42,6.15,8.26,2.72,8.86,12.5,15.3,14.87,17.8,11.5,15.82)
Disparity



#Creating a column with the sample size
n<-c(6,27,16,7,10,18,7,19,8,7,14)
n

df2 <- data.frame(ages,Disparity,n) #creates a dataframe with various columns
df2

pca_plot <- ggplot(df2) +
  geom_line(aes(x=ages, y=Disparity)) +
  geom_line(aes(x=ages, y = n, colour = "Red"))+
  scale_x_reverse (limits=c(255, 66))+
  scale_y_continuous(sec.axis = sec_axis(~ .*1, name = "Sample size")) +#adds a secondary axis
  labs(x = "Time (Mya)",
       y = "Disparity (sum of variances)")+
  theme_light()

pca_plot

#Format plot:
pca_plot  + theme (text=element_text(),
                   axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                   axis.text.y.right = element_text(size=12,color="#303030"),   
                   axis.text.y.left = element_text(size=12,color="#303030"),   
                   axis.title.x  = element_blank(), #margin is added to distance title from axis
                   axis.title.y.right = element_text(size=14,color="#303030", margin = margin(t = 0, r = 0, b = 0, l = 10)), #margin is added to distance title from axis
                   axis.title.y.left = element_text(size=14,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.ticks = element_line(colour = "#303030"),
                   panel.border = element_rect(color="#303030", fill= NA, size=0.6),
                   aspect.ratio = 0.5,
                   legend.position = "none")



#will reduce to half size, so I use double size of font
ggsave("Z_Sample size_locomotory disparity.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 16, height = 10, units = "cm",
       dpi = 300, limitsize = TRUE) 



#TRUNK LENGTH: Standard deviation and sample size (11 TIME BINS)
#--------------------------------------------------------------#

dframe1<-read.csv("MMreptiles_Trunk.csv",   #MMreptiles_Trunk.csv with 8 time bins
                  row.names = 1,
                  header=T)

dframe1

#Before plotting FILTER THE DATA THAT HAVE <NA> in a given column
dframe2  <- dframe1[!is.na(dframe1$Time_4), ]
dframe2

#Calculate the SD
dframe3<-meanSD(TRUNK ~ A, data = dframe2)
dframe3
#dframe3<-t(dframe3) #transpose = switch columns and rows
#dframe3

#assigning numerical values to the ages
ages<- c(251, 242.1, 232, 214.1, 196, 182.4, 168.8, 154.2, 122.3, 93.4, 76.1)

#Creatig a column with the SD
SD <- c(18.70557, 65.27148, 83.56136, 27.64126, 30.60752,  54.7407,  57.1060,73.81531, 111.3429, 
        83.55539, 106.7301)
SD

#Creating a column with the sample size
n<-c(10,37,26,11,14,24,7,27,23,13,24)
n

df <- data.frame(ages,SD,n) #creates a dataframe with various columns
df

pca_plot <- ggplot(df) +
  geom_line(aes(x=ages, y=SD)) +
  geom_line(aes(x=ages, y = n, colour = "Red"))+
  scale_x_reverse (limits=c(255, 66))+
  scale_y_continuous(sec.axis = sec_axis(~ .*1, name = "Sample size")) +#adds a secondary axis
  labs(x = "Time (Mya)",
       y = "Trunk length dispersion (sd)")+
  theme_light()

pca_plot


#Format plot:
pca_plot  + theme (text=element_text(),
                   axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                   axis.text.y.right = element_text(size=12,color="#303030"),   
                   axis.text.y.left = element_text(size=12,color="#303030"),   
                   axis.title.x  = element_blank(), #margin is added to distance title from axis
                   axis.title.y.right = element_text(size=14,color="#303030", margin = margin(t = 0, r = 0, b = 0, l = 10)), #margin is added to distance title from axis
                   axis.title.y.left = element_text(size=14,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.ticks = element_line(colour = "#303030"),
                   panel.border = element_rect(color="#303030", fill= NA, size=0.6),
                   aspect.ratio = 0.5,
                   legend.position = "none")



#will reduce to half size, so I use double size of font
ggsave("Z_Sample size_trunk length.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 16, height = 10, units = "cm",
       dpi = 300, limitsize = TRUE) 



################  SUPPLEMENTARY FIGURE #########################


#LOCOMOTORY DISPARITY: Disparity and sample size (13 TIME BINS)
#--------------------------------------------------------------#

#1-PCA_My Morphospace_Disparity-PManova: calculate the disparity for all taxa divided into 11 time bins


#assigning numerical values to the ages
ages<- c(251, 242.1, 232, 214.1, 196, 182.4, 168.8, 160.4, 151.1, 135, 113, 93.4, 76.1)

#Creatig a column with the median (disparity obs)
Disparity <- c(1.5,6.03,8.03,2.63,8.3,12.18,14.99,17.09,11.5,15.07,20.14, 11.2,16)
Disparity


#Creating a column with the sample size
n<-c(6,27,16,7,10,18,7,5,14,5,5,7,14)
n

df2 <- data.frame(ages,Disparity,n) #creates a dataframe with various columns
df2

pca_plot <- ggplot(df2) +
  geom_line(aes(x=ages, y=Disparity)) +
  geom_line(aes(x=ages, y = n, colour = "Red"))+
  scale_x_reverse (limits=c(255, 66))+
  scale_y_continuous(sec.axis = sec_axis(~ .*1, name = "Sample size")) +#adds a secondary axis
  labs(x = "Time (Mya)",
       y = "Disparity (sum of variances)")+
  theme_light()

pca_plot

#Format plot:
pca_plot  + theme (text=element_text(),
                   axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                   axis.text.y.right = element_text(size=12,color="#303030"),   
                   axis.text.y.left = element_text(size=12,color="#303030"),   
                   axis.title.x  = element_blank(), #margin is added to distance title from axis
                   axis.title.y.right = element_text(size=14,color="#303030", margin = margin(t = 0, r = 0, b = 0, l = 10)), #margin is added to distance title from axis
                   axis.title.y.left = element_text(size=14,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.ticks = element_line(colour = "#303030"),
                   panel.border = element_rect(color="#303030", fill= NA, size=0.6),
                   aspect.ratio = 0.5,
                   legend.position = "none")



#will reduce to half size, so I use double size of font
ggsave("Z_Sample size_locomotory disparity 13 time bins.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 16, height = 10, units = "cm",
       dpi = 300, limitsize = TRUE) 





#TRUNK LENGTH: Standard deviation and sample size (13 TIME BINS)
#--------------------------------------------------------------#

dframe1<-read.csv("MMreptiles_Trunk.csv",   #MMreptiles_Trunk.csv with 8 time bins
                  row.names = 1,
                  header=T)

dframe1

#Before plotting FILTER THE DATA THAT HAVE <NA> in a given column
dframe2  <- dframe1[!is.na(dframe1$Time_4), ]
dframe2

#Calculate the SD
dframe3<-meanSD(TRUNK ~ A, data = dframe2)
dframe3
#dframe3<-t(dframe3) #transpose = switch columns and rows
#dframe3

#assigning numerical values to the ages
ages<- c(251, 242.1, 232, 214.1, 196, 182.4, 168.8, 160.4, 151.1, 135, 113, 93.4, 76.1)

#Creatig a column with the SD
SD <- c(18.70557, 65.27148, 83.56136, 27.64126, 30.60752,  54.7407,  57.1060,41.9, 71.599, 119.3, 117.1086, 
        83.55539, 106.7301)
SD

#Creating a column with the sample size
n<-c(10,37,26,11,14,24,7,5,23,11, 17, 13,24)
n

df <- data.frame(ages,SD,n) #creates a dataframe with various columns
df

pca_plot <- ggplot(df) +
  geom_line(aes(x=ages, y=SD)) +
  geom_line(aes(x=ages, y = n, colour = "Red"))+
  scale_x_reverse (limits=c(255, 66))+
  scale_y_continuous(sec.axis = sec_axis(~ .*1, name = "Sample size")) +#adds a secondary axis
  labs(x = "Time (Mya)",
       y = "Trunk length dispersion (sd)")+
  theme_light()

pca_plot


#Format plot:
pca_plot  + theme (text=element_text(),
                   axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                   axis.text.y.right = element_text(size=12,color="#303030"),   
                   axis.text.y.left = element_text(size=12,color="#303030"),   
                   axis.title.x  = element_blank(), #margin is added to distance title from axis
                   axis.title.y.right = element_text(size=14,color="#303030", margin = margin(t = 0, r = 0, b = 0, l = 10)), #margin is added to distance title from axis
                   axis.title.y.left = element_text(size=14,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.ticks = element_line(colour = "#303030"),
                   panel.border = element_rect(color="#303030", fill= NA, size=0.6),
                   aspect.ratio = 0.5,
                   legend.position = "none")



#will reduce to half size, so I use double size of font
ggsave("Z_Sample size_trunk length 13 time bins.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 16, height = 10, units = "cm",
       dpi = 300, limitsize = TRUE) 





