#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
#http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining

#LOAD PACKAGES
library(FactoMineR) #R package for computing PCA
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



rm(list = ls())
#setwd("######")

# OPEN FILE
dframe1<-read.csv("FossilAqReptiles_data.csv", 
                  #row.names = 1,
                  header=T)
dframe1

#Exclude incomplete rows
dframe1<-na.exclude(dframe1)
dframe1

# Remove duplicates based on "Taxon" column
#dframe1<-dframe1[!duplicated(dframe1$Taxon),]

## Use the first column for row names
dframe1 <- data.frame(dframe1, row.names = 1)
dframe1

#This can be used to remove specific columns!
#dframe1 <- subset(dframe1, select = -LOWERJAW_ratio)
#dframe1


# VISUALISING PAIRS OF VARIABLES / COVARIATION
#----------------------------------------------------------#


cor <- cor(dframe1[1:16])   # request the matrix of correlation coefficients

#write.csv(cor(dframe1[1:16]), "Correlation_Variables.csv") # save it as a csv file
 

pairs (dframe1[1:16],  # Pairwise plots of numbered variables
       panel = panel.smooth)         #  add a wiggly fitted line


# install.packages("car")
library(car)
scatterplotMatrix(dframe1[c(1:16)])       # A matrix of scatterplots with frequency distributions
# Note that smoothers can't be fitted for some combinations here


#Plotting the CORRELATION of each variable
#----------------------------------------#

library(corrplot)

#Simple corr-plot
corrplot(cor, 
         method="color",
         tl.col="black")


# Computing the p-value of the correlation
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(dframe1[1:16])
head(p.mat[, 1:5])

#Correlation plot including significance 
corrplot(cor, 
         method="color",
         tl.col="black",
         tl.cex=0.8,
         cl.cex=0.8,
         p.mat = p.mat, 
         sig.level = 0.05,
         insig = "blank")

#Save 7x7 inches




# PERFORM PCA 
#-----------#

#ON ALL CHARACTER DATA - and selects the relevant columns from the data (first column = taxon doesn't count)
res.pca <- PCA(dframe1[,1:16], scale.unit = TRUE, graph=TRUE)


#1- OMIT TRUNK
res.pca <- PCA(dframe1[,2:16], scale.unit = TRUE)

#2-OMIT TRUNK AND NECK
res.pca <- PCA(dframe1[,3:16], scale.unit = TRUE)


#3-OMIT TRUNK, NECK AND TAIL
res.pca <- PCA(dframe1[,4:16], scale.unit = TRUE)


#4-OMIT TRUNK, NECK, TAIL AND JAW
res.pca <- PCA(dframe1[,5:16], scale.unit = TRUE)

#5-OMIT TRUNK, NECK, TAIL, JAW, BRACHIAL 
res.pca <- PCA(dframe1[,6:16], scale.unit = TRUE)

#6-OMIT TRUNK, NECK, TAIL, JAW, BRACHIAL, CRURAL
res.pca <- PCA(dframe1[,7:16], scale.unit = TRUE)


#7-OMIT TRUNK, NECK, TAIL, JAW, BRACHIAL, CRURAL, FORE-HIND
res.pca <- PCA(dframe1[,8:16], scale.unit = TRUE)


#8-OMIT TRUNK, NECK, TAIL, JAW, BRACHIAL, CRURAL, FORE-HIND, HUMFLARE
res.pca <- PCA(dframe1[,9:16], scale.unit = TRUE)


#9-OMIT TRUNK, NECK, TAIL, JAW, BRACHIAL, CRURAL,  FORE-HIND, HUMFLARE, FEMFLARE
res.pca <- PCA(dframe1[,10:16], scale.unit = TRUE)


#10-OMIT TRUNK, NECK, TAIL, JAW, BRACHIAL, CRURAL,FORE-HIND, HUMFLARE, FEMFLARE, FORERAT
res.pca <- PCA(dframe1[,11:16], scale.unit = TRUE)

#11-OMIT TRUNK, NECK, TAIL, JAW, BRACHIAL, CRURAL, FORE-HIND, HUMFLARE, FEMFLARE, FORERAT, HINDRAT
res.pca <- PCA(dframe1[,12:16], scale.unit = TRUE)

#12-OMIT TRUNK, NECK, TAIL, JAW, BRACHIAL, CRURAL, FORE-HIND, HUMFLARE, FEMFLARE, FORERAT, HINDRAT - MANUS_PES
res.pca <- PCA(dframe1[,12:15], scale.unit = TRUE)

#13-OMIT TRUNK, NECK, TAIL, JAW, BRACHIAL, CRURAL, FORE-HIND, HUMFLARE, FEMFLARE, FORERAT, HINDRAT, MANUSRAT
res.pca <- PCA(dframe1[,13:16], scale.unit = TRUE)


#14-OMIT TRUNK, NECK, TAIL, JAW, BRACHIAL, CRURAL, FORE-HIND, HUMFLARE, FEMFLARE, FORERAT, HINDRAT
res.pca <- PCA(dframe1[,14:15], scale.unit = TRUE)


#15-NO HINDLIMB
dframe1 <- subset(dframe1, select = -LOWERJAW_ratio)  #remove specific columns
dframe1

res.pca <- PCA(dframe1[,1:15], scale.unit = TRUE)


#16-NO TANYSTROPHEIDS
dframe1<- subset(dframe1,Group_1 %in% c("Ichthyosauromorpha", "Mosasauroidea", "Pantestudines", "Rhynchocephalia", "Sauropterygia","Thalattosauria", "Thalattosuchia", "Saurosphargidae"))
dframe1

res.pca <- PCA(dframe1[,1:16], scale.unit = TRUE, graph=TRUE)



#############################

# EIGENVALUES
#------------#

eigenvalues <- res.pca$eig #Plot showing the contribution of each axis to the overall variation
head(eigenvalues[, 1:2])

barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        xlab = "Principal Components",
        ylab = "Percentage of variance",
        col ="gray")

eigenvalues

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 30)) #Ploting eigenvalues with % (SCREE PLOT)

get_eig(res.pca) #Another way for extracting the eigenvalues


#############################

#Extract the coordinates for the variables
#----------------------------------------#
varcoord <- res.pca$var$coord
varcoord


#Extract the contribution of each variable
#----------------------------------------#
varcontrib <- res.pca$var$contrib
varcontrib


# Extract the results for variables
#---------------------------------#
var <- get_pca_var(res.pca)
var


#############################

# VARIABLE CONTRIBUTION: plots showing CHARACTERS AFFECTING THE PCA AXES 
#------------------------------------------------------------------------#

#http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining

#Basic plot of variable contribution
#----------------------------------#

plot(res.pca, choix = "var", axes=c(1,2))
plot(res.pca, choix = "var", axes=c(3,2))

#Variable contribution from factoextra (FIGURE)
#--------------------------------------------#

#PCA1-2
p<-fviz_pca_var(res.pca, 
                col.var="contrib", 
                labelsize = 4, repel = TRUE,
                gradient.cols = c("lightgrey","#c20000")) + # By default, variables/individuals are represented on dimensions 1 and 2
                xlab("PCA 1") + ylab("PCA 2") + ggtitle("") +
                labs(color='Contribution') 
p 

p + theme (text=element_text(),
           axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
           axis.text.y = element_text(size=12,color="#303030"),
           axis.title.x  = element_text(size=14,color="#303030", margin = margin(t = 10, r = 0, b = 0, l = 0)), #margin is added to distance title from axis
           axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 1, b = 0, l = 0)), #margin is added to distance title from axis
           panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank(),
           axis.ticks = element_line(colour = "#303030"),
           panel.border = element_rect(color="#303030", fill= NA, size=0.6),
           aspect.ratio = 1,
           legend.position= c(1.11, 0.8),
           legend.title = element_blank())

#legend.title = element_text( c=("legend")))
#legend.position = "none")

#will reduce to half size, so I use double size of font
ggsave("PCA 1-2 contribution.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, height = 14, units = "cm",
       dpi = 300, limitsize = TRUE) 




#PCA3-2
p<-fviz_pca_var(res.pca, 
                col.var="contrib", axes = c(3, 2),
                labelsize = 4, repel = TRUE,
                gradient.cols = c("lightgrey","#c20000")) + # By default, variables/individuals are represented on dimensions 1 and 2
                xlab("PCA 3") + ylab("PCA 2") + ggtitle("")  +
                labs(color='Contribution') 
p 


p + theme (text=element_text(),
           axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
           axis.text.y = element_text(size=12,color="#303030"),
           axis.title.x  = element_text(size=14,color="#303030", margin = margin(t = 10, r = 0, b = 0, l = 0)), #margin is added to distance title from axis
           axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 1, b = 0, l = 0)), #margin is added to distance title from axis
           panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank(),
           axis.ticks = element_line(colour = "#303030"),
           panel.border = element_rect(color="#303030", fill= NA, size=0.6),
           aspect.ratio = 1,
           legend.position= c(1.11, 0.8),
           legend.title = element_blank())



#will reduce to half size, so I use double size of font
ggsave("PCA 3-2 contribution.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, height = 14, units = "cm",
       dpi = 300, limitsize = TRUE) 



#PCA4-5
p<-fviz_pca_var(res.pca, 
                col.var="contrib", axes = c(5, 4),
                gradient.cols = c("lightgrey","#c20000")) + # By default, variables/individuals are represented on dimensions 1 and 2
                xlab("PCA 5") + ylab("PCA 4") + ggtitle("")   +
                labs(color='Contribution') 
p 

p + theme (text=element_text(),
           axis.text.x = element_text(size=20, color="#303030", hjust=0.5),
           axis.text.y = element_text(size=20,color="#303030"),
           axis.title.x  = element_text(size=24,color="#303030", margin = margin(t = 10, r = 0, b = 0, l = 0)), #margin is added to distance title from axis
           axis.title.y  = element_text(size=24,color="#303030", margin = margin(t = 0, r = 1, b = 0, l = 0)), #margin is added to distance title from axis
           panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank(),
           axis.ticks = element_line(colour = "#303030"),
           panel.border = element_rect(color="#303030", fill= NA, size=0.6),
           aspect.ratio = 1,
           legend.position= c(0.1, 0.86))


ggsave("PCA 5-4 contribution.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 18, height = 18, units = "cm",
       dpi = 300, limitsize = TRUE) 



#Plotting the CORRELATION of each variable contribution to the PCAs
#-----------------------------------------------------------------#

library(corrplot)

corrplot(var$contrib, is.corr=FALSE)
corrplot(var$cos2, is.corr=FALSE)

#Correlation for Dimensions 1-5 (FIGURE)
corrplot(var$contrib, is.corr=FALSE,
         method = "color",
         tl.col="black", 
         tl.cex=1.2,
         cl.cex=0.5,
         cl.pos = "b",
         cl.lim = c(0, 43),
         col= colorRampPalette(c("grey","white","#c20000"))(100) )
         #tl.cex=seq(.5, 2, length.out=ncol(var$contrib)))


#Save 7x9 inches


# Show contributions of each variable to the principal components
#---------------------------------------------------------------#
head(var$contrib, 16) #the number indicates how many variables will appear in the graph


# Identify the MOST SIGNIFICANTLY ASSOCIATED VARIABLES with a given PCA. Uses dimension description function dimdesc()
#--------------------------------------------------------------------#
res.desc <- dimdesc(res.pca, axes = c(1,2), proba = 0.05)
res.desc

# Description of dimension 1
res.desc$Dim.1
# Description of dimension 1
res.desc$Dim.2


# Barplots of variable contributions
#---------------------------------#
fviz_contrib(res.pca, choice = "var", axes = 1, top = 16) #variable contribution to PC1
fviz_contrib(res.pca, choice = "var", axes = 2, top = 16) #variable contribution to PC2
fviz_contrib(res.pca, choice = "var", axes = 3, top = 16) #variable contribution to PC3


############################
############################


# INDIVIDUALS CONTRIBUTIONS
#-------------------------#

ind <- get_pca_ind(res.pca)
ind

# Extract oordinates of individuals: normal PCA
head(ind$coord)
# Extract Cos 2: Quality of individuals: how much a given PCA (the importance of a principal component for a given point)
head(ind$cos2)
# Contrib: Contributions of individuals to each PCA
head(ind$contrib, 10) #the number indicates how many individuals will be shown in the table


#Simple plot of individuals
fviz_pca_ind(res.pca)

#Plot individual contributions coloring and controlling size of dot according to cos2 or contrib
#---------------------------------------------------------------------------------------------#
fviz_pca_ind(res.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

fviz_pca_ind(res.pca, col.ind = "contrib",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
repel = TRUE # Avoid text overlapping (slow if many points)
)

fviz_pca_ind(res.pca, col.ind = "cos2", geom = "point")+
  scale_color_gradient2(low="white", mid="blue",high="red", midpoint=0.5)+
  theme_minimal()


fviz_pca_ind(res.pca, col.ind = "contrib", geom = "point")+
            scale_color_gradient2(low="white", mid="blue",high="red", midpoint=1.5)+
            theme_minimal()

 
fviz_pca_ind(res.pca, geom = "point", pointsize = "contrib", 
              pointshape = 21, fill = "#E7B800") +
              theme_minimal()


# Barplot of the Total contribution
#---------------------------------#
fviz_contrib(res.pca, choice = "ind", axes = 1:2)# on PC1 and PC2
fviz_contrib(res.pca, choice = "ind", axes = 3:2)# on PC3 and PC2


#############################

# VISUALISE MORPHOSPACE: MAKE A SIMPLE MORPHOSPACE PLOT OF AXES 1 TO 4
plot(res.pca, choix = "ind", cex=0.5, axes=c(1,2))
plot(res.pca, choix = "ind", cex=0.5, axes=c(2,3))
plot(res.pca, choix = "ind", cex=0.5, axes=c(3,4))


#############################

# RENAME THE PC SCORES (MORPHOSPACE COORDINATES)
pc.scores <- res.pca$ind$coord # shows the coordinates for each individual  
pc.scores

# pc.scores <- read.table("pc.scores.final.txt", row.names=1, header=T) #
pc.scores <- as.matrix(pc.scores) 
pc.scores # shows the morphospace location coordinates of each taxa 

#Save these scores in a separate file, that will be used later to create the file: pcscores-dup.csv, 
#duplicating some taxa (Taxon / Taxon_d) that spans various time bins.


###############################

# MAKE PRETTY PLOT
#----------------#

# create lists of taxa grouped per CLADE (Column 18)
Sauropterygia_list <- row.names(dframe1)[dframe1[,17] == 'Sauropterygia']  # change names of groups to match your data
Sauropterygia_list

Thalattosauria_list <- row.names(dframe1)[dframe1[,17] == 'Thalattosauria']  # change names of groups to match your data
Thalattosauria_list

Saurosphargidae_list <- row.names(dframe1)[dframe1[,17] == 'Saurosphargidae']  # change names of groups to match your data
Saurosphargidae_list

Tanystropheidae_list <- row.names(dframe1)[dframe1[,17] == 'Tanystropheidae']  # change names of groups to match your data
Tanystropheidae_list

Thalattosuchia_list <- row.names(dframe1)[dframe1[,17] == 'Thalattosuchia']  # change names of groups to match your data
Thalattosuchia_list

Mosasauroidea_list <- row.names(dframe1)[dframe1[,17] == 'Mosasauroidea']  # change names of groups to match your data
Mosasauroidea_list

Ichthyosauromorpha_list <- row.names(dframe1)[dframe1[,17] == 'Ichthyosauromorpha']  # change names of groups to match your data
Ichthyosauromorpha_list

Pantestudines_list <- row.names(dframe1)[dframe1[,17] == 'Pantestudines']  # change names of groups to match your data
Pantestudines_list

Rhynchocephalia_list <- row.names(dframe1)[dframe1[,17] == 'Rhynchocephalia']  # change names of groups to match your data
Rhynchocephalia_list


# PCA 1-2
# MAKE EMPTY PLOT PC1 vs PC2

plot.new()
plot(pc.scores[,1], pc.scores[,2], pch=19, lwd=2, col="transparent", xlab=paste("PC1", " (", round(eigenvalues[, 2][1], 2), "%)", sep=""), ylab=paste("PC2", " (", round(eigenvalues[, 2][2], 2), "%)", sep=""))
#abline(v = -3, col = "gray99", lwd = 10000) # make gray background
abline(v = -3, col = "white", lwd = 10000) # make white background
abline(v=c(0), lwd=1, col="gray80", lty=2)
abline(h=c(0), lwd=1, col="gray80",lty=2) # plot centre of space 
box(lwd=1) # add box 
par(new=TRUE) # this creates a new layer - so you can keep on plotting in same space


# DRAW HULLS

Plot_ConvexHull<-function(xcoord, ycoord, lcolor, bgcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
  polygon(x = xcoord[hpts], y =  ycoord[hpts], col = adjustcolor(bgcolor, alpha.f = 0.27) , border = 0)
  
}  

# HULLS PER CLADE

Plot_ConvexHull(xcoord=pc.scores[Sauropterygia_list,][,1], ycoord=pc.scores[Sauropterygia_list,][,2], lcolor="#19D27E", bgcolor="#19D27E")
Plot_ConvexHull(xcoord=pc.scores[Thalattosauria_list,][,1], ycoord=pc.scores[Thalattosauria_list,][,2], lcolor="#FFD300", bgcolor="#FFD300")
Plot_ConvexHull(xcoord=pc.scores[Saurosphargidae_list,][,1], ycoord=pc.scores[Saurosphargidae_list,][,2], lcolor="#003366", bgcolor="darkgrey")
Plot_ConvexHull(xcoord=pc.scores[Thalattosuchia_list,][,1], ycoord=pc.scores[Thalattosuchia_list,][,2], lcolor="blue", bgcolor="blue")
Plot_ConvexHull(xcoord=pc.scores[Tanystropheidae_list,][,1], ycoord=pc.scores[Tanystropheidae_list,][,2], lcolor="#3AE1FF", bgcolor="#3AE1FF")
Plot_ConvexHull(xcoord=pc.scores[Mosasauroidea_list,][,1], ycoord=pc.scores[Mosasauroidea_list,][,2], lcolor="#FF8300", bgcolor="#FF8300")
Plot_ConvexHull(xcoord=pc.scores[Ichthyosauromorpha_list,][,1], ycoord=pc.scores[Ichthyosauromorpha_list,][,2], lcolor="#FF2626", bgcolor="#FF2626")
Plot_ConvexHull(xcoord=pc.scores[Pantestudines_list,][,1], ycoord=pc.scores[Pantestudines_list,][,2], lcolor="#009ccc", bgcolor="#009ccc") 
Plot_ConvexHull(xcoord=pc.scores[Rhynchocephalia_list,][,1], ycoord=pc.scores[Rhynchocephalia_list,][,2], lcolor="Purple", bgcolor="#c061ff")


# DRAW POINTS PER CLADE

points(pc.scores[Sauropterygia_list,][,1], pc.scores[Sauropterygia_list,][,2], pch=21, col="black", bg="#19D27E", cex=1)
points(pc.scores[Thalattosauria_list,][,1], pc.scores[Thalattosauria_list,][,2], pch=25, col="black", bg="#FFD300", cex=1)
points(pc.scores[Saurosphargidae_list,][,1],pc.scores[Saurosphargidae_list,][,2], pch=21, col="black", bg="darkgrey", cex=1)
points(pc.scores[Thalattosuchia_list,][,1], pc.scores[Thalattosuchia_list,][,2], pch=22, col="black", bg="blue", cex=1)
points(pc.scores[Tanystropheidae_list,][,1], pc.scores[Tanystropheidae_list,][,2], pch=24, col="black", bg="#c4f4ff", cex=1)
points(pc.scores[Mosasauroidea_list,][,1], pc.scores[Mosasauroidea_list,][,2], pch=23, col="black", bg="#FF8300", cex=1)
points(pc.scores[Ichthyosauromorpha_list,][,1],pc.scores[Ichthyosauromorpha_list,][,2], pch=24, col="black", bg="#FF2626", cex=1)
points(pc.scores[Pantestudines_list,][,1], pc.scores[Pantestudines_list,][,2], pch=23, col="black", bg="#0083ab", cex=1)
points(pc.scores[Rhynchocephalia_list,][,1],pc.scores[Rhynchocephalia_list,][,2], pch=25, col="black", bg="#c061ff", cex=1)


#for presentation save in 6.2 x 6.2


# MAKE SIMPLE PLOT WITH LABELS
text(pc.scores[,1], pc.scores[,2], labels=row.names(pc.scores),pos=4,cex=0.5,col="#565656")




# PCA 2-3
# MAKE EMPTY PLOT- PC2 vs PC3

plot.new()
plot( pc.scores[,3], pc.scores[,2], pch=19, lwd=2, col="transparent", xlab=paste("PC3", " (", round(eigenvalues[, 2][3], 2), "%)", sep=""), ylab=paste("PC2", " (", round(eigenvalues[, 2][2], 2), "%)", sep=""))
#abline(v = -3, col = "gray99", lwd = 10000) # make gray background
abline(v = -3, col = "white", lwd = 10000) # make white background
abline(v=c(0), lwd=1, col="gray80", lty=2)
abline(h=c(0), lwd=1, col="gray80",lty=2) # plot centre of space 
box(lwd=1) # add box 
par(new=TRUE) # this creates a new layer - so you can keep on plotting in same space

# DRAW HULLS

Plot_ConvexHull<-function(xcoord, ycoord, lcolor, bgcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
  polygon(x = xcoord[hpts], y =  ycoord[hpts], col = adjustcolor(bgcolor, alpha.f = 0.27) , border = 0)
  
}  


# HULLS PER CLADE

Plot_ConvexHull(xcoord=pc.scores[Sauropterygia_list,][,3], ycoord=pc.scores[Sauropterygia_list,][,2], lcolor="#19D27E", bgcolor="#19D27E")
Plot_ConvexHull(xcoord=pc.scores[Thalattosauria_list,][,3], ycoord=pc.scores[Thalattosauria_list,][,2], lcolor="#FFD300", bgcolor="#FFD300")
Plot_ConvexHull(xcoord=pc.scores[Saurosphargidae_list,][,3], ycoord=pc.scores[Saurosphargidae_list,][,2], lcolor="#003366", bgcolor="darkgrey")
Plot_ConvexHull(xcoord=pc.scores[Thalattosuchia_list,][,3], ycoord=pc.scores[Thalattosuchia_list,][,2], lcolor="blue", bgcolor="blue")
Plot_ConvexHull(xcoord=pc.scores[Tanystropheidae_list,][,3], ycoord=pc.scores[Tanystropheidae_list,][,2], lcolor="#3AE1FF", bgcolor="#3AE1FF")
Plot_ConvexHull(xcoord=pc.scores[Mosasauroidea_list,][,3], ycoord=pc.scores[Mosasauroidea_list,][,2], lcolor="#FF8300", bgcolor="#FF8300")
Plot_ConvexHull(xcoord=pc.scores[Ichthyosauromorpha_list,][,3], ycoord=pc.scores[Ichthyosauromorpha_list,][,2], lcolor="#FF2626", bgcolor="#FF2626")
Plot_ConvexHull(xcoord=pc.scores[Pantestudines_list,][,3], ycoord=pc.scores[Pantestudines_list,][,2], lcolor="#009ccc", bgcolor="#009ccc") 
Plot_ConvexHull(xcoord=pc.scores[Rhynchocephalia_list,][,3], ycoord=pc.scores[Rhynchocephalia_list,][,2], lcolor="Purple", bgcolor="#c061ff")


# DRAW POINTS PER CLADE

points(pc.scores[Sauropterygia_list,][,3], pc.scores[Sauropterygia_list,][,2], pch=21, col="black", bg="#19D27E", cex=1)
points(pc.scores[Thalattosauria_list,][,3], pc.scores[Thalattosauria_list,][,2], pch=25, col="black", bg="#FFD300", cex=1)
points(pc.scores[Saurosphargidae_list,][,3],pc.scores[Saurosphargidae_list,][,2], pch=21, col="black", bg="darkgrey", cex=1)
points(pc.scores[Thalattosuchia_list,][,3], pc.scores[Thalattosuchia_list,][,2], pch=22, col="black", bg="blue", cex=1)
points(pc.scores[Tanystropheidae_list,][,3], pc.scores[Tanystropheidae_list,][,2], pch=24, col="black", bg="#c4f4ff", cex=1)
points(pc.scores[Mosasauroidea_list,][,3], pc.scores[Mosasauroidea_list,][,2], pch=23, col="black", bg="#FF8300", cex=1)
points(pc.scores[Ichthyosauromorpha_list,][,3],pc.scores[Ichthyosauromorpha_list,][,2], pch=24, col="black", bg="#FF2626", cex=1)
points(pc.scores[Pantestudines_list,][,3], pc.scores[Pantestudines_list,][,2], pch=23, col="black", bg="#0083ab", cex=1)
points(pc.scores[Rhynchocephalia_list,][,3],pc.scores[Rhynchocephalia_list,][,2], pch=25, col="black", bg="#c061ff", cex=1)


#for presentation save in 6.2 x 6.2


# MAKE SIMPLE PLOT WITH LABELS
text(pc.scores[,3], pc.scores[,2], labels=row.names(pc.scores),pos=4,cex=0.5,col="#565656")



########################
########################
#Trunk
pca_plot <-ggplot(dframe1, aes(x=Group_1, y=Trunk, fill=Group_1)) +  # Change for each variable
  geom_boxplot(width=0.5,
               lwd=0.3,
               outlier.alpha = 0.2,
               alpha=0.4) + 
  geom_jitter(width = 0.2, alpha=0.25)+
  scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","#3944bc"))+
  labs(x = "",
       y = "Trunk (cm)")+
  theme_light()

pca_plot

#Lower jaw
pca_plot <-ggplot(dframe1, aes(x=Group_1, y=Lower.jaw.ratio, fill=Group_1)) +  # Change for each variable
  geom_boxplot(width=0.5,
               lwd=0.3,
               outlier.alpha = 0.2,
               alpha=0.4) + 
  geom_jitter(width = 0.2, alpha=0.25)+
  scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","#3944bc"))+
  labs(x = "",
       y = "Lower jaw ratio")+
  theme_light()

pca_plot


#Neck ratio
pca_plot <-ggplot(dframe1, aes(x=Group_1, y=Neck.ratio, fill=Group_1)) +  # Change for each variable
  geom_boxplot(width=0.5,
               lwd=0.3,
               outlier.alpha = 0.2,
               alpha=0.4) + 
  geom_jitter(width = 0.2, alpha=0.25)+
  scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","#3944bc"))+
  labs(x = "",
       y = "Neck ratio")+
  theme_light()

pca_plot


#Tail ratio
pca_plot <-ggplot(dframe1, aes(x=Group_1, y=Tail.ratio, fill=Group_1)) +  # Change for each variable
  geom_boxplot(width=0.5,
               lwd=0.3,
               outlier.alpha = 0.2,
               alpha=0.4) + 
  geom_jitter(width = 0.2, alpha=0.25)+
  scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","#3944bc"))+
  labs(x = "",
       y = "Tail ratio")+
  theme_light()

pca_plot

#Tail ratio
pca_plot <-ggplot(dframe1, aes(x=Group_1, y=Tail.ratio, fill=Group_1)) +  # Change for each variable
  geom_boxplot(width=0.5,
               lwd=0.3,
               outlier.alpha = 0.2,
               alpha=0.4) + 
  geom_jitter(width = 0.2, alpha=0.25)+
  scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","#3944bc"))+
  labs(x = "",
       y = "Tail ratio")+
  theme_light()

pca_plot

#Brachial
pca_plot <-ggplot(dframe1, aes(x=Group_1, y=Brachial.ind., fill=Group_1)) +  # Change for each variable
  geom_boxplot(width=0.5,
               lwd=0.3,
               outlier.alpha = 0.2,
               alpha=0.4) + 
  geom_jitter(width = 0.2, alpha=0.25)+
  scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","#3944bc"))+
  labs(x = "",
       y = "Brachial index")+
  theme_light()

pca_plot

#Crural
pca_plot <-ggplot(dframe1, aes(x=Group_1, y=Crural.ind., fill=Group_1)) +  # Change for each variable
  geom_boxplot(width=0.5,
               lwd=0.3,
               outlier.alpha = 0.2,
               alpha=0.4) + 
  geom_jitter(width = 0.2, alpha=0.25)+
  scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","#3944bc"))+
  labs(x = "",
       y = "Crural index")+
  theme_light()

pca_plot


#Fore/hind ratio
pca_plot <-ggplot(dframe1, aes(x=Group_1, y=Fore.hind.ratio, fill=Group_1)) +  # Change for each variable
  geom_boxplot(width=0.5,
               lwd=0.3,
               outlier.alpha = 0.2,
               alpha=0.4) + 
  geom_jitter(width = 0.2, alpha=0.25)+
  scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","#3944bc"))+
  labs(x = "",
       y = "Fore/hind ratio")+
  theme_light()

pca_plot



#Format plot:
pca_plot  + theme (text=element_text(),
                   axis.text.x = element_text(size=12, color="#303030", angle= 45, hjust=1),
                   axis.text.y = element_text(size=12,color="#303030"),
                   axis.title.x  = element_blank(),
                   axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 10, b = 0, l = 0)), #margin is added to distance title from axis
                   strip.text = element_text(size=12,color="#303030"), #Change size of facet title
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.ticks = element_line(colour = "#303030"),
                   panel.border = element_rect(color="#303030", fill= NA, size=0.5),
                   aspect.ratio = 0.5,
                   legend.position = "none")


#will reduce to half size, so I use double size of font
ggsave("Variables_Forehind-per group.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 14, height = 10, units = "cm",
       dpi = 300, limitsize = TRUE) 

#This takes forever...


#######################
######################

rm(list = ls())
#setwd("######")

# OPEN FILE
dframe1<-read.csv("FossilAqReptiles_data.csv", 
                  #row.names = 1,
                  header=T)
dframe1

#Exclude incomplete rows
dframe1<-na.exclude(dframe1)
dframe1

# Remove duplicates based on "Taxon" column
#dframe1<-dframe1[!duplicated(dframe1$Taxon),]

## Use the first column for row names
dframe1 <- data.frame(dframe1, row.names = 1)
dframe1



#Make plot with dots and polygons
pca_plot <- ggplot(dframe1,aes(x = Dim.1, y = Dim.2, size=Trunk, fill = Group_1, colour=Group_1, shape=Group_1
)) + #I need to specify color here so that ggMarginal can take it
  geom_point() +
  #geom_polygon(data = hull_cyl, alpha = 0.2, linetype=0)+      #0 is the size of the polygon line (before I had size=0)
  scale_colour_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","blue"))+
  scale_fill_manual(values= c("#FF2626","#FF8300","#009ccc","#c061ff","#19D27E","darkgrey","#3AE1FF", "#FFD300","blue"))+
  scale_shape_manual(values= c(24,23,21,25,21,24,21,25,22))+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #horizontal reference line in O
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #vertical reference line in O
  labs(x = "PC1",
       y = "PC2")+
  theme_light()

pca_plot

pca_plot <- pca_plot + geom_point(color="black",stroke=0.3) #add points in layers! otherwise it messes up the colors of histograms, stroke indicates thickness of line

pca_plot


#Format plot:
pca_plot <- pca_plot  + theme (text=element_text(),
                               axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                               axis.text.y = element_text(size=12,color="#303030"),
                               axis.title.x  = element_text(size=14,color="#303030", margin = margin(t = 8, r = 0, b = 0, l = 0)), #margin is added to distance title from axis
                               axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 8, b = 0, l = 0)), #margin is added to distance title from axis
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
ggsave("Variables PCA 1-2 PES AR.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 9, height = 9, units = "cm",
       dpi = 300, limitsize = TRUE) 



#To save LEGEND


#Format plot:
pca_plot  + theme (text=element_text(),
                               axis.text.x = element_text(size=12, color="#303030", hjust=0.5),
                               axis.text.y = element_text(size=12,color="#303030"),
                               axis.title.x  = element_text(size=14,color="#303030", margin = margin(t = 8, r = 0, b = 0, l = 0)), #margin is added to distance title from axis
                               axis.title.y  = element_text(size=14,color="#303030", margin = margin(t = 0, r = 8, b = 0, l = 0)), #margin is added to distance title from axis
                               strip.text = element_text(size=12,color="#303030"), #Change size of facet title
                               panel.grid.major.x = element_blank(),
                               panel.grid.minor.x = element_blank(),
                               panel.grid.major.y = element_blank(),
                               panel.grid.minor.y = element_blank(),
                               axis.ticks = element_line(colour = "#303030"),
                               panel.border = element_rect(color="#303030", fill= NA, size=0.6),
                               aspect.ratio = 1,
                               legend.justification=c(1,0), 
                               legend.position=c(0.95, 0.05),  #legend.position=c(0.95, 0.05),
                               legend.background = element_blank())



#I will reduce to 1/2 size, so I use font x2
ggsave("Variables PCA 1-2 TRUNK-legend.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 9, height = 9, units = "cm",
       dpi = 300, limitsize = TRUE) 
 
