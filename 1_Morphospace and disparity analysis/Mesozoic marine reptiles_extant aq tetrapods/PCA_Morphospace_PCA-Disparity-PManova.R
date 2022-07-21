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
dframe1<-read.csv("FossilandExtantAqTet_data.csv", 
                  row.names = 1,
                  header=T)


dframe1

#Exclude incomplete rows: ONLY DO THIS FOR THE ANALYSIS WITH FULL-DATASET
dframe1<-na.exclude(dframe1)
dframe1


#This can be used to remove specific columns!
#dframe1 <- subset(dframe1, select = -LOWERJAW_ratio)
#dframe1

################################

# VISUALISING THE DATA: PAIRS OF VARIABLES

cor(dframe1[1:16])   # request the matrix of correlation coefficients

write.csv(cor(dframe1[1:16]), "Correlation_Variables.csv") # save it as a csv file


pairs (dframe1[1:16],  # Pairwise plots of numbered variables
       panel = panel.smooth)         #  add a wiggly fitted line



#1. PCA WITH ALL RATIOS (no cetaceans)
#---------------------------------#


# PERFORM PCA ON CHARACTER DATA - and selects the relevant columns from the data (first column = taxon doesn't count)
res.pca <- PCA(dframe1[,1:16], scale.unit = TRUE)



# THIS CODE MAKES A PLOT SHOWING THE CONTRIBUTION OF EACH PC AXIS TO OVERALL VARIATION
eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        xlab = "Principal Components",
        ylab = "Percentage of variance",
        col ="gray")

eigenvalues

write.csv((eigenvalues), "Eigenvalues.csv") # save it as a csv file



#Variables contributions
fviz_pca_var(res.pca, col.var="contrib") # By default, variables/individuals are represented on dimensions 1 and 2
fviz_pca_var(res.pca, col.var="contrib", axes = c(3, 2)) # specify different axis
fviz_pca_var(res.pca, col.var="contrib", axes = c(3, 4)) # specify different axis


# Simple PCA plot with labels
plot(res.pca, choix = "ind", cex=0.5, axes=c(1,2))
plot(res.pca, choix = "ind", cex=0.5, axes=c(2,3))


#############################

# RENAME THE PC SCORES (MORPHOSPACE COORDINATES)
pc.scores <- res.pca$ind$coord
pc.scores

# pc.scores <- read.table("pc.scores.final.txt", row.names=1, header=T) #
pc.scores <- as.matrix(pc.scores) 
pc.scores # shows the morphospace location coordinates of each taxa 


###############################

#FOR THE DISPARITY ANALYSIS:
#Save these scores in a separate file, that will be used later to create the file: pcscores.csv,
#write.csv(pc.scores, "pcscores.csv")
#then duplicate the taxa that spans various time bins: Taxon / Taxon_d


###############################

#FOR THE PERMANOVA:
# to make a csv file of the pcscores (the PC coordinates of each taxa) to which you can add the clade and time-bin group in new columns
#write.csv(pc.scores, "FossilAqreptiles_Output.csv")
# then add the group and time columns (I did this manually)


##############################
##############################

# MAKE PRETTY PLOT

# create lists of extant taxa grouped per LOCOMOTION and extinct taxa grouped per GROUP (Column 20)
Sauropterygia_list <- row.names(dframe1)[dframe1[,20] == 'Sauropterygia']  # change names of groups to match your data
Sauropterygia_list

Thalattosauria_list <- row.names(dframe1)[dframe1[,20] == 'Thalattosauria']  # change names of groups to match your data
Thalattosauria_list

Saurosphargidae_list <- row.names(dframe1)[dframe1[,20] == 'Saurosphargidae']  # change names of groups to match your data
Saurosphargidae_list

Tanystropheidae_list <- row.names(dframe1)[dframe1[,20] == 'Tanystropheidae']  # change names of groups to match your data
Tanystropheidae_list

Thalattosuchia_list <- row.names(dframe1)[dframe1[,20] == 'Thalattosuchia']  # change names of groups to match your data
Thalattosuchia_list

Mosasauroidea_list <- row.names(dframe1)[dframe1[,20] == 'Mosasauroidea']  # change names of groups to match your data
Mosasauroidea_list

Ichthyosauromorpha_list <- row.names(dframe1)[dframe1[,20] == 'Ichthyosauromorpha']  # change names of groups to match your data
Ichthyosauromorpha_list

Pantestudines_list <- row.names(dframe1)[dframe1[,20] == 'Pantestudines']  # change names of groups to match your data
Pantestudines_list

Rhynchocephalia_list <- row.names(dframe1)[dframe1[,20] == 'Rhynchocephalia']  # change names of groups to match your data
Rhynchocephalia_list

Forelimb_oscillation_list <- row.names(dframe1)[dframe1[,19] == 'Forelimb oscillation']  # change names of groups to match your data
Forelimb_oscillation_list

Rowing_list <- row.names(dframe1)[dframe1[,19] == 'Rowing']  # change names of groups to match your data
Rowing_list

Paddling_list <- row.names(dframe1)[dframe1[,19] == 'Paddling']  # change names of groups to match your data
Paddling_list

Undulation_list <- row.names(dframe1)[dframe1[,19] == 'Undulation']  # change names of groups to match your data
Undulation_list

#Hindlimb_oscillation_list <- row.names(dframe1)[dframe1[,18] == 'Hindlimb_oscillation']  # change names of groups to match your data
#Hindlimb_oscillation_list



# PCA 1-2
# MAKE EMPTY PLOT PC1 vs PC2

plot.new()
plot(pc.scores[,1], pc.scores[,2], 
     pch=19, 
     lwd=2, 
     col="transparent", 
     xlab=paste("PC1", " (", round(eigenvalues[, 2][1], 2), "%)", sep=""), 
     ylab=paste("PC2", " (", round(eigenvalues[, 2][2], 2), "%)", sep=""))
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

Plot_ConvexHull(xcoord=pc.scores[Sauropterygia_list,][,1], ycoord=pc.scores[Sauropterygia_list,][,2], lcolor="#19D27E", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Thalattosauria_list,][,1], ycoord=pc.scores[Thalattosauria_list,][,2], lcolor="#FFD300", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Saurosphargidae_list,][,1], ycoord=pc.scores[Saurosphargidae_list,][,2], lcolor="#003366", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Thalattosuchia_list,][,1], ycoord=pc.scores[Thalattosuchia_list,][,2], lcolor="blue", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Tanystropheidae_list,][,1], ycoord=pc.scores[Tanystropheidae_list,][,2], lcolor="#3AE1FF", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Mosasauroidea_list,][,1], ycoord=pc.scores[Mosasauroidea_list,][,2], lcolor="#FF8300", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Ichthyosauromorpha_list,][,1], ycoord=pc.scores[Ichthyosauromorpha_list,][,2], lcolor="#FF2626", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Pantestudines_list,][,1], ycoord=pc.scores[Pantestudines_list,][,2], lcolor="#009ccc", bgcolor=F) 
Plot_ConvexHull(xcoord=pc.scores[Rhynchocephalia_list,][,1], ycoord=pc.scores[Rhynchocephalia_list,][,2], lcolor="Purple", bgcolor=F)


Plot_ConvexHull(xcoord=pc.scores[Undulation_list,][,1], ycoord=pc.scores[Undulation_list,][,2], lcolor="orange", bgcolor="orange")
#Plot_ConvexHull(xcoord=pc.scores[Hindlimb_oscillation_list,][,1], ycoord=pc.scores[Hindlimb_oscillation_list,][,2], lcolor="orange", bgcolor="orange")
Plot_ConvexHull(xcoord=pc.scores[Paddling_list,][,1], ycoord=pc.scores[Paddling_list,][,2], lcolor="lightblue", bgcolor="lightblue")
Plot_ConvexHull(xcoord=pc.scores[Rowing_list,][,1], ycoord=pc.scores[Rowing_list,][,2], lcolor="#96b1d2", bgcolor="#96b1d2")
Plot_ConvexHull(xcoord=pc.scores[Forelimb_oscillation_list,][,1], ycoord=pc.scores[Forelimb_oscillation_list,][,2], lcolor="#3340BD", bgcolor="#3340BD")



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

#points(pc.scores[Terrestrial_list,][,1], pc.scores[Terrestrial_list,][,2], pch=21, col="grey", bg="white", cex=1)
points(pc.scores[Undulation_list,][,1], pc.scores[Undulation_list,][,2], pch=21, col="orange", bg="orange", cex=1)
#points(pc.scores[Hindlimb_oscillation_list,][,1], pc.scores[Hindlimb_oscillation_list,][,2], pch=22, col="orange", bg="orange", cex=1)
points(pc.scores[Paddling_list,][,1],pc.scores[Paddling_list,][,2], pch=24, col="lightblue", bg="lightblue", cex=1)
points(pc.scores[Rowing_list,][,1], pc.scores[Rowing_list,][,2], pch=23, col="#96b1d2", bg="#96b1d2", cex=1)
points(pc.scores[Forelimb_oscillation_list,][,1],pc.scores[Forelimb_oscillation_list,][,2], pch=25, col="#3340BD", bg="#3340BD", cex=1)



#for presentation save in 6.2 x 6.2


# MAKE SIMPLE PLOT WITH LABELS
text(pc.scores[,1], pc.scores[,2], labels=row.names(pc.scores),pos=4,cex=0.5,col="#565656")



# PCA 2-3
# MAKE EMPTY PLOT PC1 vs PC2

plot.new()
plot(pc.scores[,3], pc.scores[,2], 
     pch=19, 
     lwd=2, 
     col="transparent", 
     xlab=paste("PC3", " (", round(eigenvalues[, 2][3], 2), "%)", sep=""), 
     ylab=paste("PC2", " (", round(eigenvalues[, 2][2], 2), "%)", sep=""))
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

Plot_ConvexHull(xcoord=pc.scores[Sauropterygia_list,][,3], ycoord=pc.scores[Sauropterygia_list,][,2], lcolor="#19D27E", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Thalattosauria_list,][,3], ycoord=pc.scores[Thalattosauria_list,][,2], lcolor="#FFD300", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Saurosphargidae_list,][,3], ycoord=pc.scores[Saurosphargidae_list,][,2], lcolor="#003366", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Thalattosuchia_list,][,3], ycoord=pc.scores[Thalattosuchia_list,][,2], lcolor="blue", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Tanystropheidae_list,][,3], ycoord=pc.scores[Tanystropheidae_list,][,2], lcolor="#3AE1FF", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Mosasauroidea_list,][,3], ycoord=pc.scores[Mosasauroidea_list,][,2], lcolor="#FF8300", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Ichthyosauromorpha_list,][,3], ycoord=pc.scores[Ichthyosauromorpha_list,][,2], lcolor="#FF2626", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Pantestudines_list,][,3], ycoord=pc.scores[Pantestudines_list,][,2], lcolor="#009ccc", bgcolor=F) 
Plot_ConvexHull(xcoord=pc.scores[Rhynchocephalia_list,][,3], ycoord=pc.scores[Rhynchocephalia_list,][,2], lcolor="Purple", bgcolor=F)


Plot_ConvexHull(xcoord=pc.scores[Undulation_list,][,3], ycoord=pc.scores[Undulation_list,][,2], lcolor="orange", bgcolor="orange")
#Plot_ConvexHull(xcoord=pc.scores[Hindlimb_oscillation_list,][,1], ycoord=pc.scores[Hindlimb_oscillation_list,][,2], lcolor="orange", bgcolor="orange")
Plot_ConvexHull(xcoord=pc.scores[Paddling_list,][,3], ycoord=pc.scores[Paddling_list,][,2], lcolor="lightblue", bgcolor="lightblue")
Plot_ConvexHull(xcoord=pc.scores[Rowing_list,][,3], ycoord=pc.scores[Rowing_list,][,2], lcolor="#96b1d2", bgcolor="#96b1d2")
Plot_ConvexHull(xcoord=pc.scores[Forelimb_oscillation_list,][,3], ycoord=pc.scores[Forelimb_oscillation_list,][,2], lcolor="#3340BD", bgcolor="#3340BD")



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

#points(pc.scores[Terrestrial_list,][,3], pc.scores[Terrestrial_list,][,2], pch=21, col="grey", bg="white", cex=1)
points(pc.scores[Undulation_list,][,3], pc.scores[Undulation_list,][,2], pch=21, col="orange", bg="orange", cex=1)
#points(pc.scores[Hindlimb_oscillation_list,][,3], pc.scores[Hindlimb_oscillation_list,][,2], pch=22, col="orange", bg="orange", cex=1)
points(pc.scores[Paddling_list,][,3],pc.scores[Paddling_list,][,2], pch=24, col="lightblue", bg="lightblue", cex=1)
points(pc.scores[Rowing_list,][,3], pc.scores[Rowing_list,][,2], pch=23, col="#96b1d2", bg="#96b1d2", cex=1)
points(pc.scores[Forelimb_oscillation_list,][,3],pc.scores[Forelimb_oscillation_list,][,2], pch=25, col="#3340BD", bg="#3340BD", cex=1)



#for presentation save in 6.2 x 6.2


# MAKE SIMPLE PLOT WITH LABELS
text(pc.scores[,1], pc.scores[,2], labels=row.names(pc.scores),pos=4,cex=0.5,col="#565656")




########################
#######################


#2. PCA WITH ALL RATIOS (no cetaceans)
#------------------------------------#


#Open file again (don't eliminate the rows of the cetaceans!)
dframe1<-read.csv("FossilandExtantAqTet_data.csv", 
                  row.names = 1,
                  header=T)
dframe1


# PERFORM PCA ON CHARACTER DATA - and selects the relevant columns from the data (first column = taxon doesn't count)
res.pca <- PCA(dframe1[,1:9], scale.unit = TRUE)



# THIS CODE MAKES A PLOT SHOWING THE CONTRIBUTION OF EACH PC AXIS TO OVERALL VARIATION
eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        xlab = "Principal Components",
        ylab = "Percentage of variance",
        col ="gray")

eigenvalues

write.csv((eigenvalues), "Eigenvalues.csv") # save it as a csv file



#Variables contributions
fviz_pca_var(res.pca, col.var="contrib") # By default, variables/individuals are represented on dimensions 1 and 2
#fviz_pca_var(res.pca, col.var="contrib", axes = c(3, 2)) # specify different axis

# Simple PCA plot with labels
plot(res.pca, choix = "ind", cex=0.5, axes=c(1,2))


#############################

# RENAME THE PC SCORES (MORPHOSPACE COORDINATES)
pc.scores <- res.pca$ind$coord
pc.scores

# pc.scores <- read.table("pc.scores.final.txt", row.names=1, header=T) #
pc.scores <- as.matrix(pc.scores) 
pc.scores # shows the morphospace location coordinates of each taxa 


###############################

#FOR THE DISPARITY ANALYSIS:
#Save these scores in a separate file, that will be used later to create the file: pcscores.csv,
#write.csv(pc.scores, "pcscores.csv")
#then duplicate the taxa that spans various time bins: Taxon / Taxon_d


###############################

#FOR THE PERMANOVA:
# to make a csv file of the pcscores (the PC coordinates of each taxa) to which you can add the clade and time-bin group in new columns
#write.csv(pc.scores, "FossilAqreptiles_Output.csv")
# then add the group and time columns (I did this manually)


##############################
##############################

# MAKE PRETTY PLOT

# create lists of extant taxa grouped per LOCOMOTION and extinct taxa grouped per GROUP (Column 20)
Sauropterygia_list <- row.names(dframe1)[dframe1[,20] == 'Sauropterygia']  # change names of groups to match your data
Sauropterygia_list

Thalattosauria_list <- row.names(dframe1)[dframe1[,20] == 'Thalattosauria']  # change names of groups to match your data
Thalattosauria_list

Saurosphargidae_list <- row.names(dframe1)[dframe1[,20] == 'Saurosphargidae']  # change names of groups to match your data
Saurosphargidae_list

Tanystropheidae_list <- row.names(dframe1)[dframe1[,20] == 'Tanystropheidae']  # change names of groups to match your data
Tanystropheidae_list

Thalattosuchia_list <- row.names(dframe1)[dframe1[,20] == 'Thalattosuchia']  # change names of groups to match your data
Thalattosuchia_list

Mosasauroidea_list <- row.names(dframe1)[dframe1[,20] == 'Mosasauroidea']  # change names of groups to match your data
Mosasauroidea_list

Ichthyosauromorpha_list <- row.names(dframe1)[dframe1[,20] == 'Ichthyosauromorpha']  # change names of groups to match your data
Ichthyosauromorpha_list

Pantestudines_list <- row.names(dframe1)[dframe1[,20] == 'Pantestudines']  # change names of groups to match your data
Pantestudines_list

Rhynchocephalia_list <- row.names(dframe1)[dframe1[,20] == 'Rhynchocephalia']  # change names of groups to match your data
Rhynchocephalia_list

Forelimb_oscillation_list <- row.names(dframe1)[dframe1[,19] == 'Forelimb oscillation']  # change names of groups to match your data
Forelimb_oscillation_list

Rowing_list <- row.names(dframe1)[dframe1[,19] == 'Rowing']  # change names of groups to match your data
Rowing_list

Paddling_list <- row.names(dframe1)[dframe1[,19] == 'Paddling']  # change names of groups to match your data
Paddling_list

Undulation_list <- row.names(dframe1)[dframe1[,19] == 'Undulation']  # change names of groups to match your data
Undulation_list

Caudal_oscillation_list <- row.names(dframe1)[dframe1[,19] == 'Caudal oscillation']  # change names of groups to match your data
Caudal_oscillation_list


#Hindlimb_oscillation_list <- row.names(dframe1)[dframe1[,18] == 'Hindlimb_oscillation']  # change names of groups to match your data
#Hindlimb_oscillation_list



# PCA 1-2
# MAKE EMPTY PLOT PC1 vs PC2

plot.new()
plot(pc.scores[,1], pc.scores[,2], 
     pch=19, 
     lwd=2, 
     col="transparent", 
     xlab=paste("PC1", " (", round(eigenvalues[, 2][1], 2), "%)", sep=""), 
     ylab=paste("PC2", " (", round(eigenvalues[, 2][2], 2), "%)", sep=""))
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

Plot_ConvexHull(xcoord=pc.scores[Sauropterygia_list,][,1], ycoord=pc.scores[Sauropterygia_list,][,2], lcolor="#19D27E", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Thalattosauria_list,][,1], ycoord=pc.scores[Thalattosauria_list,][,2], lcolor="#FFD300", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Saurosphargidae_list,][,1], ycoord=pc.scores[Saurosphargidae_list,][,2], lcolor="#003366", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Thalattosuchia_list,][,1], ycoord=pc.scores[Thalattosuchia_list,][,2], lcolor="blue", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Tanystropheidae_list,][,1], ycoord=pc.scores[Tanystropheidae_list,][,2], lcolor="#3AE1FF", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Mosasauroidea_list,][,1], ycoord=pc.scores[Mosasauroidea_list,][,2], lcolor="#FF8300", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Ichthyosauromorpha_list,][,1], ycoord=pc.scores[Ichthyosauromorpha_list,][,2], lcolor="#FF2626", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Pantestudines_list,][,1], ycoord=pc.scores[Pantestudines_list,][,2], lcolor="#009ccc", bgcolor=F) 
Plot_ConvexHull(xcoord=pc.scores[Rhynchocephalia_list,][,1], ycoord=pc.scores[Rhynchocephalia_list,][,2], lcolor="Purple", bgcolor=F)


Plot_ConvexHull(xcoord=pc.scores[Undulation_list,][,1], ycoord=pc.scores[Undulation_list,][,2], lcolor="orange", bgcolor="orange")
#Plot_ConvexHull(xcoord=pc.scores[Hindlimb_oscillation_list,][,1], ycoord=pc.scores[Hindlimb_oscillation_list,][,2], lcolor="orange", bgcolor="orange")
Plot_ConvexHull(xcoord=pc.scores[Paddling_list,][,1], ycoord=pc.scores[Paddling_list,][,2], lcolor="lightblue", bgcolor="lightblue")
Plot_ConvexHull(xcoord=pc.scores[Rowing_list,][,1], ycoord=pc.scores[Rowing_list,][,2], lcolor="#96b1d2", bgcolor="#96b1d2")
Plot_ConvexHull(xcoord=pc.scores[Forelimb_oscillation_list,][,1], ycoord=pc.scores[Forelimb_oscillation_list,][,2], lcolor="#3340BD", bgcolor="#3340BD")
Plot_ConvexHull(xcoord=pc.scores[Caudal_oscillation_list,][,1], ycoord=pc.scores[Caudal_oscillation_list,][,2], lcolor="darkred", bgcolor="darkred")



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

points(pc.scores[Undulation_list,][,1], pc.scores[Undulation_list,][,2], pch=21, col="orange", bg="orange", cex=1)
#points(pc.scores[Hindlimb_oscillation_list,][,1], pc.scores[Hindlimb_oscillation_list,][,2], pch=22, col="orange", bg="orange", cex=1)
points(pc.scores[Paddling_list,][,1],pc.scores[Paddling_list,][,2], pch=24, col="lightblue", bg="lightblue", cex=1)
points(pc.scores[Rowing_list,][,1], pc.scores[Rowing_list,][,2], pch=23, col="#96b1d2", bg="#96b1d2", cex=1)
points(pc.scores[Forelimb_oscillation_list,][,1],pc.scores[Forelimb_oscillation_list,][,2], pch=25, col="#3340BD", bg="#3340BD", cex=1)
points(pc.scores[Caudal_oscillation_list,][,1],pc.scores[Caudal_oscillation_list,][,2], pch=25, col="darkred", bg="darkred", cex=1)



#for presentation save in 6.2 x 6.2


# MAKE SIMPLE PLOT WITH LABELS
text(pc.scores[,1], pc.scores[,2], labels=row.names(pc.scores),pos=4,cex=0.5,col="#565656")



# DISPARITY BETWEEN CLADES 

pc.scores

clade.groups <- list(Sauropterygia_list,Thalattosauria_list,Thalattosuchia_list,Mosasauroidea_list,Ichthyosauromorpha_list,Pantestudines_list,Rhynchocephalia_list
)
names(clade.groups) <- c("Sauropterygia", "Thalattosauria","Thalattosuchia","Mosasauroidea", "Ichthyosauromorpha", "Pantestudines", "Rhynchocephalia")
clade.groups

clade.disparity.data <- dispRity.per.group(pc.scores, group = clade.groups, metric = c(sum, variances))

clade.disparity.data #this is a dispRity object
summary(clade.disparity.data)

#Nice plots:
#plot.new()
plot(clade.disparity.data, ylab="Disparity (sum of variances)", 
     col=c("#19D27E","#FFD300","#3944bc","#FF8300","#FF2626","#009ccc","#c061ff"), 
     box(lwd=0.6),
     par(mar=c(8,5,1,1)), #sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. The default is 4
     cex.axis=1, # size of font in axis labels
     cex.lab=1.2, # size of font in axis title
     ylim=c(0,15), # limits of axis
     xaxt = "n",         #this eliminates the x labels (because there is no way to rotate them using plot function)
     xlab='')            #this eliminates the x labels

## Draw the x-axis labels
text(x = 1:7, #7 labels
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c("Sauropterygia", "Thalattosauria", "Thalattosuchia", "Mosasauroidea", "Ichthyosauromorpha", "Pantestudines", "Rhynchocephalia"), 
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.99,
     ## Increase label size.
     cex = 1)

#for presentation save 7x8

#for thesis figure save 4.5 x 5.5 inches (=height of 14 cm). Then reduce to 7 cm

#######################

# Test disparity using t-test
#https://rdrr.io/cran/dispRity/man/test.dispRity.html: TESTING disPrity hypotheses


ttest <-test.dispRity(clade.disparity.data, 
              test = t.test,
              comparisons = "pairwise",
              correction = "bonferroni") #corrected

ttest
write.csv(ttest, "ttest-results-clade-.csv")



###############################
###############################

# DISPARITY BETWEEN 8 TIME BINS
#It uses the dataset with duplicates


# This is just useful to see how many bins/groups there are: Assess time bins (Time 2)
groups_bins <- with( dframe1,
                     list(
                       groups = unique(Group_1),
                       bins   = unique(Time_2)
                     ))

sapply(groups_bins, length)


##some taxa span various time-bins, use diferent dataframe for this analysis
##Read the pc.scores dataset with duplicated taxa with different names (Taxon / Taxon_d)
#MAKE SURE THIS FILE IS UPDATED:
#pcscores-dup.csv is a file created by duplicating some taxa in the original pc.scores
pc.scores <- read.csv("pcscores-dup.csv", row.names=1, header=T) 

pc.scores <- as.matrix(pc.scores) 
pc.scores


#Create list of taxa per time bin (column 19)
#This list needs to come from the file with taxa duplicated with different names (Taxon / Taxon_d)
#so that it matches the same names as in the new pc.scores
dframe2<-read.csv("FossilAqReptiles_data_dup-2.csv", 
                  row.names = 1,
                  header=T)
dframe2

#Exclude incomplete rows
dframe2<-na.exclude(dframe2)
dframe2


EarlyTriassic_list <- row.names(dframe2)[dframe2[,19] == 'Early_Triassic']  # change names of groups to match your data
EarlyTriassic_list
MidTriassic_list <- row.names(dframe2)[dframe2[,19] == 'Mid_Triassic']  # change names of groups to match your data
MidTriassic_list
LateTriassic_list <- row.names(dframe2)[dframe2[,19] == 'Late_Triassic']  # change names of groups to match your data
LateTriassic_list
EarlyJurassic_list <- row.names(dframe2)[dframe2[,19] == 'Early_Jurassic']  # change names of groups to match your data
EarlyJurassic_list
MidJurassic_list <- row.names(dframe2)[dframe2[,19] == 'Mid_Jurassic']  # change names of groups to match your data
MidJurassic_list
LateJurassic_list <- row.names(dframe2)[dframe2[,19] == 'Late_Jurassic']  # change names of groups to match your data
LateJurassic_list
EarlyCretaceous_list <- row.names(dframe2)[dframe2[,19] == 'Early_Cretaceous']  # change names of groups to match your data
EarlyCretaceous_list
LateCretaceous_list <- row.names(dframe2)[dframe2[,19] == 'Late_Cretaceous']  # change names of groups to match your data
LateCretaceous_list


#Create data grouped by time bins
time.groups <- list(EarlyTriassic_list,MidTriassic_list,LateTriassic_list,EarlyJurassic_list,MidJurassic_list,LateJurassic_list,EarlyCretaceous_list,LateCretaceous_list)
names(time.groups) <- c("Early Triassic", "Middle Triassic","Late Triassic","Early Jurassic", "Middle Jurassic", "Late Jurassic", "Early Cretaceous", "Late Cretaceous")

time.groups

time.disparity.data <- dispRity.per.group(pc.scores, group = time.groups, metric = c(sum, variances))

time.disparity.data 
summary(time.disparity.data)

#Nice plot:
plot(time.disparity.data, ylab="Disparity (sum of variances)", 
     col=c("#7C4EF9", "#A788FC", "#D6CDF9","#42A5E8","#82E3F9","#BEF3FC","#7FD33D","#C5F99D"), #col=c("#9B6BFF", "#B99AFF", "#D8C7FF","#008080","#8EE5EE","#C6E8E9","#76BA1B","#ACDF87"), 
     box(lwd=0.6),
     par(mar=c(7,5,1,1)), #sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. The default is 4
     cex.axis=1, # size of font in axis labels
     cex.lab=1.2, # size of font in axis title
     ylim=c(0,27), # limits of axis
     xaxt = "n",         #this eliminates the x labels (because there is no way to rotate them using plot function)
     xlab='')            #this eliminates the x labels

## Draw the x-axis labels
text(x = 1:8,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c("Early Triassic", "Middle Triassic","Late Triassic","Early Jurassic", "Middle Jurassic", "Late Jurassic", "Early Cretaceous", "Late Cretaceous"), 
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.99,
     ## Increase label size.
     cex = 1)

#for presentation save 8x8

#for thesis figure save 5 x 5.5 inches (=height of 14 cm). Then reduce to 7 cm



##########################

# Test disparity using t-test
#https://rdrr.io/cran/dispRity/man/test.dispRity.html: TESTING disPrity hypotheses


ttest <-test.dispRity(time.disparity.data , 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "bonferroni") #corrected

ttest
write.csv(ttest, "ttest-results-time-corr.csv")



ttest <-test.dispRity(time.disparity.data , 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "none") #uncorrected

ttest
write.csv(ttest, "ttest-results-time-uncorr.csv")



###############################
###############################


# DISPARITY BETWEEN FIRST 10 MYA AND THE REST (INDUAN-ANISIAN - LADINIAN-CRETACEOUS)


# Assess time bins (Time 3)
groups_bins <- with( dframe1,
                     list(
                       groups = unique(Group_1),
                       bins   = unique(Time_3)
                     ))

sapply(groups_bins, length)


#Create list of taxa per broad time bin (column 20)
InduanAnisian_list <- row.names(dframe2)[dframe2[,20] == 'Induan-Anisian']  # change names of groups to match your data
InduanAnisian_list
LadinianEndCretaceous_list <- row.names(dframe2)[dframe2[,20] == 'Ladinian-EndCretaceous']  # change names of groups to match your data
LadinianEndCretaceous_list

#Create data grouped by time bins
time10Mya.groups <- list(InduanAnisian_list,LadinianEndCretaceous_list)
names(time10Mya.groups) <- c("Induan-Anisian", "Ladinian-KT")
time10Mya.groups

time10Mya.disparity.data <- dispRity.per.group(pc.scores, group = time10Mya.groups, metric = c(sum, variances))

time10Mya.disparity.data 
summary(time10Mya.disparity.data)

#Nicer plots:
plot(time10Mya.disparity.data, ylab="Disparity (sum of variances)", 
     col=c("grey","grey"), 
     box(lwd=1),
     par(mar=c(7,5,1,1)), #sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. The default is 4
     cex.axis=1, # size of font in axis labels
     cex.lab=1.2, # size of font in axis title
     ylim=c(0,27), # limits of axis
     xaxt = "n",         #this eliminates the x labels (because there is no way to rotate them using plot function)
     xlab='')            #this eliminates the x labels

## Draw the x-axis labels
text(x = 1:2,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c("Induan-Anisian", "Ladinian-KT"),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.99,
     ## Increase label size.
     cex = 1)


#for presentation save 4x8

#for thesis figure save 3 x 5.5 inches (=height of 14 cm). Then reduce to 7 cm


###################
###################


# DISPARITY BETWEEN BROAD TIME BINS (TRIASSIC - JURASSIC- CRETACEOUS) 


#Create list of taxa per broad time bin (column 17)
Triassic_list <- row.names(dframe2)[dframe2[,17] == 'Triassic']  # change names of groups to match your data
Triassic_list
Jurassic_list <- row.names(dframe2)[dframe2[,17] == 'Jurassic']  # change names of groups to match your data
Jurassic_list
Cretaceous_list <- row.names(dframe2)[dframe2[,17] == 'Cretaceous']  # change names of groups to match your data
Cretaceous_list


#Create data grouped by time bins
timebroad.groups <- list(Triassic_list,Jurassic_list,Cretaceous_list)
names(timebroad.groups) <- c("Triassic", "Jurassic", "Cretaceous")
timebroad.groups

timebroad.disparity.data <- dispRity.per.group(pc.scores, group = timebroad.groups, metric = c(sum, variances)) # CHECK: this is a wrapper that includes all steps (bootstrapping and rarefaction) ??????

timebroad.disparity.data 
summary(timebroad.disparity.data)

#Nicer plots:
#plot.new()
plot(timebroad.disparity.data, ylab="Disparity (sum of variances)", 
     col=c("#7C4EF9","#42A5E8","#7FD33D"), #col=c("#9B6BFF","#8EE5EE","#ACDF87"), 
     box(lwd=0.6),
     par(mar=c(7,5,1,1)), #sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. The default is 4
     cex.axis=1, # size of font in axis labels
     cex.lab=1.2, # size of font in axis title
     ylim=c(0,27), # limits of axis
     xaxt = "n",         #this eliminates the x labels (because there is no way to rotate them using plot function)
     xlab='')            #this eliminates the x labels

## Draw the x-axis labels
text(x = 1:3,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c("Triassic", "Jurassic", "Cretaceous"), 
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.99,
     ## Increase label size.
     cex = 1)


#for presentation save 6x8

#for thesis figure save 4 x 5.5 inches (=height of 14 cm). Then reduce to 7 cm


# Test disparity using t-test
#############################

#https://rdrr.io/cran/dispRity/man/test.dispRity.html: TESTING disPrity hypotheses


ttest <-test.dispRity(timebroad.disparity.data , 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "bonferroni") #corrected

ttest
write.csv(ttest, "ttest-results-timebroad-corr.csv")



ttest <-test.dispRity(timebroad.disparity.data , 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "none") #uncorrected

ttest
write.csv(ttest, "ttest-results-timebroad-uncorr.csv")



########################
########################


#PMANOVA: PAIRWISE NON PARAMETRIC MANOVA BETWEEN CLADES

# Make a csv file of the pcscores (the PC coordinates of each taxa) to which you can add the clade and time-bin group in new columns
# this has to be done at the begining (right after calculating the pc.scores without any duplicates)

# reading the pscore file (no duplicates for PMANOVA between clades)
pcscore.data <- read.csv("FossilAqreptiles_Output.csv", row.names=1, header=T)
scores <- pcscore.data[,1:5] # just the pc scores


pairwise.adonis <- function(x, factors, sim.method = 'euclidean', p.adjust.m ='fdr')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method, permutations = 999);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}


# permanova comparing morphospace BY CLADES
cladescore.groups <- as.vector(pcscore.data$Group_1)      # groups as vector
names(cladescore.groups) <- row.names(scores)
PERMANOVA.results <- pairwise.adonis(scores, factors=cladescore.groups, sim.method = 'euclidean', p.adjust.m ='fdr')
PERMANOVA.results                                         # view results table
write.csv(PERMANOVA.results, "Permanova by clades.csv")   # save results as a csv

# if p.adjusted is less than 0.05, then the groups have significantly different centroids (positions in morphospace)
# if p.adjusted values are close to significance threshold (0.05) - increase permutations = 999 in the function to permutations = 9999 (it will take longer to run)



# WATCH! For the next 3 PMANOVAs use the file with duplicated taxa that expands various time-bins

#The file: FossilAqreptiles_dup_Output has been created by duplicating PCA scores for taxa that is in various time bins
#the duplicates (Taxon_d) have been created manually, with the same scores than the respective Taxon


# permanova comparing morphospace BY 8 TIME BINS

#load data and, before doing the analysis and FILTER THE DATA THAT HAVE <NA> in column "Time_2"
pcscore.data <- read.csv("FossilAqreptiles_dup_Output.csv", row.names=1, header=T)
pcscore.data  <- pcscore.data[!is.na(pcscore.data$Time_2), ] #Before doing the analysis, FILTER THE DATA THAT HAVE <NA> in column "Time_2"
pcscore.data
scores <- pcscore.data[,1:5] # just the pc scores


timescore.groups <- as.vector(pcscore.data$Time_2)                  # groups as vector
timescore.groups
names(timescore.groups) <- row.names(scores)
PERMANOVA.results <- pairwise.adonis(scores, factors=timescore.groups, sim.method = 'euclidean', p.adjust.m ='fdr')
PERMANOVA.results                                                   # view results table
write.csv(PERMANOVA.results, "Permanova by time bins.csv")          # save results as a csv

# if p.adjusted is less than 0.05, then the groups have significantly different centroids (positions in morphospace)
# if p.adjusted values are close to significance threshold (0.05) - increase permutations = 999 in the function to permutations = 9999 (it will take longer to run)



# permanova comparing morphospace BY BROAD TIME BINS

#load data again #Before doing the analysis and FILTER THE DATA THAT HAVE <NA> in column "Time"
pcscore.data <- read.csv("FossilAqreptiles_dup_Output.csv", row.names=1, header=T)
pcscore.data  <- pcscore.data[!is.na(pcscore.data$Time), ]
pcscore.data
scores <- pcscore.data[,1:5] # just the pc scores


broadtimescore.groups <- as.vector(pcscore.data$Time)               # groups as vector
names(broadtimescore.groups) <- row.names(scores)
PERMANOVA.results <- pairwise.adonis(scores, factors=broadtimescore.groups, sim.method = 'euclidean', p.adjust.m ='fdr')
PERMANOVA.results                                                   # view results table
write.csv(PERMANOVA.results, "Permanova by time bins-broad.csv")    # save results as a csv

# if p.adjusted is less than 0.05, then the groups have significantly different centroids (positions in morphospace)
# if p.adjusted values are close to significance threshold (0.05) - increase permutations = 999 in the function to permutations = 9999 (it will take longer to run)



# permanova comparing morphospace First 10Mya-the rest

#load data again #Before doing the analysis and FILTER THE DATA THAT HAVE <NA> in column "Time_3"
pcscore.data <- read.csv("FossilAqreptiles_dup_Output.csv", row.names=1, header=T)
pcscore.data  <- pcscore.data[!is.na(pcscore.data$Time_3), ]
pcscore.data
scores <- pcscore.data[,1:5] # just the pc scores

time10Myascore.groups <- as.vector(pcscore.data$Time_3)            # groups as vector
names(time10Myascore.groups) <- row.names(scores)
PERMANOVA.results <- pairwise.adonis(scores, factors=time10Myascore.groups, sim.method = 'euclidean', p.adjust.m ='fdr')
PERMANOVA.results                                                  # view results table
write.csv(PERMANOVA.results, "Permanova by time 10Mya-Rest.csv")   # save results as a csv

# if p.adjusted is less than 0.05, then the groups have significantly different centroids (positions in morphospace)
# if p.adjusted values are close to significance threshold (0.05) - increase permutations = 999 in the function to permutations = 9999 (it will take longer to run)



#######################
#######################

# HEXBIN PLOT (2D DENSITY DISTRIBUTION)

library(hexbin)


# FOR THIS ANALYSIS read AGAIN the pscore file from the data with no duplicates!
pcscore.data <- read.csv("FossilAqreptiles_Output.csv", row.names=1, header=T)
scores <- pcscore.data[,1:5] # just the pc scores


###HEXBIN PCA 1-2 (R-based code)
# a file was read into pcscore.data earlier for the permanova so now I'm just reducing that file into an oject that only contains what's necessary for this
hexbin.scores <- pcscore.data[,1:2] # just the scores for PC 1 and 2

# now I'm creating a hexbin object
hexbin.scores.object <- hexbin(hexbin.scores, xbins = 10, shape = 1,
                               #xbnds = range(-4, 7), ybnds = range(-3.2, 4),
                               xlab = "PC1", ylab = "PC2", IDs = FALSE)

#EXBIN plot with BTC scale
gplot.hexbin(hexbin.scores.object, style = "colorscale", legend = 1.2, lcex = 1,
             minarea = 0.04, maxarea = 0.8, mincnt = 0, maxcnt = 9,
             trans = NULL, inv = NULL, colorcut = seq(0, 1, length = min(10, 10)),
             border = NULL, density = NULL, pen = NULL,
             colramp = function(n) BTC     (n, beg=170, end=1), #`beg' and `end' must be numbers in the interval [1,256] 
             xlab = "PC1", ylab = "PC2", main = "", newpage = TRUE,
             type = c("p", "l", "n"), xaxt = c("s", "n"), yaxt = c("s", "n"),
             clip = "on", verbose = getOption("verbose"))

##Is there a way to customise the use of colors?


###HEXBIN PCA 2-3  (R-based code)
# a file was read into pcscore.data earlier for the permanova so now I'm just reducing that file into an oject that only contains what's necessary for this
hexbin.scores <- pcscore.data[,3:2] # just the scores for PC 1 and 2

# now I'm creating a hexbin object
hexbin.scores.object <- hexbin(hexbin.scores, xbins = 10, shape = 1,
                               #xbnds = range(-4, 7), ybnds = range(-3.2, 4),
                               xlab = "PC3", ylab = "PC2", IDs = FALSE)

#EXBIN plot with BTC scale
gplot.hexbin(hexbin.scores.object, style = "colorscale", legend = 1.2, lcex = 1,
             minarea = 0.04, maxarea = 0.8, mincnt = 0, maxcnt = 9,
             trans = NULL, inv = NULL, colorcut = seq(0, 1, length = min(10, 10)),
             border = NULL, density = NULL, pen = NULL,
             colramp = function(n) BTC     (n, beg=170, end=1), #`beg' and `end' must be numbers in the interval [1,256] 
             xlab = "PC3", ylab = "PC2", main = "", newpage = TRUE,
             type = c("p", "l", "n"), xaxt = c("s", "n"), yaxt = c("s", "n"),
             clip = "on", verbose = getOption("verbose"))



# KERNEL DENSITY ESTIMATION (CONTOUR MAP) no need if using hex bin plot
library(MASS)
?kde2d

# The n= option can be changed and it will change the concentration of the contour rings
dens <- kde2d(pc.scores[,1], pc.scores[,2], n=50, lims = c(range(pc.scores[,1]), range((pc.scores[,2]))))


# to plot it
plot(pc.scores[,1], pc.scores[,2], pch=19, lwd=2, col="transparent", xlab="PC1", ylab="PC2")
abline(v=c(0), lwd=1, col="gray80", lty=2)
abline(h=c(0), lwd=1, col="gray80",lty=2)
box(lwd=1.5) # add box
k <- 25 # this can also be modified to chagne concentration of points
contour(dens, drawlabels=FALSE, nlevels=k, add=TRUE, col="gray40") # you can change the col= option - to any valid colour


# to do the same but with navy colour
plot.new()
plot(pc.scores[,1], pc.scores[,2], pch=19, lwd=2, col="transparent", xlab="PC1", ylab="PC2")
abline(v=c(0), lwd=1, col="gray80", lty=2)
abline(h=c(0), lwd=1, col="gray80",lty=2)
box(lwd=1.5) # add box
k <- 25 # this can also be modified to chagne concentration of points
contour(dens, drawlabels=FALSE, nlevels=k, add=TRUE, col="navy") # you can change the col= option - to any valid colour





