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



#TO ANALYSE ONLY THE AMNIOTES:
#Create a subset of the data only with the group categories that I want to plot
#dframe1<- subset(dframe1,Group_2 %in% c("Amniotes"))
#pca_aq_reptiles2

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
res.pca <- PCA(dframe1[,1:16], 
               scale.unit = TRUE)


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

#FOR THE GGPLOTS AND LABELS
# to make a csv file of the pcscores (the PC coordinates of each taxa) with clades, time-bins and the reduced name list
#write.csv(pc.scores, "FossilandExtantAqTet_Output.csv")
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
par(mar = c(3.1, 3.1, 1, 1),             #sets the margins of the white area around the plot (bottom, left, top, right)
    mgp=c(2.1,0.5,0))                  #sets the distance between axis title and axis labels from the axis   
plot(-pc.scores[,1], pc.scores[,2],   #if it doesn't need reversing remove the - sign from the pc.scores[,1]!!
     pch=19, 
     lwd=3, 
     col="transparent", 
     cex.axis=1.05,
     cex.lab=1.22,
     tcl = -0.3,      #length of tick marks as a fraction of the height of a line of text (default is 0.5)
     las=1,           #style of the axis labels (1 is always horizontal)
     xlab=paste("PC1", " (", round(eigenvalues[, 2][1], 2), "%)", sep=""), 
     ylab=paste("PC2", " (", round(eigenvalues[, 2][2], 2), "%)", sep=""))
     abline(v=c(0), lwd=1, col="gray80", lty=5)
     abline(h=c(0), lwd=1, col="gray80",lty=5) # plot centre of space 
     box(lwd=0.5) # add box 
par(new=TRUE) # this creates a new layer - so you can keep on plotting in same space


# DRAW HULLS

Plot_ConvexHull<-function(xcoord, ycoord, lcolor, bgcolor){
  hpts <- chull(x = -xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(-xcoord[hpts], ycoord[hpts], col = lcolor)
  polygon(x = -xcoord[hpts], y =  ycoord[hpts], col = adjustcolor(bgcolor, alpha.f = 0.20) , border = 0)
  
}  

#If it doesn't need reversing:
#Plot_ConvexHull<-function(xcoord, ycoord, lcolor, bgcolor){
#  hpts <- chull(x = xcoord, y = ycoord)
#  hpts <- c(hpts, hpts[1])
#  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
#  polygon(x = xcoord[hpts], y =  ycoord[hpts], col = adjustcolor(bgcolor, alpha.f = 0.20) , border = 0)
  
#}  


# HULLS PER CLADE

Plot_ConvexHull(xcoord=pc.scores[Sauropterygia_list,][,1], ycoord=pc.scores[Sauropterygia_list,][,2], lcolor=F, bgcolor="#19D27E")
Plot_ConvexHull(xcoord=pc.scores[Thalattosauria_list,][,1], ycoord=pc.scores[Thalattosauria_list,][,2], lcolor=F, bgcolor="#FFD300")
Plot_ConvexHull(xcoord=pc.scores[Saurosphargidae_list,][,1], ycoord=pc.scores[Saurosphargidae_list,][,2], lcolor=F, bgcolor="#003366")
Plot_ConvexHull(xcoord=pc.scores[Thalattosuchia_list,][,1], ycoord=pc.scores[Thalattosuchia_list,][,2], lcolor=F, bgcolor="#3944BC")
Plot_ConvexHull(xcoord=pc.scores[Tanystropheidae_list,][,1], ycoord=pc.scores[Tanystropheidae_list,][,2], lcolor=F, bgcolor="#3AE1FF")
Plot_ConvexHull(xcoord=pc.scores[Mosasauroidea_list,][,1], ycoord=pc.scores[Mosasauroidea_list,][,2], lcolor=F, bgcolor="#FF8300")
Plot_ConvexHull(xcoord=pc.scores[Ichthyosauromorpha_list,][,1], ycoord=pc.scores[Ichthyosauromorpha_list,][,2], lcolor=F, bgcolor="#FF2626")
Plot_ConvexHull(xcoord=pc.scores[Pantestudines_list,][,1], ycoord=pc.scores[Pantestudines_list,][,2], lcolor=F, bgcolor="#009ccc") 
Plot_ConvexHull(xcoord=pc.scores[Rhynchocephalia_list,][,1], ycoord=pc.scores[Rhynchocephalia_list,][,2], lcolor=F, bgcolor="Purple")


Plot_ConvexHull(xcoord=pc.scores[Paddling_list,][,1], ycoord=pc.scores[Paddling_list,][,2], lcolor="grey", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Undulation_list,][,1], ycoord=pc.scores[Undulation_list,][,2], lcolor="grey", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Rowing_list,][,1], ycoord=pc.scores[Rowing_list,][,2], lcolor="grey", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Forelimb_oscillation_list,][,1], ycoord=pc.scores[Forelimb_oscillation_list,][,2], lcolor="grey", bgcolor=F)



# DRAW POINTS PER CLADE #(if there is no need for reversing remove the sign- from pc.scores)

points(-pc.scores[Sauropterygia_list,][,1], pc.scores[Sauropterygia_list,][,2], pch=21, col="#19D27E", bg="#19D27E", cex=0.7)
points(-pc.scores[Thalattosauria_list,][,1], pc.scores[Thalattosauria_list,][,2], pch=25, col="#FFD300", bg="#FFD300", cex=0.7)
points(-pc.scores[Saurosphargidae_list,][,1],pc.scores[Saurosphargidae_list,][,2], pch=21, col="darkgrey", bg="darkgrey", cex=0.7)
points(-pc.scores[Thalattosuchia_list,][,1], pc.scores[Thalattosuchia_list,][,2], pch=22, col="#3944BC", bg="#3944BC", cex=0.7)
points(-pc.scores[Tanystropheidae_list,][,1], pc.scores[Tanystropheidae_list,][,2], pch=24, col="#c4f4ff", bg="#c4f4ff", cex=0.7)
points(-pc.scores[Mosasauroidea_list,][,1], pc.scores[Mosasauroidea_list,][,2], pch=23, col="#FF8300", bg="#FF8300", cex=0.7)
points(-pc.scores[Ichthyosauromorpha_list,][,1],pc.scores[Ichthyosauromorpha_list,][,2], pch=24, col="#FF2626", bg="#FF2626", cex=0.7)
points(-pc.scores[Pantestudines_list,][,1], pc.scores[Pantestudines_list,][,2], pch=23, col="#0083ab", bg="#0083ab", cex=0.7)
points(-pc.scores[Rhynchocephalia_list,][,1],pc.scores[Rhynchocephalia_list,][,2], pch=25, col="#c061ff", bg="#c061ff", cex=0.7)


points(-pc.scores[Undulation_list,][,1], pc.scores[Undulation_list,][,2], pch=21, col="black", bg="#878787", cex=1, lwd = 0.7)
points(-pc.scores[Paddling_list,][,1],pc.scores[Paddling_list,][,2], pch=24, col="black", bg="white", cex=1,lwd = 0.7)
points(-pc.scores[Rowing_list,][,1], pc.scores[Rowing_list,][,2], pch=23, col="black", bg="lightgrey", cex=1,lwd = 0.7)  
points(-pc.scores[Forelimb_oscillation_list,][,1],pc.scores[Forelimb_oscillation_list,][,2], pch=22, col="black", bg="black", cex=1,lwd = 0.7)


#for presentation save in 7x7
#for figure save 5.2 x 5.2 (aprox 13.2 cm) and then scale to 6.6 cm

# MAKE SIMPLE PLOT WITH LABELS
text(-pc.scores[,1], pc.scores[,2], labels=row.names(pc.scores),pos=4,cex=0.3,col="#565656")


########################
#######################


#2. PCA WITH REDUCED RATIOS (+ cetaceans)
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

#write.csv((eigenvalues), "Eigenvalues.csv") # save it as a csv file


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

#FOR THE GGPLOTS AND LABELS
# to make a csv file of the pcscores (the PC coordinates of each taxa) with clades, time-bins and the reduced name list
write.csv(pc.scores, "FossilandExtantAqTet_Output-reduced.csv")
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


# MAKE EMPTY PLOT
plot.new()
par(mar = c(3.1, 3.1, 1, 1),             #sets the margins of the white area around the plot (bottom, left, top, right)
    mgp=c(2.1,0.5,0))                  #sets the distance between axis title and axis labels from the axis   
plot(pc.scores[,1], pc.scores[,2],   #if it doesn't need reversing remove the - sign from the pc.scores[,1]!!
     pch=19, 
     lwd=3, 
     col="transparent", 
     cex.axis=1.05,
     cex.lab=1.22,
     tcl = -0.3,      #length of tick marks as a fraction of the height of a line of text (default is 0.5)
     las=1,           #style of the axis labels (1 is always horizontal)
     xlab=paste("PC1", " (", round(eigenvalues[, 2][1], 2), "%)", sep=""), 
     ylab=paste("PC2", " (", round(eigenvalues[, 2][2], 2), "%)", sep=""))
abline(v=c(0), lwd=1, col="gray80", lty=5)
abline(h=c(0), lwd=1, col="gray80",lty=5) # plot centre of space 
box(lwd=0.5) # add box 
par(new=TRUE) # this creates a new layer - so you can keep on plotting in same space


# DRAW HULLS

Plot_ConvexHull<-function(xcoord, ycoord, lcolor, bgcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
  polygon(x = xcoord[hpts], y =  ycoord[hpts], col = adjustcolor(bgcolor, alpha.f = 0.20) , border = 0)
  
}  

# HULLS PER CLADE
Plot_ConvexHull(xcoord=pc.scores[Sauropterygia_list,][,1], ycoord=pc.scores[Sauropterygia_list,][,2], lcolor=F, bgcolor="#19D27E")
Plot_ConvexHull(xcoord=pc.scores[Thalattosauria_list,][,1], ycoord=pc.scores[Thalattosauria_list,][,2], lcolor=F, bgcolor="#FFD300")
Plot_ConvexHull(xcoord=pc.scores[Saurosphargidae_list,][,1], ycoord=pc.scores[Saurosphargidae_list,][,2], lcolor=F, bgcolor="#003366")
Plot_ConvexHull(xcoord=pc.scores[Thalattosuchia_list,][,1], ycoord=pc.scores[Thalattosuchia_list,][,2], lcolor=F, bgcolor="#3944BC")
Plot_ConvexHull(xcoord=pc.scores[Tanystropheidae_list,][,1], ycoord=pc.scores[Tanystropheidae_list,][,2], lcolor=F, bgcolor="#3AE1FF")
Plot_ConvexHull(xcoord=pc.scores[Mosasauroidea_list,][,1], ycoord=pc.scores[Mosasauroidea_list,][,2], lcolor=F, bgcolor="#FF8300")
Plot_ConvexHull(xcoord=pc.scores[Ichthyosauromorpha_list,][,1], ycoord=pc.scores[Ichthyosauromorpha_list,][,2], lcolor=F, bgcolor="#FF2626")
Plot_ConvexHull(xcoord=pc.scores[Pantestudines_list,][,1], ycoord=pc.scores[Pantestudines_list,][,2], lcolor=F, bgcolor="#009ccc") 
Plot_ConvexHull(xcoord=pc.scores[Rhynchocephalia_list,][,1], ycoord=pc.scores[Rhynchocephalia_list,][,2], lcolor=F, bgcolor="Purple")


# DRAW POINTS PER CLADE
points(pc.scores[Sauropterygia_list,][,1], pc.scores[Sauropterygia_list,][,2], pch=21, col="#19D27E", bg="#19D27E", cex=0.7)
points(pc.scores[Thalattosauria_list,][,1], pc.scores[Thalattosauria_list,][,2], pch=25, col="#FFD300", bg="#FFD300", cex=0.7)
points(pc.scores[Saurosphargidae_list,][,1],pc.scores[Saurosphargidae_list,][,2], pch=21, col="darkgrey", bg="darkgrey", cex=0.7)
points(pc.scores[Thalattosuchia_list,][,1], pc.scores[Thalattosuchia_list,][,2], pch=22, col="#3944BC", bg="#3944BC", cex=0.7)
points(pc.scores[Tanystropheidae_list,][,1], pc.scores[Tanystropheidae_list,][,2], pch=24, col="#c4f4ff", bg="#c4f4ff", cex=0.7)
points(pc.scores[Mosasauroidea_list,][,1], pc.scores[Mosasauroidea_list,][,2], pch=23, col="#FF8300", bg="#FF8300", cex=0.7)
points(pc.scores[Ichthyosauromorpha_list,][,1],pc.scores[Ichthyosauromorpha_list,][,2], pch=24, col="#FF2626", bg="#FF2626", cex=0.7)
points(pc.scores[Pantestudines_list,][,1], pc.scores[Pantestudines_list,][,2], pch=23, col="#0083ab", bg="#0083ab", cex=0.7)
points(pc.scores[Rhynchocephalia_list,][,1],pc.scores[Rhynchocephalia_list,][,2], pch=25, col="#c061ff", bg="#c061ff", cex=0.7)


# DRAW HULLS FOR LOCOMOTORY TYPES
Plot_ConvexHull(xcoord=pc.scores[Paddling_list,][,1], ycoord=pc.scores[Paddling_list,][,2], lcolor="lightgrey", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Undulation_list,][,1], ycoord=pc.scores[Undulation_list,][,2], lcolor="grey", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Rowing_list,][,1], ycoord=pc.scores[Rowing_list,][,2], lcolor="grey", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Forelimb_oscillation_list,][,1], ycoord=pc.scores[Forelimb_oscillation_list,][,2], lcolor="grey", bgcolor=F)
Plot_ConvexHull(xcoord=pc.scores[Caudal_oscillation_list,][,1], ycoord=pc.scores[Caudal_oscillation_list,][,2], lcolor="grey", bgcolor=F)



# DRAW POINTS FOR LOCOMOTORY TYPES
points(pc.scores[Undulation_list,][,1], pc.scores[Undulation_list,][,2], pch=21, col="black", bg="#878787", cex=1, lwd = 0.7)
points(pc.scores[Paddling_list,][,1],pc.scores[Paddling_list,][,2], pch=24, col="black", bg="white", cex=1,lwd = 0.7)
points(pc.scores[Rowing_list,][,1], pc.scores[Rowing_list,][,2], pch=23, col="black", bg="lightgrey", cex=1,lwd = 0.7)  
points(pc.scores[Forelimb_oscillation_list,][,1],pc.scores[Forelimb_oscillation_list,][,2], pch=22, col="black", bg="black", cex=1,lwd = 0.7)
points(pc.scores[Caudal_oscillation_list,][,1],pc.scores[Caudal_oscillation_list,][,2], pch=21, col="black", bg="black", cex=1, lwd = 0.7)


#print legend for locomotory types:

legend(3, 2, 
       legend=c("undulators", "paddlers", "rowers","underwater fliers", "caudal oscillators"), 
       cex=0.9,
       col=c("black", "black", "black", "black", "black"),
       pch =c(1, 2, 5, 15, 15)) 
     # appropriate legend, the end numbers are the symbols, and the first numbers are the location of the box in terms of the top left corner on the x and y axes



#for presentation save in 7x7
#for figure save 5.2 x 5.2 (aprox 13.2 cm) and then scale to 6.6 cm


# MAKE SIMPLE PLOT WITH LABELS
text(pc.scores[,1], pc.scores[,2], labels=row.names(pc.scores),pos=4,cex=0.3,col="#565656")


#########################
########################


##PCA1-2_ALL RATIOS, WITH NUMBERED LABELS (FIGURE)
#------------------------------------------------#

rm(list = ls())
#setwd("######")

# OPEN FILE
dframe1<-read.csv("FossilandExtantAqTet_Output.csv", 
                  row.names = 1,
                  header=T)
dframe1

Locomotion1 <- factor(dframe1$Locomotion, 
                      levels=c("Ichthyosauromorpha", "Mosasauroidea", "Pantestudines", "Rhynchocephalia", "Sauropterygia","Saurosphargidae","Thalattosauria", "Thalattosuchia", "Tanystropheidae",
                               "Paddling", "Rowing","Undulation", "Forelimb oscillation"),
                      labels = c("Ichthyosauromorpha", "Mosasauroidea", "Pantestudines", "Rhynchocephalia", "Sauropterygia","Saurosphargidae","Thalattosauria", "Thalattosuchia", "Tanystropheidae",
                                 "paddlers", "rowers","undulators", "underwater fliers"))

Locomotion1


pca_plot <- ggplot(dframe1, aes(x = -Dim.1, y = Dim.2,fill = Locomotion1, colour=Locomotion1, shape=Locomotion1))+
  geom_point(size=2, colour="black", stroke=0.4)+
  xlim(-6.5,5.2)+
  scale_colour_manual(values= c("#FF2626","#FF8300","#0083ab","#c061ff","#19D27E","darkgrey","#FFD300","#3944BC", "#add8f0","#303030","#303030","#303030","#303030")) +
  scale_fill_manual(values= c("#FF2626","#FF8300","#0083ab","#c061ff","#19D27E","darkgrey","#FFD300", "#3944BC", "#add8f0","white","lightgrey","#878787","black"))+
  scale_shape_manual(values= c(24,23,21,25,21,24,25,22,21,24,23,21,22))+  
  #scale_shape_manual(values= c(3,3,3,3,3,3,3,3,3,24,23,21,22))+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #horizontal reference line in O
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #vertical reference line in O
  labs(x = "PC1 (25.2%)",
       y = "PC2 (21.75%)")+
  #ggrepel::geom_text_repel(label = dframe1$Name, size=3)+     #package to avoid overlap for crowded labels
  #ggrepel::geom_text_repel(aes(label = Name),
  #                         data = subset (dframe1, Time_2 < 2),
  #                         size=2.8,
  #                         segment.size=0.2,     #thickness of label thickness
  #                         segment.alpha=0.5,    #alpha for labels
  #)+
  geom_text(aes(label=Num, color=Locomotion1), size=3.5, hjust=-0.4)+

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
ggsave("PCA 1-2_All ratios-labels_2.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 16, height = 16, units = "cm",
       dpi = 300, limitsize = TRUE) 


#Zoom of the right-inferior quadrant
pca_plot <- ggplot(dframe1, aes(x = -Dim.1, y = Dim.2,fill = Locomotion1, colour=Locomotion1, shape=Locomotion1))+
  geom_point(size=2.5, stroke=0.4, color="black")+
  xlim(-0.7,2.7)+
  ylim(-2.2,0.64)+
  scale_colour_manual(values= c("#FF2626","#FF8300","#0083ab","#c061ff","#19D27E","darkgrey","#FFD300","#3944BC", "#add8f0","#303030","#303030","#303030","#303030")) +
  scale_fill_manual(values= c("#FF2626","#FF8300","#0083ab","#c061ff","#19D27E","darkgrey","#FFD300", "#3944BC", "#add8f0","white","lightgrey","#878787","black"))+
  scale_shape_manual(values= c(24,23,21,25,21,24,25,22,21,24,23,21,22))+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #horizontal reference line in O
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #vertical reference line in O
  labs(x = "PC1 (25.2%)",
       y = "PC2 (21.75%)")+
  #ggrepel::geom_text_repel(label = dframe1$Name, size=3)+     #package to avoid overlap for crowded labels
  #ggrepel::geom_text_repel(aes(label = Name),
  #                         data = subset (dframe1, Time_2 < 2),
  #                         size=2.8,
  #                         segment.size=0.2,     #thickness of label thickness
  #                         segment.alpha=0.5,    #alpha for labels
  #)+
  #geom_text(aes(label=Num, color=Locomotion1), size=5, hjust=-0.2)+
  ggrepel::geom_text_repel(aes(label=Num,color=Locomotion1), 
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
                   aspect.ratio = 1.08,
                   legend.position = "none")


#will reduce to half size, so I use double size of font
ggsave("PCA 1-2_All ratios-labels_zoom.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 13, height = 13 , units = "cm",
       dpi = 300, limitsize = TRUE) 



######################
######################


##PCA1-2_REDUCED RATIOS -added cetaceans-, WITH NUMBERED LABELS (FIGURE)
#----------------------------------------------------------------------#

rm(list = ls())
#setwd("######")

# OPEN FILE
dframe1<-read.csv("FossilandExtantAqTet_Output-reduced.csv", 
                  row.names = 1,
                  header=T)


dframe1


Locomotion1 <- factor(dframe1$Locomotion, 
                      levels=c("Ichthyosauromorpha", "Mosasauroidea", "Pantestudines", "Rhynchocephalia", "Sauropterygia","Saurosphargidae","Thalattosauria", "Thalattosuchia", "Tanystropheidae",
                               "Paddling", "Rowing","Undulation", "Forelimb oscillation", "Caudal oscillation"),
                      labels = c("Ichthyosauromorpha", "Mosasauroidea", "Pantestudines", "Rhynchocephalia", "Sauropterygia","Saurosphargidae","Thalattosauria", "Thalattosuchia", "Tanystropheidae",
                                 "paddlers", "rowers","undulators", "underwater fliers", "caudal oscillators"))

Locomotion1


pca_plot <- ggplot(dframe1, aes(x = Dim.1, y = Dim.2, fill = Locomotion1, colour=Locomotion1, shape=Locomotion1))+
  geom_point(size=2, stroke=0.4,colour="black")+
  scale_colour_manual(values= c("#FF2626","#FF8300","#0083ab","#c061ff","#19D27E","darkgrey","#FFD300","#3944BC", "#add8f0","#303030","#303030","#303030","#303030","#303030")) +
  scale_fill_manual(values= c("#FF2626","#FF8300","#0083ab","#c061ff","#19D27E","darkgrey","#FFD300", "#3944BC", "#add8f0","white","lightgrey","#878787","black","black"))+  
  scale_y_continuous (labels = scales::number_format(accuracy = 0.1))+
  scale_shape_manual(values= c(24,23,21,25,21,24,25,22,21,24,23,21,22,21))+ 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #horizontal reference line in O
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #vertical reference line in O
  labs(x = "PC1 (31.48%)",
       y = "PC2 (21.57%)")+
  geom_text(aes(label=Num, color=Locomotion1), size=3.5, hjust=-0.4
  )+
  
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
ggsave("PCA 1-2_Reduced ratios-labels_2.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 16, height = 16, units = "cm",
       dpi = 300, limitsize = TRUE) 

#ggsave("LEGEND2.pdf", plot = last_plot(), device = "pdf", path = NULL,
#       scale = 1, width = 15, height = 15, units = "cm",
#       dpi = 300, limitsize = TRUE) 




#Zoom of the left-inferior quadrant
pca_plot <- ggplot(dframe1, aes(x = Dim.1, y = Dim.2, fill = Locomotion1, colour=Locomotion1, shape=Locomotion1))+
  geom_point(size=2.5, stroke=0.4, color="black")+
  xlim(-2.7,1.3)+
  ylim(-2.2,1.8)+
  scale_colour_manual(values= c("#FF2626","#FF8300","#0083ab","#c061ff","#19D27E","darkgrey","#FFD300","#3944BC", "#add8f0","#303030","#303030","#303030","#303030","#303030")) +
  scale_fill_manual(values= c("#FF2626","#FF8300","#0083ab","#c061ff","#19D27E","darkgrey","#FFD300", "#3944BC", "#add8f0","white","lightgrey","#878787","black","black"))+  
  scale_shape_manual(values= c(24,23,21,25,21,24,25,22,21,24,23,21,22,21))+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #horizontal reference line in O
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, col="gray80")+  #vertical reference line in O
  labs(x = "PC1 (25.2%)",
       y = "PC2 (21.75%)")+
  ggrepel::geom_text_repel(aes(label=Num,color=Locomotion1), 
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
                   aspect.ratio = 4/3.5,
                   legend.position = "none")


#will reduce to half size, so I use double size of font
ggsave("PCA 1-2_Reduced ratios-labels_zoom.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, height = 13, units = "cm",
       dpi = 300, limitsize = TRUE) 


#############################

# VARIABLE CONTRIBUTION: plots showing CHARACTERS AFFECTING THE PCA AXES 
#------------------------------------------------------------------------#

#http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining



#Variable contribution from factoextra (FIGURE)
#--------------------------------------------#

#PCA1-2 All ratios
p<-fviz_pca_var(res.pca, 
                col.var="contrib", 
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
           legend.position= c(1.1, 0.82),
           legend.title=element_text(size=10))

#legend.title = element_text( c=("legend")))
#legend.position = "none")

#will reduce to half size, so I use double size of font
ggsave("PCA 1-2_All ratios-contribution.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 17.5, height = 14, units = "cm",
       dpi = 300, limitsize = TRUE) 





#PCA1-2 Reduced ratios
p<-fviz_pca_var(res.pca, 
                col.var="contrib", 
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
           legend.position= c(1.1, 0.82),
           legend.title=element_text(size=10))

#legend.title = element_text( c=("legend")))
#legend.position = "none")

#will reduce to half size, so I use double size of font
ggsave("PCA 1-2_Reduced ratios-contribution.pdf", plot = last_plot(), device = "pdf", path = NULL,
       scale = 1, width = 17.5, height = 14, units = "cm",
       dpi = 300, limitsize = TRUE) 

