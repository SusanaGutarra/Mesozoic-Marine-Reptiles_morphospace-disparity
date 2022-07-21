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
dframe1<-dframe1[!duplicated(dframe1$Taxon),]

## Use the first column for row names
dframe1 <- data.frame(dframe1, row.names = 1)
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

## there are some very strong patterns going on


#################################


# PERFORM PCA ON CHARACTER DATA - and selects the relevant columns from the data (first column = taxon doesn't count)

res.pca <- PCA(dframe1[,1:16], scale.unit = TRUE)

res.pca


# THIS CODE MAKES A PLOT SHOWING THE CONTRIBUTION OF EACH PC AXIS TO OVERALL VARIATION
eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        xlab = "Principal Components",
        ylab = "Percentage of variance",
        col ="gray")

eigenvalues

#write.csv((eigenvalues), "Eigenvalues.csv") # save it as a csv file


# PLOT SHOWING THE CHARACTERS AFFECTING THE PCA AXES 
plot(res.pca, choix = "var", axes=c(1,2))
plot(res.pca, choix = "var", axes=c(3,2))


#############################

# RENAME THE PC SCORES (MORPHOSPACE COORDINATES)
pc.scores <- res.pca$ind$coord #we can specify how many axis to show by adding [,1:5] 
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

# create lists of taxa grouped per CLADE (Column 18)
Sauropterygia_list <- row.names(dframe1)[dframe1[,18] == 'Sauropterygia']  # change names of groups to match your data
Sauropterygia_list

Thalattosauria_list <- row.names(dframe1)[dframe1[,18] == 'Thalattosauria']  # change names of groups to match your data
Thalattosauria_list

Saurosphargidae_list <- row.names(dframe1)[dframe1[,18] == 'Saurosphargidae']  # change names of groups to match your data
Saurosphargidae_list

Tanystropheidae_list <- row.names(dframe1)[dframe1[,18] == 'Tanystropheidae']  # change names of groups to match your data
Tanystropheidae_list

Thalattosuchia_list <- row.names(dframe1)[dframe1[,18] == 'Thalattosuchia']  # change names of groups to match your data
Thalattosuchia_list

Mosasauroidea_list <- row.names(dframe1)[dframe1[,18] == 'Mosasauroidea']  # change names of groups to match your data
Mosasauroidea_list

Ichthyosauromorpha_list <- row.names(dframe1)[dframe1[,18] == 'Ichthyosauromorpha']  # change names of groups to match your data
Ichthyosauromorpha_list

Pantestudines_list <- row.names(dframe1)[dframe1[,18] == 'Pantestudines']  # change names of groups to match your data
Pantestudines_list

Rhynchocephalia_list <- row.names(dframe1)[dframe1[,18] == 'Rhynchocephalia']  # change names of groups to match your data
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
#######################



# DISPARITY BETWEEN CLADES 
#------------------------#

#making groups
clade.groups <- list(Sauropterygia_list,Thalattosauria_list,Thalattosuchia_list,Mosasauroidea_list,Ichthyosauromorpha_list,Pantestudines_list,Rhynchocephalia_list
)
names(clade.groups) <- c("Sauropterygia", "Thalattosauria","Thalattosuchia","Mosasauroidea", "Ichthyosauromorpha", "Pantestudines", "Rhynchocephalia")
clade.groups

#disparity calculation using "dispRity.per.group"wrapper function, with 100 bootstraps default
cladedisparity <- dispRity.per.group(pc.scores,  
                                     group = clade.groups, 
                                     metric = c(sum, variances)) #for sum of ranges use c(sum, ranges)


#adding boot.matrix function for a finer control over bootstrap and rarefaction
#boot <- custom.subsets(pc.scores,
#                       group = clade.groups) %>%
#        boot.matrix(bootstraps = 500,            
#                       rarefaction = TRUE)       
#boot 

#disparity calculation
#cladedisparity <- dispRity(boot,                        #produces a dispRity object
#                           metric = c(sum, variances)) #for sum of ranges use c(sum, ranges)


summary(cladedisparity)

#Nice plots:
plot.new()
plot(cladedisparity, ylab="Disparity (sum of variances)", 
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



# DISPARITY BETWEEN FIRST 10 MYA AND THE REST (INDUAN-ANISIAN - LADINIAN-CRETACEOUS)
#----------------------------------------------------------------------------------#


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
     cex.axis=1.2, # size of font in axis labels
     cex.lab=1.3, # size of font in axis title
     ylim=c(0,25), # limits of axis
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
     cex = 1.3)


#for presentation save 4x8

#for thesis figure save 2.74 x 6 inches, then reduce


###################
###################


# DISPARITY BETWEEN BROAD TIME BINS (TRIASSIC - JURASSIC- CRETACEOUS) 
#-------------------------------------------------------------------#


#Create list of taxa per broad time bin (column 17)
Triassic_list <- row.names(dframe2)[dframe2[,18] == 'Triassic']  # change names of groups to match your data
Triassic_list
Jurassic_list <- row.names(dframe2)[dframe2[,18] == 'Jurassic']  # change names of groups to match your data
Jurassic_list
Cretaceous_list <- row.names(dframe2)[dframe2[,18] == 'Cretaceous']  # change names of groups to match your data
Cretaceous_list


#Create data grouped by time bins
timebroad.groups <- list(Triassic_list,Jurassic_list,Cretaceous_list)
names(timebroad.groups) <- c("Triassic", "Jurassic", "Cretaceous")
timebroad.groups

timebroad.disparity.data <- dispRity.per.group(pc.scores, group = timebroad.groups, metric = c(sum, variances))

timebroad.disparity.data 
summary(timebroad.disparity.data)

#Nicer plots:
#plot.new()
plot(timebroad.disparity.data, ylab="Disparity (sum of variances)", 
     #col=c("#7C4EF9","#42A5E8","#7FD33D"), #col=c("#9B6BFF","#8EE5EE","#ACDF87"), 
     box(lwd=0.6),
     par(mar=c(7,5,1,1)), #sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. The default is 4
     cex.axis=1.2, # size of font in axis labels
     cex.lab=1.3, # size of font in axis title
     ylim=c(0,25), # limits of axis
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
     cex = 1.3)


#for presentation save 6x8

#for thesis figure save 2.9 x 6 inches

# Test disparity using t-test
#############################

#https://rdrr.io/cran/dispRity/man/test.dispRity.html: TESTING disPrity hypotheses


ttest <-test.dispRity(timebroad.disparity.data , 
                      test = t.test,
                      #comparisons = "sequential",
                      correction = "bonferroni") #corrected



ttest
write.csv(ttest, "ttest-results-timebroad-corr.csv")



ttest <-test.dispRity(timebroad.disparity.data , 
                      test = t.test,
                      #comparisons = "sequential",
                      correction = "none") #uncorrected

ttest
write.csv(ttest, "ttest-results-timebroad-uncorr.csv")

########################
########################


# DISPARITY BETWEEN 13 TIME BINS
#---------------------------------#


##some taxa span various time-bins
##Read the pc.scores dataset with duplicated taxa with different names (Taxon / Taxon_d)
#MAKE SURE THIS FILE IS UPDATED:
#pcscores-dup.csv is a file created by duplicating some taxa in the original pc.scores
pc.scores <- read.csv("pcscores-dup.csv", 
                      row.names=1, 
                      header=T) 

pc.scores <- as.matrix(pc.scores) 
pc.scores



#Create list of taxa per time bin (column 19)
#This list needs to come from the file with taxa duplicated with different names (Taxon / Taxon_d)
#so that it matches the same names as in the new pc.scores
dframe2<-read.csv("FossilAqReptiles_data_dup.csv", 
                  row.names = 1,
                  header=T)
dframe2

#Exclude incomplete rows
dframe2<-na.exclude(dframe2)
dframe2


#Assess number of bins/groups there are in Time_4 
#(there should be 14, including the category "DUP", not to be accounted for)
groups_bins <- with( dframe2,
                     list(
                       groups = unique(Group_1),
                       bins   = unique(Time_4)
                     ))

sapply(groups_bins, length)


IndOle_list <- row.names(dframe2)[dframe2[,21] == 'Ind-Ole']  # change names of groups to match your data
IndOle_list
AniLad_list <- row.names(dframe2)[dframe2[,21] == 'Ani-Lad']  # change names of groups to match your data
AniLad_list
Car_list <- row.names(dframe2)[dframe2[,21] == 'Car']  # change names of groups to match your data
Car_list
NorRha_list <- row.names(dframe2)[dframe2[,21] == 'Nor-Rha']  # change names of groups to match your data
NorRha_list
HetSin_list <- row.names(dframe2)[dframe2[,21] == 'Het-Sin']  # change names of groups to match your data
HetSin_list
PliToa_list <- row.names(dframe2)[dframe2[,21] == 'Pli-Toa']  # change names of groups to match your data
PliToa_list
AalCall_list <- row.names(dframe2)[dframe2[,21] == 'Aal-Call']  # change names of groups to match your data
AalCall_list
Oxf_list <- row.names(dframe2)[dframe2[,21] == 'Oxf']  # change names of groups to match your data
Oxf_list
KimTith_list <- row.names(dframe2)[dframe2[,21] == 'Kim-Tith']  # change names of groups to match your data
KimTith_list
BerBar_list <- row.names(dframe2)[dframe2[,21] == 'Ber-Bar']  # change names of groups to match your data
BerBar_list
AptAlb_list <- row.names(dframe2)[dframe2[,21] == 'Apt-Alb']  # change names of groups to match your data
AptAlb_list
CenCon_list <- row.names(dframe2)[dframe2[,21] == 'Cen-Con']  # change names of groups to match your data
CenCon_list
SanMaas_list <- row.names(dframe2)[dframe2[,21] == 'San-Maas']  # change names of groups to match your data
SanMaas_list



#Create data grouped by time bins
time.groups <- list(IndOle_list,AniLad_list,Car_list,NorRha_list,HetSin_list,PliToa_list,AalCall_list,Oxf_list,
                    KimTith_list,BerBar_list,AptAlb_list,CenCon_list,SanMaas_list)
names(time.groups) <- c("Ind-Ole", "Ani-Lad","Car","Nor-Rha", "Het-Sin", "Pli-Toa", "Aal-Call", "Oxf",
                        "Kim-Tith", "Ber-Bar","Apt-Alb","Cen-Con", "San-Maas")

time.groups


#Calculate disparity:
timedisparity <- dispRity.per.group(pc.scores, 
                                    group = time.groups, 
                                    metric = c(sum, variances))

timedisparity
summary(timedisparity)



#Not using wrapper function
timedisparity  <-dispRity(boot.matrix(custom.subsets(pc.scores,  
                                                     group = time.groups),
                                      bootstraps = 500),
                          metric = c(sum, variances))

timedisparity 
summary(timedisparity)



#Nice plot:
plot.new ()
plot(timedisparity, ylab="Disparity (sum of variances)", 
     col=c("#7C4EF9", "#A788FC", "#D6CDF9","#D6CDF9","#42A5E8","#42A5E8","#82E3F9","#BEF3FC","#BEF3FC","#7FD33D","#7FD33D","#C5F99D","#C5F99D"), 
     par(mar=c(7,5,1,1)), #sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. The default is 4
     cex.axis=1, # size of font in axis labels
     cex.lab=1.2, # size of font in axis title
     ylim=c(0,30), # limits of axis
     xaxt = "n",         #this eliminates the x labels (because there is no way to rotate them using plot function)
     xlab='')            #this eliminates the x labels

## Draw the x-axis labels
text(x = 1:13,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c("Ind-Ole", "Ani-Lad","Car","Nor-Rha", "Het-Sin", "Pli-Toa", "Aal-Call", "Oxf",
                "Kim-Tith", "Ber-Bar","Apt-Alb","Cen-Con", "San-Maas"), 
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

####################
####################

#1-Nice plot WITH NUMBERS:
plot.new ()

plot(timedisparity, ylab="Disparity (sum of variances)", 
     col=c("#7C4EF9", "#A788FC", "#D6CDF9","#42A5E8","#82E3F9","#BEF3FC","#7FD33D","#C5F99D"),
     at = c(250, 242.1, 232, 214.1, 196.05, 182.5, 168.8, 160.4,151.15, 135, 112.75, 93.4, 76.15),  #Assigns numerical numbers to the categoric time bins. Ben recommended at = c(252, 247, 220, 190, 170, 155, 125, 80) 
     box(lwd=0.6),
     par(mar=c(7,5,1,1)), #sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. The default is 4
     cex.axis=1, # size of font in axis labels
     cex.lab=1.2, # size of font in axis title
     ylim=c(0,28), # limits of axis
     xaxt = "n",       #this eliminates the x labels (because there is no way to rotate them using plot function)
     xlab="Time (Ma)")            #this eliminates the x labels



## 2- using Geoscale plots
#https://cran.r-project.org/web/packages/geoscale/geoscale.pdf
#https://rdrr.io/cran/geoscale/man/geoscaleBox.html

require(geoscale)

## Converting the data into a list
data_obs <- extract.dispRity(timedisparity, observed = TRUE)
data_distribution <- extract.dispRity(timedisparity, observed = FALSE)
## Removing one list level
data_distribution <- unlist(data_distribution, recursive = FALSE)
data_obs <- as.vector(data_obs)


#assigning numerical values to the ages
ages <- c(251, 242.1, 232, 214.1, 196, 182.4, 168.8, 160.4, 151.1, 135, 113, 93.4, 76.1)
ages



## Plotting the results distribution: THIS PIECE OF CODE IS A BIT BUGGY
geoscaleBox(data_distribution, ages, 
            boxes = "Epoch", 
            data.lim = c(2,28),  #limits of the y axis
            age.lim= c(251,66),  #limits of the x axis
            units=c("Period"),   #time units shown in geological scale
            box.width=3,         #width of the boxplot in millions of years
            tick.scale = 50,     #Time scale resolution.
            label="Disparity (sum of variances)",   #y axis label
            cex.ts=1,            #Size of time scale labels.
            cex.age =0.6,        #Size of labels in axis            
            ts.col=F,            #Whether to include colours in the timescale.
            ts.width=0.2)    +    #thickness of the time scale
  box(lwd=0.8)


#Save at 6.8 x 4 inches, then reduce to 5.5 cm tall



##########################

# Test disparity using t-test
#---------------------------#


#https://rdrr.io/cran/dispRity/man/test.dispRity.html: TESTING disPrity hypotheses


ttest <-test.dispRity(timedisparity, 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "bonferroni") #corrected
ttest
write.csv(ttest, "ttest-results-13time-corrected.csv")



ttest <-test.dispRity(timedisparity, 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "none") #uncorrected

ttest
write.csv(ttest, "ttest-results-13time-uncorr.csv")



##########################

# Test disparity using t-test
#---------------------------#


#https://rdrr.io/cran/dispRity/man/test.dispRity.html: TESTING disPrity hypotheses


ttest <-test.dispRity(timedisparity, 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "bonferroni") #corrected
ttest
write.csv(ttest, "ttest-results-13time-corrected.csv")



ttest <-test.dispRity(timedisparity, 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "none") #uncorrected

ttest
write.csv(ttest, "ttest-results-13time-uncorr.csv")



########################
########################


# DISPARITY BETWEEN 11 TIME BINS
#---------------------------------#

##some taxa span various time-bins
##Read the pc.scores dataset with duplicated taxa with different names (Taxon / Taxon_d)
#MAKE SURE THIS FILE IS UPDATED:
#pcscores-dup.csv is a file created by duplicating some taxa in the original pc.scores
pc.scores <- read.csv("pcscores-dup.csv", 
                      row.names=1, 
                      header=T) 

pc.scores <- as.matrix(pc.scores) 
pc.scores



#Create list of taxa per time bin (column 19)
#This list needs to come from the file with taxa duplicated with different names (Taxon / Taxon_d)
#so that it matches the same names as in the new pc.scores
dframe2<-read.csv("FossilAqReptiles_data_dup.csv", 
                  row.names = 1,
                  header=T)
dframe2

#Exclude incomplete rows
dframe2<-na.exclude(dframe2)
dframe2




#Assess number of bins/groups there are in Time_5 
#(there should be 12, including the category "DUP", not to be accounted for)
groups_bins <- with( dframe2,
                     list(
                       groups = unique(Group_1),
                       bins   = unique(Time_5)
                     ))

sapply(groups_bins, length)


IndOle_list <- row.names(dframe2)[dframe2[,22] == 'Ind-Ole']  # change names of groups to match your data
IndOle_list
AniLad_list <- row.names(dframe2)[dframe2[,22] == 'Ani-Lad']  # change names of groups to match your data
AniLad_list
Car_list <- row.names(dframe2)[dframe2[,22] == 'Car']  # change names of groups to match your data
Car_list
NorRha_list <- row.names(dframe2)[dframe2[,22] == 'Nor-Rha']  # change names of groups to match your data
NorRha_list
HetSin_list <- row.names(dframe2)[dframe2[,22] == 'Het-Sin']  # change names of groups to match your data
HetSin_list
PliToa_list <- row.names(dframe2)[dframe2[,22] == 'Pli-Toa']  # change names of groups to match your data
PliToa_list
AalCall_list <- row.names(dframe2)[dframe2[,22] == 'Aal-Call']  # change names of groups to match your data
AalCall_list
OxfTith_list <- row.names(dframe2)[dframe2[,22] == 'Oxf-Tith']  # change names of groups to match your data
OxfTith_list
BerAlb_list <- row.names(dframe2)[dframe2[,22] == 'Ber-Alb']  # change names of groups to match your data
BerAlb_list
CenCon_list <- row.names(dframe2)[dframe2[,22] == 'Cen-Con']  # change names of groups to match your data
CenCon_list
SanMaas_list <- row.names(dframe2)[dframe2[,22] == 'San-Maas']  # change names of groups to match your data
SanMaas_list



#Create data grouped by time bins
time.groups <- list(IndOle_list,AniLad_list,Car_list,NorRha_list,HetSin_list,PliToa_list,AalCall_list,OxfTith_list,
                    BerAlb_list,CenCon_list,SanMaas_list)
names(time.groups) <- c("Ind-Ole", "Ani-Lad","Car","Nor-Rha", "Het-Sin", "Pli-Toa", "Aal-Call", "Oxf-Tith",
                        "Ber-Alb","Cen-Con", "San-Maas")

time.groups

timedisparity <- dispRity.per.group(pc.scores, 
                                    group = time.groups, 
                                    metric = c(sum, variances))

timedisparity
summary(timedisparity)

#Nice plot:
plot.new ()
plot(timedisparity, ylab="Disparity (sum of variances)", 
     col=c("#7C4EF9", "#A788FC", "#D6CDF9","#D6CDF9","#42A5E8","#42A5E8","#82E3F9","#BEF3FC","#7FD33D","#C5F99D","#C5F99D"), 
     par(mar=c(7,5,1,1)), #sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. The default is 4
     cex.axis=1, # size of font in axis labels
     cex.lab=1.2, # size of font in axis title
     ylim=c(0,30), # limits of axis
     xaxt = "n",         #this eliminates the x labels (because there is no way to rotate them using plot function)
     xlab='')            #this eliminates the x labels

## Draw the x-axis labels
text(x = 1:11,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c("Ind-Ole", "Ani-Lad","Car","Nor-Rha", "Het-Sin", "Pli-Toa", "Aal-Call", "Oxf-Tith",
                "Ber-Alb","Cen-Con", "San-Maas"), 
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

####################
####################

#1-Nice plot WITH NUMBERS:
plot.new ()

plot(timedisparity, ylab="Disparity (sum of variances)", 
     col=c("#7C4EF9", "#A788FC", "#D6CDF9","#42A5E8","#82E3F9","#BEF3FC","#7FD33D","#C5F99D"),
     at = c(250, 242.1, 232, 214.1, 196.05, 182.5, 168.8, 154.3, 122.8, 93.4, 76.15),  #Assigns numerical numbers to the categoric time bins. Ben recommended at = c(252, 247, 220, 190, 170, 155, 125, 80) 
     box(lwd=0.6),
     par(mar=c(7,5,1,1)), #sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. The default is 4
     cex.axis=1, # size of font in axis labels
     cex.lab=1.2, # size of font in axis title
     ylim=c(0,28), # limits of axis
     xaxt = "n",       #this eliminates the x labels (because there is no way to rotate them using plot function)
     xlab="Time (Ma)")            #this eliminates the x labels



## 2- using Geoscale plots
#https://cran.r-project.org/web/packages/geoscale/geoscale.pdf
#https://rdrr.io/cran/geoscale/man/geoscaleBox.html

require(geoscale)

## Converting the data into a list
data_obs <- extract.dispRity(timedisparity, observed = TRUE)
data_distribution <- extract.dispRity(timedisparity, observed = FALSE)
## Removing one list level
data_distribution <- unlist(data_distribution, recursive = FALSE)
data_obs <- as.vector(data_obs)


#assigning numerical values to the ages
ages <- c(250, 242.1, 232, 214.1, 196.05, 182.5, 168.8, 154.3, 122.8, 93.4, 76.15)
ages


## Plotting the results distribution: THIS PIECE OF CODE IS A BIT BUGGY
geoscaleBox(data_distribution, ages, 
            boxes = "Stage", 
            data.lim = c(2,30),  #limits of the y axis
            age.lim= c(251,66),  #limits of the x axis
            units=c("Period"),   #time units shown in geological scale
            box.width=3,         #width of the boxplot in millions of years
            tick.scale = 50,     #Time scale resolution.
            label="Disparity (sum of variances)",   #y axis label
            cex.ts=1,            #Size of time scale labels.
            cex.age =0.6,        #Size of labels in axis            
            ts.col=F,            #Whether to include colours in the timescale.
            ts.width=0.2)    +    #thickness of the time scale
  box(lwd=0.8)


#Save at 6.8 x 4 inches, then reduce to 5.5 cm tall


##########################

# Test disparity using t-test
#---------------------------#


#https://rdrr.io/cran/dispRity/man/test.dispRity.html: TESTING disPrity hypotheses


ttest <-test.dispRity(timedisparity, 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "bonferroni") #corrected
ttest
write.csv(ttest, "ttest-results-11time-corrected.csv")



ttest <-test.dispRity(timedisparity, 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "none") #uncorrected

ttest
write.csv(ttest, "ttest-results-11time-uncorr.csv")



#PMANOVA: PAIRWISE NON PARAMETRIC MANOVA BETWEEN CLADES
#----------------------------------------------------#

# Make a csv file of the pcscores (the PC coordinates of each taxa) to which you can add the clade and time-bin group in new columns
# this has to be done at the begining (right after calculating the pc.scores without any duplicates)

# reading the pscore file (no duplicates for PMANOVA between clades)
pcscore.data <- read.csv("FossilAqreptiles_Output.csv", row.names=1, header=T)
scores <- pcscore.data[,1:5] # just the pc scores


pairwise.adonis <- function(x, factors, 
                            sim.method = 'euclidean', 
                            p.adjust.m ='bonferroni')     #'fdr'another adjustment method
                            {
                            library(vegan)
                            co = combn(unique(factors),2)
                            pairs = c()
                            F.Model =c()
                            R2 = c()
                            p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , 
                  method =sim.method, 
                  permutations = 999);
                  pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
                  F.Model =c(F.Model,ad$aov.tab[1,4]);
                  R2 = c(R2,ad$aov.tab[1,5]);
                  p.value = c(p.value,ad$aov.tab[1,6])
                  }
                  p.adjusted = p.adjust(p.value,
                                        method=p.adjust.m)
                  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
                  return(pairw.res)
                  }


# permanova comparing morphospace BY CLADES
cladescore.groups <- as.vector(pcscore.data$Group_1)      # groups as vector
names(cladescore.groups) <- row.names(scores)

PERMANOVA.results <- pairwise.adonis(scores, 
                                     factors=cladescore.groups, 
                                     sim.method = 'euclidean', #Similarity method from vegdist, default is 'bray'
                                     p.adjust.m ='bonferroni') # p-values are multiplied by the number of comparisons     

#another adjustment method is "fdr" that controls the false discovery rate

PERMANOVA.results                                         # view results table
write.csv(PERMANOVA.results, "Permanova by clades.csv")   # save results as a csv


# if p.adjusted is less than 0.05, then the groups have significantly different centroids (positions in morphospace)
# if p.adjusted values are close to significance threshold (0.05) - increase permutations = 999 in the function to permutations = 9999 (it will take longer to run)



#####################################


# For the next PMANOVA analyses use the file with duplicated taxa that expands various time-bins

#The file: FossilAqreptiles_dup_Output has been created by duplicating PCA scores for taxa that is in various time bins
#the duplicates (Taxon_d) have been created manually, with the same scores than the respective Taxon



# PERMANOVA FOR 11 TIME BINS
#-------------------------#

#load data and, before doing the analysis and FILTER THE DATA THAT HAVE <NA> in column "Time_2"
pcscore.data <- read.csv("FossilAqreptiles_dup_Output.csv", 
                         row.names=1, 
                         header=T)

pcscore.data  <- pcscore.data[!is.na(pcscore.data$Time_5), ] #Before doing the analysis, FILTER THE DATA THAT HAVE <NA> in column "Time_2"
pcscore.data
scores <- pcscore.data[,1:5] # just the pc scores


timescore.groups <- as.vector(pcscore.data$Time_5)                  # groups as vector
timescore.groups
names(timescore.groups) <- row.names(scores)

PERMANOVA.results <- pairwise.adonis(scores, 
                                     factors=timescore.groups, 
                                     sim.method = 'euclidean', 
                                     p.adjust.m ='bonferroni') # p-values are multiplied by the number of comparisons     

#another adjustment method is "fdr" that controls the false discovery rate

PERMANOVA.results                                                   # view results table
write.csv(PERMANOVA.results, "Permanova by 11 time bins.csv")          # save results as a csv

# if p.adjusted is less than 0.05, then the groups have significantly different centroids (positions in morphospace)
# if p.adjusted values are close to significance threshold (0.05) - increase permutations = 999 in the function to permutations = 9999 (it will take longer to run)




# PERMANOVA FOR 13 TIME BINS
#-------------------------#

#load data and, before doing the analysis and FILTER THE DATA THAT HAVE <NA> in column "Time_2"
pcscore.data <- read.csv("FossilAqreptiles_dup_Output.csv", 
                         row.names=1, 
                         header=T)

pcscore.data  <- pcscore.data[!is.na(pcscore.data$Time_4), ] #Before doing the analysis, FILTER THE DATA THAT HAVE <NA> in column "Time_2"
pcscore.data
scores <- pcscore.data[,1:5] # just the pc scores


timescore.groups <- as.vector(pcscore.data$Time_4)                  # groups as vector
timescore.groups
names(timescore.groups) <- row.names(scores)

PERMANOVA.results <- pairwise.adonis(scores, 
                                     factors=timescore.groups, 
                                     sim.method = 'euclidean', 
                                     p.adjust.m ='bonferroni') # p-values are multiplied by the number of comparisons     

#another adjustment method is "fdr" that controls the false discovery rate

PERMANOVA.results                                                   # view results table
write.csv(PERMANOVA.results, "Permanova by 13 time bins.csv")          # save results as a csv

# if p.adjusted is less than 0.05, then the groups have significantly different centroids (positions in morphospace)
# if p.adjusted values are close to significance threshold (0.05) - increase permutations = 999 in the function to permutations = 9999 (it will take longer to run)






# PERMANOVA FOR 10Mya VERSUS the rest
#----------------------------#

#load data again #Before doing the analysis and FILTER THE DATA THAT HAVE <NA> in column "Time_3"
pcscore.data <- read.csv("FossilAqreptiles_dup_Output.csv", 
                         row.names=1, 
                         header=T)
pcscore.data  <- pcscore.data[!is.na(pcscore.data$Time_3), ]
pcscore.data
scores <- pcscore.data[,1:5] # just the pc scores

time10Myascore.groups <- as.vector(pcscore.data$Time_3)            # groups as vector
names(time10Myascore.groups) <- row.names(scores)

PERMANOVA.results <- pairwise.adonis(scores, 
                                     factors=time10Myascore.groups, 
                                     sim.method = 'euclidean', 
                                     p.adjust.m ='bonferroni') # p-values are multiplied by the number of comparisons     

#another adjustment method is "fdr" that controls the false discovery rate

PERMANOVA.results                                                  # view results table
write.csv(PERMANOVA.results, "Permanova by time 10Mya-Rest.csv")   # save results as a csv

# if p.adjusted is less than 0.05, then the groups have significantly different centroids (positions in morphospace)
# if p.adjusted values are close to significance threshold (0.05) - increase permutations = 999 in the function to permutations = 9999 (it will take longer to run)




# PERMANOVA FOR BROAD TIME BINS
#--------------------------------------------------#

#load data again #Before doing the analysis and FILTER THE DATA THAT HAVE <NA> in column "Time"
pcscore.data <- read.csv("FossilAqreptiles_dup_Output.csv", 
                         row.names=1, 
                         header=T)

pcscore.data  <- pcscore.data[!is.na(pcscore.data$Time), ]
pcscore.data
scores <- pcscore.data[,1:5] # just the pc scores


broadtimescore.groups <- as.vector(pcscore.data$Time)               # groups as vector
names(broadtimescore.groups) <- row.names(scores)

PERMANOVA.results <- pairwise.adonis(scores, 
                                     factors=broadtimescore.groups, 
                                     sim.method = 'euclidean', 
                                     p.adjust.m ='bonferroni') # p-values are multiplied by the number of comparisons     

#another adjustment method is "fdr" that controls the false discovery rate

PERMANOVA.results                                                   # view results table
write.csv(PERMANOVA.results, "Permanova by time bins-broad.csv")    # save results as a csv

# if p.adjusted is less than 0.05, then the groups have significantly different centroids (positions in morphospace)
# if p.adjusted values are close to significance threshold (0.05) - increase permutations = 999 in the function to permutations = 9999 (it will take longer to run)





# PERMANOVA FOR 10Mya VERSUS the rest
#----------------------------#

#load data again #Before doing the analysis and FILTER THE DATA THAT HAVE <NA> in column "Time_3"
pcscore.data <- read.csv("FossilAqreptiles_dup_Output.csv", 
                         row.names=1, 
                         header=T)
pcscore.data  <- pcscore.data[!is.na(pcscore.data$Time_3), ]
pcscore.data
scores <- pcscore.data[,1:5] # just the pc scores

time10Myascore.groups <- as.vector(pcscore.data$Time_3)            # groups as vector
names(time10Myascore.groups) <- row.names(scores)

PERMANOVA.results <- pairwise.adonis(scores, 
                                     factors=time10Myascore.groups, 
                                     sim.method = 'euclidean', 
                                     p.adjust.m ='bonferroni') # p-values are multiplied by the number of comparisons     

#another adjustment method is "fdr" that controls the false discovery rate

PERMANOVA.results                                                  # view results table
write.csv(PERMANOVA.results, "Permanova by time 10Mya-Rest.csv")   # save results as a csv

# if p.adjusted is less than 0.05, then the groups have significantly different centroids (positions in morphospace)
# if p.adjusted values are close to significance threshold (0.05) - increase permutations = 999 in the function to permutations = 9999 (it will take longer to run)



