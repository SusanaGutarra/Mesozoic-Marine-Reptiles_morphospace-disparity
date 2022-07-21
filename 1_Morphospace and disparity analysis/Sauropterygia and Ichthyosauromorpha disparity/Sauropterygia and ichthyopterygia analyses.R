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
require(geoscale)

rm(list = ls())
#setwd("######")


#DISPARITY FOR SAUROPTERYGIA
#---------------------------#

# 7 TIME BINS
#------------#

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

#Create a subset of the data
dframe2<- subset(dframe2,Group_1 %in% c("Sauropterygia"))
dframe2



#Assess number of bins/groups there are in Time_2 
#(there should be 8, no Sauropterygia in Early Triassic. Including the category "DUP" not to be accounted for)
groups_bins <- with( dframe2,
                     list(
                       groups = unique(Group_1),
                       bins   = unique(Time_2)
                     ))

sapply(groups_bins, length)


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
time.groups <- list(MidTriassic_list,LateTriassic_list,EarlyJurassic_list,MidJurassic_list,LateJurassic_list,EarlyCretaceous_list,LateCretaceous_list)
names(time.groups) <- c("Middle Triassic","Late Triassic","Early Jurassic", "Middle Jurassic", "Late Jurassic", "Early Cretaceous", "Late Cretaceous")

time.groups

timedisparity <- dispRity.per.group(pc.scores, 
                                    group = time.groups, 
                                    metric = c(sum, variances))

timedisparity
summary(timedisparity)

#Nice plot:
plot.new ()
plot(timedisparity, ylab="Disparity (sum of variances)", 
     col=c("#A788FC", "#D6CDF9","#42A5E8","#82E3F9","#BEF3FC","#7FD33D","#C5F99D"), 
     box(lwd=0.6),
     par(mar=c(7,5,1,1)), #sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. The default is 4
     cex.axis=1, # size of font in axis labels
     cex.lab=1.2, # size of font in axis title
     ylim=c(0,15), # limits of axis
     xaxt = "n",         #this eliminates the x labels (because there is no way to rotate them using plot function)
     xlab='')            #this eliminates the x labels

## Draw the x-axis labels
text(x = 1:7,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c("Middle Triassic","Late Triassic","Early Jurassic", "Middle Jurassic", "Late Jurassic", "Early Cretaceous", "Late Cretaceous"), 
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.99,
     ## Increase label size.
     cex = 1)

#for presentation save 8x8

#for thesis figure save 5 x 4 inches

####################
####################

#1-Nice plot WITH NUMBERS:
plot.new ()

plot(timedisparity, ylab="Disparity (sum of variances)", 
     col=c("#A788FC", "#D6CDF9","#42A5E8","#82E3F9","#BEF3FC","#7FD33D","#C5F99D"),
     at = c(242.1, 219.1, 187.7, 168.8, 154.2, 125,  86),  #Assigns numerical numbers to the categoric time bins. Ben recommended at = c(252, 247, 220, 190, 170, 155, 125, 80) 
     box(lwd=0.2),
     par(mar=c(7,5,1,1)), #sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. The default is 4
     cex.axis=1.2, # size of font in axis labels
     cex.lab=1.3, # size of font in axis title
     ylim=c(0,10), # limits of axis
     xlim=c(251,66), # limits of axis
     xaxt = "n",       #this eliminates the x labels (because there is no way to rotate them using plot function)
     xlab="Time (Ma)")            #this eliminates the x labels




## 2 Nice plot with numbers using Geoscale plots
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
ages <- c(242.1, 219.1, 187.7, 168.8, 154.2, 125,  86)
ages


## Plotting the results distribution: THIS PIECE OF CODE IS A BIT BUGGY
geoscaleBox(data_distribution, ages, 
            boxes = "Epoch", 
            data.lim = c(0,10),  #limits of the y axis
            age.lim= c(251,66),  #limits of the x axis
            units=c("Period"),   #time units shown in geological scale
            box.width=3,         #width of the boxplot in millions of years
            tick.scale = 50,     #Time scale resolution.
            label="Disparity (sum of variances)",   #y axis label
            cex.ts=1,            #Size of time scale labels.
            cex.age =0.6,        #Size of labels in axis  
            cex.pt = 0.7,
            urotate=90,
            ts.col=F,            #Whether to include colours in the timescale.
            ts.width=0.2)    +    #thickness of the time scale
  box(lwd=0.8)


#Save at 5.8 x 3.5 inches



##########################

# Test disparity using t-test
#---------------------------#


#https://rdrr.io/cran/dispRity/man/test.dispRity.html: TESTING disPrity hypotheses


ttest <-test.dispRity(timedisparity, 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "bonferroni") #corrected
ttest
write.csv(ttest, "ttest-results-time-corrected-SAUROPTERYGIA.csv")


ttest <-test.dispRity(timedisparity, 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "none") #uncorrected

ttest
write.csv(ttest, "ttest-results-time-uncorr-SAUROPTERYGIA.csv")

###############################

# BROAD TIME BINS
#----------------#

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

#Create a subset of the data
dframe2<- subset(dframe2,Group_1 %in% c("Sauropterygia"))
dframe2



#Assess number of bins/groups there are in Time
#(there should be 4, no Sauropterygia in Early Triassic. Including the category "DUP" not to be accounted for)
groups_bins <- with( dframe2,
                     list(
                       groups = unique(Group_1),
                       bins   = unique(Time)
                     ))

sapply(groups_bins, length)


Triassic_list <- row.names(dframe2)[dframe2[,18] == 'Triassic']  # change names of groups to match your data
Triassic_list
Jurassic_list <- row.names(dframe2)[dframe2[,18] == 'Jurassic']  # change names of groups to match your data
Jurassic_list
Cretaceous_list <- row.names(dframe2)[dframe2[,18] == 'Cretaceous']  # change names of groups to match your data
Cretaceous_list


#Create data grouped by time bins
time.groups <- list(Triassic_list,Jurassic_list,Cretaceous_list)
names(time.groups) <- c("Triassic","Jurassic", "Cretaceous")

time.groups

timedisparity <- dispRity.per.group(pc.scores, 
                                    group = time.groups, 
                                    metric = c(sum, variances))

timedisparity
summary(timedisparity)

#Nice plot:
plot.new ()
plot(timedisparity, ylab="Disparity (sum of variances)", 
     #col=c("#A788FC","#42A5E8","#7FD33D","#C5F99D"), 
     box(lwd=0.6),
     par(mar=c(1,5,1,1)), #sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. The default is 4
     cex.axis=1.2, # size of font in axis labels
     cex.lab=1.3, # size of font in axis title
     ylim=c(0,10), # limits of axis
     xaxt = "n",         #this eliminates the x labels (because there is no way to rotate them using plot function)
     xlab='')            #this eliminates the x labels

## Draw the x-axis labels
text(x = 1:3,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c("Triassic","Jurassic", "Cretaceous"), 
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.99,
     ## Increase label size.
     cex = 1.2)

#for thesis figure save 2.73 x 4 inches

##########################

# Test disparity using t-test
#---------------------------#


#https://rdrr.io/cran/dispRity/man/test.dispRity.html: TESTING disPrity hypotheses


ttest <-test.dispRity(timedisparity, 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "bonferroni") #corrected
ttest
write.csv(ttest, "ttest-results-time-corrected-SAUROPTERYGIAbroad.csv")


ttest <-test.dispRity(timedisparity, 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "none") #uncorrected

ttest
write.csv(ttest, "ttest-results-time-uncorr-SAUROPTERYGIAbroad.csv")


###############################


#DISPARITY FOR ICHTHYOSAUROMORPHA
#-------------------------------#

# 5 TIME BINS
#------------#

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

#Create a subset of the data
dframe2<- subset(dframe2,Group_1 %in% c("Ichthyosauromorpha"))
dframe2



#Assess number of bins/groups there are in Time_2 
#(there should be 8, no Ichthyosauromorpha in Late Cretaceous. Including the category "DUP" not to be accounted for)
groups_bins <- with( dframe2,
                     list(
                       groups = unique(Group_1),
                       bins   = unique(Time_2)
                     ))

sapply(groups_bins, length)

EarlyTriassic_list <- row.names(dframe2)[dframe2[,19] == 'Early_Triassic']  # change names of groups to match your data
EarlyTriassic_list
MidTriassic_list <- row.names(dframe2)[dframe2[,19] == 'Mid_Triassic']  # change names of groups to match your data
MidTriassic_list
LateTriassic_list <- row.names(dframe2)[dframe2[,19] == 'Late_Triassic']  # change names of groups to match your data
LateTriassic_list
EarlyJurassic_list <- row.names(dframe2)[dframe2[,19] == 'Early_Jurassic']  # change names of groups to match your data
EarlyJurassic_list
EarlyJurassic_list <- row.names(dframe2)[dframe2[,19] == 'Early_Jurassic']  # change names of groups to match your data
EarlyJurassic_list
LateJurassic_list <- row.names(dframe2)[dframe2[,19] == 'Late_Jurassic']  # change names of groups to match your data
LateJurassic_list
EarlyCretaceous_list <- row.names(dframe2)[dframe2[,19] == 'Early_Cretaceous']  # change names of groups to match your data
EarlyCretaceous_list
LateJurassicEarlyCretaceous_list<- append(LateJurassic_list,EarlyCretaceous_list)
LateJurassicEarlyCretaceous_list


#Create data grouped by time bins

#with Late Jurassic and Early Cretaceous separate
#time.groups <- list(EarlyTriassic_list, MidTriassic_list,LateTriassic_list,EarlyJurassic_list,LateJurassic_list,EarlyCretaceous_list)
#names(time.groups) <- c("Early Triassic","Middle Triassic","Late Triassic","Early Jurassic", "Late Jurassic", "Early_Cretaceous")

#with Late Jurassic and Early Cretaceous together
time.groups <- list(EarlyTriassic_list, MidTriassic_list,LateTriassic_list,EarlyJurassic_list,LateJurassicEarlyCretaceous_list)
names(time.groups) <- c("Early Triassic","Middle Triassic","Late Triassic","Early Jurassic", "Late Jurassic-Early Cretaceous")

time.groups

timedisparity <- dispRity.per.group(pc.scores, 
                                    group = time.groups, 
                                    metric = c(sum, variances))

timedisparity
summary(timedisparity)

#Nice plot:
plot.new ()
plot(timedisparity, ylab="Disparity (sum of variances)", 
     col=c("#7C4EF9","#A788FC", "#D6CDF9","#42A5E8","#BEF3FC"), 
     box(lwd=0.6),
     par(mar=c(7,5,1,1)), #sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. The default is 4
     cex.axis=1, # size of font in axis labels
     cex.lab=1.2, # size of font in axis title
     ylim=c(0,6), # limits of axis
     xaxt = "n",         #this eliminates the x labels (because there is no way to rotate them using plot function)
     xlab='')            #this eliminates the x labels

## Draw the x-axis labels
text(x = 1:5,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c("Early Triassic","Middle Triassic","Late Triassic","Early Jurassic", "Late Jurassic-Early Cretaceous"), 
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.99,
     ## Increase label size.
     cex = 1)


#for thesis figure save 4 x 3 inches

####################
####################

## Nice plot with numbers using Geoscale plots

## Converting the data into a list
data_obs <- extract.dispRity(timedisparity, observed = TRUE)
data_distribution <- extract.dispRity(timedisparity, observed = FALSE)
## Removing one list level
data_distribution <- unlist(data_distribution, recursive = FALSE)
data_obs <- as.vector(data_obs)


#assigning numerical values to the ages

#with Late Jurassic separate
#ages <- c(250, 242.1, 219.1, 187.7, 158, 129)
#ages

#Late Jurassic and Early cretaceous together
ages <- c(250, 242.1, 219.1, 187.7, 132)
ages



## Plotting the results distribution: THIS PIECE OF CODE IS A BIT BUGGY
geoscaleBox(data_distribution, ages, 
            boxes = "Epoch", 
            data.lim = c(0,6),  #limits of the y axis
            age.lim= c(251,66),  #limits of the x axis
            units=c("Period"),   #time units shown in geological scale
            box.width=3,         #width of the boxplot in millions of years
            tick.scale = 50,     #Time scale resolution.
            label="Disparity (sum of variances)",   #y axis label
            cex.ts=1,            #Size of time scale labels.
            cex.age =0.6,        #Size of labels in axis    
            cex.pt = 0.7,
            ts.col=F,            #Whether to include colours in the timescale.
            ts.width=0.2)    +    #thickness of the time scale
  box(lwd=0.8)


#Save at 5.8 x 3.5 inches, then reduce to 5.5 cm tall



##########################

# Test disparity using t-test
#---------------------------#


#https://rdrr.io/cran/dispRity/man/test.dispRity.html: TESTING disPrity hypotheses


ttest <-test.dispRity(timedisparity, 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "bonferroni") #corrected
ttest
write.csv(ttest, "ttest-results-time-corrected-ICHTHYOSAUROMORPHA.csv")


ttest <-test.dispRity(timedisparity, 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "none") #uncorrected

ttest
write.csv(ttest, "ttest-results-time-uncorr-ICHTHYOSAUROMORPHA.csv")


# BROAD TIME BINS
#----------------#

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

#Create a subset of the data
dframe2<- subset(dframe2,Group_1 %in% c("Ichthyosauromorpha"))
dframe2



#Assess number of bins/groups there are in Time
#(there should be 4, no Sauropterygia in Early Triassic. Including the category "DUP" not to be accounted for)
groups_bins <- with( dframe2,
                     list(
                       groups = unique(Group_1),
                       bins   = unique(Time)
                     ))

sapply(groups_bins, length)


Triassic_list <- row.names(dframe2)[dframe2[,18] == 'Triassic']  # change names of groups to match your data
Triassic_list
Jurassic_list <- row.names(dframe2)[dframe2[,18] == 'Jurassic']  # change names of groups to match your data
Jurassic_list
Cretaceous_list <- row.names(dframe2)[dframe2[,18] == 'Cretaceous']  # change names of groups to match your data
Cretaceous_list
JurassicCretaceous_list <- append(Jurassic_list,Cretaceous_list)
JurassicCretaceous_list

#Create data grouped by time bins
time.groups <- list(Triassic_list,JurassicCretaceous_list)
names(time.groups) <- c("Triassic","Jurassic-Cretaceous")

time.groups

timedisparity <- dispRity.per.group(pc.scores, 
                                    group = time.groups, 
                                    metric = c(sum, variances))

timedisparity
summary(timedisparity)

#Nice plot:
plot.new ()
plot(timedisparity, ylab="Disparity (sum of variances)", 
     #col=c("#A788FC","#42A5E8","#7FD33D","#C5F99D"), 
     box(lwd=0.6),
     par(mar=c(2,5,1,1)), #sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. The default is 4
     cex.axis=1.2, # size of font in axis labels
     cex.lab=1.3, # size of font in axis title
     ylim=c(0,10), # limits of axis
     xaxt = "n",         #this eliminates the x labels (because there is no way to rotate them using plot function)
     xlab='')            #this eliminates the x labels

## Draw the x-axis labels
text(x = 1:2,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c("Triassic","Jurassic-Cretaceous"), 
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.99,
     ## Increase label size.
     cex = 1.2)

#for thesis figure save 2.38 x 4.19 inches

##########################

# Test disparity using t-test
#---------------------------#


#https://rdrr.io/cran/dispRity/man/test.dispRity.html: TESTING disPrity hypotheses


ttest <-test.dispRity(timedisparity, 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "bonferroni") #corrected
ttest
write.csv(ttest, "ttest-results-time-corrected-SAUROPTERYGIAbroad.csv")


ttest <-test.dispRity(timedisparity, 
                      test = t.test,
                      comparisons = "sequential",
                      correction = "none") #uncorrected

ttest
write.csv(ttest, "ttest-results-time-uncorr-SAUROPTERYGIAbroad.csv")


