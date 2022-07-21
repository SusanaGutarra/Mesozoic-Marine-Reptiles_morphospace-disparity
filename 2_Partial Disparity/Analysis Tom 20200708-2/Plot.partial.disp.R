# makes function for nice scale bar - made by Ben
xAxisPeriods <- function (cex = 1.0) {
  # Plots boxes with geological period labels along the x-axis.
  # Only for periods of the Mesozoic so far.
  #
  # Args:
  #   cex: character expansion of labels
  #
  # Returns:
  #   In plots, an x-axis with boxes for each geological period, labels, and
  #   time on the x-axis.
  # get x-axis extent
  axis.bottom <- par()$usr[3] - 0.1 * diff(par()$usr[3:4])
  axis.top <- par()$usr[3]
  # plot x-axis
  axis(1, pos = axis.bottom)
  # plot period boxes
  # Triassic
  rect(min(par()$usr[1], 252.17), axis.bottom, 201.3, axis.top,
       col = "white", border = TRUE, xpd = NA)
  text(mean(c(min(par()$usr[1], 252.17), 201.3)),
       mean(c(axis.bottom, axis.top)),
       "TRIASSIC",
       adj = c(0.5, 0.5), cex = cex, xpd = NA)
  # Jurassic
  rect(201.3, axis.bottom, 145, axis.top,
       col = "white", border = TRUE, xpd = NA)
  text(mean(c(201.3, 145)),
       mean(c(axis.bottom, axis.top)),
       "JURASSIC",
       adj = c(0.5, 0.5), cex = cex, xpd = NA)
  # Cretaceous
  rect(145, axis.bottom, max(par()$usr[2], 66.0), axis.top,
       col = "white", border = TRUE, xpd = NA)
  text(mean(c(145, max(par()$usr[2], 66.0))),
       mean(c(axis.bottom, axis.top)),
       "CRETACEOUS",
       adj = c(0.5, 0.5), cex = cex, xpd = NA)
  # tidy up box around
  rect(par()$usr[1], axis.bottom, par()$usr[2], axis.top,
       col = NA, border = TRUE, xpd = NA)
  # tidy up box around
  rect(par()$usr[1], axis.bottom, par()$usr[2], axis.top,
       col = NA, border = TRUE, xpd = NA)
}



# plot partial disparity data from MATLAB MDA analysis
library(scales)

pd.data <- read.table ("PD_RESULTS_R.txt", header = T, row.names = 1)
pd.data 

midpoints <-  c(251, 242, 220, 190, 168, 155, 125, 83)

nbins <- length(midpoints)
ngroups <- ncol(pd.data)

#~~~~~~~~~~~
# Begin extracting data
#~~~~~~~~~~~

#The calculations below will require editing
#You will have to change the number of time bins (rows, 16 here)
#And the number of groups (columns, 6 here) to suit your data

relative.disparity = array (dim = c(nbins, ngroups)) # This creates an empty array of the proper proportions
colnames(pd.data)

pd.data[is.na(pd.data)] <- 0

#the below script claculate the relative disparities/diversities ready  for plotting
for (i in 1:nrow(pd.data))
  
{
  
  relative.disparity [i, ] = c(pd.data[i, "Ichthyosauromorpha"], pd.data[i, "Sauropterygia"], pd.data[i, "Thalattosauria"], pd.data[i, "Other"], pd.data[i, "Pantestudines"], pd.data[i, "Thalattosuchia"], pd.data[i, "Rhynchocephalia"], pd.data[i, "Mosasauroidea"])
  
  relative.disparity[i, ] = cumsum(relative.disparity [i, ])
  
  relative.disparity [i, ] = relative.disparity [i, ]/relative.disparity [i, ngroups]
  
}

colnames(pd.data)
colnames(relative.disparity) <- colnames(pd.data)

# Set up plot area

x.limits <- c(252, 66)
y.limits <- c(0.02, .98) ## will leave tiny gaps at top and bottom - for aesthetics
plot (midpoints, y = c(rep(0.8,length(midpoints))), xlim = x.limits, ylim = y.limits, xaxt='n', col = "transparent", xlab = "", ylab = paste("Partial disparity" , sep=""))
abline(v=201.3, lwd=2, col="gray70")
abline(v=145, lwd=2, col="gray70")
abline(v=66, lwd=2, col="gray70")
box(lwd=0.6)

colnames(pd.data)

coordinates.x <- c(midpoints[1], midpoints, midpoints [nbins])
coordinates.y <- c(0, relative.disparity[, "Mosasauroidea"], 0)
polygon (coordinates.x, coordinates.y, col = "#FF8300", border = NA)

coordinates.x <- c(midpoints[1], midpoints, midpoints [nbins])
coordinates.y <- c(0, relative.disparity[, "Rhynchocephalia"], 0)
polygon (coordinates.x, coordinates.y, col = "#c061ff", border = NA)

coordinates.x <- c(midpoints[1], midpoints, midpoints [nbins])
coordinates.y <- c(0, relative.disparity[, "Thalattosuchia"], 0)
polygon (coordinates.x, coordinates.y, col = "blue", border = NA)

coordinates.x <- c(midpoints[1], midpoints, midpoints [nbins])
coordinates.y <- c(0, relative.disparity[, "Pantestudines"], 0)
polygon (coordinates.x, coordinates.y, col = "#009ccc", border = NA)

coordinates.x <- c(midpoints[1], midpoints, midpoints [nbins])
coordinates.y <- c(0, relative.disparity[, "Other"], 0)
polygon (coordinates.x, coordinates.y, col = "gray70", border = NA)

coordinates.x <- c(midpoints[1], midpoints, midpoints [nbins])
coordinates.y <- c(0, relative.disparity[, "Thalattosauria"], 0)
polygon (coordinates.x, coordinates.y, col = "#FFD300", border = NA)

coordinates.x <- c(midpoints[1], midpoints,midpoints [nbins])
coordinates.y <- c(0, relative.disparity[, "Sauropterygia"], 0)
polygon (coordinates.x, coordinates.y, col = "#19D27E", border = NA)

coordinates.x <- c(midpoints[1], midpoints, midpoints [nbins])
coordinates.y <- c(0, relative.disparity[, "Ichthyosauromorpha"], 0)
polygon (coordinates.x, coordinates.y, col = "#FF2626", border = NA)


lines(midpoints, relative.disparity[,1], lty=1, lwd=0.8,col = "white")
lines(midpoints, relative.disparity[,2], lty=1, lwd=0.8,col = "white")
lines(midpoints, relative.disparity[,3], lty=1, lwd=0.8,col = "white")
lines(midpoints, relative.disparity[,4], lty=1, lwd=0.8,col = "white")
lines(midpoints, relative.disparity[,5], lty=1, lwd=0.8,col = "white")
lines(midpoints, relative.disparity[,6], lty=1, lwd=0.8,col = "white")
lines(midpoints, relative.disparity[,7], lty=1, lwd=0.8,col = "white")
lines(midpoints, relative.disparity[,8], lty=1, lwd=0.8,col = "white")

xAxisPeriods()

#segments(201.3, 0, 201.3, 2, col=rgb (0, 0, 0, 0.25),lwd=2)
#segments(145, 0, 145, 2, col=rgb (0, 0, 0, 0.25),lwd=2)
#segments(66, 0, 66, 2, col=rgb (0, 0, 0, 0.25),lwd=2)

dev.copy(pdf,"Partial_disparity_new.pdf", width=6, height=3.5)
dev.off()

plot.new ()
