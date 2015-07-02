# Date:         January 12, 2014
# Author:       Johanna Wren 
# Email:        jwren@hawaii.edu
# Purpose:      Script to test simulation data using daily locations. Make sure currents are loading correctly and not plotting over land. 
#------------------------------------------------

#------------------------------------------------
# Setup working directory
#------------------------------------------------
mainDir <- 'DispersalModel/' ## Change to fit user settings ##
setwd(mainDir)

#------------------------------------------------
# Load libraries
#------------------------------------------------
library(marmap)

#------------------------------------------------
# Define Variables
#------------------------------------------------
### Change these to suit user (sites numbered lower than 250 are outside of geographic domanin of hawaii bathymetry file loaded below so won't show on map) ####
infile <- 'Dispersal_python_8kmHA_test_site_630.txt'  # file containing daily output data from python DispersalModel
data(hawaii)                 # Loads bathymetry for Hawaii

#------------------------------------------------
# Load data
#------------------------------------------------
Sites <- read.table(infile, header=F)  # loads the output file. Can take some time depending on the size of the file.
head(Sites)
Sites[,4] <- Sites[,4]-360   # Changes longitude to degrees West (originally in deg E) for plotting
col.na
sites <- cbind(Sites[,4]-360, Sites[,5])  # makes subset of lat and lon, and converts lon to degrees West (from deg E)
head(sites)

#------------------------------------------------
# Plot bathymetry of MHI and red points for dispersal data
#------------------------------------------------
# # This plots all locations as points, useful to see if locations are located over land. Commetn out if looking at trajectoreis. 
# plot(hawaii, 
#      deep = c(-5000, -200, 0), shallow = c(-200, 0, 0),
#      col = c("grey", "blue", "black"), step = c(1000, 200, 1),
#      lty = c(1, 1, 1), lwd = c(0.6, 0.6, 1.2),
#      draw=c(FALSE, FALSE, FALSE))
# points(sites, pch = 21, col = "red", cex = 0.1)

#------------------------------------------------
# Plot individual particle trajectories and pauses after each
#------------------------------------------------
#### Change maxPLD to whatever was used in your run ########
maxPLD=15                    # enter max PLD so the script knows how many days one trajectory is. 
temp <- which(Sites[,2]==1)  # This pulls out currentday=1 (so day of dispersal) for each particle
for (i in 1:length(temp)) {  # loops through trajectories and plots them as a line, one at a time
  plot(hawaii, 
       deep = c(-5000, -200, 0), shallow = c(-200, 0, 0),
       col = c("grey", "blue", "black"), step = c(1000, 200, 1),
       lty = c(1, 1, 1), lwd = c(0.6, 0.6, 1.2),
       draw=c(FALSE, FALSE, FALSE))
  lines(Sites[temp[i]:(temp[i]+maxPLD-1),4], Sites[temp[i]:(temp[i]+maxPLD-1),5], col='red', pch=19, cex=20, lwd=3)
  readline(prompt = "Please press <Return> for next trajectory")
}
