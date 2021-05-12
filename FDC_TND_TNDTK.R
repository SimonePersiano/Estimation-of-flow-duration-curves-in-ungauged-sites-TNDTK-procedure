#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Authors: Dr. Simone Persiano, Dr. Alessio Pugliese
# Acknowledgments: Prof. Attilio Castellarin
# 
# The present code performs the following main operations:
# (1) Extraction of the POR-FDC from the daily streamflow series observed at the given sites
# (2) Computation of the total negative deviation (TND)
# (3) Application of total negative deviation top-kriging (TNDTK) for computing POR-FDCs at 
#     the given (ungauged) target sites
#
# The required inputs are:
# (1) Daily streamflow series observed for different river cross-sections
# (2) Shapefiles of the catchment boundaries associated with the gauged river cross-sections
# (3) Shapefiles of the catchment boundaries associated with the target ungauged river cross-sections 
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Importing required packacges
library(rtop)
library(rgeos)
library(rgdal)
#
# Setting path of the working directory
path <- "/media/simone/dati/ALFFA_project_Eurac_Bolzano/ScriptsZenodo"
setwd(path)
#
# Reading shapefile of the gauged catchment boundaries
GaugedCatchments <- readOGR("GaugedCatchments/GaugedCatchments.shp",stringsAsFactors=FALSE)
View(as.data.frame(GaugedCatchments)) # Visualizing attribute table
#
# Reading shapefile of the gauged stations (sites)
GaugedStations <- readOGR("GaugedStations/GaugedStations.shp",stringsAsFactors=FALSE)
View(as.data.frame(GaugedStations)) # Visualizing attribute table
#
# Reading shapefile of catchment boundaries for ungauged target sites
UngaugedCatchments <- readOGR("UngaugedCatchments/UngaugedCatchments.shp",stringsAsFactors=FALSE)
View(as.data.frame(UngaugedCatchments)) # Visualizing attribute table
#
# Plotting
plot(GaugedCatchments) # Visualizing watershed boundaries
plot(UngaugedCatchments,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5),add=TRUE) # Visualizing watershed boundaries
plot(GaugedStations,pch=24,bg="green",cex=1.2,add=TRUE) # Adding station locations to the previous plot
#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
# (1) Extraction of the POR-FDC from the daily streamflow series observed at the given gauged sites 
# Creating a list with FDCs for the gauged sites
FDC <- list(Q=list(),d=list())
for (i in 1:length(GaugedCatchments$Cod)) # For each gauged site...
{
  code <- as.character(GaugedCatchments$Cod[i])
  file_i <- list.files(path=paste(path,"/StreamflowData/",sep=""),
                       pattern=code)
  dataQ_i <- read.csv(paste(path,"/StreamflowData/",file_i,sep=""),
                      header=TRUE,stringsAsFactors=FALSE,
                      encoding="latin1",sep=",")
  # Building POR-FDC
  FDC$Q[[code]] <- sort(dataQ_i$Q,decreasing=TRUE)
  FDC$d[[code]] <- 1:length(FDC$Q[[code]])/length(FDC$Q[[code]]+1) # Weibull plotting position  
}
#
# Computing Mean Annual Flow (MAF) 
MAF <- c()
for (i in 1:length(GaugedCatchments$Cod))
{
  MAF[i] <- mean(FDC$Q[[i]])
}
GaugedCatchments$MAF <- MAF
#
# Re-computing catchment areas
GaugedCatchments$Area_km2 <- as.numeric(gArea(GaugedCatchments,byid=TRUE)/(1000*1000))
#
#
# Function for plotting logarithmic sequences (useful for the plot) 
lseq <- function(from,to,length.out)
{
  exp(seq(log(from), log(to), length.out = length.out))
}
#
# Plot with the different dimensionless FDCs
#
plot(FDC$d[[1]],FDC$Q[[1]]/mean(FDC$Q[[1]]), 
     type="l",log="y",
     ylim=c(0.001,100),xlim=c(0,1),
     xlab="Duration [-]", 
     ylab="Dimensionless streamflow [-]",
     col="black",axes=FALSE)
grid(NULL,NULL, lty = 3, col = "lightgray")
abline(v=seq(0,1,0.1), lty = 3, col = "lightgray")
abline(h=lseq(0.01,100,5), lty = 3, col = "lightgray")
axis(1,at=seq(0,1,0.1),labels = paste(seq(0,1,0.1)))
axis(2,at=lseq(0.01,100,5),labels = paste(lseq(0.01,100,5)))
for (i in 2:length(GaugedCatchments$Cod))
{
  lines(FDC$d[[i]],FDC$Q[[i]]/mean(FDC$Q[[i]]),col="blue")
}
#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
# (2) Computation of the total negative deviation (TND)
# TND is computed using a Gaussian transformation of the duration axis
#
library(pracma) # Needed for the tnd function
source("tnd.R") # Running function for computing TND
#
years <- rep(10,length(GaugedCatchments$Cod))
min.year <- min(years) # Minimum record length for maximum duration to be used in TND computations
maxd <- min.year*365/(min.year*365+1) # lower bound for maximum durations
#
TND <- c()
for (i in 1:length(GaugedCatchments$Cod))
{
  TND[i] <- fdc.tnd(FDC$Q[[i]]/MAF[i],norm=TRUE,maxd=maxd,logq=FALSE)
}
#
# Adding a column (field) with the computed TND values
GaugedCatchments$TND <- TND
#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
# (3) Application of total negative deviation top-kriging (TNDTK) for computing the POR-FDC at 
#     a given (ungauged) target site
#
# (a) Scaling relationship between MAF and drainage area (for gauged catchments) 
plot(GaugedCatchments$Area_km2,GaugedCatchments$MAF,log="xy",
     xlab=expression("Drainage area [km"^2*"]"),
     ylab=expression("MAF [m"^3*"/s]"))
#
loglin.mod <- lm(log(GaugedCatchments$MAF)~log(GaugedCatchments$Area_km2))
c1 <- exp(loglin.mod$coefficients[1])
c2 <- loglin.mod$coefficients[2]
curve(c1*x^c2,add=TRUE,col="black")
#
# * * * 
# (b) Top-kriging prediction of TND
set.seed(1)
vic <- 6 # Neighbourhood (no. of considered catchments)
# Create an rtop object
rtopObj <- createRtopObject(observations=GaugedCatchments,
                            predictionLocations=UngaugedCatchments,
                            formulaString=TND~1,       # Variable to krige
                            params=list(gDist=TRUE,      # How to compute distances between areas
                                        rresol=500,      # Catchment area discretization
                                        nmax=vic,        # No. of neighbours
                                        wlim=1,          # Upper limit for the norm of the weights in kriging
                                        debug.level=0,   # No additional outputs required 
                                        partialOverlap=TRUE))  # For inaccurate catchment boundaries
# Build the empirical variogram
rtopObj <- rtopVariogram(rtopObj)     # Create variogram for the rtop object
rtopObj <- rtopFitVariogram(rtopObj)  # Fit theoretical variogram 
rtopObj <- checkVario(rtopObj)        # Check variogram fitting 
# Perform the predicton
rtopObj <- rtopKrige(rtopObj,wret=TRUE)
#
TND.pred <- rtopObj$predictions$var1.pred  # TND values estimated with TK
weights <- rtopObj$weight                  # Geostatistical TK weights
#
# Adding a column (field) with the predicted TND values
UngaugedCatchments$TND.pred <- TND.pred
#
# * * * 
# (c) Top-kriging prediction of MAF
GaugedCatchments$obs <- GaugedCatchments$MAF/(GaugedCatchments$Area_km2^c2)
# Create an rtop object
rtopObj.MAF <- createRtopObject(observations=GaugedCatchments,
                                predictionLocations=UngaugedCatchments,
                                formulaString=obs~1,       # Variable to krige
                                params=list(gDist=TRUE,      # How to compute distances between areas
                                            rresol=500,      # Catchment area discretization
                                            nmax=vic,        # No. of neighbours
                                            wlim=1,          # Upper limit for the norm of the weights in kriging
                                            debug.level=0,   # No additional outputs required 
                                            partialOverlap=TRUE))  # For inaccurate catchment boundaries
# Build the empirical variogram
rtopObj.MAF <- rtopVariogram(rtopObj.MAF)     # Create variogram for the rtop object
rtopObj.MAF <- rtopFitVariogram(rtopObj.MAF)  # Fit theoretical variogram 
rtopObj.MAF <- checkVario(rtopObj.MAF)        # Check variogram fitting 
# Perform the predicton
rtopObj.MAF <- rtopKrige(rtopObj.MAF)        
#
# MAF values estimated with TK
MAF.pred <- rtopObj.MAF$predictions$var1.pred*(UngaugedCatchments$Are_km2^c2)
##
UngaugedCatchments$MAF.pred <- MAF.pred
#
# * * * 
# (d) Computation of FDCs for ungauged catchments
#
# Resampling observed FDCs in view of the estimation for target sites
source("resample_FDC.R") # Running function for resampling FDCs
np <- 100 # number of sampling points (for the FDC)
sam <- pnorm(seq(qnorm(1-maxd),qnorm(maxd),(qnorm(maxd)-qnorm(1-maxd))/(np-1))) # sampling points
duration_percentage <- round(sam*100,2)
#
# Empirical FDCs (dimensionless) for gauged catchments
FDC.obs.dimless <- matrix(NA,nrow=np,ncol=length(GaugedCatchments$Cod)) 
colnames(FDC.obs.dimless) <- GaugedCatchments$Cod
rownames(FDC.obs.dimless) <- paste(duration_percentage)
for(j in 1:length(GaugedCatchments$Cod)) # Resampling...
{
  resamp <- resample.FDC(FDC$Q[[j]]/GaugedCatchments$MAF[j],sam=sam,norm=FALSE)  
  FDC.obs.dimless[,j] <- resamp$y
}
#
# Final predictions
# Predicted FDCs (dimensionless)
FDC.pred.dimless <- FDC.obs.dimless%*%t(weights)  ### estimates
colnames(FDC.pred.dimless) <- UngaugedCatchments$Locatin
rownames(FDC.pred.dimless) <- paste(duration_percentage)
#
# Predicted FDCs (dimensional)
FDC.pred <- t(t(FDC.pred.dimless)*MAF.pred)
#
# 
# Plotting final results:
# * Blue lines represent observed dimensional FDCs
# * Red lines represent predicted dimensional FDCs (at ungauged target sites)
#
plot(FDC$d[[1]],FDC$Q[[1]]/mean(FDC$Q[[1]]), 
     type="l",log="y",
     ylim=c(0.01,1000),xlim=c(0,1),
     xlab="Duration [-]", 
     ylab=expression("Dimensional streamflow [m"^3*"/s]"),
     col="black",axes=FALSE)
grid(NULL,NULL, lty = 3, col = "lightgray")
abline(v=seq(0,1,0.1), lty = 3, col = "lightgray")
abline(h=lseq(0.01,1000,6), lty = 3, col = "lightgray")
axis(1,at=seq(0,1,0.1),labels = paste(seq(0,1,0.1)))
axis(2,at=lseq(0.01,1000,6),labels = paste(lseq(0.01,1000,6)))
for (i in 2:length(GaugedCatchments$Cod))
{
  lines(FDC$d[[i]],FDC$Q[[i]],col="blue",lwd=1.5)
}
for (j in 1:length(UngaugedCatchments$Locatin))
{
  lines(as.numeric(rownames(FDC.pred))/100,as.numeric(FDC.pred[,j]),col="red",lwd=3)
}
