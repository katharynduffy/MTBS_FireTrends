

require(raster)
require(rgdal)
require(maptools)
require(dplyr)
require(tidyr)
require(raptr)
require(geosphere)
require(rgeos)


#I'm going to comment on your ANALYSIS
setwd("I:/Research/Shares/maxnp_lab/RESEARCH/McEvoy2/National Forest Foundation/SW Carbon Methodology/R Projects/Data")


# This is a script for creating the sampling frame to be used in the SW Forest Restoration Carbon Methodology. This script, when complete, should be combined with the rest of the project script in 'SW MTBS DRAFT.R'


# The first criteria of the porject area is that it be within a specfic pyrome. Import the national pyrome shapefile, identify the pyrome of interest, and project the shapefile into the correct crs.

US_Pyrome <- readOGR(dsn = "E:/Spatial Data/National FSim/Pyromes", layer = "Pyrome_20150605")
head(US_Pyrome)
US_Pyrome[39,] # Find the pyrome of interest

# Pull out the pyrome of interest and create new SpatialPolygonsDataFrame
Pyrome <- US_Pyrome[US_Pyrome$NAME == "Arizona/New Mexico Mountains (Mogollon Rim)",]
plot(Pyrome)
crs(Pyrome) # Notice how Pyrome is in Lat/Long crs. We Want it in the appropriate UTM.

Pyrome <- spTransform(Pyrome, "+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs")
crs(Pyrome) # check the crs again and see that it has changed to the appropriate UTM projection

acres <- gArea(Pyrome)/4046.856
acres
# Now, identify all the areas within the pyrome that are in the appropriate fire regime group, national vegetation class, and vegetation condition class.

#Fire Regime Group
FRG <- raster("./Inputs/us_140frg.tif")
plot(FRG) # Notice how the extent downloaded from LANDFIRE is quite large. INCLUDE description of LANDFIRE FRG classes
FRG

#National Vegetation Class
NVC <- raster("./Inputs/us_200nvc.tif")
plot(NVC)
NVC

# Vegetation Condition Class
VCC <- raster("./Inputs/us_140vcc.tif")
plot(VCC)
VCC




# Function starts with the downloaded rasters (FRG, NVC, VCC), projects them into the appropriate crs, crops them to the smaller extent of the Pyrome, and extracts only those pixels with values inside the pyrome

#The input rasters need to be projected into the same crs as Pyrome, and they all need to be in the same extent as one another for later calculations. In order to do this, start by creating a blank raster set to the spatial extent of Pyrome and at the desired resolution.

Template <- blank.raster(Pyrome, 360) #Create a blank raster which has the same extent as Pyrome and the desired resolution
crs(Template)
crs(Pyrome) # notice how the crs of our blank raster, 'Template,' is blank.
crs(Template) <- crs(Pyrome) # set the crs of the blank raster to the same as Pyrome so that the blank  raster can serve as a template extent for later steps
crs(Template)
Proj <- crs(Template) # create an object which defines the projection to be used throughout the analysis

plot(Template)
plot(Pyrome, add = TRUE)


# Using an iterative function like the one described below, project all input rasters into the same crs, resolution, and extent
inputs <- list.files("./Inputs/", pattern = "[.]tif$", full.names = TRUE)
inputs

AlignInputs <- function(filename){
  r <- raster(filename)
  rpro <- projectRaster(r, Template, method = 'ngb', filename = paste0("./Scratch/rpro", basename(filename)), overwrite = TRUE)
  rmask <- mask(rpro, Pyrome, filename = paste0("./Scratch/rmask", basename(filename), overwrite =TRUE))
}

lapply(inputs, FUN = AlignInputs)



FRGmask <- raster("./Scratch/rmaskus_140frg.grd")
plot(FRGmask)
NVCmask <- raster("./Scratch/rmaskus_200nvc.grd")
plot(NVCmask)
VCCmask <- raster("./Scratch/rmaskus_140vcc.grd")
plot(VCCmask)



FRGsample <- FRGmask %in% c(1,3) # all pixels in Fire Regime Groups 1 and 3 will have a value of 1 in the resulting raster. Other fire regimes will have a value of 0
FRGsample
freq(FRGsample)

NVCsample <- NVCmask %in% c(6075, 6076)# All pixels in Southern Rocky Mountain Ponderosa Pine Forest and Southern Rocky Mountain Ponderosa Pine Woodland will have a value of 1 while all other vegetation types will have a value of 0
NVCsample
freq(NVCsample)

VCCsample <- VCCmask %in% 4:6
VCCsample
freq(VCCsample)

SampleStack <- stack(FRGsample, NVCsample, VCCsample)

SampleFrame <- overlay(SampleStack, fun=function(x,y,z){x*y*z}) # double check this. Should it be (x,y,z){x*y*z}
plot(SampleFrame)
plot(Pyrome, add = TRUE)
freq(SampleFrame)
# See that Sample Frame includes only values of 0 and 1. The area (acres) in the Sample Frame can be calculated by converting pixel number to area.


#END SAMPLE FRAME CREATION




## ------------------------------------------------------------------------------------------------------


#BEGIN MTBS ANALYSIS


# Extract and unzip downloaded MTBS severity mosaics
mtbs <- list.files("./MTBS Original/composite_data/MTBS_BSmosaics/", pattern = "[.]zip", recursive = TRUE, full.names = TRUE)
files

for (i in 1:length(mtbs)){
  unzip(mtbs[i])
}


#list raster files to be analyzed
years <- list.files("./AZ MTBS/", pattern = "[.]tif$", full.names = TRUE) # created a lookup table with year and raster file in advance


# create an extent that matches the extent of the Sample Frame. This will be used to bring all mtbs rasters into the same dimensions and extent.
extnt <- extent(SampleFrame)



## Create a function which will (1)project all MTBS rasters into the same crs as the Sample Frame; (2) extract only those pixels which are within the pyrome, and (3) create a new raster for each year in which only those pixels coincident with the Sample Frame retain their fire severity value. All other pixels - where eithe no fire occurred, or fire occurred but not in the Sample Frame - are assigned a value of zero.

MTBSalign <- function (filename){
  r <- raster(filename)
  rpro <- projectRaster(r, Template, method = 'ngb',filename = paste0("./Scratch/", basename(filename)))
  rmask <- mask(rpro, Pyrome, filename = paste0("./Scratch/mask", basename(filename)))
  overlay(rmask, SampleFrame, FUN = function(x,y){return(x*y)}, filename = "./Algebra Outputs/", basename(filename))
}

# Make sure you create a file for the results of the function. Example "Algebra Outputs" as written above.
#Apply that function to each MTBS raster in the list 'years.'
lapply(years, FUN = MTBSalign )




#############################################################################

#STEP 3. Calculate the 33-year average acres burned for each severity class.

Severity <- list.files("./Algebra Outputs/", pattern = "[.]tif$", full.names = TRUE)

SeverityAnalysis <- function(filename){
  r <- raster(filename)
  tbl <- dplyr::tbl_df(freq(r))
  return(tbl)
}

for (i in 1:length(Severity)){
  tbl <-SeverityAnalysis(filename = Severity[i])
  write.table(tbl, file = paste0("./SeverityTables/", tools::file_path_sans_ext(basename(Severity[i]))))
}


# Now merge all resulting tables into a single dataframe

file_list <- list.files("./SeverityTables/", full.names = TRUE) #these are tables to be joined
for (i in 1:length(file_list)){
  if (!exists("dataset")){
    dataset <- read.table(file_list[i])
  }
  if (exists("dataset")){
    temp_dataset <- read.table(file_list[i])
    dataset <- left_join(dataset, temp_dataset, by  = "value")
    rm(temp_dataset)
  }
}


colnames(dataset) <- c("Value", 1984, "Dup", seq(1985,2016,1))
head(dataset)
dataset$Dup <- NULL
dataset
dataset_t <- t(dataset)
dataset_t





SeverityAnalysis("./Algebra Outputs/Algebra Outputsmtbs_AZ_1984.tif")


mlow <- mean(Severity$Low)
mlow
mMod <- mean(Severity$Moderate)
mMod
mHigh <- mean(Severity$High)
mHigh

sdlow <- sd(Severity$Low)
sdMod <- sd(Severity$Moderate)
sdHigh <- sd(Severity$High)


MeanFire <- data.frame(
  Variable = c("Total", "Mean", "StdDev"),
  Low  = c(sum(Severity$Low), mlow, sdlow),
  Moderate = c(sum(Severity$Moderate), mMod, sdMod),
  High = c(sum(Severity$High),mHigh, sdHigh)
)

#Data frame has to be converted to matrix before plotting barplot
MeanFire.m <- as.matrix(MeanFire[2:4])
rownames(MeanFire.m) <- MeanFire[,1]
x <- barplot(MeanFire.m[2,], ylim = c(0,12000), main = "Average Annual Acres Burned by Severity Class", ylab = "Acres")
y <- as.matrix(MeanFire[2,2:4])
text(x,y+500,labels=as.character(y))
