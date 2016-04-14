library(sp)
library(rgdal)
library(raster)
source("C:\\Users\\johnwithrow\\Desktop\\Projects\\R_Libraries\\spatialRoutinesInR.R")

plgOuterFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF"
plgOuterFnShort <- "BHNFa"
numSamplePoints <- 750
minDistance <- 2500
propInLocations <- 0.5
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\BHNF_bYs_2007_2010_aYs_2011_2012_jDs_200_275_ZScore_Analysis2Zs_mean.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2012.shp"
plgDisturbanceFnShort <- "IDS2012"

fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
