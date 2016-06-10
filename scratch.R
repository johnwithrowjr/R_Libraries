library(sp)
library(rgdal)
library(raster)
source("C:\\Users\\johnwithrow\\Desktop\\Projects\\R_Libraries\\JRWBase.R")
source("C:\\Users\\johnwithrow\\Desktop\\Projects\\R_Libraries\\spatialRoutinesInR.R")



plgOuterFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\BHNF\\AOI.shp"
plgOuterFnShort <- "AOI"
numSamplePoints <- 100000
minDistance <- 0

rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2012_2013_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2013_Medium_Buffer.shp"
plgDisturbanceFnShort <- "IDS2013_Medium_Buffer"
analysis01 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS_Mgmt2013_Dissolve.shp"
plgDisturbanceFnShort <- "IDS_Mgmt2013_Dissolve"
analysis01b <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2012_2013_jDs_200_275_Landsat_NDVI_no7_TDDlm2.tif"
analysis02 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2013_2014_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2014_Medium_Buffer.shp"
plgDisturbanceFnShort <- "IDS2014_Medium_Buffer"
analysis03 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2013_2014_jDs_200_275_Landsat_NDVI_no7_TDDlm2.tif"
analysis04 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2013_2014_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2013_14_Medium_Buffer.shp"
plgDisturbanceFnShort <- "IDS2013_14_Medium_Buffer"
analysis03c <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2013_2014_jDs_200_275_Landsat_NDVI_no7_TDDlm2.tif"
analysis04c <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2014_2015_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2015_Medium_Buffer.shp"
plgDisturbanceFnShort <- "IDS2015_Medium_Buffer"
analysis05 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2014_2015_jDs_200_275_Landsat_NDVI_no7_TDDlm2.tif"
analysis06 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2014_2015_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2013_15_Medium_Buffer.shp"
plgDisturbanceFnShort <- "IDS2013_15_Medium_Buffer"
analysis05c <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2014_2015_jDs_200_275_Landsat_NDVI_no7_TDDlm2.tif"
analysis06c <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

plgOuterFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\BHNF\\AOI_StudyArea.shp"
plgOuterFnShort <- "AOI_StudyArea"
numSamplePoints <- 20000
minDistance <- 0

rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2012_2013_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2013_Medium_Buffer.shp"
plgDisturbanceFnShort <- "IDS2013_Medium_Buffer"
analysis07 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2012_2013_jDs_200_275_Landsat_NDVI_no7_TDDlm2.tif"
analysis08 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2013_2014_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2014_Medium_Buffer.shp"
plgDisturbanceFnShort <- "IDS2014_Medium_Buffer"
analysis09 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2013_2014_jDs_200_275_Landsat_NDVI_no7_TDDlm2.tif"
analysis10 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2014_2015_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2015_Medium_Buffer.shp"
plgDisturbanceFnShort <- "IDS2015_Medium_Buffer"
analysis11 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2014_2015_jDs_200_275_Landsat_NDVI_no7_TDDlm2.tif"
analysis12 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)


plgOuterFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\PA\\AOI_2014.shp"
plgOuterFnShort <- "AOI_2014"
numSamplePoints <- 100000
minDistance <- 0

rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160512\\PA_BL2010to2012_2013and2014_landsat_ZSCORE_Continuous.tif"

plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDSandBuffers\\SDEVW_DMG_2014_DECODE.shp"
plgDisturbanceFnShort <- "SDEVW_DMG_2014_DECODE"
analysis13 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDSandBuffers\\SDEVW_DMG_2014_DECODE_buf200m.shp"
plgDisturbanceFnShort <- "SDEVW_DMG_2014_DECODE_buf200m"
analysis14 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDSandBuffers\\SDEVW_DMG_2014_DECODE_buf400m.shp"
plgDisturbanceFnShort <- "SDEVW_DMG_2014_DECODE_buf400m"
analysis15 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDSandBuffers\\SDEVW_DMG_2014_DECODE_buf600m.shp"
plgDisturbanceFnShort <- "SDEVW_DMG_2014_DECODE_buf600m"
analysis16 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDSandBuffers\\SDEVW_DMG_2014_DECODE_buf800m.shp"
plgDisturbanceFnShort <- "SDEVW_DMG_2014_DECODE_buf800m"
analysis17 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDSandBuffers\\SDEVW_DMG_2014_DECODE_buf1000m.shp"
plgDisturbanceFnShort <- "SDEVW_DMG_2014_DECODE_buf1000m"
analysis18 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)


plgOuterFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\PA\\AOI_2015.shp"
plgOuterFnShort <- "AOI_2015"
numSamplePoints <- 100000
minDistance <- 0

rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160512\\PA_BL2011to2013_2014and2015_landsat_ZSCORE_Continuous.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDSandBuffers\\SDEVW_DMG_2015_DECODE.shp"
plgDisturbanceFnShort <- "SDEVW_DMG_2015_DECODE"
analysis19 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDSandBuffers\\SDEVW_DMG_2015_DECODE_buf200m.shp"
plgDisturbanceFnShort <- "SDEVW_DMG_2015_DECODE_buf200m"
analysis20 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDSandBuffers\\SDEVW_DMG_2015_DECODE_buf400m.shp"
plgDisturbanceFnShort <- "SDEVW_DMG_2015_DECODE_buf400m"
analysis21 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDSandBuffers\\SDEVW_DMG_2015_DECODE_buf600m.shp"
plgDisturbanceFnShort <- "SDEVW_DMG_2015_DECODE_buf600m"
analysis22 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDSandBuffers\\SDEVW_DMG_2015_DECODE_buf800m.shp"
plgDisturbanceFnShort <- "SDEVW_DMG_2015_DECODE_buf800m"
analysis23 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDSandBuffers\\SDEVW_DMG_2015_DECODE_buf1000m.shp"
plgDisturbanceFnShort <- "SDEVW_DMG_2015_DECODE_buf1000m"
analysis24 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)


plgOuterFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\CA_AOI.shp"
plgOuterFnShort <- "CA_AOI"
numSamplePoints <- 1000
minDistance <- 0
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\CA\\KNF_PM2015_landsat_TDD_ContinuousLM2NDVI_EL5.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\CA\\pandoraMoth2014_15_IDS_Alb.shp"
plgDisturbanceFnShort <- "pandoraMoth2014_15_IDS_Alb"
analysis25 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)









propInLocations <- 0.5

plgOuterFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\BHNF"
plgOuterFnShort <- "AOI"
numSamplePoints <- 100000
minDistance <- 0
propInLocations <- 0.5
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2014_2015_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2015.shp"
plgDisturbanceFnShort <- "IDS2015"

fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
minDistance <- 0
numSamplePoints <- 250000

rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2012_2013_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS_Mgmt2013_Dissolve.shp"
plgDisturbanceFnShort <- "IDS_Mgmt2013_Dissolve"
analysis01b <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2012_2013_jDs_200_275_Landsat_NDVI_no7_TDDlm2.tif"
analysis2013_no7 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2013_2014_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS_Mgmt2014_Dissolve.shp"
plgDisturbanceFnShort <- "IDS_Mgmt2014_Dissolve"
analysis2014 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2013_2014_jDs_200_275_Landsat_NDVI_no7_TDDlm2.tif"
analysis2014_no7 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2014_2015_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS_Mgmt2015_Dissolve.shp"
plgDisturbanceFnShort <- "IDS_Mgmt2015_Dissolve"
analysis2015 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2014_2015_jDs_200_275_LandsatNDVI_no7_TDDlm2.tif"
analysis2015_no7 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)


rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2012_2013_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2013.shp"
plgDisturbanceFnShort <- "IDS2013"
analysis2013r <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2012_2013_jDs_200_275_Landsat_NDVI_no7_TDDlm2.tif"
analysis2013r_no7 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2013_2014_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2014.shp"
plgDisturbanceFnShort <- "IDS2014"
analysis2014r <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2013_2014_jDs_200_275_Landsat_NDVI_no7_TDDlm2.tif"
analysis2014r_no7 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2014_2015_jDs_200_275_Landsat_NDVI_TDDlm2.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2015.shp"
plgDisturbanceFnShort <- "IDS2015"
analysis2015r <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160505\\BHNF_bYs_2003_2008_aYs_2014_2015_jDs_200_275_LandsatNDVI_no7_TDDlm2.tif"
analysis2015r_no7 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)


rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160325\\BHNF_bYs_2006_2009_aYs_2010_2011_jDs_200_275_ZScore.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2011.shp"
plgDisturbanceFnShort <- "IDS2011"
analysis2011old1 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160325\\BHNF_bYs_2006_2009_aYs_2010_2011_jDs_200_275_ZScore_Analysis2Zs_mean.tif"
analysis2011old2 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160325\\BHNF_bYs_2007_2010_aYs_2011_2012_jDs_200_275_ZScore.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2012.shp"
plgDisturbanceFnShort <- "IDS2012"
analysis2012old1 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160325\\BHNF_bYs_2007_2010_aYs_2011_2012_jDs_200_275_ZScore_Analysis2Zs_mean.tif"
analysis2012old2 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160325\\BHNF_bYs_2008_2011_aYs_2012_2013_jDs_200_275_ZScore.tif"
plgDisturbanceFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\IDS2013.shp"
plgDisturbanceFnShort <- "IDS2013"
analysis2013old1 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
rstXFn <- "D:\\Projects\\GIS_Projects\\20160310 ORS BHNF\\Models_20160325\\BHNF_bYs_2008_2011_aYs_2012_2013_jDs_200_275_ZScore_Analysis2Zs_mean.tif"
analysis2013old2 <- fullAnalysis(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)

write.csv(analysis2013$probabilities[[1]],"c:\\Users\\johnwithrow\\Desktop\\test.csv")
write.csv(analysis2013$probabilities[[1]],"c:\\Users\\johnwithrow\\Desktop\\Prob2013.csv")
write.csv(analysis2014$probabilities[[1]],"c:\\Users\\johnwithrow\\Desktop\\Prob2014.csv")
write.csv(analysis2015$probabilities[[1]],"c:\\Users\\johnwithrow\\Desktop\\Prob2015.csv")
write.csv(analysis2013r$probabilities[[1]],"c:\\Users\\johnwithrow\\Desktop\\Prob2013r.csv")
write.csv(analysis2014r$probabilities[[1]],"c:\\Users\\johnwithrow\\Desktop\\Prob2014r.csv")
write.csv(analysis2015r$probabilities[[1]],"c:\\Users\\johnwithrow\\Desktop\\Prob2015r.csv")
write.csv(analysis2013_no7$probabilities[[1]],"c:\\Users\\johnwithrow\\Desktop\\Prob2013_no7.csv")
write.csv(analysis2014_no7$probabilities[[1]],"c:\\Users\\johnwithrow\\Desktop\\Prob2014_no7.csv")
write.csv(analysis2015_no7$probabilities[[1]],"c:\\Users\\johnwithrow\\Desktop\\Prob2015_no7.csv")
write.csv(analysis2013r_no7$probabilities[[1]],"c:\\Users\\johnwithrow\\Desktop\\Prob2013r_no7.csv")
write.csv(analysis2014r_no7$probabilities[[1]],"c:\\Users\\johnwithrow\\Desktop\\Prob2014r_no7.csv")
write.csv(analysis2015r_no7$probabilities[[1]],"c:\\Users\\johnwithrow\\Desktop\\Prob2015r_no7.csv")

write.csv(analysis2013$roc[[1]],"c:\\Users\\johnwithrow\\Desktop\\Roc2013.csv")
write.csv(analysis2014$roc[[1]],"c:\\Users\\johnwithrow\\Desktop\\Roc2014.csv")
write.csv(analysis2015$roc[[1]],"c:\\Users\\johnwithrow\\Desktop\\Roc2015.csv")
write.csv(analysis2013r$roc[[1]],"c:\\Users\\johnwithrow\\Desktop\\Roc2013r.csv")
write.csv(analysis2014r$roc[[1]],"c:\\Users\\johnwithrow\\Desktop\\Roc2014r.csv")
write.csv(analysis2015r$roc[[1]],"c:\\Users\\johnwithrow\\Desktop\\Roc2015r.csv")
write.csv(analysis2013_no7$roc[[1]],"c:\\Users\\johnwithrow\\Desktop\\Roc2013_no7.csv")
write.csv(analysis2014_no7$roc[[1]],"c:\\Users\\johnwithrow\\Desktop\\Roc2014_no7.csv")
write.csv(analysis2015_no7$roc[[1]],"c:\\Users\\johnwithrow\\Desktop\\Roc2015_no7.csv")
write.csv(analysis2013r_no7$roc[[1]],"c:\\Users\\johnwithrow\\Desktop\\Roc2013r_no7.csv")
write.csv(analysis2014r_no7$roc[[1]],"c:\\Users\\johnwithrow\\Desktop\\Roc2014r_no7.csv")
write.csv(analysis2015r_no7$roc[[1]],"c:\\Users\\johnwithrow\\Desktop\\Roc2015r_no7.csv")

write.csv(analysis2013$samplesizes[[1]],"c:\\Users\\johnwithrow\\Desktop\\Samplesizes2013.csv")
write.csv(analysis2014$samplesizes[[1]],"c:\\Users\\johnwithrow\\Desktop\\Samplesizes2014.csv")
write.csv(analysis2015$samplesizes[[1]],"c:\\Users\\johnwithrow\\Desktop\\Samplesizes2015.csv")
write.csv(analysis2013r$samplesizes[[1]],"c:\\Users\\johnwithrow\\Desktop\\Samplesizes2013r.csv")
write.csv(analysis2014r$samplesizes[[1]],"c:\\Users\\johnwithrow\\Desktop\\Samplesizes2014r.csv")
write.csv(analysis2015r$samplesizes[[1]],"c:\\Users\\johnwithrow\\Desktop\\Samplesizes2015r.csv")
write.csv(analysis2013_no7$samplesizes[[1]],"c:\\Users\\johnwithrow\\Desktop\\Samplesizes2013_no7.csv")
write.csv(analysis2014_no7$samplesizes[[1]],"c:\\Users\\johnwithrow\\Desktop\\Samplesizes2014_no7.csv")
write.csv(analysis2015_no7$samplesizes[[1]],"c:\\Users\\johnwithrow\\Desktop\\Samplesizes2015_no7.csv")
write.csv(analysis2013r_no7$samplesizes[[1]],"c:\\Users\\johnwithrow\\Desktop\\Samplesizes2013r_no7.csv")
write.csv(analysis2014r_no7$samplesizes[[1]],"c:\\Users\\johnwithrow\\Desktop\\Samplesizes2014r_no7.csv")
write.csv(analysis2015r_no7$samplesizes[[1]],"c:\\Users\\johnwithrow\\Desktop\\Samplesizes2015r_no7.csv")

plotProbabilities(analysis2011old1,lblTitle="Validation Results - 2011")
plotROC(analysis2011old1,lblTitle="ROC Results - 2011")

plotProbabilities(analysis2011old2,lblTitle="Validation Results - 2011 (Mean)")
plotROC(analysis2011old2,lblTitle="ROC Results - 2011 (Mean)")

plotProbabilities(analysis2012old1,lblTitle="Validation Results - 2012")
plotROC(analysis2012old1,lblTitle="ROC Results - 2012")

plotProbabilities(analysis2012old2,lblTitle="Validation Results - 2012 (Mean)")
plotROC(analysis2012old2,lblTitle="ROC Results - 2012 (Mean)")

plotProbabilities(analysis2013old1,lblTitle="Validation Results - 2013")
plotROC(analysis2013old1,lblTitle="ROC Results - 2013")

plotProbabilities(analysis2013old2,lblTitle="Validation Results - 2013 (Mean)")
plotROC(analysis2013old2,lblTitle="ROC Results - 2013 (Mean)")