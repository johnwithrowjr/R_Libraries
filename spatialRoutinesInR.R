library(sp)
library(raster)
dist <- function(x1,y1,x2,y2) {sqrt((x1-x2)^2+(y1-y2)^2)}
pointInPolygon <- function(x,y,plgX)
{
	result <- FALSE
	for (i in 1:length(plgX@polygons))
	{
		for (j in 1:length(plgX@polygons[[i]]@Polygons))
		{
			crds <- plgX@polygons[[i]]@Polygons[[j]]@coords
			result <- result + point.in.polygon(x,y,crds[,1],crds[,2])
		}
	}
	#print (result)
	result
}
generateSamplePointsInPolygonWithMinDistance <- function(numPts,outerPolygon,minDistance=0,maxIter=10000,rstAssist=NULL)
{
	numIter <- 0
	samplepoints.x <- numeric(0)
	samplepoints.y <- numeric(0)
	while (length(samplepoints.x) < numPts & numIter < maxIter)
	{
		numIter <- numIter + 1
		x <- runif(1,outerPolygon@bbox[1,1],outerPolygon@bbox[1,2])
		y <- runif(1,outerPolygon@bbox[2,1],outerPolygon@bbox[2,2])
		if (pointInPolygon(x,y,outerPolygon)>0)
		{
			isolated = TRUE
			if (length(samplepoints.x) > 0)
			{
				for (i in 1:length(samplepoints.x))
				{
					if (dist(x,y,samplepoints.x[i],samplepoints.y[i]) < minDistance)
					{
						isolated = FALSE
					}
				}
			}
			if (isolated)
			{
				samplepoints.x <- c(samplepoints.x,x)
				samplepoints.y <- c(samplepoints.y,y)
				print (paste("numIter = ",numIter,", numPts = ",length(samplepoints.x)," x=",x," y=",y,sep=""))
			}
		}
	}
	result <- cbind(samplepoints.x,samplepoints.y)
	colnames(result) <- c("x","y")
	result
}

generateSamplePointsWithORSBalanceAndMinDistance <- function(numPtsTotal,rstORS,plgOuter,breaks,minDistance=0,seed=NULL)
{
  if (!is.null(seed)) {set.seed(seed)}
  dfORSPts <- cbind(coordinates(rstORS),rstORS@data[[1]])
  colnames(dfORSPts) <- c("x","y","z")
  dfORSPts <- dfORSPts[!is.na(dfORSPts[,3]),]
  m <- length(breaks) - 1
  result <- matrix(nrow=0,ncol=4)
  numPts <- ceiling(numPtsTotal/m)
  for (i in 1:m)
  {
    dfORSPtsTemp <- dfORSPts[dfORSpts[,3]>=breaks[i] & dfORSPts[,3] < breaks[i+1],]
    totalrows <- dim(dfORSPtsTemp)[1]
    dfORSPtsTemp <- dfORSPtsTemp[sample(totalrows,totalrows),]
    j <- 0
    resultTemp <- matrix(nrow=0,ncol=4)
    while (j < numPts & k < totalrows)
    {
      isolated <- TRUE
      if (j > 0)
      {
        for (p in 1:j)
        {
          if (dist(resultTemp[p,1],resultTemp[p,2],dfORSPtsTemp[k,1],dfORSPtsTemp[k,2])<minDistance) { isolated <- FALSE }
        }
      }
      if (isolated) {resultTemp <- rbind(resultTemp,cbind(dfORSPtsTemp[k,],m))}
      k <- k + 1
    }
    result <- rbind(result,resultTemp)
  }
  result
}

generateSamplePointsInIsolatedLocationsWithMinDistance <- function(numPtsTotal,rstLocations,plgOuter,rstGuide,propInLocations=0.5,minDistance=0,maxIter=10000)
{
	#rstTemp <- as(calc(stack(rstLocations,as(rstGuide,"RasterLayer")),fun=function(x){x[1]*x[2]}),"SpatialGridDataFrame")
	#rstGuide <- as(rstGuide,"SpatialGridDataFrame")
	#dfGuide <- rstGuide@data[[1]][!is.na(rstGuide@data[[1]])]
	print("Converting to a List of Points...")
	dfLocationsPts <- cbind(coordinates(rstLocations),!is.na(rstLocations@data[[1]]))
	#print(dfLocations)
	dfLocationsPtsInside <- dfLocationsPts[dfLocationsPts[,3]==1,1:2]
	dfLocationsPtsOutside <- dfLocationsPts[dfLocationsPts[,3]==0,1:2]
	print("Randomizing...")
	dfLocationsPtsInside <- dfLocationsPtsInside[sample(1:dim(dfLocationsPtsInside)[1]),]
	dfLocationsPtsOutside <- dfLocationsPtsOutside[sample(1:dim(dfLocationsPtsOutside)[1]),]
	#pI <- (dim(dfLocationsPtsInside)[1])/(dim(dfLocationsPtsInside)[1]+dim(dfLocationsPtsOutside)[1])
	samplepoints.x <- numeric(0)
	samplepoints.y <- numeric(0)
	samplepoints.z <- numeric(0)
	numIter <- 0
	numPts <- floor(numPtsTotal/2.0)
	m <- sample(1:dim(dfLocationsPtsInside)[1])
	while (length(samplepoints.x) < numPts & numIter < maxIter)
	{
		numIter <- numIter + 1
		x <- dfLocationsPtsInside[m[numIter],1]
		y <- dfLocationsPtsInside[m[numIter],2]
		if (pointInPolygon(x,y,plgOuter))
		{
			isolated = TRUE
			if (length(samplepoints.x) > 0)
			{
				if (minDistance > 0)
				{
					for (i in 1:length(samplepoints.x))
					{
						if (dist(x,y,samplepoints.x[i],samplepoints.y[i]) < minDistance)
						{
							isolated = FALSE
						}
					}
				}
			}
			if (isolated)
			{
				samplepoints.x <- c(samplepoints.x,x)
				samplepoints.y <- c(samplepoints.y,y)
				samplepoints.z <- c(samplepoints.z,1)
				if (isDivisibleBy(length(samplepoints.x),100)) {print (paste("numIter = ",numIter,", numPts = ",length(samplepoints.x)," x=",x," y=",y,sep=""))}
			}
		}
	}
	numIter <- 0
	m <- sample(1:dim(dfLocationsPtsOutside)[1])
	numptsTotal <- 2 * length(samplepoints.x)
	while (length(samplepoints.x) < numPtsTotal & numIter < maxIter)
	{
		numIter <- numIter + 1
		x <- dfLocationsPtsOutside[m[numIter],1]
		y <- dfLocationsPtsOutside[m[numIter],2]
		if (pointInPolygon(x,y,plgOuter))
		{
			isolated = TRUE
			if (length(samplepoints.x) > 0)
			{
				if (minDistance > 0)
				{
					for (i in 1:length(samplepoints.x))
					{
						if (dist(x,y,samplepoints.x[i],samplepoints.y[i]) < minDistance)
						{
							isolated = FALSE
						}
					}
				}
			}
			if (isolated)
			{
				samplepoints.x <- c(samplepoints.x,x)
				samplepoints.y <- c(samplepoints.y,y)
				samplepoints.z <- c(samplepoints.z,0)
				if (isDivisibleBy(length(samplepoints.x),100)) {print (paste("numIter = ",numIter,", numPts = ",length(samplepoints.x)," x=",x," y=",y,sep=""))}
			}
		}
	}
	result <- cbind(samplepoints.x,samplepoints.y,samplepoints.z)
	colnames(result) <- c("x","y","z")
	#pI <- (dim(dfLocationsPtsInside[,])[1])/(dim(dfLocationsPtsInside)[1]+dim(dfLocationsPtsOutside)[1])
	result
}

sampleVariogramOfRaster <- function(rstX, maxDistance, numBreaks, maxIter = 10000)
{
	if (!(class(rstX)=="SpatialGridDataFrame")) { stop("Object of type SpatialGridDataFrame needed") }
	rstXcol <- cbind(coordinates(rstX),rstX@data[[1]])
	colnames(rstXcol) <- c("x","y","z")
	rstXcol <- subset(rstXcol,!is.na(rstXcol[,3]))
	n <- dim(rstXcol)[1]
	breaks <- cbind((0:(numBreaks-1))*maxDistance/numBreaks,(1:numBreaks)*maxDistance/numBreaks)
	print(breaks)
	pts <- rep(0,numBreaks)
	smpPts <- NULL
	print("Getting samples...")
	s1 <- sample(1:n,maxIter)
	s2 <- sample(1:n,maxIter)
	print("Placing into bins...")
	k <- 0
	maxIterDec <- as.integer(maxIter/10)
	for (i in 1:maxIter)
	{
		if (i > (maxIterDec*(k+1))) {
			k <- k + 1
			print (paste(10*k,"% Done...",sep=""))
		}
		z <- dist(rstXcol[s1[i],1],rstXcol[s1[i],2],rstXcol[s2[i],1],rstXcol[s2[i],2])
		#print(z)
		if (z<maxDistance)
		{
			smpPts <- rbind(smpPts,c(as.integer(z*numBreaks/maxDistance)+1,rstXcol[s1[i],3],rstXcol[s2[i],3]))
		}
	}
	print("Calculating correlations...")
	smpPtsCorr <- NULL
	for (i in 1:numBreaks)
	{
		temp <- subset(smpPts,smpPts[,1]==i)
		#print(temp)
		temp.cor <- NULL
		if (dim(temp)[1] < 3)
		{
			temp.cor <- list(estimate=NA, conf.int=c(NA,NA))
		} else {
			temp.cor <- cor.test(temp[,2],temp[,3])
		}
		smpPtsCorr <- rbind(smpPtsCorr,c(temp.cor$estimate,temp.cor$conf.int[1],temp.cor$conf.int[2]))
		#print(smpPtsCorr)
	}
	result <- cbind(rowMeans(breaks),breaks,smpPtsCorr)
	colnames(result) <- c("xmid","xlow","xhigh","corr","corr_low_95","corr_high_95")
	result <- as.data.frame(result)
	plot(result$xmid,result$corr,ylim=c(-1,1),xlab="Distance",ylab="Pearson Correlation Coefficient")
	abline(0,0)
	lines(result$xmid,result$corr,lwd=2)
	lines(result$xmid,result$corr_low_95,lwd=1,lty=3)
	lines(result$xmid,result$corr_high_95,lwd=1,lty=3)
	title("Sample Correllogram")
	result
}

simpleRaster <- function()
{
	df = data.frame(z = c(1:6,NA,8,9),
    xc = c(1,1,1,2,2,2,3,3,3),
    yc = c(rep(c(0, 1.5, 3),3)))
	coordinates(df) = ~xc+yc
	gridded(df) = TRUE
	df = as(df, "SpatialGridDataFrame") # to full grid
	df
}

valuesOfRasterWithinPolygon <- function(rstX,plgX,inside=TRUE)
{
	rstXcol <- cbind(coordinates(rstX),rstX@data[[1]])
	colnames(rstXcol) <- c("x","y","z")
	x <- numeric(0)
	z <- numeric(0)
	for (i in 1:dim(rstXcol)[1])
	{
		if (inside == as.boolean(pointInPolygon(rstXcol[i,1],rstXcol[i,2],plgX)>0)) 
		{
			x <- c(x,i)
			z <- c(z,rstXcol[i,3])
		}
	}
	cbind(x,z)
}

valuesOfRasterAtPoints <- function(rstX,smpPts)
{
	matX <- as.matrix(rstX@data[[1]])
	dim(matX) <- rstX@grid@cells.dim
	smpPts[,1] <- floor((smpPts[,1] - rstX@bbox[1,1])/rstX@grid@cellsize[1]) + 1
	smpPts[,2] <- rstX@grid@cells.dim[2] - floor((smpPts[,2] - rstX@bbox[2,1])/rstX@grid@cellsize[2])
	n <- dim(smpPts)[1]
	result <- numeric(n)
	for (i in 1:n)
	{
		if (smpPts[i,1]>0 & smpPts[i,2]>0 & smpPts[i,1]<=dim(matX)[1] & smpPts[i,2]<=dim(matX)[2])
		{
			result[i] <- matX[smpPts[i,1],smpPts[i,2]]
		} else {
			result[i] <- NA
			print (paste(str(smpPts[i,])," is outside ",str(rstX@bbox)))
		}
	}
	result
}

confirmSameGrid <- function(rstX,rstY,sampleProp=0.1)
{
	# This function takes two parameters, both of type SpatialGridDataFrame {sp}.  It checks that both grids have the same
	# bounding boxes, projection strings, first and last grid coordinates, as well as the same corresponding grid coordinates
	# from a random sample of 10% of the grid.  TRUE is returned if all those tests are true.  Otherwise it returns FALSE.
	if (!(class(rstX)=="SpatialGridDataFrame")) { stop("Object of type SpatialGridDataFrame needed for first parameter.") }
	if (!(class(rstY)=="SpatialGridDataFrame")) { stop("Object of type SpatialGridDataFrame needed for second parameter.") }
	if (!all(rstX@bbox==rstY@bbox)) { return(-1) }
	prjX <- rstX@proj4string@projargs
	prjY <- rstY@proj4string@projargs
	if ((is.na(prjX) | is.na(prjY)) & !(is.na(prjX)&is.na(prjY))) { return(-2) }
	if (!(is.na(prjX)) & !(is.na(prjY)) & !(rstX@proj4string@projargs==rstY@proj4string@projargs)) { return(-3) }
	if (!(coordinates(rstX)[1,1]==coordinates(rstY)[1,1])) { return(-4) }
	if (!(coordinates(rstX)[1,2]==coordinates(rstY)[1,2])) { return(-5) }
	if (!(coordinates(rstX)[dim(coordinates(rstX))[1],1]==coordinates(rstY)[dim(coordinates(rstY))[1],1])) { return(-6) }
	if (!(coordinates(rstX)[dim(coordinates(rstX))[1],2]==coordinates(rstY)[dim(coordinates(rstY))[1],2])) { return(-7) }
	m <- ceiling(dim(coordinates(rstX))[1]*sampleProp)
	smpPts <- sample(dim(coordinates(rstX))[1],m)
	for (i in 1:m)
	{
		if (!(coordinates(rstX)[i,1]==coordinates(rstY)[i,1])) { return(-8) }
		if (!(coordinates(rstX)[i,2]==coordinates(rstY)[i,2])) { return(-9) }
	}
	1
}

timeSequenceGrids <- function(gridSeries,dateSeries)
{
	if (length(gridSeries)!=length(dateSeries)) { stop("GridSeries and dateSeries must be of the same length.") }
	n <- length(gridSeries)
	if (n > 1)
	{
		for (i in 1:(n-1))
		{
			if (!confirmSameGrid(gridSeries[[i]],gridSeries[[i+1]])) { stop(paste("Grids ",i," and ",(i+1)," are of different projections or cell sizes.",sep="")) }
		}
	}
	result <- coordinates(gridSeries[[1]])
	for (i in 1:n)
	{
		result <- cbind(result,gridSeries[[i]]@data)
	}
	colnames(result) <- c("x","y",dateSeries)
	result
}

addPolygon <- function(x,plgX,plgName)
{
	n <- dim(x)[1]
	result <- numeric(n)
	for (i in 1:n)
	{
		result[i] <- pointInPolygon(x[i,1],x[i,2],plgX)
	}
	result <- cbind(x,result)
	colnames(result)[dim(result)[2]] <- plgName
	result
}

analysisFirstOrder <- function(x,rstX,rstLocations,params,polys)
{
	x <- x[complete.cases(x),]
	header <- list(length(params))
	pI <- list(length(params))
	pX <- list(length(params))
	analyses <- list(length(polys))
	samplesizes <- list(length(polys))
	names(header) <- colnames(x[,params])
	names(analyses) <- colnames(x[,polys])
	for (i in 1:length(header))
	{
		#header[[i]] <- histogramBreaksSturges(rstX@data[!is.na(rstX@data)],length(x[,params[i]]))
		header[[i]] <- histogramBreaksSturges(x[,params[i]],length(x[,params[i]]))
		pX[[i]] <- hist(rstX@data[!is.na(rstX@data) & rstX@data>header[[i]]$xMin & rstX@data<header[[i]]$xMax],breaks=(c(header[[i]]$bins[1,1],header[[i]]$bins[,2])),plot=FALSE)$density
		isWithinRange <- !is.na(rstX@data) & rstX@data>header[[i]]$xMin & rstX@data<header[[i]]$xMax
		pI[[i]] <- length(rstLocations@data[!is.na(rstLocations@data) & isWithinRange]) / length(rstLocations@data[isWithinRange])
	}
	for (j in 1:length(polys))
	{
		buildTable <- matrix(nrow=0,ncol=16)
		colnames(buildTable) <- c("XMid","Successes","Failures","PmaxL(I|X)","Pmean(I|X)","P95-(I|X)","P95+(I|X)","P99-(I|X)","P99+(I|X)","P(X)","PmaxL(X|I)","Pmean(X|I)","P95-(X|I)","P95+(X|I)","P99-(X|I)","P99+(X|I)")
		for (i in 1:length(params))
		{
			for (k in 1:header[[i]]$nBins)
			{
				ss <- x[(x[,params[i]]>=header[[i]]$bins[k,1] & x[,params[i]]<header[[i]]$bins[k,2] & x[,polys[j]]==1),]
				if (is.matrix(ss))
				{
					s <- dim(ss)[1]
				} else {
					s <- 1
				}
				ff <- x[(x[,params[i]]>=header[[i]]$bins[k,1] & x[,params[i]]<header[[i]]$bins[k,2] & x[,polys[j]]==0),]
				if (is.matrix(ff))
				{
					f <- dim(ff)[1]
				} else {
					f <- 1
				}
				#print(header[[i]]$bins[k,1])
				buildTable <- rbind(buildTable,c(header[[i]]$bins[k,3],s,f,betaSuccessMode(s,f),betaSuccessMean(s,f),betaSuccessQuantiles(c(0.025,0.975,0.005,0.995),s,f),pX[[i]][k],betaSuccessMode(s,f)*pX[[i]][k]/pI[[i]],betaSuccessMean(s,f)*pX[[i]][k]/pI[[i]],betaSuccessQuantiles(c(0.025,0.975,0.005,0.995),s,f)*pX[[i]][k]/pI[[i]]))
			}
		}
		analyses[[j]] <- buildTable
		plot(buildTable[,1],buildTable[,5],ylim=c(0,1),xlab="ORS Value",ylab="Proportion of Samples in IDS Polygons")
		lines(buildTable[,1],buildTable[,5])
		lines(buildTable[,1],buildTable[,6],lty=2)
		lines(buildTable[,1],buildTable[,7],lty=2)
		lines(buildTable[,1],buildTable[,8],lty=3)
		lines(buildTable[,1],buildTable[,9],lty=3)
		title("ORS / IDS Comparison Results")
		n <- as.matrix(rep(0,dim(buildTable)[1]^2))
		dim(n) <- c(dim(buildTable)[1],dim(buildTable)[1])
		for (p in 1:(dim(n)[1]-1))
		{
			for (q in (p+1):(dim(n)[1]))
			{
				n[p,q] <- observationsRequiredForBinomialDifferentiation(buildTable[p,5],buildTable[q,5],0.05)
			}
		}
		samplesizes[[j]] <- n
	}
	list(x=x,pI=pI,header=header,pX=pX,analyses=analyses,samplesizes=samplesizes)
}

fullAnalysis <- function(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,propInLocations=0.5,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
{
	print("Reading in Outer Polygons...")
	plgOuter <- readOGR(plgOuterFn,plgOuterFnShort)
	print("Reading in Disturbance Polygons...")
	plgDisturbance <- readOGR(plgDisturbanceFn,plgDisturbanceFnShort)
	print("Reading in Raster Layer...")
	rstX <- readGDAL(rstXFn)
	rstX <- as(rstX,"SpatialGridDataFrame")
	print("Creating Rasterized Polygons...")
	rstDisturbance <- raster(rstX,layer=1,values=FALSE)
	rstDisturbance <- rasterize(plgDisturbance,rstDisturbance,colnames(plgDisturbance@data)[1])
	rstDisturbance <- as(rstDisturbance,"SpatialGridDataFrame")
	samplePts <- generateSamplePointsInIsolatedLocationsWithMinDistance(numSamplePoints,rstDisturbance,plgOuter,rstX,propInLocations,minDistance,1000000)
	samplePts <- cbind(samplePts,valuesOfRasterAtPoints(rstX,samplePts))
	colnames(samplePts)[4] <- "rstX"
 	analysisFirstOrder(samplePts,rstX,rstDisturbance,4,3)
}
