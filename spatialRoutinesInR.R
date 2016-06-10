library(sp)
library(rgdal)
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
			if (!(plgX@polygons[[i]]@Polygons[[j]]@hole)) { result <- result + point.in.polygon(x,y,crds[,1],crds[,2]) }
		}
	}
	#print (result)
	result
}
generateSamplePointsInPolygonWithMinDistance <- function(numPts,outerPolygon,minDistance=0,maxIter=10000,seed=NULL)
{
	if (!is.null(seed)) {set.seed(seed)}
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

generateSamplePointsInStratifiedLocationsWithMinDistance <- function(numPtsTotal,rstLocations,plgOuter,minDistance=0,maxIter=10000)
{
	#rstTemp <- as(calc(stack(rstLocations,as(rstGuide,"RasterLayer")),fun=function(x){x[1]*x[2]}),"SpatialGridDataFrame")
	#rstGuide <- as(rstGuide,"SpatialGridDataFrame")
	#dfGuide <- rstGuide@data[[1]][!is.na(rstGuide@data[[1]])]
	print("Converting to a List of Points...")
	dfLocationsPts <- cbind(coordinates(rstLocations),!is.na(rstLocations@data[[1]]))
	print(dfLocationsPts)
	dfLocationsPtsInside <- dfLocationsPts[dfLocationsPts[,3]==1,1:2]
	dfLocationsPtsOutside <- dfLocationsPts[dfLocationsPts[,3]==0,1:2]
	print(dfLocationsPtsInside)
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
		#print(paste(length(samplepoints.x),numPts,sep=" : "))
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

exportSamplePoints <- function(workspace,fn,analObj,projObj)
{
	pts <- analObj$x
	dimnames(pts)[[1]] <- seq(1,dim(pts)[1])
	pts <- as.data.frame(pts)
	temp <- SpatialPointsDataFrame(pts[,1:2],data=pts,proj4string=CRS(proj4string(projObj)))
	#coordinates(temp) <- pts[,1:2]
	#temp <- SpatialPointsDataFrame(pts[,c(1,2)],pts[,3])
	#print("Got this far")
	#proj4string(temp) <- proj4string(projObj)
	#print("got this far")
	writeOGR(obj=temp,dsn=paste(workspace,fn,".shp",sep=""),layer=fn,driver="ESRI Shapefile")
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

analysisFirstOrderOld <- function(x,rstX,rstLocations,params,polys)
{
	x <- x[complete.cases(x),]
	header <- list(length(params))
	N <- list(length(params))
	nI <- list(length(params))
	nNotI <- list(length(params))
	pI <- list(length(params))
	pNotI <- list(length(params))
	probabilities <- list(length(polys))
	roc <- list(length(polys))
	samplesizes <- list(length(polys))
	optMax <- list(length(polys))
	optPt <- list(length(polys))
	optMtx <- list(length(polys))
	kappa <- list(length(polys))
	names(header) <- colnames(x[,params])
	names(probabilities) <- colnames(x[,polys])
	names(roc) <- colnames(x[,polys])
	for (i in 1:length(header))
	{
		header[[i]] <- histogramBreaksSturges(x[,params[i]],length(x[,params[i]]))
		isWithinRange <- !is.na(rstX@data) & rstX@data>header[[i]]$xMin & rstX@data<header[[i]]$xMax
		nI[[i]] <- length(rstLocations@data[!is.na(rstLocations@data) & isWithinRange])
		N[[i]] <- length(rstLocations@data[isWithinRange])
		nNotI[[i]] <- N[[i]] - nI[[i]]
		pI[[i]] <- nI[[i]] / N[[i]]
		pNotI[[i]] <- 1 - pI[[i]]
	}
	for (j in 1:length(polys))
	{
		for (i in 1:length(params))
		{
			buildTable <- matrix(nrow=0,ncol=12)
			colnames(buildTable) <- c("Xmin","XMid","Xmax","nXI","nnotXI","Pmean(X|I)","p(I)","nXnotI","nnotXnotI","Pmean(X|-I)","P(-I)","P(I|X)")
			for (k in 1:header[[i]]$nBins)
			{
				nXI       <- numRows(x[(x[,params[i]]>=header[[i]]$bins[k,1] & x[,params[i]]<header[[i]]$bins[k,2] & x[,polys[j]]==1),])
				nnotXI    <- numRows(x[((x[,params[i]]<header[[i]]$bins[k,1] | x[,params[i]]>=header[[i]]$bins[k,2]) & x[,polys[j]]==1),])
				nXnotI    <- numRows(x[(x[,params[i]]>=header[[i]]$bins[k,1] & x[,params[i]]<header[[i]]$bins[k,2] & x[,polys[j]]==0),])
				nnotXnotI <- numRows(x[((x[,params[i]]<header[[i]]$bins[k,1] | x[,params[i]]>=header[[i]]$bins[k,2]) & x[,polys[j]]==0),])
				pXgivenI <- betaSuccessMean(nXI,nnotXI)
				#pI <- betaSuccessMean(nI[[i]],nNotI[[i]])
				pXgivennotI <- betaSuccessMean(nXnotI,nnotXnotI)
				#pnotI <- 1 - pI
				pIgivenX <- pXgivenI*pI[[i]]/(pXgivenI*pI[[i]]+pXgivennotI*pNotI[[i]])
				buildTable <- rbind(buildTable,c(header[[i]]$bins[k,1],header[[i]]$bins[k,3],header[[i]]$bins[k,2],nXI,nnotXI,pXgivenI,pI[[i]],nXnotI,nnotXnotI,pXgivennotI,pNotI[[i]],pIgivenX))
			}

			np <- 1000
			dx <- (header[[i]]$xMax - header[[i]]$xMin) / np
			xMin <- header[[i]]$xMin
			buildTable2 <- matrix(nrow=0,ncol=3)
			colnames(buildTable2) <- c("x","tp","fp")
			optMax <- 0
			optPt <- 0
			for (p in 1:np)
			{
				Xstar <- xMin + dx*p
				print(Xstar)
				nXI       <- numRows(x[x[,params[i]]<Xstar & x[,polys[j]]==1,])
				nnotXI    <- numRows(x[x[,params[i]]>=Xstar & x[,polys[j]]==1,])
				nXnotI    <- numRows(x[x[,params[i]]<Xstar & x[,polys[j]]==0,])
				nnotXnotI <- numRows(x[x[,params[i]]>=Xstar & x[,polys[j]]==0,])
				nXI_nI <- nXI/numRows(x[x[,polys[j]]==1,])
				nXnotI_nnotI <- nXnotI/numRows(x[x[,polys[j]]==0,])
				ss <- nXI_nI + (1 - nXnotI_nnotI)
				if (ss > optMax) { optMax <- ss ; optPt <- Xstar ; fnXI <- nXI ; fnnotXI <- nnotXI ; fnXnotI <- nXnotI ; fnnotXnotI <- nnotXnotI }
				buildTable2 <- rbind(buildTable2,c(Xstar,nXI_nI,nXnotI_nnotI))
			}
			buildTable3 <- unique(buildTable2[,c(3,2)])
			print(buildTable3)
			AUC <- 0
			np <- dim(buildTable3)[1]
			if (buildTable3[1,1] > 0) { AUC <- AUC + buildTable3[1,2]*(buildTable3[1,1])/2 }
			for (p in 2:(np-1))
			{
				AUC <- AUC + buildTable3[p,2]*(buildTable3[p+1,1]-buildTable3[p-1,1])/2
			}
			if (buildTable3[np,1] > 0) { AUC <- AUC + buildTable3[np,2]*(1-buildTable3[np,1])/2 }
		}
		probabilities[[i]] <- buildTable
		roc[[i]] <- buildTable2
		optMax[[i]] <- optMax
		optPt[[i]] <- optPt
		optMtx[[i]] <- as.matrix(c(fnXI/N[[i]],fnnotXI/N[[i]],fnXnotI/N[[i]],fnnotXnotI/N[[i]]),nrow=2)
		pX <- fnXI/(fnXI+fnnotXI)*pI[[i]] + fnXnotI/(fnXnotI+fnnotXnotI)*(1-pI[[i]])
		pO <- (fnXI/N[[i]])*pI[[i]] + (fnnotXnotI/N[[i]])*(1-pI[[i]])
		pE <- pX*pI[[i]] + (1-pX)*(1-pI[[i]])
		kappa[[i]] <- (pO-pE)/(1-pE)
		samplesizes <- as.matrix(rep(0,dim(buildTable)[1]^2))
		dim(samplesizes) <- c(dim(buildTable)[1],dim(buildTable)[1])
		for (p in 1:(dim(samplesizes)[1]-1))
		{
			print(paste(p,dim(samplesizes)[1],sep=" "))
			for (q in (p+1):(dim(samplesizes)[1]))
			{
				samplesizes[p,q] <- observationsRequiredForBinomialDifferentiation(buildTable[p,12],buildTable[q,12],0.05)
			}
		}
	}
	list(x=x,pI=pI,header=header,probabilities=probabilities,roc=roc,AUC=AUC,optMax=optMax,optPt=optPt,optMtx=optMtx,kappa=kappa,samplesizes=samplesizes,buildTable2=buildTable2)
}

analysisFirstOrder <- function(x,rstX,rstLocations,param,poly)
{
	x <- x[complete.cases(x),]
	x <- as.matrix(x)
	x[,param] <- as.numeric(x[,param])
	x[,poly] <- as.numeric(x[,poly])
	header <- histogramBreaksSturges(x[,param],length(x[,param]))
	isWithinRange <- !is.na(rstX@data) & rstX@data>header$xMin & rstX@data<header$xMax
	nI <- length(rstLocations@data[!is.na(rstLocations@data) & isWithinRange])
	N <- length(rstLocations@data[isWithinRange])
	nNotI <- N - nI
	pI <- nI / N
	pNotI <- 1 - pI
	buildTable <- matrix(nrow=0,ncol=12)
	colnames(buildTable) <- c("Xmin","XMid","Xmax","nXI","nnotXI","Pmean(X|I)","p(I)","nXnotI","nnotXnotI","Pmean(X|-I)","P(-I)","P(I|X)")
	for (k in 1:header$nBins)
	{
		nXI       <- numRows(x[(x[,param]>=header$bins[k,1] & x[,param]<header$bins[k,2] & x[,poly]==1),])
		nnotXI    <- numRows(x[((x[,param]<header$bins[k,1] | x[,param]>=header$bins[k,2]) & x[,poly]==1),])
		nXnotI    <- numRows(x[(x[,param]>=header$bins[k,1] & x[,param]<header$bins[k,2] & x[,poly]==0),])
		nnotXnotI <- numRows(x[((x[,param]<header$bins[k,1] | x[,param]>=header$bins[k,2]) & x[,poly]==0),])
		pXgivenI <- betaSuccessMean(nXI,nnotXI)
		pXgivennotI <- betaSuccessMean(nXnotI,nnotXnotI)
		pIgivenX <- pXgivenI*pI/(pXgivenI*pI+pXgivennotI*pNotI)
		buildTable <- rbind(buildTable,c(header$bins[k,1],header$bins[k,3],header$bins[k,2],nXI,nnotXI,pXgivenI,pI,nXnotI,nnotXnotI,pXgivennotI,pNotI,pIgivenX))
	}
	np <- 1000
	dx <- (header$xMax - header$xMin) / np
	xMin <- header$xMin
	buildTable2 <- matrix(nrow=0,ncol=3)
	colnames(buildTable2) <- c("x","tp","fp")
	optMax <- 0
	optPt <- 0
	for (p in 1:np)
	{
		Xstar <- xMin + dx*p
		print(Xstar)
		nXI       <- numRows(x[x[,param]<Xstar & x[,poly]==1,])
		nnotXI    <- numRows(x[x[,param]>=Xstar & x[,poly]==1,])
		nXnotI    <- numRows(x[x[,param]<Xstar & x[,poly]==0,])
		nnotXnotI <- numRows(x[x[,param]>=Xstar & x[,poly]==0,])
		nXI_nI <- nXI/numRows(x[x[,poly]==1,])
		nXnotI_nnotI <- nXnotI/numRows(x[x[,poly]==0,])
		ss <- nXI_nI + (1 - nXnotI_nnotI)
		if (ss > optMax) { optMax <- ss ; optPt <- Xstar ; fnXI <- nXI ; fnnotXI <- nnotXI ; fnXnotI <- nXnotI ; fnnotXnotI <- nnotXnotI }
		buildTable2 <- rbind(buildTable2,c(Xstar,nXI_nI,nXnotI_nnotI))
		print(paste(Xstar,nXI_nI,nXnotI_nnotI,sep=" : "))
	}
	buildTable3 <- unique(buildTable2[,c(3,2)])
	print(buildTable3)
	AUC <- 0
	np <- dim(buildTable3)[1]
	if (buildTable3[1,1] > 0) { AUC <- AUC + buildTable3[1,2]*(buildTable3[1,1])/2 }
	for (p in 2:(np-1)) { AUC <- AUC + buildTable3[p,2]*(buildTable3[p+1,1]-buildTable3[p-1,1])/2 }
	if (buildTable3[np,1] > 0) { AUC <- AUC + buildTable3[np,2]*(1-buildTable3[np,1])/2 }
	probabilities <- buildTable
	roc <- buildTable2
	optMax <- optMax
	optPt <- optPt
	optMtx <- as.matrix(c(fnXI/N,fnnotXI/N,fnXnotI/N,fnnotXnotI/N),nrow=2)
	pX <- fnXI/(fnXI+fnnotXI)*pI + fnXnotI/(fnXnotI+fnnotXnotI)*(1-pI)
	pO <- (fnXI/N)*pI + (fnnotXnotI/N)*(1-pI)
	pE <- pX*pI + (1-pX)*(1-pI)
	kappa <- (pO-pE)/(1-pE)
	samplesizes <- as.matrix(rep(0,dim(buildTable)[1]^2))
	dim(samplesizes) <- c(dim(buildTable)[1],dim(buildTable)[1])
	for (p in 1:(dim(samplesizes)[1]-1))
	{
		print(paste(p,dim(samplesizes)[1],sep=" "))
		for (q in (p+1):(dim(samplesizes)[1]))
		{
			samplesizes[p,q] <- observationsRequiredForBinomialDifferentiation(buildTable[p,12],buildTable[q,12],0.05)
		}
	}
	list(x=x,pI=pI,header=header,probabilities=probabilities,roc=roc,AUC=AUC,optMax=optMax,optPt=optPt,optMtx=optMtx,kappa=kappa,samplesizes=samplesizes,buildTable2=buildTable2)
}

plotROC <- function(analysis,lblTitle="ROC Curve")
{
	plot(analysis$roc[,3],analysis$roc[,2],xlab="False Positives",ylab="True Positives",main=lblTitle)
	abline(0,1)
	text(.8,.2,paste("AUC=",analysis$AUC,sep=""))
}

plotROC2 <- function(analysis,lblTitle="ROC Curve")
{
	plot(analysis$roc[[1]][,3],analysis$roc[[1]][,2],xlab="False Positives",ylab="True Positives",main=lblTitle)
	abline(0,1)
	text(.8,.2,paste("AUC=",analysis$AUC[[1]],sep=""))
}

plotProbabilities <- function(analysis,lblX="ORS Value",lblY="Proportion of Samples in IDS Polygons",lblTitle="ORS / IDS Comparison Results")
{
	plot(analysis$probabilities[[1]][,1],analysis$probabilities[[1]][,12],ylim=c(0,1),xlab=lblX,ylab=lblY)
	lines(analysis$probabilities[[1]][,1],analysis$probabilities[[1]][,12])
	title(lblTitle)
}

plotProbabilities2 <- function(analysis,lblX="ORS Value",lblY="Proportion of Samples in IDS Polygons",lblTitle="ORS / IDS Comparison Results")
{
	plot(analysis$probabilities[,1],analysis$probabilities[,12],ylim=c(0,1),xlab=lblX,ylab=lblY)
	lines(analysis$probabilities[,1],analysis$probabilities[,12])
	title(lblTitle)
}

createMask <- function(rstX)
{
	rstY <- raster(rstX,layer=1,values=TRUE)
	rstY <- as(rstY,"SpatialGridDataFrame")
	for (i in 1:length(rstX@data))
	{
		if (!is.na(rstX@data[i,1]))
		{
			if (rstX@data[i,1] > 0)
			{
				rstY@data[i,1] <- 1
			} else {
				rstY@data[i,1] <- 0
			}
		} else {
			rstY@data[i,1] <- 0
		}
	}
	rstY
}

fullAnalysis <- function(plgOuterFn,plgOuterFnShort,numSamplePoints,minDistance,rstXFn,plgDisturbanceFn,plgDisturbanceFnShort)
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
	rstOuter <- raster(rstX,layer=1,values=FALSE)
	rstDisturbance <- rasterize(plgDisturbance,rstDisturbance,colnames(plgDisturbance@data)[1])
	rstDisturbance <- as(rstDisturbance,"SpatialGridDataFrame")
	rstOuter <- rasterize(plgOuter,rstOuter,colnames(plgOuter@data)[1])
	rstOuter <- as(rstOuter,"SpatialGridDataFrame")
	rstOuter <- createMask(rstOuter)
	#rstDisturbance@data <- as.data.frame(as.vector(rstOuter@data[,1]) * as.vector(rstDisturbance@data[,1]))
	samplePts <- generateSamplePointsInStratifiedLocationsWithMinDistance(numSamplePoints,rstDisturbance,plgOuter,minDistance,1000000)
	samplePts <- cbind(samplePts,valuesOfRasterAtPoints(rstX,samplePts))
	colnames(samplePts)[4] <- "rstX"
 	analysisFirstOrder(samplePts,rstX,rstDisturbance,4,3)
}
