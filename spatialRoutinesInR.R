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
			#print(paste("x=",x))
			#print(paste("y=",y))
			#print(paste("crds_x=",str(crds[,1])))
			#print(paste("crds_y=",str(crds[,2])))
			#print(point.in.polygon(x,y,crds[,1],crds[,2]))
			result <- result + point.in.polygon(x,y,crds[,1],crds[,2])
		}
	}
	#print (result)
	result
}
generateSamplePointsInPolygonWithMinDistance <- function(numPts,outerPolygon,minDistance=0,maxIter=10000)
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
	#print(s1[1:10])
	#print(s2[1:10])
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
			#print (paste("Got one... ",as.integer(z*numBreaks/maxDistance),rstXcol[s1,3],rstXcol[s2,3],sep=" "))
			smpPts <- rbind(smpPts,c(as.integer(z*numBreaks/maxDistance)+1,rstXcol[s1[i],3],rstXcol[s2[i],3]))
			#print(smpPts)
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
	#image(df["z"])
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
	matX <- as.matrix(rstX@data)
	dim(matX) <- rstX@grid@cells.dim
	smpPts[,1] <- floor((smpPts[,1] - rstX@bbox[1,1])/rstX@grid@cellsize[1]) + 1
	smpPts[,2] <- floor((smpPts[,2] - rstX@bbox[2,1])/rstX@grid@cellsize[2]) + 1
	n <- dim(smpPts)[1]
	result <- numeric(n)
	for (i in 1:n)
	{
		if (smpPts[i,1]>0 & smpPts[i,2]>0 & smpPts[i,1]<=dim(matX)[1] & smpPts[i,2]<=dim(matX)[2])
		{
			result[i] <- matX[smpPts[i,1],smpPts[i,2]]
		} else {
			result[i] <- NA
			print (paste(smpPts[i,]," is outside ",rstX@bbox))
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

analysisFirstOrder <- function(x,params,polys)
{
	header <- list(length(params))
	names(header) <- colnames(x[,params])
	for (i in 1:length(header))
	{
		header[[i]] <- histogramBreaksSturges(x[params[i]])
	}
	analyses <- list(length(polys))
	names(analyses) <- colnames(x[,polys])
	for (j in 1:length(polys))
	{
		buildTable <- matrix(nrow=0,ncol=7)
		colnames(buildTable) <- c("Successes","Failures","Pmean","P95-","P95+","P99-","P99+")
		for (i in 1:length(params))
		{
			s <- dim(x[x[,params[i]]>=header[[i]]$xMin & x[,params[i]]<header[[i]]$xMax & x[,polys[j]]==1,])[1]
			f <- dim(x[x[,params[i]]>=header[[i]]$xMin & x[,params[i]]<header[[i]]$xMax & x[,polys[j]]==0,])[1]
			buildTable <- rbind(rowTable,c(s,f,betaSuccessMean(s,f),betaSuccessQuantiles(c(0.025,0.975,0.005,0.995),s,f)))
		}
		analyses[[j]] <- buildTable
	}
	list(header=header,analyses=analyses)
}







setClass("rasterSequenceSample",representation(boundingPolygon="SpatialPolygonsDataFrame",minDistance="numeric",samplePoints="matrix",results="matrix"))
setClass("rasterObject",representation(source="character",raster="SpatialGridDataFrame",variogramList="list",minDistance="numeric"),prototype(source="",raster=simpleRaster(),variogramList=as.list(x=NULL),minDistance=9999))
setClass("rasterSequence",representation(numRasters="numeric",rasterObjects="list",rasterSequenceSample="rasterSequenceSample"))