createVectorDivergence <- function(rstX,rstY,betamax=-1)
{
	matX <- as.matrix(rstX@data)
	dim(matX) <- rstX@grid@cells.dim
	n <- dim(matX)
	matY <- as.matrix(rstY@data)
	dim(matY) <- rstY@grid@cells.dim
	if (dim(matX) != dim(matY)) { stop("The two rasters must have identical dimensions") }
	matXx <- as.matrix(rep(0),n[1]*n[2]))
	dim(matXx) <- n
	matXy <- as.matrix(rep(0),n[1]*n[2]))
	dim(matXy) <- n
	z <- as.matrix(rep(complex(real=0,imaginary=0),n[1]*n[2]))
	dim(z) <- n	
	for (i in 2:(n[1]-1))
	{
		for (j in 2:(n[2]-1))
		{
			matXx[i,j] <- matX[i,j+1] - matX[i,j-1]
			matXy[i,j] <- matX[i+1,j] - matX[i-1,j]
			for (p in 1:betamax)
			{
				for (q in 1:betamax)
				{
					z[i,j] <- z[i,j] + complex(real=matY[i,j]/(matXx[i,j]-p/q*matXy[i,j]),imaginary=matY[i,j]/(matXy[i,j]-q/p*matXx[i,j]))
				}
			}
		}
	}
	vectorField()
	z
}