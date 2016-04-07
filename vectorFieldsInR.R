createVectorDivergence <- function(rstX)
{
	matX <- as.matrix(rstX@data)
	dim(matX) <- x@grid@cells.dim
	n <- dim(matX)
	z <- rep(complex(real=0,imaginary=0),n[1]*n[2])
	dim(z) <- n
	for (i in 1:(n[1]-1))
	{
		for (j in 1:(n[2]-1))
		{
			print(paste(i,j,sep=" "))
			if (i==1 & j==1) { z[i,j] = 0 }
			if ( i>1 & j==1) { z[i,j] = z[i-1,j] + complex(imaginary=(matX[i+1,j]-matX[i-1,j])) }
			if (i==1 & j>1 ) { z[i,j] = z[i,j-1] + complex(real=(matX[i,j+1]-matX[i,j+1])) }
			if ( i>1 & j>1 )
			{
				z[i,j] = 0.5*z[i,j-1] + 0.5*z[i-1,j] + complex(real=(matX[i,j+1]-matX[i,j+1]),imaginary=(matX[i+1,j]-matX[i-1,j]))
			}
		}
	}
	z
}