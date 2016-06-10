rmLike <- function(expression)
{
	print("Removing...")
	print(ls(pattern=expression))
	rm(list=ls(pattern=expression))
}

histogramBreaksSturges <- function(x,n)
{
	xMin <- min(x)
	xMax <- max(x)
	#n <- length(x)
	nBins <- ceiling(log2(n)+1)
	print(class(x))
	print(class(xMin))
	dx <- (xMax-xMin)/nBins
	bins <- cbind(xMin + dx*(0:(nBins-1)),xMin+dx*(1:nBins),xMin+dx*(0:(nBins-1) + 0.5))
	colnames(bins) <- c("xMin","xMax","xMid")
	list(xMin=xMin,xMax=xMax,dx=dx,n=n,nBins=nBins,bins=bins)
}

betaSuccessProbDistribution <- function(successes=0, failures=0, step=0.01) 
{
	dbeta(seq(0,1,step),shape1=successes+1,shape2=failures+1)
}

betaSuccessQuantiles <- function(p=c(0.025,0.975),successes=0,failures=0)
{
	qbeta(p,shape1=successes+1,shape2=failures+1)
}

betaSuccessMean <- function(successes=0,failures=0)
{
	(successes+1)/(successes+failures+2)
}

betaSuccessMode <- function(successes=0,failures=0)
{
	#if (successes + failures == 0) { stop("No mode for b(0,0)")}
	if (successes + failures == 0)
	{
		result <- NA
	} else {
		result <- successes/(successes+failures)
	}
	result
}

betaSuccessMedian <- function(successes=0,failures=0)
{
	betaSuccessQuantiles(p=0.5,successes,failures)
}

isDivisibleBy <- function(x,y)
{
	x %% y == 0
}

numRows <- function(x)
{
	if (is.matrix(x))
	{
		y <- dim(x)[1]
	} else {
		y <- 1
	}
	y
}

observationsRequiredForBinomialDifferentiation <- function(p1,p2,alpha=0.05,countAsInf=1000000)
{
	if (p1<0 | p1>1 | p2<0 | p2>1) { stop("Probabilities need to be between zero and one.") }
	if (p1==p2) 
	{ 
		n <- Inf 
  	} else {
		pa <- min(p1,p2)
		pb <- max(p1,p2)
		n <- floor(1/(pb-pa))
		while ((1-pbinom(floor(pb*n),n,pa))>alpha) 
		{ 
			n <- n + 1 
			if (n > countAsInf) 
			{ 
				n <- Inf
				break
			}
		}
		#print(paste(n,pa,pb,(1-pbinom(floor(pb*n),n,pa)),sep=" : "))
	}
	n
}