atanc <- function(x)
{
	if (!is.complex(x)) { stop("Argument is not complex.") }
	atan2(Im(x),Re(x))
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

histogramBreaksSturges <- function(x,n)
{
	xMin <- min(x)
	xMax <- max(x)
	#n <- length(x)
	nBins <- ceiling(log2(n)+1)
	binWidth <- (xMax-xMin)/nBins
	bins <- cbind(xMin + binWidth*(0:(nBins-1)),xMin+binWidth*(1:nBins),xMin+binWidth*(0:(nBins-1) + 0.5))
	colnames(bins) <- c("xMin","xMax","xMid")
	list(xMin=xMin,xMax=xMax,n=n,nBins=nBins,bins=bins)
}

isDivisibleBy <- function(x,y)
{
	x %% y == 0
}

mag <- function(x)
{
	if (!is.complex(x)) { stop("Argument is not complex.") }
	sqrt(Im(x)^2+Re(x)^2)
}

observationsRequiredForBinomialDifferentiation <- function(p1,p2,alpha=0.05)
{
	if (p1<0 | p1>1 | p2<0 | p2>1) { stop("Probabilities need to be between zero and one.") }
	if (p1==p2) 
	{ 
		n <- Inf 
	} else {
		pa <- min(p1,p2)
		pb <- max(p1,p2)
		n <- 1
		while ((1-pbinom(floor(pb*n),n,pa))>alpha) 
		{ 
			n <- n + 1 
			print(paste(n,(1-pbinom(floor(pb*n),n,pa)),sep=" : "))
		}
	}
	n
}

rmLike <- function(expression)
{
	print("Removing...")
	print(ls(pattern=expression))
	rm(list=ls(pattern=expression))
}