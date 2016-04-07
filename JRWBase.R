rmLike <- function(expression)
{
	print("Removing...")
	print(ls(pattern=expression))
	rm(list=ls(pattern=expression))
}