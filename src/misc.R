scaleFactor <- function(x,y){ 
	scaleSetting <- list()
	scaleSetting$x <- c(min(x),max(x)-min(x))
	scaleSetting$y <- c(min(y),max(y)-min(y))
	return(scaleSetting)
}

applyScale <- function(x,y,scaleSetting){
	return(list(
		sx = (x - scaleSetting$x[1])/scaleSetting$x[2],
		sy = (y - scaleSetting$y[1])/scaleSetting$y[2])
	)
}