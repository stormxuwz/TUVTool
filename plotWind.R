library(zoo)
library(dygraphs)

pere_wind <- read.table("/Users/WenzhaoXu/Developer/Triaxus/input/ldtm4h2011.txt",skip=2,sep="")
muskegon_wind <- read.table("/Users/WenzhaoXu/Developer/Triaxus/input/mkgm4h2011.txt",skip=2,sep="")


transData <- function(dataSet){
	dataSet$time <- strptime(sprintf("%d-%d-%d %d:%d",dataSet$V1,dataSet$V2,dataSet$V3,dataSet$V4,dataSet$V5),format = "%Y-%m-%d %H:%M",tz="UTC")
	names(dataSet)[6:8] <- c("WDIR","WSPD","GST")
	return(zoo(dataSet[,c("WDIR","WSPD","GST")],order.by = dataSet$time))
}

pere_wind <- transData(pere_wind)
muskegon_wind <- transData(muskegon_wind)

dygraph(window(pere_wind,start = as.POSIXct("2011-06-28 11:00:00",tz="UTC"),end = as.POSIXct("2011-06-28 23:59:00",tz="UTC")))
dygraph(window(muskegon_wind,start = as.POSIXct("2011-06-27 18:00:00",tz="UTC"),end = as.POSIXct("2011-06-28 02:00:00",tz="UTC")))
