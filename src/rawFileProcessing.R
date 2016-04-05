

filePreproessing <- function(data){
	useColumns <- c(1,3,4,5,6,7,8,9,20,25:59,78:204)
	names(data)[1] <- "time_second"
	data <- data[1:7057,]
	# data$Seabird_time <- strptime(data$BBE_time,format="%Y-%m-%d-%h-%m-%s")
	data$BBE_time <- paste("2011",data[,36],data[,37],data[,38],data[,39],data[,40],sep="-")
	data$BBE_time <- strptime(data$BBE_time,format="%Y-%m-%d-%H-%M-%S")
	# data  <- data[,useColumns]
}


fileProcessing <- function(){
	Man_fileNameList <- c("Manitowoc_01.csv","Manitowoc_02.csv","Manitowoc_03.csv","Manitowoc_04.csv","Manitowoc_05.csv")
	Mus_fileNameList <- c("Muskegon River_01.csv","Muskegon River_02.csv","Muskegon River_03.csv","Muskegon River_04.csv","Muskegon River_03a.csv","Muskegon River_04a.csv")
	Pere_fileNameList <- c("Pere Marquette River_01.csv","Pere Marquette River_02.csv","Pere Marquette River_03.csv","Pere Marquette River_04.csv","Pere Marquette River_05.csv","Pere Marquette River_06.csv")


	for(f in Man_fileNameList){
		# csvFilePreprocessing(f)
		anchorPointFinder(f)
	}


	for(f in Mus_fileNameList){
		# csvFilePreprocessing(f)
		anchorPointFinder(f)
	}


	for(f in Pere_fileNameList){
		# csvFilePreprocessing(f)
		anchorPointFinder(f)
	}
}


csvFilePreprocessing <- function(csvFileName,newFile=csvFileName){
	data <- read.csv(csvFileName,skip=58)
	names(data)[1] <- "timePassSecond"
	names(data)[78:81] <- c("Zdens","Zug","OvrSzd","Zug_less_5")
	
	lastIndex <- min(which(as.character(data$timePassSecond)=="# MEP particles=            "))-1
	print(lastIndex)
	available <- 1:lastIndex
	data <- data[available,]
	newdata <- data[,c(1,3:9,20,25:35,41:52,55:59,78:203)]
	newdata$BBE_time <- paste("2011",data[,36],data[,37],data[,38],data[,39],data[,40],sep="-")
	newdata$Phyto_time <- paste(data[,53],data[,54],sep=" ")
	print(head(newdata[,1:10]))
	write.csv(newdata,newFile,row.names=FALSE)
}	



anchorPointFinder <- function(csvFileName,range=NULL,auto=TRUE){
	data <- read.csv(csvFileName)
	data$n <- 1:nrow(data)
	
	if(!is.null(range)){
		data <- data[range[1]:range[2],]
	}

	seabird_minDepth_index_max <- which.max(data$depth)
	BBE_minDepth_index_max <- which.max(data[,"depth.1"])


	seabird_minDepth_index_min <- which.min(data$depth)
	BBE_minDepth_index_min <- which.min(data[,"depth.1"])


	plot(-depth~n,data,type="l")
	lines(-depth.1 ~ n,data,col="green")
	
	points(-data$depth[seabird_minDepth_index_min]~data$n[seabird_minDepth_index_min],col="red")
	points(-data[BBE_minDepth_index_min,"depth.1"]~data$n[BBE_minDepth_index_min],col="red")
	
	points(-data$depth[seabird_minDepth_index_max]~data$n[seabird_minDepth_index_max],col="blue")
	points(-data[BBE_minDepth_index_max,"depth.1"]~data$n[BBE_minDepth_index_max],col="blue")
	
	print("Red")
	print(data$n[seabird_minDepth_index_min])
	print(data$n[BBE_minDepth_index_min])

	print("Blue")
	print(data$n[seabird_minDepth_index_max])
	print(data$n[BBE_minDepth_index_max])

}



tmpPlot <- function(){
	oldResult <- "./testFile/allTriaxus.rds"
	allTriaxus <- readRDS(oldResult)
	myTriaxus <- allTriaxus[[1]]

	rgldata <- myTriaxus@resultData[,c("DO","longitude","latitude","depth")]

	rgldata <- arrange(rgldata,depth,latitude)
	x <- rgldata$depth
	y <- rgldata$latitude
	z <- rgldata$longitude
	varVal <- rgldata[,"DO"]


	coloPal <- colorNumeric(topo.colors(10), varVal)
	valColor <- coloPal(varVal)

	open3d()
	plot3d(x,y,z,alpha=0,xlab="x",ylab="y",zlab="z",box=F)
	par3d(windowRect = c(10, 10, 600, 600))
  	# x = seq(0, 1, len = nrow(z)), y = seq(0, 1, len = ncol(z)), z,
  	
	x_tmp <- unique(x)
	y_tmp <- unique(y)
	z_tmp <- t(matrix(z,length(y_tmp),length(x_tmp)))
	
	varVal_tmp <- t(matrix(varVal,length(y_tmp),length(x_tmp)))
	col_tmp <- t(matrix(valColor,length(y_tmp),length(x_tmp)))
	bgplot3d( suppressWarnings ( image.plot( legend.only=TRUE, zlim= range(varVal,na.rm=TRUE))))#legend

	persp3d(x_tmp,y_tmp,z_tmp,col=col_tmp,add=T,xlab="depth",ylab="Lat",zlab="Long")

}

checkCutoff <- function(fileName,cutoffIndex,new=FALSE){
	data <- read.csv(fileName)
	print(nrow(data))
	data <- data[cutoffIndex,]
	plot(-depth~c(1:nrow(data)),data=data,type="l")
	lines(-depth.1~c(1:nrow(data)),data=data,col="blue")

	if(new){
		newFileName <- paste("./testFile/new_",basename(fileName),sep="")
		write.csv(data,newFileName)
	}
}

discarded <- function(){
	

# checkCutoff("./testFile/newFile/Manitowoc_01.csv",800:7000)
# checkCutoff("./testFile/newFile/Manitowoc_02.csv",3170:9352)
# checkCutoff("./testFile/newFile/Manitowoc_03.csv",3140:6415)
# checkCutoff("./testFile/newFile/Manitowoc_04.csv")
# checkCutoff("./testFile/newFile/Manitowoc_05.csv")

# checkCutoff("./testFile/newFile/Muskegon River_01.csv",1200:7200)
# checkCutoff("./testFile/newFile/Muskegon River_02.csv",1880:8200)
# checkCutoff("./testFile/newFile/Muskegon River_03a.csv",150:6745)
# checkCutoff("./testFile/newFile/Muskegon River_04a.csv",1400:6614)

# checkCutoff("./testFile/Pere Marquette River_01.csv",200:7450)
# checkCutoff("./testFile/Pere Marquette River_02.csv",350:4750)
# checkCutoff("./testFile/Pere Marquette River_03.csv",2200:8637)
# checkCutoff("./testFile/Pere Marquette River_04.csv",1:8299)
# checkCutoff("./testFile/Pere Marquette River_05.csv",3030:7999)
# checkCutoff("./testFile/Pere Marquette River_06.csv",360:5100)
}