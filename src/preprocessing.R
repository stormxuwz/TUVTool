library(dplyr)
library(lubridate)



alignment <- function(timePassSecond,BBE_time,seabirdIndex,bbeIndex){
	anchorBBETime <- BBE_time[bbeIndex]
	anchorSeabirdTime <- timePassSecond[seabirdIndex]

	newSeabirdTimeDiff <- timePassSecond-anchorSeabirdTime
	newSeabirdTime <- anchorBBETime+newSeabirdTimeDiff
	# newSeabirdTime <- BBE_time[bbeIndex]

	# timeAdjust <- (seq(1:length(BBE_time))-seabirdIndex)*0.5
	# newSeabirdTime <- newSeabirdTime+timeAdjust
	return(newSeabirdTime)
}


preprocessing <- function(Triaxus){
	seabirdIndex <- Triaxus@SeabirdAnchorIndex
	bbeIndex <- Triaxus@BBEAnchorIndex
	
	Triaxus@rawData <- initProcess(Triaxus@rawData)
	BBE_name <- Triaxus@config$BBE_name
	Seabird_name <- Triaxus@config$Seabird_name

	# print(summary(Triaxus@rawData))
	geoData <- Triaxus@rawData[,c("latitude","longitude","distance","UTC")]
	Seabird_data <- Triaxus@rawData[,c(Seabird_name,"Seabird_depth")]
	
	if("Zug" %in% names(Triaxus@rawData)){
		LOPC_data <- Triaxus@rawData[,c("Zug","Zdens")] # for back compatibility 
	}else{
		LOPC_data <- Triaxus@rawData[,Triaxus@config$LOPC_name]
	}
	Triaxus@cleanData <- cbind(geoData,Seabird_data,LOPC_data)
	
	n <- nrow(Triaxus@cleanData)
	Triaxus@cleanData[,BBE_name] <- NA

	if(seabirdIndex>0 & bbeIndex>0){
		print("Aligning BBE Data")
		BBE_data <- Triaxus@rawData[,c(BBE_name)]
		BBE_time <- strptime(Triaxus@rawData$BBE_time,format="%Y-%m-%d-%H-%M-%S")
		newSeabirdTime <- alignment(Triaxus@rawData$timePassSecond,BBE_time,seabirdIndex,bbeIndex)
		Triaxus@cleanData[,BBE_name] <- linearAlignment(newSeabirdTime,BBE_time,BBE_data,rule=1)
	}
	
	# if(length())
	# startIndex <- ifelse(Triaxus@Seabird_cutoff[1]<1,1,Triaxus@Seabird_cutoff[1])
	# endIndex <- ifelse(Triaxus@Seabird_cutoff[2]<1,n,Triaxus@Seabird_cutoff[2])

	# Triaxus@cleanData <- Triaxus@cleanData[startIndex:endIndex,]
	Triaxus@cleanData <- Triaxus@cleanData[Triaxus@Seabird_cutoff,]
	
	Seabird_separation <- depth_separation(Triaxus@cleanData$Seabird_depth)
	Seabird_nodes <- Seabird_separation[[2]]

	finalStartNodes <- Seabird_nodes[2]
	finalEndNodes <- Seabird_nodes[length(Seabird_nodes)-1]
	Triaxus@cleanData$direction <- Seabird_separation[[1]]
	Triaxus@cleanData$depth <- Triaxus@cleanData$Seabird_depth
	
	Triaxus@cleanData <- Triaxus@cleanData[finalStartNodes:finalEndNodes,]

	# Adding salinity and density estimation.
	tmpData <- Triaxus@cleanData[,c("latitude","longitude","conductivity")]
	tmpData$temperature <- Triaxus@cleanData$Seabird_temperature
	# densityData <- cal_density(tmpData)
	# Triaxus@cleanData <- cbind(Triaxus@cleanData,densityData)

	
	Triaxus@numCycle <- (length(Seabird_separation[[2]])-1)/2-1
	return(Triaxus)
}



initProcess <- function(dataSet){
	# dataSet <- na.omit(dataSet)
	print(summary(dataSet))
	dataSet <- dplyr::select(dataSet,-c(latitude,longitude))
	dataSet <- dplyr::rename(dataSet,BBE_depth=depth.1,Seabird_depth=depth,Seabird_temperature=temp,
		BBE_temperature=temp.1,distance=Distance,DO=DO.43.mg.L,DOsat=DO43...sat,latitude=DDLat,longitude=DDLong,conductivity = cond)
	
	if("Zdens" %in% names(dataSet)){
		dataSet$Zdens <- as.numeric(as.character(dataSet$Zdens))
	}

	if("Zug" %in% names(dataSet)){
		dataSet$Zug <- as.numeric(as.character(dataSet$Zug))
	}
	return(dataSet)
}



linearAlignment <- function(targetIndex,sourceIndex,sourceDataSet,rule=2){
	targetDataSet <- list()
	for(var in names(sourceDataSet))
		targetDataSet[[var]] <- approx(x=sourceIndex,y=sourceDataSet[,var],xout=targetIndex,rule=rule)$y

	return(data.frame(targetDataSet))
}



repeatLast <- function(x){
	n <- length(x)
	return(c(x,x[n]))
}

depth_separation <- function(depth,...){
	n <- length(depth)
	depth_smooth <- -smooth.spline(depth,spar=0.01)$y
	N <- length(depth_smooth)
	if(N!=n){
		warning("Spline generates different number of data")
	}

	diff_D <- diff(depth_smooth)
	D1 <- diff_D[1:(N-1)]
    D2 <- diff_D[2:N]

    transition_points <- c(1,(which(D1*D2<0)+1),N)

    directionVector <- c()
    previousNode <- 1

    for(i in transition_points[-1]){
    	if(depth_smooth[i]>depth_smooth[i-1]){# this is a peak point
    		directionVector <- c(directionVector,rep(1,i-previousNode))
    	}
    	else{# this is a valley point
    		directionVector <- c(directionVector,rep(-1,i-previousNode))
    	}
    	previousNode <- i
    }
    directionVector <- repeatLast(directionVector)
    return(list(directionVector,transition_points))
}	

init_filter <- function(x,varName){
	  threshold <- list() 
    threshold["DOsat"] <- 130
    threshold["DO"] <- 18
    threshold["Seabird_temperature"] <- 40
    threshold["Zdens"] <- quantile(x,0.995,na.rm=T)
    threshold["Zug"] <- quantile(x,0.995,na.rm=T)
    
    for(zooplanktonBinName in names(config$ZugBins)){
    	threshold[[paste0("Zdens_",zooplanktonBinName)]] <- quantile(x,0.995,na.rm=T)
    	threshold[[paste0("Zug_",zooplanktonBinName)]] <- quantile(x,0.995,na.rm=T)
    }
    
    if(varName %in% names(threshold))
      outlier <- (x>threshold[varName])
    else
      outlier <- rep(0,length(x))
    return(outlier)
}

spatialOutlier <- function(spData,x_y_ratio,nbRange,threshold){
	require(spdep)
	distance_range <- range(spData$distance)
	depth_range <- range(spData$depth)

	xycoords <- cbind(spData$distance*x_y_ratio, spData$depth)
	nb <- dnearneigh(xycoords,0,nbRange)

	hset <- c() # the distance to the median
	outlier <- rep(0,nrow(spData))
	avail <- outlier
	for(i in 1:length(nb)){
		neighborPoints <- nb[[i]]
		if(length(neighborPoints)>10 & !is.na(spData[i,"res"])){
			hset=c(hset,spData[i,"res"]-median(spData[neighborPoints,"res"],na.rm=TRUE))
			avail[i] <- 1
		}
	}

	robustmean <- median(hset,na.rm=T)
	robustsd <- sd(hset,na.rm=T)
	hset <- (hset-robustmean)/robustsd # normalize 
	outlier[avail==1] <- ifelse(abs(hset)>threshold,1,0)
	return(outlier)
}


cal_density <- function(data){
	# data is a data frame that contains latitude,longitude, temperature, conductivity and pressure
	require(oce)
	salinity <- swSCTp(conductivity=data$conductivity/1000,temperature=data$temperature,pressure=data$pressure,"mS/cm")
	density <- swRho(salinity,data$temperature,data$pressure,longitude = data$longitude,latitude = data$latitude)
	return(data.frame(salinity = salinity, density = density))
}