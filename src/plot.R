library(ggplot2)
library(ggmap)
library(gridExtra)
plot_raw <- function(Triaxus,var){
	visdata <- Triaxus@cleanData[,c(var,"distance","depth")]
	colorRange <- range(c(range(Triaxus@cleanData[,var],na.rm = TRUE),range(Triaxus@resultData[,var],na.rm = TRUE)))

	names(visdata)[1] <- "value"
	p <- ggplot(visdata)+geom_point(aes(distance,-depth,color=value))+scale_color_gradientn(colours = topo.colors(10),name=var,limits=colorRange)
	p <- p+xlab("Distance (km)")+ylab("Depth (m)")
	return(p)
}

plot_2d <- function(Triaxus,var){
	visdata <- Triaxus@resultData[,c(var,"distance","depth")]
	colorRange <- range(c(range(Triaxus@cleanData[,var],na.rm = TRUE),range(Triaxus@resultData[,var],na.rm = TRUE)))
	print(colorRange)
	names(visdata)[1] <- "value"
	p <- ggplot(visdata)+geom_tile(aes(distance,-depth,fill=value))+scale_fill_gradientn(colours = topo.colors(10),name=var,limits=colorRange)
	p <- p+geom_point(aes(distance,-depth),data=Triaxus@cleanData,alpha=0.5)
	p <- p+xlab("Distance (km)")+ylab("Depth (m)")
	return(p)
}

plot_hotspot <- function(Triaxus,var){
	visdata <- Triaxus@hotspotData[,c(var,"distance","depth")]
	names(visdata)[1] <- "value"
	p <- ggplot(visdata)+geom_tile(aes(distance,-depth,fill=value))+scale_fill_gradient2(low="blue4",high="red",mid="white",midpoint=0,name=var)
	p <- p+xlab("Distance (km)")+ylab("Depth (m)")
	return(p)
}


plot_3d_base <- function(dataList,isFactor=FALSE,hotspot=FALSE,...){
	# dataList is a list where each element is a data frame with "val","longitude","latitude","depth"
	additionalArgs <- list(...)
	if(is.null(additionalArgs$transmatrix)){
		transmatrix=matrix(c(0.01585434,-0.68476230,-0.72859418, 0.00000000,0.4700032, 0.6482660,-0.5990396,0.0000000,0.8825220,-0.3329444,0.3321184,0,0,0,0,1),nrow=4,ncol=4)	
	}else{
		transmatrix=list(...)$transmatrix
	}

	allVar <- c()
	geoLocation <- data.frame()
	
	for(subData in dataList){
		allVar <- c(allVar,subData[,"val"])
		geoLocation <- rbind(geoLocation,subData[,c("longitude","latitude")])
	}
	
	print(summary(allVar))
	if(isFactor){
			if(hotspot){
				coloPal <- colorFactor(palette = c("blue","white","Red"),as.factor(c(-1,0,1)),na.color =NA)
			}
			else{
				coloPal <- colorFactor(topo.colors(5), allVar,na.color = NA)
			}
			# print("coloPal is a factor coloPal")
		}else{
			varRange <- range(allVar,na.rm=TRUE)
			coloPal <- colorNumeric(topo.colors(10), c(floor(varRange[1]),ceiling(varRange[2])))
	}

	open3d()
	for(i in 1:length(dataList)){
		rgldata <- dataList[[i]]

		x_tmp <- unique(rgldata$depth)
		y_tmp <- rgldata$latitude[1:(nrow(rgldata)/length(x_tmp))]
		
		# rgldata <- dplyr::arrange(rgldata,depth,latitude)

		x <- rgldata$depth
		y <- rgldata$latitude
		z <- rgldata$longitude

		varVal <- rgldata[,"val"]
		naIndex <- is.na(varVal)
		if(isFactor){
			varVal <- as.factor(varVal)
		}

		if(i==1){
		  plot3d(x,y,z,alpha=0,xlab="",ylab="",zlab="",box=F)
		  aspect3d(0.8,1,0.5)
		}
		
		# print(summary(varVal))
		valColor <- coloPal(varVal)
		# z[naIndex] <- NA
		z_tmp <- t(matrix(z,length(y_tmp),length(x_tmp)))
		col_tmp <- t(matrix(valColor,length(y_tmp),length(x_tmp)))
		surface3d(x_tmp,y_tmp,z_tmp,col=col_tmp,add=T,xlab="depth",ylab="Lat",zlab="Long")
		# persp3d(x_tmp,y_tmp,z_tmp,col=col_tmp,add=T,xlab="depth",ylab="Lat",zlab="Long")
	}
	rgl.viewpoint(zoom = .8)
	par3d(userMatrix=transmatrix,windowRect=c(720,0,1920,1200),cex=0.0001)


	# Create the map cap
	bbox <- ggmap::make_bbox(geoLocation$longitude,geoLocation$latitude,f=0.02)
	triaxusPath.map <- get_map(location = bbox, maptype="watercolor", source="stamen",zoom=9)
	coor=attr(triaxusPath.map,"bb")
	numY=dim(triaxusPath.map)[1]
	numX=dim(triaxusPath.map)[2]

	byX=(coor$ur.lon-coor$ll.lon)/(numX-1)
	byY=(coor$ur.lat-coor$ll.lat)/(numY-1)

	long=seq(from=coor$ll.lon,to=coor$ur.lon,by=byX)
	lat=seq(from=coor$ll.lat,to=coor$ur.lat,by=byY)
	cols=as.matrix(triaxusPath.map)
	
	geox=rep(0,length(long))  # Plot geo map at the surface 
	geoy=lat
	geoz=matrix(rep(c(long)),nrow=length(long),ncol=length(lat))
	geocol=pracma::rot90(cols,k=-1)

	surface3d(geox,geoy,geoz,col=geocol,add=T,lit=F,alpha=0.3)

	if(!is.null(additionalArgs$riverMouth)){
		# plot the dot of riverMonth
		# riverMonth: list(river1,river2),river1=c(c1,c2)
		print("print riverMouth")
		for(i in 1:length(riverMouth))
    		plot3d(0,additionalArgs$riverMouth[[i]][2],additionalArgs$riverMouth[[i]][1],add=T,size=10,col="red")
	}
}

plot_3d_value <- function(allTriaxus,var,...){
	dataList <- list()

	for(myTriaxus in allTriaxus){
		dataList[[myTriaxus@pathName]] <- myTriaxus@resultData[,c(var,"longitude","latitude","depth")]
		names(dataList[[myTriaxus@pathName]])[1] <- "val"
	}

	plot_3d_base(dataList,FALSE,FALSE,...)
}

plot_3d_hotspot <- function(allTriaxus,var,...){
	dataList <- list()
	for(myTriaxus in allTriaxus){
		rawSub <- myTriaxus@hotspotData[,c(var,"longitude","latitude","depth")]
		rawHotspot <- rawSub[,var]
		rawSub[,var] <- 0
		rawSub[,var] <- ifelse(rawHotspot<(-3),-1,rawSub[,var])
		rawSub[,var] <- ifelse(rawHotspot>3,1,rawSub[,var])
		rawSub[,var] <- as.factor(rawSub[,var])
		dataList[[myTriaxus@pathName]] <- rawSub
		names(dataList[[myTriaxus@pathName]])[1] <- "val"
	}
	plot_3d_base(dataList,TRUE,TRUE,...)
}

plot_3d_clustering <- function(allTriaxus,K,...){
	clusterName <- paste("cluster_",K,sep="")
	dataList <- list()
	for(myTriaxus in allTriaxus){
		dataList[[myTriaxus@pathName]] <- myTriaxus@grid[,c("longitude","latitude","depth")]
		dataList[[myTriaxus@pathName]][,"val"] <- as.factor(myTriaxus@clusteringResults$clusteringIndex[[clusterName]])
	}

	plot_3d_base(dataList,TRUE,FALSE,...)
}

plot_boxplot <- function(allTriaxus,variables,K){
	totalData <- data.frame()
	clusterName <- paste("cluster_",K,sep="")

	for(myTriaxus in allTriaxus){
		mydata <- myTriaxus@resultData[,variables]
		mydata$cluster <- myTriaxus@clusteringResults$clusteringIndex[[clusterName]]
		totalData <- rbind(totalData,mydata)
	}

	totalData$cluster <- factor(totalData$cluster)
	plot_theme=theme(axis.title.x=element_blank(),axis.text.x=element_text(size=6.5,face="bold"),axis.title.y=element_text(size=6.5,face="bold"),
      axis.text.y=element_text(size=6.5,face="bold"),plot.margin = unit(c(0.1,0.1,0.1,0),"cm"))

	plot.list <- lapply(variables, function(x){
		subData <- na.omit(totalData[,c("cluster",x)]);qplot(x = subData$cluster, y=subData[, x])+ylab(config$varUnit[x])+geom_boxplot(outlier.size=0.5,na.rm=TRUE)+plot_theme
	})

	args.list <- c(plot.list,list(nrow=2,ncol=length(plot.list)/2))
	do.call(grid.arrange, args.list)
}

