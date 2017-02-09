#### This is the main part for krig Interpolation, with following functions


crossValidataion <- function(data,K=5,R=100,int_method, det_method, plotid, maxdist){
	require(cvTools)
	n <- nrow(data)
	cvIndex <- cvFolds(n,K = K,R = R)
	cvFoldIndex <- cvIndex$which

	RMS_error <- matrix(-99, nrow = R, ncol = K)
	max_error <- matrix(-99,nrow = R, ncol = K)
	dataRange <- matrix(-99,nrow = R, ncol = K)
	meanError <- matrix(-99, nrow = R, ncol = K)
	ratioError <- matrix(-99, nrow = R, ncol = K)
	dataMean <- matrix(-99, nrow = R, ncol = K)
	dataMean_whole <- matrix(-99, nrow = R, ncol = K)
	dataRange_whole <- matrix(-99, nrow = R, ncol = K)
	dataMedian <- matrix(-99, nrow = R, ncol = K)
	dataMedian_whole <- matrix(-99, nrow = R, ncol = K)

	for(i in 1:R){
		for(fold in 1:K){
			print(paste("cv progress:",i,fold))
			testIndex <- cvIndex$subsets[cvFoldIndex==fold,i]
			trainIndex <- c(1:n)[-testIndex]

			trainData <- data[trainIndex,]
			testData <- data[testIndex,]
			testData$availableIndex <- TRUE
			capture.output(cvRes <- interpolation_sub(trainData,testData[,c("distance","depth","availableIndex")],int_method,det_method,plotid = paste(plotid,"_cvwhole",sep=""),maxdist = maxdist))
			
			dataRange[i,fold] <- diff(range(testData[,1]),na.rm=TRUE)
			RMS_error[i,fold] <- sqrt(mean((cvRes-testData[,1])^2,na.rm = TRUE))
			max_error[i,fold] <- max(abs(cvRes-testData[,1]),na.rm = TRUE)
			meanError[i,fold] <- mean(cvRes-testData[,1], na.rm = TRUE)
			ratioError[i, fold] <-mean(abs(cvRes-testData[,1]/testData[,1]), na.rm = TRUE) 
			
			dataMean[i, fold] <- mean(testData[,1],na.rm = TRUE)
			dataMedian[i,fold]  <- median(testData[,1],na.rm = TRUE)

			dataMean_whole[i,fold] <- mean(data[,1],na.rm=TRUE)
			dataRange_whole[i,fold] <- diff(range(data[,1],na.rm=TRUE))
			dataMedian_whole[i, fold] <- median(data[,1],na.rm = TRUE)
			# print(summary(cvRes-testData[,1]))
		}
	}
	print(RMS_error/dataRange)
	print(max_error/dataRange)
	# print("CV ERROR:")
	# sleep(10000)
	# print(mean(RMS_error))
	# print(mean(max_error))
	return(list(cvIndex = cvIndex, 
		RMS_error = RMS_error, 
		dataRange = dataRange, 
		max_error = max_error,
		dataMean = dataMean,
		dataMean_whole = dataMean_whole,
		dataRange_whole = dataRange_whole))
}


interpolation_main <- function(Triaxus,int_method,det_method=NULL,cv=FALSE,...){
	
	# create the grid
	Triaxus@grid <- createGrid(
		x=Triaxus@cleanData$distance,
		y=Triaxus@cleanData$depth,
		dx=Triaxus@config$gridSize[1],
		dy=Triaxus@config$gridSize[2],
		longitude=Triaxus@cleanData$longitude,
		latitude=Triaxus@cleanData$latitude)

	
	Triaxus@resultData <- Triaxus@grid

	# maxdist <- ifelse(config$krigingRange/Triaxus@numCycle>0.15,config$krigingRange/Triaxus@numCycle,0.15)
	# maxdist <- 0.33
	# maxdist <- config$maxdist
	# maxdist <- NULL
	maxdist <- min(max(config$maxdist,3/Triaxus@numCycle),1)

	for(var in Triaxus@config$interestVar){
		print(var)
		if(Triaxus@separate == FALSE){
			# | var == "BAT"
			print("no separating")
			plotid <- paste(Triaxus@pathName,var,sep="_")
			outlierIndex <- init_filter(Triaxus@cleanData[,var],var)
			pred <- interpolation_sub(Triaxus@cleanData[outlierIndex<1,c(var,"distance","depth")],
				Triaxus@grid,int_method,det_method,plotid = paste(plotid,"_whole",sep=""),maxdist = maxdist)
			Triaxus@resultData[,var] <- ifelse(pred>0,pred,0)
		}else{
			plotid <- paste(Triaxus@pathName,var,sep="_")
	  		outlierIndex <- init_filter(Triaxus@cleanData[,var],var)
	    	
	    	downData <- subset(Triaxus@cleanData[outlierIndex<1,],direction==-1)[,c(var,"distance","depth")]
	    	upData <- subset(Triaxus@cleanData[outlierIndex<1,],direction==1)[,c(var,"distance","depth")]
	    	
	    	if(cv & config$meta){
				dir.create(file.path(config$outputFolder, "meta/cv"), showWarnings = FALSE)
	    		print("starting cross validation")
	    		# K is the fold of the CV, R is the times for CV
	    		saveRDS(crossValidataion(downData,K=3,R=10,int_method,det_method,plotid = paste(plotid,"_down",sep=""),maxdist = maxdist),
	    			file = paste(config$outputFolder,"/meta/cv/",plotid,"_cv_down.rds",sep=""))
	    		saveRDS(crossValidataion(upData,K=3,R=10,int_method,det_method,plotid = paste(plotid,"_up",sep=""),maxdist = maxdist),
	    			file = paste(config$outputFolder,"/meta/cv/",plotid,"_cv_up.rds",sep=""))
	    	}
	    	# else{
    		pred_down <- interpolation_sub(downData,Triaxus@grid,int_method,det_method,plotid = paste(plotid,"_down",sep=""),
    			maxdist = maxdist)
    		pred_up <- interpolation_sub(upData,Triaxus@grid,int_method,det_method,plotid = paste(plotid,"_up",sep=""),
    			maxdist = maxdist)
    		finalPrediction <- (pred_down+pred_up)/2.0
    		Triaxus@resultData[,var] <- ifelse(finalPrediction>0,finalPrediction,0)
	    	# }
		}		
	}
	# Triaxus@grid$available <- !as.logical(apply(is.na(Triaxus@resultData[,Triaxus@config$interestVar]),1,sum)) # update
	return(Triaxus)
}


interpolation_sub <- function(dataSet,grid,int_method,det_method,...){
	pred <- rep(NA,nrow(grid))
    availableIndex <- grid$available
    # print(list(...))
    dataSet <- na.omit(dataSet)
    if(nrow(dataSet)==0){
    	print("no variable")
      return(pred)
    }

    if(int_method=="tps"){
        pred[availableIndex] <- interpolation_tps(dataSet,grid[availableIndex,],...)
        return(pred)
    }

    detrending_result <- detrending(dataSet,grid[availableIndex,],method=det_method)
    dataSet$res <- detrending_result$res

    # save the TPS trend data
    if(config$meta){
		dir.create(file.path(config$outputFolder, "meta/detrending"), showWarnings = FALSE)
    	saveRDS(list(dataset = dataSet,df = detrending_result$df),
    	file = paste0(config$outputFolder,"/meta/detrending/",list(...)$plotid,"_trend.rds"))
    }
    
    outlierIndex <- spatialOutlier(dataSet,1,0.75,4)
    print(paste("outlier Num:",sum(outlierIndex)))
    dataSet <- dataSet[outlierIndex<1,]
	
	if(int_method=="krige"){
  		pred_res <- interpolation_krig(dataSet,grid[availableIndex,],idw=FALSE,...)
	} 
	else if(int_method=="idw"){
  		pred_res <- interpolation_krig(dataSet,grid[availableIndex,],idw=TRUE,...)
	}
	else{
		stop("wrong method")
	}

    pred[availableIndex] <- detrending_result$trendSurface+pred_res

    return(pred)
}



createGrid <- function(x,y,dx,dy,longitude,latitude){
    require(dismo)
		xRange <- range(x)
		yRange <- range(y)
		# yRange <- quantile(y,c(0.05,0.8))
		startDepth=yRange[1]-yRange[1]%%dy+dy
		grid <- expand.grid(distance=seq(from=xRange[1],to=xRange[2],by=dx),
			depth=seq(from=startDepth,to=yRange[2],by=dy))

		grid$latitude <- approx(y=latitude,x,xout=grid$distance,rule=2)$y
		grid$longitude <- approx(y=longitude,x,xout=grid$distance,rule=2)$y

		# do convex hull
    convexHullModel<-convHull(data.frame(distance=x,depth=y))
    grid$available <- as.logical(predict(convexHullModel,grid[,c("distance","depth")]))
		return(grid)
}


interpolation_tps<-function(spData,grid,...){
  	require(fields)
	  
  	subSample <- seq(1,nrow(spData),by=1)
  	tpsModel <- Tps(spData[subSample,c("distance","depth")],spData[subSample,1])
  	# print(summary(tpsModel))
  	print("error")
  	print(max(tpsModel$residuals/spData[subSample,1]))
  	
  	pred <- predict(tpsModel,grid[,c("distance","depth")])

  	spData$res <- tpsModel$residuals
  	
  	plotid <- list(...)$plotid

  	if(config$meta){
  		dir.create(file.path(config$outputFolder, "meta/tps"), showWarnings = FALSE)
  		saveRDS(list(tpsModel,spData),file = paste0(config$outputFolder,"/meta/tps/",plotid,"_tpsModel.rds"))
  		pdf(file=paste(config$outputFolder,"/meta/tps/",plotid,"_tpsInt.pdf",sep=""))
	  	print(qplot(distance,-depth,data=spData[subSample,],color=spData[subSample,1]))
		print(qplot(distance,-depth,data=spData[subSample,],color=tpsModel$residuals)+scale_colour_gradient2(low="red",high="blue",mid="white",midpoint=0,name = "Residual"))
	  	print(qplot(distance,-depth,data=spData[subSample,],color=tpsModel$residuals/spData[subSample,1])+scale_colour_gradient2(low="red",high="blue",mid="white",midpoint=0,name = "Residual Ratio"))
	  	dev.off()
  	}

  	return(pred)
}


interpolation_krig <- function(spData,grid,idw,...){
	require(gstat)
 	var_name <- names(spData)[1]
 	# print(spData$distance)
 	print(range(spData$distance))

 	scaleSetting <- scaleFactor(spData$distance,spData$depth)
    scaledCoord <- applyScale(grid$distance,grid$depth,scaleSetting)
    grid$scaled_x <- scaledCoord$sx
    grid$scaled_y <- scaledCoord$sy
    
    scaledCoord <- applyScale(spData$distance,spData$depth,scaleSetting)
    spData$scaled_distance <- scaledCoord$sx
    spData$scaled_depth <- scaledCoord$sy

    coordinates(spData)=~scaled_distance+scaled_depth
    # coordinates(grid)=~scaled_x+scaled_y
 	
 	spData <- remove.duplicates(spData,zero=0.005,remove.second=TRUE) 	# remove too near points
 	# spData=remove.duplicates(spData,zero=0.01,remove.second=TRUE) 	# remove too near points
 	
 	var_formu=as.formula(paste("res","~1"))
 	# print(list(...)$maxdist)
 	if(idw==TRUE){
 	    model_gstat <- gstat(NULL,id=var_name,formula=var_formu,data=spData,maxdist=list(...)$maxdist,nmax = 100,nmin = 20)
 	}else{
 	    model_gstat <- gstat(NULL,id=var_name,formula=var_formu,data=spData,maxdist=list(...)$maxdist,nmax = 100,nmin = 20)

		vgmFitting <- variogram_fitting(model_gstat,list(...)$plotid)

		# vgmFitting[[2]] is the scaling factor
		model_gstat$data[[1]]$data@coords[,2] <- model_gstat$data[[1]]$data@coords[,2]*vgmFitting[[2]]
		model_gstat<-gstat(model_gstat,id=var_name,model=vgmFitting[[1]])
		grid$scaled_y <- grid$scaled_y*vgmFitting[[2]]
	}
 	 coordinates(grid)=~scaled_x+scaled_y

 	 pred<-predict(model_gstat,grid,debug.level=-1)
 	 print(summary(spData$res))
 	 print(summary(pred[[1]]))
 	 return(pred[[1]])
}


variogram_fitting <- function(g,plotid){
	var_name <- names(g$data)[[1]]
	localRange <- g$data[[1]]$maxdist
	print(paste("localRange:",localRange))
	g0 <- g
	Y0 <- g0$data[[1]]$data@coords[,2]
	kgmodel <- config$model

	optimFunc <- function(K,g){
		# function to do optimization
		g$data[[1]]$data@coords[,2] <- Y0*K

		# variogram	
		v <- variogram(g,cressie=T,alpha=c(0,90),cutoff=localRange,width=localRange/15)
		v_0 <- subset(v,dir.hor == 0 & np>10)  # 0 direction 
		v_90 <- subset(v,dir.hor == 90 & np>10) # 90 direction
		
		if(nrow(v_0)<1 | nrow(v_90)<1)
			return(Inf)

		v_0_model <- fit.variogram(v_0,vgm(NA,kgmodel,NA,v_0$gamma[1]),fit.method=7) # fit the 0 direction
		v_90_model <- fit.variogram(v_90,vgm(NA,kgmodel,NA,v_0$gamma[1]),fit.method=7)  # fit the 90 direction

		if( (v_0_model$range[2])<0 | (v_90_model$range[2])<0){
			print("bad WLS fit, fix nugget = 0 to fit")
			v_0_model <- fit.variogram(v_0,vgm(NA,kgmodel,NA,0),fit.method=7,fit.sill = c(FALSE,TRUE))
			v_90_model <- fit.variogram(v_90,vgm(NA,kgmodel,NA,0),fit.method=7, fit.sill = c(FALSE,TRUE))
			if((v_0_model$range[2])<0 | (v_90_model$range[2])<0){
				v_0_model <- fit.variogram(v_0,vgm(NA,kgmodel,NA,0),fit.method=7,fit.sill = c(FALSE,FALSE))
				v_90_model <- fit.variogram(v_90,vgm(NA,kgmodel,NA,0),fit.method=7, fit.sill = c(FALSE,FALSE))
			}
		}

		# RMSE of the variogram models
		dist_vector <- unique(v$dist)
		a <- mean((variogramLine(v_0_model, dist_vector = dist_vector)$gamma-variogramLine(v_90_model, dist_vector = dist_vector)$gamma)^2)
		return(a)
	}

	### change optimK tol from 0.1 to 0.05
	optimK <- optimize(optimFunc,config$K,g=g0,tol = 0.05)$minimum[1]
	
	print(paste("bestK:",optimK))
	print("optim Finish")

	g0$data[[1]]$data@coords[,2] <- Y0*optimK
	
	v <- variogram(g0,cressie=T,cutoff=localRange,width=localRange/15)
	v_model<- fit.variogram(v,vgm(NA,kgmodel,NA,v$gamma[1]),fit.method=7)

	if(v_model$range[2]<0){
		# if the range is not positive, set to positive
		v_model <- fit.variogram(v,vgm(NA,kgmodel,NA,0),fit.method=7, fit.sill = c(FALSE,TRUE))
		if(v_model$range[2]<0){
			v_model <- fit.variogram(v,vgm(NA,kgmodel,NA,0),fit.method=7, fit.sill = c(FALSE,FALSE))
		}
		# v_model$range[2] <- 0.02
	}

	print(v_model)
	
	if(config$meta){
		dir.create(file.path(config$outputFolder, "meta/variogram"), showWarnings = FALSE)
		bestModel <- kgmodel
		v_4direction <- variogram(g0,cressie=T,alpha=c(0,45,90,135),cutoff=localRange,width=localRange/15)
		v_2direction <- variogram(g0,cressie=T,alpha=c(0,90),cutoff=localRange,width=localRange/15)

		v_0 <- subset(v_2direction,dir.hor == 0)
		v_90 <- subset(v_2direction,dir.hor == 90)

		pdf(file=paste(config$outputFolder,"/meta/variogram/",plotid,"_krig_meta.pdf",sep=""))
	   	par(mfrow=c(2,2))
	   	print(plot(v_4direction,v_model,main=paste("4 Direction Variogram in Kriging Range"),pl=T))
	   	print(plot(v_2direction,v_model,main=paste("2 Direction Variogram in Kriging Range"),pl=T))
	    print(plot(v,v_model,main=paste("omnidirection Variogram in kriging range, model=",bestModel,"range = ",localRange),pl=T))
	    print(qplot(g0$data[[1]]$data@coords[,1],-g0$data[[1]]$data@coords[,2],colour=g$data[[1]]$data$res)+scale_colour_gradient2(low="red",high="blue",mid="white",midpoint=0,name="Residuals")+xlab("Adjusted Distance")+ylab(paste("Adjusted depth","Ratio:",optimK))+coord_fixed())
	    print(qplot(g0$data[[1]]$data@coords[,1],-g0$data[[1]]$data@coords[,2],colour=g$data[[1]]$data$res)+scale_colour_gradient2(low="red",high="blue",mid="white",midpoint=0,name="Residuals")+xlab("Adjusted Distance")+ylab(paste("Adjusted depth","Ratio:",optimK)))
	    dev.off()
	}
	
	return(list(v_model,optimK))
}


detrending <- function(dataSet,grid,method){
  require(fields,quietly = TRUE)
  # print(head(dataSet))
  var_name <- names(dataSet)[1]
  num <- nrow(grid)
  df <- -1
  
  if(method=="none"){
      return(list(trendSurface=rep(0,num),res=dataSet[,1],df = df))
  }

  if(method=="tps"){
      index <- seq(from = 1, to = nrow(dataSet), by = 2)
      # df <- max(5,as.integer(length(index)/200))
      df <- config$tpsDf
      trend_model <- Tps(dataSet[index,c("distance","depth")],dataSet[index,1],df = df) #df=config$tpsDf
      pred_orig <- predict(trend_model,dataSet[,c("distance","depth")])
  }
  else if(method=="loess"){
      var_formu <- as.formula(paste(var_name,"~distance+depth"))
      trend_model <- loess(var_formu,data=dataSet,span=0.25)
      pred_orig <- predict(trend_model,dataSet[,2:3])
  }else if(method == "linear"){
  	  var_formu <- as.formula(paste(var_name,"~I(distance)+I(depth)+I(distance*depth)+I(distance^2)+I(depth^2)"))
  	  trend_model <- lm(var_formu, data = dataSet)
  	  pred_orig <- predict(trend_model,dataSet[,2:3])
  }

  trendSurface <- predict(trend_model,grid[,c("distance","depth")])
  res <- dataSet[,1]-pred_orig
  return(list(trendSurface=c(trendSurface),res=res,df = df))
}




# interpolation_idw <- function(spData,grid){
#   # grid is the grid
#   # spData is the raw data, with first variable as the interest
#   require(gstat)
#   scaleSetting <- scaleFactor(spData$distance,spData$depth)
#   scaledCoord <- applyScale(grid$distance,grid$depth,scaleSetting)
#   grid$scaled_x <- scaledCoord$sx
#   grid$scaled_y <- scaledCoord$sy
  
#   scaledCoord <- applyScale(spData$distance,spData$depth,scaleSetting)
#   spData$scaled_distance0 <- scaledCoord$sx
#   spData$scaled_depth <- scaledCoord$sy

#   coordinates(spData)=~scaled_distance+scaled_depth
#   coordinates(grid)=~scaled_x+scaled_y

#   var_name <- names(spData)[1]
#   formu <- as.formula(paste(var_name,"~1",sep=""))
#   pred <- idw(formu, spData, grid)

# }

