#### This is the main part for krig Interpolation, with following functions
## 1: kriging Interpolation
## 2: GA variogram fit
require("GA")
interpolation_main <- function(Triaxus,int_method,det_method=NULL){
	Triaxus@grid <- createGrid(
		x=Triaxus@cleanData$distance,
		y=Triaxus@cleanData$depth,
		dx=Triaxus@config$gridSize[1],
		dy=Triaxus@config$gridSize[2],
		longitude=Triaxus@cleanData$longitude,
		latitude=Triaxus@cleanData$latitude)

	Triaxus@resultData <- Triaxus@grid

	# maxdist <- ifelse(config$krigingRange/Triaxus@numCycle>0.15,config$krigingRange/Triaxus@numCycle,0.15)
	maxdist <- 0.33
	# maxdist <- NULL

	if(Triaxus@seperate == FALSE){
		print("no seperating")
		for(var in Triaxus@config$interestVar){
			plotid <- paste(Triaxus@pathName,var,sep="_")
			outlierIndex <- init_filter(Triaxus@cleanData[,var],var)
			pred <- interpolation_sub(Triaxus@cleanData[outlierIndex<1,c(var,"distance","depth")],Triaxus@grid,int_method,det_method,plotid = paste(plotid,"_whole",sep=""),maxdist = maxdist)
			Triaxus@resultData[,var] <- ifelse(pred>0,pred,0)
		}
	}else{
		for(var in Triaxus@config$interestVar){
			print(var)
			plotid <- paste(Triaxus@pathName,var,sep="_")
	  		outlierIndex <- init_filter(Triaxus@cleanData[,var],var)
	    	pred_down <- interpolation_sub(subset(Triaxus@cleanData[outlierIndex<1,],direction==-1)[,c(var,"distance","depth")],Triaxus@grid,int_method,det_method,plotid = paste(plotid,"_down",sep=""),maxdist = maxdist)
	    	pred_up <- interpolation_sub(subset(Triaxus@cleanData[outlierIndex<1,],direction==1)[,c(var,"distance","depth")],Triaxus@grid,int_method,det_method,plotid = paste(plotid,"_up",sep=""),maxdist = maxdist)
	    	
	    	# pred <- apply(c(1:length(pred_down)),
	    	# 	function(x){
	    	# 		if(is.na(pred_down[i]))
	    	# 			return(pred_up*2)

	    	# 	)

	    	finalPrediction <- (pred_down+pred_up)/2.0
	    	Triaxus@resultData[,var] <- ifelse(finalPrediction>0,finalPrediction,0)
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

  	plotid <- list(...)$plotid
  	pdf(file=paste("~/Developer/Triaxus/output/meta/",plotid,"_tpsInt.pdf",sep=""))
  	# qplot(distance,depth,data=spData[subSample,],color=pred)
  	print(qplot(distance,depth,data=spData[subSample,],color=spData[subSample,1]))
	print(qplot(distance,depth,data=spData[subSample,],color=tpsModel$residuals)+scale_colour_gradient2(low="red",high="blue",mid="white",midpoint=0))
  	print(qplot(distance,depth,data=spData[subSample,],color=tpsModel$residuals/spData[subSample,1])+scale_colour_gradient2(low="red",high="blue",mid="white",midpoint=0))
  	dev.off()
  	return(pred)
}


interpolation_krig <- function(spData,grid,idw,...){
	require(gstat)
 	var_name <- names(spData)[1]
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
 	
 	spData=remove.duplicates(spData,zero=0.005,remove.second=TRUE) 	# remove too near points
 	var_formu=as.formula(paste("res","~1"))
 	# print(list(...)$maxdist)
 	if(idw==TRUE){
 	    model_gstat=gstat(NULL,id=var_name,formula=var_formu,data=spData,maxdist=list(...)$maxdist,nmax = 100,nmin = 20)
 	}else{
 	    model_gstat=gstat(NULL,id=var_name,formula=var_formu,data=spData,maxdist=list(...)$maxdist,nmax = 100, nmin = 20)

		vgmFitting=variogram_fitting(model_gstat,list(...)$plotid)
		model_gstat$data[[1]]$data@coords[,2] <- model_gstat$data[[1]]$data@coords[,2]*vgmFitting[[2]]
		grid$scaled_y <- grid$scaled_y*vgmFitting[[2]]
		model_gstat<-gstat(model_gstat,id=var_name,model=vgmFitting[[1]])
	}
 	 coordinates(grid)=~scaled_x+scaled_y
 	 # print(grid)
 	 pred<-predict(model_gstat,grid,debug.level=-1)	
 	 return(pred[[1]])
}


variogram_fitting <- function(g,plotid){
	var_name <- names(g$data)[[1]]
	localRange <- g$data[[1]]$maxdist
	print(paste("localRange:",localRange))
	g0 <- g
	Y0=g0$data[[1]]$data@coords[,2];
	
	optimFunc <- function(K,g){
		g$data[[1]]$data@coords[,2] <- Y0*K
		v <- variogram(g,cressie=T,alpha=c(0,90),cutoff=localRange,width=localRange/15)
		v_0 <- subset(v,dir.hor == 0 & np>10)
		v_90 <- subset(v,dir.hor == 90 & np>10)
		# v_0 <- subset(v,dir.hor == 0)
		# v_90 <- subset(v,dir.hor == 90)

		if(nrow(v_0)<1 | nrow(v_90)<1)
			return(Inf)
		# print(v)
		# v_0_model <- fit.variogram(v_0,vgm(NA,"Exp",NA,v_0$gamma[1]),fit.method=7)
		v_0_model <- fit.variogram(v_0,vgm(NA,"Gau",NA,v_0$gamma[1]),fit.method=7)

		# v_90_model <- fit.variogram(v_90,vgm(NA,"Exp",NA,v_0$gamma[1]),fit.method=7)
		v_90_model <- fit.variogram(v_90,vgm(NA,"Gau",NA,v_0$gamma[1]),fit.method=7)
		# if(v_0_model$range[2]<0 | v_90_model$range[2]<0)
			# return(Inf)
		# v_model <- fit.variogram(v_90,vgm(NA,"Gau",NA,v$gamma[1]),fit.method=7)
		# a <- max((variogramLine(v_model,dist_vector = v_0$dist)$gamma - v_0$gamma)^2)+max((variogramLine(v_model,dist_vector = v_90$dist)$gamma - v_90$gamma)^2)
		
		if( (v_0_model$range[2])<0 | (v_90_model$range[2])<0){
			print("bad WLS fit, switch to OLS")
			# print(K)
			# print(v_0)
			# print(v_90)
			v_0_model <- fit.variogram(v_0,vgm(NA,"Gau",NA,v_0$gamma[1]),fit.method=0)
			v_90_model <- fit.variogram(v_90,vgm(NA,"Gau",NA,v_0$gamma[1]),fit.method=0)
		}

		a <- mean((variogramLine(v_0_model, dist_vector = v_0$dist)$gamma-variogramLine(v_90_model, dist_vector = v_0$dist)$gamma)^2)
		# print(a)
		# v_0_model_gau2 <- fit.variogram(v_0,vgm(max(v$gamma)*2,"Gau",localRange*2,v_0$gamma[1]),fit.method=7)
		# v_90_model_gau2 <- fit.variogram(v_90,vgm(max(v$gamma)*2,"Gau",localRange*2,v_0$gamma[1]),fit.method=7)


		# a <- abs(sum(v_0_model$psill)/sum(v_90_model$psill)-1)+abs(v_0_model$range[2]/v_90_model$range[2]-1)+abs(v_0_model$psill[1]/v_90_model$psill[1]-1)
		# b = abs(sum(v_0_model_exp$psill)/sum(v_90_model_exp$psill)-1)+abs(v_0_model_exp$range[2]/v_90_model_exp$range[2]-1)+abs(v_0_model_exp$psill[1]/v_90_model_exp$psill[1]-1)
		# a2 <- abs(sum(v_0_model_gau2$psill)/sum(v_90_model_gau2$psill)-1)+abs(v_0_model_gau2$range[2]/v_90_model_gau2$range[2]-1)+abs(v_0_model_gau2$psill[1]/v_90_model_gau2$psill[1]-1)
		# return(min(a,a2))
		# a <- a*(1+(K-1)/4)
		return(a)
	}

	# miniError=c()
	# for(K in seq(1,5,0.2)){
		# print(K)
		# miniError=c(miniError,optimFunc(K,g))
	# }
	# optimK <- seq(0.2,5,0.2)[which.min(miniError)]
	
	# optimK <- optim(1, optimFunc,gr=NULL,g0,
 #      method = c("Nelder-Mead"),
 #      lower = 0.2, upper = 10,
 #      control = list(ndeps=0.5), 
 #      hessian = FALSE)$par[1]
	# for(i in ){
		# g$data[[1]]$data@coords[,2] <- Y0*K
		# v <- variogram(g,cressie=T,alpha=c(0,90),cutoff=localRange,width=localRange/15)
		# v_0 <- subset(v,dir.hor == 0)
		# v_90 <- subset(v,dir.hor == 90)
		# if(v_0)
	# }

	optimK <- optimize(optimFunc,config$K,g=g0,tol = 0.1)$minimum[1]
	# GA works also good maybe
	# optimK <- ga(type = "real-valued",fitness = function(x,g) -optimFunc(x,g0), g=g,min = c(1),max = c(5),popSize = 20,maxiter = 10)@solution[1]
	# K_optins <- seq(1,6,0.2)
	# optimK <- K_optins[which.min(sapply(K_optins,optimFunc,g=g0))]


	print(paste("bestK:",optimK))
	# print(optimK_normal)
	print("optim Finish")
	g0$data[[1]]$data@coords[,2] <- Y0*optimK
	
	v <- variogram(g0,cressie=T,cutoff=localRange,width=localRange/15)
	# v_model<- fit.variogram(v,vgm(NA,"Exp",NA,v$gamma[1]),fit.method=7)
	v_model<- fit.variogram(v,vgm(NA,"Gau",NA,v$gamma[1]),fit.method=7)


	if(v_model$range[2]<0){
		v_model$range[2] <- 0.02
	}
	# v_model_gau2 <- fit.variogram(v,vgm(max(v$gamma)*2,"Gau",localRange,v$gamma[1]),fit.method=7)
	
  
  	# v_model_gau_err <- attr(v_model_gau,"SSErr")
  	# v_model_exp_err <- attr(v_model_exp,"SSErr")
  	# v_model_gau_err2 <- attr(v_model_gau2,"SSErr")

  	v_4direction <- variogram(g0,cressie=T,alpha=c(0,45,90,135),cutoff=localRange,width=localRange/15)
	v_2direction <- variogram(g0,cressie=T,alpha=c(0,90),cutoff=localRange,width=localRange/15)

	
	v_0 <- subset(v_2direction,dir.hor == 0)
	v_90 <- subset(v_2direction,dir.hor == 90)

	# print(v_model_gau_err)
	# print(v_model_exp)
	# if(v_model_gau_err>v_model_exp_err){
	# 	v_model <- v_model_exp
	# }else{
	# 	v_model <- v_model_gau
	# }
	# # v_model <- ifelse(v_model_gau_err>v_model_exp_err,v_model_exp,v_model_gau) # This will raise warnings, don't write like this
	# bestModel <- ifelse(v_model_gau_err>v_model_gau_err2,"Exp","Gau")

	# v_model <- v_model_gau
	bestModel <- "Gau"

	#pdf(file=paste("~/Developer/Triaxus/output/variogram/",plotid,"_krig_meta.pdf",sep=""))
  #  par(mfrow=c(2,2))
  #  print(plot(v_4direction,v_model,main=paste("4 Direction Variogram in Kriging Range"),pl=T))
  #  print(plot(v_2direction,v_model,main=paste("2 Direction Variogram in Kriging Range gau"),pl=T))
    # print(plot(v_2direction,v_model_gau2,main=paste("2 Direction Variogram in Kriging Range gau2"),pl=T))

  #  print(plot(v,v_model,main=paste("omnidirection Variogram in kriging range, model=",bestModel),pl=T))

    # print(plot(v_0,v_0_model_sph,main=paste("horizon fit individually"),pl=T))
    # print(plot(v_90,v_90_model_sph,main=paste("vertical fit individually"),pl=T))

  #  print(qplot(g0$data[[1]]$data@coords[,1],-g0$data[[1]]$data@coords[,2],colour=g$data[[1]]$data$res)+scale_colour_gradient2(low="red",high="blue",mid="white",midpoint=0,name="Residuals")+xlab("Adjusted Distance")+ylab(paste("Adjusted depth","Ratio:",optimK))+coord_fixed())
  #	dev.off()
  # print(optimK)
	return(list(v_model,optimK))
}


detrending <- function(dataSet,grid,method){
  require(fields,,quietly = TRUE)
  # print(head(dataSet))
  var_name=names(dataSet)[1]
  num=nrow(grid)
  # print(paste("Detrending on",var_name," ,using",type))

  if(method=="none"){
      return(list(trendSurface=rep(0,num),res=dataSet[,1]))
  }

  if(method=="tps"){
      index=seq(from = 1, to = nrow(dataSet), by = 2)
      trend_model=Tps(dataSet[index,c("distance","depth")],dataSet[index,1],df=config$tpsDf)
      pred_orig=predict(trend_model,dataSet[,c("distance","depth")])
  }
  else if(method=="linear"){
      var_formu=as.formula(paste(var_name,"~distance+depth"))
      trend_model=loess(var_formu,data=dataSet,span=0.25)
      pred_orig=predict(trend_model,dataSet[,2:3])
  }

  trendSurface=predict(trend_model,grid[,c("distance","depth")])
  # print(trendSurface)
  res=dataSet[,1]-pred_orig
  return(list(trendSurface=c(trendSurface),res=res))
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

