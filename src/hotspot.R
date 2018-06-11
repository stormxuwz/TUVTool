######## This function is to give local G values to find the local hotspot and coldspot
hotspot <- function(mydata,nbrange=0.75,x_y_ratio=1){
  library(spdep,quietly = TRUE)
  library(gridExtra,,quietly = TRUE)

  GAvailIndex <- !is.na(mydata[,1])
  mydata <- mydata[GAvailIndex,]
  xycoords <- cbind(mydata$distance*x_y_ratio, mydata$depth)
  nb <- dnearneigh(xycoords,0,nbrange)
  G <- localG(as.numeric(mydata[,1]), nb2listw(nb, style="B"))
  G <- as.numeric(G)

  return(list(G=G,GAvailIndex=GAvailIndex))
}

hotspot_main <- function(Triaxus){
    availableIndex <- Triaxus@grid$available
    Triaxus@hotspotData <- Triaxus@grid

    for(var in Triaxus@config$interestVar){
        subData <- Triaxus@resultData[,c(var,"distance","depth")]
        
        hotspot_index <- rep(NA,nrow(subData))
        if(sum(!is.na(subData[,var]))>0){
          print("starting hotspot")
          hotspotModel <- hotspot(subData,nbrange=Triaxus@config$nbrange,x_y_ratio=Triaxus@config$depth_distance_ratio)
          hotspot_index[hotspotModel$GAvailIndex] <-  hotspotModel$G
        }
        Triaxus@hotspotData[,var] <- hotspot_index
    }

    print("hotspot analysis finished")
    return(Triaxus)
}

