##### Clustering model
require(cluster)

clustering <- function(allTriaxus,variableForClustering,Ks){
    dataSet <- data.frame()
    totalAvail <- c()

    for(Triaxus in allTriaxus){
      # available <- Triaxus@grid$available
      available <- !as.logical(apply(is.na(Triaxus@resultData[,variableForClustering]),1,sum))

      totalAvail <- c(totalAvail,available)
      dataSet <- rbind(dataSet,Triaxus@resultData[available,variableForClustering])
    }

    clusteringIndexList <- list()
    avg_silWidthList <- list()

    for(K in Ks){
      clusterName <- paste("cluster_",K,sep="")
      clusterModel <- clustering_sub(dataSet,K)

      clusterIndex <- rep(NA,length(totalAvail))
      clusterIndex[totalAvail] <- clusterModel$clusterIndex

      avg_silWidthList[[clusterName]] <- clusterModel$avg_silWidth
      clusteringIndexList[[clusterName]] <- clusterIndex
    }
    
    return(list(clusteringIndexList = clusteringIndexList, avg_silWidthList = avg_silWidthList))
}


clustering_sub <- function(dataSet,K){
  # dataSet is the data frame with clustering variables as columns
  print(K)
  dataSet <- scale(dataSet)
  dataSet <- round(dataSet,2)
  kmcluster=kmeans(dataSet,centers=K,nstart = 30,iter.max=40)
  dissE <- daisy(dataSet)
  sk <- silhouette(kmcluster$cluster, dissE^2)
  avgSw <- summary(sk)$avg.width
  return(list(clusterIndex = kmcluster$cluster,avg_silWidth = avgSw))
}

clustering_main <- function(allTriaxus,variableForClustering,Ks=c(2)){
    clusteringResult <- clustering(allTriaxus,variableForClustering,Ks)
    endingPoint <- 0


    for(i in 1:length(allTriaxus)){
      
      startingPoint <- endingPoint+1

      sampleNum <- nrow(allTriaxus[[i]]@resultData)

      #available <- !as.logical(apply(is.na(allTriaxus[[i]]@resultData[,variableForClustering]),1,sum))
      
      clusterIndexList <- list()
      silhouetteList <- list()
          
      endingPoint <- startingPoint+sampleNum-1
      for(j in 1:length(clusteringResult$avg_silWidthList)){

          clusterName <- names(clusteringResult$avg_silWidthList)[j]
          
          clusterIndexList[[clusterName]] <- clusteringResult$clusteringIndexList[[j]][startingPoint:endingPoint]
          silhouetteList[[clusterName]] <- clusteringResult$avg_silWidthList[[j]]

      }

      allTriaxus[[i]]@clusteringResults$variable <- variableForClustering
      allTriaxus[[i]]@clusteringResults$clusteringIndex <- clusterIndexList
      allTriaxus[[i]]@clusteringResults$silhouetteList <- silhouetteList

    }
      

    return(allTriaxus)
}


