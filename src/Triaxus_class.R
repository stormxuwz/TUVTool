# source("./src/rawFileParser.R")
require(methods)
setClass(
	"Base_Triaxus",
	slots = c(
		pathName="character",
		fileName = "character",
		config = "list",
		rawData = "data.frame",
		cleanData = "data.frame",
		grid = "data.frame",
		resultData = "data.frame",
		hotspotData = "data.frame",
		numCycle = "numeric",
		Seabird_cutoff="numeric",
		BBE_cutoff="numeric",
		clusteringResults = "list",
		SeabirdAnchorIndex="numeric",
		BBEAnchorIndex="numeric",
		separate = "logical",
		separationNode = "numeric"
	)
)

setMethod("initialize","Base_Triaxus", 
	function(.Object,config,fileName,SeabirdAnchorIndex=0,BBEAnchorIndex=SeabirdAnchorIndex,Seabird_cutoff,realName=NULL,separate= TRUE,rawData = NULL){
		
		.Object@fileName <- fileName
		.Object@config <- config
		
		if(is.null(realName)){
			.Object@pathName <- substr(basename(fileName), 1, nchar(basename(fileName)) - 4)
			FileExtension <- tools::file_ext(fileName)
		}
		else{
			.Object@pathName <- substr(basename(realName), 1, nchar(basename(realName)) - 4)
			FileExtension <- tools::file_ext(realName)
		}
		
		if(is.null(rawData)){
			if(!is.null(fileName)){
				if(FileExtension=="csv"){
					.Object@rawData <- read.csv(fileName)
				}else if(FileExtension=="dat"){
					.Object@rawData <- readingRawFile(fileName)
				}else{
					return(NULL)
				}
			}
		}else{
			.Object@rawData <- rawData
		}

		.Object@Seabird_cutoff <- Seabird_cutoff
		.Object@SeabirdAnchorIndex <- SeabirdAnchorIndex
		.Object@BBEAnchorIndex <- BBEAnchorIndex
		.Object@separate <- separate
		# .Object@BBE_cutoff <- BBE_cutoff
		if(length(Seabird_cutoff)%%2==0){
			cutoff <- c()
			for(i in seq(1,length(Seabird_cutoff),by=2)){
				if(Seabird_cutoff[i]<1){
					Seabird_cutoff[i]  <- 1
				}

				if(Seabird_cutoff[i]>nrow(.Object@rawData)){
					Seabird_cutoff[i]  <- 1
				}

				if(Seabird_cutoff[i+1]<1){
					Seabird_cutoff[i+1] <- nrow(.Object@rawData)
				}

				if(Seabird_cutoff[i+1]>nrow(.Object@rawData)){
					Seabird_cutoff[i+1] <- nrow(.Object@rawData)
				}
				cutoff <- c(cutoff,Seabird_cutoff[i]:Seabird_cutoff[i+1])
			}
			.Object@Seabird_cutoff <- cutoff
			# print(cutoff)
		}
		return(.Object)
	}
)

