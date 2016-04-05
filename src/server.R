rm(list=ls())
library(shiny)
library(RMySQL)
library(dygraphs)
library(zoo)
library(reshape2)
library(leaflet)
library(RColorBrewer)
library(rgl)
library(rglwidget)
library(ggplot2)
source("plot.R")
library(plotly)
library(dplyr)
source("config.R")
source("Triaxus_class.R")
source("preprocessing.R")
source("interpolation.R")
source("hotspot.R")
source("clustering.R")
source("misc.R")


# color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
# 	scale = (length(lut)-1)/(max-min)
#     dev.new(width=1.75, height=5)
#     plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
#     axis(2, ticks, las=1)
# for (i in 1:(length(lut)-1)) {
# y = (i-1)/scale + min
#     	rect(0,y,10,y+1/scale, col=lut[i], border=NA)
#     }	
# }



# colorBarPlot <- function(value,color,title=""){
# 	# value and color are already sorted from min to max and is 31 levels by default. nticks = 5

# 	scale <- (max(value)-min(value))/10

# 	plot(c(0,10), c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)

# 	axis(2, (seq(min(value),max(value),length.out=5)-min(value))/scale, las=1)
	
# 	for (i in 1:10) {
# 		y = (i-1)/scale + min
#     	rect(0,y,10,y+1/scale, col=lut[i], border=NA)
#     }	
# }


# x <- seq(-10, 10, length = 30)
# y <- x
# f <- function(x, y) { r <- sqrt(x^2 + y^2); 10 * sin(r)/r }
# z <- outer(x, y, f)
# z[is.na(z)] <- 1

options(rgl.useNULL=TRUE)
options(shiny.maxRequestSize=30*1024^2) 



shinyServer(function(input,output,session)
{
	output$threeDMap <- renderRglwidget({
		# progress <- shiny::Progress$new()
    	# Make sure it closes when we exit this reactive, even if there's an error
    	# on.exit(progress$close())
    	# progress$set(message = "Making plot", value = 0)

		allTriaxus <- visResultData()
		if(is.null(allTriaxus)){
			return()
		}
		if(input$hotspot3d){
			plot_3d_hotspot(allTriaxus,input$varToVis)
		}else{
			plot_3d_value(allTriaxus,input$varToVis)
		}
		print("finished plotting")
		rglwidget()   # This step is slow...
	})
	
	output$colorBar <- renderPlot({
		allTriaxus <- visResultData()
		if(is.null(allTriaxus)){
			return(NULL)
		}
		if(input$hotspot3d){
			# coloPal <- colorFactor(palette = c("blue","white","Red"),as.factor(c(-1,0,1)))
			# colorBarValue <- as.factor(c(-1,0,1))
			# colorBarColor <- coloPal(colorBarValue)
			# colorBarPlot(colorBarValue,colorBarColor,nticks=3)
			return(NULL)
		}
		else{
			v <- c()
			for(myTriaxus in allTriaxus){
				v <- c(v,myTriaxus@resultData[,input$varToVis])
			}
			varRange <- range(v,na.rm=TRUE)
			coloPal <- colorNumeric(topo.colors(10),c(floor(varRange[1]),ceiling(varRange[2])))

			colorBarValue <- seq(floor(varRange[1]),ceiling(varRange[2]),length.out=100)
			colorBarColor <- coloPal(colorBarValue)
			colorBarPlot(colorBarValue,colorBarColor)
		}
	})

	output$dygraph <- renderDygraph({
		myTriaxus <- isolate(readTriaxus())
		dataPlot <- myTriaxus@rawData[,c("Seabird_depth","BBE_depth")]
		n <- nrow(dataPlot)
		dataPlot <- zoo(dataPlot,order.by=1:n)
		print(dataPlot)
	    dygraph(dataPlot, main = "Path segments")
  	})

	
	clusteredTriaxus <- reactive({
		clusterButton()
		allTriaxus <- isolate(visResultData())

		if(is.null(allTriaxus)){
			print("no data")
			return(NULL)
		}
		
		clusterName <- paste("cluster_",isolate(input$clusteringNum),sep="")

		if(is.null(allTriaxus[[1]]@clusteringResults$silhouetteList[[clusterName]])){
			progress <- shiny::Progress$new()
	    	# Make sure it closes when we exit this reactive, even if there's an error
	    	on.exit(progress$close())
	    	progress$set(message = "Start Clustering", value = 0)
			
			variableForClustering <- isolate(input$varForClustering)
			K <- isolate(as.integer(input$clusteringNum))
			allTriaxus <- clustering_main(allTriaxus,variableForClustering,Ks=c(K))
		}
		return(allTriaxus)
	})

	output$clustering <- renderRglwidget({
		allTriaxus <- clusteredTriaxus()
		if(is.null(allTriaxus)){
			return(NULL)
		}
		plot_3d_clustering(allTriaxus,isolate(input$clusteringNum))
		rglwidget()
	})

	output$boxplot <- renderPlot({
		allTriaxus <- clusteredTriaxus()
		if(is.null(allTriaxus)){
			return(NULL)
		}
		K <- isolate(as.integer(input$clusteringNum))
		print(K)
		plot_boxplot(allTriaxus,isolate(input$varForClustering),K)
	})


	# output$clustering <- renderRglwidget({
	# 	clusterButton()
	# 	allTriaxus <- isolate(visResultData())

	# 	if(is.null(allTriaxus)){
	# 		return()
	# 	}
		
	# 	clusterName <- paste("cluster_",isolate(input$clusteringNum),sep="")

	# 	if(is.null(allTriaxus[[1]]@clusteringResults$silhouetteList[[clusterName]])){
	# 		progress <- shiny::Progress$new()
	#     	# Make sure it closes when we exit this reactive, even if there's an error
	#     	on.exit(progress$close())
	#     	progress$set(message = "Start Clustering", value = 0)
			
	# 		variableForClustering <- isolate(input$varForClustering)
	# 		K <- isolate(as.integer(input$clusteringNum))
	# 		allTriaxus <- clustering_main(allTriaxus,variableForClustering,Ks=c(K))
	# 	}

	# 	plot_3d_clustering(allTriaxus,isolate(input$clusteringNum))
	# 	rglwidget()
	# })

	# observe({
	# 	clusterButton()
	# 	allTriaxus <- isolate(visResultData())

	# 	if(is.null(allTriaxus)){
	# 		return()
	# 	}
		
	# 	clusterName <- paste("cluster_",isolate(input$clusteringNum),sep="")

	# 	if(is.null(allTriaxus[[1]]@clusteringResults$silhouetteList[[clusterName]])){
	# 		progress <- shiny::Progress$new()
	#     	# Make sure it closes when we exit this reactive, even if there's an error
	#     	on.exit(progress$close())
	#     	progress$set(message = "Start Clustering", value = 0)
			
	# 		variableForClustering <- isolate(input$varForClustering)
	# 		K <- isolate(as.integer(input$clusteringNum))
	# 		allTriaxus <- clustering_main(allTriaxus,variableForClustering,Ks=c(K))
	# 	}

	# 	output$clustering <- renderRglwidget({
	# 		plot_3d_clustering(allTriaxus,isolate(input$clusteringNum))
	# 	})
	# 	output$boxplot <- renderPlot({
	# 		plot_boxplot(allTriaxus,input$variableForClustering,K)
	# 	})
	# })
	# geoData <<-NULL
	# latRange<<-NULL
	# longRange<<-NULL
	# rawData<<-NULL
	# allTriaxus <<-
	gobutton <- eventReactive(input$goButton, {})
	readButton <- eventReactive(input$readButton, {})
	clusterButton <- eventReactive(input$startClustering, {})

	visResultData <- reactive({
		if(is.null(input$resultFile$datapath[1])){
			return(NULL)
		}

		allTriaxus <- readRDS(input$resultFile$datapath[1])
		allNames <- names(allTriaxus)
		names(allNames) <- allNames
		varNames <- allTriaxus[[1]]@config$interestVar
		names(varNames) <- varNames
		updateSelectInput(session, "pathToVis",
			choices = as.list(allNames),selected=allNames[1]
		)

		updateSelectInput(session, "varToVis",
			choices = as.list(varNames),selected=varNames[1]
		)
		updateSelectizeInput(session, 'varForClustering', choices = config$interestVar, selected= NULL, server = FALSE)
		return(allTriaxus)
	})

	output$rawVis <- renderPlot({
		allTriaxus <- visResultData()
		if(is.null(allTriaxus)){
			return(NULL)
		}
		myTriaxus <- allTriaxus[[input$pathToVis]]

		if(is.null(myTriaxus)){
			return(NULL)
		}
		plot_raw(myTriaxus,input$varToVis)
	})

	
	
	output$rawExplor <- renderPlot({
		mydata <- rawData()
		mydata <- mydata[,c(input$calRawVar,"Seabird_depth")]
		mydata$id <- 1:nrow(mydata)
		names(mydata)[1]="value"
		print(names(mydata))
		qplot(id, -Seabird_depth, data=mydata,color=value)
	})


	output$twoD_Vis <- renderPlot({
		if(is.null(input$resultFile)){
			return()
		}

		allTriaxus <- visResultData()
		myTriaxus <- allTriaxus[[input$pathToVis]]

		if(input$hotspot2d){
			plot_hotspot(myTriaxus,input$varToVis)
		}else{
			plot_2d(myTriaxus,input$varToVis)
		}
		
	})


	readTriaxus <- reactive({
		realName <- input$newFile$name[1]
		progress <- Progress$new()
		myTriaxus <- new("Base_Triaxus",config,input$newFile$datapath[1],input$Seabird_anchor_index,
			input$BBE_anchor_index,c(input$startIndex,input$endIndex),realName)
		progress$set(message = "Preprocessing", value = 0)
		myTriaxus <- preprocessing(myTriaxus)
		on.exit(progress$close())
		print("preprocessing finished")
		return(myTriaxus)
	})

	readPreviousTriaxus <- reactive({
		if(is.null(input$PreviousFile$datapath[1])){
			return(NULL)
		}
		else{
			print("reading new file")
			return(readRDS(input$PreviousFile$datapath[1]))
		}
	})
	
	observe({
		readButton()
		myTriaxus <- isolate(readTriaxus())
		oldTriaxus <- isolate(readPreviousTriaxus())

		geoData <- myTriaxus@cleanData[,c("latitude","longitude","distance","UTC")]
		latRange <- range(geoData$latitude)
		longRange <- range(geoData$longitude)
		leafletProxy("calMap", data =geoData) %>% clearShapes() %>%addPolylines(lng=~longitude,lat=~latitude,color="red")  %>% fitBounds(longRange[1]-0.1, latRange[1]-0.1,longRange[2]+0.1, latRange[2]+0.1)

		if(!is.null(oldTriaxus)){
			print("oldTriaxus exist")	
			for(myTriaxus in oldTriaxus){
				previousGeoData <- myTriaxus@cleanData[,c("latitude","longitude","distance")]
				leafletProxy("calMap", data =previousGeoData) %>% addPolylines(lng=~longitude,lat=~latitude,color="yellow")
			}
		}
		updateTabsetPanel(session, "mapFilePage", selected = "Map")
	})


	readRawFile <- reactive({
		if(is.null(input$newFile$datapath[1])){
			return(NULL)
		}
		rawData <- read.csv(input$newFile$datapath[1])
		rawData$n <- 1:nrow(rawData)
		rawData <- rawData[,c("depth","depth.1","n")]
		names(rawData)[2] <- "depth_1"
		rawData		
	})


	output$rawData_check <- renderPlot({
		rawData <- readRawFile()
		if(!is.null(rawData)){
			startIndex <- ifelse(input$startIndex_check<1,1,input$startIndex_check)
			endIndex <- ifelse(input$endIndex_check<1,nrow(rawData),input$endIndex_check)

			subData <- rawData[startIndex:endIndex,]
			# print(summary(rawData))
			seabird_maxDepth_index <- which.max(subData$depth)
			BBE_maxDepth_index <- which.max(subData$depth_1)

			seabirdIndex <- subData$n[seabird_maxDepth_index]
			BBEIndex <- subData$n[BBE_maxDepth_index]

			plot(-depth~n,data=subData,type="l")
			lines(-depth_1~n,data=subData,col="blue")
			points(-subData$depth[seabird_maxDepth_index]~seabirdIndex,col="red")
			points(-subData$depth_1[BBE_maxDepth_index]~BBEIndex,col="red")
			
			# print(seabirdIndex)
			suggestion <- paste("Seabird Index:",seabirdIndex,"/FlorProb Index:",BBEIndex)
			updateTextInput(session,"suggestedAnchorIndex",value = suggestion)
		}
	})


	# Calculate 
	observe({
		gobutton()
		progress <- Progress$new()

		progress$set(message = "Reading file", value = 0)
		myTriaxus <- isolate(readTriaxus())

		progress$set(message = "Interpolation", value = 0)
		myTriaxus <- interpolation_main(myTriaxus,int_method=isolate(input$InterpolationMethod),det_method=isolate(input$detrendingMethod)) %>% hotspot_main()
		
		progress$set(message = "Reading Previous Results", value = 0)

		# if(is.null(isolate(input$PreviousFile$datapath[1]))){
		# 	allTriaxus=list()
		# }
		# else{
		# 	allTriaxus <- readRDS(isolate(input$PreviousFile$datapath[1]))
		# }
		allTriaxus <- isolate(readPreviousTriaxus())

		if(is.null(allTriaxus)){
			allTriaxus=list()
		}

		allTriaxus[[myTriaxus@pathName]] <- myTriaxus
		newTriaxusFileName <- paste("allTriaxus_upto_",myTriaxus@pathName,sep="")
		saveRDS(allTriaxus,paste("./testFile/",newTriaxusFileName,".rds",sep=""))
		on.exit(progress$close())

		
	})
	
	


	output$calRaw <- renderPlot({
		if(!is.null(rawData)){
			visData <- rawData[,c("Seabird_depth",input$calRawVar)]
			names(visData)[2]="value"
			ggplot(visData)+geom_point(aes(1:nrow(visData),-Seabird_depth,color=value))
		}
	})

	output$calMap <- renderLeaflet({		
       	leaflet("calMap") %>% clearShapes() %>% addTiles()
    })
	
	output$visMap <- renderLeaflet({
		print("visMap initialize")
       	leaflet("visMap") %>% clearShapes() %>% addTiles()
       	# print(geoData)
    })
})