rm(list=ls())
source("./src/config.R")
source("./src/Triaxus_class.R")
source("./src/preprocessing.R")
source("./src/interpolation.R")
source("./src/hotspot.R")
source("./src/clustering.R")
source("./src/plot.R")
source("./src/misc.R")
library(RColorBrewer)
library(leaflet)
library(ggplot2)
options(rgl.useNULL=FALSE)
library(rgl)

outputFolder <- "../output/"


main<- function(newfile,oldResult,seabirdIndex,bbeIndex,seabird_cutoff,newResult=FALSE){
	myTriaxus <- new("Base_Triaxus",config,newfile,seabirdIndex,bbeIndex,seabird_cutoff)
	print(newfile)
	myTriaxus <- preprocessing(myTriaxus) %>% interpolation_main(int_method="krige",det_method="tps") %>% hotspot_main() 
  	
  	if(newResult | !file.exists(oldResult)) {
		print("creating new File")
		allTriaxus <- list()
		allTriaxus[[myTriaxus@pathName]] <- myTriaxus
		saveRDS(allTriaxus,oldResult)
	}
	else{
		allTriaxus <- readRDS(oldResult)
		allTriaxus[[myTriaxus@pathName]] <- myTriaxus
		saveRDS(allTriaxus,oldResult)
	}
	# 
}

calculate <- function(){
	oldResult <- "../output/Manitowoc.rds"
	main("../input/Manitowoc_01.csv",oldResult,409,371,c(850,-1),TRUE)
	main("../input/Manitowoc_02.csv",oldResult,3106,3067,c(3120,9352))
	main("../input/Manitowoc_03.csv",oldResult,6409,6370,c(3120,6360))
	main("../input/Manitowoc_04.csv",oldResult,244,207,c(1,3311))
	main("../input/Manitowoc_05.csv",oldResult,0,0,c(1930,8630))
	allTriaxus <- clustering_main(readRDS(oldResult),c("Seabird_temperature","Spec.Cond"),Ks=c(2,3,4,5,6))
	saveRDS(allTriaxus,"../output/Manitowoc_Results.rds")
	rm(allTriaxus)

	oldResult <- "../output/Muskegon.rds"
	main("../input/Muskegon River_01.csv",oldResult,7267,7228,c(1930,7200),TRUE)
	main("../input/Muskegon River_02.csv",oldResult,2017,1980,c(1880,-1))
	main("../input/Muskegon River_03a.csv",oldResult,120,81,c(150,6745))
	main("../input/Muskegon River_04a.csv",oldResult,377,338,c(1400,6614))
	allTriaxus <- clustering_main(readRDS(oldResult),c("Seabird_temperature","Spec.Cond"),Ks=c(2,3,4,5,6))
	saveRDS(allTriaxus,"../output/Muskegon_Results.rds")
	rm(allTriaxus)

	oldResult <- "../output/Pere.rds"
	main("../input/Pere Marquette River_01.csv",oldResult,5828,5789,c(200,3260,3780,7450),TRUE)
	main("../input/Pere Marquette River_02.csv",oldResult,513,473,c(350,4705))
	main("../input/Pere Marquette River_03.csv",oldResult,2138,2100,c(2200,8637))
	main("../input/Pere Marquette River_04.csv",oldResult,833,794,c(2000,8299))
	main("../input/Pere Marquette River_05.csv",oldResult,6806,6769,c(3030,7999))
	main("../input/Pere Marquette River_06.csv",oldResult,3654,3615,c(360,5100))
	allTriaxus <- clustering_main(readRDS(oldResult),c("Seabird_temperature","Spec.Cond"),Ks=c(2,3,4,5,6))
	saveRDS(allTriaxus,"../output/Pere_Results.rds")
	rm(allTriaxus)
}

paper_plot <- function(){
	# Manitowoc
	transmatrix=matrix(c(0.01585434,-0.68476230,-0.72859418, 0.00000000,0.4700032, 0.6482660,-0.5990396,0.0000000,0.8825220,-0.3329444,0.3321184,0,0,0,0,1),nrow=4,ncol=4)
	riverMouth <- list(c(-87.5600,44.1425),c(-87.6436,44.0900))
	#plot_3d_clustering(readRDS("../output/Manitowoc_Results.rds"),3,transmatrix=transmatrix,riverMouth=riverMouth)
	paperPlotSub("../output/Manitowoc_Results.rds",transmatrix,riverMouth,3)


	transmatrix=matrix(c(0.0670390948653221 , -0.821753919124603 , -0.565885365009308 , 0 , -0.435517132282257 , 0.486176699399948 , -0.757599711418152 , 0 , 0.897680580615997 , 0.297241151332855 , -0.325294703245163 , 0 , 0 , 0 , 0 , 1),nrow=4,ncol=4)
	riverMouth <- list(c(-86.33,43.22))
	#plot_3d_clustering(readRDS("../output/Muskegon_Results.rds"),3,transmatrix=transmatrix,riverMouth=riverMouth)
	paperPlotSub("../output/Muskegon_Results.rds",transmatrix,riverMouth,3)

	# Pere
	transmatrix=matrix(c(-0.02324756 , -0.68113685 , -0.731787025928497 , 0 , -0.583894312381744 , 0.603414714336395 , -0.543100714683533 , 0 , 0.811497092247009 , 0.41466024518013 , -0.411739528179169 , 0 , 0 , 0 , 0 , 1),nrow=4,ncol=4)
	riverMouth <- list(c(-86.46,43.95))
	#plot_3d_clustering(readRDS("../output/Pere_Results.rds"),2,transmatrix=transmatrix,riverMouth=riverMouth)
	paperPlotSub("../output/Pere_Results.rds",transmatrix,riverMouth,2)
}

paperPlotSub <- function(triaxusFile,transmatrix,riverMouth,K){
	locationName <- strsplit(basename(triaxusFile),"_")[[1]][1]
	print(locationName)
	allTriaxus <- readRDS(triaxusFile)
	for(var in config$interestVar){
		print(var)

		plot_3d_value(allTriaxus,var,transmatrix=transmatrix,riverMouth=riverMouth)
		rgl.snapshot(paste(outputFolder,locationName,var,"value.png",sep="_"))
    	rgl.close()

		# postscript()
		plot_3d_hotspot(allTriaxus,var,transmatrix=transmatrix,riverMouth=riverMouth)
		rgl.snapshot(paste(outputFolder,locationName,var,"hotspot.png",sep="_"))
    	rgl.close()
	}

	# print silhouetter width

	write.csv(x=as.data.frame(allTriaxus[[1]]@clusteringResults$silhouetteList),file=paste(outputFolder,locationName,"clusterSil.csv",sep="_"))
	# plot clustering analysis
	plot_3d_clustering(allTriaxus,K,transmatrix=transmatrix,riverMouth=riverMouth)
	rgl.snapshot(paste(outputFolder,locationName,"cluster.png",sep="_"))
    rgl.close()
	
	# print boxplot
	pdf(paste(outputFolder,locationName,K,"boxplot.pdf",sep="_"),height=2.75,width=3.25)
	plot_boxplot(allTriaxus,config$interestVar,K)
	dev.off()
}

# calculate()
paper_plot()
