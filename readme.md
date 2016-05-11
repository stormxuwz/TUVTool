# Towed Undulating Vehicles Analyzing Tool
	
This folder creates a tool to analyze Towed undulating data, especially for Triaxus TUV data.

It can 
* Perform automated kriging analysis
* Perform hotspot analysis and cluster analysis
* Visualize the results in 2D and 3D

---
To use:

	source("./src/config.R")  # configuration list
	source("./src/Triaxus_class.R") # Triaxus base class
	source("./src/preprocessing.R")
	source("./src/interpolation.R")
	source("./src/hotspot.R")
	source("./src/clustering.R")
	source("./src/plot.R")
	source("./src/misc.R")
	source("./src/rawFileParser.R")
	library(RColorBrewer)
	library(leaflet)
	library(ggplot2)
	options(rgl.useNULL=FALSE)
	library(rgl)
	myTriaxus <- new("Base_Triaxus",config,newfile,seabirdIndex,bbeIndex,seabird_cutoff,separate=config$separate,rawData = NULL,realName=NULL)
	
where 

*	config: list defined in config.R
*	newFile: newfile name, can take .dat LOPC rawfile or csv file (sample file can be downloaded here [link](https://www.dropbox.com/sh/jezabaryohpfdnf/AABCzXVb0AVOPhRmDVgX538-a?dl=0))
*	seabirdIndex, bbeIndex:  the index of seabird data and bbe data that represents the same valley points. Need to determine manually 
*	seabird_cutoff:  a vector contains the start and end index for the
*	separate: whether to separate upcast and downcast before analyzing
*	rawData: bind data to the class rather than reading the data. If set, no reading file is performed
*	realName: for web app only since the file uploaded may not preserve the original name
	
To interpolate:

	myTriaxus <- preprocessing(myTriaxus) %>% interpolation_main(int_method="krige",det_method="tps") %>% hotspot_main() 
	
"%>%" is the pipeline operation. The preprocessing, interpolation_main and hotspot_main will take Triaxus class as the first input.


	allTriaxus <- list(myTriaxus)
	
	# To clustering
	# variableForClustering: which variable considered in cluster analysis
	# Ks: how many clusters to find
	allTriaxus <- clustering_main(allTriaxus, variableForClustering= c("Seabird_temperature","Spec.Cond"),Ks=c(2,3,4,5,6))
		
	# To plot boxplot in each cluster
	plot_boxplot (allTriaxus,variables = c("Seabird_Temperature"))
	
	# To plot 3D values, 3D plot take allTriaxus (a list of Triaxus_base class) as input
	plot_3d_value(allTriaxus,var="Seabird_Temperature",)
	plot_3d_hotspot (allTriaxus,var="Seabird_Temperature")
	plot_3d_clustering (allTriaxus,K=2)

	# To plot 2D values, 2D plot take Triaxus_base class and variable name
	plot_2d(myTriaxus,var= "Seabird_Temperature")
	plot_hotspot (myTriaxus,var= "Seabird_Temperature")
	plot_raw (myTriaxus,var= "Seabird_Temperature")
