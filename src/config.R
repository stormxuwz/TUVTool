require("leaflet")
config <- list()

#####
# install.packages(c("leaflet","ggplot2","sp","gstat","RColorBrewer","rgl","dplyr","ggmap","lubridate","dismo","rglwidget","shiny","reshape2","dygraphs","pracma"))
#####


###########################
### Data File settings ####
###########################

config$BBE_name <- c("total","green","bluegreen","diatom")
config$Seabird_name <- c("Seabird_temperature","DO","DOsat","conductivity","BAT","Spec.Cond")
config$LOPC_name <- c("Zdens","Zug", "Zug_small", "Zug_medium", "Zug_large","Zdens_small","Zdens_medium","Zdens_large")
config$varUnit <- c(Seabird_temperature="Temperature (C)",
	DO = "DO (mg/L)",
	Spec.Cond = "Spec Cond (uS/cm)",
	DOsat = "DO Sat (%)",
	total = "Chl (ug/L)",
	Zdens = "Zooplankton Density (1/m^3)",
	Zug = "Zooplankton Biomass (ug/m^3)",
	BAT = "Beam Attenuation Coeff (1/m)",
	density = "Density (kg/m^3)",
	Zdens_small = "Small Zooplankton Density (1/m^3)",
	Zug_small = "Small Zooplankton Biomass (ug/m^3)",
	Zdens_medium = "Medium Zooplankton Density (1/m^3)",
	Zug_medium = "Medium Zooplankton Biomass (ug/m^3)",
	Zdens_large = "Large Zooplankton Density (1/m^3)",
	Zug_large = "Large Zooplankton Biomass (ug/m^3)"
)
config$factorColor <- colorFactor(c("blue4","white","red","blue","yellow","green","aquamarine","darkorange3","darkorchid4","lightpink1"),c(-1:8))

######################
### User settings ####
######################
# the hotspot neighbor size
config$depth_distance_ratio <- 1
config$nbrange <- 0.75

config$ZugBins = list(all=10:128,small=10:50, medium=51:100, large=101:128)

# Interpolation 
config$maxdist <- 0.33
config$separate <- TRUE
config$tpsDf <- 10  # tps detrending results
config$K <- c(0.3,5)  # y axis scale factor range after scaled to [0,1]
config$model <- "Sph"
config$gridSize <- c(dx=0.2,dy=0.25) # grid size dx in KM unit and dy in m unit

config$interestVar <- c("Seabird_temperature","Spec.Cond","DO","DOsat","total","BAT","Zdens","Zug") # config$interestVar <- c("Zdens") for test
# config$interestVar <- c("BAT") # for test only

# output folder
config$outputFolder <- "~/Developer/Triaxus/output/"

#plot meta 
config$meta <- FALSE
