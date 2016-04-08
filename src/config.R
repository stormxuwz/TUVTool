require("leaflet")
config <- list()

config$depth_distance_ratio <- 1
config$nbrange <- 0.75
config$BBE_name <- c("total","green","bluegreen","diatom")
config$Seabird_name <- c("Seabird_temperature","DO","DOsat","cond","BAT","Spec.Cond")
config$gridSize=c(dx=0.2,dy=0.25)
config$interestVar <- c("Seabird_temperature","Spec.Cond","DO","DOsat","total","BAT","Zdens","Zug")
# config$interestVar <- c("Zdens","Seabird_temperature")
config$LOPC_name <- c("Zdens","Zug")
config$varUnit <- c(
	Seabird_temperature="Temperature (C)",
	DO = "DO (mg/L)",
	Spec.Cond = "Spec Cond (uS/cm)",
	DOsat = "DO Sat (%)",
	total = "Chl (ug/L)",
	Zdens = "Zooplankton Density (1/m^3)",
	Zug = "Zooplankton Biomass (ug/m^3)",
	BAT = "Beam Attenuation Coeff (1/m)"
)

config$factorColor <- colorFactor(c("blue4","white","red","blue","yellow","green","aquamarine","darkorange3","darkorchid4","lightpink1"),c(-1:8))