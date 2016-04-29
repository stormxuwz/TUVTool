library(shiny)
library(dygraphs)
library(leaflet)
library(RColorBrewer)
library(rgl)
library(rglwidget)
library(plotly)
source("config.R")

shinyUI(
	navbarPage("TRIAXUS",
		### Calculation Pannel
		tabPanel("Calculate",
			sidebarLayout(
				sidebarPanel(
					fileInput("newFile","Select New File"),
					fileInput("PreviousFile","Select Previous Result File"),	

					selectInput("detrendingMethod", 
  						label = h5("Detrending Method"), 
        				choices = list("Thin Plate Spline" = "tps", "Linear Regression" = "linear"), 
        				selected = "lm"),

					selectInput("InterpolationMethod", 
  						label = h5("Interpolation Method"), 
        				choices = list("Thin Plate Spline" = "tps", "Kriging" = "krige","IDW"="idw"), 
        				selected = "idw"),
					
					helpText(h5("Enter the time anchor index of sensors")),
					column(6,
						numericInput("Seabird_anchor_index", 
  						label = h6("Seabird"), 
        				value = 0)
					),

					column(6,
						numericInput("BBE_anchor_index", 
  						label = h6("FluoroProbe"), 
  						value = 0)
					),

					helpText(h5("Enter the range of data to analyze")),
					column(6,
						numericInput("startIndex", 
  						label = h6("From"), 
        				value = 0)
					),
						
					column(6,
						numericInput("endIndex", 
  						label = h6("To"), 
  						value = 0)
					),
        			

					hr(),
					actionButton("readButton", "Preprocessing"),
					actionButton("goButton", "Calculate")
				),

				mainPanel(
					tabsetPanel(
						id = "mapFilePage",
						tabPanel("Input Check", 
							# selectInput("varToVis", 
				            	# label = h5("Variable"), 
				            	# choices = list("DO" = "DO"), 
				            	# selected = "DO"),
							
							column(6,
								numericInput("startIndex_check", 
  								label = h6("From"), 
        						value = 0)
							),
						
							column(6,
								numericInput("endIndex_check", 
  								label = h6("To"), 
  								value = 100)
							),
							textInput("suggestedAnchorIndex","Suggested Anchor Index"),
							plotOutput("rawData_check")
						),
						tabPanel("Map", leafletOutput("calMap"))
					)
				)
			)
		),

		### Visualization Panel
		tabPanel("Visualize",
			sidebarPanel(
				fileInput("resultFile","Select the result file to visualize"),
				selectInput("varToVis", 
				            label = h5("Variable"), 
				            choices = list("DO" = "DO"), 
				            selected = "DO"),
				selectInput("pathToVis",
							label = h5("Path to Visualize"), 
				            choices = list(), 
				            selected = "tmp")

       			# leafletOutput("visMap")
			),
			mainPanel(
				tabsetPanel(
					# tabPanel("visMap",leafletOutput("visMap")),
					tabPanel("Raw", 
						plotOutput("rawVis")
					),
					tabPanel("2D Interpolation", 
						plotOutput("twoD_Vis"),
						checkboxInput("hotspot2d","hotspot")),

					tabPanel("3D Visualize",
						column(8,
   							rglwidgetOutput("threeDMap"),
							checkboxInput("hotspot3d","hotspot")
  						),
						column(4,
							plotOutput("colorBar")
						)
					),
					tabPanel("Clustering", 
						column(6,	
							selectizeInput("varForClustering", label=h5("Selected variable for clustering"), choices=NULL, selected = NULL, multiple = TRUE, options = NULL)
						),
						# sliderInput("clusteringNum", "clusteringNum", min = 2, max = 5, value = 3, step= 1),
						column(6,
							selectInput("clusteringNum", label = h5("Cluster Num"), choices = list("2"=2,"3"=3,"4"=4,"5"=5,"6"=6), selected = 2)
						),
						actionButton("startClustering", "Start Clustering"),

						column(9,
							rglwidgetOutput("clustering")
						),
						column(3,
							plotOutput("colorBarCluster")
						),
						# column(6,
							plotOutput("boxplot")
						# )
					)
					# tabPanel("Clustering",rglwidgetOutput("threeDMap2"),plotOutput("statistics"))
				)
				
			)
		)
	)
)