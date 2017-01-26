# source("airpollutionstats.R")
# source("converter.R")
# source("hourlyaverage.R")
# source("maxhour.R")
# source("roll8avg.R")
source('pipeline.R')

pipeline.load_TCEQ_libs()

dataframe <- load.file_as_dataframe("examples/raw_data/TCEQ/2007_raw.csv", maximum_rows=100000)

stats <- TCEQ.generate_air_pollution_statistics(dataframe, region=NA, pollutants=NA, AQS_codes=NA, stats=NA)

# 
# shinyServer(function(input, output) {
#   output$smap <- renderLeaflet({
#     (smap <- leaflet() %>% addTiles()) 
#     smap %>% setView(lng = -99.9018, lat = 31.9686, zoom = 5) # set centre and extent of map
#     addMarkers(smap,lng=latlong$lng, lat=latlong$lat)
#   })
#   
#   output$hmm <- renderText({
#     latitude=c(latitude,input$smap_marker_click[3])
#     longitude=c(longitude,input$smap_marker_click[4])
#     paste(latitude)
#   })
#   
#   observeEvent(input$do, {
#     output$table1 <- renderTable({ 
#       PStats=airpollutionstats(data=data,region=NA,pollutants=input$Pollutants,AQS.Code=NA,stats=input$Statistics)
#       PStats[[1]]
#     })
#     output$table2 <- renderTable({ 
#       PStats=airpollutionstats(data=data,region=NA,pollutants=input$Pollutants,AQS.Code=NA,stats=input$Statistics)
#       PStats[[2]]
#     })
#     output$table3 <- renderTable({ 
#       PStats=airpollutionstats(data=data,region=NA,pollutants=input$Pollutants,AQS.Code=NA,stats=input$Statistics)
#       PStats[[3]]
#     })
#     output$table4 <- renderTable({ 
#       PStats=airpollutionstats(data=data,region=NA,pollutants=input$Pollutants,AQS.Code=NA,stats=input$Statistics)
#       PStats[[4]]
#     })
#   })
#   
# })
# 
# library(shiny)
# library(leaflet)
# shinyUI(fluidPage(
#   titlePanel("Pollution App"),
#   
#   sidebarLayout(
#     sidebarPanel(
#       helpText("Please Select the data you are interested in"),
#       
#       checkboxGroupInput('Statistics', 
#                          label = h4("Statistics"), 
#                          choices = list("Hourly Averages" = 1, 
#                                         "Daily Hourly Max" = 'maxhr', "Rolling 8 Hour Averages" = 'eighthrrolling', "8 Hour Max" = 'eighthrmax'),
#                          selected = 1),
#       checkboxGroupInput('Pollutants', 
#                          label = h4("Pollutants"), 
#                          choices = list("CO" = "CO", 
#                                         "SO2" = "SO2", "NO" = "NO", "NO2" = "NO2", "NOx" = "NOx", "O3" = "O3", "PM" = "PM", "Temperature" = "Temperature", "Wind Speed" = "Wind Speed", "Wind Direction" = "Wind Direction", "Humidity" = "Humidity", "Solar Radiation" = "Solar Radiation"),
#                          selected = 1),
#       
#       dateRangeInput("dates", label = h3("Date range")),
#       actionButton("do", "Generate Statistics")
#     ),
#     
#     mainPanel(
#       tabsetPanel(
#         tabPanel("Sites Map", leafletOutput("smap")), 
#         tabPanel("Hourly Averages", tableOutput("table1")),
#         tabPanel("Daily Hourly Max", tableOutput("table2")),
#         tabPanel("Rolling 8 Hour Averages", tableOutput("table3")),
#         tabPanel("Daily 8 Hour Maxes", tableOutput("table4")),
#         tabPanel("Let's See if this is working", textOutput("hmm"))
#       )
#     )   
#   )
# ))