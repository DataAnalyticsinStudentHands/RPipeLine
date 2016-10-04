library(shiny)
library(leaflet)
data=read.csv("data/2007_raw.csv",nrows=100000)
latlong=read.csv("data/latlongAQS.csv")
source("airpollutionstats.R")
source("converter.R")
source("hourlyaverage.R")
source("maxhour.R")
source("roll8avg.R")
latitude=vector()
longitude=vector()
shinyServer(function(input, output) {
  output$smap <- renderLeaflet({
    (smap <- leaflet() %>% addTiles()) 
    smap %>% setView(lng = -99.9018, lat = 31.9686, zoom = 5) # set centre and extent of map
    addMarkers(smap,lng=latlong$lng, lat=latlong$lat)
  })

  output$hmm <- renderText({
      latitude=c(latitude,input$smap_marker_click[3])
      longitude=c(longitude,input$smap_marker_click[4])
      paste(latitude)
  })

  observeEvent(input$do, {
    output$table1 <- renderTable({ 
      PStats=airpollutionstats(data=data,region=NA,pollutants=input$Pollutants,AQS.Code=NA,stats=input$Statistics)
      PStats[[1]]
    })
    output$table2 <- renderTable({ 
      PStats=airpollutionstats(data=data,region=NA,pollutants=input$Pollutants,AQS.Code=NA,stats=input$Statistics)
      PStats[[2]]
    })
    output$table3 <- renderTable({ 
      PStats=airpollutionstats(data=data,region=NA,pollutants=input$Pollutants,AQS.Code=NA,stats=input$Statistics)
      PStats[[3]]
    })
    output$table4 <- renderTable({ 
      PStats=airpollutionstats(data=data,region=NA,pollutants=input$Pollutants,AQS.Code=NA,stats=input$Statistics)
      PStats[[4]]
    })
  })

})