library(shiny)
library(leaflet)
shinyUI(fluidPage(
  titlePanel("Pollution App"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Please Select the data you are interested in"),
      
      checkboxGroupInput('Statistics', 
                         label = h4("Statistics"), 
                         choices = list("Hourly Averages" = 1, 
                                        "Daily Hourly Max" = 'maxhr', "Rolling 8 Hour Averages" = 'eighthrrolling', "8 Hour Max" = 'eighthrmax'),
                         selected = 1),
      checkboxGroupInput('Pollutants', 
                         label = h4("Pollutants"), 
                         choices = list("CO" = "CO", 
                                        "SO2" = "SO2", "NO" = "NO", "NO2" = "NO2", "NOx" = "NOx", "O3" = "O3", "PM" = "PM", "Temperature" = "Temperature", "Wind Speed" = "Wind Speed", "Wind Direction" = "Wind Direction", "Humidity" = "Humidity", "Solar Radiation" = "Solar Radiation"),
                         selected = 1),

      dateRangeInput("dates", label = h3("Date range")),
      actionButton("do", "Generate Statistics")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Sites Map", leafletOutput("smap")), 
        tabPanel("Hourly Averages", tableOutput("table1")),
        tabPanel("Daily Hourly Max", tableOutput("table2")),
        tabPanel("Rolling 8 Hour Averages", tableOutput("table3")),
        tabPanel("Daily 8 Hour Maxes", tableOutput("table4")),
        tabPanel("Let's See if this is working", textOutput("hmm"))
      )
    )   
  )
))