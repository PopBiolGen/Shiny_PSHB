
#### UI object ####
## (Can put this in separate script, then 'source' here in app.R)

library(shiny)
library(bslib)
library(leaflet)
library(dplyr)
library(lubridate)
library(maps)
library(ozmaps)
library(sf)
library(zoo)
library(ggplot2)
library(grid)
library(gridExtra)
library(markdown)
library(knitr)
#sf_oz <- subset(ozmap("country"))

source("src/TPCFunctions.R")
source("src/modelFunctions_tyears.R")
source("src/growthFun.R")
source("src/plotFun.R")
source("src/collapse_monthsFun.R")

rmdfiles <- c("Documentation.rmd")
sapply(rmdfiles, knit, quiet = T)

# Define UI ----
ui <- page_navbar(title = "PSHB Survey Planner", # Separate tab with Readme
    
    nav_panel(title = "Survey planner",
              
              layout_columns(
                col_widths = c(5,7),
  
 card(span("Select a location",
            style = "font-size:20px"),
       leafletOutput("map"),
  
  sliderInput("weeks",
              label = span("Select number of survey weeks",
                           style = "font-size:20px"),
              min = 1, max = 52, value = 52),
  
  
  span(textOutput("selected_values"), style = "font-size:10px") # Disaply selected coords
  ),
  
  card(plotOutput("plot"))
 
 )
),

nav_panel(title = "Read me",
          
          fluidRow(
            div(withMathJax(includeMarkdown("Documentation.md")), style = "font-size: 15px;"),
            
            
            # you can add input selectors here as needed
            img(src = "PBG_Curtin_Logo.png", width="100%")
          )
)
  
  )


# Define server logic ----

server <- function(input, output, session) {
  
  ### Initial Map
  output$map <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%  # OpenStreetMap tiles
      setView(lng = 133, lat = -25, zoom = 3)  # Australia-ish
  })
  
  ### Store clicked coordinates
  clicked_point <- reactiveVal(NULL)
  
  ### When user clicks on the map
  observeEvent(input$map_click, {
    
    lat <- input$map_click$lat
    lng <- input$map_click$lng
    
    clicked_point(c(lat = lat, lon = lng)) # Save Coords
    
    # Update map: clear old markers, add new one
    leafletProxy("map") %>%
      clearMarkers() %>%
      addMarkers(lng = lng, lat = lat)
  })

  
  # Print values (checking working)
  output$selected_values <- renderText({
    req(clicked_point())
    
    paste(
      "Latitude =",  round(clicked_point()[["lat"]], 4),
      "| Longitude =", round(clicked_point()[["lon"]], 4))
  })
  
  output$plot <- renderPlot({
    
    req(clicked_point(), input$weeks)
    
    locLat <- clicked_point()[["lat"]]
    locLong <- clicked_point()[["lon"]]
    surv_weeks <- input$weeks
    
    plot_fun(locLat,
             locLong,
             surv_weeks)
    
  })
  
  
}


# Run the app ----
shinyApp(ui = ui, server = server)