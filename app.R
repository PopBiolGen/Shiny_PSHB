
#### PSHB Survey Spp ####

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
library(shinycssloaders)

source("src/TPCFunctions.R")
source("src/modelFunctions_tyears.R")
source("src/growthFun.R")
source("src/plotFun.R")
source("src/collapse_monthsFun.R")

#rmdfiles <- c("Documentation.rmd")
#sapply(rmdfiles, knit, quiet = T)

# Define UI ----
ui <- page_navbar(title = "PSHB Survey Planner", # Separate tab with Readme
    
    nav_panel(title = "Survey planner",
              
              layout_columns(
                col_widths = c(5,7),
  
 card(span("Select a location",
            style = "font-size:20px"),
       leafletOutput("map"),
  
  radioButtons("survey_mode",
               label = span("Survey type", style = "font-size:20px"),
               choices = c("Survey period" = "period",
                           "Multiple surveys" = "multiple"),
               selected = "period",
               inline = TRUE),

  conditionalPanel(
    condition = "input.survey_mode == 'period'",
    sliderInput("weeks",
                label = span("Number of survey weeks", style = "font-size:20px"),
                min = 1, max = 52, value = 10)
  ),

  conditionalPanel(
    condition = "input.survey_mode == 'multiple'",
    sliderInput("n_surveys",
                label = span("Number of surveys per year", style = "font-size:20px"),
                min = 1, max = 12, value = 4)
  ),
  
  
  span(textOutput("selected_values"), style = "font-size:10px") # Disaply selected coords
  ),
  
  card(withSpinner(plotOutput("plot")))
 
 )
),

nav_panel(title = "Read me",
          
          fluidRow(
            div(withMathJax(includeMarkdown("Documentation.md")), style = "font-size: 17px;"),
            
            
            # you can add input selectors here as needed
            tags$a(img(src = "PBG_Curtin_Logo.png", width="800px"),
            href="https://popbiolgenomics.org/")
          )
)
  
  )


# Define server logic ----

server <- function(input, output, session) {
  
  ### Initial Map
  output$map <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%  # OpenStreetMap tiles
      setView(lng = 135, lat = -26, zoom = 3.5)  %>% # Centre on Australia
      addMarkers(lng = 115.830, lat = -31.963)
  })
  
  ### Store clicked coordinates as reactive value
  clicked_point <- reactiveVal(c(lat = -31.963, # Start on King's Park, Perth
                                 lon = 115.830))
  
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

  
  # Print coord values
  output$selected_values <- renderText({
    req(clicked_point())
    
    paste(
      "Latitude =",  round(clicked_point()[["lat"]], 4),
      "| Longitude =", round(clicked_point()[["lon"]], 4))
  })
  
  # Expensive reactive: fetch climate data + run population model
  # Keyed on location only — not re-run when survey settings change
  grow_data <- eventReactive(clicked_point(), {
    tryCatch(
      list(ok = TRUE,
           data = growth_loc_fun(clicked_point()[["lat"]],
                                 clicked_point()[["lon"]])),
      error = function(e) {
        if (grepl("HTTP \\(502\\)", e$message))
          return(list(ok = FALSE, data = NULL))
        stop(e)
      }
    )
  })

  # Cheap reactive: re-renders whenever survey settings change, using cached data
  output$plot <- renderPlot({

    result <- grow_data()

    validate(
      need(result$ok,
           "Climate data from SILO may be unavailable between 11am and 1pm (Brisbane time) each Wednesday and Thursday to allow for essential system maintenance")
    )

    plot_fun(
      grow_data  = result$data,
      surv_weeks = input$weeks,
      mode       = input$survey_mode,
      n_surveys  = input$n_surveys
    )

  })
  
  
}


# Run the app ----
shinyApp(ui = ui, server = server)