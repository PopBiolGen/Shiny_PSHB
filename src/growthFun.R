##################################

#### Calculate Daily Growth Rate at Location

growth_loc_fun <- function(locLat, locLong){
  
  yearSim <- run_year(locLat = locLat, locLong = locLong, # Run original single pop model (using coords)
                      make_plot = FALSE)
  
  population_data <- yearSim$popDat # Population size over 3 years
  # Focus on pre-oviposition stage
  pre_adult_vec <- population_data[2,]
  pre_adult_growth_rate <- diff(log(pre_adult_vec)) # Daily growth rate (change in log population size)
  
  # Calculate cumulative population growth over 30 day window
  # I.e.extended periods of growth are weighted higher than single days of very high growth
  total.t <- ncol(population_data) # 3 years = 1098 days
  
  # Create new data frame of cumulative growth rates
  PAgrow_cum <- data.frame(t = c(1:total.t),
                           growth = c(NA, # No growth on t = 1
                                      pre_adult_growth_rate))
  
  PAgrow_cum <- PAgrow_cum %>% # Sum growth over rolling 30 day window
    mutate(cum_grow = rollapply(PAgrow_cum$growth, 
                                width = 30, # number of previous days in window
                                FUN = sum,
                                align = "right", fill = NA)) %>%
    filter(t >= 367 & # Subset year 2
             t <= 732) %>%
    
    mutate(t = t - 366, # convert t to 1-366
           date = as.Date(t - 1, # Add dates (starts t = 0)
                          origin = "2000-01-01"), # Use leap year as dummy date
           date_dm = format(date, "%d-%b")) # Convert to day/month
  
  return(PAgrow_cum)
  
}

# Other stages:
#juv_vec <- population_data[1,]
#ad_vec <- population_data[3,]
#juv_growth_rate <- diff(log(juv_vec))
#ad_growth_rate <- diff(log(ad_vec))
#juv_growth_rate_smoothed <- stats::filter(juv_growth_rate, rep(1/window_size, window_size), sides = 2)
#ad_growth_rate_smoothed <- stats::filter(ad_growth_rate, rep(1/window_size, window_size), sides = 2)
# Rolling average to smooth growth rate:
#window_size = 5 
#pre_adult_growth_rate_smoothed <- stats::filter(pre_adult_growth_rate, rep(1/window_size, window_size), sides = 2) # Smooth growth rate using moving average


######


