
months_fun <- function(locLat, locLong, surv_weeks){
  
  PAgrow_cum <- growth_loc_fun(locLat, locLong) # Calc cumulative growth values
  
  ### Subset highest growth rates, according to number of survey weeks
  top_grow <- PAgrow_cum %>%
    slice_max(order_by = cum_grow, # highest cumulative growth days, for number of survey weeks
              n = ifelse(surv_weeks < 52,
                         surv_weeks * 7, # weeks to days
                         366)) %>% # Assume if 52 weeks = 366 days (not 364)
    mutate(month = format(date, "%b"))

months_ok <- top_grow %>%
  mutate(month = droplevels(month(date, label = TRUE, abbr = TRUE))) %>%
  count(month) %>%
  filter(n >= 7) %>%
  arrange(month) %>%
  pull(month)

month_string <- collapse_months(months_ok)

return(month_string)

}
