#### E.g. Plot function ####

# How to display weeks/months to survey?
# https://www.royfrancis.com/calendar-plot-shiny-app-and-dynamic-ui/
# Fill out calendar
# Summarise months (with x or more days in them)

plot_fun <- function(locLat, locLong, surv_weeks){
  
  PAgrow_cum <- growth_loc_fun(locLat, locLong) # Calc cumulative growth values
  
  ### Subset highest growth rates, according to number of survey weeks
  top_grow <- PAgrow_cum %>%
    slice_max(order_by = cum_grow, # highest cumulative growth days, for number of survey weeks
            n = ifelse(surv_weeks < 52,
                       surv_weeks * 7, # weeks to days
                       366)) %>% # Assume if 52 weeks = 366 days (not 364)
  mutate(month = format(date, "%b"))

# Create ordered list of months
  months_ok <- top_grow %>%
    mutate(month = droplevels(month(date, label = TRUE, abbr = TRUE))) %>%
    count(month) %>%
    filter(n >= ifelse(surv_weeks > 2,
                       7,
                       2)) %>%
    arrange(month) %>%
    pull(month)
    
growth_plot <- ggplot()+
  
  scale_x_date(date_labels = "%b",
               date_breaks = "1 month",
               limits = c(min(PAgrow_cum$date)-1,
                          max(PAgrow_cum$date)+1),
               name=NULL) + # want to remove 2nd 01-Jan
  
  scale_y_continuous(name = "PSHB growth\n(30-day sum)\n",
                     limits = c(min(PAgrow_cum$cum_grow) - 0.05,
                          max(PAgrow_cum$cum_grow) + 0.05)) +
  
  geom_rect(data=top_grow, inherit.aes=FALSE,
            aes(xmin= c(date - 0.5),
                xmax=c(date + 0.5), 
                ymin= min(PAgrow_cum$cum_grow) - 0.05, 
                ymax= max(PAgrow_cum$cum_grow) + 0.05), 
            fill="lightblue", alpha=0.8)+
  
  geom_line(data = PAgrow_cum,
            aes(x = date, y = cum_grow),
            lwd=1.8)+
  
 # ggtitle(label = collapse_months(months_ok)) + # Call collapse function (lists months as title)
  
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(hjust = 0.05, size = 12),
        axis.title.y = element_text(size = 16))

# CALENDAR PLOT
month_data <- data.frame(
  month_num = 1:12,
  month_abb = month.abb,
  row = rep(3:1, each = 4),
  col = rep(1:4, 3)
)
df <- data.frame(month_abb = c(unique(top_grow$month)),
                 surv = "y")
month_data <- merge(month_data, df, by='month_abb', all.x=T)
month_data[is.na(month_data)] <- "n" 
month_data$surv <- factor(month_data$surv, levels = c("y", "n")) # make sure y always coloured blue

cal_plot <- ggplot(month_data, aes(x = col, y = row)) +
#  geom_tile(color = "white", size = 1, fill = "lightgrey") +  # Create squares
  geom_tile(data = na.omit(month_data),
            aes(fill = surv), 
            color = "white", size = 1) +
  geom_text(aes(label = month_abb), color = "white", size = 6) + # Add labels
  scale_fill_manual(values = c("lightblue", "lightgrey")) + # Colour gradient
  theme_void() + # Clean background
  theme(legend.position = "none") # Remove legend


plot_grid <- grid.arrange(cal_plot, growth_plot,
             nrow=2,
             heights = c(2, 1),
             top = textGrob(paste("Best months for PSHB surveys:",
                                  collapse_months(months_ok)),
                            gp=gpar(fontsize=20))) 
                                   # font = 3,
                                   # fontface = 'bold.italic',
                                    #col = "lightblue")))

return(plot_grid)
}