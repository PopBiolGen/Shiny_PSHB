#### Plot function ####

# Find the optimal start offset for n equally-spaced surveys
# Returns the day numbers (t = 1:366) of the optimal survey days
find_optimal_surveys <- function(grow_df, n_surveys) {
  spacing   <- 366 / n_surveys
  n_starts  <- max(1, floor(spacing))

  scores <- sapply(seq_len(n_starts), function(start) {
    days <- round(start + (seq_len(n_surveys) - 1) * spacing)
    days <- pmin(pmax(days, 1), 366)
    sum(grow_df$cum_grow[grow_df$t %in% days], na.rm = TRUE)
  })

  best_start <- which.max(scores)
  days <- round(best_start + (seq_len(n_surveys) - 1) * spacing)
  pmin(pmax(days, 1), 366)
}


plot_fun <- function(grow_data, surv_weeks = 10, mode = "period", n_surveys = 4) {

  PAgrow_cum <- grow_data

  # ---- Shared axis / theme settings ----
  x_scale <- scale_x_date(
    date_labels = "%b",
    date_breaks = "1 month",
    limits = c(min(PAgrow_cum$date) - 1, max(PAgrow_cum$date) + 1),
    name = NULL
  )

  y_scale <- scale_y_continuous(
    name   = "PSHB growth\n(30-day sum)\n",
    limits = c(min(PAgrow_cum$cum_grow, na.rm = TRUE) - 0.05,
               max(PAgrow_cum$cum_grow, na.rm = TRUE) + 0.05)
  )

  base_theme <- theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x  = element_text(hjust = 0.05, size = 12),
    axis.title.y = element_text(size = 16)
  )

  # ---- Survey period mode (original behaviour) ----
  if (mode == "period") {

    top_grow <- PAgrow_cum |>
      slice_max(order_by = cum_grow,
                n = ifelse(surv_weeks < 52, surv_weeks * 7, 366)) |>
      mutate(month = format(date, "%b"))

    months_ok <- top_grow |>
      mutate(month = droplevels(month(date, label = TRUE, abbr = TRUE))) |>
      count(month) |>
      filter(n >= ifelse(surv_weeks > 2, 7, 2)) |>
      arrange(month) |>
      pull(month)

    growth_plot <- ggplot() +
      x_scale + y_scale +
      geom_rect(data = top_grow, inherit.aes = FALSE,
                aes(xmin = date - 0.5, xmax = date + 0.5,
                    ymin = min(PAgrow_cum$cum_grow, na.rm = TRUE) - 0.05,
                    ymax = max(PAgrow_cum$cum_grow, na.rm = TRUE) + 0.05),
                fill = "lightblue", alpha = 0.8) +
      geom_line(data = PAgrow_cum, aes(x = date, y = cum_grow), lwd = 1.8) +
      base_theme

    survey_months <- unique(top_grow$month)
    plot_title    <- paste("Best months for PSHB surveys:", collapse_months(months_ok))

  # ---- Multiple surveys mode ----
  } else {

    survey_days <- find_optimal_surveys(PAgrow_cum, n_surveys)
    survey_data <- PAgrow_cum[PAgrow_cum$t %in% survey_days, ]

    growth_plot <- ggplot() +
      x_scale + y_scale +
      geom_line(data = PAgrow_cum, aes(x = date, y = cum_grow), lwd = 1.8) +
      geom_vline(data = survey_data, aes(xintercept = date),
                 color = "steelblue", lwd = 1.2, linetype = "dashed") +
      geom_point(data = survey_data, aes(x = date, y = cum_grow),
                 color = "steelblue", size = 4) +
      base_theme

    survey_months <- unique(format(survey_data$date, "%b"))
    plot_title    <- paste("Optimal survey dates:",
                           paste(format(survey_data$date, "%d %b"), collapse = " | "))
  }

  # ---- Calendar grid (shared) ----
  month_data <- data.frame(
    month_num = 1:12,
    month_abb = month.abb,
    row = rep(3:1, each = 4),
    col = rep(1:4, 3)
  )

  df         <- data.frame(month_abb = survey_months, surv = "y")
  month_data <- merge(month_data, df, by = "month_abb", all.x = TRUE)
  month_data[is.na(month_data)] <- "n"
  month_data$surv <- factor(month_data$surv, levels = c("y", "n"))

  cal_plot <- ggplot(month_data, aes(x = col, y = row)) +
    geom_tile(data = na.omit(month_data),
              aes(fill = surv),
              color = "white", size = 1) +
    geom_text(aes(label = month_abb), color = "white", size = 6) +
    scale_fill_manual(values = c("lightblue", "lightgrey")) +
    theme_void() +
    theme(legend.position = "none")

  # ---- Combine ----
  # Convert to grobs and align column widths so plot panels are horizontally flush
  g_cal    <- ggplotGrob(cal_plot)
  g_growth <- ggplotGrob(growth_plot)

  max_widths      <- grid::unit.pmax(g_cal$widths, g_growth$widths)
  g_cal$widths    <- max_widths
  g_growth$widths <- max_widths

  plot_grid <- grid.arrange(
    g_cal, g_growth,
    nrow    = 2,
    heights = c(1, 2),
    top     = textGrob(plot_title, gp = gpar(fontsize = 18))
  )

  return(plot_grid)
}
