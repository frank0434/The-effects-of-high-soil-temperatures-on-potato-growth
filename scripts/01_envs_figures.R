source("scripts/00_setup.R")
# Load weather data 

daily_weather_data <- fread("data/daily_weather.csv")
# rain  -------------------------------------------------------------------

p1 <- ggplot(daily_weather_data, aes(x = DAP)) +
  geom_col(aes(y = Rain), alpha = 0.7, fill = "black") +
  geom_line(aes(y = Tmean, color = "Tmean"), linewidth = lw - 0.3, show.legend = FALSE) +
  geom_line(aes(y = Tmax, color = "Tmax"), linewidth = lw - 0.3, linetype = "dashed") +
  geom_line(aes(y = Tmin, color = "Tmin"), linewidth = lw - 0.3, linetype = "dotted") +
  scale_y_continuous(
    name = "Temperature (°C)",expand = c(0,0),
    sec.axis = sec_axis(~ ., name = "Rainfall (mm)"), 
    limits = c(0, 40)
  ) +
  scale_x_continuous(expand = c(0,0),limits = c(1,120)) +
  geom_rect(
    data = treatment_rect_dt,
    aes(xmin =  DAP_min, xmax = DAP_max, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = colors_temp, guide = "none") +
  scale_color_manual(values = c("Tmean" = "black", "Tmax" = "black", "Tmin" = "black")) +
  labs(color = "Legend") +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  facet_grid(~Season) +
  theme(# strip.background = element_blank(),
        legend.box.margin = margin(b = -5, unit = "mm"),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.x = unit(5, "mm"))

# Radiation plot
p2 <- ggplot(daily_weather_data, aes(x = DAP, y = IRRAD)) +
  geom_line(color = "grey50", alpha = 0.8, linewidth = lw) +
  geom_line(aes(y = Cumulative_Radiation/80), color = "black", linewidth = lw - 0.3) +
  scale_y_continuous(
    name = expression(atop(Daily~Radiation, (MJ~m^{-2}))),
    sec.axis = sec_axis(~ . * 80, name = expression(atop(Cumulative~Radiation, (MJ~m^{-2})))), 
    limits = c(0, 35), expand = c(0,0)
  ) + 
  scale_x_continuous(name = "Days After Planting", expand = c(0,0), limits = c(0,120)) +
  geom_rect(
    data = treatment_rect_dt,
    aes(xmin =  DAP_min, xmax = DAP_max, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15,
    inherit.aes = FALSE) +
  scale_fill_manual(values = colors_temp, guide = "none") +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  theme(panel.spacing.x = unit(5, "mm"),
        panel.grid = element_blank())+
  facet_grid(~Season)

weather_plot <- (p1 +
                   theme(axis.text.x = element_blank(),
                         axis.title.x = element_blank()))/ (
                           p2 + 
                             theme(strip.text = element_blank(),
                                   strip.background = element_blank())) + 
  plot_layout(guides = "collect", axes = "collect", heights = c(0.95, 1)) & 
  theme(legend.position = "top")
print(weather_plot)
# treatment period weather metrics ----------------------------------------
# Get DAP ranges for each treatment period

trt_ranges <- treatments_DAP_arrows[
  variable %in% c("2024_HE", "2024_HTI", "2025_HTI", "2025_HB"),
  .(DAP_start = DAP[1], DAP_end = DAP[2]), by = variable
][, Season := as.integer(substr(variable, 1, 4))]

# Non-equi join: match daily_weather_data rows falling within each period
air_periods <- daily_weather_data[trt_ranges, on = .(Season, DAP >= DAP_start, DAP <= DAP_end), nomatch = NULL]
# Mean temperature per treatment period
mean_temps <- air_periods[, .(Mean_Temp = round(mean(Tmean, na.rm = TRUE), 1),
                              sd_Temp = sd(Tmean), 
                              Rain = sum(Rain) |> round(),
                              Rain_mean = sum(Rain)/.N,
                              Days = .N), by = .(Treatment = variable)]
print(mean_temps)
rad_periods <- daily_weather_data[trt_ranges, on = .(Season, DAP >= DAP_start, DAP <= DAP_end), nomatch = NULL]
sum_rad_periods <- rad_periods[, .(rad_mean = mean(IRRAD) |> round(),
                                   rad_sum= sum(IRRAD) |> round()), by = .(Treatment = variable)]
weather_metric_periods <- mean_temps[sum_rad_periods, on = "Treatment"]

# soil and treatment verification -----------------------------------------

# Summary stats

soil_temperature_daily <- soil_temperature_daily[!(Season == 2025 & 
                                                     Treatment == "HB" & 
                                                     PlotID == 8)]
stats <- get_env_summary_stats(soil_temperature_daily, soil_moisture)
print(stats)
#  what is the proportion of time above 22 C in the 20 cm soil temperature before the heat events?
supra_threshold <- soil_temperature_daily[Season == 2024 & DAP > 39 & DAP < 46 & Depth == 20][, 
                        supra_threshold := mean_temp > 22, 
                       by = Treatment]
air_temp_daily[Season == 2024 & DAP > 39 & DAP < 46][, 
                        supra_threshold := mean_temp > 22, 
                       by = Treatment]
supra_threshold[supra_threshold == TRUE, .N/4, by = Treatment]

# soil temperature time series ------------------------------------------------------
soil_temperature_daily <- fread("data/soil_temperature_daily.csv")
# calculate mean and sd time series by Treatment/Depth/Season
soil_temp_ts <- soil_temperature_daily[, .(
    mean_temp = mean(mean_temp, na.rm = TRUE),
    sd_temp = sd(mean_temp, na.rm = TRUE)
    ), by = .(Season, Depth, Treatment, Date, DAP)
    ][order(Season, Depth, Treatment, Date, DAP)]

# 2-day moving average for mean and upper/lower bounds
soil_temp_ts[, `:=`(
  mean_ma = frollmean(mean_temp, n = 2, align = "right"),
  upper_ma = frollmean(mean_temp + sd_temp, n = 2, align = "right"),
  lower_ma = frollmean(mean_temp - sd_temp, n = 2, align = "right")
), by = .(Season, Depth, Treatment)]
soil_temp_ts[, Depth := paste0(Depth, " cm")]
fsize <- 12
# Plot mean +/- SD  -----------------
panel_tag_soil <- data.table(Season = c(2024, 2025), Depth = c("20 cm", "20 cm"),
                             label = c("a", "b"))
# calculate mean and sd time series by Treatment/Depth/Season
# !!!!! PlotID 8 in season2025 had a defected cable
soil_temp_ts <- soil_temperature_daily[, .(
  mean_temp = mean(mean_temp, na.rm = TRUE),
  sd_temp = sd(mean_temp, na.rm = TRUE)
), by = .(Season, Depth, Treatment, Date, DAP)
][order(Season, Depth, Treatment, Date, DAP)]

# 2-day moving average for mean and upper/lower bounds
soil_temp_ts[, `:=`(
  mean_ma = frollmean(mean_temp, n = 2, align = "right"),
  upper_ma = frollmean(mean_temp + sd_temp, n = 2, align = "right"),
  lower_ma = frollmean(mean_temp - sd_temp, n = 2, align = "right")
), by = .(Season, Depth, Treatment)]
soil_temp_ts[, Depth := paste0(Depth, " cm")]
fsize <- 12
# Plot mean +/- SD  -----------------
panel_tag_soil <- data.table(Season = c(2024, 2025), Depth = c("20 cm", "20 cm"),
                             label = c("a", "b"))
# Keep AC over the whole season; keep each heated treatment only inside its own
# heating window
soil_temp_20 <- soil_temp_ts[Depth != "40 cm"]
soil_temp_sd_dt <- rbind(
  soil_temp_20[Treatment == "AC"],
  soil_temp_20[treatment_rect_dt, on = .(Season, Treatment), nomatch = NULL
               ][DAP >= DAP_min - 2 & DAP <= DAP_max + 2
                 ][, c("DAP_min", "DAP_max") := NULL]
)
# figure for the paper
soil_temp_sd_p <- soil_temp_sd_dt |>
  ggplot(aes(DAP, mean_temp)) +
  geom_rect(
    data = treatment_rect_dt,
    aes(xmin =  DAP_min, xmax = DAP_max, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15,
    inherit.aes = FALSE
  ) +
  geom_hline(yintercept = c(14, 22), color = "#e74c3c")+
  geom_ribbon(aes(ymin = mean_temp - sd_temp, 
                  ymax = mean_temp + sd_temp,
                  fill = Treatment),
              alpha = 0.5) +
  geom_line(aes(linetype = Treatment, color = Treatment), linewidth = lw - .3) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 120)) +
  scale_y_continuous(limits = c(10, 35))+
  facet_grid( ~ Season) +
  scale_linetype_manual(
    values = linetype_temp,
    breaks = legend_order
  ) +
  scale_color_manual(
    values = colors_temp,
    breaks = legend_order
  ) +
  scale_fill_manual(
    values = colors_temp,
    breaks = legend_order
  ) +
  geom_text(data = panel_tag_soil,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.5,    vjust = 1.3,    size = 5,  inherit.aes = FALSE) +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  theme(legend.position = "top",
        panel.spacing.x = unit(5, "mm"),
        # axis.text.x = element_blank(),
        # axis.title.x =  element_blank(),
        panel.grid = element_blank(),
        legend.key.width = unit(15, "mm"),
      legend.title = element_text(margin = margin(r = 15))) +
  labs(y = "Soil temperature (°C)", x = "Days after planting")
soil_temp_sd_p
# canopy temperature as a verification of non stressed canopy -------------
canopy_temp <- fread("data/canopy_temperature.csv")
flightime <- fread("data/flightime.csv")
# Add date as character for labeling
canopy_temp[, Date_label := format(Date, "%b %d")]
canopy_temp_sum <- canopy_temp[, .(mean = mean(canopy_temp),
                                   sd = sd(canopy_temp)), 
                               by = .(Season, Date, DAP, Treatment)]
# visualise the canopy temperature time series 
canopy_temp_sum |>
  ggplot(aes(DAP, mean, color = Treatment)) +
  geom_point(size = ps, alpha = 0.6, position = position_dodge(width = 5)) +
    geom_rect(
    data = treatment_rect_dt,
    aes(xmin =  DAP_min, xmax = DAP_max, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 3,
                position = position_dodge(width = 5)) +
  facet_grid(~Season) +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  theme(legend.position = "top",
        strip.background = element_blank()) +
  labs(x = "Days After Planting",
       y = "Canopy Temperature (°C)",
       color = "Treatment") +
  scale_y_continuous(
    name = "Surface canopy temperature (°C)",
    limits = c(0, 35),
    expand = c(0, 0)
  ) +
  scale_shape_manual(values = point_shape, breaks = c("AC", "HE", "HTI", "HB")) +
  scale_color_manual(values = colors_temp, breaks = c("AC", "HE", "HTI", "HB")) +
  scale_fill_manual(values = colors_temp, breaks = c("AC", "HE", "HTI", "HB")) +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  theme(panel.grid = element_blank()) +
  labs(x = "Days after planting")


