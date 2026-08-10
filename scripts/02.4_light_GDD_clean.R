# ==============================================================================
# Canopy light interception (fPAR) on a growing-degree-day axis
#
# Figures  fig_fPAR_GDD      treatment-level fitted curves, both seasons
#          figS_LI_GDD_24    per-plot fitted vs observed, 2024
#          figS_LI_GDD_25    per-plot fitted vs observed, 2025
#
# Stats    curve parameters (vmax, t1, tm1, t2, te, delta_t) and total
#          intercepted radiation: lm(y ~ Block + Treatment) per season,
#          Type II SS, Tukey-adjusted pairwise contrasts on the emmeans.
#
# Inputs   data/LI_cleaned.csv     cleaned per-plot fPAR  
#          data/daily_weather.csv  daily Tmean and IRRAD
#
# Self-contained: reads only the two CSVs above, no {targets} store required.
# ==============================================================================
source("scripts/functions.R")   # constants, palettes, curve fitting

# One season, one response: Block as a fixed blocking factor, no interaction.
# Type II SS because 2025 is unbalanced (one HB plot lost); for the balanced
# 2024 season Type II is identical to Type I/III.

fit_treatment_lm <- function(data, response) {
  data  <- droplevels(data)
  model <- lm(stats::as.formula(paste(response, "~ Block + Treatment")), data = data)
  list(model   = model,
       anova   = as.data.table(car::Anova(model, type = "II"), keep.rownames = "term"),
       emmeans = emmeans::emmeans(model, specs = ~ Treatment))
}

# ==============================================================================
# PART 1: Cumulative simple GDD from emergence (Tbase = 2°C)
# ==============================================================================
daily_weather <- fread("data/daily_weather.csv")
daily_weather[, Date := as.Date(Date)]
daily_weather[, DAE := fifelse(
  Season == 2024,
  as.integer(Date - emergence_date$`2024`),
  as.integer(Date - emergence_date$`2025`)
)]

setorder(daily_weather, Season, Date)
daily_gdd_join <- daily_weather[DAE >= 0]
daily_gdd_join[, cumGDD_DAE := cumsum(Tmean - 2), by = Season]
daily_gdd_join <- daily_gdd_join[, .(Season, Date, DAP, cumGDD_DAE, IRRAD)]

# ==============================================================================
# PART 2: Cleaned LI data on the GDD axis
# ==============================================================================
li_clean <- fread("data/LI_cleaned.csv")
li_clean[, Date := as.Date(Date)]

li_gdd_sum <- merge(
  li_clean,
  daily_gdd_join[, .(Season, Date, DAP, cumGDD_DAE)],
  by = c("Season", "Date", "DAP"), all.x = TRUE
)

# te_fixed: cumulative GDD at the leaf removal date, used to cap each fit
leaf_removal_gdd <- data.table(
  Season = c(2024L, 2025L),
  te_gdd = c(
    daily_gdd_join[Date == leaf_removal_dates$`2024`, cumGDD_DAE][1],
    daily_gdd_join[Date == leaf_removal_dates$`2025`, cumGDD_DAE][1]
  )
)
li_gdd_sum <- merge(li_gdd_sum, leaf_removal_gdd, by = "Season", all.x = TRUE)

li_gdd_sum[, ':='(
  Season    = factor(Season),
  Treatment = factor(Treatment, levels = legend_order),
  Block     = factor(Block),
  PlotID    = factor(PlotID)
)]
setorder(li_gdd_sum, Season, Treatment, Block, PlotID)

# Treatment (heating) windows converted from DAP to GDD units
rect_gdd_min <- daily_gdd_join[
  treatment_rect_dt, on = c("Season", "DAP == DAP_min")
][, .(Season, Treatment, DAP_min = DAP, DAP_max, xmin = cumGDD_DAE)]
rect_gdd <- daily_gdd_join[
  rect_gdd_min, on = c("Season", "DAP == DAP_max")
][, .(Season, Treatment, xmin, xmax = cumGDD_DAE)]

max_gdd <- li_gdd_sum[, ceiling(max(cumGDD_DAE, na.rm = TRUE))]

# ==============================================================================
# PART 3: Per-plot curve fitting on the GDD axis
# ==============================================================================
li_gdd_fits <- li_gdd_sum[
  LI > 0 & is.finite(cumGDD_DAE),
  .(data = list(.SD)),
  by = .(Season, Treatment, Block, PlotID, te_gdd)
]
setorder(li_gdd_fits, Season, Treatment, Block, PlotID)

li_gdd_fits[, starts := mapply(
  function(x, te) optimize_start_values_combined(x, y_col = "LI", x_col = "cumGDD_DAE", tmax = te),
  data, te_gdd, SIMPLIFY = FALSE
)]
li_gdd_fits[, fit := mapply(
  function(data, start, te)
    fit_combined(data, y_col = "LI", x_col = "cumGDD_DAE", te_cap = te, start_vals = start),
  data, starts, te_gdd, SIMPLIFY = FALSE
)]
li_gdd_fits[, r2 := sapply(fit, r2_nls)]

li_gdd_params <- li_gdd_fits[
  , .(sapply(fit, extract_combined_params, simplify = FALSE)),
  by = .(Season, Treatment, Block, PlotID)
][, unlist(V1, recursive = FALSE), by = .(Season, Treatment, Block, PlotID)]
li_gdd_params[, delta_t := t2 - t1]

# Per-plot prediction grid (step = 1 so axis break values are present)
prediction_grid_gdd <- generate_time_grid(
  li_gdd_sum[!is.na(cumGDD_DAE)], end_time = max_gdd, step = 1
)
setorder(prediction_grid_gdd, Season, Block, PlotID, cumGDD_DAE)

predictions_gdd <- merge(
  prediction_grid_gdd,
  li_gdd_params[, .(Season, Treatment, Block, PlotID, tm1, t1, t2, te, vmax)],
  by = c("Season", "Treatment", "Block", "PlotID"), all.x = TRUE
)
predictions_gdd[, prediction := combined_equation(cumGDD_DAE, tm1, t1, t2, te, vmax)]
predictions_gdd[!is.finite(prediction) | prediction < 0, prediction := 0]

# ==============================================================================
# PART 4: Statistical testing of the curve parameters
# ==============================================================================
canopy_param_cols <- c("vmax", "t1", "tm1", "t2", "te", "delta_t")

li_gdd_param_models_24 <- setNames(
  lapply(canopy_param_cols, function(p) fit_treatment_lm(li_gdd_params[Season == 2024], p)),
  canopy_param_cols
)
li_gdd_param_models_25 <- setNames(
  lapply(canopy_param_cols, function(p) fit_treatment_lm(li_gdd_params[Season == 2025], p)),
  canopy_param_cols
)

# Type II ANOVA tables
rbindlist(list(
  rbindlist(lapply(canopy_param_cols, function(p)
    li_gdd_param_models_24[[p]]$anova[, .(Season = 2024L, Parameter = p, term,
                                          `Sum Sq`, Df, `F value`, `Pr(>F)`)])),
  rbindlist(lapply(canopy_param_cols, function(p)
    li_gdd_param_models_25[[p]]$anova[, .(Season = 2025L, Parameter = p, term,
                                          `Sum Sq`, Df, `F value`, `Pr(>F)`)]))
)) |> fwrite("results/LI_GDD_parameters_anova_typeII.csv", bom = TRUE)

param_gdd_pair_24 <- rbindlist(lapply(canopy_param_cols, function(p)
  summarise_pairwise_table(li_gdd_param_models_24[[p]]$emmeans, 2024, "LI_GDD_nls", p)
))
param_gdd_pair_25 <- rbindlist(lapply(canopy_param_cols, function(p)
  summarise_pairwise_table(li_gdd_param_models_25[[p]]$emmeans, 2025, "LI_GDD_nls", p)
))
rbindlist(
  list(param_gdd_pair_24[, .(contrast, p.value, Season, Parameter, Significant)],
       param_gdd_pair_25[, .(contrast, p.value, Season, Parameter, Significant)]),
  use.names = TRUE, fill = TRUE
) |> fwrite("results/LI_GDD_parameters_pairwise_comparisons.csv", bom = TRUE)

# Parameter summary table (GDD units)
li_gdd_params_summary <- melt(
  li_gdd_params[, .(Season, Treatment, Block, PlotID, tm1, t1, t2, te, vmax, delta_t)],
  id.vars = c("Season", "Treatment", "Block", "PlotID")
)[, .(mean = mean(value), se = sd(value) / sqrt(.N)),
  by = .(Season, Treatment, variable)]

li_gdd_params_summary[variable == "vmax",
  formatted := sprintf("%.2f ± %.2f", mean, se)]
li_gdd_params_summary[variable != "vmax",
  formatted := sprintf("%.0f ± %.1f", mean, se)]

li_gdd_params_summary |>
  dcast(Season + Treatment ~ variable, value.var = "formatted") |>
  fwrite("results/LI_GDD_parameters_summary.csv", bom = TRUE)

# ==============================================================================
# PART 5: Treatment-level fits for plotting
# ==============================================================================
li_gdd_trt_fits <- li_gdd_sum[
  LI > 0 & is.finite(cumGDD_DAE),
  .(data = list(.SD)),
  by = .(Season, Treatment, te_gdd)
]
setorder(li_gdd_trt_fits, Season, Treatment)

li_gdd_trt_fits[, starts := mapply(
  function(x, te) optimize_start_values_combined(x, y_col = "LI", x_col = "cumGDD_DAE", tmax = te),
  data, te_gdd, SIMPLIFY = FALSE
)]
li_gdd_trt_fits[, fit := mapply(
  function(data, start, te)
    fit_combined(data, y_col = "LI", x_col = "cumGDD_DAE", te_cap = te, start_vals = start),
  data, starts, te_gdd, SIMPLIFY = FALSE
)]
li_gdd_trt_params <- li_gdd_trt_fits[
  , .(sapply(fit, extract_combined_params, simplify = FALSE)),
  by = .(Season, Treatment)
][, unlist(V1, recursive = FALSE), by = .(Season, Treatment)]

li_gdd_trt_params[, curve_dat := mapply(function(tm1, t1, t2, te, vmax) {
  time <- seq(1, max_gdd, 1)   # step = 1 ensures axis break values are present
  data.table(cumGDD_DAE = time,
             predictions = combined_equation(time, tm1, t1, t2, te, vmax))
}, tm1, t1, t2, te, vmax, SIMPLIFY = FALSE)]

LI_gdd_curve_dt <- li_gdd_trt_params[
  , unlist(curve_dat, recursive = FALSE), by = .(Season, Treatment)
][predictions >= 0][, ':='(Season = as.integer(as.character(Season)))]

# Secondary DAP axis: nearest observed day to each GDD break
gdd_axis_breaks <- seq(0, max_gdd, by = 200)
dap_axis_gdd <- rbindlist(lapply(gdd_axis_breaks[gdd_axis_breaks > 0], function(b) {
  daily_gdd_join[, .SD[which.min(abs(cumGDD_DAE - b))], by = Season][
    , .(Season, cumGDD_DAE = b, DAP)]
}))
emergence_dap_dt <- data.table(
  Season     = c(2024L, 2025L),
  DAP        = c(as.integer(emergence_date$`2024` - planting_dates$`2024`),
                 as.integer(emergence_date$`2025` - planting_dates$`2025`)),
  cumGDD_DAE = 0L
)
dap_axis_gdd <- rbindlist(list(dap_axis_gdd, emergence_dap_dt),
                          use.names = TRUE, fill = TRUE)

# Observed treatment means for the points and error bars
li_gdd_summary <- li_gdd_sum[
  , .(mean_LI = mean(LI, na.rm = TRUE),
      se_LI   = sd(LI, na.rm = TRUE) / sqrt(.N)),
  by = .(Season, Treatment, cumGDD_DAE)
][mean_LI > 0 & is.finite(cumGDD_DAE)]

# ==============================================================================
# PART 6: fig_fPAR_GDD
#
# 2025 t2 arrows mark the onset of canopy decline, one per treatment with a
# Tukey letter. The head is pinned to the plotted curve — x is the treatment-
# level fitted t2, y is that curve evaluated there — so only the tail needs
# setting.
# ==============================================================================
t2_letters_25 <- as.data.table(
  multcomp::cld(li_gdd_param_models_25$t2$emmeans, Letters = letters,
                adjust = "tukey", quiet = TRUE)
)[, .(Treatment = as.character(Treatment), letter = trimws(.group))]

# tail of each arrow, where the label sits — the only coordinates to tune
t2_arrow_tail <- data.table(
  Treatment = c("AC",  "HB", "HTI"),
  x_tail    = c(1200, 1200, 1200),
  y_tail    = c(1.05,    1, 0.95),
  lab_nudge = c(0, 0, 0),   # label offset from the tail x
  lab_hjust = c(0, 0, 0)
)

t2_arrows_25 <- li_gdd_trt_params[
  Season == "2025",
  .(Treatment = as.character(Treatment),
    x_head = t2,
    y_head = combined_equation(t2, tm1, t1, t2, te, vmax))
][t2_letters_25, on = "Treatment"][t2_arrow_tail, on = "Treatment"]

t2_arrows_25[, ':='(
  Season    = 2025L,
  Treatment = factor(Treatment, levels = legend_order),
  label     = paste(Treatment, letter)
)]

fig_fPAR_GDD <- LI_gdd_curve_dt |>
  ggplot(aes(cumGDD_DAE, predictions, colour = Treatment)) +
  # heating period rectangles
  geom_rect(
    data = rect_gdd,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15, inherit.aes = FALSE
  ) +
  geom_segment(
    data = t2_arrows_25,
    aes(x = x_head, xend = x_tail, y = y_head, yend = y_tail, colour = Treatment),
    arrow = arrow(length = unit(0.5, "cm"), type = "open", ends = "first"),
    linewidth = lw - 0.3, inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_text(
    data = t2_arrows_25,
    aes(x = x_tail + lab_nudge, y = y_tail, label = label,
        hjust = lab_hjust, colour = Treatment),
    vjust = 0.3, inherit.aes = FALSE, show.legend = FALSE,
    size = floor(fontsize / ggplot2::.pt) - 1
  ) +
  # treatment mean curves
  geom_line(aes(linetype = Treatment), linewidth = lw - .3) +
  # observed means with SE
  geom_point(data = li_gdd_summary,
             aes(x = cumGDD_DAE, y = mean_LI), size = 1, alpha = 0.5) +
  geom_errorbar(data = li_gdd_summary,
                aes(x = cumGDD_DAE, ymin = mean_LI - se_LI, ymax = mean_LI + se_LI),
                width = max_gdd * 0.01, alpha = 0.5, inherit.aes = FALSE) +
  # DAP annotation under the GDD breaks
  geom_text(
    data = dap_axis_gdd,
    aes(x = cumGDD_DAE, y = 0, label = paste0("(", DAP, ")")),
    inherit.aes = FALSE,
    vjust = 3.2, size = floor(fontsize / ggplot2::.pt) - 1
  ) +
  facet_grid(. ~ Season) +
  coord_cartesian(ylim = c(0, 1.1), clip = "off") +
  scale_color_manual(name = "Treatment", values = colors_temp, breaks = legend_order) +
  scale_fill_manual(values = colors_temp, breaks = legend_order, guide = "none") +
  scale_linetype_manual(name = "Treatment", values = linetype_temp, breaks = legend_order) +
  scale_x_continuous(
    name   = "Growing degree days from emergence (GDD, Tbase = 2°C)\n(Days after planting)",
    limits = c(0, max_gdd + 50),
    breaks = gdd_axis_breaks
  ) +
  scale_y_continuous(name = "Fraction of light intercepted",
                     breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  theme(
    legend.position    = "top",
    legend.key.size    = unit(15, "mm"),
    legend.background  = element_blank(),
    legend.title       = element_text(margin = margin(r = 6, unit = "mm")),
    panel.grid         = element_blank(),
    panel.spacing.x    = unit(1, "mm"),
    plot.margin        = margin(-5, 1, 5, 1, "mm"),
    axis.title.x       = element_text(vjust = -3)
  )
save_plot(fig_fPAR_GDD, "fig_fPAR_GDD", height = 5)

# ==============================================================================
# PART 7: figS_LI_GDD_24 / _25 — per-plot fitted vs observed
# ==============================================================================
li_gdd_r2 <- li_gdd_fits[, .(Season, Treatment, Block, PlotID, r2)
][, ':='(label = sprintf("R² = %.2f", r2), x_pos = 150, y_pos = 0.1)]

figS_LI_GDD_24 <- copy(predictions_gdd[Season == 2024])[cumGDD_DAE >= te, prediction := NA] |>
  ggplot(aes(cumGDD_DAE, prediction, color = Treatment)) +
  geom_rect(
    data = rect_gdd[Season == 2024],
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15, inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_line() +
  geom_point(data = li_gdd_sum[Season == "2024"], aes(cumGDD_DAE, LI)) +
  geom_text(
    data = li_gdd_r2[Season == "2024"],
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE, hjust = 0, size = floor(fontsize / ggplot2::.pt) - 1
  ) +
  scale_color_manual(values = colors_temp, breaks = legend_order) +
  scale_fill_manual(values = colors_temp, breaks = legend_order, guide = "none") +
  scale_x_continuous(
    name   = "Growing degree days from emergence (GDD, Tbase = 2°C)",
    limits = c(0, max_gdd), breaks = seq(0, max_gdd, 400)
  ) +
  scale_y_continuous(
    name   = "Fraction of light intercepted",
    limits = c(0, 1), breaks = seq(0, 1, 0.2)
  ) +
  facet_wrap(~PlotID, ncol = 3) +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  theme(
    legend.position   = "top",
    legend.key.size   = unit(15, "mm"),
    legend.background = element_blank(),
    panel.grid        = element_blank(),
    panel.spacing.x   = unit(1, "mm"),
    plot.margin       = margin(1, 1, 0, 1, "mm")
  )
save_plot(figS_LI_GDD_24, "figS_LI_GDD_24", height = 6, width = 5)

figS_LI_GDD_25 <- copy(predictions_gdd[Season == 2025])[cumGDD_DAE >= te, prediction := NA] |>
  ggplot(aes(cumGDD_DAE, prediction, color = Treatment)) +
  geom_rect(
    data = rect_gdd[Season == 2025],
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15, inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_line() +
  geom_point(data = li_gdd_sum[Season == "2025"], aes(cumGDD_DAE, LI)) +
  geom_text(
    data = li_gdd_r2[Season == "2025"],
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE, hjust = 0, size = floor(fontsize / ggplot2::.pt) - 1
  ) +
  scale_color_manual(values = colors_temp, breaks = legend_order) +
  scale_fill_manual(values = colors_temp, breaks = legend_order, guide = "none") +
  scale_x_continuous(
    name   = "Growing degree days from emergence (GDD, Tbase = 2°C)",
    limits = c(0, max_gdd), breaks = seq(0, max_gdd, 400)
  ) +
  scale_y_continuous(
    name   = "Fraction of light intercepted",
    limits = c(0, 1), breaks = seq(0, 1, 0.2)
  ) +
  facet_wrap(~PlotID, ncol = 3) +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  theme(
    legend.position   = "top",
    legend.key.size   = unit(15, "mm"),
    legend.background = element_blank(),
    panel.grid        = element_blank(),
    panel.spacing.x   = unit(1, "mm"),
    plot.margin       = margin(0, 1, 0, 1, "mm")
  )
save_plot(figS_LI_GDD_25, "figS_LI_GDD_25", height = 6, width = 5)

# ==============================================================================
# PART 8: Total intercepted radiation per treatment
# Modelled daily LI × observed daily IRRAD, summed from emergence to fitted te.
# ==============================================================================
daily_rad_gdd <- daily_gdd_join[!is.na(IRRAD)]

integrate_radiation <- function(params) {
  rbindlist(lapply(seq_len(nrow(params)), function(i) {
    p  <- params[i]
    dg <- daily_rad_gdd[Season == p$Season & cumGDD_DAE <= p$te]
    li <- pmax(0, combined_equation(dg$cumGDD_DAE, p$tm1, p$t1, p$t2, p$te, p$vmax))
    cbind(p[, .SD, .SDcols = intersect(c("Season", "Treatment", "Block", "PlotID"), names(p))],
          data.table(total_MJ_m2 = sum(dg$IRRAD * li, na.rm = TRUE), n_days = nrow(dg)))
  }))
}

trt_params_int <- copy(li_gdd_trt_params)[
  , .(Season = as.integer(as.character(Season)), Treatment, tm1, t1, t2, te, vmax)
]
total_rad_intcp <- integrate_radiation(trt_params_int)
total_rad_intcp[, total_MJ_m2 := round(total_MJ_m2)]
print(total_rad_intcp)
fwrite(total_rad_intcp, "results/LI_total_intercepted_radiation.csv", bom = TRUE)

# Per-plot totals for the significance test
plot_params_int <- copy(li_gdd_params)[
  , .(Season = as.integer(as.character(Season)), Treatment, Block, PlotID,
      tm1, t1, t2, te, vmax)
]
plot_rad_intcp <- integrate_radiation(plot_params_int)
plot_rad_intcp[, ':='(
  Season    = factor(Season),
  Treatment = factor(Treatment, levels = legend_order),
  Block     = factor(Block)
)]

for (s in c("2024", "2025")) {
  cat("\n===", s, "===\n")
  L <- fit_treatment_lm(plot_rad_intcp[Season == s], "total_MJ_m2")
  print(L$anova)
  cat("Shapiro-Wilk on residuals:\n")
  print(shapiro.test(residuals(L$model)))
  cat("Pairwise (Tukey):\n")
  print(pairs(L$emmeans, adjust = "tukey"))
}

plot_rad_intcp[, .(mean_total = mean(total_MJ_m2),
                   se_total   = sd(total_MJ_m2) / sqrt(.N)),
               by = .(Season, Treatment)
][, formatted := sprintf("%.1f ± %.1f", mean_total, se_total)] |>
  fwrite("results/LI_total_intercepted_radiation_summary.csv", bom = TRUE)

# ==============================================================================
# PART 9: t2 and te in GDD and in calendar days (AC vs HTI)
# Reports how many calendar days earlier the HTI canopy senesced.
# ==============================================================================
t2_te_gdd_summary <- li_gdd_params[
  Treatment %in% c("AC", "HTI"),
  .(mean_t2 = mean(t2), se_t2 = sd(t2) / sqrt(.N),
    mean_te = mean(te), se_te = sd(te) / sqrt(.N)),
  by = .(Season = as.integer(as.character(Season)), Treatment)
]
t2_te_gdd_wide <- dcast(t2_te_gdd_summary, Season ~ Treatment,
                        value.var = c("mean_t2", "se_t2", "mean_te", "se_te"))
t2_te_gdd_wide[, ':='(
  gdd_diff_t2    = mean_t2_AC - mean_t2_HTI,
  gdd_diff_t2_se = sqrt(se_t2_AC^2 + se_t2_HTI^2),
  gdd_diff_te    = mean_te_AC - mean_te_HTI,
  gdd_diff_te_se = sqrt(se_te_AC^2 + se_te_HTI^2)
)]
cat("\n=== t2 and te GDD: mean ± SE per treatment and AC - HTI differences ===\n")
print(t2_te_gdd_wide[, .(Season, mean_t2_AC, se_t2_AC, mean_t2_HTI, se_t2_HTI,
                         mean_te_AC, se_te_AC, mean_te_HTI, se_te_HTI,
                         gdd_diff_t2, gdd_diff_t2_se, gdd_diff_te, gdd_diff_te_se)])
