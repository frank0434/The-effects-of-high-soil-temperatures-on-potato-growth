# ==============================================================================
# Constants and functions for 02.4_light_GDD_clean.R
#
# Self-contained: this is the only file that script sources. Extracted from
# 00_setup.R, 00_equations_and_models.R, functions_nlme.R,
# functions_logistic_fitting.R and functions_visualise.R.
# ==============================================================================
library(data.table)
library(ggplot2)
library(emmeans)   # emmeans(), pairs() methods
library(DEoptim)   # global search for curve starting values
# also called via ::  — car, multcomp, minpack.lm

# ==============================================================================
# EXPERIMENTAL CONSTANTS
# ==============================================================================
planting_dates     <- list(`2024` = as.Date("2024-05-13"),
                           `2025` = as.Date("2025-05-07"))
emergence_date     <- list(`2024` = as.Date("2024-06-07"),
                           `2025` = as.Date("2025-05-23"))
leaf_removal_dates <- list(`2024` = as.Date("2024-09-02"),
                           `2025` = as.Date("2025-08-18"))

# Heating windows, converted to days after planting for the plot rectangles
treatment_rect_dt <- data.table(
  Season    = c(2024L, 2024L, 2025L, 2025L),
  Treatment = c("HE", "HTI", "HTI", "HB"),
  start     = as.Date(c("2024-06-07", "2024-06-28", "2025-06-13", "2025-07-12")),
  end       = as.Date(c("2024-06-21", "2024-07-15", "2025-06-29", "2025-07-28"))
)
treatment_rect_dt[, ':='(
  DAP_min = as.integer(start - planting_dates[[as.character(Season)]]),
  DAP_max = as.integer(end   - planting_dates[[as.character(Season)]])
), by = Season]
treatment_rect_dt[, c("start", "end") := NULL]

# ==============================================================================
# PLOTTING THEME
# ==============================================================================
legend_order  <- c("AC", "HE", "HTI", "HB")
colors_temp   <- c("AC" = "#000000", "HE" = "#314856",
                   "HTI" = "#8B1C1C", "HB" = "#fd8700")
linetype_temp <- c("AC" = "solid", "HE" = "dashed",
                   "HTI" = "solid", "HB" = "dotdash")
lw       <- 0.8
fontsize <- 12

save_plot <- function(plot, name, width = 7, height = 5, dpi = 300, format = "svg") {
  filename <- paste0("figures/", name, ".", format)
  ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
}

# ==============================================================================
# CANOPY COVER EQUATION
# ==============================================================================

#' Combined canopy cover equation with plateau
#'
#' Three phases: logistic expansion (up to t1), plateau (t1 to t2) and
#' senescence decay (t2 to te).
#'
#' @param td_sum numeric vector - cumulative thermal time
#' @param tm1 numeric - thermal time at the inflection of canopy expansion
#' @param t1 numeric - thermal time at the end of canopy expansion
#' @param t2 numeric - thermal time at the end of the plateau
#' @param te numeric - thermal time at crop termination
#' @param vmax numeric - maximum canopy cover / light interception
#'
#' @return numeric vector of predicted values
combined_equation <- function(td_sum, tm1, t1, t2, te, vmax) {
  ifelse(
    td_sum < t1,
    # Growth phase
    vmax * (1 + (t1 - td_sum) / (t1 - tm1)) * (td_sum / t1) ^ (t1 / (t1 - tm1)),
    ifelse(
      td_sum <= t2,
      # Plateau phase
      vmax,
      # Senescence phase
      vmax * ((te - td_sum) / (te - t2)) * ((td_sum + t1 - t2) / t1) ^ (t1 / (te - t2))
    )
  )
}

# ==============================================================================
# CURVE FITTING
# ==============================================================================

#' Find starting values for the combined equation with DEoptim
#'
#' The combined equation is badly behaved from arbitrary starts, so a global
#' search supplies nlsLM with a feasible starting point. Bounds keep the four
#' break points ordered and separated by at least min_gap.
#'
#' @param clean_data data.table with the response and time columns
#' @param y_col character - response column name
#' @param x_col character - time column name
#' @param tmax numeric - upper bound for te (cannot exceed leaf removal)
#' @param seed numeric - optional seed; derived from the data when NULL
#'
#' @return list with tm1, t1, t2, te, vmax; NULL if fewer than 5 observations
optimize_start_values_combined <- function(clean_data, y_col, x_col, tmax, seed = NULL) {

  required_cols <- c(y_col, x_col)
  if (!all(required_cols %in% names(clean_data))) {
    stop(sprintf("Missing required columns: %s",
                 paste(setdiff(required_cols, names(clean_data)), collapse = ", ")))
  }

  dat <- copy(clean_data)[!is.na(get(y_col)) & !is.na(get(x_col))]
  if (nrow(dat) < 5) {
    return(NULL)
  }

  setorderv(dat, x_col)

  if (is.null(seed)) {
    seed <- as.integer(abs(sum(dat[[y_col]], na.rm = TRUE) * 1000) %% 10000)
  }
  set.seed(seed)

  observed_vmax <- max(dat[[y_col]], na.rm = TRUE)
  if (!is.finite(observed_vmax)) {
    return(NULL)
  }

  x_min <- min(dat[[x_col]], na.rm = TRUE)
  x_max <- max(dat[[x_col]], na.rm = TRUE)
  x_range <- x_max - x_min
  if (!is.finite(x_range) || x_range <= 0) {
    return(NULL)
  }

  min_gap <- max(5, x_range * 0.1)

  tm1_lower <- max(0, x_min)
  tm1_upper <- max(tm1_lower + min_gap, x_min + x_range * 0.35)

  t1_lower <- tm1_lower + min_gap
  t1_upper <- max(t1_lower + min_gap, x_min + x_range * 0.6)

  t2_lower <- t1_lower + min_gap
  t2_upper <- max(t2_lower + min_gap, x_min + x_range * 0.9)

  te_lower <- t2_lower + min_gap
  te_upper <- tmax

  vmax_lower <- max(0.5 * observed_vmax, 0.5)
  vmax_upper <- max(vmax_lower + 1e-3, 1)

  rss_combined <- function(params, data, y_col, x_col, vmax_lower, vmax_upper, min_gap) {
    tm1  <- params[1]
    t1   <- params[2]
    t2   <- params[3]
    te   <- params[4]
    vmax <- params[5]

    if (!is.finite(tm1) || !is.finite(t1) || !is.finite(t2) || !is.finite(te) || !is.finite(vmax) ||
        t1 <= (tm1 + min_gap) || t2 <= (t1 + min_gap) || te <= (t2 + min_gap) ||
        vmax < vmax_lower || vmax > vmax_upper) {
      return(1e12)
    }

    pred <- tryCatch(
      combined_equation(data[[x_col]], tm1 = tm1, t1 = t1, t2 = t2, te = te, vmax = vmax),
      error = function(e) rep(NA_real_, nrow(data))
    )

    if (any(!is.finite(pred))) {
      return(1e12)
    }

    rss <- sum((data[[y_col]] - pred)^2, na.rm = TRUE)
    if (!is.finite(rss)) 1e12 else rss
  }

  de_result <- DEoptim(
    fn = rss_combined,
    lower = c(tm1_lower, t1_lower, t2_lower, te_lower, vmax_lower),
    upper = c(tm1_upper, t1_upper, t2_upper, te_upper, vmax_upper),
    control = DEoptim.control(itermax = 400, F = 0.8, CR = 0.9, trace = FALSE),
    data = dat,
    y_col = y_col,
    x_col = x_col,
    vmax_lower = vmax_lower,
    vmax_upper = vmax_upper,
    min_gap = min_gap
  )

  list(
    tm1  = unname(de_result$optim$bestmem[1]),
    t1   = unname(de_result$optim$bestmem[2]),
    t2   = unname(de_result$optim$bestmem[3]),
    te   = unname(de_result$optim$bestmem[4]),
    vmax = unname(de_result$optim$bestmem[5])
  )
}

#' Fit the combined equation to one plot or one treatment
#'
#' Bounded Levenberg-Marquardt fit. Bounds mirror those used for the starting
#' value search but with a smaller minimum gap, so nlsLM can still move.
#'
#' @param data data.table for a single fitting unit
#' @param y_col character - response column name
#' @param x_col character - time column name
#' @param te_cap numeric - upper bound for te (thermal time at leaf removal)
#' @param start_vals list from optimize_start_values_combined()
#'
#' @return nlsLM object, or NULL if the fit failed
fit_combined <- function(data, y_col, x_col, te_cap, start_vals) {
  vmax_obs <- max(data[[y_col]], na.rm = TRUE)
  x_min <- min(data[[x_col]], na.rm = TRUE)
  x_max <- max(data[[x_col]], na.rm = TRUE)
  x_range <- x_max - x_min
  min_gap <- max(0.5, x_range * 0.03)

  tm1_lower <- max(0, x_min)
  tm1_upper <- max(tm1_lower + min_gap, x_min + x_range * 0.35)

  t1_lower <- tm1_lower + min_gap
  t1_upper <- max(t1_lower + min_gap, x_min + x_range * 0.6)

  t2_lower <- t1_lower + min_gap
  t2_upper <- max(t2_lower + min_gap, x_min + x_range * 0.9)

  te_lower <- t2_lower + min_gap
  te_upper <- te_cap

  vmax_lower <- max(0.5 * vmax_obs, 0.5)
  vmax_upper <- max(vmax_lower + 1e-3, 1)

  model_formula <- stats::as.formula(
    sprintf("%s ~ combined_equation(%s, tm1, t1, t2, te, vmax)", y_col, x_col)
  )

  tryCatch(
    minpack.lm::nlsLM(
      formula = model_formula,
      data = data,
      start = start_vals,
      lower = c(tm1_lower, t1_lower, t2_lower, te_lower, vmax_lower),
      upper = c(tm1_upper, t1_upper, t2_upper, te_upper, vmax_upper),
      control = minpack.lm::nls.lm.control(maxiter = 200)
    ),
    error = function(e) NULL
  )
}

#' R2 of an nls fit
#'
#' @param fit nls/nlsLM object or NULL
#' @return numeric R2, NA if the fit is NULL or the response has no variance
r2_nls <- function(fit) {
  if (is.null(fit)) {
    return(NA_real_)
  }

  y    <- as.numeric(fit$m$lhs())
  yhat <- as.numeric(stats::fitted(fit))

  sst <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  sse <- sum((y - yhat)^2, na.rm = TRUE)

  if (!is.finite(sst) || sst == 0) {
    return(NA_real_)
  }

  1 - sse / sst
}

#' Pull the five combined-equation parameters out of a fit
#'
#' @param fit nls/nlsLM object or NULL
#' @return one-row data.table (all NA if the fit is NULL)
extract_combined_params <- function(fit) {
  if (is.null(fit)) {
    return(data.table(tm1 = NA_real_, t1 = NA_real_, t2 = NA_real_,
                      te = NA_real_, vmax = NA_real_))
  }

  est <- coef(fit)

  data.table(
    tm1  = unname(est[["tm1"]]),
    t1   = unname(est[["t1"]]),
    t2   = unname(est[["t2"]]),
    te   = unname(est[["te"]]),
    vmax = unname(est[["vmax"]])
  )
}

#' Regular prediction grid spanning 0 to end_time for every plot
#'
#' @param data data.table carrying the Season/Treatment/Block/PlotID factors
#' @param end_time numeric - upper end of the grid
#' @param step numeric - grid spacing
#'
#' @return data.table with one cumGDD_DAE column per plot
generate_time_grid <- function(data, end_time = 60, step = 0.01) {
  time_grid <- data[, .(cumGDD_DAE = seq(0, end_time, step)),
                    by = .(Season, Treatment, Block, PlotID)]
  time_grid[, `:=`(
    Treatment = factor(Treatment, levels = levels(Treatment)),
    Block     = factor(Block, levels = levels(Block)),
    PlotID    = factor(PlotID, levels = levels(PlotID)))]

  setorder(time_grid, Season, Block, PlotID, cumGDD_DAE)
  time_grid
}

# ==============================================================================
# SIGNIFICANCE TABLES
# ==============================================================================
significance_label <- function(p_value) {
  ifelse(
    is.na(p_value),
    NA_character_,
    ifelse(p_value < 0.001, "**", ifelse(p_value < 0.05, "*", "ns"))
  )
}

#' Tukey-adjusted pairwise contrasts as a tidy table
#'
#' @param emm emmGrid from emmeans()
#' @param season integer - season label to attach
#' @param model_name character - model label to attach
#' @param parameter character - response label to attach
#'
#' @return data.table of contrasts with a Significant column
summarise_pairwise_table <- function(emm, season, model_name, parameter) {
  out <- as.data.table(summary(pairs(emm, adjust = "tukey"), infer = c(TRUE, TRUE)))
  out[, `:=`(
    Season      = season,
    Model       = model_name,
    Parameter   = parameter,
    Significant = significance_label(p.value)
  )]
  out
}
