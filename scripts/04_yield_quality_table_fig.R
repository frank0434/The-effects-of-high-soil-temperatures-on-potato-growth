# ==============================================================================
# Final yield and quality #
#
# Inputs  data/tuber_yield_final_harvest.csv   per-plot fresh yield at final harvest
#         data/tuber_avt_combined.csv          per-plot tuber count/weight
#         data/tuber_w_pct_combined.csv        per-plot marketable/secondary %
#         data/dmc_final_harvest.csv           per-plot dry matter content
#         data/weibull_fit_results_2025.csv    per-plot Weibull shape/scale
#         data/frying_colour_data.csv          per-plot frying colour index
# ==============================================================================
source("scripts/00_setup.R")
source("scripts/functions.R")
# ==============================================================================
# HELPERS (mean +/- SE formatting and Tukey CLD significance letters)
# ==============================================================================

fmt_ms <- function(m, s, ltr = "", dm = 0, ds = 1) {
  if (is.na(m) || is.na(s)) return(NA_character_)
  sprintf(paste0("%.", dm, "f \u00b1 %.", ds, "f%s"), m, s, ltr)
}

# Fit lm(var ~ Treatment [+ Block]), return cld letter DT (empty letters if NS)
fit_letters <- function(dt, var_col) {
  sub <- dt[!is.na(get(var_col))]
  sub[, Treatment := droplevels(factor(Treatment, levels = c("AC", "HE", "HTI", "HB")))]
  sub[, Block := factor(Block)]
  empty <- data.table(Treatment = character(), letter = character())
  if (nrow(sub) < 4 || length(unique(sub$Treatment)) < 2) return(empty)

  has_blk <- length(unique(sub$Block)) > 1
  fml <- paste(var_col, if (has_blk) "~ Treatment + Block" else "~ Treatment")
  mod <- tryCatch(lm(as.formula(fml), data = sub), error = function(e) NULL)
  if (is.null(mod)) return(empty)

  a_res <- tryCatch(car::Anova(mod, type = "II"), error = function(e) NULL)
  trt_p <- if (!is.null(a_res)) a_res[["Pr(>F)"]][1] else NA_real_
  if (is.na(trt_p) || trt_p >= 0.05) {
    return(data.table(Treatment = as.character(unique(sub$Treatment)), letter = ""))
  }

  emm <- tryCatch(emmeans(mod, ~ Treatment), error = function(e) NULL)
  if (is.null(emm)) return(data.table(Treatment = as.character(unique(sub$Treatment)), letter = ""))
  cld_res <- tryCatch(
    as.data.table(multcomp::cld(emm, Letters = letters, adjust = "tukey", quiet = TRUE)),
    error = function(e) NULL)
  if (is.null(cld_res))
    return(data.table(Treatment = as.character(unique(sub$Treatment)), letter = ""))
  cld_res[, .(Treatment = as.character(Treatment), letter = trimws(.group))]
}

# Summarise var per season: raw per-plot mean +/- SE + significance letters
summarise_var <- function(dt, var_col, seasons, excl_plot_2025 = 8L, dm = 0, ds = 1) {
  rbindlist(lapply(seasons, function(s) {
    sub <- copy(dt[Season == s & !is.na(get(var_col))])
    if (s == 2025L && !is.null(excl_plot_2025))
      sub <- sub[as.integer(PlotID) != excl_plot_2025]
    if (nrow(sub) == 0) return(NULL)
    ltr <- fit_letters(sub, var_col)
    st  <- sub[, .(Season = s,
                    mean = mean(get(var_col), na.rm = TRUE),
                    se   = sd(get(var_col),  na.rm = TRUE) / sqrt(.N)),
               by = Treatment]
    st[ltr, letter := i.letter, on = "Treatment"]
    st[is.na(letter), letter := ""]
    st[, formatted := mapply(fmt_ms, mean, se, letter,
                              MoreArgs = list(dm = dm, ds = ds))]
  }))
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

yield_dt <- fread("data/tuber_yield_final_harvest.csv")
yield_dt[, `:=`(Season = as.integer(Season), Block = factor(Block), PlotID = as.integer(PlotID))]
yield_dt[, fresh_yield_kg_m2 := tuber_fw_kg_ha / 10000]

avt_dt <- fread("data/tuber_avt_combined.csv")
avt_dt[, avt := weight_of_tubers / number_of_tubers]

avt_dt[, tuber_n_per_m2 := number_of_tubers / (number_of_plants * 0.75 * 0.3)]
avt_dt[, `:=`(Season = as.integer(Season), Block = factor(Block), PlotID = as.integer(PlotID))]

mkt_dt <- fread("data/tuber_w_pct_combined.csv")
mkt_dt[, `:=`(Season = as.integer(Season), Block = factor(Block), PlotID = as.integer(PlotID))]

dmc_dt <- fread("data/dmc_final_harvest.csv")
dmc_dt[, tuber_dmc := tuber_dmc * 100]
dmc_dt[, `:=`(Season = as.integer(Season), Block = factor(Block), PlotID = as.integer(PlotID))]

weib_dt <- fread("data/weibull_fit_results_2025.csv")
weib_dt[, `:=`(Season = 2025L, Block = factor(Block), PlotID = as.integer(PlotID))]
weib_dt <- weib_dt[!is.na(scale)]

fry_dt <- fread("data/frying_colour_data.csv")
fry_dt[, `:=`(Season = as.integer(Season), Block = factor(Block), PlotID = as.integer(PlotID))]
fry_dt <- fry_dt[!is.na(Fry.index)]

# ==============================================================================
# TABLE 5: final yield and quality
# ==============================================================================

yield_s   <- summarise_var(yield_dt, "fresh_yield_kg_m2", c(2024L, 2025L), dm = 2, ds = 3)
tuber_n_s <- summarise_var(avt_dt,   "tuber_n_per_m2",    c(2024L, 2025L), dm = 0, ds = 0)
avt_s     <- summarise_var(avt_dt,   "avt",               c(2024L, 2025L), dm = 0, ds = 1)
dmc_s     <- summarise_var(dmc_dt,   "tuber_dmc",         c(2024L, 2025L), dm = 0, ds = 1)
fry_s     <- summarise_var(fry_dt,   "Fry.index",         c(2024L, 2025L), dm = 1, ds = 2)

get_f <- function(src, seas, trt) {
  v <- src[Season == seas & Treatment == trt, formatted]
  if (length(v) == 0 || all(is.na(v))) return("\u2014")
  v[1]
}

t5_grid <- rbindlist(list(
  data.table(Season = 2024L, Treatment = c("AC", "HE", "HTI")),
  data.table(Season = 2025L, Treatment = c("AC", "HTI", "HB"))
))

table5 <- t5_grid[, .(
  Season, Treatment,
  `Fresh tuber yield (kg m-2)` = get_f(yield_s,   Season, Treatment),
  `Tuber number per m-2`       = get_f(tuber_n_s, Season, Treatment),
  `Average tuber weight (g)`   = get_f(avt_s,     Season, Treatment),
  `Dry matter content (%)`     = get_f(dmc_s,     Season, Treatment),
  `Frying colour index`        = get_f(fry_s,     Season, Treatment)
), by = seq_len(nrow(t5_grid))][, -1]

fwrite(table5, "results/table5_yield_quality.csv", bom = TRUE)

# ==============================================================================
# FIGURE 7: marketable yield, deformed tubers, Weibull scale
# ==============================================================================

# --- Weibull scale model + emmeans CI (2025 only) ---
scale_model_2025 <- lm(scale ~ Treatment, data = weib_dt)
scale_emm <- emmeans(scale_model_2025, ~ Treatment)
scale_means <- as.data.table(summary(scale_emm, infer = c(TRUE, TRUE)))
scale_means[, `:=`(Treatment = factor(Treatment, levels = legend_order),
                   Season = factor("2025", levels = c("2024", "2025")))]
scale_letters <- fit_letters(weib_dt, "scale")
scale_letters[, Treatment := factor(Treatment, levels = legend_order)]
scale_letters <- scale_letters[scale_means, on = "Treatment"]

# --- Market %/deformed % (both seasons; exclude 2025 plot 8, as in the models above) ---
fig_pct <- droplevels(mkt_dt[!(Season == 2025 & PlotID == 8)])
fig_pct[, Treatment := factor(Treatment, levels = legend_order)]
fig_pct[, Season := factor(Season)]

market_letters <- rbindlist(lapply(c("2024", "2025"), function(s) {
  fit_letters(fig_pct[Season == s], "market_pct")[, Season := s]
}))
sec_letters <- rbindlist(lapply(c("2024", "2025"), function(s) {
  fit_letters(fig_pct[Season == s], "Secondary_pct")[, Season := s]
}))
market_letters[, Treatment := factor(Treatment, levels = legend_order)]
sec_letters[, Treatment := factor(Treatment, levels = legend_order)]

# Shared y-axis range for Market and Secondary (both are %, directly comparable)
shared_pct_ymax <- max(fig_pct$market_pct, fig_pct$Secondary_pct, na.rm = TRUE) * 1.15

# --- Panel A: Marketable yield (%) ---
fig_market <- fig_pct |>
  ggplot(aes(x = Treatment, y = market_pct)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot(width = 0.5, outlier.shape = 16) +
  stat_summary(fun = "mean", geom = "point", color = "red", size = ps) +
  geom_text(data = market_letters, aes(x = Treatment, y = shared_pct_ymax * 0.92, label = letter),
            size = 5, inherit.aes = FALSE, family = "Times New Roman") +
  facet_wrap(~ Season, scales = "free_x") +
  scale_y_continuous(name = "Marketable yield \n(%)", limits = c(0, shared_pct_ymax), expand = c(0, 0)) +
  labs(x = NULL, tag = "A") +
  theme_bw(base_size = fontsize) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        plot.tag = element_text(size = 14),
        plot.tag.position = c(0, 1),
        text = element_text(size = fontsize, family = "Times New Roman"))

# --- Panel B: Deformed tubers (%) ---
fig_sec <- fig_pct |>
  ggplot(aes(x = Treatment, y = Secondary_pct)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot(width = 0.5, outlier.shape = 16) +
  stat_summary(fun = "mean", geom = "point", color = "red", size = ps) +
  geom_text(data = sec_letters, aes(x = Treatment, y = shared_pct_ymax * 0.92, label = letter),
            size = 5, inherit.aes = FALSE, family = "Times New Roman") +
  facet_wrap(~ Season, scales = "free_x") +
  scale_y_continuous(name = "Deformed \ntuber yield (%)", limits = c(0, shared_pct_ymax), expand = c(0, 0)) +
  labs(x = NULL, tag = "B") +
  theme_bw(base_size = fontsize) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        text = element_text(size = fontsize, family = "Times New Roman"),
        plot.tag = element_text(size = 14),
        plot.tag.position = c(0, 1))

# --- Panel C: Weibull scale — 2025 only; 2024 facet left empty ("Not applicable") ---
weibull_blank <- data.table(Season = factor("2024", levels = c("2024", "2025")),
                            Treatment = factor(c("AC", "HE", "HTI"), levels = legend_order),
                            emmean = NA_real_)
na_note <- data.table(Season = factor("2024", levels = c("2024", "2025")),
                      Treatment = factor("HE", levels = legend_order),
                      emmean = mean(scale_means$emmean))

fig_weibull <- ggplot(scale_means, aes(x = Treatment, y = emmean)) +
  geom_blank(data = weibull_blank) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, linewidth = lw) +
  geom_text(data = scale_letters, aes(x = Treatment, y = upper.CL + 1, label = letter),
            size = 5, inherit.aes = FALSE, family = "Times New Roman") +
  geom_text(data = na_note, aes(label = "Not applicable"),
            size = 5, fontface = "italic", color = "grey40", inherit.aes = TRUE) +
  facet_wrap(~ Season, scales = "free_x") +
  labs(x = NULL, tag = "C") +
  ylab("Weibull scale \n(mm)") +
  theme_bw(base_size = fontsize) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        text = element_text(size = fontsize, family = "Times New Roman"),
        plot.tag = element_text(size = 14),
        plot.tag.position = c(0, 1))

# --- Compose: stacked vertically, season shown only on top panel ---
fig7 <- (fig_market / plot_spacer() / fig_sec / plot_spacer() / fig_weibull) +
  plot_layout(heights = c(1, -0.1, 1, -0.1, 1)) &
  theme(legend.position = "none")

print(fig7)
