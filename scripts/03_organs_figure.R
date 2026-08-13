# ==============================================================================
# Dry matter accumulation by organ — timeseries figure (Fig. 3 / p_organs_anno)
#
# Figure  p_organs_anno   leaf / stem / tuber dry weight vs DAP, by season,
#                         with treatment-period shading and significance
#                         asterisks at DAP 56/83 (2025)
#
# Input   data/drymatter_organs_long.csv   tidy long-format dry matter data
#                                           (Season, Block, Treatment, PlotID,
#                                           Date, DAP, Organ, value [g/m^2])
# ==============================================================================
source("scripts/00_setup.R")
library(ggh4x)


dma_long <- fread("data/drymatter_organs_long.csv")

# ── Summary statistics (mean ± SE by Season / Treatment / Organ / DAP) ───────

dma_organs <- dma_long[Organ %in% c("leaf", "stem", "tuber", "dead_canopy", "leaf_stem")]
dma_summary <- dma_organs[, .(
  mean = mean(value, na.rm = TRUE),
  se   = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value)))
), by = .(Season, Treatment, Organ, DAP)]

# Insert NA row to break the HE line after its final measurement (DAP 42)
dma_summary_plot <- rbindlist(list(
  dma_summary,
  data.table(Season = 2024, Treatment = "HE", DAP = 63,
             Organ  = c("leaf", "stem", "tuber", "leaf_stem"),
             mean = NA_real_, se = NA_real_)
), use.names = TRUE, fill = TRUE)

dma_summary_plot[, Organ := sub("^(.)", "\\U\\1", Organ, perl = TRUE)]
dma_summary_plot[Organ == "Leaf_stem", Organ := "Leaf:Stem"]
dma_summary_plot[, Organ := factor(Organ, levels = c("Leaf", "Stem", "Leaf:Stem", "Tuber", "Dead_canopy"))]

# ── Figure ────────────────────────────────────────────────────────────────────

# Asterisks mark sampling points with at least one significant pairwise
# difference (Tukey HSD, alpha = 0.05); `pair` records which contrasts drive
# each asterisk — used verbatim in the figure caption.
significance_annotations_25 <- data.table(
  Organ  = c("Stem", "Leaf", "Tuber"),
  Season = 2025,
  x_pos  = c(56,     83,     83),
  y_pos  = c(150,    220,    1300),
  label  = "*",
  pair   = c("HTI > AC (p = 0.040) and HTI > HB (p = 0.036)",
             "AC > HB (p = 0.026)",
             "AC > HTI (p = 0.028) and AC > HB (p = 0.005)")
)
significance_annotations_25[, Organ := factor(Organ, levels = c("Leaf", "Stem", "Tuber"))]

panel_tags <- data.table(
  Organ  = factor(c("Leaf", "Leaf", "Stem", "Stem", "Tuber", "Tuber"),
                  levels = c("Leaf", "Stem", "Tuber")),
  Season = c(2024, 2025, 2024, 2025, 2024, 2025),
  tag    = c("a",  "b",  "c",  "d",  "e",  "f")
)

p_organs <- dma_summary_plot[Organ != "Leaf:Stem" & Organ != "Dead_canopy"] |>
  ggplot(aes(DAP, mean, linetype = Treatment, shape = Treatment)) +
  geom_errorbar(aes(x = DAP, y = mean, ymin = mean - se, ymax = mean + se),
                width = ps, inherit.aes = FALSE) +
  geom_rect(
    data = treatments_DAP_arrows[
      , .(.N, xmin = min(DAP), xmax = max(DAP)), by = .(Season, Treatment)],
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = colors_temp, guide = "none") +
  geom_line(aes(color = Treatment), linewidth = lw - 0.3) +
  geom_point(size = ps, alpha = 0.7) +
  scale_color_manual(values = colors_temp,     limits = legend_order) +
  scale_linetype_manual(values = linetype_temp, limits = legend_order) +
  scale_shape_manual(values = point_shape,     limits = legend_order) +
  facet_grid(Organ ~ Season, scales = "free_y") +
  scale_x_continuous(name = "Days after planting", expand = c(0, 0), limits = c(0, 120)) +
  scale_y_continuous(name = expression(Dry ~ weight ~ (g ~ m^{-2}))) +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  theme(
    legend.position  = "top",
    panel.spacing.y  = unit(0, "mm"),
    legend.key.width = grid::unit(15, "mm"),
    legend.title     = element_text(margin = margin(r = 15)),
    panel.grid       = element_blank(),
    panel.spacing.x  = unit(7, "mm")
  )

p_organs_anno <- p_organs +
  geom_text(data = panel_tags,
            aes(x = -Inf, y = Inf, label = tag),
            hjust = -0.5, vjust = 1.3, size = 5, fontface = "bold",
            inherit.aes = FALSE) +
  geom_text(data = significance_annotations_25,
            aes(x = x_pos, y = y_pos, label = label),
            size = 8, color = "black", inherit.aes = FALSE) +
  ggh4x::facetted_pos_scales(y = list(
    Organ == "Leaf"  ~ scale_y_continuous(limits = c(0, 250)),
    Organ == "Stem"  ~ scale_y_continuous(limits = c(0, 250)),
    Organ == "Tuber" ~ scale_y_continuous(limits = c(0, 1500))
  )) +
  ylab(expression(Dry ~ weight ~ (g ~ m^{-2})))

p_organs_anno
