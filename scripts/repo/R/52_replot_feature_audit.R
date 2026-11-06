# ================================================================
# 52_replot_feature_audit.R
# ----------------------------------------------------------------
# Re-draws two figures of the feature audit (03) at paper size from
# its CSV outputs, so the slow audit does not have to be re-run:
#
#   P5_RESULTS/figures/03_feature_audit/feature_set_comparison.{png,pdf}
#   plots/drop_cluster_importance_multimetric.{png,pdf}
#   plots/drop_cluster_importance_mining_f1.{png,pdf}
#
# Inputs (PATHS$out_validation/tables):
#   feature_set_cv_summary.csv, drop_cluster_importance.csv
# ================================================================

source("R/00_config.R")
pacman::p_load(dplyr, tidyr, readr, ggplot2)

out_dir <- PATHS$out_validation
tab_dir <- file.path(out_dir, "tables")
plt_dir <- results_fig_dir("03_feature_audit")   # same folder as 03

BASE <- 16    # base text size for both figures

SCOPE_LABELS <- c(national = "National", ZONE_LABELS)
scope_fac <- function(s) factor(s, levels = names(SCOPE_LABELS),
                                labels = unname(SCOPE_LABELS))

theme_pub <- function(base = BASE) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
          strip.background = element_rect(fill = "grey96", colour = NA),
          strip.text       = element_text(face = "bold", size = base),
          axis.text        = element_text(size = base - 2, colour = "grey15"),
          axis.title       = element_text(size = base),
          legend.text      = element_text(size = base - 1),
          legend.title     = element_text(size = base),
          plot.background  = element_rect(fill = "white", colour = NA))
}
png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
save_both <- function(fig, path, width, height) {
  fig <- fig + labs(title = NULL, subtitle = NULL)   # captions go in LaTeX
  ggsave(path, fig, device = png_dev, width = width, height = height,
         units = "in", dpi = 600, bg = "white")
  ggsave(sub("\\.png$", ".pdf", path), fig, device = pdf_dev,
         width = width, height = height, units = "in", bg = "white")
  cat("  ->", path, "(+ .pdf)\n")
}

# ================================================================
# 1. FEATURE-SET COMPARISON
# ================================================================
cat("[1] Feature-set comparison\n")
set_summary <- read_csv(file.path(tab_dir, "feature_set_cv_summary.csv"),
                        show_col_types = FALSE)

cfg_levels <- c("all_candidates", "literature_union", "top_k", "decorrelated",
                "nested_top_k", "nested_decorrelated")
cfg_labels <- c(all_candidates      = "All 27 bands",
                literature_union    = "Literature bands",
                top_k               = "Top-10",
                decorrelated        = "Decorrelated\ntop-10",
                nested_top_k        = "Top-10\n(nested)",
                nested_decorrelated = "Decorrelated\ntop-10 (nested)")

cmp_plot <- bind_rows(
  set_summary |> transmute(config, metric = "Accuracy", mean = mean_accuracy,
                           lo = ci95_lower, hi = ci95_upper, selection_on),
  set_summary |> transmute(config, metric = "Mining F1", mean = mean_mining_f1,
                           lo = pmax(mean_mining_f1 - sd_mining_f1, 0),
                           hi = pmin(mean_mining_f1 + sd_mining_f1, 1), selection_on)
) |>
  mutate(config = factor(config, levels = cfg_levels))

fig_cfg <- ggplot(cmp_plot, aes(config, mean, fill = selection_on)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", mean), y = hi), vjust = -0.5, size = 4.8) +
  facet_wrap(~ metric, nrow = 1) +
  scale_x_discrete(labels = cfg_labels) +
  scale_fill_manual(values = c("fixed a priori"        = "#BBBBBB",
                               "all data (optimistic)" = "#CCBB44",
                               "training folds only"   = "#4477AA"),
                    name = "Bands selected on") +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = NULL, y = "Spatial CV score\n(accuracy: mean \u00b1 95% CI; F1: mean \u00b1 SD)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = BASE - 2),
        legend.position = "bottom",
        legend.key.size = unit(0.7, "cm"))
save_both(fig_cfg, file.path(plt_dir, "feature_set_comparison.png"), 13, 7.5)

# ================================================================
# 2. DROP-FAMILY IMPORTANCE
# ================================================================
cat("[2] Drop-family importance\n")
cluster_importance <- read_csv(file.path(tab_dir, "drop_cluster_importance.csv"),
                               show_col_types = FALSE)

# same family labels as in 03
make_family_label <- function(features_str) {
  m <- sort(strsplit(features_str, ",\\s*")[[1]])
  n_m <- length(m)
  if (n_m == 1) {
    pretty1 <- c(IRON_OXIDE = "Iron oxide", WATER_PROB = "Water probability",
                 NDVI_ANOMALY = "NDVI anomaly", NIR_corr = "NIR correlation",
                 NIR_STD_DEV = "NIR std dev", NIR_contrast = "NIR contrast",
                 NIR_var = "NIR variance", ELEVATION = "Elevation",
                 ASPECT = "Aspect", SLOPE = "Slope", TRI = "TRI",
                 CLAY_INDEX = "Clay index", BSI = "BSI", MNDWI = "MNDWI",
                 NDBI = "NDBI", ENTROPY = "Entropy")
    return(unname(if (m %in% names(pretty1)) pretty1[m] else m))
  }
  exact_themes <- list(
    "Bare soil (BSI+CLAY)"   = c("BSI", "CLAY_INDEX"),
    "Topography (slope+TRI)" = c("SLOPE", "TRI"),
    "SWIR + built-up"        = c("SWIR1", "SWIR2", "AWEI_NSH", "NDBI"))
  for (lab in names(exact_themes)) if (setequal(m, exact_themes[[lab]])) return(lab)
  subset_themes <- list(
    "NIR texture" = c("NIR_STD_DEV", "NIR_contrast", "NIR_var", "NIR_corr"))
  for (lab in names(subset_themes)) if (all(m %in% subset_themes[[lab]])) return(lab)
  tag_rules <- list(Visible = c("BLUE", "GREEN", "RED"), NIR = "NIR",
                    SWIR = c("SWIR1", "SWIR2"),
                    greenness = c("NDVI", "NDWI", "NDVI_ANOMALY"),
                    `water idx.` = c("MNDWI", "AWEI_NSH", "AWEI_SH", "WATER_PROB"),
                    color = c("CI", "BI"))
  tags <- names(tag_rules)[vapply(tag_rules, function(t) any(m %in% t), logical(1))]
  if (length(tags)) return(sprintf("%s (n=%d)", paste(tags, collapse = " + "), n_m))
  if (n_m <= 3) return(paste(m, collapse = " + "))
  sprintf("%s + %d more", paste(m[1:2], collapse = " + "), n_m - 2)
}

plot_data <- cluster_importance |>
  mutate(family_label = vapply(features, make_family_label, character(1))) |>
  filter(!is.na(delta_mining_f1)) |>
  group_by(zone) |>
  mutate(family_label = make.unique(family_label, sep = " #")) |>
  ungroup() |>
  arrange(match(zone, names(SCOPE_LABELS)), delta_mining_f1) |>
  mutate(key = factor(paste(zone, family_label, sep = " | "),
                      levels = unique(paste(zone, family_label, sep = " | "))))
strip_scope <- function(b) sub(".*\\| ", "", b)
n_scopes <- length(unique(plot_data$zone))

theme_bars <- function(base = BASE) {
  theme_pub(base) +
    theme(panel.grid.major.y = element_blank(),
          strip.text.y = element_text(angle = 0, hjust = 0, size = base),
          axis.text.y  = element_text(size = base - 3),
          panel.spacing = unit(0.7, "lines"))
}

# (a) multimetric
plot_multi <- plot_data |>
  select(zone, key, delta_accuracy, delta_mining_f1, delta_log_loss) |>
  pivot_longer(starts_with("delta_"), names_to = "metric", values_to = "delta") |>
  mutate(metric = factor(recode(metric, delta_accuracy = "Accuracy",
                                delta_mining_f1 = "Mining F1",
                                delta_log_loss = "Log-loss"),
                         levels = c("Accuracy", "Mining F1", "Log-loss")))
fig_multi <- ggplot(plot_multi, aes(y = key, x = delta, fill = delta > 0)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey20") +
  geom_col(colour = "grey15", linewidth = 0.2, width = 0.75) +
  facet_grid(rows = vars(scope_fac(zone)), cols = vars(metric),
             scales = "free", space = "free_y") +
  scale_y_discrete(labels = strip_scope) +
  scale_x_continuous(n.breaks = 4, expand = expansion(mult = 0.08)) +
  scale_fill_manual(values = c(`TRUE` = "#4477AA", `FALSE` = "#BBBBBB"), guide = "none") +
  labs(y = NULL,
       x = "\u0394 metric (baseline \u2212 without family; log-loss: without \u2212 baseline)") +
  theme_bars()
save_both(fig_multi, file.path(plt_dir, "drop_cluster_importance_multimetric.png"),
          15, max(8, 0.36 * nrow(plot_data) + 0.6 * n_scopes + 2))

# (b) mining F1 only
lim <- max(abs(plot_data$delta_mining_f1), 1e-3, na.rm = TRUE) * 1.35
fig_mF1 <- ggplot(plot_data, aes(y = key, x = delta_mining_f1,
                                 fill = delta_mining_f1 > 0)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey20") +
  geom_col(colour = "grey15", linewidth = 0.25, width = 0.75) +
  geom_text(aes(label = sprintf("%+.3f", delta_mining_f1),
                hjust = ifelse(delta_mining_f1 >= 0, -0.15, 1.15)),
            size = 4, colour = "grey15") +
  facet_grid(rows = vars(scope_fac(zone)), scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = strip_scope) +
  scale_x_continuous(limits = c(-lim, lim), expand = expansion(mult = 0.02)) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#BBBBBB"), guide = "none") +
  labs(y = NULL, x = "\u0394 Mining F1 (baseline \u2212 without family)") +
  theme_bars()
save_both(fig_mF1, file.path(plt_dir, "drop_cluster_importance_mining_f1.png"),
          11, max(7, 0.36 * nrow(plot_data) + 0.6 * n_scopes + 2))