# ================================================================
# 40_performance_evaluation.R
# ----------------------------------------------------------------
# Performance of the NATIONAL model (10) from its spatial-CV OOF
# predictions. Headline numbers are national; the same OOF predictions
# are also split by latitude zone for reporting (the model itself does
# not use zones).
#
# Outputs (PATHS$out_unconstrained):
#   metrics/tab_cv_summary.csv          national + per zone
#   metrics/tab_per_class_metrics.csv   precision / recall / F1
#   metrics/tab_auc.csv                 one-vs-rest AUC, DeLong 95% CI
#   metrics/tab_calibration.csv         reliability bins + ECE per class
#   metrics/tab_feature_importance.csv  permutation importance (CV fold models)
#   metrics/tab_power_mde.csv           minimum detectable effect
#   metrics/tab_power_mcnemar.csv       power of the 41 comparisons (if run)
# Figures (P5_RESULTS/figures/40_performance[/appendix]):
#   fig_roc_national.png
#   figures/main/fig_confusion_national.png (+ .pdf)
#   figures/main/fig_f1_bars.png
#   figures/main/fig_feature_importance.png
#   figures/main/fig_calibration.png
#   figures/appendix/fig_roc_grid_zones.png
#   figures/appendix/fig_confusion_zones.png (+ .pdf)
# ================================================================

source("R/00_config.R")
source("R/01_utils.R")

pacman::p_load(dplyr, tidyr, readr, ggplot2, patchwork, pROC, ranger, pwr)

# ----------------------------------------------------------------
# 1. LOAD
# ----------------------------------------------------------------

p1_dir           <- PATHS$out_unconstrained
metrics_dir      <- file.path(p1_dir, "metrics")
fig_main_dir     <- results_fig_dir("40_performance")
fig_appendix_dir <- results_fig_dir("40_performance/appendix")
for (d in c(metrics_dir, fig_main_dir, fig_appendix_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

cat("STEP 1 - Loading inputs\n")

national_model <- load_national_model()
oof_df <- readr::read_csv(file.path(metrics_dir, "oof_predictions.csv"),
                          show_col_types = FALSE)
stopifnot(all(c("zone", "fold", "y_true", "y_pred",
                "prob_water", "prob_mining", "prob_veg") %in% names(oof_df)))
prob_cols <- c("prob_water", "prob_mining", "prob_veg")

SCOPES       <- c("national", ZONE_ORDER)
SCOPE_LABELS <- c(national = "National", ZONE_LABELS)
scope_df     <- function(s) if (s == "national") oof_df else oof_df[oof_df$zone == s, ]
scope_fac    <- function(s) factor(SCOPE_LABELS[s], levels = unname(SCOPE_LABELS))

cat(sprintf("  OOF rows: %d; bands: %s\n", nrow(oof_df),
            paste(national_model$features, collapse = ", ")))
cat(sprintf("  CV buffer %.0f m, block %.0f m\n",
            national_model$buffer_used_m, national_model$block_size_m))

# ----------------------------------------------------------------
# 2. CV SUMMARY (national + per zone)
# ----------------------------------------------------------------

cat("\nSTEP 2 - CV summary\n")

cv_summary_tbl <- bind_rows(lapply(SCOPES, function(s) {
  d  <- scope_df(s)
  fm <- fold_metrics(d$y_true, d$y_pred, d$fold, as.matrix(d[, prob_cols]))
  sm <- summarise_fold_metrics(fm)
  if (is.null(sm)) return(NULL)
  data.frame(scope = s, n_folds = sm$n_folds, n_test_total = sm$n_test,
             accuracy = sm$mean_accuracy, acc_sd = sm$sd_accuracy,
             ci95_lo = sm$ci95_lower, ci95_hi = sm$ci95_upper,
             kappa = sm$mean_kappa, kappa_sd = sm$sd_kappa,
             mining_f1 = sm$mean_mining_f1, mining_f1_sd = sm$sd_mining_f1,
             log_loss = sm$mean_log_loss)
}))
write_csv(cv_summary_tbl, file.path(metrics_dir, "tab_cv_summary.csv"))
print(cv_summary_tbl, digits = 3, row.names = FALSE)

# ----------------------------------------------------------------
# 3. PER-CLASS PRECISION / RECALL / F1 (CV mean +/- SD and pooled)
# ----------------------------------------------------------------

cat("\nSTEP 3 - Per-class precision/recall/F1\n")

per_class_tbl <- bind_rows(lapply(SCOPES, function(s) {
  d <- scope_df(s)
  fold_prfs <- bind_rows(lapply(split(d, d$fold), function(fd) {
    if (nrow(fd) < 10) return(NULL)
    per_class_prf(fd$y_true, fd$y_pred)$prf |> mutate(fold = fd$fold[1])
  }))
  if (!nrow(fold_prfs)) return(NULL)
  cv_agg <- fold_prfs |>
    group_by(class_id, class) |>
    summarise(precision_cv_mean = mean(precision, na.rm = TRUE),
              precision_cv_sd   = sd(precision,   na.rm = TRUE),
              recall_cv_mean    = mean(recall,    na.rm = TRUE),
              recall_cv_sd      = sd(recall,      na.rm = TRUE),
              f1_cv_mean        = mean(f1,        na.rm = TRUE),
              f1_cv_sd          = sd(f1,          na.rm = TRUE),
              n_folds           = sum(!is.na(f1)),
              .groups = "drop")
  pooled <- per_class_prf(d$y_true, d$y_pred)$prf |>
    rename(precision_pooled = precision, recall_pooled = recall, f1_pooled = f1)
  cv_agg |> left_join(pooled, by = c("class_id", "class")) |>
    mutate(scope = s, n_true = as.integer(table(factor(d$y_true, 1:3)))) |>
    select(scope, class_id, class, n_true, everything())
}))
write_csv(per_class_tbl, file.path(metrics_dir, "tab_per_class_metrics.csv"))

# ----------------------------------------------------------------
# 4. ROC + AUC (one-vs-rest, DeLong 95% CI)
# ----------------------------------------------------------------

cat("\nSTEP 4 - ROC / AUC\n")

roc_scope <- function(s) {
  d <- scope_df(s)
  if (nrow(d) < 20) return(NULL)
  Filter(Negate(is.null), lapply(1:3, function(cls) {
    y_bin <- as.integer(d$y_true == cls)
    p_cls <- d[[prob_cols[cls]]]
    if (length(unique(y_bin)) < 2 || sd(p_cls, na.rm = TRUE) < 1e-8) return(NULL)
    roc_obj <- pROC::roc(response = y_bin, predictor = p_cls,
                         levels = c(0, 1), direction = "<", quiet = TRUE)
    ci <- tryCatch(as.numeric(pROC::ci.auc(roc_obj, method = "delong")),
                   error = function(e) c(NA, as.numeric(roc_obj$auc), NA))
    list(curve = data.frame(scope = s, class = CLASS_NAMES[cls],
                            fpr = 1 - roc_obj$specificities,
                            tpr = roc_obj$sensitivities),
         stat  = data.frame(scope = s, class = CLASS_NAMES[cls],
                            n_pos = sum(y_bin), auc = as.numeric(roc_obj$auc),
                            auc_ci_lo = ci[1], auc_ci_hi = ci[3]))
  }))
}
roc_nested <- lapply(SCOPES, roc_scope)
roc_curves <- bind_rows(lapply(roc_nested, function(z) bind_rows(lapply(z, `[[`, "curve"))))
auc_tbl    <- bind_rows(lapply(roc_nested, function(z) bind_rows(lapply(z, `[[`, "stat"))))
write_csv(auc_tbl, file.path(metrics_dir, "tab_auc.csv"))
print(auc_tbl, digits = 3, row.names = FALSE)

# ----------------------------------------------------------------
# 5. CALIBRATION (the class probabilities feed the simulations in 30)
# ----------------------------------------------------------------

cat("\nSTEP 5 - Probability calibration (national)\n")

N_BINS <- 10
calib_tbl <- bind_rows(lapply(1:3, function(cls) {
  p   <- oof_df[[prob_cols[cls]]]
  obs <- as.integer(oof_df$y_true == cls)
  bin <- cut(p, breaks = seq(0, 1, length.out = N_BINS + 1), include.lowest = TRUE)
  data.frame(class = CLASS_NAMES[cls], p, obs, bin) |>
    group_by(class, bin) |>
    summarise(n = n(), mean_pred = mean(p), obs_freq = mean(obs), .groups = "drop") |>
    mutate(ece_contrib = n / sum(n) * abs(mean_pred - obs_freq))
}))
ece_tbl <- calib_tbl |> group_by(class) |>
  summarise(ece = sum(ece_contrib), .groups = "drop") |>
  mutate(brier = vapply(match(class, CLASS_NAMES), function(cls)
    mean((oof_df[[prob_cols[cls]]] - (oof_df$y_true == cls))^2), numeric(1)))
write_csv(calib_tbl |> left_join(ece_tbl, by = "class"),
          file.path(metrics_dir, "tab_calibration.csv"))
print(ece_tbl, digits = 3, row.names = FALSE)

# ----------------------------------------------------------------
# 6. FEATURE IMPORTANCE (permutation, CV fold models)
# ----------------------------------------------------------------
# ranger's permutation importance is computed on each model's
# out-of-bag samples, not on the spatial test fold.

cat("\nSTEP 6 - Feature importance\n")

fold_models <- Filter(Negate(is.null), national_model$cv$fold_models)
imp_models  <- if (length(fold_models)) fold_models else national_model$models
feats <- national_model$features
imp_mat <- sapply(imp_models, function(m) ranger::importance(m)[feats])
imp_mat <- matrix(imp_mat, nrow = length(feats), dimnames = list(feats, NULL))
importance_tbl <- data.frame(feature  = feats,
                             imp_mean = rowMeans(imp_mat, na.rm = TRUE),
                             imp_sd   = apply(imp_mat, 1, sd, na.rm = TRUE),
                             n_models = ncol(imp_mat)) |>
  arrange(desc(imp_mean))
write_csv(importance_tbl, file.path(metrics_dir, "tab_feature_importance.csv"))

# ----------------------------------------------------------------
# 7. POWER ANALYSIS
# ----------------------------------------------------------------
# Minimum detectable accuracy difference (two-proportion test, Cohen's
# h) at 80% power and alpha = 0.05, plus n needed for 3/5/10 pp.

cat("\nSTEP 7a - Minimum detectable effect\n")

power_mde <- function(s, power = 0.80, alpha = 0.05) {
  d <- scope_df(s)
  if (nrow(d) < 30) return(NULL)
  p1 <- mean(d$y_pred == d$y_true); n <- nrow(d)
  mde <- tryCatch({
    h <- pwr::pwr.2p.test(n = n, sig.level = alpha, power = power)$h
    sin(asin(sqrt(p1)) + h / 2)^2 - p1
  }, error = function(e) NA_real_)
  n_for <- function(delta) tryCatch({
    h <- pwr::ES.h(min(p1 + delta, 0.999), p1)
    ceiling(pwr::pwr.2p.test(h = h, sig.level = alpha, power = power)$n)
  }, error = function(e) NA_real_)
  data.frame(scope = s, n_test = n, observed_accuracy = p1,
             mde_delta_80pct_power = mde,
             n_needed_for_3pp = n_for(0.03), n_needed_for_5pp = n_for(0.05),
             n_needed_for_10pp = n_for(0.10))
}
power_mde_tbl <- bind_rows(lapply(SCOPES, power_mde))
write_csv(power_mde_tbl, file.path(metrics_dir, "tab_power_mde.csv"))
print(power_mde_tbl, digits = 3, row.names = FALSE)

cat("\nSTEP 7b - McNemar power (from 41)\n")
mcnemar_path <- file.path(PATHS$out_validation, "tables", "mcnemar_tests.csv")
if (file.exists(mcnemar_path)) {
  mcnemar_in <- readr::read_csv(mcnemar_path, show_col_types = FALSE)
  
  # Achieved power at the observed discordance rate, assuming a 65:35
  # split of discordant pairs.
  mcnemar_power_sim <- function(n_samples, disc_rate, p_a_wins = 0.65,
                                alpha = 0.05, n_sim = 2000, seed = 42) {
    set.seed(seed)
    mean(vapply(seq_len(n_sim), function(i) {
      n_disc <- rbinom(1, n_samples, disc_rate)
      if (n_disc < 6) return(FALSE)
      a <- rbinom(1, n_disc, p_a_wins)
      pv <- tryCatch(mcnemar.test(matrix(c(0L, a, n_disc - a, 0L), 2),
                                  correct = FALSE)$p.value,
                     error = function(e) NA_real_)
      isTRUE(pv < alpha)
    }, logical(1)))
  }
  n_for_power <- function(disc_rate, target = 0.80,
                          grid = c(500, 1000, 2000, 4000, 8000, 16000, 32000)) {
    pw  <- vapply(grid, function(n) mcnemar_power_sim(n, disc_rate), numeric(1))
    idx <- which(pw >= target)
    if (length(idx) == 0) NA_integer_ else as.integer(grid[min(idx)])
  }
  mcnemar_power_tbl <- mcnemar_in |>
    mutate(n_discordant     = a_only_correct + b_only_correct,
           discordance_rate = n_discordant / n_paired,
           achieved_power_at_65_35 = vapply(seq_len(n()), function(i)
             mcnemar_power_sim(n_paired[i], discordance_rate[i]), numeric(1)),
           n_for_80pct_power = vapply(seq_len(n()), function(i)
             n_for_power(discordance_rate[i]), integer(1))) |>
    select(comparison, n_paired, n_discordant, discordance_rate, p_value,
           significant_005, achieved_power_at_65_35, n_for_80pct_power)
  write_csv(mcnemar_power_tbl, file.path(metrics_dir, "tab_power_mcnemar.csv"))
  print(mcnemar_power_tbl, digits = 3, row.names = FALSE)
  cat("  Low power means non-significance is not evidence of equivalence.\n")
} else {
  cat(sprintf("  Skipped: run 41 first (expected %s)\n", mcnemar_path))
}

# ================================================================
# 8. FIGURES
# ================================================================

theme_pub <- function(base = 10) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
          strip.background = element_rect(fill = "grey96", colour = NA),
          strip.text       = element_text(face = "bold"),
          plot.title       = element_text(face = "bold", size = base + 1))
}
png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

# --- ROC: national (main) -------------------------------------------
cat("\nSTEP 8a - ROC figures\n")
auc_lbl <- auc_tbl |> mutate(legend_lbl = sprintf("%s (AUC = %.2f)", class, auc))
roc_lbl <- roc_curves |> left_join(auc_lbl |> select(scope, class, legend_lbl),
                                   by = c("scope", "class"))

nat_lbls <- auc_lbl |> filter(scope == "national") |>
  arrange(factor(class, levels = CLASS_NAMES))
nat_roc  <- roc_lbl |> filter(scope == "national") |>
  mutate(legend_lbl = factor(legend_lbl, levels = nat_lbls$legend_lbl))
fig_roc <- ggplot(nat_roc, aes(fpr, tpr, colour = legend_lbl)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey80",
              linetype = "dashed", linewidth = 0.3) +
  geom_path(linewidth = 0.9) +
  scale_colour_manual(values = setNames(CLASS_COLS[nat_lbls$class], nat_lbls$legend_lbl),
                      name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  coord_fixed() +
  labs(x = "False positive rate", y = "True positive rate",
       title = "National model, spatial CV") +
  theme_pub() +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.85),
                                         colour = "grey80", linewidth = 0.2),
        legend.text = element_text(size = 8))
ggsave(file.path(fig_main_dir, "fig_roc_national.png"),
       fig_roc, width = 4.2, height = 4.2, dpi = 600, bg = "white")

# --- ROC: zones x classes (appendix) --------------------------------
zone_roc <- roc_lbl |> filter(scope != "national") |>
  mutate(scope_lbl = scope_fac(scope), class = factor(class, levels = CLASS_NAMES))
zone_ann <- auc_lbl |> filter(scope != "national") |>
  mutate(scope_lbl = scope_fac(scope), class = factor(class, levels = CLASS_NAMES),
         txt = sprintf("AUC = %.2f", auc))
fig_grid <- ggplot(zone_roc, aes(fpr, tpr)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey80",
              linetype = "dashed", linewidth = 0.3) +
  geom_path(aes(colour = class), linewidth = 0.8, show.legend = FALSE) +
  geom_text(data = zone_ann, aes(x = 0.95, y = 0.05, label = txt),
            hjust = 1, vjust = 0, size = 2.8, colour = "grey30") +
  scale_colour_manual(values = CLASS_COLS) +
  facet_grid(scope_lbl ~ class, switch = "y") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  coord_fixed() +
  labs(x = "False positive rate", y = "True positive rate",
       title = "National model, OOF predictions split by latitude zone") +
  theme_pub() + theme(panel.spacing = unit(0.35, "lines"))
ggsave(file.path(fig_appendix_dir, "fig_roc_grid_zones.png"),
       fig_grid, width = 7.5, height = 7.8, dpi = 600, bg = "white")

# --- Confusion matrices ---------------------------------------------
cat("STEP 8b - Confusion matrices\n")

conf_table <- function(s) {
  d  <- scope_df(s)
  yt <- factor(CLASS_NAMES[d$y_true], levels = CLASS_NAMES)
  yp <- factor(CLASS_NAMES[d$y_pred], levels = CLASS_NAMES)
  out <- as.data.frame(table(True = yt, Predicted = yp), stringsAsFactors = FALSE) |>
    rename(Count = Freq) |>
    group_by(True) |>
    mutate(pct = if (sum(Count) > 0) 100 * Count / sum(Count) else NA_real_) |>   # per true class
    ungroup() |>
    mutate(True      = factor(True, levels = CLASS_NAMES),
           Predicted = factor(Predicted, levels = rev(CLASS_NAMES)))
  attr(out, "n")   <- nrow(d)
  attr(out, "acc") <- if (nrow(d)) mean(d$y_true == d$y_pred) else NA_real_
  out
}

make_conf_panel <- function(s, base = 13, show_y = TRUE) {
  df <- conf_table(s)
  if (all(df$Count == 0)) return(patchwork::plot_spacer())
  diag_df <- df[as.character(df$True) == as.character(df$Predicted), ]
  ggplot(df, aes(True, Predicted, fill = pct)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_tile(data = diag_df, fill = NA, colour = "grey10", linewidth = 0.6,
              width = 0.97, height = 0.97) +
    geom_text(aes(label = sprintf("%d", Count), colour = pct > 55),
              size = base * 0.42, fontface = "bold", vjust = -0.15, show.legend = FALSE) +
    geom_text(aes(label = ifelse(is.na(pct), "", sprintf("(%.0f%%)", pct)),
                  colour = pct > 55),
              size = base * 0.32, vjust = 1.25, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"),
                        na.value = "grey15") +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b", limits = c(0, 100),
                        na.value = "grey92", breaks = seq(0, 100, 25),
                        name = "Recall (% of true class)") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "True class", y = if (show_y) "Predicted class" else NULL,
         title = SCOPE_LABELS[[s]],
         subtitle = sprintf("n = %s,  accuracy = %.2f",
                            format(attr(df, "n"), big.mark = ","), attr(df, "acc"))) +
    theme_pub(base = base) +
    theme(aspect.ratio = 1, panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey30", linewidth = 0.4),
          axis.ticks = element_blank(),
          axis.text = element_text(size = base - 1, colour = "grey10"),
          axis.text.y = element_text(angle = 90, hjust = 0.5),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(size = base - 2, colour = "grey30", hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_text(size = base - 2, vjust = 0.85),
          legend.text = element_text(size = base - 3),
          legend.key.width = unit(1.4, "cm"), legend.key.height = unit(0.35, "cm"))
}

conf_nat <- file.path(fig_main_dir, "fig_confusion_national.png")
fig_conf_nat <- make_conf_panel("national", base = 14)
ggsave(conf_nat, fig_conf_nat, device = png_dev, width = 5.6, height = 6.4,
       units = "in", dpi = 600, bg = "white")
ggsave(sub("\\.png$", ".pdf", conf_nat), fig_conf_nat, device = pdf_dev,
       width = 5.6, height = 6.4, units = "in", bg = "white")

fig_conf_z <- Reduce(`|`, lapply(seq_along(ZONE_ORDER), function(i)
  make_conf_panel(ZONE_ORDER[i], show_y = (i == 1)))) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
conf_z <- file.path(fig_appendix_dir, "fig_confusion_zones.png")
ggsave(conf_z, fig_conf_z, device = png_dev, width = 12, height = 5.4,
       units = "in", dpi = 600, bg = "white")
ggsave(sub("\\.png$", ".pdf", conf_z), fig_conf_z, device = pdf_dev,
       width = 12, height = 5.4, units = "in", bg = "white")

for (s in SCOPES) {       # sanity print: diagonal recall
  d     <- conf_table(s)
  diagd <- d[as.character(d$True) == as.character(d$Predicted), ]
  r     <- diagd$pct[match(CLASS_NAMES, as.character(diagd$True))]
  cat(sprintf("  [%-17s] recall W/M/V = %.0f%% / %.0f%% / %.0f%%  (acc %.2f, n %d)\n",
              s, r[1], r[2], r[3], attr(d, "acc"), attr(d, "n")))
}

# --- F1 bars (national + zones) -------------------------------------
cat("STEP 8c - F1 bars\n")
f1_plot_df <- per_class_tbl |>
  mutate(scope_lbl = scope_fac(scope), class = factor(class, levels = CLASS_NAMES))
fig_f1 <- ggplot(f1_plot_df, aes(class, f1_cv_mean, fill = class)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(f1_cv_mean - f1_cv_sd, 0),
                    ymax = pmin(f1_cv_mean + f1_cv_sd, 1)),
                width = 0.2, linewidth = 0.4, na.rm = TRUE) +
  # values printed inside the bar foot so they never collide with the error bars
  geom_text(aes(y = 0.04, label = sprintf("%.2f", f1_cv_mean)), vjust = 0,
            size = 3.4, colour = "white", fontface = "bold") +
  scale_fill_manual(values = CLASS_COLS, guide = "none") +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  facet_wrap(~ scope_lbl, nrow = 1) +
  labs(x = NULL, y = "F1 (CV mean \u00b1 SD across folds)") +
  theme_pub()
ggsave(file.path(fig_main_dir, "fig_f1_bars.png"),
       fig_f1, width = 11, height = 3.4, dpi = 600, bg = "white")

# --- Feature importance ---------------------------------------------
cat("STEP 8d - Feature importance\n")
fig_imp <- ggplot(importance_tbl,
                  aes(imp_mean, reorder(feature, imp_mean))) +
  geom_col(fill = "#4477AA", alpha = 0.85, width = 0.75) +
  geom_errorbar(aes(xmin = pmax(imp_mean - imp_sd, 0), xmax = imp_mean + imp_sd),
                width = 0.25, linewidth = 0.3, colour = "grey30",
                orientation = "y", na.rm = TRUE) +
  labs(x = sprintf("Permutation importance (mean \u00b1 SD, %d CV fold models)",
                   ncol(imp_mat)), y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_blank())
ggsave(file.path(fig_main_dir, "fig_feature_importance.png"),
       fig_imp, width = 6, height = max(3, 0.3 * nrow(importance_tbl) + 1),
       dpi = 600, bg = "white")

# --- Calibration ----------------------------------------------------
cat("STEP 8e - Calibration\n")
ece_lbl  <- setNames(sprintf("%s (ECE = %.3f)", ece_tbl$class, ece_tbl$ece),
                     ece_tbl$class)[CLASS_NAMES]
cal_plot <- calib_tbl |>
  mutate(class_lbl = factor(ece_lbl[class], levels = unname(ece_lbl)))
fig_cal <- ggplot(cal_plot, aes(mean_pred, obs_freq, colour = class)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linetype = "dashed") +
  geom_line(linewidth = 0.6) +
  geom_point(aes(size = n)) +
  scale_colour_manual(values = CLASS_COLS, guide = "none") +
  scale_size_area(max_size = 3.5, name = "points per bin") +
  facet_wrap(~ class_lbl, nrow = 1) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean predicted probability", y = "Observed frequency",
       title = "Reliability of OOF class probabilities (national)") +
  theme_pub() + theme(legend.position = "bottom")
ggsave(file.path(fig_main_dir, "fig_calibration.png"),
       fig_cal, width = 10, height = 4, dpi = 600, bg = "white")

# ================================================================
# 9. LATEX-READY ROWS
# ================================================================

cat("\n\n== LaTeX rows for the performance table ==\n\n")
latex_table <- cv_summary_tbl |>
  left_join(per_class_tbl |> select(scope, class, f1_cv_mean) |>
              pivot_wider(names_from = class, values_from = f1_cv_mean),
            by = "scope")
for (i in seq_len(nrow(latex_table))) {
  r <- latex_table[i, ]
  cat(sprintf("%-22s & %d & %.3f $\\pm$ %.3f & [%.2f, %.2f] & %.2f & %.2f & %.2f \\\\\n",
              SCOPE_LABELS[[r$scope]], as.integer(r$n_test_total),
              r$accuracy, r$acc_sd, r$ci95_lo, r$ci95_hi,
              r$Water, r$Mining, r$Vegetation))
  if (r$scope == "national") cat("\\midrule\n")
}

cat("\nDone.\n")
cat("  CSVs    ->", metrics_dir, "\n")
cat("  Figures ->", fig_main_dir, "\n          ", fig_appendix_dir, "\n")