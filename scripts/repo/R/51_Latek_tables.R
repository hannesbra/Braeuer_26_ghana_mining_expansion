# ================================================================
# 51_Latek_tables.R
# ----------------------------------------------------------------
# Writes booktabs .tex tables from the CSV outputs of 03, 40 and 41.
# Inputs are named explicitly (no filename guessing).
#
# Tables (P5_RESULTS/tables/):
#   table_rf_performance.tex       national model: national + per zone (40)
#   table_approach_comparison.tex  P0 vs M1/M2/M3/B3 (41)
#   table_band_sets.tex            national band sets under spatial CV (03)
#   table_feature_importance.tex   permutation importance, top 15 (40)
#   table_power_analysis.tex       minimum detectable effect (40)
#   table_power_mcnemar.tex        power of the paired tests (40, if run)
# ================================================================

source("R/00_config.R")
pacman::p_load(dplyr, readr, tidyr)

out_tex <- results_tab_dir()   # P5_RESULTS/tables
p1_metrics <- file.path(PATHS$out_unconstrained, "metrics")

SRC <- list(
  cv_summary = file.path(p1_metrics, "tab_cv_summary.csv"),                  # 40
  per_class  = file.path(p1_metrics, "tab_per_class_metrics.csv"),           # 40
  importance = file.path(p1_metrics, "tab_feature_importance.csv"),          # 40
  power_mde  = file.path(p1_metrics, "tab_power_mde.csv"),                   # 40
  power_mcn  = file.path(p1_metrics, "tab_power_mcnemar.csv"),               # 40
  approaches = file.path(PATHS$out_validation, "tables", "master_results_long.csv"),  # 41
  band_sets  = file.path(PATHS$out_validation, "tables", "feature_set_cv_summary.csv") # 03
)

SCOPE_LABELS <- c(national = "National", ZONE_LABELS)
APPROACH_LABELS <- c(
  P0_national_all   = "P0 National, all 27 bands",
  M1_national_top10 = "M1 National, top-10 bands",
  M2_zone_specific  = "M2 Zone models, literature bands",
  M3_zone_all       = "M3 Zone models, all 27 bands",
  B3_random_cv      = "B3 P0 under random CV"
)
BANDSET_LABELS <- c(
  all_candidates      = "All candidate bands",
  literature_union    = "Literature union",
  top_k               = "Top-10 by importance",
  decorrelated        = "Decorrelated top-10",
  nested_top_k        = "Top-10, nested selection",
  nested_decorrelated = "Decorrelated top-10, nested selection"
)

cat("================================================================\n")
cat(" LaTeX TABLE EXPORT\n")
cat("================================================================\n")

# ----------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------

tex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([&%$#_{}])", "\\\\\\1", x)
  x <- gsub("~", "\\\\textasciitilde{}", x)
  x
}

# Column-wise formatting (numeric columns keep their type until here).
fmt_col <- function(x, digits) {
  if (is.numeric(x)) {
    is_int <- all(is.na(x) | abs(x - round(x)) < 1e-9)
    out <- if (is_int) formatC(round(x), format = "d", big.mark = ",")
    else formatC(x, format = "f", digits = digits)
    out[is.na(x)] <- "--"
    out
  } else {
    out <- as.character(x)
    out[is.na(out)] <- "--"
    out
  }
}

# `header` = LaTeX column headings (already escaped/typeset).
# Character cells are escaped unless listed in `raw_cols`.
write_table <- function(df, header, caption, label, file, align = NULL,
                        digits = 3, raw_cols = character(0)) {
  if (is.null(df) || !nrow(df)) { cat(sprintf("  [skip] %s: no rows\n", file)); return(invisible()) }
  stopifnot(length(header) == ncol(df))
  cells <- lapply(names(df), function(nm) {
    v <- fmt_col(df[[nm]], digits)
    if (!is.numeric(df[[nm]]) && !nm %in% raw_cols) v <- tex_escape(v)
    v
  })
  body <- do.call(paste, c(cells, sep = " & "))
  if (is.null(align)) align <- paste0("l", strrep("r", ncol(df) - 1))
  tex <- c("\\begin{table}[!ht]", "\\centering", "\\small",
           sprintf("\\caption{%s}", caption), sprintf("\\label{tab:%s}", label),
           sprintf("\\begin{tabular}{%s}", align), "\\toprule",
           paste0(paste(header, collapse = " & "), " \\\\"), "\\midrule",
           paste0(body, " \\\\"), "\\bottomrule", "\\end{tabular}", "\\end{table}")
  out_p <- file.path(out_tex, file)
  writeLines(tex, out_p)
  cat(sprintf("  -> %s\n", out_p))
}

read_src <- function(key) {
  p <- SRC[[key]]
  if (!file.exists(p)) { cat(sprintf("  [missing] %s (%s)\n", key, p)); return(NULL) }
  read_csv(p, show_col_types = FALSE)
}

# ----------------------------------------------------------------
# 1. NATIONAL MODEL PERFORMANCE (40)
# ----------------------------------------------------------------
cat("\n[1] RF performance\n")
cv <- read_src("cv_summary"); pc <- read_src("per_class")
if (!is.null(cv) && !is.null(pc)) {
  f1 <- pc |> select(scope, class, f1_cv_mean) |>
    pivot_wider(names_from = class, values_from = f1_cv_mean)
  perf <- cv |> left_join(f1, by = "scope") |>
    mutate(scope = unname(SCOPE_LABELS[scope]),
           acc   = sprintf("%.3f $\\pm$ %.3f", accuracy, acc_sd),
           ci    = sprintf("[%.3f, %.3f]", ci95_lo, ci95_hi)) |>
    transmute(scope, n_test_total, acc, ci, kappa, Water, Mining, Vegetation)
  write_table(perf,
              header = c("", "$n$", "Accuracy", "95\\% CI", "$\\kappa$",
                         "$F_1$ Water", "$F_1$ Mining", "$F_1$ Vegetation"),
              caption = paste(
                "Performance of the national random forest (27 bands) under spatial",
                "5-fold cross-validation with a variogram-derived buffer. Zone rows",
                "split the same out-of-fold predictions by latitude zone; the model",
                "does not use zones. Accuracy: mean $\\pm$ SD across folds."),
              label = "rf_performance", file = "table_rf_performance.tex",
              raw_cols = c("acc", "ci"))
}

# ----------------------------------------------------------------
# 2. APPROACH COMPARISON (41)
# ----------------------------------------------------------------
cat("\n[2] Approach comparison\n")
ap <- read_src("approaches")
if (!is.null(ap)) {
  ap_w <- ap |>
    filter(metric %in% c("accuracy", "kappa", "mining_f1", "water_f1",
                         "vegetation_f1", "log_loss")) |>
    select(approach, metric, mean) |>
    pivot_wider(names_from = metric, values_from = mean) |>
    left_join(ap |> filter(metric == "accuracy") |>
                select(approach, ci_lo, ci_hi), by = "approach") |>
    mutate(ci = sprintf("[%.3f, %.3f]", ci_lo, ci_hi),
           approach = factor(approach, levels = names(APPROACH_LABELS))) |>
    arrange(approach) |>
    mutate(approach = unname(APPROACH_LABELS[as.character(approach)])) |>
    select(approach, accuracy, ci, kappa, water_f1, mining_f1, vegetation_f1, log_loss)
  write_table(ap_w,
              header = c("Approach", "Accuracy", "95\\% CI", "$\\kappa$",
                         "$F_1$ Water", "$F_1$ Mining", "$F_1$ Vegetation", "Log-loss"),
              caption = paste(
                "Modelling approaches under identical spatial cross-validation",
                "(shared folds and buffer). P0 is the model used for all maps.",
                "B3 repeats P0 with random folds to show the optimism of",
                "non-spatial validation."),
              label = "approach_comparison", file = "table_approach_comparison.tex",
              raw_cols = "ci")
}

# ----------------------------------------------------------------
# 3. BAND SETS (03)
# ----------------------------------------------------------------
cat("\n[3] Band sets\n")
bs <- read_src("band_sets")
if (!is.null(bs)) {
  bs_t <- bs |>
    mutate(config = factor(config, levels = names(BANDSET_LABELS))) |>
    arrange(config) |>
    transmute(config = unname(BANDSET_LABELS[as.character(config)]),
              n_bands, selection_on, mean_accuracy,
              ci = sprintf("[%.3f, %.3f]", ci95_lower, ci95_upper),
              mean_kappa, mean_mining_f1, mean_log_loss)
  write_table(bs_t,
              header = c("Band set", "Bands", "Selected on", "Accuracy", "95\\% CI",
                         "$\\kappa$", "$F_1$ Mining", "Log-loss"),
              caption = paste(
                "National band sets under spatial cross-validation (shared folds).",
                "Nested rows repeat the band selection inside each training fold and",
                "give an unbiased estimate of the selection procedure. The final",
                "model uses all candidate bands."),
              label = "band_sets", file = "table_band_sets.tex",
              align = "lrlrrrrr", raw_cols = "ci")
}

# ----------------------------------------------------------------
# 4. FEATURE IMPORTANCE (40)
# ----------------------------------------------------------------
cat("\n[4] Feature importance\n")
fi <- read_src("importance")
if (!is.null(fi)) {
  fi_t <- fi |> arrange(desc(imp_mean)) |> head(15) |>
    transmute(rank = row_number(), feature, imp_mean, imp_sd)
  write_table(fi_t,
              header = c("Rank", "Band", "Importance (mean)", "SD"),
              caption = paste(
                "Top 15 bands of the national model by permutation importance",
                "(mean and SD across the spatial-CV fold models; computed on each",
                "model's out-of-bag samples)."),
              label = "feature_importance", file = "table_feature_importance.tex",
              align = "rlrr", digits = 4)
}

# ----------------------------------------------------------------
# 5. POWER (40)
# ----------------------------------------------------------------
cat("\n[5] Power analysis\n")
pw <- read_src("power_mde")
if (!is.null(pw)) {
  pw_t <- pw |> transmute(scope = unname(SCOPE_LABELS[scope]), n_test,
                          observed_accuracy, mde_delta_80pct_power,
                          n_needed_for_3pp, n_needed_for_5pp, n_needed_for_10pp)
  write_table(pw_t,
              header = c("", "$n$", "Accuracy", "MDE (80\\% power)",
                         "$n$ for 3 pp", "$n$ for 5 pp", "$n$ for 10 pp"),
              caption = paste(
                "Minimum detectable accuracy difference at 80\\% power",
                "($\\alpha = 0.05$, two-proportion test) and the test-set size",
                "needed to detect differences of 3, 5 and 10 percentage points."),
              label = "power_analysis", file = "table_power_analysis.tex")
}
pm <- read_src("power_mcn")
if (!is.null(pm)) {
  pm_t <- pm |> transmute(comparison, n_paired, n_discordant, p_value,
                          achieved_power_at_65_35, n_for_80pct_power)
  write_table(pm_t,
              header = c("Comparison", "$n$", "Discordant", "$p$",
                         "Power (65:35)", "$n$ for 80\\% power"),
              caption = paste(
                "Power of the paired McNemar tests at the observed discordance rate,",
                "assuming a 65:35 split of discordant pairs."),
              label = "power_mcnemar", file = "table_power_mcnemar.tex")
}

cat("\n================================================================\n")
cat("Include with \\input{tables/<file>}; needs \\usepackage{booktabs}\n")
cat("================================================================\n")