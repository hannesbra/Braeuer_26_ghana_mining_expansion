# ================================================================
# 41_approach_comparison.R
# ----------------------------------------------------------------
# Compares the main national approach (P0) with benchmarks under
# identical spatial CV (shared national folds and buffer), so the OOF
# predictions are paired and McNemar's test applies.
#
#   P0  National model, all 27 bands (config)           <- main text
#   M1  National model, decorrelated top-10 from 03      <- band selection
#   M2  Zone models, zone-specific literature bands      <- appendix
#       (the former main approach)
#   M3  Zone models, P0's 27 bands                       <- appendix
#       (separates the effect of zoning from the effect of the bands)
#   B3  P0 under RANDOM CV (no blocks, no buffer)        <- optimism audit
#
# Zone models are trained only on points of their own latitude zone
# and predict only that zone.
#
# P0 uses the same seeds and settings as 10, so its OOF predictions
# reproduce metrics/oof_predictions.csv (checked below).
#
# Outputs (PATHS$out_validation):
#   tables/master_results_long.csv, master_results.csv
#   tables/results_by_zone.csv          approach x latitude zone
#   tables/mcnemar_tests.csv            paired tests (read by 40)
#   tables/optimism_audit.csv           B3 random - P0 spatial
#   tables/variable_importance_zone_models.csv
#   tables/benchmark_oof.csv
#   P5_RESULTS/figures/41_approach_comparison/approach_comparison.png
#   models/benchmarks.rds
# ================================================================

source("R/00_config.R")
source("R/01_utils.R")

pacman::p_load(terra, sf, dplyr, tidyr, readr, ranger, blockCV, ggplot2, parallel)

terraOptions(memfrac = 0.6, tempdir = tempdir())
set.seed(CV_CONFIG$seed)
start_time <- Sys.time()

out_dir <- PATHS$out_validation
ensure_output_dirs(out_dir, subdirs = c("tables", "models"))
plot_dir <- results_fig_dir("41_approach_comparison")

APPROACH_LABELS <- c(
  P0_national_all      = "P0 National, all 27 bands",
  M1_national_top10    = "M1 National, top-10 bands (03)",
  M2_zone_specific     = "M2 Zone models, literature bands",
  M3_zone_all          = "M3 Zone models, all 27 bands",
  B3_random_cv         = "B3 P0 under random CV"
)

cat("================================================================\n")
cat("41 - APPROACH COMPARISON: P0 vs M1 / M2 / M3 / B3\n")
cat("================================================================\n")

# ================================================================
# 1. DATA, BAND SETS AND FOLDS
# ================================================================

cat("\n[1/6] Training data, band sets and folds\n")

df <- prepare_training_data()
y  <- as.integer(df$landcover_3class)
n  <- nrow(df)

p0_features <- read_selected_features()$features        # locked 27 bands
if (!file.exists(NATIONAL_CONFIG$selected_file))
  stop("Top-10 set from 03 not found: ", NATIONAL_CONFIG$selected_file, call. = FALSE)
m1_sel      <- read_csv(NATIONAL_CONFIG$selected_file, show_col_types = FALSE)
m1_features <- m1_sel$feature[order(m1_sel$order)]
# Zone models drop bands that are constant within a zone
zone_specific <- lapply(REGIONAL_FEATURES, function(f) intersect(f, p0_features))
zone_all      <- setNames(rep(list(p0_features), length(ZONE_ORDER)), ZONE_ORDER)

cat(sprintf("  P0 bands (%d): %s\n", length(p0_features),
            paste(p0_features, collapse = ", ")))
cat(sprintf("  M1 bands (%d, %s rule): %s\n", length(m1_features), m1_sel$rule[1],
            paste(m1_features, collapse = ", ")))
cat(sprintf("  M2 bands per zone: %s\n", paste(lengths(zone_specific), collapse = " / ")))

folds    <- get_national_folds(df)
fold_ids <- folds$fold_ids

# Random folds for B3 (same fold sizes, no spatial structure, no buffer)
set.seed(CV_CONFIG$seed)
random_fold_ids <- sample(rep(seq_len(CV_CONFIG$n_folds), length.out = n))
random_masks    <- lapply(seq_len(CV_CONFIG$n_folds),
                          function(k) which(random_fold_ids != k))

cat("  rows per zone: ",
    paste(sprintf("%s = %d", ZONE_ORDER, as.integer(table(factor(df$zone, ZONE_ORDER)))),
          collapse = ", "), "\n", sep = "")

# ================================================================
# 2. CV ENGINE
# ================================================================
# approach "national": one model per fold on all training rows.
# approach "zone":     one model per fold and zone, trained and
#                      evaluated on that zone's rows only.

run_cv <- function(approach, features, fold_ids, masks, label) {
  prob <- matrix(NA_real_, n, 3, dimnames = list(NULL, CLASS_LEVELS))
  for (k in seq_len(CV_CONFIG$n_folds)) {
    te_all <- which(fold_ids == k)
    tr_all <- masks[[k]]
    if (approach == "national") {
      if (length(te_all) < 10 || length(tr_all) < 30 ||
          length(unique(y[tr_all])) < 2) next
      m <- fit_rf_weighted(df[tr_all, features, drop = FALSE], y[tr_all],
                           seed = CV_CONFIG$seed + k, importance = "none")
      prob[te_all, ] <- predict_probs(m, df[te_all, features, drop = FALSE])
    } else if (approach == "zone") {
      for (zi in seq_along(ZONE_ORDER)) {
        z  <- ZONE_ORDER[zi]
        te <- te_all[df$zone[te_all] == z]
        tr <- tr_all[df$zone[tr_all] == z]
        f  <- features[[z]]
        if (!length(te) || length(tr) < 30 || length(unique(y[tr])) < 2 ||
            length(f) < 2) next
        m <- fit_rf_weighted(df[tr, f, drop = FALSE], y[tr],
                             seed = CV_CONFIG$seed + 1000 * zi + k,
                             importance = "none")
        prob[te, ] <- predict_probs(m, df[te, f, drop = FALSE])
      }
    } else stop("Unknown approach: ", approach)
  }
  pred <- rep(NA_integer_, n)
  ok   <- !is.na(prob[, 1])
  pred[ok] <- argmax_class(prob[ok, , drop = FALSE])
  fm <- fold_metrics(y, pred, fold_ids, prob)
  s  <- summarise_fold_metrics(fm)
  cat(sprintf("  %-22s acc %.3f [%.3f-%.3f]  kappa %.3f  mining F1 %.3f  (predicted %d/%d)\n",
              label, s$mean_accuracy, s$ci95_lower, s$ci95_upper, s$mean_kappa,
              s$mean_mining_f1, sum(ok), n))
  list(approach = label, prob = prob, pred = pred, fold_ids = fold_ids,
       fold_details = fm, summary = s)
}

# ================================================================
# 3. RUN
# ================================================================

cat("\n[2/6] Spatial CV per approach\n")

B <- list(
  P0_national_all   = run_cv("national", p0_features,   fold_ids, folds$buffer_masks,
                             "P0_national_all"),
  M1_national_top10 = run_cv("national", m1_features,   fold_ids, folds$buffer_masks,
                             "M1_national_top10"),
  M2_zone_specific  = run_cv("zone",     zone_specific, fold_ids, folds$buffer_masks,
                             "M2_zone_specific"),
  M3_zone_all       = run_cv("zone",     zone_all,      fold_ids, folds$buffer_masks,
                             "M3_zone_all"),
  B3_random_cv      = run_cv("national", p0_features,   random_fold_ids, random_masks,
                             "B3_random_cv")
)

# ---- Consistency with 10 ------------------------------------------
oof10_path <- file.path(PATHS$out_unconstrained, "metrics", "oof_predictions.csv")
if (file.exists(oof10_path)) {
  oof10 <- read_csv(oof10_path, show_col_types = FALSE)
  j <- match(oof10$row_id, df$row_id)
  agree <- mean(oof10$y_pred == B$P0_national_all$pred[j], na.rm = TRUE)
  cat(sprintf("\n  P0 vs 10's OOF predictions: %.1f%% identical%s\n", 100 * agree,
              if (agree < 0.999) "  (<100%: 10 was run with another band set or folds; rerun 10)" else ""))
}

# ================================================================
# 4. RESULTS TABLES
# ================================================================

cat("\n[3/6] Results tables\n")

summary_long <- function(bench) {
  s <- bench$summary
  if (is.null(s)) return(NULL)
  data.frame(approach = bench$approach,
             metric   = c("accuracy", "kappa", "water_f1", "mining_f1",
                          "vegetation_f1", "log_loss"),
             mean     = c(s$mean_accuracy, s$mean_kappa, s$mean_water_f1,
                          s$mean_mining_f1, s$mean_vegetation_f1, s$mean_log_loss),
             sd       = c(s$sd_accuracy, s$sd_kappa, NA, s$sd_mining_f1, NA, NA),
             ci_lo    = c(s$ci95_lower, NA, NA, NA, NA, NA),
             ci_hi    = c(s$ci95_upper, NA, NA, NA, NA, NA),
             n_folds  = s$n_folds,
             n_test   = s$n_test,
             stringsAsFactors = FALSE)
}

master_long <- bind_rows(lapply(B, summary_long))
write_csv(master_long, file.path(out_dir, "tables", "master_results_long.csv"))

master_wide <- master_long |>
  mutate(display = ifelse(!is.na(ci_lo),
                          sprintf("%.3f [%.3f-%.3f]", mean, ci_lo, ci_hi),
                          sprintf("%.3f", mean))) |>
  select(approach, metric, display) |>
  pivot_wider(names_from = approach, values_from = display)
write_csv(master_wide, file.path(out_dir, "tables", "master_results.csv"))
print(master_wide, n = 20)

# ---- By latitude zone (same OOF predictions, split for reporting) ---
by_zone <- bind_rows(lapply(B, function(b) bind_rows(lapply(ZONE_ORDER, function(z) {
  i <- df$zone == z
  s <- summarise_fold_metrics(fold_metrics(y[i], b$pred[i], b$fold_ids[i],
                                           b$prob[i, , drop = FALSE]))
  if (is.null(s)) return(NULL)
  data.frame(approach = b$approach, zone = z, n_test = s$n_test,
             accuracy = s$mean_accuracy, ci95_lower = s$ci95_lower,
             ci95_upper = s$ci95_upper, kappa = s$mean_kappa,
             mining_f1 = s$mean_mining_f1, water_f1 = s$mean_water_f1,
             vegetation_f1 = s$mean_vegetation_f1, log_loss = s$mean_log_loss)
}))))
write_csv(by_zone, file.path(out_dir, "tables", "results_by_zone.csv"))
cat("\n  Accuracy / mining F1 by latitude zone:\n")
print(by_zone |> select(approach, zone, n_test, accuracy, mining_f1),
      digits = 3, row.names = FALSE)

# ---- Optimism audit (random - spatial CV, same bands) --------------
optimism <- master_long |>
  filter(approach %in% c("P0_national_all", "B3_random_cv")) |>
  select(approach, metric, mean) |>
  pivot_wider(names_from = approach, values_from = mean) |>
  mutate(optimism_bias = B3_random_cv - P0_national_all)
write_csv(optimism, file.path(out_dir, "tables", "optimism_audit.csv"))
cat("\n  Optimism audit (B3 random - P0 spatial):\n")
print(optimism, digits = 3, row.names = FALSE)

# ================================================================
# 5. McNEMAR TESTS (paired, shared spatial folds)
# ================================================================

cat("\n[4/6] McNemar tests\n")

local({     # sanity check: A clearly better than B
  truth <- rep(1L, 100)
  pa <- truth; pa[1:5]  <- 2L
  pb <- truth; pb[1:40] <- 2L
  chk <- run_mcnemar(pa, pb, truth, "A", "B")
  stopifnot(chk$better == "A", chk$a_only_correct == 35, chk$b_only_correct == 0)
})

pairs <- list(c("P0_national_all", "M1_national_top10"),
              c("P0_national_all", "M2_zone_specific"),
              c("P0_national_all", "M3_zone_all"),
              c("M3_zone_all",     "M2_zone_specific"))
short <- function(a) sub("_.*$", "", a)
mcnemar_df <- bind_rows(lapply(pairs, function(p)
  run_mcnemar(B[[p[1]]]$pred, B[[p[2]]]$pred, y, short(p[1]), short(p[2]))))
write_csv(mcnemar_df, file.path(out_dir, "tables", "mcnemar_tests.csv"))
print(mcnemar_df |> select(comparison, n_paired, a_only_correct, b_only_correct,
                           p_value, test_used, better),
      row.names = FALSE)

# ================================================================
# 6. IMPORTANCE OF THE ZONE MODELS (appendix)
# ================================================================

cat("\n[5/6] Permutation importance of the zone models (M2, all zone data)\n")

importance_rows <- lapply(seq_along(ZONE_ORDER), function(zi) {
  z  <- ZONE_ORDER[zi]
  i  <- which(df$zone == z)
  f  <- zone_specific[[z]]
  if (length(i) < 30 || length(unique(y[i])) < 2) return(NULL)
  m  <- fit_rf_weighted(df[i, f, drop = FALSE], y[i],
                        seed = CV_CONFIG$seed + 1000 * zi,
                        threads = max(1L, detectCores() - 1L))
  imp <- data.frame(zone = z, feature = names(m$variable.importance),
                    importance = as.numeric(m$variable.importance)) |>
    arrange(desc(importance))
  cat(sprintf("  %-18s top 5: %s\n", z, paste(head(imp$feature, 5), collapse = ", ")))
  imp
})
write_csv(bind_rows(importance_rows),
          file.path(out_dir, "tables", "variable_importance_zone_models.csv"))

# ================================================================
# 7. FIGURE, OOF AND OBJECTS
# ================================================================

cat("\n[6/6] Figure and saved objects\n")

plot_df <- bind_rows(
  master_long |> filter(metric == "accuracy") |>
    transmute(approach, zone = "all", metric = "Accuracy", mean,
              lo = ci_lo, hi = ci_hi),
  master_long |> filter(metric == "mining_f1") |>
    transmute(approach, zone = "all", metric = "Mining F1", mean,
              lo = pmax(mean - sd, 0), hi = pmin(mean + sd, 1))
) |>
  mutate(approach = factor(approach, levels = names(APPROACH_LABELS),
                           labels = unname(APPROACH_LABELS)))

fig <- ggplot(plot_df, aes(mean, approach, fill = grepl("^P0", approach))) +
  geom_col(width = 0.65, colour = "grey20", linewidth = 0.2) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, linewidth = 0.4,
                orientation = "y", na.rm = TRUE) +
  geom_text(aes(label = sprintf("%.3f", mean), x = pmax(hi, mean, na.rm = TRUE)),
            hjust = -0.2, size = 3) +
  facet_wrap(~ metric, nrow = 1) +
  scale_x_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#BBBBBB"), guide = "none") +
  labs(x = "Spatial CV (accuracy: mean ± 95% CI; mining F1: mean ± SD)", y = NULL,
       title = "Modelling approaches under identical spatial cross-validation",
       subtitle = sprintf("Shared folds, buffer %.1f km (B3: random folds, no buffer), n = %d",
                          folds$buffer_m / 1000, n)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey96", colour = NA),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))
ggsave(file.path(plot_dir, "approach_comparison.png"), fig,
       width = 10, height = 3.6, dpi = 600, bg = "white")

oof_long <- bind_rows(lapply(B, function(b)
  data.frame(approach = b$approach, row_id = df$row_id, zone = df$zone,
             fold = b$fold_ids, y_true = y, y_pred = b$pred,
             prob_water = b$prob[, 1], prob_mining = b$prob[, 2],
             prob_veg = b$prob[, 3])))
write_csv(oof_long, file.path(out_dir, "tables", "benchmark_oof.csv"))

saveRDS(list(benchmarks      = lapply(B, function(b) b[c("approach", "pred", "prob",
                                                         "fold_details", "summary")]),
             features        = list(P0 = p0_features, M1 = m1_features,
                                    M2 = zone_specific, M3 = zone_all),
             folds           = folds[c("fold_ids", "buffer_m", "block_size")],
             random_fold_ids = random_fold_ids,
             row_id          = df$row_id),
        file.path(out_dir, "models", "benchmarks.rds"))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("\n================================================================\n")
cat(sprintf("Approach comparison complete in %.1f min\n", as.numeric(elapsed)))
cat("================================================================\n")