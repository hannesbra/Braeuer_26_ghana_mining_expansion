# ================================================================
# 10_Training_Prediction.R
# ----------------------------------------------------------------
# MAIN classification: one NATIONAL random forest on all 27 usable
# candidate bands (locked in NATIONAL_CONFIG$features). Spatial block
# CV on the shared national folds (variogram buffer), a bootstrapped
# ensemble on all training points, then yearly rasters 2000-2019:
# class, normalised entropy and 3-band class probabilities.
#
# 03_feature_importance_validation. wrote the shared fold cache used here. With the same folds, band
# order and seeds, the CV here should reproduce 03's "all_candidates"
# result (accuracy 0.870). The zone-specific approach is evaluated
# only as a benchmark in 41.
#
# Outputs (PATHS$out_unconstrained):
#   models/national_model.rds         ensemble, CV fold models, buffer diagnostics
#   models/national_features.csv      band set used (audit trail)
#   metrics/oof_predictions.csv       OOF predictions (zone column for reporting)
#   metrics/cv_fold_details.csv       per-fold metrics
#   metrics/cv_summary_national.csv   fold summary
#   metrics/cv_summary_by_zone.csv    same OOF predictions split by latitude zone
#   statistics/yearly_class_stats_national.csv
# Yearly rasters (PATHS$national_predictions):
#   prediction_{year}.tif, uncertainty_{year}.tif
#   probabilities/prob_{year}.tif
# ================================================================

source("R/00_config.R")
source("R/01_utils.R")

pacman::p_load(terra, sf, dplyr, readr, ranger, progress, parallel, jsonlite)

terraOptions(memfrac = 0.6, tempdir = tempdir())
set.seed(CV_CONFIG$seed)
start_time <- Sys.time()

cfg      <- NATIONAL_CONFIG
out_dir  <- PATHS$out_unconstrained
ensure_output_dirs(out_dir)
pred_dir <- PATHS$national_predictions
prob_dir <- PATHS$national_probabilities
dir.create(prob_dir, showWarnings = FALSE, recursive = TRUE)   # creates both

# Layer index of the probability stack must equal the class id
stopifnot(identical(unname(CLASS_IDS), 1:3))
N_THREADS <- max(1L, detectCores() - 1L)

cat("================================================================\n")
cat("P1 - NATIONAL MODEL: TRAIN & PREDICT\n")
cat("================================================================\n")
cat(sprintf("Model/metrics: %s\nPredictions:   %s\n\n", out_dir, pred_dir))

# ================================================================
# 1. FEATURE SET
# ================================================================

cat("[1/6] Feature set\n")
fs       <- read_selected_features(cfg)
features <- fs$features
cat(sprintf("  source  : %s\n", fs$source))
cat(sprintf("  features: %d  (%s)\n", length(features),
            paste(features, collapse = ", ")))
write_csv(data.frame(order = seq_along(features), feature = features,
                     source = fs$source),
          file.path(out_dir, "models", "national_features.csv"))

# ================================================================
# 2. TRAINING DATA
# ================================================================

cat("\n[2/6] Training data\n")
if (normalizePath(cfg$train_composite, mustWork = FALSE) !=
    normalizePath(sprintf(cfg$year_composite_fmt, PREDICT_CONFIG$base_year),
                  mustWork = FALSE))
  cat(sprintf("  NOTE: training composite (%s) is not the %d prediction composite; ",
              basename(cfg$train_composite), PREDICT_CONFIG$base_year),
      "both must have identical band processing.\n", sep = "")

# Same rows as 41/03 (NA dropped over all candidate bands), so the
# shared fold cache applies.
training_df <- prepare_training_data(cfg$train_composite)
cat("  class counts:\n");  print(table(CLASS_NAMES[training_df$landcover_3class]))
cat("  latitude zones (reporting only):\n"); print(table(training_df$zone))

X     <- training_df[, features, drop = FALSE]
y_int <- as.integer(training_df$landcover_3class)
n     <- nrow(training_df)

cell_ha <- prod(terra::res(terra::rast(cfg$train_composite))) / 1e4   # UTM grid

# ================================================================
# 3. SPATIAL CV (shared national folds)
# ================================================================

cat("\n[3/6] Spatial cross-validation\n")
folds    <- get_national_folds(training_df)
fold_ids <- folds$fold_ids

oof_prob    <- matrix(NA_real_, n, 3, dimnames = list(NULL, CLASS_LEVELS))
fold_models <- vector("list", CV_CONFIG$n_folds)
n_train     <- setNames(rep(NA_integer_, CV_CONFIG$n_folds), seq_len(CV_CONFIG$n_folds))

for (k in seq_len(CV_CONFIG$n_folds)) {
  te <- which(fold_ids == k)
  tr <- folds$buffer_masks[[k]]
  if (length(te) < 10 || length(tr) < 30 || length(unique(y_int[tr])) < 2) {
    message(sprintf("  fold %d skipped (train = %d, test = %d)", k, length(tr), length(te)))
    next
  }
  # Same seed/settings as cv_scope() in 03 and P0 in 41
  m <- fit_rf_weighted(X[tr, , drop = FALSE], y_int[tr],
                       seed = CV_CONFIG$seed + k,
                       sample_frac = cfg$cv_sample_frac)
  oof_prob[te, ]   <- predict_probs(m, X[te, , drop = FALSE])
  fold_models[[k]] <- m
  n_train[k]       <- length(tr)
  cat(sprintf("  fold %d: train %d (after %.0f m buffer), test %d\n",
              k, length(tr), folds$buffer_m, length(te)))
}

has_pred <- !is.na(oof_prob[, 1])
oof_pred <- rep(NA_integer_, n)
oof_pred[has_pred] <- argmax_class(oof_prob[has_pred, , drop = FALSE])

oof_df <- data.frame(
  row_id      = training_df$row_id,
  zone        = training_df$zone,
  fold        = as.integer(fold_ids),
  x_utm       = training_df$x_utm,
  y_utm       = training_df$y_utm,
  y_true      = y_int,
  y_pred      = oof_pred,
  prob_water  = oof_prob[, 1],
  prob_mining = oof_prob[, 2],
  prob_veg    = oof_prob[, 3],
  stringsAsFactors = FALSE
)[has_pred, , drop = FALSE]
write_csv(oof_df, file.path(out_dir, "metrics", "oof_predictions.csv"))

prob_cols    <- as.matrix(oof_df[, c("prob_water", "prob_mining", "prob_veg")])
fold_details <- fold_metrics(oof_df$y_true, oof_df$y_pred, oof_df$fold, prob_cols)
fold_details$n_train <- n_train[as.character(fold_details$fold)]
cv_summary   <- summarise_fold_metrics(fold_details)
cv_by_zone   <- bind_rows(lapply(ZONE_ORDER, function(z) {
  i <- oof_df$zone == z
  s <- summarise_fold_metrics(fold_metrics(oof_df$y_true[i], oof_df$y_pred[i],
                                           oof_df$fold[i], prob_cols[i, , drop = FALSE]))
  if (is.null(s)) NULL else cbind(zone = z, s)
}))

write_csv(fold_details, file.path(out_dir, "metrics", "cv_fold_details.csv"))
write_csv(cv_summary,   file.path(out_dir, "metrics", "cv_summary_national.csv"))
write_csv(cv_by_zone,   file.path(out_dir, "metrics", "cv_summary_by_zone.csv"))

cat(sprintf("\n  National CV: acc %.3f (95%% CI %.3f-%.3f), kappa %.3f, mining F1 %.3f\n",
            cv_summary$mean_accuracy, cv_summary$ci95_lower, cv_summary$ci95_upper,
            cv_summary$mean_kappa, cv_summary$mean_mining_f1))
cat("  Same OOF predictions by latitude zone:\n")
print(cv_by_zone[, c("zone", "n_test", "n_folds", "mean_accuracy",
                     "mean_kappa", "mean_mining_f1")],
      digits = 3, row.names = FALSE)

# ================================================================
# 4. FINAL ENSEMBLE (all training points)
# ================================================================

cat(sprintf("\n[4/6] Ensemble of %d bootstrapped ranger models\n", CV_CONFIG$n_models))

ensemble <- vector("list", CV_CONFIG$n_models)
for (m in seq_len(CV_CONFIG$n_models)) {
  set.seed(CV_CONFIG$seed + 100 * m)
  boot <- sample(seq_len(n), replace = TRUE)
  ensemble[[m]] <- fit_rf_weighted(X[boot, , drop = FALSE], y_int[boot],
                                   seed = CV_CONFIG$seed + 100 * m,
                                   sample_frac = cfg$ens_sample_frac,
                                   threads = N_THREADS)
}

national_model <- list(
  approach         = "national",
  models           = ensemble,
  features         = features,
  feature_source   = fs$source,
  class_levels     = CLASS_LEVELS,
  class_counts     = table(factor(y_int, levels = CLASS_LEVELS)),
  n_train          = n,
  train_composite  = cfg$train_composite,
  rf_settings      = cfg[c("num_trees", "min_node_size",
                           "cv_sample_frac", "ens_sample_frac")],
  # buffer / CV diagnostics (used by 40 and 50)
  buffer_strategy  = CV_CONFIG$buffer_strategy,
  buffer_used_m    = folds$buffer_m,
  vgm_info         = folds$vgm_info,
  variogram_target = folds$variogram_target,
  block_size_m     = folds$block_size,
  cv = list(fold_ids        = fold_ids,
            row_id          = training_df$row_id,
            fold_models     = fold_models,
            fold_details    = fold_details,
            summary         = cv_summary,
            summary_by_zone = cv_by_zone),
  created = Sys.time()
)
# Warn if the previous model used other bands: "missing" mode would
# keep its rasters for the years already on disk.
prev_features <- if (file.exists(cfg$model_path))
  tryCatch(readRDS(cfg$model_path)$features, error = function(e) NULL) else NULL
model_changed <- !is.null(prev_features) && !identical(prev_features, features)
saveRDS(national_model, cfg$model_path)
cat(sprintf("  saved: %s\n", cfg$model_path))

# ================================================================
# 5. PREDICT ALL YEARS
# ================================================================
# recompute_mode (NATIONAL_CONFIG):
#   "missing" - only years not yet on disk in PATHS$national_predictions
#   "all"     - regenerate every year
#   "none"    - skip prediction (training-only run)

existing_files <- list.files(pred_dir, pattern = "^prediction_\\d{4}\\.tif$")
existing_years <- as.integer(gsub("prediction_(\\d{4})\\.tif", "\\1", existing_files))
years_to_do    <- switch(cfg$recompute_mode,
                         missing = setdiff(PREDICT_CONFIG$predict_years, existing_years),
                         all     = PREDICT_CONFIG$predict_years,
                         none    = integer(0),
                         stop("Unknown recompute_mode: ", cfg$recompute_mode))

if (model_changed && cfg$recompute_mode == "missing" && length(existing_years))
  warning("The band set differs from the previous model, but ", length(existing_years),
          " years already on disk are kept (recompute_mode = 'missing'). ",
          "Set recompute_mode = 'all' to regenerate them.", call. = FALSE)

cat(sprintf("\n[5/6] Prediction: %d years already on disk, %d to process\n",
            length(existing_years), length(years_to_do)))

stats_rows <- list()

if (length(years_to_do) > 0) {
  pb <- progress_bar$new(format = "  Year :year [:bar] :percent  eta :eta",
                         total = length(years_to_do), width = 70, clear = FALSE)
  for (yr in years_to_do) {
    pb$tick(tokens = list(year = yr))
    cp <- sprintf(cfg$year_composite_fmt, yr)
    if (!file.exists(cp)) { warning("Composite missing for ", yr); next }
    
    tryCatch({
      comp <- load_composite(cp, bands = features)
      res  <- predict_ensemble_raster(comp, national_model$models, features)
      
      writeRaster(res$pred, file.path(pred_dir, sprintf("prediction_%d.tif", yr)),
                  overwrite = TRUE, datatype = "INT1U")
      writeRaster(res$uncert, file.path(pred_dir, sprintf("uncertainty_%d.tif", yr)),
                  overwrite = TRUE, datatype = "FLT4S")
      if (isTRUE(PREDICT_CONFIG$write_probabilities))
        writeRaster(res$prob, file.path(prob_dir, sprintf("prob_%d.tif", yr)),
                    overwrite = TRUE, datatype = "FLT4S")
      
      fr <- freq(res$pred)
      if (!is.null(fr) && nrow(fr) > 0)
        stats_rows[[as.character(yr)]] <- data.frame(
          year = yr, class = as.integer(fr$value),
          class_name = CLASS_NAMES[as.integer(fr$value)],
          count = fr$count, area_ha = fr$count * cell_ha)
      
      rm(comp, res); gc(verbose = FALSE)
    }, error = function(e) {
      message(sprintf("ERROR in year %d: %s", yr, conditionMessage(e)))
    })
  }
}

# ================================================================
# 6. YEARLY CLASS STATISTICS
# ================================================================
# New years are appended; rows of re-predicted years are replaced.

cat("\n[6/6] Yearly class statistics\n")
stats_path <- file.path(out_dir, "statistics", "yearly_class_stats_national.csv")
if (length(stats_rows) > 0) {
  stats_df <- bind_rows(stats_rows)
  if (file.exists(stats_path)) {
    prev <- read_csv(stats_path, show_col_types = FALSE)
    stats_df <- bind_rows(prev |> filter(!year %in% stats_df$year), stats_df)
  }
  stats_df <- stats_df |> arrange(year, class)
  write_csv(stats_df, stats_path)
  
  cat(sprintf("  -> %s\n", stats_path))
} else {
  cat("  nothing new to write\n")
}

missing_years <- setdiff(PREDICT_CONFIG$predict_years,
                         as.integer(gsub("prediction_(\\d{4})\\.tif", "\\1",
                                         list.files(pred_dir, pattern = "^prediction_\\d{4}\\.tif$"))))
if (length(missing_years) && cfg$recompute_mode != "none")
  warning("Years still missing in ", pred_dir, ": ",
          paste(missing_years, collapse = ", "), call. = FALSE)

# ================================================================
# DONE
# ================================================================

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("\n================================================================\n")
cat(sprintf("National model complete in %.1f min\n", as.numeric(elapsed)))
cat("================================================================\n")
cat("Rerun downstream: 20 -> 30 -> 40 -> 41 -> 45 -> 50 -> 51 -> 60/61\n")