# ================================================================
# 01_utils.R
# ----------------------------------------------------------------
# Helper functions. Source 00_config.R first.
#
#   A. Setup            ensure_output_dirs()
#   B. CRS              convert_lat_to_utm_y(), get_zone_boundaries()
#   C. Recoding         recode_to_3class(), assign_zone_from_latitude()
#   D. Clipping         clip_bands_to_limits(), clip_raster_to_limits()
#   E. Data I/O         load_training_df(), extract_training_features(),
#                       load_composite(), prepare_training_data()
#   F. Buffer strategy  estimate_autocorr_range(), resolve_cv_buffer()
#   G. Spatial CV       create_spatial_cv_folds(), get_national_folds()
#   H. Modelling        calculate_class_weights(), fit_ranger_rf(),
#                       fit_rf_weighted(), align_probs(), predict_probs(),
#                       argmax_class(), predict_ensemble_raster(),
#                       read_selected_features(), load_national_model()
#   I. Metrics          per_class_prf(), cohens_kappa(), multi_log_loss(),
#                       fold_metrics(), summarise_fold_metrics(),
#                       run_mcnemar()
#   J. Mining polygons  rasterize_polygons(), class_area_by_polygon()
# ================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(terra, sf, dplyr, readr, ranger, blockCV, automap)

# ================================================================
# A. SETUP
# ================================================================

ensure_output_dirs <- function(root,
                               subdirs = c("models", "predictions",
                                           "statistics", "metrics")) {
  # Figures/tables are not created here: they go to P5_RESULTS
  # (results_fig_dir(), results_tab_dir() in 00_config.R).
  dir.create(root, showWarnings = FALSE, recursive = TRUE)
  for (sd in subdirs)
    dir.create(file.path(root, sd), showWarnings = FALSE, recursive = TRUE)
  invisible(root)
}

# ================================================================
# B. CRS
# ================================================================

convert_lat_to_utm_y <- function(lat, reference_lon = GHANA_REF_LON) {
  pt     <- sf::st_sfc(sf::st_point(c(reference_lon, lat)), crs = 4326)
  pt_utm <- sf::st_transform(pt, GHANA_CRS)
  unname(sf::st_coordinates(pt_utm)[, "Y"])
}

get_zone_boundaries <- function() {
  list(
    degrees = list(northern = NORTHERN_LAT, coastal = COASTAL_LAT),
    utm     = list(northern = convert_lat_to_utm_y(NORTHERN_LAT),
                   coastal  = convert_lat_to_utm_y(COASTAL_LAT))
  )
}

# ================================================================
# C. RECODING
# ================================================================

recode_to_3class <- function(codes) {
  codes <- as.character(codes)
  out   <- RAW_TO_3CLASS[codes]
  out[is.na(codes)]               <- NA_integer_
  out[!is.na(codes) & is.na(out)] <- 3L     # default -> Vegetation
  unname(as.integer(out))
}

assign_zone_from_latitude <- function(lat) {
  ifelse(lat >= NORTHERN_LAT, "northern_savannah",
         ifelse(lat >= COASTAL_LAT, "forest_transition", "coastal"))
}

# ================================================================
# D. BAND CLIPPING
# ================================================================

clip_bands_to_limits <- function(df, limits = BAND_LIMITS_GLOBAL) {
  for (bn in names(limits)) {
    if (bn %in% names(df)) {
      lims     <- limits[[bn]]
      df[[bn]] <- pmax(lims[1], pmin(lims[2], df[[bn]]))
    }
  }
  df
}

clip_raster_to_limits <- function(r, limits = BAND_LIMITS_GLOBAL) {
  for (bn in names(limits)) {
    if (bn %in% names(r)) {
      lims    <- limits[[bn]]
      r[[bn]] <- terra::clamp(r[[bn]], lower = lims[1], upper = lims[2],
                              values = TRUE)
    }
  }
  r
}

# ================================================================
# E. DATA I/O
# ================================================================

# Load the training CSV, recode to 3 classes, assign latitude zones.
load_training_df <- function(path = PATHS$training_data) {
  raw <- readr::read_csv(path, show_col_types = FALSE)
  
  lat_col <- intersect(c("latitude",  "Latitude",  "LAT", "lat"), names(raw))[1]
  lon_col <- intersect(c("longitude", "Longitude", "LON", "lon"), names(raw))[1]
  lc_col  <- intersect(c("landcover_code", "LandcoverCode", "landcover"),
                       names(raw))[1]
  if (is.na(lat_col) || is.na(lon_col) || is.na(lc_col))
    stop("Training CSV is missing latitude/longitude/landcover columns. ",
         "Found: ", paste(names(raw), collapse = ", "))
  
  df <- data.frame(
    row_id           = seq_len(nrow(raw)),
    latitude         = as.numeric(raw[[lat_col]]),
    longitude        = as.numeric(raw[[lon_col]]),
    landcover_code   = as.integer(raw[[lc_col]]),
    stringsAsFactors = FALSE
  )
  df$landcover_3class <- recode_to_3class(df$landcover_code)
  df$zone             <- assign_zone_from_latitude(df$latitude)
  df[!is.na(df$landcover_3class) & !is.na(df$latitude) &
       !is.na(df$longitude), , drop = FALSE]
}

# Extract raster values at the training points and bind them on.
extract_training_features <- function(training_df, raster, crs = GHANA_CRS) {
  training_sf <- sf::st_as_sf(training_df,
                              coords = c("longitude", "latitude"),
                              crs = 4326) |>
    sf::st_transform(crs)
  vals <- terra::extract(raster, training_sf, ID = FALSE)
  out  <- cbind(training_df, vals)
  out  <- clip_bands_to_limits(out)
  xy   <- sf::st_coordinates(training_sf)
  out$x_utm <- xy[, 1]
  out$y_utm <- xy[, 2]
  out
}

# Open a composite, name its bands and clip the requested ones.
# Stops if the layer count does not match BAND_NAMES, because names
# could then not be assigned safely.
load_composite <- function(path, bands = BAND_NAMES) {
  if (!file.exists(path)) stop("Composite not found: ", path, call. = FALSE)
  r <- terra::rast(path)
  if (terra::nlyr(r) != length(BAND_NAMES))
    stop(sprintf("%s has %d layers, BAND_NAMES has %d.",
                 basename(path), terra::nlyr(r), length(BAND_NAMES)),
         call. = FALSE)
  names(r) <- BAND_NAMES
  clip_raster_to_limits(r[[bands]])
}

# Training table used by 03, 10, 41, 44 and 50. Rows with NA in any of the
# `required` bands are dropped, so every script sees the same rows
# (and the shared fold cache stays valid).
prepare_training_data <- function(composite_path = NATIONAL_CONFIG$train_composite,
                                  required       = CANDIDATE_BANDS,
                                  verbose        = TRUE) {
  comp <- load_composite(composite_path)
  df   <- load_training_df() |>
    extract_training_features(raster = comp, crs = GHANA_CRS)
  keep <- stats::complete.cases(df[, required, drop = FALSE])
  if (verbose)
    cat(sprintf("  training points: %d (dropped %d with NA in candidate bands)\n",
                sum(keep), sum(!keep)))
  df[keep, , drop = FALSE]
}

# ================================================================
# F. BUFFER STRATEGY
# ================================================================

# Spatial autocorrelation range of `variable` over the given points.
# Returns list(range_m = buffer distance, ...). Falls back to twice the
# mean nearest-neighbour distance if the variogram fit fails.
estimate_autocorr_range <- function(
    data_sf, variable,
    cutoff_frac = CV_CONFIG$variogram_cutoff_frac,
    min_buffer  = CV_CONFIG$buffer_min_m,
    max_buffer  = CV_CONFIG$buffer_max_m,
    verbose     = TRUE
) {
  stopifnot(variable %in% names(data_sf))
  data_sf <- data_sf[!is.na(data_sf[[variable]]), ]
  if (nrow(data_sf) < 30) {
    if (verbose) message("    [vgm] <30 usable points; using min_buffer.")
    return(list(range_m = min_buffer, source = "fallback_min_buffer",
                fit = NULL, emp = NULL))
  }
  
  bbox     <- sf::st_bbox(data_sf)
  extent_m <- max(bbox[3] - bbox[1], bbox[4] - bbox[2])
  cutoff   <- min(max_buffer * 2.5, extent_m * cutoff_frac)
  
  vgm_fit <- tryCatch(
    automap::autofitVariogram(
      stats::as.formula(paste(variable, "~ 1")),
      input_data = data_sf,
      cutoff     = cutoff,
      model      = c("Sph", "Exp", "Gau", "Mat")
    ),
    error = function(e) {
      if (verbose) message("    [vgm] autofit failed: ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(vgm_fit) || is.null(vgm_fit$var_model)) {
    nn <- sf::st_nearest_feature(data_sf)
    nn_mean <- mean(as.numeric(sf::st_distance(data_sf, data_sf[nn, ],
                                               by_element = TRUE)))
    r <- max(min_buffer, min(max_buffer, round(nn_mean * 2)))
    return(list(range_m = r, source = "fallback_nn_distance",
                fit = NULL, emp = NULL))
  }
  
  model_type <- as.character(vgm_fit$var_model$model[2])
  raw_range  <- vgm_fit$var_model$range[2]
  # 95% practical range (spherical range is already practical)
  practical_range <- switch(model_type,
                            "Exp" = 3 * raw_range,
                            "Gau" = sqrt(3) * raw_range,
                            "Sph" = raw_range,
                            "Mat" = 4 * raw_range,
                            raw_range)
  r <- max(min_buffer, min(max_buffer, round(practical_range)))
  
  if (verbose)
    cat(sprintf("    [vgm] %s model; raw range %.0f m; practical %.0f m; buffer %.0f m\n",
                model_type, raw_range, practical_range, r))
  
  list(range_m = r, practical_range = practical_range, raw_range = raw_range,
       model = model_type, nugget = vgm_fit$var_model$psill[1],
       sill = sum(vgm_fit$var_model$psill), source = "variogram_fit",
       fit = vgm_fit, emp = vgm_fit$exp_var)
}

# CV buffer for a point set: fixed, or variogram of the mining
# indicator (NDVI if there are fewer than variogram_min_mining_n
# mining points).
resolve_cv_buffer <- function(data, data_sf,
                              strategy = CV_CONFIG$buffer_strategy,
                              verbose  = TRUE) {
  if (strategy == "fixed")
    return(list(buffer_m = CV_CONFIG$buffer_fixed_m, vgm_info = NULL,
                variogram_target = NA_character_))
  if (strategy != "variogram") stop("Unknown buffer strategy: ", strategy)
  
  data_sf$mining_binary <- as.numeric(data$landcover_3class == CLASS_IDS[["Mining"]])
  if ("NDVI" %in% names(data)) data_sf$NDVI <- data$NDVI
  
  target <- if (sum(data_sf$mining_binary, na.rm = TRUE) >=
                CV_CONFIG$variogram_min_mining_n) "mining_binary" else "NDVI"
  if (!target %in% names(data_sf))
    stop("Variogram target '", target, "' not found.")
  
  vgm_info <- estimate_autocorr_range(data_sf, variable = target,
                                      verbose = verbose)
  list(buffer_m = vgm_info$range_m, vgm_info = vgm_info,
       variogram_target = target)
}

# ================================================================
# G. SPATIAL CV
# ================================================================

# Spatial blocks (blockCV) plus, per fold, the training indices that
# lie outside `buffer_distance` of every test point.
create_spatial_cv_folds <- function(training_sf,
                                    n_folds         = CV_CONFIG$n_folds,
                                    buffer_distance = CV_CONFIG$buffer_min_m,
                                    block_size      = NULL,
                                    seed            = CV_CONFIG$seed) {
  set.seed(seed)
  sf::sf_use_s2(FALSE)
  
  if (is.null(block_size)) {
    bbox        <- sf::st_bbox(training_sf)
    study_range <- min(bbox[3] - bbox[1], bbox[4] - bbox[2])
    block_size  <- max(CV_CONFIG$min_block_size,
                       study_range / CV_CONFIG$block_divisor)
  }
  
  sb <- tryCatch(
    blockCV::cv_spatial(
      x = training_sf, column = "landcover_3class",
      k = n_folds, size = block_size,
      hexagon = FALSE, selection = "random",
      iteration = 100, progress = FALSE, plot = FALSE, report = FALSE
    ),
    error = function(e) {
      warning("blockCV failed; falling back to k-means on coordinates: ",
              conditionMessage(e))
      km <- stats::kmeans(sf::st_coordinates(training_sf),
                          centers = n_folds, nstart = 20)
      list(folds_ids = km$cluster)
    }
  )
  
  fold_ids <- as.integer(sb$folds_ids)
  buffer_masks <- lapply(seq_len(n_folds), function(k) {
    test_idx  <- which(fold_ids == k)
    train_idx <- which(fold_ids != k)
    if (buffer_distance <= 0 || !length(test_idx) || !length(train_idx))
      return(train_idx)
    near <- lengths(sf::st_is_within_distance(training_sf[train_idx, ],
                                              training_sf[test_idx, ],
                                              dist = buffer_distance)) > 0
    train_idx[!near]
  })
  
  list(fold_ids = fold_ids, buffer_masks = buffer_masks,
       block_size = block_size, buffer_distance = buffer_distance)
}

# One national set of spatial folds, cached on disk and reused by 03,10,
# 41, and 50 so their OOF predictions are paired. The cache is
# rebuilt if the training rows, seed, fold count or buffer strategy
# change.
get_national_folds <- function(df,
                               cache_path = CV_CONFIG$fold_cache,
                               rebuild    = FALSE,
                               verbose    = TRUE) {
  key <- list(row_id   = df$row_id,
              seed     = CV_CONFIG$seed,
              n_folds  = CV_CONFIG$n_folds,
              strategy = CV_CONFIG$buffer_strategy,
              fixed_m  = CV_CONFIG$buffer_fixed_m)
  if (!rebuild && file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (identical(cached$key, key)) {
      if (verbose)
        cat(sprintf("  shared folds loaded from cache (buffer %.0f m, block %.0f m)\n",
                    cached$buffer_m, cached$block_size))
      return(cached)
    }
    if (verbose) cat("  fold cache is out of date; rebuilding\n")
  }
  
  data_sf <- sf::st_as_sf(df, coords = c("x_utm", "y_utm"),
                          crs = GHANA_CRS, remove = FALSE)
  if (verbose) cat("  resolving CV buffer...\n")
  br <- resolve_cv_buffer(df, data_sf, verbose = verbose)
  sb <- create_spatial_cv_folds(data_sf, buffer_distance = br$buffer_m)
  
  out <- list(key              = key,
              row_id           = df$row_id,
              fold_ids         = sb$fold_ids,
              buffer_masks     = sb$buffer_masks,
              buffer_m         = br$buffer_m,
              block_size       = sb$block_size,
              vgm_info         = br$vgm_info,
              variogram_target = br$variogram_target,
              created          = Sys.time())
  dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(out, cache_path)
  if (verbose)
    cat(sprintf("  shared folds built: buffer %.0f m, block %.0f m, sizes %s\n",
                out$buffer_m, out$block_size,
                paste(table(out$fold_ids), collapse = " | ")))
  out
}

# ================================================================
# H. MODELLING
# ================================================================

# Square-root inverse-frequency class weights, in the order of the
# factor levels of y.
calculate_class_weights <- function(y) {
  freq      <- table(y)
  n_classes <- length(freq)
  w         <- sqrt(length(y) / (n_classes * freq))
  w         <- w / sum(w) * n_classes
  as.numeric(w)
}

fit_ranger_rf <- function(x, y, class_weights = NULL,
                          num_trees     = 500,
                          min_node_size = 3,
                          sample_frac   = 0.8,
                          num_threads   = 1,
                          seed          = NULL,
                          importance    = "permutation") {
  ranger::ranger(
    x               = x,
    y               = y,
    num.trees       = num_trees,
    mtry            = max(1, floor(sqrt(ncol(x)))),
    min.node.size   = min_node_size,
    sample.fraction = sample_frac,
    importance      = importance,
    probability     = TRUE,
    class.weights   = class_weights,
    num.threads     = num_threads,
    seed            = seed,
    verbose         = FALSE
  )
}

# Weighted probability forest with the project's settings. Classes
# absent from y are dropped before the weights are computed, so the
# weights always match ranger's factor levels.
fit_rf_weighted <- function(X, y, seed,
                            sample_frac   = NATIONAL_CONFIG$cv_sample_frac,
                            threads       = 1,
                            num_trees     = NATIONAL_CONFIG$num_trees,
                            min_node_size = NATIONAL_CONFIG$min_node_size,
                            importance    = "permutation") {
  y <- droplevels(factor(y, levels = CLASS_LEVELS))
  fit_ranger_rf(X, y,
                class_weights = calculate_class_weights(y),
                num_trees     = num_trees,
                min_node_size = min_node_size,
                sample_frac   = sample_frac,
                num_threads   = threads,
                seed          = seed,
                importance    = importance)
}

# Probability matrix with all three classes in fixed order (columns
# "1","2","3"). ranger only returns columns for classes present in the
# training data.
align_probs <- function(p) {
  out  <- matrix(0, nrow(p), length(CLASS_LEVELS),
                 dimnames = list(NULL, CLASS_LEVELS))
  cols <- intersect(colnames(p), CLASS_LEVELS)
  out[, cols] <- p[, cols, drop = FALSE]
  out
}

predict_probs <- function(model, X)
  align_probs(predict(model, X, num.threads = 1)$predictions)

argmax_class <- function(p)
  as.integer(CLASS_LEVELS)[max.col(p, ties.method = "first")]

# Ensemble-averaged class probabilities over a raster. Returns
#   pred   integer class (layer index == class id)
#   uncert normalised Shannon entropy in [0, 1]
#   prob   3-band probability raster
predict_ensemble_raster <- function(composite, models, features) {
  r <- composite[[features]]
  prob_fun <- function(model, data, ...)
    align_probs(predict(model, data, num.threads = 1)$predictions)
  
  prob_sum <- NULL
  for (m in models) {
    p <- terra::predict(r, m, fun = prob_fun, na.rm = TRUE)
    prob_sum <- if (is.null(prob_sum)) p else prob_sum + p
    rm(p); gc(verbose = FALSE)
  }
  prob <- prob_sum / length(models)
  names(prob) <- paste0("prob_class_", CLASS_LEVELS)
  
  cls <- terra::which.max(prob)
  names(cls) <- "class"
  ent <- -sum(prob * log(prob + 1e-10)) / log(terra::nlyr(prob))
  names(ent) <- "entropy"
  list(pred = cls, uncert = ent, prob = prob)
}

# National band set: frozen vector in NATIONAL_CONFIG$features, else
# the set written by the feature audit (03).
read_selected_features <- function(cfg = NATIONAL_CONFIG) {
  if (!is.null(cfg$features)) {
    feats <- cfg$features
    src   <- "frozen list in NATIONAL_CONFIG$features"
  } else {
    if (!file.exists(cfg$selected_file))
      stop("No band set: NATIONAL_CONFIG$features is NULL and ",
           cfg$selected_file, " is missing. Run 03 first.", call. = FALSE)
    sel   <- readr::read_csv(cfg$selected_file, show_col_types = FALSE)
    sel   <- sel[order(sel$order), , drop = FALSE]
    feats <- sel$feature
    src   <- sprintf("%s (%s rule, from 03)", basename(cfg$selected_file),
                     sel$rule[1])
  }
  feats <- unique(feats)
  bad <- intersect(feats, HOUSEKEEPING_BANDS)
  if (length(bad)) stop("Quality bands in the feature set: ",
                        paste(bad, collapse = ", "), call. = FALSE)
  miss <- setdiff(feats, BAND_NAMES)
  if (length(miss)) stop("Features not in BAND_NAMES: ",
                         paste(miss, collapse = ", "), call. = FALSE)
  unstable <- intersect(feats, TEMPORALLY_UNSTABLE_BANDS)
  if (length(unstable))
    warning("Temporally unstable bands selected (", paste(unstable, collapse = ", "),
            "): the nightlight sensor changes in 2013, which biases the ",
            "2000-2012 back-cast.", call. = FALSE)
  if (length(feats) < 3) stop("Fewer than 3 features selected.", call. = FALSE)
  list(features = feats, source = src)
}

load_national_model <- function(path = NATIONAL_CONFIG$model_path) {
  if (!file.exists(path))
    stop("National model not found: ", path, ". Run 10 first.", call. = FALSE)
  readRDS(path)
}

# ================================================================
# I. METRICS
# ================================================================

per_class_prf <- function(y_true, y_pred, levels = 1:3) {
  cm <- table(pred = factor(y_pred, levels = levels),
              true = factor(y_true, levels = levels))
  out <- data.frame(class_id = levels, class = CLASS_NAMES[levels],
                    precision = NA_real_, recall = NA_real_, f1 = NA_real_)
  for (i in seq_along(levels)) {
    tp <- cm[i, i]
    fp <- sum(cm[i, ]) - tp
    fn <- sum(cm[, i]) - tp
    out$precision[i] <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    out$recall[i]    <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    p <- out$precision[i]; r <- out$recall[i]
    out$f1[i] <- if (!is.na(p) && !is.na(r) && (p + r) > 0)
      2 * p * r / (p + r) else NA_real_
  }
  list(cm = cm, prf = out)
}

cohens_kappa <- function(cm) {
  n <- sum(cm); if (n == 0) return(NA_real_)
  po <- sum(diag(cm)) / n
  pe <- sum((rowSums(cm) / n) * (colSums(cm) / n))
  (po - pe) / (1 - pe)
}

# Multi-class negative log-likelihood. truth in 1..ncol(prob).
multi_log_loss <- function(truth, prob, eps = 1e-15) {
  p_true <- prob[cbind(seq_along(truth), as.integer(truth))]
  -mean(log(pmax(p_true, eps)), na.rm = TRUE)
}

# Per-fold metrics. Rows without a prediction are ignored; folds with
# fewer than min_n rows are dropped. prob (n x 3) adds log-loss.
fold_metrics <- function(y_true, y_pred, fold, prob = NULL, min_n = 10) {
  ok <- !is.na(y_pred) & !is.na(y_true)
  rows <- lapply(sort(unique(fold[ok])), function(k) {
    i <- which(ok & fold == k)
    if (length(i) < min_n) return(NULL)
    pc <- per_class_prf(y_true[i], y_pred[i])
    data.frame(fold          = k,
               n_test        = length(i),
               accuracy      = mean(y_pred[i] == y_true[i]),
               kappa         = cohens_kappa(pc$cm),
               water_f1      = pc$prf$f1[1],
               mining_f1     = pc$prf$f1[2],
               vegetation_f1 = pc$prf$f1[3],
               log_loss      = if (is.null(prob)) NA_real_ else
                 multi_log_loss(y_true[i], prob[i, , drop = FALSE]))
  })
  do.call(rbind, rows)
}

# Fold summary with t-distribution 95% CI for accuracy. NULL if empty.
summarise_fold_metrics <- function(fold_df) {
  if (is.null(fold_df) || nrow(fold_df) == 0) return(NULL)
  n      <- nrow(fold_df)
  t_crit <- stats::qt(0.975, df = max(n - 1, 1))
  acc_m  <- mean(fold_df$accuracy, na.rm = TRUE)
  acc_sd <- stats::sd(fold_df$accuracy, na.rm = TRUE)
  se     <- acc_sd / sqrt(n)
  data.frame(
    n_folds            = n,
    n_test             = if ("n_test" %in% names(fold_df)) sum(fold_df$n_test) else NA,
    mean_accuracy      = acc_m,
    sd_accuracy        = acc_sd,
    ci95_lower         = acc_m - t_crit * se,
    ci95_upper         = acc_m + t_crit * se,
    mean_kappa         = mean(fold_df$kappa,         na.rm = TRUE),
    sd_kappa           = stats::sd(fold_df$kappa,    na.rm = TRUE),
    mean_water_f1      = mean(fold_df$water_f1,      na.rm = TRUE),
    mean_mining_f1     = mean(fold_df$mining_f1,     na.rm = TRUE),
    sd_mining_f1       = stats::sd(fold_df$mining_f1, na.rm = TRUE),
    mean_vegetation_f1 = mean(fold_df$vegetation_f1, na.rm = TRUE),
    mean_log_loss      = if ("log_loss" %in% names(fold_df))
      mean(fold_df$log_loss, na.rm = TRUE) else NA_real_
  )
}

# Paired McNemar test on two OOF prediction vectors. Exact binomial
# test when there are few discordant pairs.
run_mcnemar <- function(pred_a, pred_b, truth, label_a, label_b,
                        exact_threshold = 25) {
  ok <- !is.na(pred_a) & !is.na(pred_b) & !is.na(truth)
  a_correct <- pred_a[ok] == truth[ok]
  b_correct <- pred_b[ok] == truth[ok]
  n_a_only  <- sum( a_correct & !b_correct)    # A right, B wrong
  n_b_only  <- sum(!a_correct &  b_correct)    # B right, A wrong
  n_disc    <- n_a_only + n_b_only
  
  if (n_disc == 0) {
    chi_sq <- NA_real_; p_val <- NA_real_; test_used <- "none (no discordant pairs)"
  } else if (n_disc < exact_threshold) {
    chi_sq <- NA_real_
    p_val  <- stats::binom.test(n_a_only, n_disc, p = 0.5)$p.value
    test_used <- "exact binomial"
  } else {
    chi_sq <- (abs(n_a_only - n_b_only) - 1)^2 / n_disc
    p_val  <- stats::pchisq(chi_sq, df = 1, lower.tail = FALSE)
    test_used <- "McNemar (continuity-corrected)"
  }
  
  data.frame(
    comparison      = paste(label_a, "vs", label_b),
    n_paired        = sum(ok),
    n_discordant    = n_disc,
    a_only_correct  = n_a_only,
    b_only_correct  = n_b_only,
    acc_a           = mean(a_correct),
    acc_b           = mean(b_correct),
    chi_squared     = chi_sq,
    p_value         = p_val,
    test_used       = test_used,
    significant_005 = !is.na(p_val) & p_val < 0.05,
    better          = if (n_a_only > n_b_only) label_a
    else if (n_b_only > n_a_only) label_b else "tie",
    stringsAsFactors = FALSE
  )
}

# ================================================================
# J. MINING POLYGONS
# ================================================================

# Rasterise the Maus (2022) polygons onto a reference grid. A cell is
# inside if its centre falls in a polygon (touches = FALSE).
rasterize_polygons <- function(polygons, reference_raster, touches = FALSE) {
  if (inherits(polygons, "character"))
    polygons <- sf::st_read(polygons, quiet = TRUE)
  polygons <- sf::st_transform(polygons, terra::crs(reference_raster))
  r <- terra::rasterize(terra::vect(polygons), reference_raster,
                        field = 1, background = 0, touches = touches)
  names(r) <- "inside"
  r
}

# Per-class area (ha) split into inside / outside a polygon mask.
class_area_by_polygon <- function(prediction, inside_mask,
                                  cell_area_ha = NULL,
                                  classes      = seq_along(CLASS_NAMES)) {
  if (is.null(cell_area_ha))
    cell_area_ha <- as.numeric(
      terra::global(terra::cellSize(prediction, unit = "ha"),
                    "mean", na.rm = TRUE)[1, 1])
  inside_bool <- inside_mask == 1L
  
  rows <- lapply(classes, function(cls) {
    cls_mask  <- prediction == cls
    n_total   <- as.numeric(terra::global(cls_mask, "sum", na.rm = TRUE)[1, 1])
    n_inside  <- as.numeric(terra::global(cls_mask & inside_bool, "sum",
                                          na.rm = TRUE)[1, 1])
    n_outside <- n_total - n_inside
    data.frame(class           = as.integer(cls),
               count_total     = as.integer(n_total),
               count_inside    = as.integer(n_inside),
               count_outside   = as.integer(n_outside),
               area_ha_total   = n_total   * cell_area_ha,
               area_ha_inside  = n_inside  * cell_area_ha,
               area_ha_outside = n_outside * cell_area_ha,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}