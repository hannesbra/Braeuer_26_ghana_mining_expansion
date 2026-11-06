# ================================================================
# Transition Constraints
# ----------------------------------------------------------------
#
# Rule: going backward in time from reference year, mining and
# water extents can only shrink. A pixel without data in a year keeps
# that year NA; the rule is carried across the gap using the last
# observed later class (same rule as the simulations in 30).
#
# Script goes back from REFERENCE_YEAR, writes constrained rasters to PATHS$out_constrained, and emits:
#   predictions/prediction_{year}_shrinking.tif    (shrunk class)
#   predictions/uncertainty_{year}_shrinking.tif   (copy of P1's)
#   statistics/constraint_application_stats.csv   (pixels changed / year)
#   statistics/area_before_after.csv              (ha per class, P1 vs P2)
#   models/shrinking_config.rds                   (audit trail)
# ================================================================

source("R/00_config.R")
source("R/01_utils.R")

pacman::p_load(terra, dplyr, readr, jsonlite)

terraOptions(memfrac = 0.6, tempdir = tempdir())
start_time <- Sys.time()

# ================================================================
# 0. Transition Constraints CONFIG
# ----------------------------------------------------------------
# ================================================================

if (!exists("Transition_CONFIG")) {
  SHRINKING_CONFIG <- list(
    reference_year = PREDICT_CONFIG$base_year,        # 2019
    
    # "if future is class X, past is allowed to be in this set"
    allowed_past_from_future = list(
      "1" = c(1L, 2L, 3L),   # future = Water       → any past
      "2" = c(1L, 2L, 3L),   # future = Mining      → any past
      "3" = c(3L)            # future = Vegetation  → past must be Vegetation
    ),
    
    # When a past pixel violates the rule, replace it with this class
    default_correction = 3L,
    
    # Written filename suffix — keeps P2 outputs visually distinct
    filename_suffix = "_shrinking"
  )
}

# ================================================================
# 1. SETUP
# ================================================================

in_dir  <- PATHS$national_predictions   # national model (10)
out_dir <- PATHS$out_constrained
ensure_output_dirs(out_dir)

cat("================================================================\n")
cat("P2 — Transition constraints POST-PROCESSING\n")
cat("================================================================\n")
cat(sprintf("P1 input :  %s\n", in_dir))
cat(sprintf("P2 output:  %s\n", out_dir))
cat(sprintf("Reference year (no constraint): %d\n",
            SHRINKING_CONFIG$reference_year))
cat("\nRules (future → allowed past):\n")
for (fc in names(SHRINKING_CONFIG$allowed_past_from_future)) {
  cat(sprintf("  future=%s → past in {%s}\n",
              CLASS_NAMES[as.integer(fc)],
              paste(CLASS_NAMES[SHRINKING_CONFIG$allowed_past_from_future[[fc]]],
                    collapse = ", ")))
}
cat("\n")

# Config as an audit trail
saveRDS(SHRINKING_CONFIG,
        file.path(out_dir, "models", "shrinking_config.rds"))
jsonlite::write_json(SHRINKING_CONFIG,
                     file.path(out_dir, "models", "shrinking_config.json"),
                     auto_unbox = TRUE, pretty = TRUE)

# ================================================================
# 2. CORE FUNCTION
# ----------------------------------------------------------------

shrink_one_year <- function(past_r, future_r, cfg) {
  # Ensure the rasters line up; the prediction engine produces identical geometry
  # across years, guarding against mismatched extents/resolutions.
  if (!compareGeom(past_r, future_r, stopOnError = FALSE)) {
    future_r <- resample(future_r, past_r, method = "near")
  }
  
  needs_correction <- NULL
  
  for (fc_str in names(cfg$allowed_past_from_future)) {
    allowed <- as.integer(cfg$allowed_past_from_future[[fc_str]])
    fc      <- as.integer(fc_str)
    
    # No-op rules (any past allowed)
    if (setequal(allowed, 1:3)) next
    
    # past ∈ allowed set — express as OR of equality checks (portable)
    past_ok <- Reduce(`|`, lapply(allowed, function(v) past_r == v))
    
    # Violation mask: future is this class, and past is not allowed
    viol <- (future_r == fc) & (!past_ok)
    # Missing future or missing past: no correction (past stays as it is)
    viol <- terra::ifel(is.na(viol), 0, viol)
    
    needs_correction <- if (is.null(needs_correction)) viol
    else needs_correction | viol
  }
  
  if (is.null(needs_correction)) {
    # All rules were no-ops (e.g. config allows anything)
    return(list(corrected = past_r, n_changed = 0L))
  }
  
  corrected <- ifel(needs_correction, cfg$default_correction, past_r)
  n_changed <- as.integer(global(needs_correction, "sum", na.rm = TRUE)[1, 1])
  
  list(corrected = corrected, n_changed = n_changed)
}

# ================================================================
# 3. REVERSE-TIME LOOP
# ================================================================

years        <- sort(PREDICT_CONFIG$predict_years)
years_desc   <- sort(years, decreasing = TRUE)
ref_year     <- SHRINKING_CONFIG$reference_year
stats_rows   <- list()
future_rast  <- NULL

cat("[1/3] Applying transition constraints (reverse chronological)\n")

for (yr in years_desc) {
  in_path <- file.path(in_dir, sprintf("prediction_%d.tif", yr))
  if (!file.exists(in_path)) {
    warning(sprintf("Missing base prediction for %d — skipping", yr))
    next
  }
  
  current     <- rast(in_path)
  out_pred    <- file.path(out_dir, "predictions",
                           sprintf("prediction_%d%s.tif",
                                   yr, SHRINKING_CONFIG$filename_suffix))
  out_uncert  <- file.path(out_dir, "predictions",
                           sprintf("uncertainty_%d%s.tif",
                                   yr, SHRINKING_CONFIG$filename_suffix))
  
  if (yr == ref_year) {
    # Reference year: no constraint, just write through
    writeRaster(current, out_pred, overwrite = TRUE,
                datatype = "INT1U", NAflag = 255)
    future_rast <- current
    stats_rows[[as.character(yr)]] <- data.frame(
      year             = yr,
      reference_year   = NA_integer_,
      n_pixels_changed = 0L,
      change_rate_pct  = 0,
      is_reference     = TRUE
    )
    cat(sprintf("  %d (reference): written through\n", yr))
  } else {
    res <- shrink_one_year(current, future_rast, SHRINKING_CONFIG)
    writeRaster(res$corrected, out_pred, overwrite = TRUE,
                datatype = "INT1U", NAflag = 255)
    
    # Rate computed against non-NA pixels of the unconstrained raster
    n_valid <- as.numeric(global(!is.na(current), "sum", na.rm = TRUE)[1, 1])
    rate <- if (n_valid > 0) 100 * res$n_changed / n_valid else NA_real_
    
    stats_rows[[as.character(yr)]] <- data.frame(
      year             = yr,
      reference_year   = ref_year,
      n_pixels_changed = res$n_changed,
      change_rate_pct  = rate,
      is_reference     = FALSE
    )
    cat(sprintf("  %d: %s pixels modified (%.2f%% of valid pixels)\n",
                yr, format(res$n_changed, big.mark = ","), rate))
    
    # Gaps: where this year has no data, keep the last observed later
    # class, so the rule carries across the gap (same rule as 30).
    future_rast <- terra::cover(res$corrected, future_rast)
    rm(res); gc(verbose = FALSE)
  }
  
  # Copy uncertainty raster unchanged
  u_in <- file.path(in_dir, sprintf("uncertainty_%d.tif", yr))
  if (file.exists(u_in)) {        # always refresh: rasters may come from a new model
    file.copy(u_in, out_uncert, overwrite = TRUE)
  }
  
  rm(current); gc(verbose = FALSE)
}

# ================================================================
# 4. CONSTRAINT-APPLICATION STATS
# ================================================================

stats_df <- do.call(rbind, stats_rows)
stats_df <- stats_df[order(stats_df$year), , drop = FALSE]
stats_path <- file.path(out_dir, "statistics", "constraint_application_stats.csv")
write_csv(stats_df, stats_path)

cat(sprintf("\n[2/3] Transition Constraint stats → %s\n", stats_path))
cat(sprintf("  Total pixels changed across %d years: %s\n",
            nrow(stats_df) - 1L,
            format(sum(stats_df$n_pixels_changed), big.mark = ",")))
cat(sprintf("  Peak change rate: %.2f%% (year %d)\n",
            max(stats_df$change_rate_pct, na.rm = TRUE),
            stats_df$year[which.max(stats_df$change_rate_pct)]))

# ================================================================
# 5. BEFORE/AFTER AREA COMPARISON
# ================================================================
# Per-year class area (ha) before and after transition constraints.

cat("\n[3/3] Computing before/after area comparison\n")

compute_class_area_ha <- function(r) {
  cell_ha <- as.numeric(global(cellSize(r, unit = "ha"),
                               "mean", na.rm = TRUE)[1, 1])
  fr <- terra::freq(r)
  if (is.null(fr) || nrow(fr) == 0) return(NULL)
  fr$area_ha <- fr$count * cell_ha
  # freq() emits column "value" for the class code in terra ≥ 1.7
  if ("value" %in% names(fr)) names(fr)[names(fr) == "value"] <- "class"
  fr[, c("class", "count", "area_ha")]
}

compare_rows <- list()
for (yr in years) {
  p1_path <- file.path(in_dir,  sprintf("prediction_%d.tif", yr))
  p2_path <- file.path(out_dir, "predictions",
                       sprintf("prediction_%d%s.tif",
                               yr, SHRINKING_CONFIG$filename_suffix))
  if (!file.exists(p1_path) || !file.exists(p2_path)) next
  
  a_before <- compute_class_area_ha(rast(p1_path))
  a_after  <- compute_class_area_ha(rast(p2_path))
  if (is.null(a_before) || is.null(a_after)) next
  
  names(a_before) <- c("class", "count_before", "area_ha_before")
  names(a_after)  <- c("class", "count_after",  "area_ha_after")
  merged <- merge(a_before, a_after, by = "class", all = TRUE)
  merged$year       <- yr
  merged$class_name <- CLASS_NAMES[merged$class]
  merged$delta_ha   <- merged$area_ha_after - merged$area_ha_before
  compare_rows[[as.character(yr)]] <- merged[, c(
    "year", "class", "class_name",
    "count_before", "area_ha_before",
    "count_after",  "area_ha_after",
    "delta_ha"
  )]
}

area_df <- do.call(rbind, compare_rows)
area_path <- file.path(out_dir, "statistics", "area_before_after.csv")
write_csv(area_df, area_path)
cat(sprintf("  Area comparison → %s\n", area_path))

# Mining-focused one-liner for the console
mining_summary <- area_df |>
  dplyr::filter(class == 2L) |>
  dplyr::arrange(year) |>
  dplyr::mutate(shrank_by_ha = round(area_ha_before - area_ha_after, 1))

cat("\n  Mining area before vs after (ha):\n")
print(mining_summary[, c("year", "area_ha_before", "area_ha_after", "shrank_by_ha")],
      row.names = FALSE, digits = 5)

# ================================================================
# 6. INTERSECTION WITH MINING POLYGONS
# ----------------------------------------------------------------
# Raw totals give us more an idea about data quality. Mining class only makes sense inside polygons
# ================================================================

poly_path <- PATHS$mining_polygons
if (!file.exists(poly_path)) {
  cat("\n[4/4] SKIPPED — no mining polygon file at ", poly_path, "\n", sep = "")
} else {
  cat("\n[4/4] Intersecting predictions with known mining polygons\n")
  
  ref_r <- rast(file.path(in_dir, sprintf("prediction_%d.tif", ref_year)))
  
  cat("  Rasterizing polygons onto prediction grid...\n")
  inside_mask <- rasterize_polygons(poly_path, ref_r)
  writeRaster(inside_mask,
              file.path(out_dir, "models", "inside_mining_polygons_mask.tif"),
              datatype = "INT1U", overwrite = TRUE)
  
  cell_ha <- as.numeric(global(cellSize(ref_r, unit = "ha"),
                               "mean", na.rm = TRUE)[1, 1])
  
  intersect_rows <- list()
  for (yr in years) {
    p1_path <- file.path(in_dir, sprintf("prediction_%d.tif", yr))
    p2_path <- file.path(out_dir, "predictions",
                         sprintf("prediction_%d%s.tif",
                                 yr, SHRINKING_CONFIG$filename_suffix))
    if (!file.exists(p1_path) || !file.exists(p2_path)) next
    
    before <- class_area_by_polygon(rast(p1_path), inside_mask, cell_ha)
    after  <- class_area_by_polygon(rast(p2_path), inside_mask, cell_ha)
    before$stage <- "unconstrained"
    after$stage  <- "shrinking"
    
    both <- rbind(before, after)
    both$year       <- yr
    both$class_name <- CLASS_NAMES[both$class]
    intersect_rows[[as.character(yr)]] <- both[, c(
      "year", "stage", "class", "class_name",
      "count_total",   "count_inside",   "count_outside",
      "area_ha_total", "area_ha_inside", "area_ha_outside"
    )]
  }
  
  intersect_df <- do.call(rbind, intersect_rows)
  intersect_path <- file.path(out_dir, "statistics",
                              "area_within_mining_polygons.csv")
  write_csv(intersect_df, intersect_path)
  cat(sprintf("  Polygon-aware area stats → %s\n", intersect_path))
  
  # ---- Mining-focused console summary ---------------------------
  mining_wide <- intersect_df |>
    dplyr::filter(class %in% c(1L, 2L)) |>            # Mining + Water (disturbed)
    dplyr::group_by(year, stage) |>
    dplyr::summarise(area_ha_inside  = sum(area_ha_inside),
                     area_ha_outside = sum(area_ha_outside), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from   = stage,
      values_from  = c(area_ha_inside, area_ha_outside),
      names_glue   = "{.value}_{stage}"
    ) |>
    dplyr::arrange(year) |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 0)))
  
  cat("\n  Mining area (ha) — split by mining-polygon overlap:\n")
  cat("    inside  = within a delineated mining polygon\n")
  cat("    outside = elsewhere\n\n")
  print(mining_wide, n = Inf)
  
  # ---- Headline numbers -----------------------------------------
  ref_row     <- mining_wide |> dplyr::filter(year == ref_year)
  earliest_yr <- min(mining_wide$year)
  early_row   <- mining_wide |> dplyr::filter(year == earliest_yr)
  
  if (nrow(ref_row) == 1 && nrow(early_row) == 1) {
    cat(sprintf(
      "\n  At reference year %d: %.0f ha inside polygons, %.0f ha outside.\n",
      ref_year,
      ref_row$area_ha_inside_unconstrained,
      ref_row$area_ha_outside_unconstrained
    ))
    cat(sprintf(
      "  At %d (unconstrained): %.0f ha inside, %.0f ha outside.\n",
      earliest_yr,
      early_row$area_ha_inside_unconstrained,
      early_row$area_ha_outside_unconstrained
    ))
    cat(sprintf(
      "  At %d (shrinking):     %.0f ha inside, %.0f ha outside.\n",
      earliest_yr,
      early_row$area_ha_inside_shrinking,
      early_row$area_ha_outside_shrinking
    ))
    
    inside_removed  <- early_row$area_ha_inside_unconstrained -
      early_row$area_ha_inside_shrinking
    outside_removed <- early_row$area_ha_outside_unconstrained -
      early_row$area_ha_outside_shrinking
    denom <- inside_removed + outside_removed
    pct_outside <- if (denom > 0) 100 * outside_removed / denom else NA
    cat(sprintf(
      "  In %d the shrinking constraint removed %.0f ha inside and %.0f ha\n  outside known mine polygons (%.1f%% of the removed area was outside).\n",
      earliest_yr, inside_removed, outside_removed, pct_outside
    ))
  }
}

# ================================================================
# DONE
# ================================================================

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("\n================================================================\n")
cat(sprintf("Transition constraints complete in %.1f min\n", as.numeric(elapsed)))
cat("================================================================\n")