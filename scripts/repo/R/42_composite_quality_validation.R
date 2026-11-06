# ================================================================
# 42 Date quality and validation
# ----------------------------------------------------------------
# Per zone × year quality diagnostics for the medoid composites.
# Looks at "could classification issues stem from input compositing rather than the model?"
#
# For each year × zone, summarises four indicators:
#   1. VALID_PIXEL_COUNT  – number of cloud-free Landsat scenes for that period
#   2. CLOUD_PROB         – mean residual cloud probability
#   3. temporal noise     – mean of NIR/RED/SWIR1 std bands
#   4. NA coverage        – fraction of pixels with no usable data
#
# Indicators gets a green / yellow / red rating against thresholds in QC_CONFIG. 
#
# Outputs (in PATHS$out_validation):
#   tables/qc_zone_year_summary.csv       extraction results, all bands
#   tables/qc_reflectance_quality.csv     saturation/extreme stats
#   tables/qc_index_quality.csv           index range/coherence
#   tables/qc_scorecard.csv               OK/CAUTION/SUSPECT per zone×year
#   P5_RESULTS/figures/42_composite_qc/qc_scorecard_heatmap.png        20-year × 3-zone overview
# ================================================================

source("R/00_config.R")
source("R/01_utils.R")

pacman::p_load(terra, sf, dplyr, tidyr, readr, ggplot2, progress)

terraOptions(memfrac = 0.6, tempdir = tempdir())
set.seed(CV_CONFIG$seed)
start_time <- Sys.time()

# ================================================================
# 0. CONFIG
# ================================================================
# Thresholds chosen from a first pass over Ghana 2000–2019 composites.
# If my scorecard ends up flagging more than ~30% of zone-years, relax thresholds or accept early Landsat 7 era 2000–2003 poorer quality

QC_CONFIG <- list(
  n_sample_pixels = 15000,        # stratified sample per run
  
  quality_bands     = c("PIXEL_COUNT", "CLOUD_PROB", "VALID_PIXEL_COUNT",
                        "NIR_STD_TEMPORAL", "RED_STD_TEMPORAL",
                        "SWIR1_STD_TEMPORAL", "DATA_MASK"),
  reflectance_bands = c("BLUE", "GREEN", "RED", "NIR", "SWIR1", "SWIR2"),
  index_bands       = c("NDVI", "BSI", "MNDWI", "NDBI"),
  
  thresholds = list(
    valid_pixel_count = list(good = 15,   moderate = 8),    # ≥
    cloud_prob_pct    = list(good = 10,   moderate = 25),   # ≤
    temporal_noise    = list(good = 0.02, moderate = 0.05), # ≤
    na_pct            = list(good = 2,    moderate = 10)    # ≤
  )
)

out_dir <- PATHS$out_validation
ensure_output_dirs(out_dir, subdirs = c("tables"))
plot_dir <- results_fig_dir("42_composite_qc")

cat("================================================================\n")
cat("COMPOSITE QC — per zone × year diagnostics\n")
cat("================================================================\n")

# ================================================================
# 1. ZONE RASTER + STRATIFIED SAMPLE
# ================================================================

cat("\n[1/4] Building zone raster + stratified sample\n")

ref_path <- file.path(PATHS$composite_dir,
                      sprintf("Mining_Composite_%d.tif",
                              PREDICT_CONFIG$base_year))
ref <- rast(ref_path)
if (nlyr(ref) >= length(BAND_NAMES))
  names(ref) <- BAND_NAMES[seq_len(nlyr(ref))]

# Boundaries in the raster's native CRS
bds <- get_zone_boundaries()
if (terra::is.lonlat(ref)) {
  north_b <- bds$degrees$northern; coast_b <- bds$degrees$coastal
  cat("  Raster CRS: geographic (degrees)\n")
} else {
  north_b <- bds$utm$northern;     coast_b <- bds$utm$coastal
  cat("  Raster CRS: projected (UTM)\n")
}

y_coords    <- init(ref[[1]], "y")
zone_raster <- ifel(y_coords >= north_b, 1L,
                    ifel(y_coords >= coast_b, 2L, 3L))
names(zone_raster) <- "zone"
rm(y_coords); gc(verbose = FALSE)

zone_id_to_key <- c("1" = "northern_savannah",
                    "2" = "forest_transition",
                    "3" = "coastal")

zf <- freq(zone_raster)
for (i in seq_len(nrow(zf))) {
  cat(sprintf("    %s: %s pixels\n",
              ZONE_LABELS[zone_id_to_key[as.character(zf$value[i])]],
              format(zf$count[i], big.mark = ",")))
}

# Stratified sample: cells + their zone IDs
sample_cells <- spatSample(zone_raster,
                           size = QC_CONFIG$n_sample_pixels,
                           method = "random", cells = TRUE, na.rm = TRUE)
sample_zones <- terra::extract(zone_raster, sample_cells$cell)$zone

cat(sprintf("  Sampled %d pixels (per zone: %s)\n",
            nrow(sample_cells),
            paste(table(sample_zones), collapse = " / ")))

# Convert numeric zone id → zone key string for cleaner downstream joins
sample_zone_keys <- zone_id_to_key[as.character(sample_zones)]

# ================================================================
# 2. EXTRACT PER YEAR × ZONE × BAND
# ================================================================

cat("\n[2/4] Extracting per-band statistics across years\n")

YEARS <- PREDICT_CONFIG$predict_years

# Helper: per-zone summary of a single band
summarise_band_by_zone <- function(values, zones, year, band) {
  rows <- list()
  for (zk in unique(zones)) {
    if (is.na(zk)) next
    z_vals <- values[zones == zk]
    z_vals <- z_vals[!is.na(z_vals)]
    if (length(z_vals) < 10) next
    
    rows[[zk]] <- data.frame(
      year     = year, zone = zk, zone_name = ZONE_LABELS[zk],
      band     = band,
      mean     = mean(z_vals),
      median   = median(z_vals),
      sd       = sd(z_vals),
      p05      = unname(quantile(z_vals, 0.05)),
      p25      = unname(quantile(z_vals, 0.25)),
      p75      = unname(quantile(z_vals, 0.75)),
      p95      = unname(quantile(z_vals, 0.95)),
      pct_zero = 100 * mean(z_vals == 0),
      n_valid  = length(z_vals),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# Helper: NA fraction inside the sampled pixels per zone
sample_na_pct_by_zone <- function(values_with_na, zones) {
  data.frame(
    zone     = unique(zones[!is.na(zones)]),
    pct_na   = sapply(unique(zones[!is.na(zones)]), function(zk)
      100 * mean(is.na(values_with_na[zones == zk]))),
    stringsAsFactors = FALSE
  )
}

quality_rows     <- list()
reflectance_rows <- list()
index_rows       <- list()
na_rows          <- list()

pb <- progress_bar$new(
  format = "  :year [:bar] :percent | :elapsed",
  total  = length(YEARS), clear = FALSE, width = 60
)

for (yr in YEARS) {
  pb$tick(tokens = list(year = yr))
  cp <- file.path(PATHS$composite_dir,
                  sprintf("Mining_Composite_%d.tif", yr))
  if (!file.exists(cp)) { warning("Missing composite: ", yr); next }
  
  comp <- rast(cp)
  names(comp) <- if (nlyr(comp) >= length(BAND_NAMES))
    BAND_NAMES[seq_len(nlyr(comp))]
  else paste0("Band_", seq_len(nlyr(comp)))
  
  # Single extract call — much faster than one per band
  used_bands <- intersect(c(QC_CONFIG$quality_bands,
                            QC_CONFIG$reflectance_bands,
                            QC_CONFIG$index_bands), names(comp))
  if (length(used_bands) == 0) { rm(comp); next }
  ext_vals <- terra::extract(comp[[used_bands]], sample_cells$cell)
  
  # ---- Quality bands ----------------------------------------
  for (b in intersect(QC_CONFIG$quality_bands, used_bands)) {
    quality_rows[[paste0(yr, "_", b)]] <-
      summarise_band_by_zone(ext_vals[[b]], sample_zone_keys, yr, b)
  }
  
  # ---- Reflectance bands  (saturation/extreme tracking) -----
  for (b in intersect(QC_CONFIG$reflectance_bands, used_bands)) {
    for (zk in unique(sample_zone_keys)) {
      if (is.na(zk)) next
      v <- ext_vals[[b]][sample_zone_keys == zk]
      v <- v[!is.na(v)]
      if (length(v) < 10) next
      reflectance_rows[[paste0(yr, "_", b, "_", zk)]] <- data.frame(
        year = yr, zone = zk, zone_name = ZONE_LABELS[zk], band = b,
        mean_refl     = mean(v),
        sd_refl       = sd(v),
        pct_saturated = 100 * mean(v >= 0.95 | v <= -0.01),
        pct_negative  = 100 * mean(v < 0),
        p01           = unname(quantile(v, 0.01)),
        p99           = unname(quantile(v, 0.99)),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # ---- Index bands  (range / IQR coherence) -----------------
  for (b in intersect(QC_CONFIG$index_bands, used_bands)) {
    for (zk in unique(sample_zone_keys)) {
      if (is.na(zk)) next
      v <- ext_vals[[b]][sample_zone_keys == zk]
      v <- v[!is.na(v)]
      if (length(v) < 10) next
      index_rows[[paste0(yr, "_", b, "_", zk)]] <- data.frame(
        year = yr, zone = zk, zone_name = ZONE_LABELS[zk], band = b,
        mean_idx         = mean(v),
        sd_idx           = sd(v),
        pct_out_of_range = 100 * mean(abs(v) > 1),
        iqr              = IQR(v),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # ---- NA coverage (computed on a NIR (key band), reported per zone)
  na_band <- if ("NIR" %in% names(comp)) "NIR" else names(comp)[1]
  na_global <- as.numeric(global(is.na(comp[[na_band]]),
                                 fun = "mean", na.rm = FALSE)[1, 1]) * 100
  for (zk in unique(sample_zone_keys)) {
    if (is.na(zk)) next
    v <- ext_vals[[na_band]][sample_zone_keys == zk]
    na_rows[[paste0(yr, "_", zk)]] <- data.frame(
      year     = yr, zone = zk, zone_name = ZONE_LABELS[zk],
      band     = "NA_COVERAGE",
      pct_na_zone   = 100 * mean(is.na(v)),
      pct_na_global = na_global,
      stringsAsFactors = FALSE
    )
  }
  
  rm(comp, ext_vals); gc(verbose = FALSE)
}

zone_year_summary <- do.call(rbind, quality_rows)
reflectance_df    <- do.call(rbind, reflectance_rows)
indices_df        <- do.call(rbind, index_rows)
na_df             <- do.call(rbind, na_rows)

write_csv(zone_year_summary, file.path(out_dir, "tables",
                                       "qc_zone_year_summary.csv"))
write_csv(reflectance_df,    file.path(out_dir, "tables",
                                       "qc_reflectance_quality.csv"))
write_csv(indices_df,        file.path(out_dir, "tables",
                                       "qc_index_quality.csv"))

# ================================================================
# 3. SCORECARD
# ================================================================

cat("\n[3/4] Building per-zone-year scorecard\n")

# Helper that converts a value + threshold into "good"/"moderate"/"poor"
classify <- function(value, thr, lower_is_better = TRUE) {
  if (is.na(value)) return(NA_character_)
  if (lower_is_better) {
    if (value <= thr$good)     return("good")
    if (value <= thr$moderate) return("moderate")
    return("poor")
  } else {
    if (value >= thr$good)     return("good")
    if (value >= thr$moderate) return("moderate")
    return("poor")
  }
}

scorecard_rows <- list()
for (yr in YEARS) {
  for (zk in ZONE_ORDER) {
    
    # 1. Valid pixel count (higher is better)
    vpc <- zone_year_summary |>
      filter(year == yr, zone == zk, band == "VALID_PIXEL_COUNT") |>
      pull(mean)
    vpc_v     <- if (length(vpc) == 1) vpc else NA_real_
    vpc_score <- classify(vpc_v,
                          QC_CONFIG$thresholds$valid_pixel_count,
                          lower_is_better = FALSE)
    
    # 2. Cloud probability (lower better)
    cld <- zone_year_summary |>
      filter(year == yr, zone == zk, band == "CLOUD_PROB") |>
      pull(mean)
    cld_v     <- if (length(cld) == 1) cld else NA_real_
    cld_score <- classify(cld_v,
                          QC_CONFIG$thresholds$cloud_prob_pct,
                          lower_is_better = TRUE)
    
    # 3. Temporal noise: mean of three std bands (lower better)
    nse <- zone_year_summary |>
      filter(year == yr, zone == zk,
             band %in% c("NIR_STD_TEMPORAL", "RED_STD_TEMPORAL",
                         "SWIR1_STD_TEMPORAL")) |>
      pull(mean)
    nse_v     <- if (length(nse) > 0) mean(nse, na.rm = TRUE) else NA_real_
    nse_score <- classify(nse_v,
                          QC_CONFIG$thresholds$temporal_noise,
                          lower_is_better = TRUE)
    
    # 4. NA coverage in this zone (lower better)
    na_pct <- na_df |>
      filter(year == yr, zone == zk) |>
      pull(pct_na_zone)
    na_v     <- if (length(na_pct) == 1) na_pct else NA_real_
    na_score <- classify(na_v, QC_CONFIG$thresholds$na_pct,
                         lower_is_better = TRUE)
    
    # Overall: 2+ poor → SUSPECT, 1 poor or 2+ moderate → CAUTION, else OK
    scs <- c(vpc_score, cld_score, nse_score, na_score)
    scs <- scs[!is.na(scs)]
    n_poor <- sum(scs == "poor")
    n_mod  <- sum(scs == "moderate")
    overall <- if (n_poor >= 2) "SUSPECT"
    else if (n_poor >= 1 || n_mod >= 2) "CAUTION"
    else "OK"
    
    scorecard_rows[[paste0(yr, "_", zk)]] <- data.frame(
      year         = yr, zone = zk, zone_name = ZONE_LABELS[zk],
      vpc          = round(vpc_v, 1),
      cloud        = round(cld_v, 1),
      noise        = round(nse_v, 4),
      na_pct       = round(na_v, 2),
      vpc_score    = vpc_score,
      cloud_score  = cld_score,
      noise_score  = nse_score,
      na_score     = na_score,
      overall      = overall,
      stringsAsFactors = FALSE
    )
  }
}
scorecard <- do.call(rbind, scorecard_rows)
write_csv(scorecard, file.path(out_dir, "tables", "qc_scorecard.csv"))

cat("  Overall distribution:\n")
print(table(scorecard$overall))

suspect <- scorecard |> filter(overall == "SUSPECT")
caution <- scorecard |> filter(overall == "CAUTION")

if (nrow(suspect) > 0) {
  cat("\n  SUSPECT zone-years (consider exclusion or re-compositing):\n")
  print(suspect |> select(year, zone_name, vpc, cloud, noise, na_pct),
        row.names = FALSE)
} else {
  cat("\n  No zone-years flagged as SUSPECT\n")
}
if (nrow(caution) > 0) {
  cat("\n  CAUTION zone-years (usable but noisier):\n")
  print(caution |> select(year, zone_name, vpc, cloud, noise, na_pct),
        row.names = FALSE)
}

# ================================================================
# 4. SCORECARD HEATMAP
# ================================================================

cat("\n[4/4] Plotting scorecard heatmap\n")

heat_df <- scorecard |>
  mutate(
    zone_name = factor(zone_name,
                       levels = unname(ZONE_LABELS[ZONE_ORDER])),
    score_val = case_when(
      overall == "OK"      ~ 1,
      overall == "CAUTION" ~ 2,
      overall == "SUSPECT" ~ 3,
      TRUE                 ~ NA_real_
    )
  )

p <- ggplot(heat_df,
            aes(x = year, y = zone_name, fill = factor(score_val))) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c("1" = "#228833", "2" = "#EE7733", "3" = "#CC3311"),
    labels = c("OK", "Caution", "Suspect"),
    name   = "Quality"
  ) +
  scale_x_continuous(breaks = seq(min(YEARS), max(YEARS), by = 2)) +
  labs(title    = "Composite quality scorecard",
       subtitle = "Indicators: valid-pixel count, cloud probability, temporal noise, NA fraction",
       x = "Year", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid       = element_blank(),
        axis.text.x      = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "qc_scorecard_heatmap.png"),
       p, width = 11, height = 3.4, dpi = 600, bg = "white")

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("\n================================================================\n")
cat(sprintf("QC complete in %.1f min\n", as.numeric(elapsed)))
cat("================================================================\n")
cat("Outputs:\n")
cat("  ", file.path(out_dir, "tables", "qc_zone_year_summary.csv"), "\n")
cat("  ", file.path(out_dir, "tables", "qc_reflectance_quality.csv"), "\n")
cat("  ", file.path(out_dir, "tables", "qc_index_quality.csv"), "\n")
cat("  ", file.path(out_dir, "tables", "qc_scorecard.csv"), "\n")
cat("  ", file.path(plot_dir, "qc_scorecard_heatmap.png"), "\n")