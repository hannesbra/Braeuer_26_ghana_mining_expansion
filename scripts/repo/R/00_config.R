# ================================================================
# 00_config.R
# ----------------------------------------------------------------
# Paths, band definitions, class mapping, CV and model settings.
# Every other script sources this file first.
#
# Main classification: one NATIONAL random forest (script 10) on all
# 27 usable candidate bands (see NATIONAL_CONFIG$features; the band
# audit in 03 supports this choice). The latitude
# zones below are kept for two purposes only:
#   * reporting accuracy per zone (40), composite QC (42), maps (50)
#   * the zone-specific benchmark in the appendix (41)
# ================================================================

# ----------------------------------------------------------------
# 1. PATHS
# ----------------------------------------------------------------
# Data live outside the repository. Set GHANA_MINING_DATA in ~/.Renviron
# (or set it in the session) to point at the data folder; the default
# is a folder called data/ in the project root.
DATA_ROOT <- Sys.getenv("GHANA_MINING_DATA", unset = "data")

PATHS <- list(
  # Inputs
  composite_dir     = file.path(DATA_ROOT, "Ghana_Tifs", "B_MERGED"),
  training_data     = file.path(DATA_ROOT, "Training_data",
                                "Mining_Training_Data_2019.csv"),
  ghana_districts   = file.path(DATA_ROOT, "district_shapefiles",
                                "geo2_gh2000_2010.shp"),
  mining_polygons   = file.path(DATA_ROOT, "Mining_polygons",
                                "Ghana_Mining_Polygons.shp"),
  # MinCom cadastre export (used by 45, 45a)
  lease_kml         = Sys.getenv("GHANA_MINING_LEASES", unset =
                                   file.path(DATA_ROOT, "Licences_Mining_Ghana", "Ghana_Mining_Repository.kml")),
  
  # Outputs - one root per pipeline stage
  out_unconstrained = file.path(DATA_ROOT, "Ghana_Tifs", "P1_UNCONSTRAINED"),
  out_constrained   = file.path(DATA_ROOT, "Ghana_Tifs", "P2_CONSTRAINED"),
  out_uncertainty   = file.path(DATA_ROOT, "Ghana_Tifs", "P3_UNCERTAINTY"),
  out_validation    = file.path(DATA_ROOT, "Ghana_Tifs", "P4_VALIDATION")
)

# Yearly rasters of the national model (10), kept apart from the old
# zone-model rasters in P1_UNCONSTRAINED/predictions. Downstream scripts
# (20, 30, 45, 45a, 50, 60, 61) read from here.
PATHS$national_predictions   <- file.path(PATHS$out_unconstrained, "national_predictions")
PATHS$national_probabilities <- file.path(PATHS$national_predictions, "probabilities")

# ----------------------------------------------------------------
# 1b. RESULTS AND FINAL PRODUCT
# ----------------------------------------------------------------
# P1-P4 hold models, rasters and intermediate tables. Everything that
# goes into the paper lives in P5_RESULTS, and the shareable data
# product for the socio-economic analysis in P6_FINAL_FOOTPRINT.
#
#   P5_RESULTS/
#     figures/<NN_topic>/   one folder per producing script (NN = script number)
#     tables/               all .tex tables
#     statistics/           key CSVs (mining+water area inside polygons,
#                           per year, per district x year)
#   P6_FINAL_FOOTPRINT/     constrained footprint per year (rasters,
#                           polygons, onset raster, README) - script 70
PATHS$results   <- file.path(DATA_ROOT, "Ghana_Tifs", "P5_RESULTS")
PATHS$figures   <- file.path(PATHS$results, "figures")
PATHS$tables    <- file.path(PATHS$results, "tables")
PATHS$key_stats <- file.path(PATHS$results, "statistics")
PATHS$final     <- file.path(DATA_ROOT, "Ghana_Tifs", "P6_FINAL_FOOTPRINT")

# Folder helpers (base R only, so scripts that source just this file
# can use them). Each returns the path and creates it if missing.
results_fig_dir <- function(topic) {
  d <- file.path(PATHS$figures, topic)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}
results_tab_dir <- function() {
  dir.create(PATHS$tables, showWarnings = FALSE, recursive = TRUE)
  PATHS$tables
}
results_stat_dir <- function() {
  dir.create(PATHS$key_stats, showWarnings = FALSE, recursive = TRUE)
  PATHS$key_stats
}

# ----------------------------------------------------------------
# 2. CRS
# ----------------------------------------------------------------
GHANA_CRS     <- "EPSG:32630"   # UTM Zone 30N
GHANA_REF_LON <- -1.5           # used for lat -> UTM conversion of zone edges

# ----------------------------------------------------------------
# 3. BAND DEFINITIONS
# ----------------------------------------------------------------
# 40 bands in export order, matches the GEE composite layout.
# NEVER reorder without regenerating composites.

BAND_NAMES <- c(
  # Reflectance (6)
  "BLUE", "GREEN", "RED", "NIR", "SWIR1", "SWIR2",
  # Spectral indices (12)
  "NDVI", "NDWI", "MNDWI", "AWEI_NSH", "AWEI_SH",
  "BSI", "NDBI", "IRON_OXIDE", "CLAY_INDEX", "CI", "BI", "WATER_PROB",
  # Texture (5)
  "ENTROPY", "NIR_STD_DEV", "NIR_contrast", "NIR_corr", "NIR_var",
  # Topography (4)
  "ELEVATION", "SLOPE", "ASPECT", "TRI",
  # Nightlight (4)
  "NIGHTLIGHT", "NIGHTLIGHT_NORM", "IS_LIT", "NIGHTLIGHT_DENSITY",
  # Quality (8)
  "PIXEL_COUNT", "CLOUD_PROB", "COMPOSITE_YEAR", "DATA_MASK",
  "NIR_STD_TEMPORAL", "RED_STD_TEMPORAL", "SWIR1_STD_TEMPORAL",
  "VALID_PIXEL_COUNT",
  # Anomaly (1, GEE-precomputed)
  "NDVI_ANOMALY"
)

# Quality/housekeeping bands: never used as predictors
HOUSEKEEPING_BANDS <- c(
  "PIXEL_COUNT", "CLOUD_PROB", "COMPOSITE_YEAR", "DATA_MASK",
  "NIR_STD_TEMPORAL", "RED_STD_TEMPORAL", "SWIR1_STD_TEMPORAL",
  "VALID_PIXEL_COUNT"
)

# Bands whose meaning changes within 2000-2019 (DMSP -> VIIRS in 2013).
# A model trained on 2019 values would read the 2000-2012 values on a
# different scale, which biases the back-cast.
TEMPORALLY_UNSTABLE_BANDS <- c("NIGHTLIGHT", "NIGHTLIGHT_NORM", "IS_LIT",
                               "NIGHTLIGHT_DENSITY")
EXCLUDE_UNSTABLE_BANDS <- TRUE

# Candidate predictors for the feature audit (03) and the national
# "all bands" benchmark (41)
CANDIDATE_BANDS <- setdiff(
  BAND_NAMES,
  c(HOUSEKEEPING_BANDS,
    if (EXCLUDE_UNSTABLE_BANDS) TEMPORALLY_UNSTABLE_BANDS)
)

BAND_LIMITS_GLOBAL <- list(
  BLUE = c(0, 0.25),       GREEN = c(0, 0.25),      RED   = c(0, 0.35),
  NIR  = c(0, 0.75),       SWIR1 = c(0, 0.75),      SWIR2 = c(0, 0.70),
  NDVI = c(-1, 1),         NDWI  = c(-1, 1),        MNDWI = c(-1, 1),
  BSI  = c(-1, 1),         NDBI  = c(-1, 1),
  AWEI_NSH = c(-10, 5),    AWEI_SH = c(-10, 5),
  CI   = c(-1, 1),         BI     = c(0, 1),        WATER_PROB = c(0, 3),
  IRON_OXIDE = c(0, 3),    CLAY_INDEX = c(-1, 1),
  ENTROPY    = c(0, 0.40), NIR_STD_DEV = c(0, 0.25),
  ELEVATION  = c(-32, 850),SLOPE = c(0, 75),
  ASPECT     = c(0, 360),  TRI   = c(0, 85),
  NDVI_ANOMALY = c(-1.5, 1.5)
)

# ----------------------------------------------------------------
# 4. ZONE-SPECIFIC FEATURE SETS  (appendix benchmark only, 41 and 03)
# ----------------------------------------------------------------
# Literature-derived sets used by the former zone-specific approach.
REGIONAL_FEATURES <- list(
  northern_savannah = c(
    "BSI", "NDVI", "SWIR1", "NIR", "MNDWI", "NDBI",
    "ENTROPY", "ELEVATION", "SLOPE", "NDVI_ANOMALY"
  ),
  forest_transition = c(
    "NDVI", "BSI", "MNDWI", "SWIR2", "NIR", "RED",
    "NIR_STD_DEV", "TRI", "AWEI_NSH", "AWEI_SH",
    "NDVI_ANOMALY", "ENTROPY"
  ),
  coastal = c(
    "MNDWI", "BSI", "AWEI_NSH", "AWEI_SH", "GREEN", "BLUE",
    "WATER_PROB", "SWIR1", "ENTROPY", "NDVI",
    "NDVI_ANOMALY", "ASPECT"
  )
)

# ----------------------------------------------------------------
# 5. ZONE DEFINITIONS  (reporting strata + appendix benchmark)
# ----------------------------------------------------------------
NORTHERN_LAT <- 8.0
COASTAL_LAT  <- 6.0

ZONE_ORDER <- c("northern_savannah", "forest_transition", "coastal")

ZONE_LABELS <- c(
  northern_savannah = "Northern Savannah",
  forest_transition = "Forest-Transition",
  coastal           = "Coastal"
)

ZONE_DISPLAY <- c(
  northern_savannah = "Northern Transition Zone",
  forest_transition = "Central Gold Belt (Forest)",
  coastal           = "SW Mining Belt (Coastal)"
)

# ----------------------------------------------------------------
# 6. CLASS SYSTEM
# ----------------------------------------------------------------
CLASS_NAMES  <- c("Water", "Mining", "Vegetation")
CLASS_IDS    <- c(Water = 1L, Mining = 2L, Vegetation = 3L)
CLASS_LEVELS <- as.character(CLASS_IDS)          # "1" "2" "3"
CLASS_COLS   <- c(Water = "#0072B2", Mining = "#D55E00", Vegetation = "#009E73")

RAW_TO_3CLASS <- c(
  "0" = 3L, "1" = 1L, "2" = 3L, "3" = 3L,
  "4" = 2L, "5" = 2L, "6" = 3L
)

# ----------------------------------------------------------------
# 7. CROSS-VALIDATION
# ----------------------------------------------------------------
# One set of spatial folds for the whole country, shared by 03, 10, 41
# and 50 (cached in fold_cache), so all OOF predictions are paired.
# Buffer = practical range of a variogram fitted to the mining
# indicator (NDVI if there are too few mining points).

CV_CONFIG <- list(
  n_folds                 = 5,
  
  # Buffer strategy ------------------------------------------------
  buffer_strategy         = "variogram",   # or "fixed"
  buffer_min_m            = 1000,          # floor when variogram range is short
  buffer_max_m            = 25000,         # ceiling
  variogram_cutoff_frac   = 0.33,          # fraction of the study extent
  variogram_min_mining_n  = 10,            # else fall back to NDVI
  buffer_fixed_m          = 5000,          # used when buffer_strategy = "fixed"
  
  # Block size ------------------------------------------------------
  min_block_size          = 5000,          # floor
  block_divisor           = 20,            # study_range / 20 -> block size
  
  # Shared fold cache -----------------------------------------------
  fold_cache              = file.path(PATHS$out_validation, "models",
                                      "cv_folds_national.rds"),
  
  # Ensemble / reproducibility -------------------------------------
  n_models                = 5,
  seed                    = 42
)

# ----------------------------------------------------------------
# 8. NATIONAL MODEL
# ----------------------------------------------------------------
NATIONAL_CONFIG <- list(
  # --- Band set ------------------------------------------------------
  # LOCKED: all 27 candidate bands (CANDIDATE_BANDS minus ENTROPY, which
  # is constant in the 2019 composite). In the feature audit (03) this
  # set beat the decorrelated top-10 under spatial CV (accuracy 0.870
  # vs 0.858, mining F1 0.784 vs 0.758, McNemar p = 0.03). The top-10
  # set stays in selected_file and is the M1 benchmark in 41.
  # Set to NULL to read selected_file instead.
  features       = c("BLUE", "GREEN", "RED", "NIR", "SWIR1", "SWIR2",
                     "NDVI", "NDWI", "MNDWI", "AWEI_NSH", "AWEI_SH",
                     "BSI", "NDBI", "IRON_OXIDE", "CLAY_INDEX", "CI", "BI",
                     "WATER_PROB",
                     "NIR_STD_DEV", "NIR_contrast", "NIR_corr", "NIR_var",
                     "ELEVATION", "SLOPE", "ASPECT", "TRI",
                     "NDVI_ANOMALY"),
  selected_file  = file.path(PATHS$out_validation, "tables",
                             "national_selected_features.csv"),
  top_k          = 10,
  selection_rule = "decorrelated",   # "decorrelated" = max. one band per
  #                                    correlation family; or "top_k"
  cluster_h      = 0.30,             # families = bands with |rho| > 0.7
  
  # --- Composites ----------------------------------------------------
  train_composite    = file.path(PATHS$composite_dir,
                                 "Mining_Composite_2019_full.tif"),
  year_composite_fmt = file.path(PATHS$composite_dir,
                                 "Mining_Composite_%d.tif"),
  
  # --- Random forest -------------------------------------------------
  num_trees       = 500,
  min_node_size   = 3,
  cv_sample_frac  = 0.8,
  ens_sample_frac = 0.6,
  
  # --- Outputs -------------------------------------------------------
  model_path    = file.path(PATHS$out_unconstrained, "models",
                            "national_model.rds"),
  
  # --- Prediction ----------------------------------------------------
  # "missing": only years not yet in PATHS$national_predictions
  # "all"    : regenerate every year (use after changing the model)
  # "none"   : train only
  recompute_mode = "missing"
)

# ----------------------------------------------------------------
# 9. PREDICTION CONFIG
# ----------------------------------------------------------------
PREDICT_CONFIG <- list(
  base_year           = 2019,
  predict_years       = 2000:2019,
  write_probabilities = TRUE   # 3-band probability rasters for 30
)