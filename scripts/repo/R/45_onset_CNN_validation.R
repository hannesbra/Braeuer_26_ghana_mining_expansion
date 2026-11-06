# ================================================================
# 45_CNN_RF_onset_consistency.R
# ----------------------------------------------------------------
# CNN-RF ONSET AGREEMENT CONSISTENCY CHECK for industrial mines,
# with Minerals Commission (MinCom) production data as external
# evidence and a full sensitivity grid.
#
# Design
#   * Industrial mines = the ANCHORS of the socio-economic analysis
#     (Mining_expansion_Ghana/00_config.R). The mining polygons are
#     mixed ASM/industrial clusters, so each industrial footprint is an
#     ANCHOR CLIP:  mining polygons  ∩  buffer(anchor, r)  ∩  Voronoi
#     cell of the anchor (neighbouring anchors such as Tarkwa /
#     Iduapriem / Nsuta or Bibiani / Awaso do not share pixels).
#     r = 5 km (main, = ANCHOR_BUF) with 3 and 10 km as sensitivity
#     (= ANCHOR_ALT).
#   * Positives: patches inside anchor footprints (anchors established
#     <= 2019). Anchors established after 2019 (Namdini) are NEGATIVE
#     CONTROLS: they must not be detected before 2020.
#   * Negatives: artisanal polygons (mining polygons outside every
#     anchor clip + gap), plus context ring, bright bare surfaces and
#     random landscape.
#   * Excluded inputs: topography, nightlights, quality bands,
#     COMPOSITE_YEAR, NDVI_ANOMALY.
#   * Patch classifier (small CNN), NOT a segmentation network.
#   * Out-of-fold inference: anchors are split into spatial groups
#     and folds; every anchor's onset is estimated by a model that
#     never saw that anchor (no memorisation of the mine itself).
#   * Thresholds: CNN tile threshold at a fixed false-positive rate. MAIN
#     = YEAR-SPECIFIC (constant false-alarm rate): each year's threshold is
#     the spec-quantile of a held-out background set scored in that year,
#     which removes sensor/composite drift (Landsat 5/7 vs 8, SLC-off).
#     The 2019-only threshold is kept as "fixed2019" for comparison.
#   * Reference extent (MAX EXTENT) of each mine = RF-P2 2019 cells classed
#     Mining or Water INSIDE the footprint (mining polygons ∩ clip).
#     Vegetation inside the polygons is NOT counted. All area series are
#     restricted to this 2019 mined mask:
#       RF  : Mining+Water cells in year t that are also mined in 2019;
#       CNN : share of the 2019-mined cells lying in detected tiles
#             (tiles weighted by their 2019-mined hectares) x ref extent.
#   * Onset rules (applied identically to CNN and RF area series):
#       T10  (MAIN) first year the area reaches max(10 ha, 10% of the 2019
#            mined extent) and stays there for patience+1 consecutive
#            observed years (establishment threshold, sustained);
#       T50  same with 50% of the 2019 mined extent (expansion timing);
#       walk backward walk from 2019 (>= k tiles detected, stop after
#            `patience` misses);
#       first  60_ rule (first year >= X ha; RF only).
#     share_first = extent in 2000 / extent in 2019 classifies each mine as
#     greenfield (<10%), expanding (10-50%) or established (>=50%).
#   * MinCom: company-level production (published from 2013) gives
#       - a hard upper bound: onset must not be after first production,
#       - operating recall: producing anchor-years must be detected,
#       - an uncensored start where production begins inside the
#         record (e.g. zero months before first output).
#
# Anchor coordinates were verified against the 2019 composite and
# Mindat / company reports and the MinCom cadastre (45a_check_anchors.R).
# Obuasi, Prestea (underground) and Awaso (bauxite) are EXCLUDED from CNN
# training and from the consistency statistics (scored, reported
# descriptively).
#
# Relies on 00_config.R and 01_utils.R (remote-sensing project).
#
# Outputs (PATHS$out_validation/cnn_rf_consistency):
#   tables/mincom_production_raw.csv, mincom_anchor_summary.csv,
#          mincom_unmatched_companies.csv
#   tables/cnn_fold_metrics.csv
#   tables/tile_scores_r{r}.csv, rf_panel_r{r}.csv
#   tables/thresholds_by_year.csv       fixed vs year-specific thresholds
#   tables/background_drift_by_year.csv background FPR at the 2019 threshold
#   tables/cnn_detected_area_by_year.csv
#   tables/relative_onset_frac_sensitivity.csv
#   tables/onset_grid_long.csv          config x anchor onsets
#   tables/agreement_grid.csv           CNN-RF consistency per config
#   tables/mincom_checks_grid.csv       external checks per config
#   tables/sensitivity_marginals.csv
#   tables/MAIN_consistency_table.csv   <- final product (per anchor)
#   tables/MAIN_consistency_summary.csv <- final product (headline)
#   tables/MAIN_consistency_table.tex
#   anchor_footprints_r{r}.gpkg
#   Figures: P5_RESULTS/figures/45_cnn_rf_onset/*.png
# ================================================================

# 00_config.R / 01_utils.R: use what is already in the session, otherwise
# look in <root>/R/ and <root>/ (root = MINING_PROJECT_ROOT or getwd()).
if (!exists("PATHS") || !exists("clip_raster_to_limits")) {
  root <- Sys.getenv("MINING_PROJECT_ROOT", unset = getwd())
  for (f in c("00_config.R", "01_utils.R")) {
    hit <- c(file.path(root, "R", f), file.path(root, f))
    hit <- hit[file.exists(hit)][1]
    if (is.na(hit)) stop("Cannot find ", f, " under ", root,
                         ". setwd() to the project root or source 00/01 first.")
    source(hit)
  }
}

pacman::p_load(terra, sf, dplyr, tidyr, readr, ggplot2, jsonlite, rvest,
               tibble, torch)
if (!torch::torch_is_installed()) {
  torch::install_torch()
  stop("torch backend installed. Restart R (RStudio: Ctrl+Shift+F10) ",
       "and run the script again.", call. = FALSE)
}

sf::sf_use_s2(FALSE)
terraOptions(memfrac = 0.6, tempdir = tempdir())
start_time <- Sys.time()

# ================================================================
# 0. CONFIG
# ================================================================

CFG <- list(
  out_dir            = file.path(PATHS$out_validation, "cnn_rf_consistency"),
  train_composite    = file.path(PATHS$composite_dir, "Mining_Composite_2019_full.tif"),
  year_composite_fmt = file.path(PATHS$composite_dir, "Mining_Composite_%d.tif"),
  years              = sort(PREDICT_CONFIG$predict_years, decreasing = TRUE),
  reference_year     = PREDICT_CONFIG$base_year,
  
  # --- Industrial anchors ----------------------------------------------
  # "embedded" = verified table below (default). "socio" = parse ANCHORS
  # from the socio-economic config (only once it has the same coordinates).
  anchor_source = "embedded",
  socio_config  = Sys.getenv("SOCIO_CONFIG", unset = NA),
  # Excluded from CNN training and from the agreement statistics (still
  # scored and reported descriptively):
  #   Obuasi  - mainly underground; small surface footprint; main lease
  #             missing from the cadastre export
  #   Prestea - mainly underground (Prestea shafts); surface area dominated
  #             by galamsey; anchor 3.5 km outside the current lease
  #   Awaso   - only bauxite mine (lateritic benches, no tailings ponds);
  #             not detected in the 5 km main spec even after the lease snap
  exclude_anchors = c("Obuasi", "Prestea", "Awaso"),
  # Anchors moved onto the largest mining polygon INSIDE their own lease
  # (within snap_max_m). Awaso sat ~300 m outside lease MC7/SF8.
  snap_to_lease = c("Awaso"),
  snap_max_m    = 5000,
  anchor_radius_main = 5000,                  # = ANCHOR_BUF
  # 0 = LEASE CLIP: mining polygons inside the company's official Mining
  # Lease(s) (MinCom cadastre export) and within lease_max_m of the anchor.
  anchor_radii       = c(3000, 5000, 10000, 0),  # ANCHOR_ALT + main + lease
  min_footprint_ha   = 5,     # below: footprint = disc around anchor point
  point_disc_m       = 1500,
  artisanal_gap_m    = 2000,  # artisanal negatives start 2 km outside clips
  
  # --- Official licences (MinCom Mining Cadastre export, 25.09.2026) ------
  lease_dir   = dirname(PATHS$lease_kml),
  lease_kml   = PATHS$lease_kml,
  lease_max_m = 10000,       # lease clip is also limited to 10 km of the anchor
  lease_fallback_m = 5000,   # no lease found -> ordinary 5 km clip
  # anchor -> regex on "owner | code" of Mining Lease records
  lease_map = c(
    # Obuasi: main lease missing from the export (only Fiankuma satellites)
    Iduapriem = "Iduapriem|ML2/187",
    Tarkwa = "Gold Fields Ghana", Damang = "Abosso", Ahafo = "Newmont Ghana Gold|RL7/2",
    Akyem = "Golden Ridge", Chirano = "Chirano", Bibiani = "Bibiani",
    Bogoso = "Blue Gold Bogoso", Prestea = "Blue Gold Bogoso|Prestea Sankofa",
    Asanko = "Asanko", Wassa = "Star \\(Wassa\\)", Edikan = "Perseus",
    Nzema = "Adamus", Namdini = "Cardinal Namdini", Nsuta = "Ghana Manganese",
    Awaso = "Ghana Bauxite"),
  # artisanal negatives must not lie inside ANY large-scale Mining Lease
  art_exclude_leases = TRUE,
  
  # --- MinCom ----------------------------------------------------------
  mincom_url = "https://www.mincom.gov.gh/industry-statistics-data/",
  mincom_csv = file.path(PATHS$out_validation, "cnn_rf_consistency",
                         "mincom_manual_export.csv"),   # fallback
  # anchor -> regex on MinCom "Company Name" (check unmatched list!)
  mincom_map = c(
    Obuasi = "Obuasi", Iduapriem = "Iduapriem", Tarkwa = "Tarkwa",
    Damang = "Damang|Abosso", Ahafo = "Ahafo", Akyem = "Akyem|Golden Ridge",
    Chirano = "Chirano", Bibiani = "Bibiani", Bogoso = "Bogoso",
    Prestea = "Prestea", Asanko = "Asanko", Wassa = "Wassa",
    Edikan = "Perseus|Edikan", Nzema = "Adamus|Nzema",
    Namdini = "Namdini|Cardinal", Nsuta = "Manganese", Awaso = "Bauxite"),
  mincom_minerals = c("Gold", "Manganese", "Bauxite"),
  mincom_late_tol = 1,        # onset may lag first production by <= 1 yr
  
  # --- CNN inputs --------------------------------------------------------
  channels = c("BLUE", "GREEN", "RED", "NIR", "SWIR1", "SWIR2",
               "NDVI", "MNDWI", "BSI", "NDBI", "IRON_OXIDE",
               "CLAY_INDEX", "ENTROPY", "NIR_STD_DEV"),
  qc_bands = c("VALID_PIXEL_COUNT", "CLOUD_PROB"),
  patch_px = 32,
  
  # --- Samples -----------------------------------------------------------
  min_pos_cover      = 0.30,
  pos_per_km2        = 4, pos_min = 5, pos_max = 150,
  n_negatives        = 4000,
  neg_mix            = c(artisanal = 0.40, context = 0.20,
                         hard = 0.20, random = 0.20),
  min_art_cover      = 0.20,
  neg_context_ring_m = 10000,
  neg_exclusion_m    = 1000,
  hard_bsi_quantile  = 0.75,
  
  # --- Folds (spatial, anchor-grouped) ------------------------------------
  k_folds       = 4,
  group_link_m  = 20000,   # anchors closer than this share a group
  fold_assign_m = 30000,   # negatives within this of an anchor join its fold
  block_m       = 25000,
  
  # --- Optimisation -------------------------------------------------------
  epochs = 40, batch_size = 64, lr = 1e-3, weight_decay = 1e-4,
  early_stop = 8,
  
  # --- Back-casting -------------------------------------------------------
  infer_min_cover   = 0.25,   # fallback tile rule (no 2019 mined cells)
  infer_min_mined   = 0.05,   # tile kept if >= 5% of it is 2019 Mining/Water
  infer_max_tiles   = 60,
  min_valid_frac    = 0.60,
  min_cluster_valid = 0.50,
  qc_excuse         = c("SUSPECT"),
  
  # --- Threshold calibration (section 6b) ---------------------------------
  calib_n          = 1500,   # held-out background patches, scored every year
  calib_hard_share = 0.5,    # share of bright bare-soil (high-BSI) patches
  calib_max_dist_m = 50000,  # background within 50 km of a positive anchor
  calib_min_n      = 200,    # fewer usable patches in a year -> 2019 threshold
  
  # --- Relative, sustained onset (T10 / T50) -------------------------------
  # persistence = patience + 1 consecutive observed years
  onset_rel = list(frac = 0.10, frac_grid = c(0.05, 0.10, 0.25), min_ha = 10),
  
  # --- MAIN specification -------------------------------------------------
  main = list(radius = 5000, spec = 0.99, k = 3L, patience = 2L, rf_ha = 10,
              thr_mode = "per_year", comparison = "CNN T10 vs RF-P2 T10"),
  
  # --- SENSITIVITY GRID ---------------------------------------------------
  grid = list(spec     = c(0.95, 0.99, 0.995),
              k        = c(1L, 2L, 3L, 4L),
              patience = c(1L, 2L, 3L),
              rf_ha    = c(5, 10, 25),
              # "per_year_cons" = max(per-year, 2019) is also available
              thr_mode = c("fixed2019", "per_year")),
  agree_tol = c(1, 2),       # years
  
  seed = CV_CONFIG$seed
)

out_dir <- CFG$out_dir
ensure_output_dirs(out_dir, subdirs = c("models", "tables"))
plot_dir <- results_fig_dir("45_cnn_rf_onset")   # figures -> P5_RESULTS
set.seed(CFG$seed); torch_manual_seed(CFG$seed)
device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
jsonlite::write_json(CFG, file.path(out_dir, "models", "config.json"),
                     auto_unbox = TRUE, pretty = TRUE, null = "null")

cat("================================================================\n")
cat("45 — CNN-RF onset consistency (anchor clips, MinCom, sensitivity)\n")
cat("================================================================\n")

# ================================================================
# 1. HELPERS
# ================================================================

open_composite <- function(path) {
  r <- rast(path)
  if (nlyr(r) != length(BAND_NAMES))
    stop(sprintf("%s has %d layers, expected %d", basename(path),
                 nlyr(r), length(BAND_NAMES)))
  names(r) <- BAND_NAMES
  r[[c(CFG$channels, CFG$qc_bands)]]
}

read_patches <- function(r, xy, size = CFG$patch_px) {
  n <- nrow(xy); nl <- nlyr(r); half <- size %/% 2
  out <- array(NA_real_, dim = c(n, nl, size, size),
               dimnames = list(NULL, names(r), NULL, NULL))
  rows <- terra::rowFromY(r, xy[, 2]); cols <- terra::colFromX(r, xy[, 1])
  terra::readStart(r); on.exit(terra::readStop(r))
  for (i in seq_len(n)) {
    if (is.na(rows[i]) || is.na(cols[i])) next
    r0 <- rows[i] - half; c0 <- cols[i] - half
    rr <- max(1, r0):min(nrow(r), r0 + size - 1)
    cc <- max(1, c0):min(ncol(r), c0 + size - 1)
    v  <- terra::readValues(r, row = rr[1], nrows = length(rr),
                            col = cc[1], ncols = length(cc), mat = TRUE)
    for (k in seq_len(nl))
      out[i, k, rr - r0 + 1, cc - c0 + 1] <-
      matrix(v[, k], nrow = length(rr), ncol = length(cc), byrow = TRUE)
  }
  out
}

clamp_channels <- function(a) {
  for (b in dimnames(a)[[2]]) {
    lim <- BAND_LIMITS_GLOBAL[[b]]
    if (!is.null(lim)) a[, b, , ] <- pmin(pmax(a[, b, , ], lim[1]), lim[2])
  }
  a
}

patch_qc <- function(a) {
  n <- dim(a)[1]
  m <- function(b) {
    if (!b %in% dimnames(a)[[2]]) return(rep(NA_real_, n))
    v <- rowMeans(matrix(a[, b, , ], nrow = n), na.rm = TRUE)
    v[is.nan(v)] <- NA_real_; v
  }
  data.frame(valid_pixel_count = m("VALID_PIXEL_COUNT"),
             cloud_prob = m("CLOUD_PROB"))
}

# raw clamped float tensor (NaN = no data) + QC
read_patch_set <- function(r, xy, chunk = 500) {
  xs <- list(); qcs <- list()
  for (s in seq(1, nrow(xy), by = chunk)) {
    e   <- min(nrow(xy), s + chunk - 1)
    arr <- read_patches(r, xy[s:e, , drop = FALSE])
    qcs[[length(qcs) + 1]] <- patch_qc(arr)
    arr <- clamp_channels(arr[, CFG$channels, , , drop = FALSE])
    xs[[length(xs) + 1]] <- torch_tensor(arr, dtype = torch_float())
    rm(arr); gc(verbose = FALSE)
  }
  list(X = torch_cat(xs, dim = 1), qc = bind_rows(qcs))
}

compute_norm <- function(X, rows) {
  mu <- sdv <- numeric(length(CFG$channels))
  for (k in seq_along(CFG$channels)) {
    v <- as_array(X[rows, k, , ])
    mu[k] <- mean(v, na.rm = TRUE); sdv[k] <- max(sd(v, na.rm = TRUE), 1e-6)
  }
  list(mean = mu, sd = sdv)
}

normalise_tensor <- function(X, norm) {
  mu  <- torch_tensor(norm$mean, dtype = torch_float())$view(c(1, -1, 1, 1))
  sdv <- torch_tensor(norm$sd,   dtype = torch_float())$view(c(1, -1, 1, 1))
  valid <- X$isnan()$logical_not()$all(dim = 2, keepdim = TRUE)$to(dtype = torch_float())
  Z <- ((X - mu) / sdv)$nan_to_num(nan = 0) * valid
  list(X = torch_cat(list(Z, valid), dim = 2),
       valid_frac = as.numeric(as_array(valid$mean(dim = c(2, 3, 4)))))
}

cover_frac <- function(squares, target_union) {
  sq <- st_sf(sid = seq_along(squares), geometry = squares)
  ix <- suppressWarnings(st_intersection(sq, st_sf(geometry = target_union)))
  out <- numeric(length(squares))
  if (nrow(ix) > 0) {
    a <- tapply(as.numeric(st_area(ix)), ix$sid, sum)
    out[as.integer(names(a))] <- a
  }
  out / as.numeric(st_area(squares))
}

safe_inter <- function(a, b) {
  x <- suppressWarnings(st_intersection(a, b))
  if (!length(x)) return(st_sfc(st_polygon(), crs = st_crs(a)))
  st_union(x)
}

auc_mw <- function(y, p) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  r <- rank(p); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Backward walk. `st` in DESCENDING year order.
walk_onset <- function(years_desc, st, patience) {
  last_det <- NA_integer_; misses <- 0; stop_year <- NA_integer_
  pending <- integer(0); isolated <- integer(0)
  for (i in seq_along(years_desc)) {
    s <- st[i]; yr <- years_desc[i]
    if (is.na(s) || s == "no_data") next
    if (s == "detected") {
      if (!is.na(last_det) && length(pending)) isolated <- c(isolated, pending)
      last_det <- yr; pending <- integer(0); misses <- 0
    } else {
      misses <- misses + 1; pending <- c(pending, yr)
      if (misses >= patience) { stop_year <- yr; break }
    }
  }
  censored <- is.na(stop_year) && !is.na(last_det)
  onset_min <- NA_integer_
  if (!is.na(last_det) && !censored) {
    mb <- years_desc[st %in% "missed" & years_desc < last_det]
    onset_min <- as.integer(max(mb) + 1)
  }
  list(onset = last_det, onset_min = onset_min, left_censored = censored,
       detected_ref = !is.na(st[1]) && st[1] == "detected",
       n_isolated_miss = length(isolated),
       isolated_miss_years = paste(isolated, collapse = ";"))
}

walk_table <- function(df, patience) {
  df |> arrange(anchor_id, desc(year)) |>
    group_by(anchor_id) |>
    group_modify(function(d, key) {
      w <- walk_onset(d$year, d$status, patience)
      as.data.frame(w)
    }) |> ungroup()
}

# Relative, sustained onset. T10 = first year the disturbed area reaches
# max(min_ha, frac x area in ref_year) and stays at/above it for `persist`
# consecutive OBSERVED years (NA years are skipped, not counted as misses).
# Scale-free (a 3,000 ha and a 300 ha mine are treated alike) and robust to
# small pre-existing disturbance (artisanal pits, misclassified bare soil)
# that trips an absolute 10-ha rule in 2000. Also returns
#   t50         same rule with half of the ref extent (expansion timing),
#   share_first extent in the first observed year / ref extent
#               (~0 greenfield, ~1 established before the study period),
#   censored    TRUE if the rule is already met in the first observed year.
# ref_ha: the mine's 2019 MINED extent (Mining+Water cells in the footprint,
# vegetation excluded). If NULL, the series' own ref_year value is used.
onset_relative <- function(year, area, frac = 0.10, min_ha = 10, persist = 3,
                           ref_year = max(year), half = 0.5, ref_ha = NULL) {
  o <- order(year); year <- year[o]; area <- area[o]
  out <- list(onset = NA_integer_, censored = FALSE, t50 = NA_integer_,
              thr_ha = NA_real_, ref_ha = NA_real_, share_first = NA_real_)
  a_ref <- if (!is.null(ref_ha)) ref_ha else area[year == ref_year]
  obs <- which(!is.na(area))
  if (!length(a_ref) || is.na(a_ref) || a_ref < min_ha || !length(obs)) return(out)
  first_sustained <- function(thr) {
    ok <- area >= thr
    for (j in seq_along(obs)) {
      win <- obs[j:min(length(obs), j + persist - 1)]
      if (all(ok[win])) return(as.integer(year[obs[j]]))
    }
    NA_integer_
  }
  thr   <- max(min_ha, frac * a_ref)
  onset <- first_sustained(thr)
  list(onset = onset,
       censored = !is.na(onset) && onset == year[obs[1]],
       t50 = first_sustained(max(min_ha, half * a_ref)),
       thr_ha = thr, ref_ha = a_ref,
       share_first = area[obs[1]] / a_ref)
}

# df must carry `ref_ha` (2019 mined extent) per anchor
rel_table <- function(df, value_col, persist, min_ha, frac = CFG$onset_rel$frac) {
  df |> arrange(anchor_id, year) |> group_by(anchor_id) |>
    group_modify(function(d, key) as.data.frame(
      onset_relative(d$year, d[[value_col]], frac = frac, min_ha = min_ha,
                     persist = persist, ref_year = CFG$reference_year,
                     ref_ha = if ("ref_ha" %in% names(d)) d$ref_ha[1] else NULL))) |>
    ungroup()
}

# First year (forward) meeting the RF area threshold, as in 60_
first_year_rule <- function(df) {
  df |> group_by(anchor_id) |>
    summarise(onset_first = { h <- sort(year[status == "detected"])
    if (length(h)) h[1] else NA_integer_ },
    .groups = "drop")
}

# ================================================================
# 2. ANCHORS
# ================================================================

cat("[1/9] Industrial anchors\n")

# Parse ONLY the ANCHORS (and EFFECT_LEAD) assignments from the
# socio-economic config, so its libraries/paths are not side-loaded.
read_anchors <- function(path) {
  ex   <- parse(path)
  pick <- function(nm) Filter(function(e) is.call(e) &&
                                identical(e[[1]], as.name("<-")) &&
                                identical(e[[2]], as.name(nm)), ex)
  env <- new.env()
  for (e in pick("EFFECT_LEAD")) eval(e, env)
  a <- pick("ANCHORS")
  if (!length(a)) stop("No ANCHORS assignment in ", path)
  eval(a[[1]], env)
  list(anchors = env$ANCHORS,
       effect_lead = if (exists("EFFECT_LEAD", env)) env$EFFECT_LEAD else 2)
}

ANCHORS_VERIFIED <- tribble(   # verified in 45a (Sept 2026)
  ~name, ~established, ~lat, ~lon,
  "Obuasi",1897,6.22030,-1.65310,  "Iduapriem",1992,5.24310,-2.04330,
  "Tarkwa",1998,5.31840,-2.01360,  "Damang",1997,5.51639,-1.84250,
  "Ahafo",2006,7.01248,-2.35604,   "Akyem",2013,6.35287,-1.01944,
  "Chirano",2005,6.31028,-2.38306, "Bibiani",1997,6.47560,-2.30500,
  "Bogoso",1990,5.52700,-2.03900,  "Prestea",1990,5.43270,-2.14280,
  "Asanko",1997,6.37060,-1.96310,  "Wassa",1998,5.49420,-1.71810,
  "Edikan",2011,5.95800,-1.89700,  "Nzema",2011,5.00000,-2.23330,
  "Namdini",2021,10.64780,-0.65920,"Nsuta",1916,5.28310,-1.96670,
  "Awaso",1941,6.23310,-2.26670)
# Notes: Obuasi = Pompora TSF; Iduapriem = Mindat Iduapriem Mine (not
# Tarkwa town); Bogoso = plant/tailings read off the 2019 composite;
# Namdini = centroid of Cardinal Namdini lease PL9/29 (cadastre);
# Wassa established 1998 (first open-pit ore), not 2005 (Golden Star
# restart); Asanko/Nkran mined 1997-2001, restarted 2015/16.

if (CFG$anchor_source == "socio" && file.exists(CFG$socio_config)) {
  ra <- read_anchors(CFG$socio_config)
  anchors <- ra$anchors
  cat("  anchors read from socio-economic config\n")
} else {
  anchors <- ANCHORS_VERIFIED
  cat("  anchors: verified embedded table\n")
}

r_train <- open_composite(CFG$train_composite)
if (terra::is.lonlat(r_train)) stop("Composite must be projected (UTM 30N).")
res_m <- terra::res(r_train)[1]; rcrs <- terra::crs(r_train)
y_min <- min(CFG$years)

anchors <- anchors |>
  select(name, established, lat, lon) |>
  mutate(anchor_id = seq_len(n()),
         role = case_when(name %in% CFG$exclude_anchors         ~ "excluded",
                          established > CFG$reference_year         ~ "negative_control",
                          TRUE                                     ~ "positive"),
         # "window" anchors: established inside the observable period
         in_window = established > y_min & established <= CFG$reference_year)
anc_sf <- st_as_sf(anchors, coords = c("lon", "lat"), crs = 4326,
                   remove = FALSE) |> st_transform(rcrs)
print(anchors |> select(name, established, role, in_window), n = Inf)

# Voronoi partition so neighbouring anchors never share a footprint
make_voronoi <- function(anc_sf) {
  envl <- st_as_sfc(st_bbox(st_buffer(st_union(anc_sf), 300000)))
  v <- st_collection_extract(st_voronoi(st_union(st_geometry(anc_sf)),
                                        envelope = envl), "POLYGON")
  v <- st_sf(geometry = v)
  v$anchor_id <- anc_sf$anchor_id[vapply(st_intersects(v, anc_sf),
                                         function(i) i[1], integer(1))]
  v[order(v$anchor_id), ]
}
vor <- make_voronoi(anc_sf)

polys <- st_read(PATHS$mining_polygons, quiet = TRUE) |>
  st_transform(rcrs) |> st_make_valid()
polys <- polys[!st_is_empty(polys), ]

# ---- Official Mining Leases (MinCom cadastre KML) --------------------
# The KML keeps attributes in <ExtendedData>, which GDAL's KML driver does
# not always expose, so it is parsed directly with xml2.
read_licences <- function(path) {
  doc <- xml2::read_xml(path); xml2::xml_ns_strip(doc)
  pm  <- xml2::xml_find_all(doc, ".//Placemark")
  rows <- lapply(pm, function(p) {
    v <- function(nm) {
      x <- xml2::xml_find_first(p, sprintf(".//Data[@name='%s']/value", nm))
      if (inherits(x, "xml_missing")) NA_character_ else xml2::xml_text(x)
    }
    cs <- xml2::xml_text(xml2::xml_find_all(p, ".//outerBoundaryIs/LinearRing/coordinates"))
    pl <- lapply(cs, function(z) {
      xy <- strsplit(strsplit(trimws(z), "\\s+")[[1]], ",")
      m  <- do.call(rbind, lapply(xy, function(q) as.numeric(q[1:2])))
      if (is.null(m) || nrow(m) < 4) return(NULL)
      if (any(m[1, ] != m[nrow(m), ])) m <- rbind(m, m[1, ])
      list(m)
    })
    pl <- Filter(Negate(is.null), pl)
    if (!length(pl)) return(NULL)
    st_sf(code = v("code"), owner = v("owner"), type = v("type"),
          status = v("status"), start_date = as.Date(v("start_date")),
          geometry = st_sfc(st_multipolygon(pl), crs = 4326))
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

leases <- NULL
if (file.exists(CFG$lease_kml)) {
  leases <- read_licences(CFG$lease_kml) |> st_transform(rcrs) |> st_make_valid()
  leases <- leases[grepl("Mining Lease", leases$type) &
                     !grepl("Restricted", leases$type), ]
  leases$key <- paste(leases$owner, "|", leases$code)
  lease_union <- st_union(leases)
  cat(sprintf("  %d large-scale Mining Lease polygons read from cadastre\n", nrow(leases)))
  lease_check <- lapply(seq_len(nrow(anc_sf)), function(i) {
    pat <- CFG$lease_map[anc_sf$name[i]]
    own <- if (is.na(pat)) leases[0, ] else leases[grepl(pat, leases$key), ]
    d   <- if (nrow(own)) min(as.numeric(st_distance(anc_sf[i, ], own))) else NA_real_
    tibble(name = anc_sf$name[i], n_leases = nrow(own),
           lease_ha = if (nrow(own)) round(sum(as.numeric(st_area(own))) / 1e4) else 0,
           dist_to_own_lease_m = round(d),
           codes = paste(own$code, collapse = "; "))
  }) |> bind_rows()
  write_csv(lease_check, file.path(out_dir, "tables", "anchor_lease_check.csv"))
  print(lease_check, n = Inf)
} else {
  cat("  lease KML not found — lease clip and lease exclusion disabled\n")
  CFG$anchor_radii <- setdiff(CFG$anchor_radii, 0)
}

# ---- Snap selected anchors onto a mined polygon inside their lease -----
if (!is.null(leases) && length(CFG$snap_to_lease)) {
  for (nm in CFG$snap_to_lease) {
    i   <- which(anc_sf$name == nm); pat <- CFG$lease_map[nm]
    if (!length(i) || is.na(pat)) next
    own  <- leases[grepl(pat, leases$key), ]
    near <- polys[lengths(st_is_within_distance(polys, anc_sf[i, ], CFG$snap_max_m)) > 0, ]
    if (!nrow(own) || !nrow(near)) { cat("  snap skipped for", nm, "\n"); next }
    pcs <- suppressWarnings(st_intersection(st_geometry(near), st_union(st_geometry(own))))
    pcs <- suppressWarnings(st_collection_extract(pcs, "POLYGON"))
    pcs <- pcs[!st_is_empty(pcs)]
    if (!length(pcs)) { cat("  no mined polygon inside the lease for", nm, "\n"); next }
    new <- st_point_on_surface(pcs[which.max(as.numeric(st_area(pcs)))])
    ll  <- st_coordinates(st_transform(new, 4326))
    cat(sprintf("  snapped %s onto its lease: lat %.5f, lon %.5f (moved %.0f m)\n",
                nm, ll[2], ll[1], as.numeric(st_distance(anc_sf[i, ], new))))
    st_geometry(anc_sf)[i] <- new
    anchors$lat[i] <- anc_sf$lat[i] <- ll[2]; anchors$lon[i] <- anc_sf$lon[i] <- ll[1]
  }
  vor <- make_voronoi(anc_sf)
}

# ---- Spatial folds by anchor group ---------------------------------
dmat <- as.numeric(st_distance(anc_sf)); dim(dmat) <- c(nrow(anc_sf), nrow(anc_sf))
grp  <- cutree(hclust(as.dist(dmat), method = "single"), h = CFG$group_link_m)
anchors$group <- anc_sf$group <- grp
gtab <- anchors |> group_by(group) |>
  summarise(n_pos = sum(role == "positive"), n = n(), .groups = "drop") |>
  arrange(desc(n_pos), desc(n))
fold_load <- setNames(rep(0, CFG$k_folds), seq_len(CFG$k_folds))
gtab$fold <- NA_integer_
for (i in seq_len(nrow(gtab))) {
  f <- as.integer(names(which.min(fold_load)))
  gtab$fold[i] <- f; fold_load[f] <- fold_load[f] + max(gtab$n_pos[i], 0.1)
}
anchors$fold <- anc_sf$fold <- gtab$fold[match(grp, gtab$group)]
cat("  anchor groups -> folds:\n")
print(anchors |> group_by(fold) |>
        summarise(anchors = paste(name, collapse = ", "), .groups = "drop"))

# ================================================================
# 3. MINERALS COMMISSION PRODUCTION DATA
# ================================================================

cat("\n[2/9] Minerals Commission production statistics\n")

get_mincom <- function() {
  cache <- file.path(out_dir, "tables", "mincom_production_raw.csv")
  if (file.exists(cache)) return(read_csv(cache, show_col_types = FALSE))
  tb <- tryCatch({
    pg   <- rvest::read_html(CFG$mincom_url)
    tabs <- rvest::html_table(pg)
    hit  <- which(vapply(tabs, function(t) any(grepl("Company", names(t))),
                         logical(1)))
    if (length(hit)) tabs[[hit[1]]] else NULL
  }, error = function(e) { message("  scrape failed: ", e$message); NULL })
  if (is.null(tb) && file.exists(CFG$mincom_csv))
    tb <- read_csv(CFG$mincom_csv, show_col_types = FALSE)
  if (is.null(tb)) return(NULL)
  write_csv(tb, cache)
  tb
}

mc_raw <- get_mincom()
mincom_ok <- !is.null(mc_raw) && nrow(mc_raw) > 0
if (!mincom_ok) {
  warning("MinCom table unavailable. Save the table from ", CFG$mincom_url,
          " as ", CFG$mincom_csv, " and re-run. MinCom checks are skipped.")
  mc_anchor <- tibble(anchor_id = anchors$anchor_id)
} else {
  names(mc_raw) <- trimws(names(mc_raw))
  mc <- mc_raw |>
    transmute(size    = trimws(`Mine Size`),
              company = trimws(`Company Name`),
              mineral = trimws(`Type of Mineral`),
              year    = as.integer(Year),
              month   = trimws(Month),
              prod    = readr::parse_number(as.character(Production))) |>
    filter(size == "Large-Scale", toupper(company) != "TOTAL PRODUCTION",
           grepl(paste(CFG$mincom_minerals, collapse = "|"), mineral)) |>
    mutate(month_i = match(tolower(month), tolower(month.name)))
  
  first_company_year <- min(mc$year[!is.na(mc$month_i)], na.rm = TRUE)
  cat(sprintf("  %d Large-Scale company rows; company-level monthly data from %d\n",
              nrow(mc), first_company_year))
  
  # yearly production per company (monthly sum; ANNUAL row if no months)
  mc_year <- mc |>
    group_by(company, year) |>
    summarise(prod = if (any(!is.na(month_i))) sum(prod[!is.na(month_i)], na.rm = TRUE)
              else sum(prod, na.rm = TRUE),
              .groups = "drop")
  
  mc_anchor <- lapply(seq_len(nrow(anchors)), function(i) {
    pat <- CFG$mincom_map[anchors$name[i]]
    if (is.na(pat)) return(NULL)
    comp <- unique(mc$company[grepl(pat, mc$company, ignore.case = TRUE)])
    if (!length(comp)) return(tibble(anchor_id = anchors$anchor_id[i],
                                     mc_company = NA_character_,
                                     mc_in_record = FALSE))
    m  <- mc |> filter(company %in% comp, !is.na(month_i)) |>
      group_by(year, month_i) |> summarise(prod = sum(prod, na.rm = TRUE),
                                           .groups = "drop") |>
      arrange(year, month_i)
    y  <- mc_year |> filter(company %in% comp) |>
      group_by(year) |> summarise(prod = sum(prod), .groups = "drop")
    pos <- which(m$prod > 0)
    first_pos_year <- if (length(pos)) m$year[pos[1]] else NA_integer_
    # censored: producing already in the very first reported month
    censored <- length(pos) > 0 && pos[1] == 1 && m$year[1] == first_company_year
    tibble(anchor_id        = anchors$anchor_id[i],
           mc_company       = paste(comp, collapse = " | "),
           mc_in_record     = TRUE,
           mc_first_report  = min(m$year),
           mc_first_prod    = first_pos_year,
           mc_start_censored= censored,
           mc_prod_years    = paste(y$year[y$prod > 0], collapse = ";"),
           mc_zero_years    = paste(y$year[y$prod <= 0], collapse = ";"))
  }) |> bind_rows()
  
  shared <- mc_anchor |> filter(!is.na(mc_company)) |>
    count(mc_company) |> filter(n > 1)
  mc_anchor$mc_shared_company <- mc_anchor$mc_company %in% shared$mc_company
  
  unmatched <- setdiff(unique(mc$company),
                       unique(unlist(strsplit(na.omit(mc_anchor$mc_company), " \\| "))))
  write_csv(tibble(company = unmatched),
            file.path(out_dir, "tables", "mincom_unmatched_companies.csv"))
  cat("  Large-Scale companies NOT mapped to an anchor (check!):\n   ",
      paste(unmatched, collapse = "\n    "), "\n")
}
mc_anchor <- left_join(anchors |> select(anchor_id, name, established, role),
                       mc_anchor, by = "anchor_id")
write_csv(mc_anchor, file.path(out_dir, "tables", "mincom_anchor_summary.csv"))

# producing anchor-years inside the study period (for operating recall)
mc_prod_long <- mc_anchor |>
  filter(!is.na(mc_prod_years), mc_prod_years != "") |>
  select(anchor_id, mc_prod_years) |>
  separate_rows(mc_prod_years, sep = ";") |>
  transmute(anchor_id, year = as.integer(mc_prod_years), producing = TRUE) |>
  filter(year %in% CFG$years)

# ================================================================
# 4. MODEL
# ================================================================

conv_block <- function(i, o) nn_sequential(
  nn_conv2d(i, o, 3, padding = 1), nn_batch_norm2d(o), nn_relu(),
  nn_conv2d(o, o, 3, padding = 1), nn_batch_norm2d(o), nn_relu())

mine_cnn <- nn_module(
  "mine_cnn",
  initialize = function(in_ch) {
    self$b1 <- conv_block(in_ch, 32); self$b2 <- conv_block(32, 64)
    self$b3 <- conv_block(64, 128)
    self$pool <- nn_max_pool2d(2); self$gap <- nn_adaptive_avg_pool2d(1)
    self$drop <- nn_dropout(0.3); self$fc <- nn_linear(128, 1)
  },
  forward = function(x) {
    x <- self$pool(self$b1(x)); x <- self$pool(self$b2(x)); x <- self$b3(x)
    x <- torch_flatten(self$gap(x), start_dim = 2)
    self$fc(self$drop(x))$squeeze(2)
  })

augment_batch <- function(x) {
  C <- x$size(2) - 1; H <- x$size(3); n <- x$size(1)
  for (j in seq_len(sample(0:3, 1))) x <- torch_flip(x$transpose(3, 4), dims = 4)
  if (runif(1) < 0.5) x <- torch_flip(x, dims = 4)
  spec <- x[, 1:C, , , drop = FALSE]
  mask <- x[, C + 1, , , drop = FALSE]
  spec <- spec * (1 + 0.10 * torch_randn(n, C, 1, 1, device = x$device)) +
    0.20 * torch_randn(n, C, 1, 1, device = x$device)
  if (runif(1) < 0.3) {                        # SLC-off-like stripes
    m <- matrix(1, H, H); m[sample(H, round(0.15 * H)), ] <- 0
    mask <- mask * torch_tensor(m, dtype = torch_float(),
                                device = x$device)$view(c(1, 1, H, H))
  }
  torch_cat(list(spec * mask, mask), dim = 2)
}

predict_probs <- function(model, X, bs = 256) {
  model$eval(); n <- X$size(1); out <- numeric(n)
  with_no_grad({
    for (s in seq(1, n, by = bs)) {
      e <- min(n, s + bs - 1)
      out[s:e] <- as.numeric(as_array(torch_sigmoid(
        model(X[s:e, , , , drop = FALSE]$to(device = device)))$cpu()))
    }
  })
  out
}

train_cnn <- function(X, y, idx_tr, idx_va, seed) {
  torch_manual_seed(seed); set.seed(seed)
  model <- mine_cnn(in_ch = length(CFG$channels) + 1)$to(device = device)
  opt   <- optim_adamw(model$parameters, lr = CFG$lr, weight_decay = CFG$weight_decay)
  pw    <- sum(y[idx_tr] == 0) / max(1, sum(y[idx_tr] == 1))
  lossf <- nn_bce_with_logits_loss(pos_weight = torch_tensor(pw, device = device))
  best <- -Inf; best_state <- NULL; since <- 0
  for (ep in seq_len(CFG$epochs)) {
    model$train()
    b_list <- split(sample(idx_tr), ceiling(seq_along(idx_tr) / CFG$batch_size))
    for (b in b_list[lengths(b_list) >= 2]) {
      xb <- augment_batch(X[b, , , , drop = FALSE]$to(device = device))
      yb <- torch_tensor(y[b], dtype = torch_float(), device = device)
      opt$zero_grad(); loss <- lossf(model(xb), yb); loss$backward(); opt$step()
    }
    auc <- auc_mw(y[idx_va], predict_probs(model, X[idx_va, , , , drop = FALSE]))
    if (!is.na(auc) && auc > best) {
      best <- auc; since <- 0
      best_state <- lapply(model$state_dict(), function(t) t$clone())
    } else if ((since <- since + 1) >= CFG$early_stop) break
  }
  if (!is.null(best_state)) model$load_state_dict(best_state)
  list(model = model, val_auc = best, epochs = ep,
       p_va = predict_probs(model, X[idx_va, , , , drop = FALSE]),
       y_va = y[idx_va])
}

# ================================================================
# 5. PER RADIUS: footprints, samples, out-of-fold CNNs, tile scores, RF
# ================================================================

qc_path <- file.path(PATHS$out_validation, "tables", "qc_scorecard.csv")
qc_card <- if (file.exists(qc_path))
  read_csv(qc_path, show_col_types = FALSE) |> select(year, zone, qc = overall) else
    tibble(year = integer(0), zone = character(0), qc = character(0))
anchors$zone <- assign_zone_from_latitude(anchors$lat)

path_p2 <- function(yr) file.path(PATHS$out_constrained, "predictions",
                                  sprintf("prediction_%d_shrinking.tif", yr))
path_p1 <- function(yr) file.path(PATHS$national_predictions,
                                  sprintf("prediction_%d.tif", yr))

build_footprints <- function(r) {
  rows <- lapply(seq_len(nrow(anc_sf)), function(i) {
    pt   <- st_geometry(anc_sf)[i]
    cell <- st_geometry(vor)[vor$anchor_id == anc_sf$anchor_id[i]]
    csrc <- "buffer"
    if (r == 0) {                                  # lease clip
      pat  <- CFG$lease_map[anc_sf$name[i]]
      own  <- if (is.na(pat)) NULL else leases[grepl(pat, leases$key), ]
      clip <- if (!is.null(own) && nrow(own))
        safe_inter(safe_inter(st_buffer(pt, CFG$lease_max_m), cell),
                   st_union(st_geometry(own)))
      else st_sfc(st_polygon(), crs = st_crs(cell))
      if (!length(clip) || as.numeric(st_area(clip)) < 1) {
        clip <- safe_inter(st_buffer(pt, CFG$lease_fallback_m), cell)
        csrc <- "buffer_fallback"
      } else csrc <- "lease"
    } else clip <- safe_inter(st_buffer(pt, r), cell)
    hit  <- lengths(st_intersects(polys, clip)) > 0
    fp   <- if (any(hit)) safe_inter(clip, st_union(st_geometry(polys)[hit]))
    else st_sfc(st_polygon(), crs = st_crs(clip))
    a    <- as.numeric(st_area(fp)) / 1e4
    src  <- "polygons"
    if (!length(a) || is.na(a) || a < CFG$min_footprint_ha) {
      fp  <- safe_inter(st_buffer(pt, CFG$point_disc_m), cell)
      src <- "point_disc"; a <- as.numeric(st_area(fp)) / 1e4
    }
    list(clip = clip, fp = fp, area = a, src = src, csrc = csrc)
  })
  fp <- st_sf(anchor_id = anc_sf$anchor_id, name = anc_sf$name,
              role = anc_sf$role, fold = anc_sf$fold,
              fp_area_ha = vapply(rows, `[[`, numeric(1), "area"),
              fp_source = vapply(rows, `[[`, character(1), "src"),
              clip_source = vapply(rows, `[[`, character(1), "csrc"),
              geometry = do.call(c, lapply(rows, `[[`, "fp")))
  clips <- do.call(c, lapply(rows, `[[`, "clip"))
  list(fp = fp, clips_u = st_union(clips))
}

build_samples <- function(FP) {
  fp <- FP$fp; half_m <- CFG$patch_px * res_m / 2
  posfp <- fp[fp$role == "positive" & fp$fp_source == "polygons", ]
  pos <- lapply(seq_len(nrow(posfp)), function(i) {
    n <- max(CFG$pos_min, min(CFG$pos_max, round(posfp$fp_area_ha[i] / 100 * CFG$pos_per_km2)))
    p <- st_sample(st_geometry(posfp)[i], size = n, type = "random")
    if (!length(p)) return(NULL)
    cv <- cover_frac(st_buffer(p, half_m, endCapStyle = "SQUARE"), st_geometry(posfp)[i])
    st_sf(label = 1L, type = "industrial", anchor_id = posfp$anchor_id[i],
          fold = posfp$fold[i], geometry = p)[cv >= CFG$min_pos_cover, ]
  })
  pos <- do.call(rbind, Filter(function(x) !is.null(x) && nrow(x) > 0, pos))
  
  n_neg <- round(CFG$n_negatives * CFG$neg_mix)
  mk <- function(g, type) if (!length(g)) NULL else
    st_sf(label = 0L, type = type, anchor_id = NA_integer_,
          fold = NA_integer_, geometry = g)
  
  cat("    positives:", nrow(pos), "\n")
  
  # artisanal: whole polygons >= gap from every clip (indexed, no union)
  far <- lengths(st_is_within_distance(polys, FP$clips_u, CFG$artisanal_gap_m)) == 0
  if (CFG$art_exclude_leases && !is.null(leases))   # not inside any industrial lease
    far <- far & lengths(st_intersects(polys, lease_union)) == 0
  art <- polys[far, ]
  ids <- table(sample(nrow(art), 3 * n_neg[["artisanal"]], replace = TRUE,
                      prob = as.numeric(st_area(art))))
  ap  <- do.call(c, lapply(names(ids), function(j)
    st_sample(st_geometry(art)[as.integer(j)], ids[[j]], exact = TRUE)))
  ac  <- pmin(cover_frac(st_buffer(ap, half_m, endCapStyle = "SQUARE"),
                         st_geometry(polys)), 1)
  ap  <- ap[ac >= CFG$min_art_cover]
  ap  <- ap[sample(length(ap), min(length(ap), n_neg[["artisanal"]]))]
  cat("    artisanal negatives:", length(ap), "\n")
  
  # context: ring around clips, drop points near any polygon
  ring <- st_difference(st_buffer(FP$clips_u, CFG$neg_context_ring_m), FP$clips_u)
  cp   <- st_sample(ring, 3 * n_neg[["context"]])
  cp   <- cp[lengths(st_is_within_distance(cp, polys, CFG$neg_exclusion_m)) == 0]
  cp   <- cp[seq_len(min(length(cp), n_neg[["context"]]))]
  cat("    context negatives:", length(cp), "\n")
  
  cs <- st_as_sf(CAND, coords = c("x", "y"), crs = rcrs)   # sampled once, see below
  cs <- cs[lengths(st_is_within_distance(cs, polys, CFG$neg_exclusion_m)) == 0, ]
  hp <- which(cs$BSI >= quantile(cs$BSI, CFG$hard_bsi_quantile, na.rm = TRUE))
  hi <- sample(hp, min(length(hp), n_neg[["hard"]]))
  ri <- sample(setdiff(seq_len(nrow(cs)), hi), min(nrow(cs) - length(hi), n_neg[["random"]]))
  cat("    hard / random negatives:", length(hi), "/", length(ri), "\n")
  
  neg <- do.call(rbind, Filter(Negate(is.null), list(
    mk(ap, "artisanal"), mk(cp, "context"),
    mk(st_geometry(cs)[hi], "hard"), mk(st_geometry(cs)[ri], "random"))))
  
  # negatives join the fold of a nearby anchor, otherwise a block-hashed fold
  nn <- st_nearest_feature(neg, anc_sf)
  dd <- as.numeric(st_distance(neg, anc_sf[nn, ], by_element = TRUE))
  xy <- st_coordinates(neg)
  hash_fold <- ((floor(xy[, 1] / CFG$block_m) * 31 +
                   floor(xy[, 2] / CFG$block_m)) %% CFG$k_folds) + 1
  neg$fold <- ifelse(dd <= CFG$fold_assign_m, anc_sf$fold[nn], hash_fold)
  
  s <- rbind(pos, neg)
  xy <- st_coordinates(s); s$x <- xy[, 1]; s$y <- xy[, 2]
  s
}

# RF area inside each footprint, restricted to the 2019 MINED MASK
# (RF-P2 2019 = Mining or Water, cell centre inside the footprint).
# Processed per anchor on small crops (fast, low memory). Returns
#   panel : anchor x year x series rf_ha (Mining+Water in year t AND mined 2019)
#   ref   : ref_ha = 2019 mined extent, fp_ha = polygon footprint area
#   masks : per-anchor 2019 mined mask (1/NA) for weighting CNN tiles
rf_panel_for <- function(fp) {
  ref_path <- if (file.exists(path_p2(CFG$reference_year))) path_p2(CFG$reference_year)
  else path_p1(CFG$reference_year)
  if (!file.exists(ref_path)) return(NULL)
  ref <- rast(ref_path)
  fpt <- st_transform(fp, crs(ref))
  cell_ha <- prod(res(ref)) / 1e4
  is_mined <- function(x) x == CLASS_IDS[["Mining"]] | x == CLASS_IDS[["Water"]]
  masks <- list(); refs <- list(); out <- list()
  for (i in seq_len(nrow(fpt))) {
    id <- fpt$anchor_id[i]; g <- fpt[i, ]
    if (st_is_empty(g) || as.numeric(st_area(g)) == 0) next
    e <- terra::intersect(ext(vect(st_buffer(g, 100))), ext(ref))
    if (is.null(e)) {             # footprint outside the RF prediction extent
      cat(sprintf("    %s: no RF prediction coverage -> skipped\n", fpt$name[i]))
      refs[[length(refs) + 1]] <- data.frame(anchor_id = id, ref_ha = NA_real_,
                                             fp_ha = NA_real_)
      next
    }
    refc   <- crop(ref, e)
    inside <- terra::rasterize(vect(g), refc, field = 1, background = NA)
    m19    <- ifel(is_mined(refc) & !is.na(inside), 1, NA)
    masks[[as.character(id)]] <- m19
    n19 <- as.numeric(global(m19, "sum", na.rm = TRUE)[1, 1])
    nfp <- as.numeric(global(inside, "sum", na.rm = TRUE)[1, 1])
    refs[[length(refs) + 1]] <- data.frame(
      anchor_id = id,
      ref_ha = ifelse(is.finite(n19), n19, 0) * cell_ha,
      fp_ha  = ifelse(is.finite(nfp), nfp, 0) * cell_ha)
    for (series in c("P1", "P2")) for (yr in CFG$years) {
      pth <- if (series == "P1") path_p1(yr) else path_p2(yr)
      if (!file.exists(pth)) next
      rr    <- crop(rast(pth), refc)
      n     <- as.numeric(global(is_mined(rr) * m19, "sum", na.rm = TRUE)[1, 1])
      valid <- as.numeric(global((!is.na(rr)) * m19, "sum", na.rm = TRUE)[1, 1])
      out[[length(out) + 1]] <- data.frame(
        anchor_id = id, year = yr, series = series,
        rf_ha = if (is.finite(valid) && valid > 0) n * cell_ha else NA_real_)
    }
  }
  list(panel = bind_rows(out), ref = bind_rows(refs), masks = masks, cell_ha = cell_ha)
}

# CNN scoring tiles. Tiles are weighted by their 2019 MINED hectares (from
# the mask above), so vegetation inside the polygons carries no weight.
# Kept: tiles with >= infer_min_mined of their area mined in 2019 (the
# infer_max_tiles most-mined tiles if there are more). Anchors without
# 2019 mined cells (negative control) fall back to polygon cover.
build_tiles <- function(fp, masks = NULL) {
  side <- CFG$patch_px * res_m; tile_ha <- side^2 / 1e4
  lapply(seq_len(nrow(fp)), function(i) {
    g  <- st_geometry(fp)[i]; id <- fp$anchor_id[i]
    if (st_is_empty(g)) return(NULL)
    gr <- st_make_grid(g, cellsize = side, square = TRUE)
    gr <- gr[lengths(st_intersects(gr, g)) > 0]
    if (!length(gr)) return(NULL)
    cv <- cover_frac(gr, g)
    m  <- masks[[as.character(id)]]
    mined <- numeric(length(gr))
    if (!is.null(m)) {
      ex  <- terra::extract(m, vect(st_transform(gr, crs(m))))
      cnt <- tapply(!is.na(ex[[2]]), ex$ID, sum)
      mined[as.integer(names(cnt))] <- as.numeric(cnt) * prod(res(m)) / 1e4
    }
    use_mined <- sum(mined) > 0
    key  <- if (use_mined) mined else cv
    keep <- if (use_mined) which(mined / tile_ha >= CFG$infer_min_mined)
    else which(cv >= CFG$infer_min_cover)
    if (!length(keep)) keep <- which.max(key)
    if (length(keep) > CFG$infer_max_tiles)
      keep <- keep[order(-key[keep])][seq_len(CFG$infer_max_tiles)]
    xy <- st_coordinates(st_centroid(gr[keep]))
    data.frame(anchor_id = id, fold = fp$fold[i], x = xy[, 1], y = xy[, 2],
               cover = cv[keep], mined_ha = mined[keep],
               tile_rule = if (use_mined) "mined2019" else "cover")
  }) |> bind_rows() |> mutate(tile_id = row_number())
}

# Background candidates for hard/random negatives do not depend on the
# radius, so they are sampled ONCE. method = "regular" reads the big
# composite row-wise; "random" does thousands of scattered block reads
# and can take very long on a Ghana-wide multiband GeoTIFF.
cat("\n  sampling background candidates (once) ...\n")
n_bg <- sum(round(CFG$n_negatives * CFG$neg_mix[c("hard", "random")]))
CAND <- spatSample(r_train[[c("BSI", "NDVI")]], size = 12 * n_bg,
                   method = "regular", na.rm = TRUE, xy = TRUE)
CAND <- CAND[sample(nrow(CAND)), ]
CAND$cid <- seq_len(nrow(CAND))
cat("  ", nrow(CAND), "candidates\n")

# Held-out CALIBRATION background for year-specific thresholds (section 6b).
# Removed from CAND so it is never used for training in any fold/radius.
cal_sf <- st_as_sf(CAND, coords = c("x", "y"), crs = rcrs, remove = FALSE)
keep <- lengths(st_is_within_distance(cal_sf, polys, CFG$neg_exclusion_m)) == 0 &
  lengths(st_is_within_distance(cal_sf, anc_sf[anc_sf$role == "positive", ],
                                CFG$calib_max_dist_m)) > 0
if (exists("lease_union")) keep <- keep & lengths(st_intersects(cal_sf, lease_union)) == 0
pool     <- st_drop_geometry(cal_sf[keep, ])
bsi_q    <- quantile(pool$BSI, CFG$hard_bsi_quantile, na.rm = TRUE)
hard_ids <- pool$cid[pool$BSI >= bsi_q]
rest_ids <- setdiff(pool$cid, hard_ids)
n_h      <- round(CFG$calib_n * CFG$calib_hard_share)
cal_ids  <- c(sample(hard_ids, min(n_h, length(hard_ids))),
              sample(rest_ids, min(CFG$calib_n - n_h, length(rest_ids))))
CALIB <- CAND[CAND$cid %in% cal_ids, ]
CALIB$type <- ifelse(CALIB$cid %in% hard_ids, "hard", "random")
CAND  <- CAND[!CAND$cid %in% cal_ids, ]
write_csv(CALIB, file.path(out_dir, "tables", "calibration_background_points.csv"))
cat(sprintf("  calibration background: %d patches (%d hard / %d random), held out\n",
            nrow(CALIB), sum(CALIB$type == "hard"), sum(CALIB$type == "random")))
rm(cal_sf, pool)

fold_metrics <- list(); val_neg <- list(); tile_all <- list(); rf_all <- list()
ref_all <- list()
models_all <- list(); norms_all <- list(); radius_of <- list()

for (r in CFG$anchor_radii) {
  tag <- if (r == 0) "lease" else sprintf("r%dkm", r / 1000)
  cat(sprintf("\n[3/9] Radius %d m — footprints & samples\n", r))
  FP <- build_footprints(r)
  st_write(FP$fp, file.path(out_dir, sprintf("anchor_footprints_%s.gpkg", tag)),
           delete_dsn = TRUE, quiet = TRUE)
  print(st_drop_geometry(FP$fp) |> select(name, role, fold, fp_area_ha, fp_source, clip_source))
  
  S <- build_samples(FP)
  write_csv(st_drop_geometry(S), file.path(out_dir, "tables", sprintf("samples_%s.csv", tag)))
  print(table(S$type, S$fold))
  
  ps   <- read_patch_set(r_train, cbind(S$x, S$y))
  Xraw <- ps$X; y <- S$label
  
  cat(sprintf("[4/9] Radius %d m — out-of-fold CNNs\n", r))
  models <- list(); norms <- list()
  pos_folds <- sort(unique(S$fold[S$label == 1]))
  for (f in seq_len(CFG$k_folds)) {
    if (!any(anchors$fold == f)) next
    # validation fold = next fold that contains positives
    cand   <- setdiff(pos_folds, f)
    vf     <- if (any(cand > f)) min(cand[cand > f]) else min(cand)
    idx_tr <- which(!S$fold %in% c(f, vf)); idx_va <- which(S$fold == vf)
    idx_te <- which(S$fold == f)
    norm   <- compute_norm(Xraw, idx_tr)
    nt     <- normalise_tensor(Xraw, norm)
    ok     <- nt$valid_frac >= CFG$min_valid_frac
    fit    <- train_cnn(nt$X, y, idx_tr[ok[idx_tr]], idx_va[ok[idx_va]],
                        seed = CFG$seed + f)
    p_te   <- predict_probs(fit$model, nt$X[idx_te[ok[idx_te]], , , , drop = FALSE])
    y_te   <- y[idx_te[ok[idx_te]]]; ty <- S$type[idx_te[ok[idx_te]]]
    thr99  <- unname(quantile(fit$p_va[fit$y_va == 0], 0.99))
    fold_metrics[[paste(tag, f)]] <- data.frame(
      radius = r, fold = f, val_auc = fit$val_auc, epochs = fit$epochs,
      test_auc = auc_mw(y_te, p_te),
      test_sens_at_spec99 = mean(p_te[y_te == 1] >= thr99),
      test_fpr_artisanal  = mean(p_te[ty == "artisanal"] >= thr99),
      test_fpr_other      = mean(p_te[y_te == 0 & ty != "artisanal"] >= thr99),
      n_test_pos = sum(y_te == 1), n_test_neg = sum(y_te == 0),
      test_anchors = paste(anchors$name[anchors$fold == f], collapse = ", "))
    val_neg[[paste(tag, f)]] <- data.frame(radius = r, fold = f,
                                           p = fit$p_va[fit$y_va == 0])
    torch_save(fit$model, file.path(out_dir, "models", sprintf("cnn_%s_fold%d.pt", tag, f)))
    models[[f]] <- fit$model; norms[[f]] <- norm
    cat(sprintf("  fold %d: val AUC %.3f | test AUC %.3f | FPR artisanal %.3f\n",
                f, fit$val_auc, fold_metrics[[paste(tag, f)]]$test_auc,
                fold_metrics[[paste(tag, f)]]$test_fpr_artisanal))
    rm(nt); gc(verbose = FALSE)
  }
  models_all[[tag]] <- models; norms_all[[tag]] <- norms; radius_of[[tag]] <- r
  rm(Xraw, ps); gc(verbose = FALSE)
  
  cat(sprintf("[5/9] Radius %d m — RF inside the 2019 mined mask of each footprint\n", r))
  rfp <- rf_panel_for(FP$fp)
  if (!is.null(rfp)) {
    rf <- rfp$panel; rf$radius <- r
    write_csv(rf, file.path(out_dir, "tables", sprintf("rf_panel_%s.csv", tag)))
    rf_all[[tag]]  <- rf
    ref_all[[tag]] <- rfp$ref |> mutate(radius = r,
                                        mined_share_2019 = ifelse(fp_ha > 0, ref_ha / fp_ha, NA))
    print(ref_all[[tag]] |> left_join(anchors |> select(anchor_id, name), by = "anchor_id") |>
            transmute(name, fp_ha = round(fp_ha), mined2019_ha = round(ref_ha),
                      mined_share = round(mined_share_2019, 2),
                      T10_threshold_ha = round(pmax(CFG$onset_rel$min_ha,
                                                    CFG$onset_rel$frac * ref_ha))))
  }
  
  cat(sprintf("[6/9] Radius %d m — scoring tiles 2019 -> 2000 (out-of-fold)\n", r))
  tiles <- build_tiles(FP$fp, if (!is.null(rfp)) rfp$masks else NULL)
  write_csv(tiles, file.path(out_dir, "tables", sprintf("tiles_%s.csv", tag)))
  ts_r  <- list()
  for (yr in CFG$years) {
    cp <- sprintf(CFG$year_composite_fmt, yr)
    if (!file.exists(cp)) { warning("Missing composite ", yr); next }
    rs <- read_patch_set(open_composite(cp), cbind(tiles$x, tiles$y))
    prob <- vf <- rep(NA_real_, nrow(tiles))
    for (f in seq_len(CFG$k_folds)) {
      ii <- which(tiles$fold == f); if (!length(ii)) next
      nt <- normalise_tensor(rs$X[ii, , , , drop = FALSE], norms[[f]])
      prob[ii] <- predict_probs(models[[f]], nt$X); vf[ii] <- nt$valid_frac
    }
    ts_r[[as.character(yr)]] <- data.frame(radius = r, year = yr,
                                           tiles[, c("tile_id", "anchor_id", "fold", "cover", "mined_ha")],
                                           prob = prob, valid_frac = vf, rs$qc)
    rm(rs); gc(verbose = FALSE)
  }
  ts_r <- bind_rows(ts_r)
  write_csv(ts_r, file.path(out_dir, "tables", sprintf("tile_scores_%s.csv", tag)))
  tile_all[[tag]] <- ts_r
}

fold_metrics <- bind_rows(fold_metrics); val_neg <- bind_rows(val_neg)
tile_all <- bind_rows(tile_all); rf_all <- bind_rows(rf_all)
ref_tab  <- bind_rows(ref_all)
write_csv(ref_tab, file.path(out_dir, "tables", "mined_extent_2019_by_footprint.csv"))
write_csv(fold_metrics, file.path(out_dir, "tables", "cnn_fold_metrics.csv"))
if (!nrow(rf_all)) stop("No RF predictions found — the consistency check needs P1/P2.")

# ================================================================
# 6b. YEAR-SPECIFIC THRESHOLD CALIBRATION (constant false-alarm rate)
# ================================================================
# The 2019 threshold is the spec-quantile of 2019 validation negatives.
# Older composites (Landsat 5/7, SLC-off stripes, fewer clear scenes) shift
# the score distribution of NON-mine land, so a fixed threshold lets more
# background through in early years and every mine looks "present since
# 2000". Here a held-out background set (never used in training, >= 1 km
# from any mining polygon, outside all industrial leases, <= 50 km from a
# positive anchor) is scored in EVERY year and
#     thr[radius, fold, year, spec] = quantile(p_background[year], spec)
# keeps the background false-positive rate at 1 - spec in each year
# (CFAR logic, Finn & Johnson 1968). Assumption: land not mined in 2019
# was not mined before (reclamation is negligible in the study period).

cat("\n[6b/9] Year-specific threshold calibration on held-out background\n")
cal_rows <- list()
for (yr in CFG$years) {
  cp <- sprintf(CFG$year_composite_fmt, yr)
  if (!file.exists(cp)) next
  rs <- read_patch_set(open_composite(cp), cbind(CALIB$x, CALIB$y))
  for (tag in names(models_all)) for (f in seq_along(models_all[[tag]])) {
    m <- models_all[[tag]][[f]]; if (is.null(m)) next
    nt <- normalise_tensor(rs$X, norms_all[[tag]][[f]])
    cal_rows[[paste(tag, f, yr)]] <- data.frame(
      radius = radius_of[[tag]], fold = f, year = yr, cid = CALIB$cid,
      type = CALIB$type, p = predict_probs(m, nt$X), valid_frac = nt$valid_frac)
    rm(nt)
  }
  rm(rs); gc(verbose = FALSE)
  cat(sprintf("  %d scored\n", yr))
}
cal_scores <- bind_rows(cal_rows) |> filter(valid_frac >= CFG$min_valid_frac)
write_csv(cal_scores, file.path(out_dir, "tables", "calibration_background_scores.csv"))

thr_fixed <- val_neg |> group_by(radius, fold) |>
  reframe(spec = CFG$grid$spec, thr = quantile(p, CFG$grid$spec))
thr_year  <- cal_scores |> group_by(radius, fold, year) |>
  reframe(spec = CFG$grid$spec, thr_y = quantile(p, CFG$grid$spec), n_cal = n()) |>
  left_join(thr_fixed, by = c("radius", "fold", "spec"))
thr_all <- bind_rows(
  tidyr::crossing(thr_fixed, year = CFG$years) |> mutate(thr_mode = "fixed2019"),
  thr_year |> mutate(thr = ifelse(n_cal >= CFG$calib_min_n, thr_y, thr),
                     thr_mode = "per_year") |>
    select(radius, fold, spec, thr, year, thr_mode),
  thr_year |> mutate(thr = ifelse(n_cal >= CFG$calib_min_n, pmax(thr_y, thr), thr),
                     thr_mode = "per_year_cons") |>
    select(radius, fold, spec, thr, year, thr_mode)) |>
  filter(thr_mode %in% CFG$grid$thr_mode)
write_csv(thr_all, file.path(out_dir, "tables", "thresholds_by_year.csv"))

# Drift diagnostic: background FPR per year at the FIXED 2019 threshold.
# Flat at 1 - spec = no drift; higher in early years = sensor/composite drift.
drift <- cal_scores |>
  left_join(thr_fixed |> filter(spec == CFG$main$spec) |> select(radius, fold, thr),
            by = c("radius", "fold")) |>
  group_by(radius, year, type) |>
  summarise(mean_p = mean(p), fpr_fixed2019 = mean(p >= thr), n = n(),
            .groups = "drop")
write_csv(drift, file.path(out_dir, "tables", "background_drift_by_year.csv"))
print(drift |> filter(radius == CFG$main$radius) |>
        select(year, type, fpr_fixed2019) |>
        pivot_wider(names_from = type, values_from = fpr_fixed2019) |>
        arrange(year), n = Inf)

# ================================================================
# 7. SENSITIVITY GRID: onsets per configuration
# ================================================================

cat("\n[7/9] Sensitivity grid\n")
tile_ha  <- (CFG$patch_px * res_m)^2 / 1e4
tile_all <- tile_all |>
  left_join(anchors |> select(anchor_id, zone), by = "anchor_id") |>
  left_join(qc_card, by = c("year", "zone"))

onset_rows <- list(); status_rows <- list(); area_rows <- list(); status_main <- NULL
for (tm in CFG$grid$thr_mode) for (r in CFG$anchor_radii) for (sp in CFG$grid$spec) {
  th <- thr_all |> filter(thr_mode == tm, radius == r, spec == sp) |>
    select(fold, year, thr)
  base <- tile_all |> filter(radius == r) |>
    left_join(th, by = c("fold", "year")) |>
    mutate(usable = !is.na(valid_frac) & valid_frac >= CFG$min_valid_frac,
           det    = usable & !is.na(thr) & prob >= thr) |>
    group_by(anchor_id, year) |>
    summarise(n_tiles = n(), n_usable = sum(usable), n_detected = sum(det),
              mined_usable = sum(mined_ha[usable]), mined_det = sum(mined_ha[det]),
              qc = first(qc), .groups = "drop") |>
    left_join(ref_tab |> filter(radius == r) |> select(anchor_id, ref_ha),
              by = "anchor_id") |>
    mutate(enough = n_usable / n_tiles >= CFG$min_cluster_valid,
           # share of the 2019-mined cells (usable tiles) that lie in detected
           # tiles, expressed in hectares of the 2019 mined extent
           cnn_frac = ifelse(enough & mined_usable > 0, mined_det / mined_usable, NA_real_),
           cnn_ha   = cnn_frac * ref_ha)
  area_rows[[paste(tm, r, sp)]] <- base |>
    transmute(anchor_id, year, thr_mode = tm, radius = r, spec = sp, cnn_frac, cnn_ha, ref_ha)
  
  rel_cnn <- bind_rows(lapply(CFG$grid$patience, function(pt)
    rel_table(base, "cnn_ha", persist = pt + 1, min_ha = CFG$onset_rel$min_ha) |>
      transmute(anchor_id, patience = pt, cnn_rel = onset, cnn_rel_cens = censored,
                cnn_t50 = t50, cnn_share_first = share_first)))
  
  for (k in CFG$grid$k) {
    cy <- base |> mutate(status = case_when(
      !enough                                  ~ "no_data",
      n_detected >= pmin(k, pmax(n_usable, 1)) ~ "detected",
      qc %in% CFG$qc_excuse                    ~ "no_data",
      TRUE                                     ~ "missed"))
    status_rows[[length(status_rows) + 1]] <- cy |>
      transmute(anchor_id, year, thr_mode = tm, radius = r, spec = sp, k = k,
                cnn_det = status == "detected")
    if (tm == CFG$main$thr_mode && r == CFG$main$radius &&
        sp == CFG$main$spec && k == CFG$main$k) status_main <- cy
    for (pt in CFG$grid$patience) {
      onset_rows[[length(onset_rows) + 1]] <- walk_table(cy, pt) |>
        transmute(anchor_id, thr_mode = tm, radius = r, spec = sp, k = k,
                  patience = pt, cnn_onset = onset, cnn_onset_min = onset_min,
                  cnn_censored = left_censored, cnn_detected_ref = detected_ref,
                  cnn_isolated_miss = n_isolated_miss) |>
        left_join(rel_cnn |> filter(patience == pt) |> select(-patience),
                  by = "anchor_id")
    }
  }
}
cnn_grid       <- bind_rows(onset_rows)
cnn_status_all <- bind_rows(status_rows)
cnn_area       <- bind_rows(area_rows)
write_csv(cnn_area, file.path(out_dir, "tables", "cnn_detected_area_by_year.csv"))

rf_rows <- list(); rf_status_all <- list()
for (r in unique(rf_all$radius)) for (s in unique(rf_all$series)) for (h in CFG$grid$rf_ha) {
  d <- rf_all |> filter(radius == r, series == s) |>
    left_join(ref_tab |> filter(radius == r) |> select(anchor_id, ref_ha), by = "anchor_id") |>
    mutate(status = case_when(is.na(rf_ha) ~ "no_data",
                              rf_ha >= h ~ "detected", TRUE ~ "missed"))
  rf_status_all[[length(rf_status_all) + 1]] <- d |>
    transmute(anchor_id, year, radius = r, series = s, rf_ha_thr = h,
              rf_det = status == "detected")
  fy <- first_year_rule(d)
  for (pt in CFG$grid$patience) {
    rl <- rel_table(d, "rf_ha", persist = pt + 1, min_ha = h) |>
      transmute(anchor_id, rf_rel = onset, rf_rel_cens = censored, rf_t50 = t50,
                rf_share_first = share_first)
    rf_rows[[length(rf_rows) + 1]] <- walk_table(d, pt) |>
      transmute(anchor_id, radius = r, series = s, rf_ha = h, patience = pt,
                rf_walk = onset, rf_walk_censored = left_censored) |>
      left_join(fy, by = "anchor_id") |> rename(rf_first = onset_first) |>
      left_join(rl, by = "anchor_id")
  }
}
rf_grid <- bind_rows(rf_rows) |>
  pivot_wider(id_cols = c(anchor_id, radius, rf_ha, patience),
              names_from = series,
              values_from = c(rf_walk, rf_walk_censored, rf_first,
                              rf_rel, rf_rel_cens, rf_t50, rf_share_first),
              names_glue = "{.value}_{tolower(series)}")
rf_status_all <- bind_rows(rf_status_all)

cfg_id <- function(radius, spec, k, patience, rf_ha, thr_mode)
  sprintf("r%g_s%g_k%d_p%d_ha%g_%s", radius / 1000, spec, k, patience, rf_ha, thr_mode)

grid_long <- cnn_grid |>
  inner_join(rf_grid, by = c("anchor_id", "radius", "patience"),
             relationship = "many-to-many") |>
  left_join(anchors |> select(anchor_id, name, established, role, in_window, fold),
            by = "anchor_id") |>
  left_join(mc_anchor |> select(anchor_id, any_of(c("mc_first_prod", "mc_start_censored",
                                                    "mc_shared_company"))),
            by = "anchor_id") |>
  left_join(ref_tab |> select(anchor_id, radius, ref_ha, fp_ha, mined_share_2019),
            by = c("anchor_id", "radius")) |>
  mutate(config = cfg_id(radius, spec, k, patience, rf_ha, thr_mode),
         t10_threshold_ha = pmax(CFG$onset_rel$min_ha, CFG$onset_rel$frac * ref_ha),
         greenfield = rf_share_first_p2 < CFG$onset_rel$frac)
for (cc in c("mc_first_prod", "mc_start_censored", "mc_shared_company"))
  if (!cc %in% names(grid_long)) grid_long[[cc]] <- NA
write_csv(grid_long, file.path(out_dir, "tables", "onset_grid_long.csv"))
cat(sprintf("  %d configurations x %d anchors\n",
            n_distinct(grid_long$config), nrow(anchors)))

# ================================================================
# 7b. CONSISTENCY METRICS
# ================================================================
# Censored onsets carry the first observed year (<= 2000). Two onsets are
# CONSISTENT if |a - b| <= tol, which counts "both <= 2000" as agreement.
# Rates are therefore also reported for WINDOW mines (established
# 2001-2019) and GREENFIELD mines (RF extent in 2000 < 10% of 2019), where
# agreement is not automatic.

consistency <- function(a, b, ca, cb, tol) {
  ok <- !is.na(a) & !is.na(b)
  d  <- a[ok] - b[ok]; unc <- !ca[ok] & !cb[ok]
  tibble(n_pairs = sum(ok),
         n_both_censored = sum(ca[ok] & cb[ok]),
         n_one_missing = sum(xor(is.na(a), is.na(b))),
         consistent_rate = if (sum(ok)) mean(abs(d) <= tol) else NA_real_,
         mae_uncensored  = if (sum(unc)) mean(abs(d[unc])) else NA_real_,
         bias_uncensored = if (sum(unc)) mean(d[unc]) else NA_real_,
         n_uncensored    = sum(unc))
}

# name = c(onset A, censor A, onset B, censor B); censor NA -> "== first year"
COMPARISONS <- list(
  "CNN T10 vs RF-P2 T10"         = c("cnn_rel",   "cnn_rel_cens", "rf_rel_p2",   "rf_rel_cens_p2"),
  "CNN walk vs RF-P2 walk"       = c("cnn_onset", "cnn_censored", "rf_walk_p2",  "rf_walk_censored_p2"),
  "CNN walk vs RF-P2 first (60_)"= c("cnn_onset", "cnn_censored", "rf_first_p2", NA),
  "CNN T10 vs RF-P1 T10"         = c("cnn_rel",   "cnn_rel_cens", "rf_rel_p1",   "rf_rel_cens_p1"))
cens_of <- function(d, onset_col, cens_col)
  if (is.na(cens_col)) d[[onset_col]] %in% y_min else d[[cens_col]] %in% TRUE

agree <- grid_long |> filter(role == "positive") |>
  group_by(config, thr_mode, radius, spec, k, patience, rf_ha) |>
  group_modify(function(d, key) {
    subsets <- list("all positive" = rep(TRUE, nrow(d)),
                    "window 2001-2019" = d$in_window %in% TRUE,
                    "greenfield (RF)" = d$greenfield %in% TRUE)
    bind_rows(lapply(names(COMPARISONS), function(cn) {
      cc <- COMPARISONS[[cn]]
      ca <- cens_of(d, cc[1], cc[2]); cb <- cens_of(d, cc[3], cc[4])
      bind_rows(lapply(names(subsets), function(sn) {
        w <- subsets[[sn]]
        bind_rows(lapply(CFG$agree_tol, function(tl)
          consistency(d[[cc[1]]][w], d[[cc[3]]][w], ca[w], cb[w], tl) |>
            mutate(comparison = cn, subset = sn, tol = tl)))
      }))
    }))
  }) |> ungroup()
write_csv(agree, file.path(out_dir, "tables", "agreement_grid.csv"))

# ---- External (MinCom + anchor) checks per configuration -----------
ext_checks <- function(onset, censored, mc_first, established) {
  has_mc <- !is.na(mc_first)
  tibble(
    n_with_mincom = sum(has_mc),
    late_vs_mincom = if (sum(has_mc)) mean(is.na(onset[has_mc]) |
                                             onset[has_mc] > mc_first[has_mc] + CFG$mincom_late_tol) else NA_real_,
    mae_vs_established = { w <- !is.na(onset) & !censored &
      established > y_min & established <= CFG$reference_year
    if (sum(w)) mean(abs(onset[w] - established[w])) else NA_real_ })
}
recall_tab <- function(det_df) {
  if (!nrow(mc_prod_long)) return(NA_real_)
  j <- mc_prod_long |> left_join(det_df, by = c("anchor_id", "year"))
  mean(j$det %in% TRUE)
}
ctrl_ids <- anchors$anchor_id[anchors$role == "negative_control"]

cfgs <- grid_long |> distinct(config, thr_mode, radius, spec, k, patience, rf_ha)
mc_rows <- vector("list", nrow(cfgs))
for (i in seq_len(nrow(cfgs))) {
  c_ <- cfgs[i, ]
  d  <- grid_long |> filter(config == c_$config, role == "positive")
  cs <- cnn_status_all |> filter(thr_mode == c_$thr_mode, radius == c_$radius,
                                 spec == c_$spec, k == c_$k) |>
    transmute(anchor_id, year, det = cnn_det)
  rs <- rf_status_all |> filter(radius == c_$radius, series == "P2",
                                rf_ha_thr == c_$rf_ha) |>
    transmute(anchor_id, year, det = rf_det)
  mfp <- d$mc_first_prod
  cnn_extra <- function(x) mutate(x, operating_recall = recall_tab(cs),
                                  ctrl_detected_pre2020 = any(cs$det[cs$anchor_id %in% ctrl_ids]))
  rf_extra  <- function(x) mutate(x, operating_recall = recall_tab(rs),
                                  ctrl_detected_pre2020 = any(rs$det[rs$anchor_id %in% ctrl_ids]))
  mc_rows[[i]] <- bind_rows(
    ext_checks(d$cnn_rel, d$cnn_rel_cens %in% TRUE, mfp, d$established) |>
      mutate(method = "CNN T10") |> cnn_extra(),
    ext_checks(d$cnn_onset, d$cnn_censored, mfp, d$established) |>
      mutate(method = "CNN walk") |> cnn_extra(),
    ext_checks(d$rf_rel_p2, d$rf_rel_cens_p2 %in% TRUE, mfp, d$established) |>
      mutate(method = "RF-P2 T10") |> rf_extra(),
    ext_checks(d$rf_walk_p2, d$rf_walk_censored_p2 %in% TRUE, mfp, d$established) |>
      mutate(method = "RF-P2 walk") |> rf_extra()) |>
    mutate(config = c_$config, thr_mode = c_$thr_mode, radius = c_$radius,
           spec = c_$spec, k = c_$k, patience = c_$patience, rf_ha = c_$rf_ha)
}
mc_checks <- bind_rows(mc_rows)
write_csv(mc_checks, file.path(out_dir, "tables", "mincom_checks_grid.csv"))

# ---- Marginal sensitivity: which choice moves the agreement? --------
headline <- agree |> filter(comparison == CFG$main$comparison, tol == 2)
marg <- bind_rows(lapply(c("thr_mode", "radius", "spec", "k", "patience", "rf_ha"), function(p)
  headline |> group_by(subset, value = as.character(.data[[p]])) |>
    summarise(mean_consistent = mean(consistent_rate, na.rm = TRUE),
              min_consistent  = min(consistent_rate, na.rm = TRUE),
              max_consistent  = max(consistent_rate, na.rm = TRUE),
              .groups = "drop") |> mutate(parameter = p)))
write_csv(marg, file.path(out_dir, "tables", "sensitivity_marginals.csv"))

# ================================================================
# 8. MAIN SPECIFICATION — FINAL PRODUCT
# ================================================================

cat("\n[8/9] Main specification consistency table\n")
M <- CFG$main
main_cfg <- cfg_id(M$radius, M$spec, M$k, M$patience, M$rf_ha, M$thr_mode)
fmt_onset <- function(o, c) ifelse(is.na(o), "not detected",
                                   ifelse(c, paste0("<=", o), as.character(o)))

main_tab <- grid_long |> filter(config == main_cfg) |>
  mutate(pattern = case_when(is.na(rf_share_first_p2)   ~ NA_character_,
                             rf_share_first_p2 < CFG$onset_rel$frac ~ "greenfield",
                             rf_share_first_p2 < 0.5    ~ "expanding",
                             TRUE                       ~ "established"),
         cnn_T10  = fmt_onset(cnn_rel,   cnn_rel_cens %in% TRUE),
         rf_T10   = fmt_onset(rf_rel_p2, rf_rel_cens_p2 %in% TRUE),
         rf_T50   = ifelse(is.na(rf_t50_p2), "--", as.character(rf_t50_p2)),
         cnn_walk = fmt_onset(cnn_onset, cnn_censored),
         rf_walk  = fmt_onset(rf_walk_p2, rf_walk_censored_p2 %in% TRUE),
         rf_first = fmt_onset(rf_first_p2, rf_first_p2 %in% y_min),
         diff_T10 = cnn_rel - rf_rel_p2,
         consistent_T10_2y = abs(diff_T10) <= 2,
         mincom_first_prod = ifelse(is.na(mc_first_prod), "not in record",
                                    ifelse(mc_start_censored %in% TRUE,
                                           paste0("<=", mc_first_prod),
                                           as.character(mc_first_prod))),
         cnn_late_vs_mincom = !is.na(mc_first_prod) &
           (is.na(cnn_rel) | cnn_rel > mc_first_prod + CFG$mincom_late_tol),
         rf_late_vs_mincom = !is.na(mc_first_prod) &
           (is.na(rf_rel_p2) | rf_rel_p2 > mc_first_prod + CFG$mincom_late_tol)) |>
  select(name, role, established, in_window, fold, pattern,
         footprint_ha = fp_ha, mined2019_ha = ref_ha, mined_share_2019,
         t10_threshold_ha,
         rf_share_2000 = rf_share_first_p2, cnn_share_2000 = cnn_share_first,
         cnn_T10, rf_T10, diff_T10, consistent_T10_2y, rf_T50,
         cnn_walk, rf_walk, rf_first, mincom_first_prod, mc_shared_company,
         cnn_late_vs_mincom, rf_late_vs_mincom, cnn_isolated_miss,
         cnn_rel, rf_rel_p2, cnn_onset, rf_walk_p2) |>
  arrange(role, established)

rng_txt <- function(x) if (all(is.na(x))) "not detected" else
  paste(range(x, na.rm = TRUE), collapse = "-")
# robustness only over VALID settings (main threshold mode; fixed2019 is
# invalidated by the drift diagnostic)
anchor_robust <- grid_long |> filter(role == "positive", thr_mode == M$thr_mode) |>
  group_by(name) |>
  summarise(share_configs_consistent_T10 = mean(abs(cnn_rel - rf_rel_p2) <= 2, na.rm = TRUE),
            share_configs_consistent_walk = mean(abs(cnn_onset - rf_walk_p2) <= 2, na.rm = TRUE),
            cnn_T10_range = rng_txt(cnn_rel),
            rf_T10_range  = rng_txt(rf_rel_p2),
            .groups = "drop")
main_tab <- left_join(main_tab, anchor_robust, by = "name")
write_csv(main_tab, file.path(out_dir, "tables", "MAIN_consistency_table.csv"))
print(main_tab |> select(name, role, established, pattern, mined2019_ha, t10_threshold_ha,
                         rf_share_2000, cnn_T10, rf_T10,
                         diff_T10, rf_T50, cnn_walk, rf_walk, mincom_first_prod,
                         share_configs_consistent_T10), n = Inf, width = Inf)

rng <- agree |> filter(thr_mode == M$thr_mode) |>
  group_by(comparison, subset, tol) |>
  summarise(grid_min    = min(consistent_rate, na.rm = TRUE),
            grid_median = median(consistent_rate, na.rm = TRUE),
            grid_max    = max(consistent_rate, na.rm = TRUE), .groups = "drop")
main_sum <- agree |> filter(config == main_cfg) |>
  left_join(rng, by = c("comparison", "subset", "tol")) |>
  cross_join(mc_checks |> filter(config == main_cfg) |>
               select(method, operating_recall, late_vs_mincom, ctrl_detected_pre2020) |>
               pivot_wider(names_from = method,
                           values_from = c(operating_recall, late_vs_mincom,
                                           ctrl_detected_pre2020)))
write_csv(main_sum, file.path(out_dir, "tables", "MAIN_consistency_summary.csv"))
print(main_sum |> filter(tol == 2) |>
        select(comparison, subset, n_pairs, n_both_censored, consistent_rate,
               grid_min, grid_median, grid_max, mae_uncensored, bias_uncensored), n = Inf)

# Grid range of the main comparison over valid (per-year threshold) settings
grid_range_main <- agree |>
  filter(thr_mode == M$thr_mode, comparison == M$comparison, tol == 2) |>
  group_by(subset) |>
  summarise(n_configs = n(),
            min    = min(consistent_rate, na.rm = TRUE),
            median = median(consistent_rate, na.rm = TRUE),
            max    = max(consistent_rate, na.rm = TRUE), .groups = "drop")
write_csv(grid_range_main, file.path(out_dir, "tables", "grid_range_main_comparison.csv"))
print(grid_range_main)

# Sensitivity of T10 to the relative threshold (main spec only)
frac_sens <- bind_rows(lapply(CFG$onset_rel$frac_grid, function(fr) {
  cm <- cnn_area |> filter(thr_mode == M$thr_mode, radius == M$radius, spec == M$spec)
  rm_ <- rf_all |> filter(radius == M$radius, series == "P2") |>
    left_join(ref_tab |> filter(radius == M$radius) |> select(anchor_id, ref_ha),
              by = "anchor_id")
  full_join(
    rel_table(cm, "cnn_ha", persist = M$patience + 1, min_ha = CFG$onset_rel$min_ha,
              frac = fr) |> transmute(anchor_id, cnn_T = onset, cnn_cens = censored),
    rel_table(rm_, "rf_ha", persist = M$patience + 1, min_ha = M$rf_ha,
              frac = fr) |> transmute(anchor_id, rf_T = onset, rf_cens = censored),
    by = "anchor_id") |> mutate(frac = fr)
})) |> left_join(anchors |> select(anchor_id, name, established, role), by = "anchor_id")
write_csv(frac_sens, file.path(out_dir, "tables", "relative_onset_frac_sensitivity.csv"))

# LaTeX appendix tables: run 45b_appendix_table.R (reads the CSVs above).


# ================================================================
# 9. PLOTS
# ================================================================

cat("\n[9/9] Plots\n")
lvl <- anchors$name[order(anchors$established)]
rad_lab <- function(x) ifelse(x == "0", "lease clip", paste(as.numeric(x) / 1000, "km clip"))

# (a) per-anchor onsets
pd <- grid_long |> filter(config == main_cfg, role != "excluded") |>
  transmute(name, established,
            `CNN T10` = cnn_rel, `RF-P2 T10` = rf_rel_p2,
            `CNN walk` = cnn_onset, `RF-P2 first >= ha` = rf_first_p2,
            MinCom = mc_first_prod) |>
  pivot_longer(-c(name, established), names_to = "source", values_to = "year") |>
  mutate(name = factor(name, levels = lvl))
p_a <- ggplot(pd, aes(year, name, colour = source, shape = source)) +
  geom_vline(xintercept = y_min, linetype = 3, colour = "grey50") +
  geom_point(data = distinct(pd, name, established), inherit.aes = FALSE,
             aes(x = pmax(established, y_min - 1), y = name),
             shape = 124, size = 4, colour = "black") +
  geom_point(size = 2.6, position = position_dodge(width = 0.6), na.rm = TRUE) +
  scale_colour_manual(values = c(`CNN T10` = "#D55E00", `RF-P2 T10` = "#0072B2",
                                 `CNN walk` = "#E69F00", `RF-P2 first >= ha` = "#56B4E9",
                                 MinCom = "#009E73")) +
  labs(x = NULL, y = NULL, colour = NULL, shape = NULL,
       title = "Onset per industrial mine (main specification)",
       subtitle = "| = established (pre-2000 shown at 1999); points at 2000 = left-censored") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(plot_dir, "MAIN_onsets_per_anchor.png"), p_a,
       width = 9, height = 6.5, dpi = 300)

# (b) sensitivity heatmap for the main comparison
hm <- agree |> filter(comparison == M$comparison, tol == 2, rf_ha == M$rf_ha,
                      thr_mode == M$thr_mode)
p_b <- ggplot(hm, aes(factor(k), factor(patience), fill = consistent_rate)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", consistent_rate)), size = 2.4) +
  facet_grid(rows = vars(subset, spec), cols = vars(radius),
             labeller = labeller(spec = function(x) paste("spec", x), radius = rad_lab)) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(x = "min. detected tiles (k)", y = "patience (T10 persistence = patience + 1)",
       fill = "consistent\n(±2 yr)",
       title = sprintf("%s across the grid (%s threshold, RF >= %g ha)",
                       M$comparison, M$thr_mode, M$rf_ha)) +
  theme_minimal(base_size = 8)
ggsave(file.path(plot_dir, "sensitivity_consistency_heatmap.png"), p_b,
       width = 10, height = 11, dpi = 300)

# (c) marginal sensitivity
p_c <- ggplot(marg, aes(value, mean_consistent, colour = subset, group = subset)) +
  geom_linerange(aes(ymin = min_consistent, ymax = max_consistent),
                 position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  facet_wrap(~ parameter, scales = "free_x", nrow = 1) +
  ylim(0, 1) +
  labs(x = NULL, y = "Consistency rate (±2 yr)", colour = NULL,
       title = sprintf("Which choice moves the agreement? (%s)", M$comparison)) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(plot_dir, "sensitivity_marginals.png"), p_c,
       width = 11, height = 4, dpi = 300)

# (d) CNN detection heatmap (main) with MinCom producing years
if (!is.null(status_main)) {
  hm2 <- status_main |> left_join(anchors |> select(anchor_id, name), by = "anchor_id") |>
    mutate(name = factor(name, levels = rev(lvl)))
  mcp <- mc_prod_long |> left_join(anchors |> select(anchor_id, name), by = "anchor_id") |>
    mutate(name = factor(name, levels = levels(hm2$name)))
  p_d <- ggplot(hm2, aes(year, name, fill = status)) +
    geom_tile(colour = "white") +
    geom_point(data = mcp, aes(year, name), inherit.aes = FALSE,
               shape = 16, size = 1.2, colour = "black") +
    scale_fill_manual(values = c(detected = "#D55E00", missed = "#E0E0E0",
                                 no_data = "#56B4E9")) +
    labs(x = NULL, y = NULL, fill = NULL,
         title = sprintf("CNN detection per anchor-year (main spec, %s threshold)", M$thr_mode),
         subtitle = "• = MinCom reports Large-Scale production that year") +
    theme_minimal(base_size = 10)
  ggsave(file.path(plot_dir, "MAIN_cnn_detection_mincom.png"), p_d,
         width = 9, height = 5.5, dpi = 300)
}

# (e) drift diagnostic: background FPR at the fixed 2019 threshold
p_e <- ggplot(drift |> filter(radius == M$radius),
              aes(year, fpr_fixed2019, colour = type)) +
  geom_hline(yintercept = 1 - M$spec, linetype = 2) +
  geom_vline(xintercept = c(2003, 2012.5), linetype = 3, colour = "grey50") +
  geom_line() + geom_point() +
  scale_colour_manual(values = c(hard = "#D55E00", random = "#0072B2")) +
  labs(x = NULL, y = "Background false-positive rate", colour = "Background",
       title = "Score drift: held-out background at the 2019 threshold",
       subtitle = "dashed = target FPR; dotted = SLC-off (2003) and Landsat 8 (2013)") +
  theme_minimal(base_size = 11)
ggsave(file.path(plot_dir, "background_drift.png"), p_e,
       width = 8, height = 4, dpi = 300)

# (f) area trajectories with T10 (CNN vs RF, main spec)
traj <- bind_rows(
  rf_all |> filter(radius == M$radius, series == "P2") |>
    transmute(anchor_id, year, method = "RF-P2", ha = rf_ha),
  cnn_area |> filter(thr_mode == M$thr_mode, radius == M$radius, spec == M$spec) |>
    transmute(anchor_id, year, method = "CNN", ha = cnn_ha)) |>
  left_join(anchors |> select(anchor_id, name, role), by = "anchor_id") |>
  filter(role != "excluded") |> mutate(name = factor(name, levels = lvl))
vl <- main_tab |> filter(role != "excluded") |>
  transmute(name = factor(name, levels = lvl), CNN = cnn_rel, `RF-P2` = rf_rel_p2) |>
  pivot_longer(-name, names_to = "method", values_to = "T10")
thr_l <- main_tab |> filter(role != "excluded") |>
  transmute(name = factor(name, levels = lvl), t10_threshold_ha)
p_f <- ggplot(traj, aes(year, ha, colour = method)) +
  geom_hline(data = thr_l, aes(yintercept = t10_threshold_ha),
             linetype = 3, colour = "grey40", na.rm = TRUE) +
  geom_line(na.rm = TRUE) +
  geom_vline(data = vl, aes(xintercept = T10, colour = method), linetype = 2, na.rm = TRUE) +
  facet_wrap(~ name, scales = "free_y") +
  scale_colour_manual(values = c(CNN = "#D55E00", `RF-P2` = "#0072B2")) +
  labs(x = NULL, y = "Area within the 2019 mined extent (ha)", colour = NULL,
       title = "Mined-area trajectories, T10 threshold (dotted) and T10 onset (dashed)") +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(plot_dir, "MAIN_area_trajectories_T10.png"), p_f,
       width = 11, height = 8, dpi = 300)

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("\n================================================================\n")
cat(sprintf("Done in %.1f min. Final product: %s\n", as.numeric(elapsed),
            file.path(out_dir, "tables", "MAIN_consistency_table.csv")))
cat("================================================================\n")