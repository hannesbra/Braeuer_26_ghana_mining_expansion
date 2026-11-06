# ================================================================
# 30_Uncertainty_Simulations.R
# ----------------------------------------------------------------
# How certain are the mining numbers? Monte Carlo simulation of the
# disturbed area (Mining + Water) inside the Maus mining polygons,
# 2000-2019, using the 5 bootstrapped forests of the national model
# (10) as the source of model uncertainty.
#
# One simulation:
#   1. For every year, draw ONE ensemble member at random (shared by
#      all pixels in that year) and use its class probabilities.
#   2. Draw one uniform number u per pixel and keep it for all years.
#      A pixel is Water if u < P(water), Mining if u < P(water) +
#      P(mining), otherwise Vegetation. Using the same u in every year
#      means a pixel only changes class when its probabilities change,
#      not because of fresh coin flips.
#      (Independent draws per year make the constrained back-cast decay
#      geometrically: 0.9^19 = 14% survival even for a pixel that is 90%
#      mining every year. Set shared_u = FALSE to reproduce that.)
#   3. Apply the transition rule of 20 backwards in time: a pixel that
#      is Vegetation in a later year is Vegetation in all earlier
#      years. Years without data for a pixel are skipped and the rule
#      carries over the gap.
#
# Two versions are run with the same member draws:
#   argmax : each pixel takes the most likely class of the drawn member.
#            Comparable to the deterministic maps of 20 (which use the
#            argmax of the ensemble mean); the spread is pure member
#            disagreement. Main uncertainty band.
#   sample : each pixel draws its class from the member's probabilities
#            (steps 2-3 above). Its mean is the probability-weighted
#            expected area. Because the forests are trained on balanced,
#            class-weighted data, vegetation pixels keep small but
#            non-zero mining probabilities, which add up over ~2 million
#            pixels; this version therefore lies above the argmax maps.
#
# The spread across simulations reflects disagreement between the
# ensemble members (plus, for "sample", pixel-level classification noise). It does
# not include spatially correlated errors shared by all members, so it
# is a lower bound on the total uncertainty.
#
# Only pixels inside the mining polygons are simulated (all area and
# onset statistics are polygon-based). Member probabilities are
# predicted for those pixels only and cached per year.
#
# Outputs (PATHS$out_uncertainty):
#   member_probs/member_probs_{year}.rds   cached member probabilities
#   statistics/mining_area_in_polygons_stochastic.csv
#        year; argmax_* columns (argmax version) and unprefixed columns
#        (sample version): mean/sd/median and 50/90/95% intervals of disturbed area,
#        mean mining and water area, ensemble-mean expected area
#        (no sampling, no constraint), deterministic P2 area (from 20)
#   statistics/simulation_area_matrices.rds  sim x year matrices + draws
#   predictions/disturbed_prob_{year}.tif    P(disturbed) per pixel
#   predictions/member_sd_{year}.tif         SD of P(disturbed) across members
#   predictions/onset_median.tif, onset_iqr.tif, onset_frac_disturbed.tif
#   models/uncertainty_config.{rds,json}
# Rasters cover the bounding box of the mining polygons.
# ================================================================

source("R/00_config.R")
source("R/01_utils.R")

pacman::p_load(terra, sf, dplyr, readr, jsonlite, progress, parallel)

terraOptions(memfrac = 0.6, tempdir = tempdir())
start_time <- Sys.time()

# ================================================================
# 0. CONFIG
# ================================================================

if (!exists("UNCERTAINTY_CONFIG")) {
  UNCERTAINTY_CONFIG <- list(
    reference_year  = PREDICT_CONFIG$base_year,
    n_sim           = 1000,
    shared_u        = TRUE,        # same uniform number per pixel in all years
    apply_shrinking_constraint = TRUE,
    member_per_year = TRUE,        # FALSE: one member for all years of a sim
    pixel_maps_from = "argmax",    # per-pixel / onset rasters: "argmax" or "sample"
    map_years       = c(2000, 2005, 2010, 2015, 2019),  # per-pixel rasters
    predict_chunk   = 400000,      # rows per ranger predict call
    seed            = CV_CONFIG$seed
  )
}
cfg <- UNCERTAINTY_CONFIG
N_THREADS <- max(1L, detectCores() - 1L)

out_dir    <- PATHS$out_uncertainty
member_dir <- file.path(out_dir, "member_probs")
ensure_output_dirs(out_dir, subdirs = c("predictions", "statistics", "models",
                                        "member_probs"))

cat("================================================================\n")
cat(" 30 - UNCERTAINTY: ensemble-draw simulation inside mining polygons\n")
cat("================================================================\n")
cat(sprintf("Simulations: %d | shared u: %s | member per year: %s | constraint: %s\n",
            cfg$n_sim, cfg$shared_u, cfg$member_per_year,
            cfg$apply_shrinking_constraint))

saveRDS(cfg, file.path(out_dir, "models", "uncertainty_config.rds"))
jsonlite::write_json(cfg, file.path(out_dir, "models", "uncertainty_config.json"),
                     auto_unbox = TRUE, pretty = TRUE)

# Probabilities are stored as raw bytes: p * 250, 255 = no data.
enc <- function(p) {
  x <- round(p * 250)
  x <- pmin(pmax(x, 0), 250)
  x[is.na(x)] <- 255
  as.raw(x)
}
dec <- function(r) {
  x <- as.integer(r)
  x[x == 255L] <- NA_integer_
  x / 250
}

# ================================================================
# 1. MODEL, GRID AND POLYGON PIXELS
# ================================================================

cat("\n[1/5] Model, grid and polygon pixels\n")

nm        <- load_national_model()
feats     <- nm$features
models    <- nm$models
n_mem     <- length(models)
model_tag <- paste(format(nm$created), n_mem, paste(feats, collapse = ","))
cat(sprintf("  national model: %d members, %d bands\n", n_mem, length(feats)))

pred_dir <- PATHS$national_predictions
years    <- PREDICT_CONFIG$predict_years
years    <- years[file.exists(sprintf(NATIONAL_CONFIG$year_composite_fmt, years))]
ref_path <- file.path(pred_dir, sprintf("prediction_%d.tif", cfg$reference_year))
if (!file.exists(ref_path)) stop("Run 10 first: ", ref_path, call. = FALSE)
ref_full <- rast(ref_path)
cell_ha  <- prod(res(ref_full)) / 1e4

polys <- sf::st_read(PATHS$mining_polygons, quiet = TRUE) |>
  sf::st_transform(crs(ref_full))
poly_mask  <- rasterize_polygons(polys, ref_full)
cells_full <- terra::cells(poly_mask, 1)[[1]]
n_pix      <- length(cells_full)

ref_crop   <- crop(ref_full, ext(vect(polys)), snap = "out")
cells_crop <- cellFromXY(ref_crop, xyFromCell(ref_full, cells_full))
stopifnot(!anyNA(cells_crop))

cat(sprintf("  polygon pixels: %s (%.0f ha); pixel = %.4f ha\n",
            format(n_pix, big.mark = ","), n_pix * cell_ha, cell_ha))
cat(sprintf("  years: %d-%d (%d)\n", min(years), max(years), length(years)))

# ================================================================
# 2. MEMBER PROBABILITIES AT POLYGON PIXELS (cached per year)
# ================================================================

cat("\n[2/5] Member probabilities (cached per year)\n")

pw_l <- pd_l <- list()       # raw matrices [n_pix x n_mem]: P(water), P(water+mining)
expected_rows <- list()
map_sd <- list()

for (yr in years) {
  cache <- file.path(member_dir, sprintf("member_probs_%d.rds", yr))
  obj   <- if (file.exists(cache)) readRDS(cache) else NULL
  if (is.null(obj) || !identical(obj$tag, model_tag) || obj$n_pix != n_pix) {
    t0   <- Sys.time()
    # Read the raw bands at the polygon cells and clip the values there
    # (same limits as load_composite, without clamping the whole country)
    comp <- rast(sprintf(NATIONAL_CONFIG$year_composite_fmt, yr))
    if (nlyr(comp) != length(BAND_NAMES))
      stop("Composite ", yr, " has ", nlyr(comp), " layers.", call. = FALSE)
    names(comp) <- BAND_NAMES
    if (!compareGeom(comp, ref_full, stopOnError = FALSE))
      stop("Composite ", yr, " is not on the prediction grid.", call. = FALSE)
    X  <- clip_bands_to_limits(terra::extract(comp[[feats]], cells_full))
    ok <- which(stats::complete.cases(X))
    pw <- matrix(as.raw(255), n_pix, n_mem)
    pd <- matrix(as.raw(255), n_pix, n_mem)
    for (s in seq(1, length(ok), by = cfg$predict_chunk)) {
      ch <- ok[s:min(s + cfg$predict_chunk - 1, length(ok))]
      Xc <- X[ch, , drop = FALSE]
      for (m in seq_len(n_mem)) {
        p <- align_probs(predict(models[[m]], Xc, num.threads = N_THREADS,
                                 verbose = FALSE)$predictions)
        pw[ch, m] <- enc(p[, "1"])
        pd[ch, m] <- enc(p[, "1"] + p[, "2"])
      }
    }
    obj <- list(tag = model_tag, year = yr, n_pix = n_pix, pw = pw, pd = pd)
    saveRDS(obj, cache)
    rm(comp, X); gc(verbose = FALSE)
    cat(sprintf("  %d predicted (%d valid pixels, %.1f min)\n", yr, length(ok),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  } else {
    cat(sprintf("  %d loaded from cache\n", yr))
  }
  pw_l[[as.character(yr)]] <- obj$pw
  pd_l[[as.character(yr)]] <- obj$pd
  
  # Ensemble mean / member SD of P(disturbed): no sampling, no constraint
  pdm  <- vapply(seq_len(n_mem), function(m) dec(obj$pd[, m]), numeric(n_pix))
  mu   <- rowMeans(pdm)
  expected_rows[[as.character(yr)]] <- data.frame(
    year = yr, expected_unconstrained_ha = sum(mu, na.rm = TRUE) * cell_ha)
  if (yr %in% cfg$map_years)
    map_sd[[as.character(yr)]] <- sqrt(pmax(rowMeans(pdm^2) - mu^2, 0) * n_mem / (n_mem - 1))
  rm(pdm, mu, obj); gc(verbose = FALSE)
}

# Sanity check: member mean should equal 10's ensemble probability raster
pr_ref <- file.path(PATHS$national_probabilities, sprintf("prob_%d.tif", cfg$reference_year))
if (file.exists(pr_ref)) {
  set.seed(cfg$seed)
  smp <- sample(n_pix, min(5000, n_pix))
  r10 <- terra::extract(rast(pr_ref), cells_full[smp])
  p10 <- r10[[1]] + r10[[2]]
  p30 <- rowMeans(vapply(seq_len(n_mem), function(m)
    dec(pd_l[[as.character(cfg$reference_year)]][smp, m]), numeric(length(smp))))
  cat(sprintf("  check vs 10's ensemble raster (%d px): max |diff| = %.4f (rounding ~0.004)\n",
              length(smp), max(abs(p10 - p30), na.rm = TRUE)))
}

# ================================================================
# 3. SIMULATION
# ================================================================

cat(sprintf("\n[3/5] Simulation (%d runs)\n", cfg$n_sim))

years_desc <- sort(years, decreasing = TRUE)
years_chr  <- sort(years)
n_years    <- length(years_desc)

set.seed(cfg$seed)
member_draw <- if (cfg$member_per_year) {
  matrix(sample.int(n_mem, cfg$n_sim * n_years, replace = TRUE), cfg$n_sim, n_years,
         dimnames = list(NULL, years_desc))
} else {
  matrix(rep(sample.int(n_mem, cfg$n_sim, replace = TRUE), n_years), cfg$n_sim, n_years,
         dimnames = list(NULL, years_desc))
}

sim_dist  <- matrix(0, cfg$n_sim, n_years, dimnames = list(NULL, years_desc))
sim_water <- sim_dist
arg_dist  <- sim_dist                  # argmax version
arg_water <- sim_dist
counts    <- matrix(0L, n_pix, n_years, dimnames = list(NULL, years_desc))

pb <- progress_bar$new(format = "  sim :current/:total [:bar] :percent | ETA :eta",
                       total = cfg$n_sim, clear = FALSE, width = 70)
for (s in seq_len(cfg$n_sim)) {
  u     <- runif(n_pix)
  alive <- rep(TRUE, n_pix)            # disturbed in every later (non-NA) year
  alive_a <- rep(TRUE, n_pix)          # same, argmax version
  for (j in seq_len(n_years)) {
    yr_s <- as.character(years_desc[j])
    m    <- member_draw[s, j]
    if (!cfg$shared_u && j > 1) u <- runif(n_pix)
    pd   <- dec(pd_l[[yr_s]][, m])
    pw   <- dec(pw_l[[yr_s]][, m])
    has  <- !is.na(pd)
    dist <- has & (u < pd)
    if (cfg$apply_shrinking_constraint) dist <- dist & alive
    wat  <- dist & (u < pw)
    alive[has] <- dist[has]
    sim_dist[s, j]  <- sum(dist) * cell_ha
    sim_water[s, j] <- sum(wat)  * cell_ha
    
    # Argmax of the drawn member (ties resolved Water > Mining > Vegetation,
    # as which.max in 10)
    pm     <- pd - pw
    pv     <- 1 - pd
    dist_a <- has & (pmax(pw, pm) >= pv)
    if (cfg$apply_shrinking_constraint) dist_a <- dist_a & alive_a
    wat_a  <- dist_a & (pw >= pm)
    alive_a[has] <- dist_a[has]
    arg_dist[s, j]  <- sum(dist_a) * cell_ha
    arg_water[s, j] <- sum(wat_a)  * cell_ha
    
    counts[, j] <- counts[, j] + if (cfg$pixel_maps_from == "argmax") dist_a else dist
  }
  pb$tick()
}
sim_mining <- sim_dist - sim_water
arg_mining <- arg_dist - arg_water

# ================================================================
# 4. AREA STATISTICS
# ================================================================

cat("\n[4/5] Area statistics\n")

q <- function(x, p) unname(stats::quantile(x, p))
area_tbl <- bind_rows(lapply(years_chr, function(yr) {
  a <- sim_dist[, as.character(yr)]
  data.frame(year = yr,
             mean_ha = mean(a), sd_ha = stats::sd(a), median_ha = stats::median(a),
             pi50_lower_ha = q(a, 0.25),  pi50_upper_ha = q(a, 0.75),
             pi90_lower_ha = q(a, 0.05),  pi90_upper_ha = q(a, 0.95),
             pi95_lower_ha = q(a, 0.025), pi95_upper_ha = q(a, 0.975),
             mean_mining_ha = mean(sim_mining[, as.character(yr)]),
             mean_water_ha  = mean(sim_water[, as.character(yr)]))
})) |>
  left_join(bind_rows(lapply(years_chr, function(yr) {
    a <- arg_dist[, as.character(yr)]
    data.frame(year = yr,
               argmax_mean_ha = mean(a), argmax_sd_ha = stats::sd(a),
               argmax_median_ha = stats::median(a),
               argmax_pi50_lower_ha = q(a, 0.25),  argmax_pi50_upper_ha = q(a, 0.75),
               argmax_pi90_lower_ha = q(a, 0.05),  argmax_pi90_upper_ha = q(a, 0.95),
               argmax_pi95_lower_ha = q(a, 0.025), argmax_pi95_upper_ha = q(a, 0.975),
               argmax_mean_mining_ha = mean(arg_mining[, as.character(yr)]),
               argmax_mean_water_ha  = mean(arg_water[, as.character(yr)]))
  })), by = "year") |>
  left_join(bind_rows(expected_rows), by = "year")

# Deterministic constrained area from 20 (same polygon rule)
p2_csv <- file.path(PATHS$out_constrained, "statistics", "area_within_mining_polygons.csv")
if (file.exists(p2_csv)) {
  det <- read_csv(p2_csv, show_col_types = FALSE) |>
    filter(class %in% c(CLASS_IDS[["Water"]], CLASS_IDS[["Mining"]])) |>
    group_by(year, stage) |>
    summarise(ha = sum(area_ha_inside), .groups = "drop") |>
    tidyr::pivot_wider(names_from = stage, values_from = ha,
                       names_prefix = "deterministic_") |>
    rename_with(~ sub("deterministic_shrinking", "deterministic_constrained_ha", .x)) |>
    rename_with(~ sub("deterministic_unconstrained", "deterministic_unconstrained_ha", .x))
  area_tbl <- left_join(area_tbl, det, by = "year")
} else {
  cat("  (20's polygon table not found; deterministic columns omitted)\n")
}

area_tbl <- area_tbl |> mutate(across(-year, ~ round(.x, 1)))
area_path <- file.path(out_dir, "statistics", "mining_area_in_polygons_stochastic.csv")
write_csv(area_tbl, area_path)
saveRDS(list(sim_disturbed = sim_dist, sim_mining = sim_mining, sim_water = sim_water,
             argmax_disturbed = arg_dist, argmax_mining = arg_mining,
             argmax_water = arg_water,
             member_draw = member_draw, years_desc = years_desc, years_chrono = years_chr,
             cell_ha = cell_ha, n_sim = cfg$n_sim, config = cfg),
        file.path(out_dir, "statistics", "simulation_area_matrices.rds"))

cat(sprintf("  -> %s\n\n", area_path))
print(area_tbl |> select(year, any_of("deterministic_constrained_ha"),
                         argmax_mean_ha, argmax_pi95_lower_ha, argmax_pi95_upper_ha,
                         sample_mean_ha = mean_ha, pi95_lower_ha, pi95_upper_ha),
      row.names = FALSE)

# ================================================================
# 5. PER-PIXEL RASTERS (polygon bounding box)
# ================================================================

cat("\n[5/5] Per-pixel rasters\n")

write_cells <- function(vals, name, path, datatype = "FLT4S") {
  v <- rep(NA_real_, ncell(ref_crop))
  v[cells_crop] <- vals
  r <- setValues(rast(ref_crop), v)
  names(r) <- name
  writeRaster(r, path, overwrite = TRUE, datatype = datatype,
              gdal = c("COMPRESS=DEFLATE"))
  invisible(path)
}

# P(disturbed) per pixel = share of simulations; columns in chronological order
F_chr <- counts[, as.character(years_chr), drop = FALSE] / cfg$n_sim
rm(counts); gc(verbose = FALSE)

for (yr in intersect(cfg$map_years, years_chr)) {
  write_cells(F_chr[, as.character(yr)], "p_disturbed",
              file.path(out_dir, "predictions", sprintf("disturbed_prob_%d.tif", yr)))
  if (!is.null(map_sd[[as.character(yr)]]))
    write_cells(map_sd[[as.character(yr)]], "member_sd",
                file.path(out_dir, "predictions", sprintf("member_sd_%d.tif", yr)))
  cat(sprintf("  %d: disturbed_prob + member_sd\n", yr))
}

# Onset: under the constraint the disturbed set only grows forward in
# time, so P(onset <= t) = P(disturbed in t). Gap years are bridged
# with a running maximum. Quantiles are conditional on the pixel being
# disturbed in at least one simulation (as before).
for (j in seq_len(ncol(F_chr))[-1]) F_chr[, j] <- pmax(F_chr[, j], F_chr[, j - 1])
ever <- F_chr[, ncol(F_chr)]
G    <- F_chr / ifelse(ever > 0, ever, NA_real_)
first_year <- function(p) {
  out <- rep(NA_real_, n_pix)
  for (j in seq_along(years_chr)) {
    hit <- is.na(out) & !is.na(G[, j]) & G[, j] >= p
    out[hit] <- years_chr[j]
  }
  out
}
onset_med <- first_year(0.50)
onset_iqr <- first_year(0.75) - first_year(0.25)
write_cells(onset_med, "onset_median", file.path(out_dir, "predictions", "onset_median.tif"))
write_cells(onset_iqr, "onset_iqr",    file.path(out_dir, "predictions", "onset_iqr.tif"))
write_cells(ever, "onset_frac_disturbed",
            file.path(out_dir, "predictions", "onset_frac_disturbed.tif"))
cat("  onset_median / onset_iqr / onset_frac_disturbed\n")

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("\n================================================================\n")
cat(sprintf("P3 complete in %.1f min\n", as.numeric(elapsed)))
cat("================================================================\n")
cat("Next: 50 (section F draws the deterministic vs probabilistic figure)\n")