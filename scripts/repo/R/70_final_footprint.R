# ================================================================
# 70_final_footprint.R
# ----------------------------------------------------------------
# FINAL DATA PRODUCT for the socio-economic analysis: the
# transition-constrained mining footprint (Mining + Water) inside the
# Maus mining polygons, one layer per year, 2000-2019.
#
# Run after 20 (constrained rasters). Uses 30 and 60 outputs if present.
#
# Outputs (PATHS$final = Ghana_Tifs/P6_FINAL_FOOTPRINT):
#   rasters/footprint_{year}.tif   uint8, polygon bounding box, 30 m UTM 30N
#                                  0 = inside polygon, not disturbed (vegetation)
#                                  1 = water (ponds, flooded pits)
#                                  2 = mining (bare, excavated)
#                                  3 = inside polygon, no observation this year
#                                  NA (255) = outside the mining polygons
#   rasters/onset_year.tif         first year a pixel is Mining or Water
#                                  (2000 = already disturbed in 2000, left-censored;
#                                   NA = never disturbed or outside polygons)
#   rasters/mining_polygon_mask.tif 1 = inside a Maus polygon (cell centre rule)
#   footprint_polygons.gpkg        layer "footprint_{year}": one feature per
#                                  connected patch (class, area_ha, district)
#                                  layer "onset": patches dissolved by onset year;
#                                  footprint in year t = features with onset_year <= t
#   tables/mining_area_by_year.csv          (copy of the key statistic, below)
#   tables/district_mining_panel.csv        (copy from 60, if it was run)
#   README.txt                     coding, CRS, provenance
#
# Key statistic (PATHS$key_stats = P5_RESULTS/statistics):
#   mining_area_in_polygons_by_year.csv
#     year, polygon / mining / water / disturbed / no-data area (ha),
#     unconstrained disturbed area (P1), new disturbance per year (onset),
#     and the 30 uncertainty band (argmax mean, 95% interval) if available.
# ================================================================

source("R/00_config.R")
source("R/01_utils.R")

pacman::p_load(terra, sf, dplyr, readr)

terraOptions(memfrac = 0.6, tempdir = tempdir())
sf::sf_use_s2(FALSE)
start_time <- Sys.time()

# ================================================================
# 0. CONFIG
# ================================================================

FINAL_CONFIG <- list(
  years          = PREDICT_CONFIG$predict_years,
  write_polygons = TRUE,
  min_patch_ha   = 0,      # drop patches smaller than this (0 = keep all;
  #                          0.45 ha = 5 Landsat pixels removes speckle)
  add_district   = TRUE,   # attach district (GEOLEVEL2, ADMIN_NAME) to patches
  gdal_opts      = c("COMPRESS=DEFLATE", "TILED=YES")
)
cfg <- FINAL_CONFIG

out_dir  <- PATHS$final
rast_dir <- file.path(out_dir, "rasters")
tab_out  <- file.path(out_dir, "tables")
for (d in c(rast_dir, tab_out)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
gpkg     <- file.path(out_dir, "footprint_polygons.gpkg")

path_p2 <- function(yr) file.path(PATHS$out_constrained, "predictions",
                                  sprintf("prediction_%d_shrinking.tif", yr))
path_p1 <- function(yr) file.path(PATHS$national_predictions,
                                  sprintf("prediction_%d.tif", yr))

years <- cfg$years[file.exists(path_p2(cfg$years))]
if (!length(years)) stop("No constrained rasters found. Run 20 first.", call. = FALSE)
years <- sort(years)

cat("================================================================\n")
cat("70 - FINAL FOOTPRINT (constrained Mining + Water inside polygons)\n")
cat("================================================================\n")
cat(sprintf("Years: %d-%d (%d)\nOutput: %s\n", min(years), max(years),
            length(years), out_dir))

# ================================================================
# 1. GRID AND POLYGON MASK (polygon bounding box)
# ================================================================

cat("\n[1/5] Polygon mask\n")
ref_full <- rast(path_p2(max(years)))
polys    <- st_read(PATHS$mining_polygons, quiet = TRUE) |>
  st_transform(crs(ref_full)) |> st_make_valid()
ref      <- crop(ref_full, ext(vect(polys)), snap = "out")
inside   <- rasterize_polygons(polys, ref)                # 0/1 on the cropped grid
in_mask  <- ifel(inside == 1L, 1L, NA)
names(in_mask) <- "inside"
cell_ha  <- prod(res(ref)) / 1e4
n_inside <- as.numeric(global(in_mask, "sum", na.rm = TRUE)[1, 1])
cat(sprintf("  %s polygon cells (%.0f ha), cell = %.4f ha\n",
            format(n_inside, big.mark = ","), n_inside * cell_ha, cell_ha))

writeRaster(in_mask, file.path(rast_dir, "mining_polygon_mask.tif"),
            overwrite = TRUE, datatype = "INT1U", NAflag = 255,
            gdal = cfg$gdal_opts)

districts <- NULL
if (cfg$add_district && file.exists(PATHS$ghana_districts))
  districts <- st_read(PATHS$ghana_districts, quiet = TRUE) |>
  st_transform(crs(ref)) |> st_make_valid() |>
  select(any_of(c("GEOLEVEL2", "ADMIN_NAME")))

# ================================================================
# 2. YEARLY FOOTPRINT RASTERS + AREA STATISTICS
# ================================================================

cat("\n[2/5] Yearly footprint rasters\n")

count_vals <- function(r, vals) {
  fr <- terra::freq(r)
  vapply(vals, function(v) {
    x <- fr$count[fr$value == v]
    if (length(x)) as.numeric(sum(x)) else 0
  }, numeric(1))
}

onset <- init(in_mask, NA)                                # empty raster on the grid
stats_rows <- list()

for (yr in years) {
  r  <- crop(rast(path_p2(yr)), ref)
  fp <- ifel(is.na(in_mask), NA,
             ifel(is.na(r), 3L,
                  ifel(r == CLASS_IDS[["Water"]], 1L,
                       ifel(r == CLASS_IDS[["Mining"]], 2L, 0L))))
  names(fp) <- sprintf("footprint_%d", yr)
  writeRaster(fp, file.path(rast_dir, sprintf("footprint_%d.tif", yr)),
              overwrite = TRUE, datatype = "INT1U", NAflag = 255,
              gdal = cfg$gdal_opts)

  # first year disturbed (rule of 20: once disturbed, never vegetation later)
  onset <- ifel(is.na(onset) & (fp == 1L | fp == 2L), yr, onset)

  n <- count_vals(fp, 0:3)
  # unconstrained (P1) for comparison
  u_ha <- NA_real_
  if (file.exists(path_p1(yr))) {
    u  <- crop(rast(path_p1(yr)), ref)
    du <- (u == CLASS_IDS[["Water"]] | u == CLASS_IDS[["Mining"]]) * in_mask
    u_ha <- as.numeric(global(du, "sum", na.rm = TRUE)[1, 1]) * cell_ha
  }
  stats_rows[[as.character(yr)]] <- data.frame(
    year = yr,
    polygon_area_ha  = n_inside * cell_ha,
    water_ha         = n[2] * cell_ha,
    mining_ha        = n[3] * cell_ha,
    disturbed_ha     = (n[2] + n[3]) * cell_ha,
    vegetation_ha    = n[1] * cell_ha,
    nodata_ha        = n[4] * cell_ha,
    valid_share      = (n_inside - n[4]) / n_inside,
    disturbed_unconstrained_ha = u_ha)
  cat(sprintf("  %d: disturbed %8.0f ha (mining %7.0f, water %6.0f), no data %5.1f%%\n",
              yr, (n[2] + n[3]) * cell_ha, n[3] * cell_ha, n[2] * cell_ha,
              100 * n[4] / n_inside))
  rm(r, fp); gc(verbose = FALSE)
}

names(onset) <- "onset_year"
writeRaster(onset, file.path(rast_dir, "onset_year.tif"), overwrite = TRUE,
            datatype = "INT2U", NAflag = 65535, gdal = cfg$gdal_opts)

# ================================================================
# 3. KEY STATISTIC: AREA INSIDE POLYGONS PER YEAR
# ================================================================

cat("\n[3/5] Key statistic\n")
onset_n <- count_vals(onset, years)
area_tbl <- bind_rows(stats_rows) |>
  mutate(new_disturbed_ha        = onset_n[match(year, years)] * cell_ha,
         cumulative_onset_ha     = cumsum(new_disturbed_ha))

# Uncertainty band from 30 (same polygon rule)
stoch_csv <- file.path(PATHS$out_uncertainty, "statistics",
                       "mining_area_in_polygons_stochastic.csv")
if (file.exists(stoch_csv)) {
  st <- read_csv(stoch_csv, show_col_types = FALSE) |>
    select(year, any_of(c("argmax_mean_ha", "argmax_sd_ha",
                          "argmax_pi95_lower_ha", "argmax_pi95_upper_ha")))
  area_tbl <- left_join(area_tbl, st, by = "year")
} else cat("  (30 not run: no uncertainty columns)\n")

area_tbl <- area_tbl |> mutate(across(ends_with("_ha"), ~ round(.x, 1)),
                               valid_share = round(valid_share, 4))
key_csv <- file.path(results_stat_dir(), "mining_area_in_polygons_by_year.csv")
write_csv(area_tbl, key_csv)
file.copy(key_csv, file.path(tab_out, "mining_area_by_year.csv"), overwrite = TRUE)
cat(sprintf("  -> %s\n", key_csv))

# Consistency check against 20 (should match to the cell)
p2_csv <- file.path(PATHS$out_constrained, "statistics", "area_within_mining_polygons.csv")
if (file.exists(p2_csv)) {
  chk <- read_csv(p2_csv, show_col_types = FALSE) |>
    filter(stage == "shrinking", class %in% c(1L, 2L)) |>
    group_by(year) |> summarise(ha20 = sum(area_ha_inside), .groups = "drop") |>
    inner_join(area_tbl |> select(year, disturbed_ha), by = "year")
  cat(sprintf("  max |difference| vs 20's table: %.1f ha\n",
              max(abs(chk$ha20 - chk$disturbed_ha))))
}

# District x year panel from 60 (same rule), copied into the final product
dp <- file.path(PATHS$key_stats, "district_mining_panel.csv")
if (file.exists(dp)) {
  file.copy(dp, file.path(tab_out, "district_mining_panel.csv"), overwrite = TRUE)
  cat("  district x year panel (60) copied\n")
} else cat("  run 60 for the district x year panel\n")

# ================================================================
# 4. POLYGONS (GeoPackage)
# ================================================================

to_patches <- function(r, value_name) {
  if (as.numeric(global(!is.na(r), "sum", na.rm = TRUE)[1, 1]) == 0) return(NULL)
  v <- as.polygons(r, dissolve = TRUE, values = TRUE, na.rm = TRUE)
  v <- disagg(v)                                   # one feature per patch
  s <- st_as_sf(v)
  names(s)[1] <- value_name
  s$area_ha <- round(as.numeric(st_area(s)) / 1e4, 3)
  s <- s[s$area_ha >= cfg$min_patch_ha, ]
  if (!is.null(districts) && nrow(s)) {
    pos <- st_point_on_surface(st_geometry(s))
    j   <- st_join(st_sf(id = seq_along(pos), geometry = pos), districts,
                   join = st_intersects, left = TRUE)
    j   <- j[!duplicated(j$id), ]
    s   <- cbind(s, st_drop_geometry(j)[, setdiff(names(j), c("id", "geometry")),
                                        drop = FALSE])
  }
  s
}

if (cfg$write_polygons) {
  cat("\n[4/5] Polygons -> ", gpkg, "\n", sep = "")
  if (file.exists(gpkg)) file.remove(gpkg)
  for (yr in years) {
    fp <- rast(file.path(rast_dir, sprintf("footprint_%d.tif", yr)))
    d  <- ifel(fp == 1L | fp == 2L, fp, NA)
    s  <- to_patches(d, "class")
    if (is.null(s) || !nrow(s)) next
    s$class_name <- c("Water", "Mining")[s$class]
    s$year       <- yr
    s$patch_id   <- seq_len(nrow(s))
    st_write(s, gpkg, layer = sprintf("footprint_%d", yr),
             delete_layer = TRUE, quiet = TRUE)
    cat(sprintf("  %d: %s patches\n", yr, format(nrow(s), big.mark = ",")))
  }
  s <- to_patches(onset, "onset_year")
  if (!is.null(s) && nrow(s)) {
    s$left_censored <- s$onset_year == min(years)
    st_write(s, gpkg, layer = "onset", delete_layer = TRUE, quiet = TRUE)
    cat(sprintf("  onset: %s patches\n", format(nrow(s), big.mark = ",")))
  }
} else cat("\n[4/5] Polygons skipped (write_polygons = FALSE)\n")

# ================================================================
# 5. README
# ================================================================

cat("\n[5/5] README\n")
nm <- tryCatch(load_national_model(), error = function(e) NULL)
readme <- c(
  "Ghana mining footprint 2000-2019 (transition-constrained)",
  "==========================================================",
  sprintf("Created: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "What it is",
  "  Pixels classed Mining or Water by the national random forest (script 10),",
  "  after the backward transition rule (script 20: a pixel that is vegetation",
  "  in a later year is vegetation in all earlier years), restricted to the",
  "  Maus et al. (2022) mining polygons (cell-centre rule).",
  if (!is.null(nm)) sprintf("  Model: %d bands, %d ensemble members, trained %s.",
                            length(nm$features), length(nm$models),
                            format(nm$created, "%Y-%m-%d")),
  "",
  "Grid",
  sprintf("  CRS %s (UTM 30N), %.0f m cells (%.4f ha), extent = polygon bounding box.",
          GHANA_CRS, res(ref)[1], cell_ha),
  "",
  "rasters/footprint_YYYY.tif (uint8)",
  "  0 = inside polygon, not disturbed   1 = water   2 = mining",
  "  3 = inside polygon, no observation   255 (NA) = outside polygons",
  "rasters/onset_year.tif (uint16)",
  sprintf("  first year classed 1 or 2; %d = disturbed in %d or earlier (left-censored);",
          min(years), min(years)),
  "  NA = never disturbed or outside polygons.",
  "  Footprint in year t (gaps filled) = onset_year <= t.",
  "footprint_polygons.gpkg",
  "  footprint_YYYY: connected patches (class, class_name, area_ha, district, year)",
  "  onset: patches dissolved by onset year (left_censored flag)",
  sprintf("  patches smaller than %.2f ha removed.", cfg$min_patch_ha),
  "tables/mining_area_by_year.csv: area per year (ha); see column names.",
  "tables/district_mining_panel.csv: district x year (if script 60 was run).",
  "",
  "Notes",
  "  Areas from the yearly rasters treat code 3 (no observation) as unknown;",
  "  the onset-based footprint carries the last state across such gaps.",
  "  Uncertainty: argmax_* columns come from the ensemble simulation (script 30)",
  "  and capture member disagreement only (a lower bound)."
)
writeLines(readme, file.path(out_dir, "README.txt"))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("\n================================================================\n")
cat(sprintf("Final footprint complete in %.1f min -> %s\n", as.numeric(elapsed), out_dir))
cat("================================================================\n")
