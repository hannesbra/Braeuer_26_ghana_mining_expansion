# ================================================================
# Creates district tables for appendix.
# ----------------------------------------------------------------
# Produces the per-district appendix tables (A3-A7):
#   * area measured WITHIN the mining polygons (Maus layer)
#   * "mining area" = Mining + Water class pixels combined
#   * BOTH unconstrained (P1) and constrained (P2 "shrinking") series
#   * annual change (on constrained) + per-district z-score of change
#
# Inputs  : PATHS$out_unconstrained/predictions/prediction_{year}.tif
#           PATHS$out_constrained/predictions/prediction_{year}_shrinking.tif
#           PATHS$mining_polygons, PATHS$ghana_districts
# Outputs : <stat_dir>/district_mining_within_polygons.csv  (all districts, both series)
#           <tab_dir>/table_appendix_<district>.tex          (one per selected district)
# ================================================================

source("R/00_config.R")
source("R/01_utils.R")

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tidyr); library(readr)
})

# ---- Config ---------------------------------------------------
DISTURBANCE_CLASSES <- c(CLASS_IDS[["Mining"]], CLASS_IDS[["Water"]])  # = c(2, 1)
Z_THRESHOLD         <- 1.5     # bold annual-change z-scores above this
years               <- PREDICT_CONFIG$predict_years

# Districts to tabulate: label = regex matched against ADMIN_NAME
DISTRICTS_FOR_TABLES <- c(
  "Tarkwa Nsuaem"  = "Tarkwa Nsuaem",
  "Amansie West"   = "Amansie West",
  "Wassa Amenfi"   = "Wassa Amenfi",
  "Upper Denkyira" = "Upper Denkyira",
  "Asutifi"        = "Asutifi"
)

stat_dir <- results_stat_dir()   # P5_RESULTS/statistics
tab_dir  <- results_tab_dir()    # P5_RESULTS/tables

# ---- Path resolvers for the two series ------------------------
path_unconstr <- function(yr)
  file.path(PATHS$national_predictions, sprintf("prediction_%d.tif", yr))
path_constr <- function(yr)
  file.path(PATHS$out_constrained, "predictions", sprintf("prediction_%d_shrinking.tif", yr))

years <- years[file.exists(vapply(years, path_constr, character(1))) &
                 file.exists(vapply(years, path_unconstr, character(1)))]
stopifnot(length(years) > 0)

# ---- Reference grid, district id raster, polygon inside-mask --
ref     <- rast(path_constr(years[1]))
cell_ha <- as.numeric(global(cellSize(ref, unit = "ha"), "mean", na.rm = TRUE)[1, 1])

districts <- st_read(PATHS$ghana_districts, quiet = TRUE) |>
  st_transform(terra::crs(ref))
districts$didx <- seq_len(nrow(districts))
dist_r <- terra::rasterize(vect(districts), ref, field = "didx")
names(dist_r) <- "didx"

inside      <- rasterize_polygons(PATHS$mining_polygons, ref)   # util from 01_utils.R
inside_bool <- inside == 1L

# ---- Disturbed (Mining+Water) area inside polygons, per district ----
disturbed_ha_by_district <- function(path) {
  r <- rast(path)
  dist_in <- (r %in% DISTURBANCE_CLASSES) & inside_bool
  z <- terra::zonal(dist_in, dist_r, fun = "sum", na.rm = TRUE)
  names(z) <- c("didx", "n")
  z$area_ha <- z$n * cell_ha
  z[, c("didx", "area_ha")]
}

rows <- list()
for (yr in years) {
  u <- disturbed_ha_by_district(path_unconstr(yr)); names(u)[2] <- "area_unconstr_ha"
  c_ <- disturbed_ha_by_district(path_constr(yr));  names(c_)[2] <- "area_constr_ha"
  m <- merge(u, c_, by = "didx"); m$year <- yr
  rows[[as.character(yr)]] <- m
  cat(sprintf("  %d done\n", yr))
}

wp <- bind_rows(rows) |>
  left_join(st_drop_geometry(districts)[, c("didx", "GEOLEVEL2", "ADMIN_NAME")],
            by = "didx") |>
  arrange(ADMIN_NAME, year) |>
  group_by(ADMIN_NAME) |>
  mutate(annual_change_ha = area_constr_ha - lag(area_constr_ha)) |>   # on constrained
  ungroup()

# per-district z-score of annual change
wp <- wp |>
  group_by(ADMIN_NAME) |>
  mutate(chg_mean = mean(annual_change_ha, na.rm = TRUE),
         chg_sd   = sd(annual_change_ha,   na.rm = TRUE),
         z_change = (annual_change_ha - chg_mean) / chg_sd) |>
  ungroup() |>
  select(GEOLEVEL2, ADMIN_NAME, year,
         area_unconstr_ha, area_constr_ha, annual_change_ha, z_change,
         chg_mean, chg_sd)

write_csv(wp, file.path(stat_dir, "district_mining_within_polygons.csv"))
cat(sprintf("\nWithin-polygon panel -> %s\n",
            file.path(stat_dir, "district_mining_within_polygons.csv")))

# ---- LaTeX appendix table per selected district ---------------
fmt <- function(x) ifelse(is.na(x), "--",
                          formatC(round(x), format = "d", big.mark = "\\,"))
fmt_z <- function(z) {
  if (is.na(z)) return("--")
  s <- sprintf("%.2f", z)
  if (z > Z_THRESHOLD) sprintf("\\textbf{%s}", s) else s
}

for (i in seq_along(DISTRICTS_FOR_TABLES)) {
  label   <- names(DISTRICTS_FOR_TABLES)[i]
  pattern <- DISTRICTS_FOR_TABLES[i]
  d <- wp |> filter(grepl(pattern, ADMIN_NAME, ignore.case = TRUE)) |> arrange(year)
  if (nrow(d) == 0) { cat(sprintf("  [skip] no district matches '%s'\n", pattern)); next }
  
  full_name <- d$ADMIN_NAME[1]
  xbar <- round(d$chg_mean[1]); ssd <- round(d$chg_sd[1])
  slug <- gsub("[^A-Za-z0-9]+", "_", label)
  
  body <- sprintf("%d & %s & %s & %s & %s \\\\",
                  d$year, fmt(d$area_unconstr_ha), fmt(d$area_constr_ha),
                  fmt(d$annual_change_ha), vapply(d$z_change, fmt_z, character(1)))
  
  tex <- c(
    "\\begin{table}[!ht]", "\\centering",
    sprintf("\\caption{Estimated mining area within polygons, %s district, %d--%d. Mining area combines mining and water class pixels. Constrained estimate applies the transition rule. z-scores computed against district mean of annual change ($\\bar{x} = %d$ ha, $s = %d$ ha). Years with $z > %.1f$ in bold.}",
            full_name, min(d$year), max(d$year), xbar, ssd, Z_THRESHOLD),
    sprintf("\\label{tab:appendix_%s}", slug),
    "\\begin{tabular}{lrrrr}", "\\toprule",
    "Year & Unconstr. (ha) & Constr. (ha) & Annual change (ha) & z-score \\\\",
    "\\midrule",
    body,
    "\\bottomrule", "\\end{tabular}", "\\end{table}"
  )
  out_p <- file.path(tab_dir, sprintf("table_appendix_%s.tex", slug))
  writeLines(tex, out_p)
  cat(sprintf("  %-16s -> %s  (2019 constr = %s ha)\n",
              label, out_p, fmt(d$area_constr_ha[d$year == max(d$year)])))
}