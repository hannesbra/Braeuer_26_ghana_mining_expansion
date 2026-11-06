# ================================================================
# District aggregation of mining changes for time and extent of mining expansion
# ----------------------------------------------------------------
# Aggregates per-year mining rasters into a DISTRICT x YEAR panel, then derives onset timing and expansion
# magnitude per district. Links to socio-economic analyses. Mining pixels inside mining polygons counted only
#
# Inputs  : PATHS$out_constrained/predictions/prediction_{year}_shrinking.tif
#           PATHS$out_unconstrained/predictions/prediction_{year}.tif
#           PATHS$mining_polygons, PATHS$ghana_districts
# Outputs : P5_RESULTS/statistics/district_mining_panel.csv    (district x year, both series + share)
#           <stat_dir>/district_mining_summary.csv  (onset, magnitude per district)
#           <tab_dir>/table_district_expansion.tex  (top expanders, booktabs)
# ================================================================

source("R/00_config.R")
source("R/01_utils.R")

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(readr)
})

# ---- Config for this step -------------------------------------
ONSET_MIN_HA <- 10   # onset = first year constrained within-polygon area >= this (ha)
MINING_CLASS <- CLASS_IDS[["Mining"]]   # = 2
WATER_CLASS  <- CLASS_IDS[["Water"]]    # = 1  (combined with mining)
years        <- PREDICT_CONFIG$predict_years            # 2000:2019

stat_dir <- results_stat_dir()   # P5_RESULTS/statistics (key results)
tab_dir  <- results_tab_dir()    # P5_RESULTS/tables

# ---- Path resolvers (need BOTH series per year) ---------------
path_constr <- function(yr)
  file.path(PATHS$out_constrained, "predictions", sprintf("prediction_%d_shrinking.tif", yr))
path_unconstr <- function(yr)
  file.path(PATHS$national_predictions, sprintf("prediction_%d.tif", yr))

years <- years[file.exists(vapply(years, path_constr,   character(1))) &
                 file.exists(vapply(years, path_unconstr, character(1)))]
stopifnot(length(years) > 0)
cat(sprintf("Aggregating %d yearly rasters (%d-%d), within polygons, Mining+Water\n",
            length(years), min(years), max(years)))

# ---- Reference grid, district id raster, polygon inside-mask --
ref <- rast(path_constr(years[1]))
cell_ha <- as.numeric(global(cellSize(ref, unit = "ha"), "mean", na.rm = TRUE)[1, 1])

districts <- st_read(PATHS$ghana_districts, quiet = TRUE) |>
  st_transform(terra::crs(ref))
districts$didx <- seq_len(nrow(districts))
dist_r <- terra::rasterize(vect(districts), ref, field = "didx")
names(dist_r) <- "didx"

inside      <- rasterize_polygons(PATHS$mining_polygons, ref)   # util from 01_utils.R
inside_bool <- inside == 1L

# District total area + mining-polygon area within each district (ha)
dist_cells <- terra::zonal(cellSize(ref, unit = "ha"), dist_r, fun = "sum", na.rm = TRUE)
names(dist_cells) <- c("didx", "district_area_ha")
poly_cells <- terra::zonal(inside, dist_r, fun = "sum", na.rm = TRUE)
names(poly_cells) <- c("didx", "poly_cells")
poly_cells$polygon_area_ha <- poly_cells$poly_cells * cell_ha

# ---- Within-polygon disturbed (Mining+Water) ha per district --
# (r == Mining | r == Water) preserves NA, so a district whose polygon
# pixels are all NA (no valid data) yields NaN -> set to NA. Districts
# with no polygons yield 0 (all inside_bool FALSE), which is correct.
disturbed_inside_ha <- function(path) {
  r  <- rast(path)
  di <- (r == MINING_CLASS | r == WATER_CLASS) & inside_bool
  z  <- terra::zonal(di, dist_r, fun = "sum", na.rm = TRUE)
  names(z) <- c("didx", "n")
  ha <- z$n * cell_ha
  ha[is.nan(ha)] <- NA
  data.frame(didx = z$didx, ha = ha)
}

panel_rows <- vector("list", length(years))
for (i in seq_along(years)) {
  yr <- years[i]
  c_ <- disturbed_inside_ha(path_constr(yr));   names(c_)[2] <- "mining_ha"
  u_ <- disturbed_inside_ha(path_unconstr(yr)); names(u_)[2] <- "mining_ha_unconstr"
  z  <- merge(c_, u_, by = "didx")
  z$year <- yr
  panel_rows[[i]] <- z[, c("didx", "year", "mining_ha", "mining_ha_unconstr")]
  cat(sprintf("  %d : %s ha (constr, in-polygon Mining+Water)\n",
              yr, format(round(sum(z$mining_ha, na.rm = TRUE)), big.mark = ",")))
  rm(c_, u_, z); gc(verbose = FALSE)
}

panel <- bind_rows(panel_rows) |>
  left_join(st_drop_geometry(districts)[, c("didx", "GEOLEVEL2", "ADMIN_NAME")], by = "didx") |>
  left_join(dist_cells, by = "didx") |>
  left_join(poly_cells[, c("didx", "polygon_area_ha")], by = "didx") |>
  mutate(mining_share      = mining_ha / district_area_ha,     # mined fraction of district
         mining_share_poly = mining_ha / polygon_area_ha) |>   # fill of known mining area
  arrange(GEOLEVEL2, year) |>
  select(GEOLEVEL2, ADMIN_NAME, year, mining_ha, mining_ha_unconstr,
         district_area_ha, polygon_area_ha, mining_share, mining_share_poly)

write_csv(panel, file.path(stat_dir, "district_mining_panel.csv"))
cat(sprintf("\nPanel -> %s  (%d district-years)\n",
            file.path(stat_dir, "district_mining_panel.csv"), nrow(panel)))

# ---- District-level onset + magnitude (constrained series) ----
y0 <- min(years); y1 <- max(years)
summary_tbl <- panel |>
  group_by(GEOLEVEL2, ADMIN_NAME, district_area_ha, polygon_area_ha) |>
  summarise(
    area_first_ha         = mining_ha[year == y0],
    area_last_ha          = mining_ha[year == y1],
    area_last_unconstr_ha = mining_ha_unconstr[year == y1],
    peak_ha   = if (all(is.na(mining_ha))) NA_real_ else max(mining_ha, na.rm = TRUE),
    peak_year = if (all(is.na(mining_ha))) NA_integer_
    else year[which.max(replace(mining_ha, is.na(mining_ha), -Inf))],
    onset_year = {
      hit <- which(mining_ha >= ONSET_MIN_HA)
      if (length(hit)) year[min(hit)] else NA_integer_
    },
    .groups = "drop"
  ) |>
  mutate(
    expansion_ha    = area_last_ha - area_first_ha,
    expansion_share = expansion_ha / district_area_ha,
    annual_exp_ha   = expansion_ha / (y1 - y0),
    ever_mined      = !is.na(onset_year)
  ) |>
  arrange(desc(expansion_ha))

write_csv(summary_tbl, file.path(stat_dir, "district_mining_summary.csv"))
cat(sprintf("Summary -> %s  (%d districts, %d with onset)\n",
            file.path(stat_dir, "district_mining_summary.csv"),
            nrow(summary_tbl), sum(summary_tbl$ever_mined)))

cat("\nTop 10 districts by mining expansion (ha, constrained, in-polygon):\n")
print(summary_tbl |>
        select(ADMIN_NAME, onset_year, area_first_ha, area_last_ha,
               expansion_ha, expansion_share) |>
        head(10), n = 10)

# ---- LaTeX table (booktabs), top expanders --------------------
topN <- head(summary_tbl, 15)
tex <- c(
  "\\begin{table}[!ht]", "\\centering",
  "\\caption{District-level mining expansion within polygons, 2000--2019 (Mining+Water, constrained).}",
  "\\label{tab:district_expansion}",
  "\\begin{tabular}{lrrrr}", "\\toprule",
  "District & Onset & 2000 (ha) & 2019 (ha) & $\\Delta$ (ha) \\\\",
  "\\midrule",
  sprintf("%s & %s & %.0f & %.0f & %.0f \\\\",
          topN$ADMIN_NAME,
          ifelse(is.na(topN$onset_year), "--", as.character(topN$onset_year)),
          topN$area_first_ha, topN$area_last_ha, topN$expansion_ha),
  "\\bottomrule", "\\end{tabular}", "\\end{table}"
)
writeLines(tex, file.path(tab_dir, "table_district_expansion.tex"))
cat(sprintf("\nLaTeX table -> %s\n",
            file.path(tab_dir, "table_district_expansion.tex")))