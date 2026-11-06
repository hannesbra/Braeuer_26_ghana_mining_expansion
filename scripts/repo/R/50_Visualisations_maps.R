# ================================================================
# 50_Visualisations_maps.R
# ----------------------------------------------------------------
# Figures for the Ghana mining analysis (national model).
#
# Produces (in P5_RESULTS/figures/50_maps; caches in P4_VALIDATION/viz_cache):
#
#   A. fig_overlay_4panel.png         Mining + water cells over 2000/2019
#                                     imagery (Tarkwa, Amansie West)
#   B. fig_onset_map.png              Median first year + IQR contours,
#                                     with hillshade base and 5 largest
#                                     cities labelled
#   C. fig_methods_panel.png          (a) national CV folds, (b) variogram
#                                     behind the CV buffer
#   D. fig_aoa.png                    Area of applicability, national model
#   E. fig_constraint_vs_gold.png     Yearly mining area: constrained
#                                     vs unconstrained, LBMA gold price,
#                                     Landsat sensor changes (dotted)
#   F. fig_det_vs_prob.png            Deterministic (P2) vs probabilistic
#                                     (P3 mean, +/- 1 SD, 95% interval)
#   G. fig_study_area.png             Study area overview: terrain backdrop
#                                     (elevation tint + hillshade), districts,
#                                     lakes, main rivers, mining polygons,
#                                     cities; legend on top
#   H. fig_industrial_mines_appendix  4 mines x 3 years panel map
#   I. fig_valid_pixel_count.png      Valid Landsat observations per
#                                     composite year and zone (mining
#                                     polygons + 5 km buffer)
#
# Colouring in A and H (no legends; set in LaTeX):
#   Water (inside and outside the polygons) = blue
#   Mining inside a Maus polygon = yellow, outside = dark red
#   H only: vegetation = green
#   A and H: no data inside the study footprint = grey; outside = white
#
# Speed: the national model and its training data are loaded only when
# C or D is switched on in RUN below.
#
# v7: A/H colour scheme with study-footprint masking (zero-reflectance
# cells outside the composites are no longer classed as water); G with
# terrain backdrop and Natural Earth Africa rivers; section switches.
# ================================================================

source("R/00_config.R")
source("R/01_utils.R")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  terra, sf, dplyr, tidyr, readr, ggplot2, patchwork,
  scales, viridisLite, automap, gstat, ggnewscale, ggspatial
)
HAS_TIDYTERRA <- requireNamespace("tidyterra", quietly = TRUE)
if (HAS_TIDYTERRA) library(tidyterra)
HAS_CAST     <- requireNamespace("CAST",     quietly = TRUE)
HAS_GGREPEL  <- requireNamespace("ggrepel",  quietly = TRUE)
if (HAS_GGREPEL) library(ggrepel)

terraOptions(memfrac = 0.6, tempdir = tempdir())
sf::sf_use_s2(FALSE)
set.seed(CV_CONFIG$seed)

# Which figures to draw (C and D load the national model: slow)
RUN <- c(A = TRUE, B = TRUE, C = TRUE, D = TRUE, E = TRUE,
         F = TRUE, G = TRUE, H = TRUE, I = TRUE)

`%||%` <- function(a, b)
  if (is.null(a) || (length(a) == 1 && is.atomic(a) && is.na(a))) b else a

# A cached CSV is reused only if it is newer than every input raster;
# after a re-prediction (10/20/30) it is rebuilt automatically.
cache_fresh <- function(cache, inputs) {
  inputs <- inputs[file.exists(inputs)]
  file.exists(cache) && length(inputs) > 0 &&
    file.mtime(cache) > max(file.mtime(inputs))
}

out_fig   <- results_fig_dir("50_maps")                       # paper figures
cache_dir <- file.path(PATHS$out_validation, "viz_cache")      # rebuildable caches
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat(" VISUALISATIONS  v7\n")
cat("================================================================\n")
cat(sprintf("Output: %s\n", out_fig))
cat(sprintf("Sections: %s\n", paste(names(RUN)[RUN], collapse = ", ")))

# ----------------------------------------------------------------
# Class colours for A and H
# ----------------------------------------------------------------
COL_WATER      <- unname(CLASS_COLS["Water"])   # "#0072B2"  water inside polygons
COL_MINING     <- "#FFD300"                     # mining inside polygons
COL_WATER_OUT  <- COL_WATER                     # water outside polygons (same blue)
COL_MINING_OUT <- "#8B0000"                     # mining outside polygons
COL_VEG        <- unname(CLASS_COLS["Vegetation"])
COL_NODATA     <- "grey70"

# Class codes shared by A and H
#   1 water in, 2 mining in, 3 water out, 4 mining out, 5 vegetation, 6 no data
CLS_LEVELS <- c("water_in", "mining_in", "water_out", "mining_out", "veg", "nodata")
CLS_COLS   <- c(water_in = COL_WATER, mining_in = COL_MINING,
                water_out = COL_WATER_OUT, mining_out = COL_MINING_OUT,
                veg = COL_VEG, nodata = COL_NODATA)

WATER_ID  <- CLASS_IDS[["Water"]]
MINING_ID <- CLASS_IDS[["Mining"]]

# ----------------------------------------------------------------
# Themes
# ----------------------------------------------------------------
theme_pub <- function(base = 10) {
  theme_bw(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
      strip.background = element_rect(fill = "grey96", colour = NA),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", size = base + 1),
      plot.subtitle    = element_text(colour = "grey30", size = base - 1)
    )
}
theme_map <- function(base = 9) {
  theme_void(base_size = base) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5,
                                      margin = margin(b = 4)),
      plot.subtitle    = element_text(hjust = 0.5, colour = "grey30",
                                      margin = margin(b = 4)),
      plot.background  = element_rect(fill = "white", colour = NA),
      legend.position  = "bottom",
      legend.key.height = unit(0.3, "cm"),
      legend.key.width  = unit(1.0, "cm"),
      legend.title      = element_text(size = base - 1),
      legend.text       = element_text(size = base - 1)
    )
}

# ----------------------------------------------------------------
# Map furniture: lat/lon grid (graticule), scale bar, north arrow
# ----------------------------------------------------------------
# theme_map_grid() draws the graticule ON TOP of the rasters
# (panel.ontop) so it stays visible over full-coverage imagery, and
# prints degree labels on the axes.
theme_map_grid <- function(base = 9, dark = FALSE, border_col = "grey20") {
  grid_col <- if (dark) scales::alpha("white", 0.45) else scales::alpha("grey30", 0.45)
  theme_bw(base_size = base) +
    theme(
      panel.ontop       = TRUE,
      panel.background  = element_rect(fill = NA, colour = NA),
      panel.grid.major  = element_line(colour = grid_col, linewidth = 0.25,
                                       linetype = "dashed"),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(fill = NA, colour = border_col,
                                       linewidth = 0.5),
      plot.background   = element_rect(fill = "white", colour = NA),
      axis.title        = element_blank(),
      axis.text         = element_text(size = base - 2, colour = "grey25"),
      axis.text.y       = element_text(angle = 90, hjust = 0.5),
      axis.ticks        = element_line(colour = "grey40", linewidth = 0.25),
      axis.ticks.length = unit(0.08, "cm"),
      plot.title        = element_text(face = "bold", hjust = 0.5,
                                       margin = margin(b = 4)),
      plot.subtitle     = element_text(hjust = 0.5, colour = "grey30",
                                       margin = margin(b = 4)),
      legend.position   = "bottom",
      legend.key.height = unit(0.3, "cm"),
      legend.key.width  = unit(1.0, "cm"),
      legend.title      = element_text(size = base - 1),
      legend.text       = element_text(size = base - 1)
    )
}

# Tidy graticule breaks (in degrees) for a UTM view window.
grat_scales <- function(xlim, ylim, n = 4, crs = GHANA_CRS) {
  xlim <- unname(as.numeric(xlim)); ylim <- unname(as.numeric(ylim))
  bb <- sf::st_bbox(c(xmin = xlim[1], ymin = ylim[1],
                      xmax = xlim[2], ymax = ylim[2]),
                    crs = sf::st_crs(crs))
  ll <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(bb), 4326))
  list(
    scale_x_continuous(breaks = pretty(c(ll[["xmin"]], ll[["xmax"]]), n = n)),
    scale_y_continuous(breaks = pretty(c(ll[["ymin"]], ll[["ymax"]]), n = n))
  )
}

# Scale bar + north arrow. `dark = TRUE` switches to white text and
# outlines for use over dark imagery; `size` scales both elements.
map_furniture <- function(dark = FALSE, scale_loc = "bl", arrow_loc = "tr",
                          size = 1, width_hint = 0.25) {
  fg <- if (dark) "white" else "black"
  list(
    ggspatial::annotation_scale(
      location = scale_loc, width_hint = width_hint, style = "bar",
      unit_category = "metric",
      bar_cols = c("black", "white"), line_col = fg, text_col = fg,
      height = unit(0.18 * size, "cm"), text_cex = 0.7 * size,
      pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm")
    ),
    ggspatial::annotation_north_arrow(
      location = arrow_loc, which_north = "true",
      height = unit(0.9 * size, "cm"), width = unit(0.9 * size, "cm"),
      pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"),
      style = ggspatial::north_arrow_fancy_orienteering(
        line_col = fg, text_col = fg,
        fill = if (dark) c("black", "white") else c("white", "black"),
        text_size = 7 * size)
    )
  )
}

# Panel letters for paper figures: large bold "(a)", "(b)"; captions
# and descriptive titles are added in LaTeX. `base` sets all text sizes.
theme_paper_panel <- function(base = 14) theme(
  plot.title    = element_text(size = base + 8, face = "bold", hjust = 0),
  plot.subtitle = element_blank(),
  plot.caption  = element_blank(),
  legend.title  = element_text(size = base + 2),
  legend.text   = element_text(size = base),
  axis.text     = element_text(size = base - 1),
  axis.title    = element_text(size = base + 2))

# ================================================================
# 0. CRS COMPATIBILITY CHECK
# ================================================================
assert_utm30n <- function(r, label) {
  if (is.character(r)) {
    if (!file.exists(r)) {
      message(sprintf("  [crs] %s: missing (%s)", label, r))
      return(invisible(FALSE))
    }
    r <- rast(r)
  }
  desc <- terra::crs(r, describe = TRUE)
  ok   <- isTRUE(desc$code == "32630")
  if (!ok) {
    warning(sprintf(
      "[crs] %s is NOT EPSG:32630 (got authority=%s, code=%s, name=%s).",
      label, desc$authority %||% "NA", desc$code %||% "NA", desc$name %||% "NA"
    ), call. = FALSE)
  } else {
    cat(sprintf("  [crs] %s: EPSG:32630 OK\n", label))
  }
  invisible(ok)
}

cat("\n[crs] Auditing inputs\n")
cat("  Configured GHANA_CRS:", GHANA_CRS, "\n")

p1_dir <- PATHS$out_unconstrained
p2_dir <- PATHS$out_constrained
p3_dir <- PATHS$out_uncertainty

comp_2019_path <- file.path(PATHS$composite_dir, "Mining_Composite_2019.tif")
comp_2000_path <- file.path(PATHS$composite_dir, "Mining_Composite_2000.tif")
pred_2019_path <- file.path(p2_dir, "predictions", "prediction_2019_shrinking.tif")
pred_2000_path <- file.path(p2_dir, "predictions", "prediction_2000_shrinking.tif")
if (!file.exists(pred_2019_path))
  pred_2019_path <- file.path(PATHS$national_predictions, "prediction_2019.tif")
if (!file.exists(pred_2000_path))
  pred_2000_path <- file.path(PATHS$national_predictions, "prediction_2000.tif")

assert_utm30n(comp_2019_path, "composite_2019")
assert_utm30n(comp_2000_path, "composite_2000")
assert_utm30n(pred_2019_path, "prediction_2019")
assert_utm30n(pred_2000_path, "prediction_2000")
assert_utm30n(file.path(p3_dir, "predictions", "onset_median.tif"), "onset_median")
assert_utm30n(file.path(p3_dir, "predictions", "onset_iqr.tif"),    "onset_iqr")

read_and_warp <- function(path, label, target_crs = GHANA_CRS) {
  if (!file.exists(path)) {
    message(sprintf("  [crs] %s: missing (%s)", label, path)); return(NULL)
  }
  v <- sf::st_read(path, quiet = TRUE)
  src <- sf::st_crs(v)$epsg %||% NA
  if (is.na(src)) {
    message(sprintf("  [crs] %s: undefined CRS - assuming WGS84", label))
    sf::st_crs(v) <- 4326
  }
  v <- sf::st_transform(v, target_crs)
  cat(sprintf("  [crs] %s: src=%s -> %s OK\n", label, src %||% "NA", target_crs))
  v
}

ghana_v       <- read_and_warp(PATHS$ghana_districts, "ghana_districts")
mining_p      <- read_and_warp(PATHS$mining_polygons, "mining_polygons")
if (!is.null(mining_p)) {
  mining_p <- sf::st_make_valid(mining_p)
  mining_p <- mining_p[!sf::st_is_empty(mining_p), ]
}
ghana_outline <- if (!is.null(ghana_v))
  suppressWarnings(sf::st_union(ghana_v)) else NULL

# ----------------------------------------------------------------
# Raster helpers for A and H
# ----------------------------------------------------------------
pred_path_for <- function(year) {
  p <- file.path(p2_dir, "predictions", sprintf("prediction_%d_shrinking.tif", year))
  if (!file.exists(p))
    p <- file.path(PATHS$national_predictions, sprintf("prediction_%d.tif", year))
  p
}
comp_path_for <- function(year) sprintf(NATIONAL_CONFIG$year_composite_fmt, year)

# RGB of a composite, cropped (reads only 3 layers)
read_rgb <- function(year, bx) {
  r <- rast(comp_path_for(year), lyrs = match(c("RED", "GREEN", "BLUE"), BAND_NAMES))
  names(r) <- c("RED", "GREEN", "BLUE")
  crop(r, bx)
}

# 1 where the composite has real data, NA elsewhere. Outside the
# composite footprint the bands are NA or exactly 0; the classifier
# turns zero reflectance into "water", so those cells must be masked.
valid_mask <- function(rgb) {
  s <- rgb$RED + rgb$GREEN + rgb$BLUE
  ifel(!is.na(s) & s > 0, 1L, NA)
}

align_to <- function(r, tmpl) {
  if (compareGeom(r, tmpl, stopOnError = FALSE)) r
  else resample(r, tmpl, method = "near")
}

bbox_sf <- function(bx) sf::st_as_sfc(sf::st_bbox(
  c(xmin = xmin(bx), ymin = ymin(bx), xmax = xmax(bx), ymax = ymax(bx)),
  crs = GHANA_CRS))

# Inside-polygon mask (1 = cell centre in a Maus polygon, 0 = outside)
# on the grid of raster `r`, using only the polygons that touch `bx_sf`.
inside_mask_for <- function(r, bx_sf) {
  if (is.null(mining_p)) return(terra::init(r[[1]], 0))
  sub <- mining_p[lengths(sf::st_intersects(mining_p, bx_sf)) > 0, ]
  if (!nrow(sub)) return(terra::init(r[[1]], 0))
  rasterize_polygons(sub, r)
}

# Class raster (codes as CLS_LEVELS). Vegetation -> veg_code
# (NA in A, so the imagery shows through). Cells inside the study
# footprint without data this year -> 6; outside the footprint -> NA.
classify_cells <- function(pred, inside, valid, study, veg_code = 5L) {
  cls <- ifel(pred == WATER_ID,  ifel(inside, 1L, 3L),
              ifel(pred == MINING_ID, ifel(inside, 2L, 4L), veg_code))
  cls <- mask(cls, valid)
  nod <- !is.na(study) & (is.na(valid) | is.na(pred))
  cls <- ifel(nod, 6L, cls)
  mask(cls, study)
}

to_utm_xy <- function(lon, lat)
  sf::st_coordinates(sf::st_transform(
    sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326), GHANA_CRS))

# ----------------------------------------------------------------
# Reference data: major Ghanaian cities (2010-2021 census,
# urban agglomeration where reported).
# ----------------------------------------------------------------
GHANA_CITIES <- data.frame(
  name = c("Accra", "Kumasi", "Tamale", "Sekondi-Takoradi", "Cape Coast",
           "Sunyani", "Obuasi", "Koforidua", "Ho",  "Tarkwa",
           "Bibiani", "Bogoso", "Wa", "Bolgatanga", "Asankrangwa"),
  lon  = c(-0.1870, -1.6244, -0.8393, -1.7554, -1.2466,
           -2.3265, -1.6700, -0.2571,  0.4719, -1.9933,
           -2.3306, -2.0167, -2.5099, -0.8516, -2.4264),
  lat  = c( 5.5600,  6.6906,  9.4035,  4.9344,  5.1054,
            7.3382,  6.2027,  6.0894,  6.6111,  5.3000,
            6.4636,  5.5722, 10.0601, 10.7856,  5.7997),
  pop  = c(2500000, 2000000, 360000, 445000, 170000,
           125000,  175000, 125000, 110000,  36000,
           35000,   27000, 100000,  75000,  20000)
)
cities_sf <- sf::st_as_sf(GHANA_CITIES, coords = c("lon", "lat"),
                          crs = 4326) |>
  sf::st_transform(GHANA_CRS)

# ----------------------------------------------------------------
# National model and its training points with CV fold ids (C, D only)
# ----------------------------------------------------------------
national_model <- NULL
train_nat      <- NULL
if (RUN[["C"]] || RUN[["D"]]) {
  national_model <- if (file.exists(NATIONAL_CONFIG$model_path))
    load_national_model() else NULL
  if (!is.null(national_model)) {
    cat("\n[model] national model + training folds\n")
    train_nat <- prepare_training_data(national_model$train_composite, verbose = FALSE)
    j <- match(train_nat$row_id, national_model$cv$row_id)
    train_nat$fold <- national_model$cv$fold_ids[j]
    train_nat <- train_nat[!is.na(j), , drop = FALSE]
    cat(sprintf("  %d training points, bands: %s\n", nrow(train_nat),
                paste(national_model$features, collapse = ", ")))
  } else {
    cat("\n[model] national model not found - sections C and D are skipped\n")
  }
}

# ================================================================
# A. 4-PANEL OVERLAY
# ================================================================
if (RUN[["A"]]) {
  cat("\n[A] 4-panel mining-overlay maps\n")
  
  SITES <- list(
    tarkwa       = list(lon = -2.024863, lat = 5.329213,
                        label = "Tarkwa  (industrial)"),
    amansie_west = list(lon = -1.835355, lat = 6.366951,
                        label = "Amansie West  (artisanal)")
  )
  ZOOM_HALF_M <- 25000
  
  draw_A <- function(site, year) {
    xy  <- to_utm_xy(site$lon, site$lat)
    bx  <- ext(xy[1] - ZOOM_HALF_M, xy[1] + ZOOM_HALF_M,
               xy[2] - ZOOM_HALF_M, xy[2] + ZOOM_HALF_M)
    bsf <- bbox_sf(bx)
    
    rgb   <- read_rgb(year, bx)
    valid <- valid_mask(rgb)
    rgb   <- clamp(rgb, lower = 0, upper = 0.4, values = TRUE)
    gray  <- mask(0.299 * rgb$RED + 0.587 * rgb$GREEN + 0.114 * rgb$BLUE, valid)
    names(gray) <- "luminance"
    qs   <- as.numeric(global(gray, fun = function(x)
      quantile(x, c(0.02, 0.98), na.rm = TRUE))[1, ])
    gray <- clamp(gray, qs[1], qs[2], values = TRUE)
    
    pred   <- align_to(crop(rast(pred_path_for(year)), bx), gray)
    study  <- align_to(valid_mask(read_rgb(2019, bx)), gray)   # composite footprint
    inside <- inside_mask_for(pred, bsf) == 1
    cls <- classify_cells(pred, inside, valid, study, veg_code = NA)
    odf <- as.data.frame(cls, xy = TRUE, na.rm = TRUE)
    names(odf)[3] <- "type"
    odf$type <- factor(CLS_LEVELS[odf$type], levels = CLS_LEVELS)
    
    maus <- suppressWarnings(sf::st_intersection(mining_p, bsf))
    
    g <- ggplot()
    g <- if (HAS_TIDYTERRA) g + tidyterra::geom_spatraster(data = gray) else
      g + geom_raster(data = as.data.frame(gray, xy = TRUE, na.rm = TRUE),
                      aes(x, y, fill = luminance))
    g <- g + scale_fill_gradient(low = "black", high = "white", limits = qs,
                                 na.value = NA, guide = "none") +
      ggnewscale::new_scale_fill()
    if (nrow(odf))
      g <- g + geom_raster(data = odf, aes(x, y, fill = type)) +
      scale_fill_manual(values = CLS_COLS, drop = FALSE, guide = "none")
    if (nrow(maus))
      g <- g + geom_sf(data = maus, fill = scales::alpha("white", 0.18),
                       colour = scales::alpha("white", 0.65), linewidth = 0.25)
    xl <- c(xmin(bx), xmax(bx)); yl <- c(ymin(bx), ymax(bx))
    g + map_furniture(dark = TRUE) + grat_scales(xl, yl) +
      coord_sf(crs = GHANA_CRS, datum = sf::st_crs(4326),
               xlim = xl, ylim = yl, expand = FALSE) +
      labs(title = sprintf("%s  -  %d", site$label, year)) +
      theme_map_grid(dark = TRUE) + theme(legend.position = "none")
  }
  
  pa <- list()
  for (s in names(SITES)) for (yr in c(2000, 2019)) {
    cat(sprintf("  building %s_%d\n", s, yr))
    pa[[paste0(s, "_", yr)]] <- tryCatch(draw_A(SITES[[s]], yr), error = function(e) {
      message("    failed: ", e$message); NULL })
  }
  if (length(Filter(Negate(is.null), pa)) == 4) {
    fig_A <- (pa$tarkwa_2000 | pa$tarkwa_2019) /
      (pa$amansie_west_2000 | pa$amansie_west_2019)
    ggsave(file.path(out_fig, "fig_overlay_4panel.png"), fig_A,
           width = 10, height = 10, dpi = 500, bg = "white")
    cat(sprintf("  -> %s\n", file.path(out_fig, "fig_overlay_4panel.png")))
  } else cat("  Skipped: not all 4 panels could be built.\n")
}

# ================================================================
# B. ONSET MAP
# ================================================================
onset_median_path <- file.path(p3_dir, "predictions", "onset_median.tif")
onset_iqr_path    <- file.path(p3_dir, "predictions", "onset_iqr.tif")

if (RUN[["B"]]) {
  cat("\n[B] Onset map\n")
  if (file.exists(onset_median_path) && file.exists(onset_iqr_path)) {
    
    onset_med <- rast(onset_median_path)
    onset_iqr <- rast(onset_iqr_path)
    frac_path <- file.path(p3_dir, "predictions", "onset_frac_disturbed.tif")
    ONSET_MIN_FRAC <- 0.5
    if (file.exists(frac_path)) {
      frac <- rast(frac_path)
      onset_med <- mask(onset_med, frac >= ONSET_MIN_FRAC, maskvalues = c(0, NA))
      onset_iqr <- mask(onset_iqr, frac >= ONSET_MIN_FRAC, maskvalues = c(0, NA))
    }
    
    agg_fac <- max(1, floor(min(nrow(onset_med), ncol(onset_med)) / 1500))
    if (agg_fac > 1) {
      onset_med_p <- aggregate(onset_med, fact = agg_fac,
                               fun = "median", na.rm = TRUE)
      onset_iqr_p <- aggregate(onset_iqr, fact = agg_fac,
                               fun = "median", na.rm = TRUE)
    } else {
      onset_med_p <- onset_med; onset_iqr_p <- onset_iqr
    }
    
    valid_cells <- which(!is.na(values(onset_med_p)[, 1]))
    if (length(valid_cells) > 0) {
      pts <- terra::xyFromCell(onset_med_p, valid_cells)
      qx <- as.numeric(quantile(pts[, 1], probs = c(0.01, 0.99), na.rm = TRUE))
      qy <- as.numeric(quantile(pts[, 2], probs = c(0.01, 0.99), na.rm = TRUE))
      pad_m <- 5000
      bx_view <- ext(qx[1] - pad_m, qx[2] + pad_m,
                     qy[1] - pad_m, qy[2] + pad_m)
      cat(sprintf("  view  width = %.1f km, height = %.1f km\n",
                  (xmax(bx_view) - xmin(bx_view)) / 1000,
                  (ymax(bx_view) - ymin(bx_view)) / 1000))
    } else {
      bx_view <- ext(onset_med_p)
    }
    
    onset_med_c <- crop(onset_med_p, bx_view)
    onset_iqr_c <- crop(onset_iqr_p, bx_view)
    
    hillshade <- NULL
    if (file.exists(comp_2019_path)) {
      cat("  building hillshade\n")
      elev_lyr <- which(BAND_NAMES == "ELEVATION")
      elev <- tryCatch(
        rast(comp_2019_path, lyrs = elev_lyr) |> crop(bx_view),
        error = function(e) { message("    elevation crop failed: ", e$message); NULL }
      )
      if (!is.null(elev)) {
        elev <- aggregate(elev, fact = 10, fun = "mean", na.rm = TRUE)
        slope  <- terrain(elev, "slope",  unit = "radians")
        aspect <- terrain(elev, "aspect", unit = "radians")
        hillshade <- shade(slope, aspect, angle = 35, direction = 315)
        names(hillshade) <- "hillshade"
      }
    }
    
    view_box_geom <- sf::st_as_sfc(sf::st_bbox(c(
      xmin = xmin(bx_view), ymin = ymin(bx_view),
      xmax = xmax(bx_view), ymax = ymax(bx_view)
    ), crs = GHANA_CRS))
    districts_view <- if (!is.null(ghana_v))
      suppressWarnings(sf::st_intersection(ghana_v, view_box_geom)) else NULL
    
    in_view <- lengths(sf::st_intersects(cities_sf, view_box_geom)) > 0
    cities_view <- cities_sf[in_view, ]
    if (nrow(cities_view) > 5)
      cities_view <- cities_view[order(-cities_view$pop), ][1:5, ]
    cat(sprintf("  cities: %s\n",
                if (nrow(cities_view) > 0) paste(cities_view$name, collapse = ", ")
                else "(none in view)"))
    
    yr_min <- 2000; yr_max <- 2019
    
    fig_B_main <- ggplot()
    if (!is.null(ghana_outline))
      fig_B_main <- fig_B_main +
      geom_sf(data = ghana_outline, fill = "grey94", colour = "grey55",
              linewidth = 0.25)
    
    if (!is.null(hillshade) && HAS_TIDYTERRA) {
      fig_B_main <- fig_B_main +
        tidyterra::geom_spatraster(data = hillshade, alpha = 0.45) +
        scale_fill_gradient(low = "grey20", high = "grey98",
                            na.value = NA, guide = "none") +
        ggnewscale::new_scale_fill()
    }
    
    if (!is.null(districts_view) && nrow(districts_view) > 0)
      fig_B_main <- fig_B_main +
      geom_sf(data = districts_view, fill = NA,
              colour = scales::alpha("grey45", 0.55),
              linewidth = 0.18)
    
    if (HAS_TIDYTERRA) {
      fig_B_main <- fig_B_main +
        tidyterra::geom_spatraster(data = onset_med_c) +
        scale_fill_viridis_c(
          option = "C", direction = 1,
          limits = c(yr_min, yr_max),
          na.value = "transparent",
          breaks = seq(yr_min, yr_max, 5),
          labels = function(b) ifelse(b == yr_min, paste0("\u2264", yr_min), b),
          name = "Median onset year",
          guide = guide_colourbar(title.position = "top", title.hjust = 0.5)
        ) +
        tidyterra::geom_spatraster_contour(
          data = onset_iqr_c,
          breaks = c(2.5, 5, 10),
          colour = "white", alpha = 0.85, linewidth = 0.4
        )
    } else {
      df_med <- as.data.frame(onset_med_c, xy = TRUE, na.rm = TRUE)
      names(df_med)[3] <- "year"
      fig_B_main <- fig_B_main +
        geom_raster(data = df_med, aes(x, y, fill = year)) +
        scale_fill_viridis_c(
          option = "C",
          limits = c(yr_min, yr_max),
          breaks = seq(yr_min, yr_max, 5)
        )
    }
    
    if (nrow(cities_view) > 0) {
      fig_B_main <- fig_B_main +
        geom_sf(data = cities_view, shape = 21, size = 2.0, stroke = 0.5,
                fill = "white", colour = "black")
      if (HAS_GGREPEL) {
        fig_B_main <- fig_B_main +
          ggrepel::geom_label_repel(
            data = cities_view,
            aes(label = name, geometry = geometry),
            stat = "sf_coordinates",
            size = 2.7, fontface = "bold", colour = "black",
            fill = scales::alpha("white", 0.85), label.size = 0,
            box.padding = 0.4, point.padding = 0.25,
            min.segment.length = 0.1,
            segment.colour = "black", segment.size = 0.25,
            seed = 42
          )
      } else {
        fig_B_main <- fig_B_main +
          geom_sf_text(data = cities_view, aes(label = name),
                       size = 2.7, fontface = "bold", colour = "black",
                       nudge_x = 4000, nudge_y = 4000)
      }
    }
    
    xl_B <- c(xmin(bx_view), xmax(bx_view)); yl_B <- c(ymin(bx_view), ymax(bx_view))
    fig_B_main <- fig_B_main +
      map_furniture(dark = FALSE, scale_loc = "bl", arrow_loc = "tr",
                    width_hint = 0.2) +
      grat_scales(xl_B, yl_B, n = 5) +
      coord_sf(
        xlim = xl_B, ylim = yl_B,
        crs  = GHANA_CRS, datum = sf::st_crs(4326), expand = FALSE
      ) +
      theme_map_grid(base = 10) +
      theme(legend.position  = "bottom",
            legend.key.width = unit(2.5, "cm"))
    
    ggsave(file.path(out_fig, "fig_onset_map.png"),
           fig_B_main, width = 9, height = 7, dpi = 500, bg = "white")
    cat(sprintf("  -> %s\n", file.path(out_fig, "fig_onset_map.png")))
    
  } else {
    cat("  Skipped: onset rasters missing.\n")
  }
}

# ================================================================
# C. METHODS FIGURE  (national folds, variogram)
# ================================================================
if (RUN[["C"]]) {
  cat("\n[C] Methods figure\n")
  
  C_BASE <- 16   # base text size for both panels
  
  if (is.null(national_model) || is.null(train_nat)) {
    cat("  Skipped: national model missing (run 10).\n")
  } else {
    zb_utm   <- get_zone_boundaries()$utm
    train_sf <- sf::st_as_sf(train_nat, coords = c("x_utm", "y_utm"),
                             crs = GHANA_CRS, remove = FALSE)
    train_sf$class <- factor(CLASS_NAMES[train_sf$landcover_3class], levels = CLASS_NAMES)
    train_sf$fold  <- factor(train_sf$fold, levels = seq_len(CV_CONFIG$n_folds))
    
    data_geoms <- sf::st_geometry(train_sf)
    if (!is.null(mining_p)) data_geoms <- c(data_geoms, sf::st_geometry(mining_p))
    data_bbox <- sf::st_bbox(data_geoms)
    pad_m  <- 30000
    zoom_c <- list(xlim = unname(c(data_bbox["xmin"] - pad_m, data_bbox["xmax"] + pad_m)),
                   ylim = unname(c(data_bbox["ymin"] - pad_m, data_bbox["ymax"] + pad_m)))
    
    base_c <- function() {
      g <- ggplot()
      if (!is.null(ghana_outline))
        g <- g + geom_sf(data = ghana_outline, fill = "grey95", colour = "grey50",
                         linewidth = 0.3)
      for (yval in c(zb_utm$northern, zb_utm$coastal))
        g <- g + annotate("segment", x = zoom_c$xlim[1], xend = zoom_c$xlim[2],
                          y = yval, yend = yval, colour = "grey40",
                          linetype = "dashed", linewidth = 0.3)
      g
    }
    frame_c <- function(g, title) g +
      map_furniture(dark = FALSE, scale_loc = "bl", arrow_loc = "tr",
                    width_hint = 0.2, size = 1.4) +
      grat_scales(zoom_c$xlim, zoom_c$ylim, n = 4) +
      coord_sf(xlim = zoom_c$xlim, ylim = zoom_c$ylim, crs = GHANA_CRS,
               datum = sf::st_crs(4326), expand = FALSE) +
      labs(title = title) +
      theme_map_grid(base = C_BASE) + theme(legend.position = "right")
    
    # (a) national spatial folds + buffer around fold-1 test points
    cat("  C2 spatial folds\n")
    buf_m    <- national_model$buffer_used_m
    buf_poly <- sf::st_union(sf::st_buffer(train_sf[train_sf$fold == "1", ], buf_m))
    fold_pal <- setNames(scales::hue_pal()(CV_CONFIG$n_folds),
                         as.character(seq_len(CV_CONFIG$n_folds)))
    fig_C2 <- frame_c(
      base_c() +
        geom_sf(data = buf_poly, fill = "grey75", colour = "grey55",
                linewidth = 0.2, alpha = 0.5) +
        geom_sf(data = train_sf, aes(colour = fold), size = 0.8, alpha = 0.85) +
        scale_colour_manual(values = fold_pal, name = "Fold", drop = FALSE,
                            guide = guide_legend(override.aes = list(size = 5))),
      "(a)") + theme_paper_panel(C_BASE) +
      theme(legend.key.size = unit(0.9, "cm"))
    cat(sprintf("  folds: %.1f km buffer around fold 1 shown in grey; block %.1f km\n",
                buf_m / 1000, national_model$block_size_m / 1000))
    
    # (b) variogram behind the buffer
    cat("  C3 variogram\n")
    vgm <- national_model$vgm_info
    fig_C3 <- if (is.null(vgm) || is.null(vgm$fit)) {
      ggplot() + annotate("text", x = 0.5, y = 0.5, size = 5,
                          label = sprintf("no variogram fit (%s)", vgm$source %||% "fixed buffer")) +
        theme_void()
    } else {
      emp   <- as.data.frame(vgm$emp)[, c("dist", "gamma", "np")]
      d_seq <- seq(0, max(emp$dist) * 1.05, length.out = 200)
      gline <- gstat::variogramLine(vgm$fit$var_model, dist_vector = d_seq)
      rng   <- vgm$practical_range %||% buf_m
      ggplot() +
        geom_point(data = emp, aes(dist / 1000, gamma, size = np),
                   colour = "grey25", alpha = 0.75) +
        geom_line(data = gline, aes(dist / 1000, gamma), colour = "#0072B2", linewidth = 0.8) +
        geom_vline(xintercept = rng / 1000, colour = "#D55E00",
                   linetype = "dashed", linewidth = 0.6) +
        annotate("text", x = rng / 1000, y = vgm$sill * 0.05,
                 label = sprintf("%.1f km", rng / 1000),
                 hjust = -0.1, vjust = 0, size = 6.5, colour = "#D55E00") +
        scale_size_continuous(range = c(1, 4.5), name = "n pairs",
                              labels = scales::label_comma()) +
        labs(x = "Lag distance (km)", y = "Semivariance") +
        theme_pub(base = C_BASE) + theme(legend.position = "bottom")
    }
    fig_C3 <- fig_C3 + labs(title = "(b)") + theme_paper_panel(C_BASE)
    if (!is.null(vgm$fit))
      cat(sprintf("  variogram: %s model, target %s, nugget %.3f, sill %.3f, buffer %.1f km\n",
                  vgm$model %||% "NA", national_model$variogram_target %||% "NA",
                  vgm$nugget %||% NA, vgm$sill %||% NA, buf_m / 1000))
    
    fig_C <- (fig_C2 | fig_C3) + plot_layout(widths = c(1, 1.1))
    ggsave(file.path(out_fig, "fig_methods_panel.png"),
           fig_C, width = 14, height = 7.5, dpi = 500, bg = "white")
    cat(sprintf("  -> %s\n", file.path(out_fig, "fig_methods_panel.png")))
  }
}

# ================================================================
# D. AREA OF APPLICABILITY  (national model, whole country)
# ================================================================
if (RUN[["D"]]) {
  cat("\n[D] Area of applicability (national model)\n")
  
  AOA_AGG_FACT <- 10
  
  if (is.null(national_model) || is.null(train_nat)) {
    cat("  Skipped: national model missing (run 10).\n")
  } else if (!HAS_CAST) {
    cat("  Skipped: CAST not installed (install.packages('CAST')).\n")
  } else {
    feats <- national_model$features
    aoa_cache <- file.path(cache_dir, "aoa_cache.rds")
    aoa_inputs <- c(NATIONAL_CONFIG$model_path, national_model$train_composite,
                    CV_CONFIG$fold_cache)
    feats_k <- feats
    if (cache_fresh(aoa_cache, aoa_inputs)) {
      cc <- readRDS(aoa_cache)
      aoa_di <- terra::unwrap(cc$di); aoa_mask <- terra::unwrap(cc$mask)
      threshold <- cc$threshold; thr_src <- cc$thr_src
      cat(sprintf("  loaded AOA cache (%s)\n", basename(aoa_cache)))
    } else {
      comp  <- load_composite(national_model$train_composite, bands = feats)
      comp_agg <- aggregate(comp, fact = AOA_AGG_FACT, fun = "mean", na.rm = TRUE)
      rm(comp); gc(verbose = FALSE)
      
      train_x <- train_nat[, feats, drop = FALSE]
      sdv     <- vapply(train_x, stats::sd, numeric(1), na.rm = TRUE)
      feats_k <- feats[is.finite(sdv) & sdv > 1e-8]
      if (length(feats_k) < length(feats))
        cat(sprintf("  dropped constant band(s): %s\n",
                    paste(setdiff(feats, feats_k), collapse = ", ")))
      
      imp <- rowMeans(sapply(national_model$models, function(m)
        m$variable.importance[feats_k]), na.rm = TRUE)
      imp[!is.finite(imp) | imp < 0] <- 0
      if (sum(imp) <= 0) imp[] <- 1
      weight_df <- as.data.frame(t(imp)); names(weight_df) <- feats_k
      
      cv_list <- split(seq_len(nrow(train_nat)), train_nat$fold)
      thr_src <- "spatial-CV"
      tdi <- tryCatch(
        CAST::trainDI(train = train_x[, feats_k, drop = FALSE], weight = weight_df,
                      variables = feats_k, CVtest = cv_list),
        error = function(e) { message("    trainDI (CV) failed: ", conditionMessage(e)); NULL })
      if (is.null(tdi)) {
        thr_src <- "all-training"
        tdi <- CAST::trainDI(train = train_x[, feats_k, drop = FALSE], weight = weight_df,
                             variables = feats_k)
      }
      aoa_out   <- CAST::aoa(newdata = comp_agg[[feats_k]], trainDI = tdi)
      aoa_di    <- aoa_out$DI;  names(aoa_di) <- "DI"
      aoa_mask  <- aoa_out$AOA
      if (!is.null(ghana_outline)) {
        gv       <- terra::vect(sf::st_transform(ghana_outline, terra::crs(aoa_di)))
        aoa_di   <- terra::mask(aoa_di, gv)
        aoa_mask <- terra::mask(aoa_mask, gv)
      }
      threshold <- tdi$threshold
      saveRDS(list(di = terra::wrap(aoa_di), mask = terra::wrap(aoa_mask),
                   threshold = threshold, thr_src = thr_src), aoa_cache)
    }
    cat(sprintf("  AOA threshold (%s DI): %.4f\n", thr_src, threshold))
    
    aoa_vals   <- terra::values(aoa_mask)[, 1]
    ok         <- !is.na(aoa_vals)
    pct_within <- 100 * mean(aoa_vals[ok] == 1)
    pct_poly_within <- NA_real_
    if (!is.null(mining_p)) {
      pm   <- rasterize_polygons(mining_p, aoa_mask)
      psel <- terra::values(pm)[, 1] == 1L & ok
      if (any(psel)) pct_poly_within <- 100 * mean(aoa_vals[psel] == 1)
    }
    yv <- terra::init(aoa_mask, "y")
    zb_utm <- get_zone_boundaries()$utm
    zone_v <- terra::values(terra::ifel(yv >= zb_utm$northern, 1L,
                                        terra::ifel(yv >= zb_utm$coastal, 2L, 3L)))[, 1]
    aoa_zone <- data.frame(
      scope = c("national", ZONE_ORDER),
      pct_within_aoa = c(pct_within, vapply(1:3, function(z)
        100 * mean(aoa_vals[ok & zone_v == z] == 1), numeric(1))),
      pct_mining_polygons_within_aoa = c(pct_poly_within, NA, NA, NA),
      threshold = threshold, threshold_source = thr_src)
    write.csv(aoa_zone, file.path(results_stat_dir(), "aoa_summary.csv"), row.names = FALSE)
    print(aoa_zone, digits = 3, row.names = FALSE)
    
    bx_view  <- ext(trim(aoa_di))
    view_sfc <- sf::st_as_sfc(sf::st_bbox(c(xmin = xmin(bx_view), ymin = ymin(bx_view),
                                            xmax = xmax(bx_view), ymax = ymax(bx_view)),
                                          crs = GHANA_CRS))
    outline_clip <- if (is.null(ghana_outline)) NULL else
      suppressWarnings(sf::st_intersection(sf::st_make_valid(ghana_outline), view_sfc))
    maus_clip <- if (is.null(mining_p)) NULL else
      suppressWarnings(sf::st_intersection(sf::st_make_valid(mining_p), view_sfc))
    
    base_d <- function() {
      g <- ggplot()
      if (!is.null(outline_clip) && length(outline_clip) > 0)
        g <- g + geom_sf(data = outline_clip, fill = "grey95", colour = "grey60",
                         linewidth = 0.25)
      g
    }
    add_maus <- function(g, col, lwd) {
      if (!is.null(maus_clip) && nrow(maus_clip) > 0)
        g <- g + geom_sf(data = maus_clip, fill = NA, colour = col, linewidth = lwd)
      g
    }
    xl_D <- c(xmin(bx_view), xmax(bx_view)); yl_D <- c(ymin(bx_view), ymax(bx_view))
    frame_d <- function(g) g +
      map_furniture(dark = FALSE, scale_loc = "bl", arrow_loc = "tr", width_hint = 0.2) +
      grat_scales(xl_D, yl_D, n = 5) +
      coord_sf(xlim = xl_D, ylim = yl_D, crs = GHANA_CRS,
               datum = sf::st_crs(4326), expand = FALSE) +
      theme_map_grid(base = 10) + theme(legend.key.width = unit(2.0, "cm"))
    
    hi <- max(as.numeric(quantile(terra::values(aoa_di)[, 1], 0.98, na.rm = TRUE)),
              threshold * 1.15)
    di_clip  <- clamp(aoa_di, lower = 0, upper = hi, values = TRUE)
    di_guide <- guide_colourbar(title.position = "top", title.hjust = 0.5)
    figA <- base_d()
    if (HAS_TIDYTERRA) {
      figA <- figA +
        tidyterra::geom_spatraster(data = di_clip) +
        scale_fill_viridis_c(name = "Dissimilarity index (DI)", option = "rocket",
                             direction = -1, na.value = "transparent",
                             limits = c(0, hi), guide = di_guide) +
        tidyterra::geom_spatraster_contour(data = aoa_di, breaks = threshold,
                                           colour = "#00E5FF", linewidth = 0.5)
    } else {
      dA <- as.data.frame(di_clip, xy = TRUE, na.rm = TRUE); names(dA)[3] <- "DI"
      figA <- figA + geom_raster(data = dA, aes(x, y, fill = DI)) +
        scale_fill_viridis_c(name = "Dissimilarity index (DI)", option = "rocket",
                             direction = -1, limits = c(0, hi), guide = di_guide)
    }
    figA <- frame_d(add_maus(figA, "grey20", 0.12)) +
      labs(title = "(a)") + theme_paper_panel(12)
    
    dB <- as.data.frame(aoa_mask, xy = TRUE, na.rm = TRUE); names(dB)[3] <- "aoa"
    dB$AOA <- factor(ifelse(dB$aoa == 1, "Within AOA", "Outside AOA (extrapolation)"),
                     levels = c("Within AOA", "Outside AOA (extrapolation)"))
    subB <- sprintf("%.1f%% of Ghana within AOA", pct_within)
    if (!is.na(pct_poly_within))
      subB <- sprintf("%s | %.1f%% of mining polygons", subB, pct_poly_within)
    figB <- frame_d(add_maus(base_d() + geom_raster(data = dB, aes(x, y, fill = AOA)),
                             "grey15", 0.15)) +
      scale_fill_manual(values = c("Within AOA" = "#2C7FB8",
                                   "Outside AOA (extrapolation)" = "#E34A33"),
                        name = NULL, drop = FALSE) +
      labs(title = "(b)") + theme_paper_panel(12)
    cat(sprintf("  for the caption: cyan = AOA threshold (DI = %.3f, %s); %s\n",
                threshold, thr_src, subB))
    
    fig_D <- figA + figB + patchwork::plot_layout(ncol = 2)
    asp_D <- (ymax(bx_view) - ymin(bx_view)) / (xmax(bx_view) - xmin(bx_view))
    ggsave(file.path(out_fig, "fig_aoa.png"), fig_D, width = 12,
           height = min(14, max(5, (12 - 1.2) / 2 * asp_D + 1.6)),
           dpi = 500, bg = "white")
    cat(sprintf("  -> %s\n", file.path(out_fig, "fig_aoa.png")))
    
    dir.create(file.path(p3_dir, "predictions"), showWarnings = FALSE, recursive = TRUE)
    writeRaster(aoa_di, file.path(p3_dir, "predictions", "aoa_DI_national.tif"),
                overwrite = TRUE, datatype = "FLT4S")
    writeRaster(aoa_mask, file.path(p3_dir, "predictions", "aoa_mask_national.tif"),
                overwrite = TRUE, datatype = "INT1U")
  }
}

# ================================================================
# E. CONSTRAINT VS UNCONSTRAINED + GOLD PRICE + LANDSAT SENSORS
# ================================================================
if (RUN[["E"]]) {
  cat("\n[E] Constraint vs unconstrained + gold price (inside polygons)\n")
  
  E_BASE <- 16
  
  # Landsat sensor changes that affect the annual composites. Dotted
  # vertical line at the composite year, label next to it near the top.
  #   2003: Landsat 7 ETM+ scan-line corrector fails (May) -> SLC-off stripes
  #   2012: Landsat 5 TM decommissioned (Nov 2011) -> L7 SLC-off only
  #   2013: Landsat 8 OLI operational (April)
  SENSOR_EVENTS <- data.frame(
    year  = c(2003, 2012, 2013),
    label = c("L7 SLC-off", "L5 TM ends", "L8 OLI")
  )
  
  if (is.null(mining_p)) {
    cat("  Skipped: mining_polygons not available; cannot restrict.\n")
  } else {
    
    area_csv <- file.path(cache_dir, "yearly_mining_area_inside.csv")
    years_to_check <- PREDICT_CONFIG$predict_years
    
    e_inputs <- c(file.path(PATHS$national_predictions, sprintf("prediction_%d.tif", years_to_check)),
                  file.path(p2_dir, "predictions",
                            sprintf("prediction_%d_shrinking.tif", years_to_check)))
    if (!cache_fresh(area_csv, e_inputs)) {
      cat("  computing yearly mining area inside polygons (cached on first run)\n")
      
      ref <- rast(pred_path_for(years_to_check[1]))
      inside_ref <- rasterize_polygons(mining_p, ref) == 1L
      cat(sprintf("  polygon mask: %d cells inside\n",
                  as.integer(global(inside_ref, "sum", na.rm = TRUE)[1, 1])))
      
      rows <- list()
      for (yr in years_to_check) {
        for (strat in c("unconstrained", "shrinking")) {
          path <- if (strat == "unconstrained")
            file.path(PATHS$national_predictions, sprintf("prediction_%d.tif", yr))
          else
            file.path(p2_dir, "predictions",
                      sprintf("prediction_%d_shrinking.tif", yr))
          if (!file.exists(path)) {
            cat(sprintf("    [skip] %d / %s : missing\n", yr, strat)); next
          }
          r <- rast(path)
          mining_inside <- (r == 2L | r == 1L) & inside_ref   # Mining + Water
          n_mining <- as.numeric(global(mining_inside, "sum", na.rm = TRUE)[1, 1])
          cell_ha  <- as.numeric(global(cellSize(r, unit = "ha"),
                                        "mean", na.rm = TRUE)[1, 1])
          rows[[paste(yr, strat)]] <- data.frame(
            year = yr, strategy = strat,
            n_pixels_inside = n_mining,
            area_ha_inside  = n_mining * cell_ha)
          cat(sprintf("    %d  %-13s  %10.0f ha (inside polygons)\n",
                      yr, strat, n_mining * cell_ha))
        }
      }
      area_df <- do.call(rbind, rows)
      write.csv(area_df, area_csv, row.names = FALSE)
      cat(sprintf("  cached -> %s\n", area_csv))
    } else {
      area_df <- read.csv(area_csv)
      cat(sprintf("  loaded cache: %s (%d rows)\n", area_csv, nrow(area_df)))
    }
    
    # LBMA gold AM annual averages, USD/oz (rounded)
    gold_df <- data.frame(
      year        = 2000:2019,
      gold_usd_oz = c(279,  271,  310,  363,  410,
                      445,  604,  695,  872,  972,
                      1225, 1572, 1669, 1411, 1266,
                      1160, 1251, 1257, 1269, 1393)
    )
    gold_df <- gold_df[gold_df$year %in% area_df$year, ]
    sensor_df <- SENSOR_EVENTS[SENSOR_EVENTS$year >= min(area_df$year) &
                                 SENSOR_EVENTS$year <= max(area_df$year), ]
    
    ratio <- max(area_df$area_ha_inside, na.rm = TRUE) /
      max(gold_df$gold_usd_oz, na.rm = TRUE)
    
    strat_labels <- c(unconstrained = "Unconstrained (P1)",
                      shrinking     = "Transition-constrained (P2)")
    strat_cols   <- c(unconstrained = "#56B4E9", shrinking = "#D55E00")
    strat_lty    <- c(unconstrained = "22",      shrinking = "solid")
    
    fig_E <- ggplot() +
      # sensor changes (drawn first so the series sit on top)
      geom_vline(data = sensor_df, aes(xintercept = year),
                 colour = "grey55", linetype = "dotted", linewidth = 0.5) +
      geom_text(data = sensor_df, aes(x = year, y = Inf, label = label),
                angle = 90, hjust = 1.08, vjust = -0.5,
                size = E_BASE / 4.2, colour = "grey35") +
      geom_line(data = gold_df, aes(x = year, y = gold_usd_oz * ratio),
                colour = "#B8860B", linewidth = 0.9, linetype = "dotted") +
      geom_point(data = gold_df, aes(x = year, y = gold_usd_oz * ratio),
                 colour = "#B8860B", size = 2, shape = 18) +
      geom_line(data = area_df,
                aes(x = year, y = area_ha_inside, colour = strategy, linetype = strategy),
                linewidth = 1.1) +
      geom_point(data = area_df,
                 aes(x = year, y = area_ha_inside, colour = strategy), size = 2.2) +
      scale_y_continuous(
        name   = "Mining and water area inside polygons (ha)",
        labels = scales::label_comma(),
        # headroom at the top so the sensor labels do not cover the curves
        expand = expansion(mult = c(0.05, 0.24)),
        sec.axis = sec_axis(~ . / ratio, name = "LBMA gold price (USD / oz)",
                            labels = scales::label_comma())
      ) +
      scale_colour_manual(values = strat_cols, labels = strat_labels, name = NULL) +
      scale_linetype_manual(values = strat_lty, labels = strat_labels, name = NULL) +
      scale_x_continuous(breaks = seq(2000, 2020, 2)) +
      labs(x = NULL, title = NULL, subtitle = NULL, caption = NULL) +
      theme_pub(E_BASE) +
      theme(
        legend.position    = "bottom",
        legend.text        = element_text(size = E_BASE),
        legend.key.width   = unit(1.6, "cm"),
        axis.text          = element_text(size = E_BASE - 2),
        axis.title.y       = element_text(size = E_BASE),
        axis.title.y.right = element_text(colour = "#B8860B", size = E_BASE,
                                          margin = margin(l = 8)),
        axis.text.y.right  = element_text(colour = "#B8860B"),
        axis.line.y.right  = element_line(colour = "#B8860B"),
        axis.ticks.y.right = element_line(colour = "#B8860B")
      )
    
    ggsave(file.path(out_fig, "fig_constraint_vs_gold.png"),
           fig_E, width = 10, height = 5.6, dpi = 500, bg = "white")
    cat(sprintf("  -> %s\n", file.path(out_fig, "fig_constraint_vs_gold.png")))
  }
}

# ================================================================
# F. DETERMINISTIC VS PROBABILISTIC
# ================================================================
if (RUN[["F"]]) {
  cat("\n[F] Deterministic vs probabilistic (Mining + Water inside polygons)\n")
  
  F_BASE <- 15
  
  stoch_csv <- file.path(p3_dir, "statistics", "mining_area_in_polygons_stochastic.csv")
  if (!file.exists(stoch_csv)) {
    cat("  Skipped: run 30 first (", stoch_csv, ")\n", sep = "")
  } else {
    st <- read.csv(stoch_csv)
    p2_csv <- file.path(PATHS$out_constrained, "statistics", "area_within_mining_polygons.csv")
    if (file.exists(p2_csv)) {
      p2 <- read.csv(p2_csv)
      p2 <- p2[p2$stage == "shrinking" &
                 p2$class %in% c(CLASS_IDS[["Water"]], CLASS_IDS[["Mining"]]), ]
      det <- aggregate(area_ha_inside ~ year, data = p2, FUN = sum)
      st$deterministic_constrained_ha <- det$area_ha_inside[match(st$year, det$year)]
    }
    
    series <- list(
      list(key = "argmax", lab = "Ensemble draws (argmax of a random member per year)",
           col = "#0072B2", pre = "argmax_"),
      list(key = "sample", lab = "Probability-weighted (class drawn from probabilities)",
           col = "#8F7AB8", pre = ""))
    series <- Filter(function(s) paste0(s$pre, "mean_ha") %in% names(st), series)
    
    long <- do.call(rbind, lapply(series, function(s) data.frame(
      year = st$year, series = s$lab,
      mean = st[[paste0(s$pre, "mean_ha")]], sd = st[[paste0(s$pre, "sd_ha")]],
      lo95 = st[[paste0(s$pre, "pi95_lower_ha")]], hi95 = st[[paste0(s$pre, "pi95_upper_ha")]])))
    lab_det <- "Deterministic (argmax of ensemble mean, constrained)"
    cols <- c(setNames(vapply(series, `[[`, "", "col"), vapply(series, `[[`, "", "lab")),
              setNames("#D55E00", lab_det))
    
    fig_F <- ggplot() +
      geom_ribbon(data = long, aes(year, ymin = lo95, ymax = hi95, fill = series),
                  alpha = 0.18) +
      geom_ribbon(data = long, aes(year, ymin = pmax(mean - sd, 0), ymax = mean + sd,
                                   fill = series), alpha = 0.40) +
      geom_line(data = long, aes(year, mean, colour = series), linewidth = 1.1) +
      geom_point(data = long, aes(year, mean, colour = series), size = 2) +
      geom_line(data = st, aes(year, deterministic_constrained_ha, colour = lab_det),
                linewidth = 1.1, linetype = "22") +
      geom_point(data = st, aes(year, deterministic_constrained_ha, colour = lab_det),
                 shape = 17, size = 2.5) +
      scale_colour_manual(values = cols, breaks = names(cols), name = NULL) +
      scale_fill_manual(values = cols, guide = "none") +
      scale_y_continuous(name = "Mining and water area inside polygons (ha)",
                         labels = scales::label_comma(), limits = c(0, NA),
                         expand = expansion(mult = c(0, 0.04))) +
      scale_x_continuous(name = "Year", breaks = seq(2000, 2020, 2)) +
      guides(colour = guide_legend(ncol = 1,
                                   override.aes = list(linewidth = 1.6, size = 3.5))) +
      labs(title = NULL, subtitle = NULL, caption = NULL) +
      theme_pub(F_BASE) +
      theme(legend.position  = "bottom",
            legend.text      = element_text(size = F_BASE),
            legend.key.width = unit(1.8, "cm"),
            legend.key.height = unit(0.7, "cm"),
            axis.text        = element_text(size = F_BASE - 2))
    
    f_png <- file.path(out_fig, "fig_det_vs_prob.png")
    ggsave(f_png, fig_F, width = 10, height = 6.5, dpi = 500, bg = "white")
    ggsave(sub("\\.png$", ".pdf", f_png), fig_F, width = 10, height = 6.5, bg = "white",
           device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf)
    cat(sprintf("  -> %s (+ .pdf)\n", f_png))
  }
}

# ================================================================
# G. STUDY-AREA OVERVIEW MAP
# ----------------------------------------------------------------
# Terrain backdrop (muted elevation tint + light hillshade), districts,
# lakes, main rivers, mining polygons and major cities; legend on top.
# ================================================================
if (RUN[["G"]]) {
  cat("\n[G] Study-area overview map\n")
  
  G_BASE <- 15
  
  if (is.null(ghana_v) || is.null(mining_p)) {
    cat("  Skipped: districts or mining polygons missing.\n")
  } else {
    
    pacman::p_load(rnaturalearth)
    
    cat("  loading water bodies (Natural Earth)\n")
    fetch_ne <- function(type) tryCatch({
      x <- rnaturalearth::ne_download(scale = 10, type = type, category = "physical",
                                      returnclass = "sf")
      suppressWarnings(sf::st_intersection(sf::st_transform(x, GHANA_CRS), ghana_outline))
    }, error = function(e) { message("    ", type, ": ", e$message); NULL })
    lakes_g <- fetch_ne("lakes")
    # Main rivers: the global layer only carries the Volta system in Ghana;
    # the Natural Earth Africa supplement adds Pra, Ankobra, Tano, Offin,
    # Birim etc. Both are merged and drawn as one "Rivers" layer.
    rivers_g <- lapply(c("rivers_lake_centerlines", "rivers_africa"), function(tp) {
      x <- fetch_ne(tp)
      if (is.null(x) || !nrow(x)) return(NULL)
      x <- suppressWarnings(sf::st_collection_extract(x, "LINESTRING"))
      x["name"] <- if ("name" %in% names(x)) x$name else NA_character_
      x[, "name"]
    })
    rivers_g <- do.call(rbind, Filter(Negate(is.null), rivers_g))
    if (!is.null(rivers_g)) cat(sprintf("  rivers: %d segments (%s)\n", nrow(rivers_g),
                                        paste(head(sort(unique(na.omit(rivers_g$name))), 12),
                                              collapse = ", ")))
    
    cities_g <- cities_sf[order(-cities_sf$pop), ][1:8, ]
    
    # ---- Terrain backdrop: muted elevation tint + light hillshade ----
    # DEM: SRTM 30 arc-sec (~1 km) via geodata (cached in viz_cache);
    # fallback: ELEVATION band of the 2019 training composite.
    G_DEM_RES_M  <- 400     # plotting resolution (m)
    G_HS_EXAGG   <- 4       # vertical exaggeration for the hillshade (Ghana is flat)
    G_HS_ALPHA   <- 0.35    # hillshade strength
    G_TINTS      <- c("#EEF2E6", "#DDE6D3", "#CCD6C1", "#BFC5B3", "#B0AFA2")  # grey-green, no orange
    
    dem <- tryCatch({
      pacman::p_load(geodata)
      geodata::elevation_30s(country = "GHA", path = cache_dir)
    }, error = function(e) { message("    geodata DEM: ", e$message); NULL })
    if (is.null(dem) && file.exists(NATIONAL_CONFIG$train_composite)) {
      cat("  DEM fallback: ELEVATION band of the training composite\n")
      dem <- rast(NATIONAL_CONFIG$train_composite, lyrs = which(BAND_NAMES == "ELEVATION"))
      dem <- aggregate(dem, fact = max(1, round(G_DEM_RES_M / res(dem)[1])),
                       fun = "mean", na.rm = TRUE)
    }
    hs <- NULL
    if (!is.null(dem)) {
      tmpl <- rast(ext(vect(ghana_outline)), resolution = G_DEM_RES_M, crs = GHANA_CRS)
      dem  <- project(dem, tmpl, method = "bilinear")
      dem  <- mask(dem, vect(ghana_outline))
      names(dem) <- "elev"
      tr <- terrain(dem * G_HS_EXAGG, c("slope", "aspect"), unit = "radians")
      hs <- shade(tr$slope, tr$aspect, angle = 40, direction = 315)
      hs <- mask(hs, vect(ghana_outline)); names(hs) <- "hs"
      cat(sprintf("  DEM: %.0f-%.0f m\n", global(dem, "min", na.rm = TRUE)[1, 1],
                  global(dem, "max", na.rm = TRUE)[1, 1]))
    } else cat("  no DEM available - plain background\n")
    
    add_rast <- function(g, r, ...) {
      if (HAS_TIDYTERRA) return(g + tidyterra::geom_spatraster(data = r, ...))
      d <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(d)[3] <- "v"
      g + geom_raster(data = d, aes(x, y, fill = v), ...)
    }
    
    fill_vals <- c("Mining polygons" = scales::alpha("#D55E00", 0.75),
                   "Lakes"           = scales::alpha("#A6CEE3", 0.85))
    fig_G <- ggplot() +
      geom_sf(data = ghana_outline, fill = "grey97", colour = NA)
    if (!is.null(dem)) {
      fig_G <- add_rast(fig_G, dem) +
        scale_fill_gradientn(
          colours = G_TINTS, na.value = NA, name = "Elevation (m)",
          guide = guide_colourbar(order = 3, direction = "horizontal",
                                  title.position = "top", title.hjust = 0.5,
                                  barwidth = unit(5, "cm"), barheight = unit(0.35, "cm"))) +
        ggnewscale::new_scale_fill()
      fig_G <- add_rast(fig_G, hs, alpha = G_HS_ALPHA) +
        scale_fill_gradient(low = "grey15", high = "white", na.value = NA, guide = "none") +
        ggnewscale::new_scale_fill()
    }
    fig_G <- fig_G +
      geom_sf(data = ghana_v, fill = NA, colour = scales::alpha("white", 0.7), linewidth = 0.2)
    
    if (!is.null(lakes_g) && nrow(lakes_g) > 0)
      fig_G <- fig_G + geom_sf(data = lakes_g, aes(fill = "Lakes"),
                               colour = "#4682B4", linewidth = 0.3)
    if (!is.null(rivers_g) && nrow(rivers_g) > 0)
      fig_G <- fig_G + geom_sf(data = rivers_g, aes(colour = "Rivers"), linewidth = 0.4)
    
    fig_G <- fig_G +
      geom_sf(data = mining_p, aes(fill = "Mining polygons"), colour = "#A63603",
              linewidth = 0.3) +
      scale_fill_manual(values = fill_vals, breaks = names(fill_vals), name = NULL,
                        guide = guide_legend(order = 1, override.aes =
                                               list(colour = c("#A63603", "#4682B4")))) +
      scale_colour_manual(values = c(Rivers = "#4682B4"), name = NULL,
                          guide = guide_legend(order = 2,
                                               override.aes = list(linewidth = 1.2))) +
      geom_sf(data = ghana_outline, fill = NA, colour = "grey20", linewidth = 0.5) +
      geom_sf(data = cities_g, shape = 21, fill = "white",
              colour = "black", size = 3, stroke = 0.7)
    
    if (HAS_GGREPEL) {
      fig_G <- fig_G +
        ggrepel::geom_label_repel(
          data = cities_g,
          aes(label = name, geometry = geometry),
          stat = "sf_coordinates",
          size = 4.6, fontface = "bold",
          fill = scales::alpha("white", 0.85), label.size = 0,
          box.padding = 0.5, point.padding = 0.3,
          min.segment.length = 0.1,
          segment.colour = "black", segment.size = 0.3,
          seed = 42
        )
    } else {
      fig_G <- fig_G +
        geom_sf_text(data = cities_g, aes(label = name),
                     size = 4.6, fontface = "bold",
                     nudge_x = 8000, nudge_y = 8000)
    }
    
    fig_G <- fig_G +
      ggspatial::annotation_scale(
        location = "br", width_hint = 0.20,
        pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm"),
        height = unit(0.3, "cm"), style = "bar", text_cex = 1.2
      ) +
      ggspatial::annotation_north_arrow(
        location = "tr", which_north = "true",
        pad_x = unit(0.4, "cm"), pad_y = unit(0.4, "cm"),
        style = ggspatial::north_arrow_minimal(line_width = 1.2, text_size = 11),
        height = unit(1.4, "cm"), width = unit(1.4, "cm")
      ) +
      coord_sf(crs = GHANA_CRS, expand = FALSE) +
      theme_map(G_BASE) +
      theme(legend.position  = "top",
            legend.box       = "horizontal",
            legend.title     = element_text(size = G_BASE),
            legend.text      = element_text(size = G_BASE),
            legend.key.size  = unit(0.9, "cm"),
            legend.key.width = unit(1.0, "cm"),
            legend.spacing.x = unit(0.4, "cm"),
            legend.margin    = margin(b = 6))
    
    ggsave(file.path(out_fig, "fig_study_area.png"),
           fig_G, width = 10, height = 11.5, dpi = 500, bg = "white")
    cat(sprintf("  -> %s\n", file.path(out_fig, "fig_study_area.png")))
  }
}

# ================================================================
# H. INDUSTRIAL MINES OVER TIME
# ----------------------------------------------------------------
# 4 large industrial open-pit mines in 2000 / 2010 / 2019, constrained
# map with the A/H colour scheme; polygon outlines in pink. Year and
# [PRE]/[OPS] label below each map; no legend (set in LaTeX).
# ================================================================
if (RUN[["H"]]) {
  cat("\n[H] Industrial mines appendix figure\n")
  
  INDUSTRIAL_SITES <- list(
    tarkwa = list(lon = -1.9933, lat = 5.3000, start_year = 1998,
                  label = "Tarkwa Gold Mine\n(Gold Fields, 1998)"),
    ahafo  = list(lon = -2.3500, lat = 7.0500, start_year = 2006,
                  label = "Ahafo (Kenyasi) Gold Mine\n(Newmont, 2006)"),
    akyem  = list(lon = -1.0000, lat = 6.3000, start_year = 2013,
                  label = "Akyem Gold Mine\n(Newmont / Zijin, 2013)"),
    damang = list(lon = -1.8500, lat = 5.5100, start_year = 1997,
                  label = "Damang Gold Mine\n(Gold Fields, 1997)")
  )
  YEARS_SHOW   <- c(2000, 2010, 2019)
  SITE_HALF_M  <- 12000    # 24 km panels
  ROW_LABEL_SZ <- 4.4
  H_LABEL_SZ   <- 10       # [PRE]/[OPS] label below each map
  POLY_COL     <- "#FF1493"
  
  build_site_panel <- function(site, year) {
    if (!file.exists(pred_path_for(year))) return(patchwork::plot_spacer())
    xy  <- to_utm_xy(site$lon, site$lat)
    bx  <- ext(xy[1] - SITE_HALF_M, xy[1] + SITE_HALF_M,
               xy[2] - SITE_HALF_M, xy[2] + SITE_HALF_M)
    bsf <- bbox_sf(bx)
    
    pred   <- crop(rast(pred_path_for(year)), bx)
    study  <- align_to(valid_mask(read_rgb(2019, bx)), pred)   # composite footprint
    valid  <- align_to(valid_mask(read_rgb(year, bx)), pred)   # data in this year
    inside <- inside_mask_for(pred, bsf) == 1
    cls    <- classify_cells(pred, inside, valid, study, veg_code = 5L)
    
    fr   <- terra::freq(cls)
    n_of <- function(v) sum(fr$count[fr$value %in% v])
    pct  <- 100 * n_of(1:4) / max(1, n_of(1:5))    # mining + water / valid cells
    
    df <- as.data.frame(cls, xy = TRUE, na.rm = TRUE)
    names(df)[3] <- "class"
    df$class <- factor(CLS_LEVELS[df$class], levels = CLS_LEVELS)
    maus <- suppressWarnings(sf::st_intersection(mining_p, bsf))
    
    is_ops  <- year >= site$start_year
    tag_col <- if (is_ops) "#B22222" else "#228B22"
    tag_lbl <- if (is_ops) "[OPS]"   else "[PRE]"
    
    g <- ggplot() + geom_raster(data = df, aes(x, y, fill = class))
    if (nrow(maus))
      g <- g + geom_sf(data = maus, colour = POLY_COL, fill = NA, linewidth = 0.45)
    g + scale_fill_manual(values = CLS_COLS, drop = FALSE, guide = "none") +
      map_furniture(dark = FALSE, size = 0.65, width_hint = 0.3) +
      grat_scales(c(xmin(bx), xmax(bx)), c(ymin(bx), ymax(bx)), n = 3) +
      coord_sf(crs = GHANA_CRS, datum = sf::st_crs(4326),
               xlim = c(xmin(bx), xmax(bx)), ylim = c(ymin(bx), ymax(bx)),
               expand = FALSE) +
      labs(caption = sprintf("%d %s  Mining + water: %.1f%%", year, tag_lbl, pct)) +
      theme_map_grid(base = 8, border_col = tag_col) +
      theme(
        plot.caption = element_text(colour = tag_col, face = "bold",
                                    size = H_LABEL_SZ, hjust = 0.5,
                                    margin = margin(t = 4)),
        plot.caption.position = "panel",
        panel.border = element_rect(fill = NA, colour = tag_col, linewidth = 0.8),
        axis.text    = element_text(size = 6, colour = "grey25"),
        legend.position = "none",
        plot.margin  = margin(2, 3, 4, 2)
      )
  }
  
  make_row_label <- function(text) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = text, angle = 90,
               size = ROW_LABEL_SZ, fontface = "bold", lineheight = 0.95,
               hjust = 0.5, vjust = 0.5) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      theme_void() +
      theme(plot.margin      = margin(2, 0, 2, 2),
            panel.background = element_rect(fill = "white", colour = NA))
  }
  
  fig_H_rows <- lapply(names(INDUSTRIAL_SITES), function(s) {
    ps <- lapply(YEARS_SHOW, function(yr) {
      cat(sprintf("  building %s_%d\n", s, yr))
      tryCatch(build_site_panel(INDUSTRIAL_SITES[[s]], yr), error = function(e) {
        message("    failed: ", e$message); patchwork::plot_spacer() })
    })
    # widths are set per row (applying them to the outer stack has no
    # effect on nested rows)
    (make_row_label(INDUSTRIAL_SITES[[s]]$label) | ps[[1]] | ps[[2]] | ps[[3]]) +
      plot_layout(widths = c(0.13, 1, 1, 1))
  })
  
  fig_H <- patchwork::wrap_plots(fig_H_rows, ncol = 1)
  ggsave(file.path(out_fig, "fig_industrial_mines_appendix.png"),
         fig_H, width = 10, height = 3.5 * length(fig_H_rows), dpi = 500, bg = "white")
  cat(sprintf("  -> %s\n", file.path(out_fig, "fig_industrial_mines_appendix.png")))
}

# ================================================================
# I. VALID PIXEL COUNT PER COMPOSITE YEAR
# ----------------------------------------------------------------
# Number of valid (cloud-free) Landsat observations behind each annual
# composite, sampled inside the mining polygons + 5 km buffer, per
# latitude zone: mean and 5th-95th percentile band. Thresholds as in
# the composite QC (42): good >= 15, moderate >= 8; dotted line = 5.
# No title (caption in LaTeX). The sample table is cached and rebuilt
# when a composite changes.
# ================================================================
if (RUN[["I"]]) {
  cat("\n[I] Valid pixel count per composite year\n")
  
  I_BASE     <- 15
  I_N_SAMPLE <- 15000
  I_BUFFER_M <- 5000
  
  vpc_csv   <- file.path(cache_dir, "valid_pixel_count_by_zone_year.csv")
  vpc_years <- PREDICT_CONFIG$predict_years
  vpc_paths <- sprintf(NATIONAL_CONFIG$year_composite_fmt, vpc_years)
  keep      <- file.exists(vpc_paths)
  vpc_years <- vpc_years[keep]; vpc_paths <- vpc_paths[keep]
  
  if (is.null(mining_p) || !length(vpc_paths)) {
    cat("  Skipped: mining polygons or composites missing.\n")
  } else {
    if (!cache_fresh(vpc_csv, vpc_paths)) {
      cat(sprintf("  sampling %d points in polygons + %.0f km buffer\n",
                  I_N_SAMPLE, I_BUFFER_M / 1000))
      set.seed(CV_CONFIG$seed)
      study <- sf::st_buffer(sf::st_union(sf::st_make_valid(mining_p)), I_BUFFER_M)
      pts   <- sf::st_sample(study, size = I_N_SAMPLE)
      lat   <- sf::st_coordinates(sf::st_transform(pts, 4326))[, 2]
      zone  <- assign_zone_from_latitude(lat)
      xy    <- sf::st_coordinates(pts)
      vpc_lyr <- which(BAND_NAMES == "VALID_PIXEL_COUNT")
      
      vpc_df <- dplyr::bind_rows(lapply(seq_along(vpc_paths), function(i) {
        r <- rast(vpc_paths[i], lyrs = vpc_lyr)
        v <- terra::extract(r, xy)[, 1]
        cat(sprintf("    %d\n", vpc_years[i]))
        data.frame(year = vpc_years[i], zone = zone, vpc = v) |>
          dplyr::filter(!is.na(vpc)) |>
          dplyr::group_by(year, zone) |>
          dplyr::summarise(mean = mean(vpc),
                           p05  = unname(quantile(vpc, 0.05)),
                           p95  = unname(quantile(vpc, 0.95)),
                           n    = dplyr::n(), .groups = "drop")
      }))
      write.csv(vpc_df, vpc_csv, row.names = FALSE)
      cat(sprintf("  cached -> %s\n", vpc_csv))
    } else {
      vpc_df <- read.csv(vpc_csv)
      cat(sprintf("  loaded cache: %s\n", vpc_csv))
    }
    
    vpc_df$zone <- factor(vpc_df$zone, levels = ZONE_ORDER,
                          labels = unname(ZONE_LABELS[ZONE_ORDER]))
    zone_cols <- setNames(c("#C46A1E", "#2E8B2E", "#4682B4"),
                          unname(ZONE_LABELS[ZONE_ORDER]))
    thr <- data.frame(y   = c(15, 8, 5),
                      lab = c("good (\u226515)", "moderate (\u22658)", "poor (<5)"),
                      lty = c("solid", "dashed", "dotted"))
    
    fig_I <- ggplot(vpc_df, aes(year, mean, colour = zone, fill = zone)) +
      geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.15, colour = NA) +
      geom_hline(data = thr, aes(yintercept = y, linetype = lty),
                 colour = "grey40", linewidth = 0.6, inherit.aes = FALSE) +
      geom_text(data = thr, aes(x = min(vpc_df$year), y = y, label = lab),
                inherit.aes = FALSE, hjust = 0, vjust = -0.4,
                size = I_BASE / 3.6, colour = "grey30") +
      geom_line(linewidth = 1.1) +
      geom_point(size = 2.4) +
      scale_linetype_identity() +
      scale_colour_manual(values = zone_cols, name = "Zone") +
      scale_fill_manual(values = zone_cols, name = "Zone") +
      scale_x_continuous(breaks = seq(2000, 2020, 2)) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.04))) +
      labs(x = "Year", y = "Valid pixel count",
           title = NULL, subtitle = NULL, caption = NULL) +
      theme_pub(I_BASE) +
      theme(legend.position  = "bottom",
            legend.title     = element_text(size = I_BASE),
            legend.text      = element_text(size = I_BASE),
            legend.key.width = unit(1.4, "cm"),
            axis.text        = element_text(size = I_BASE - 2))
    
    i_png <- file.path(out_fig, "fig_valid_pixel_count.png")
    ggsave(i_png, fig_I, width = 10, height = 5.5, dpi = 500, bg = "white")
    ggsave(sub("\\.png$", ".pdf", i_png), fig_I, width = 10, height = 5.5, bg = "white",
           device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf)
    cat(sprintf("  -> %s (+ .pdf)\n", i_png))
    
    cat("\n  Mean valid pixel count, selected years (for the appendix table):\n")
    print(tidyr::pivot_wider(
      vpc_df[vpc_df$year %in% c(2000, 2001, 2005, 2010, 2013, 2015, 2019),
             c("year", "zone", "mean")],
      names_from = zone, values_from = mean), digits = 3)
  }
}

cat("\nDone.\n")