# ================================================================
# 45a_check_anchors.R
# Do the industrial anchors sit on mining sites and inside the
# company's official Mining Lease? Run after 00_config.R and 01_utils.R.
#
# Output (PATHS$out_validation/cnn_rf_consistency/anchor_check):
#   anchor_check_panels.png  (-> P5_RESULTS/figures/45_cnn_rf_onset) one panel per anchor, 2019 false colour
#                            (SWIR1-NIR-RED: bare/mined = pink-magenta,
#                            vegetation = green, water = dark blue)
#                            yellow  = mining polygons (Maus)
#                            orange  = company's Mining Lease(s) (cadastre)
#                            white dashed = 3 / 5 km circles
#                            red fill = 2019 Mining+Water cells inside the
#                                       5 km footprint (the mine's max extent)
#                            cyan +  = anchor, red x = largest polygon nearby
#   anchor_check.csv         polygon diagnostics + lease match per anchor
#   anchor_mined_extent_2019.csv
#                            per anchor and footprint (5 km clip, lease clip):
#                            polygon footprint ha, 2019 Mining+Water ha (max
#                            extent; vegetation in the polygons NOT counted),
#                            vegetation ha, mined share, T10 threshold ha
# ================================================================

pacman::p_load(terra, sf, dplyr, tibble, readr, xml2)
sf::sf_use_s2(FALSE)

out <- file.path(PATHS$out_validation, "cnn_rf_consistency", "anchor_check")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

LEASE_KML <- PATHS$lease_kml

# Same settings as script 45
EXCLUDE_ANCHORS <- c("Obuasi", "Prestea", "Awaso")  # not used for CNN training / statistics
SNAP_TO_LEASE   <- c("Awaso")              # moved onto mined polygon inside own lease
SNAP_MAX_M      <- 5000
CLIP_M          <- 5000    # main footprint: polygons within 5 km (Voronoi-split)
LEASE_MAX_M     <- 10000   # lease clip: polygons in own lease within 10 km
T10_FRAC        <- 0.10    # establishment threshold = 10% of 2019 mined extent
T10_MIN_HA      <- 10

PRED_2019 <- file.path(PATHS$out_constrained, "predictions", "prediction_2019_shrinking.tif")
if (!file.exists(PRED_2019))
  PRED_2019 <- file.path(PATHS$national_predictions, "prediction_2019.tif")

# Verified anchors (identical to ANCHORS_VERIFIED in script 45)
ANCHORS <- tribble(
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

# anchor -> regex on "owner | code" (same as CFG$lease_map in script 45)
LEASE_MAP <- c(
  Iduapriem = "Iduapriem|ML2/187", Tarkwa = "Gold Fields Ghana", Damang = "Abosso",
  Ahafo = "Newmont Ghana Gold|RL7/2", Akyem = "Golden Ridge", Chirano = "Chirano",
  Bibiani = "Bibiani", Bogoso = "Blue Gold Bogoso",
  Prestea = "Blue Gold Bogoso|Prestea Sankofa", Asanko = "Asanko",
  Wassa = "Star \\(Wassa\\)", Edikan = "Perseus", Nzema = "Adamus",
  Namdini = "Cardinal Namdini", Nsuta = "Ghana Manganese", Awaso = "Ghana Bauxite",
  Obuasi = "Anglogold Ashanti \\(Ghana\\)")   # shown here; not used as clip in 45

# ---- Data ----------------------------------------------------------
comp <- rast(file.path(PATHS$composite_dir, "Mining_Composite_2019_full.tif"))
names(comp) <- BAND_NAMES
fc <- comp[[c("SWIR1", "NIR", "RED")]]

anc <- st_as_sf(ANCHORS, coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  st_transform(crs(comp))
polys <- st_read(PATHS$mining_polygons, quiet = TRUE) |>
  st_transform(crs(comp)) |> st_make_valid()
polys$area_ha <- as.numeric(st_area(polys)) / 1e4

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
if (file.exists(LEASE_KML)) {
  leases <- read_licences(LEASE_KML) |> st_transform(crs(comp)) |> st_make_valid()
  leases <- leases[grepl("Mining Lease", leases$type) & !grepl("Restricted", leases$type), ]
  leases$key <- paste(leases$owner, "|", leases$code)
  cat(nrow(leases), "large-scale Mining Lease polygons\n")
} else cat("Lease KML not found:", LEASE_KML, "\n")

# ---- Snap selected anchors onto a mined polygon inside their lease -----
anc$moved_m <- 0
if (!is.null(leases)) for (nm in SNAP_TO_LEASE) {
  i   <- which(anc$name == nm); pat <- LEASE_MAP[nm]
  own  <- leases[grepl(pat, leases$key), ]
  near <- polys[lengths(st_is_within_distance(polys, anc[i, ], SNAP_MAX_M)) > 0, ]
  if (!nrow(own) || !nrow(near)) next
  pcs <- suppressWarnings(st_intersection(st_geometry(near), st_union(st_geometry(own))))
  pcs <- suppressWarnings(st_collection_extract(pcs, "POLYGON"))
  pcs <- pcs[!st_is_empty(pcs)]
  if (!length(pcs)) next
  new <- st_point_on_surface(pcs[which.max(as.numeric(st_area(pcs)))])
  ll  <- st_coordinates(st_transform(new, 4326))
  anc$moved_m[i] <- round(as.numeric(st_distance(anc[i, ], new)))
  st_geometry(anc)[i] <- new
  anc$lat[i] <- ll[2]; anc$lon[i] <- ll[1]
  cat(sprintf("snapped %s: lat %.5f, lon %.5f (moved %d m)\n", nm, ll[2], ll[1], anc$moved_m[i]))
}
anc$cnn_role <- ifelse(anc$name %in% EXCLUDE_ANCHORS, "excluded",
                       ifelse(anc$established > 2019, "negative_control", "positive"))

# ---- 2019 mined extent per footprint (the T10 reference) ----------
# Footprint = mining polygons ∩ clip ∩ Voronoi cell of the anchor (as in 45).
# Max extent = RF 2019 cells classed Mining or Water inside the footprint;
# Vegetation cells inside the polygons are reported but NOT counted.
safe_inter <- function(a, b) {
  x <- suppressWarnings(st_intersection(a, b))
  if (!length(x)) return(st_sfc(st_polygon(), crs = st_crs(a)))
  st_union(x)
}
envl <- st_as_sfc(st_bbox(st_buffer(st_union(anc), 300000)))
vor  <- st_collection_extract(st_voronoi(st_union(st_geometry(anc)), envelope = envl), "POLYGON")
vor_id <- vapply(st_intersects(vor, anc), function(i) i[1], integer(1))

pred19 <- if (file.exists(PRED_2019)) rast(PRED_2019) else NULL
if (is.null(pred19)) cat("2019 prediction not found — mined extent skipped\n")

footprint <- function(i, type = c("clip", "lease")) {
  type <- match.arg(type)
  pt   <- st_geometry(anc)[i]
  cell <- vor[vor_id == i]
  clip <- if (type == "clip") safe_inter(st_buffer(pt, CLIP_M), cell) else {
    pat <- LEASE_MAP[anc$name[i]]
    own <- if (is.null(leases) || is.na(pat)) NULL else leases[grepl(pat, leases$key), ]
    if (is.null(own) || !nrow(own)) return(NULL)
    safe_inter(safe_inter(st_buffer(pt, LEASE_MAX_M), cell), st_union(st_geometry(own)))
  }
  hit <- lengths(st_intersects(polys, clip)) > 0
  if (!any(hit)) return(NULL)
  fp <- safe_inter(clip, st_union(st_geometry(polys)[hit]))
  if (!length(fp) || as.numeric(st_area(fp)) == 0) NULL else fp
}

extent_19 <- function(fp) {
  if (is.null(fp) || is.null(pred19)) return(c(fp_ha = NA, mined_ha = NA, veg_ha = NA))
  fpt <- st_transform(st_sf(geometry = fp), crs(pred19))
  e   <- terra::intersect(ext(vect(st_buffer(fpt, 100))), ext(pred19))
  if (is.null(e)) return(c(fp_ha = NA, mined_ha = NA, veg_ha = NA))   # outside RF extent
  r   <- crop(pred19, e)
  ins <- terra::rasterize(vect(fpt), r, field = 1, background = NA)
  ch  <- prod(res(r)) / 1e4
  cnt <- function(x) { v <- as.numeric(global(x, "sum", na.rm = TRUE)[1, 1]); ifelse(is.finite(v), v, 0) }
  c(fp_ha    = cnt(ins) * ch,
    mined_ha = cnt((r == CLASS_IDS[["Mining"]] | r == CLASS_IDS[["Water"]]) * ins) * ch,
    veg_ha   = cnt((r == CLASS_IDS[["Vegetation"]]) * ins) * ch)
}

fp_clip <- lapply(seq_len(nrow(anc)), footprint, type = "clip")
mined_tab <- lapply(seq_len(nrow(anc)), function(i) {
  bind_rows(lapply(c("clip", "lease"), function(tp) {
    fp <- if (tp == "clip") fp_clip[[i]] else footprint(i, "lease")
    e  <- extent_19(fp)
    tibble(name = anc$name[i], cnn_role = anc$cnn_role[i],
           footprint = ifelse(tp == "clip", sprintf("%g km clip", CLIP_M / 1000), "lease clip"),
           fp_ha = round(e[["fp_ha"]]), mined2019_ha = round(e[["mined_ha"]]),
           veg2019_ha = round(e[["veg_ha"]]),
           mined_share = round(e[["mined_ha"]] / e[["fp_ha"]], 2),
           T10_threshold_ha = round(pmax(T10_MIN_HA, T10_FRAC * e[["mined_ha"]])))
  }))
}) |> bind_rows()
write_csv(mined_tab, file.path(out, "anchor_mined_extent_2019.csv"))
cat("\n2019 mined extent (Mining+Water in footprint; vegetation excluded):\n")
print(as.data.frame(mined_tab))

# ---- Diagnostics per anchor --------------------------------------
diag <- lapply(seq_len(nrow(anc)), function(i) {
  pt  <- anc[i, ]
  own <- if (!is.null(leases) && !is.na(LEASE_MAP[pt$name]))
    leases[grepl(LEASE_MAP[pt$name], leases$key), ] else NULL
  lease_cols <- tibble(
    n_own_leases = if (is.null(own)) 0L else nrow(own),
    in_own_lease = !is.null(own) && nrow(own) > 0 &&
      any(lengths(st_intersects(pt, own)) > 0),
    dist_own_lease_m = if (is.null(own) || !nrow(own)) NA_real_ else
      round(min(as.numeric(st_distance(pt, own)))),
    own_lease_codes = if (is.null(own)) "" else paste(own$code, collapse = "; "))
  
  p10 <- polys[lengths(st_is_within_distance(polys, pt, 10000)) > 0, ]
  if (!nrow(p10))
    return(bind_cols(tibble(name = pt$name, cnn_role = pt$cnn_role, lat = pt$lat,
                            lon = pt$lon, moved_m = pt$moved_m,
                            inside_polygon = FALSE, dist_nearest_m = NA,
                            ha_1km = 0, ha_3km = 0, ha_5km = 0, largest_poly_ha = NA,
                            largest_dist_m = NA, sugg_lat = NA, sugg_lon = NA), lease_cols))
  d   <- as.numeric(st_distance(pt, p10))
  big <- p10[which.max(p10$area_ha), ]
  sg  <- st_coordinates(st_transform(st_point_on_surface(st_geometry(big)), 4326))
  bind_cols(tibble(
    name = pt$name, cnn_role = pt$cnn_role, lat = pt$lat, lon = pt$lon,
    moved_m = pt$moved_m,
    inside_polygon  = any(d == 0),
    dist_nearest_m  = round(min(d)),
    # NOTE: full area of polygons TOUCHING the radius, not area inside it
    ha_1km = round(sum(p10$area_ha[d <= 1000])),
    ha_3km = round(sum(p10$area_ha[d <= 3000])),
    ha_5km = round(sum(p10$area_ha[d <= 5000])),
    largest_poly_ha = round(big$area_ha),
    largest_dist_m  = round(as.numeric(st_distance(pt, big))),
    sugg_lat = round(sg[2], 5), sugg_lon = round(sg[1], 5)), lease_cols)
}) |> bind_rows() |>
  mutate(flag = case_when(is.na(dist_nearest_m) ~ "NO POLYGON within 10 km",
                          dist_nearest_m > 1000 | ha_3km < 50 ~ "CHECK",
                          TRUE ~ "ok"),
         lease_flag = case_when(is.null(leases)            ~ "no lease data",
                                in_own_lease               ~ "in lease",
                                !is.na(dist_own_lease_m) &
                                  dist_own_lease_m <= 1000 ~ "near lease (<1 km)",
                                n_own_leases == 0          ~ "no lease in export",
                                TRUE                       ~ "outside lease"))
write_csv(diag, file.path(out, "anchor_check.csv"))
print(as.data.frame(diag |> select(name, cnn_role, moved_m, flag, lease_flag, dist_nearest_m, ha_3km,
                                   dist_own_lease_m, own_lease_codes)))

# ---- Panels ------------------------------------------------------
graphics.off()
png_file <- file.path(results_fig_dir("45_cnn_rf_onset"), "anchor_check_panels.png")
png(png_file, width = 4000, height = 3300, res = 200)
par(mfrow = c(4, 5), mar = c(0.5, 0.5, 2, 0.5))
for (i in seq_len(nrow(anc))) {
  pt  <- anc[i, ]
  win <- st_buffer(st_geometry(pt), 6000)
  r   <- crop(fc, ext(vect(win)))
  m_i <- mined_tab$mined2019_ha[mined_tab$name == pt$name & mined_tab$footprint != "lease clip"]
  ttl <- sprintf("%s (%d) — %s / %s — mined 2019: %s ha%s", pt$name, pt$established,
                 diag$flag[i], diag$lease_flag[i], ifelse(is.na(m_i), "NA", m_i),
                 ifelse(pt$cnn_role == "excluded", " [CNN: excluded]",
                        ifelse(pt$moved_m > 0, " [snapped]", "")))
  if (all(is.na(values(r[[1]])))) { plot.new(); title(paste(ttl, "\nno composite data")); next }
  plotRGB(r, r = 1, g = 2, b = 3, stretch = "lin", axes = FALSE, main = ttl,
          mar = c(0.5, 0.5, 2, 0.5), cex.main = 0.85)
  if (!is.null(pred19) && !is.null(fp_clip[[i]])) {
    fpt <- st_transform(st_sf(geometry = fp_clip[[i]]), crs(pred19))
    e19 <- terra::intersect(ext(vect(st_transform(st_sf(geometry = win), crs(pred19)))),
                            ext(pred19))
    m19 <- if (is.null(e19)) NULL else
      mask(ifel(crop(pred19, e19) == CLASS_IDS[["Mining"]] |
                  crop(pred19, e19) == CLASS_IDS[["Water"]], 1, NA), vect(fpt))
    if (!is.null(m19) && !all(is.na(values(m19))))
      plot(project(m19, crs(fc), method = "near"), add = TRUE, legend = FALSE,
           col = adjustcolor("red", alpha.f = 0.45))
  }
  pw <- polys[lengths(st_intersects(polys, win)) > 0, ]
  if (nrow(pw)) plot(vect(pw), add = TRUE, border = "yellow", lwd = 0.7)
  if (!is.null(leases) && !is.na(LEASE_MAP[pt$name])) {
    lw <- leases[grepl(LEASE_MAP[pt$name], leases$key) &
                   lengths(st_intersects(leases, win)) > 0, ]
    if (nrow(lw)) plot(vect(lw), add = TRUE, border = "orange", lwd = 2.2)
  }
  plot(vect(st_buffer(st_geometry(pt)[c(1, 1)], c(3000, 5000))), add = TRUE,
       border = "white", lty = 2)
  plot(vect(pt), add = TRUE, pch = 3, cex = 2.5, lwd = 3, col = "cyan")
  if (!is.na(diag$sugg_lat[i])) {
    sp <- st_transform(st_sfc(st_point(c(diag$sugg_lon[i], diag$sugg_lat[i])), crs = 4326),
                       crs(comp))
    plot(vect(sp), add = TRUE, pch = 4, cex = 2, lwd = 3, col = "red")
  }
}
dev.off()
cat("Saved:", png_file, "\n")