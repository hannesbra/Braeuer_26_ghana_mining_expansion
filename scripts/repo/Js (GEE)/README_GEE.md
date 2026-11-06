# Google Earth Engine scripts

These scripts build the annual Landsat composites and export the manually
collected training samples used in:

> Bräuer, H. (???). *Quantification of Ghana's historical mining expansion.*
> GIScience & Remote Sensing. https://doi.org/TODO

The Earth Engine part of the workflow produces two inputs for the R analysis:

1. Annual 3-year dry-season medoid composites, 2000–2019 (GeoTIFF)
2. Manually labelled training points (CSV)

| Script | Purpose | Output |
|---|---|---|
| `01_landsat_composites.js` | Builds one composite per year, 2000–2019, with spectral indices, texture and terrain bands | 20 GeoTIFFs in Google Drive folder `GEE_Ghana_Mining_Medoid` |
| `02_export_training_samples.js` | Displays reference layers for manual labelling and exports the labelled points | CSVs in Google Drive folder `GEE_Manual_Sampling`, plus Earth Engine table assets |

TODO: rename the files above if your script names differ.

---

## Requirements

- A Google Earth Engine account (free for non-commercial and research use):
  https://earthengine.google.com
- The scripts run in the Earth Engine **JavaScript Code Editor**
  (https://code.earthengine.google.com). No local installation is needed.
- About 5–10 GB of free Google Drive space for the composite exports.

### Public datasets (loaded automatically)

| Dataset | Earth Engine ID | Use |
|---|---|---|
| Landsat 5 TM C2 L2 | `LANDSAT/LT05/C02/T1_L2` | Surface reflectance |
| Landsat 7 ETM+ C2 L2 | `LANDSAT/LE07/C02/T1_L2` | Surface reflectance |
| Landsat 8 OLI C2 L2 | `LANDSAT/LC08/C02/T1_L2` | Surface reflectance |
| Landsat 9 OLI-2 C2 L2 | `LANDSAT/LC09/C02/T1_L2` | Surface reflectance (only for years after 2021; unused in the 2000–2019 run) |
| SRTM 1 arc-second | `USGS/SRTMGL1_003` | Elevation, slope, aspect, ruggedness |
| LSIB 2017 | `USDOS/LSIB_SIMPLE/2017` | Ghana boundary |
| ESA WorldCover 2020 v100 | `ESA/WorldCover/v100` | Visual guidance during labelling only |
| VIIRS DNB / DMSP-OLS | `NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG`, `NOAA/DMSP-OLS/NIGHTTIME_LIGHTS` | Night-light bands (exported but **not used** in the classifier) |

### Your own assets (must be uploaded first)

Both scripts read assets from `projects/ee-hannes23braeuer/assets/`. To
reproduce the workflow, upload these files from the data archive to your own
Earth Engine project. Then change the asset paths at the top of each script.

| Asset name in scripts | Source file | Used by |
|---|---|---|
| `ghana_polygons2019` | `data/mining_polygons_ghana_2019.gpkg` (Ghana subset of Maus et al., 2022) | both |
| `Mining_Composite_2019` | the 2019 output of script 01 | 02 |

---

## 01 – Landsat composites

### Method

- **Study area:** the 2019 mining polygons buffered by 5 km
  (`CONFIG.BUFFER_KM`). Images are filtered to Ghana, and the export is
  cropped to the bounding box of the buffered area.
- **Temporal window:** for target year *t*, all dry-season scenes
  (1 November to 30 April) from the seasons ending in *t*−1, *t* and *t*+1.
  The 2000 composite therefore draws on scenes from November 1998 to April 2001.
- **Scene filter:** scene cloud cover below 50 %.
- **Cloud masking:** uses the `QA_PIXEL` band. Cloud and cloud-shadow pixels
  are masked when the flag is set *and* their confidence is medium or high.
  Snow pixels are always masked. For Landsat 8 and 9, cirrus is masked too.
  For Landsat 7 SLC-off scenes, gap pixels are also masked.
- **Scaling:** Collection 2 Level-2 scale factor (×0.0000275 − 0.2), clamped to 0–1.
- **Cross-sensor harmonization:** OLI reflectance is transformed to ETM+ using
  the coefficients of Roy et al. (2016).
- **Sensors per year:**

  | Years | Sensors |
  |---|---|
  | ≤ 2012 | Landsat 5 + 7 |
  | 2013–2021 | Landsat 7 + 8 (harmonized) |
  | ≥ 2022 | Landsat 8 + 9 (harmonized) |

- **Compositing:** medoid (Flood, 2013). For each pixel, the composite keeps
  the actual observation whose six-band reflectance vector lies closest
  (Euclidean distance) to the per-band median.
- **Output grid:** EPSG:32630 (WGS 84 / UTM 30N), 30 m, cloud-optimized GeoTIFF.

### Exported bands

The **27 classifier features** used in the paper are marked ✔. The remaining
bands are exported for quality control only.

| Group | Band | Description | Used |
|---|---|---|---|
| Reflectance | `BLUE`, `GREEN`, `RED`, `NIR`, `SWIR1`, `SWIR2` | Medoid surface reflectance (ETM+ equivalent) | ✔ (6) |
| Vegetation | `NDVI` | (NIR − RED)/(NIR + RED); Tucker (1979) | ✔ |
| Water | `NDWI` | McFeeters (1996) | ✔ |
| | `MNDWI` | Xu (2006) | ✔ |
| | `AWEI_NSH`, `AWEI_SH` | Feyisa et al. (2014) | ✔ (2) |
| | `WATER_PROB` | 0–3 score from agreement of the two AWEI thresholds | ✔ |
| Bare soil / built-up | `BSI` | Rikimaru et al. (2002) | ✔ |
| | `NDBI` | Zha et al. (2003) | ✔ |
| Mineral / colour | `IRON_OXIDE` | RED/BLUE; Segal (1982) | ✔ |
| | `CLAY_INDEX` | (SWIR1 − SWIR2)/(SWIR1 + SWIR2) | ✔ |
| | `CI` | Coloration index; Pouget et al. (1990) | ✔ |
| | `BI` | Brightness index; Pouget et al. (1990) | ✔ |
| Texture (NIR, 3×3) | `ENTROPY`, `NIR_STD_DEV` | Local entropy and standard deviation | ✔ (2) |
| | `NIR_contrast`, `NIR_corr`, `NIR_var` | GLCM contrast, correlation, variance | ✔ (3) |
| Terrain | `ELEVATION`, `SLOPE`, `ASPECT`, `TRI` | SRTM; Farr et al. (2007) | ✔ (4) |
| Night lights | `NIGHTLIGHT`, `NIGHTLIGHT_NORM`, `IS_LIT`, `NIGHTLIGHT_DENSITY` | VIIRS / DMSP | – |
| Quality | `PIXEL_COUNT`, `VALID_PIXEL_COUNT` | Number of clear observations per pixel | – |
| | `CLOUD_PROB` | Temporal SD of BLUE (a proxy for residual contamination, not a probability) | – |
| | `NIR_STD_TEMPORAL`, `RED_STD_TEMPORAL`, `SWIR1_STD_TEMPORAL` | Temporal SD within the 3-year window | – |
| | `COMPOSITE_YEAR`, `DATA_MASK` | Target year; 1 where any observation exists | – |
| Anomaly | `NDVI_ANOMALY` | NDVI minus the mean of the 2000–2004 composites | – |

### How to run

1. Open the Code Editor and paste `01_landsat_composites.js` into a new script.
2. Change the polygon asset path (`miningPolygons2019`) to your own project.
3. Click **Run**. The Console shows image counts per year, and the map shows
   example composites for 2000, 2005, 2010 and 2019.
4. Open the **Tasks** tab and click **Run** on each of the 20 export tasks
   (`Mining_Composite_medoid_2000` … `_2019`). Each export takes roughly
   30 minutes.
5. Download the GeoTIFFs from Google Drive into TODO `data/composites/`, where
   the R scripts expect them.

---

## 02 – Training sample export

### Method

Training points were placed by hand on the 2019 composite. The 2019 RGB
composite was the primary reference. ESA WorldCover 2020 and NDVI were
displayed as visual aids only; WorldCover labels were **not** copied into the
training data. Points were stratified by distance to the 2019 mining polygons:

| Zone | Definition |
|---|---|
| 1 – `Mining_Areas` | Inside the 2019 mining polygons |
| 2 – `Buffer_1_2km` | 0–2 km outside the polygons |
| 3 – `Buffer_2_5km` | 2–5 km outside the polygons |

Labelling classes in this script:

| Code | Class |
|---|---|
| 1 | Water |
| 2 | Dense forest |
| 3 | Woodland / savanna |
| 4 | Urban |
| 5 | Bare soil |
| 6 | Agriculture |
| 0 | Other |

### How to run

1. Paste `02_export_training_samples.js` into the Code Editor and update the
   asset paths.
2. Create point **geometry imports** named `<Class>_<zone>`, for example
   `Water_1`, `Dense_Forest_2` or `Bare_Soil_3`, and place points on the map.
   Missing imports are skipped.
3. Click **Run**, then start the export tasks in the **Tasks** tab. Each run
   produces:
   - one CSV and one asset per class and zone
   - one combined CSV and asset per zone
   - one master dataset with all points (`Manual_Samples_Master_Dataset`)

The hand-placed geometries live only in the Code Editor imports. **The
exported master CSV is the canonical training dataset** and is archived at
TODO DOI. You do not need to rerun this script to reproduce the analysis.

### Output columns

| Column | Description |
|---|---|
| `longitude`, `latitude` | WGS 84 coordinates |
| `zone`, `zone_name` | Distance stratum (see above) |
| `landcover_code`, `class_name` | Label (see above) |
| `sampling_method`, `verified`, `notes` | Provenance fields |
| `timestamp` | Time of **export** (ms since 1970), not time of labelling |

---

## Full references in the paper appendix

## Licence

Code: TODO (for example MIT). Data: see the data archive.
