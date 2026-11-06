# Mining expansion in Ghana, 2000–2019

Annual Landsat composites (Google Earth Engine) are classified into Water / Mining /
Vegetation with one national random forest, constrained with a backward transition rule,
and summarised inside the Maus et al. (2022) mining polygons. The final product is a
yearly mining footprint (Mining + Water) for socio-economic analysis.

## Repository layout

```
R/
  00_config.R            paths, bands, classes, CV and model settings (sourced by all)
  01_utils.R             shared helpers
  preprocessing/         00a tile merge (GEE export) and training data visualisation,
  03_feature_importance_validation.R   band audit
  10_Training_Prediction.R             national RF, spatial CV, yearly rasters  -> P1
  20_Transition_Constraints.R          backward transition rule              -> P2
  30_Uncertainty_Simulations.R         ensemble-draw simulation              -> P3
  40-45b                                validation and benchmarks             -> P4
  50_Visualisations_maps.R, 51, 52      figures and LaTeX tables              -> P5
  60, 61                                district statistics                   -> P5
  70_final_footprint.R                  final data product                    -> P6
```

Run all scripts from the repository root (they call `source("R/00_config.R")`).

## Data layout (outside the repository)

Set `GHANA_MINING_DATA` (e.g. in `~/.Renviron`) to the data folder.

```
<GHANA_MINING_DATA>/
  Training_data/, Mining_polygons/, district_shapefiles/
  Ghana_Tifs/
    B_MERGED/              input composites (one 40-band GeoTIFF per year)
    P1_UNCONSTRAINED/      model, OOF predictions, metrics, yearly rasters
    P2_CONSTRAINED/        constrained yearly rasters + statistics
    P3_UNCERTAINTY/        simulation outputs
    P4_VALIDATION/         validation tables, fold cache, caches
    P5_RESULTS/            everything that goes into the paper
      figures/<NN_topic>/  one folder per producing script
      tables/              .tex tables
      statistics/          key CSVs: area inside polygons per year,
                           district x year panel, district summary
    P6_FINAL_FOOTPRINT/    final product (see its README.txt)
      rasters/footprint_YYYY.tif, onset_year.tif, mining_polygon_mask.tif
      footprint_polygons.gpkg
      tables/
```

## Run order

`03 → 10 → 20 → 30 → 40 → 41 → 42 → 44 → 45 → 45a → 45b → 50 → 51 → 52 → 60 → 61 → 70`

## Requirements

R ≥ 4.3 with `pacman`; packages are loaded with `pacman::p_load()` in each script
(terra, sf, ranger, blockCV, automap, gstat, CAST, torch, ggplot2, patchwork, …).
