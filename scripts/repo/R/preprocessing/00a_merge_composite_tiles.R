# Script to merge GEE Mining Composite tiles and prepare for land cover classification
# Buffered composites and 2019 Ghana-wide composite for model training

# ==============================================================================
# SETUP AND PACKAGE LOADING
# ==============================================================================

# Required packages
packages <- c("terra", "raster", "sf", "stringr", "dplyr", "ggplot2", 
              "RColorBrewer", "gridExtra", "viridis")

# Install missing packages
install_missing <- packages[!packages %in% installed.packages()[,"Package"]]
if(length(install_missing) > 0) {
  install.packages(install_missing)
}

# Load libraries
invisible(lapply(packages, library, character.only = TRUE))

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Set input and output directories
# Data folder: same convention as 00_config.R (env var GHANA_MINING_DATA)
DATA_ROOT  <- Sys.getenv("GHANA_MINING_DATA", unset = "data")
input_dir  <- file.path(DATA_ROOT, "Ghana_Tifs", "A_UNMERGED")   # GEE tiles
output_dir <- file.path(DATA_ROOT, "Ghana_Tifs", "B_MERGED")     # = PATHS$composite_dir

# Years to process
years_to_process <- c(2000:2019)

# Create output directories
dirs_to_create <- c(output_dir, 
                    file.path(output_dir, "visualizations"),
                    file.path(output_dir, "validation_reports"))

for (dir in dirs_to_create) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Created directory:", dir, "\n")
  }
}

# ==============================================================================
# BAND CONFIGURATION (Aligned with GEE composite)
# ==============================================================================

# Expected bands from GEE composite (40 bands total)
EXPECTED_BANDS <- c(
  # Base reflectance (6 bands)
  "BLUE", "GREEN", "RED", "NIR", "SWIR1", "SWIR2",
  # Spectral indices (12 bands)
  "NDVI", "NDWI", "MNDWI", "AWEI_NSH", "AWEI_SH",
  "BSI", "NDBI", "IRON_OXIDE", "CLAY_INDEX", "CI", "BI", "WATER_PROB",
  # Texture (5 bands)
  "ENTROPY", "NIR_STD_DEV", "NIR_contrast", "NIR_corr", "NIR_var",
  # Topography (4 bands)
  "ELEVATION", "SLOPE", "ASPECT", "TRI",
  # Nightlight (4 bands)
  "NIGHTLIGHT", "NIGHTLIGHT_NORM", "IS_LIT", "NIGHTLIGHT_DENSITY",
  # Quality (8 bands)
  "PIXEL_COUNT", "CLOUD_PROB", "COMPOSITE_YEAR", "DATA_MASK",
  "NIR_STD_TEMPORAL", "RED_STD_TEMPORAL", "SWIR1_STD_TEMPORAL",
  "VALID_PIXEL_COUNT",
  # Anomaly (1 band)
  "NDVI_ANOMALY"
)

# Key bands for validation visualization
KEY_VALIDATION_BANDS <- c("NDVI", "BSI", "MNDWI", "NDBI", "AWEI_NSH", "WATER_PROB")

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

# Function to extract year from filename pattern
# Pattern: Mining_Composite_medoid_YYYY-XXXXXXXXXX-XXXXXXXXXX.tif
extract_year <- function(filename) {
  # Match pattern: medoid_YYYY- where YYYY is the year
  year_match <- str_match(filename, "medoid_(\\d{4})")
  
  if (!is.na(year_match[1, 2])) {
    return(as.numeric(year_match[1, 2]))
  }
  
  # Fallback: try to find any 4-digit year between 2000-2025
  year_match <- str_extract(filename, "20[0-2][0-9]")
  return(as.numeric(year_match))
}

# Function to group files by year
group_files_by_year <- function(input_dir) {
  # List all TIFF files matching the NEW naming pattern
  # Pattern: Mining_Composite_medoid_YYYY-XXXXXXXXXX-XXXXXXXXXX.tif
  all_files <- list.files(input_dir, 
                          pattern = "Mining_Composite_medoid_\\d{4}.*\\.(tif|tiff)$", 
                          full.names = TRUE, 
                          ignore.case = TRUE)
  
  if (length(all_files) == 0) {
    # Try alternative patterns
    cat("No files found with primary pattern. Trying alternatives...\n")
    all_files <- list.files(input_dir, 
                            pattern = "Mining_Composite.*medoid.*\\.(tif|tiff)$", 
                            full.names = TRUE, 
                            ignore.case = TRUE)
  }
  
  if (length(all_files) == 0) {
    stop("No TIFF files found in the input directory!\n",
         "Expected pattern: Mining_Composite_medoid_YYYY-XXXXXXXXXX-XXXXXXXXXX.tif")
  }
  
  cat("Found", length(all_files), "total tile files\n")
  
  # Create a data frame with file info
  file_df <- data.frame(
    filepath = all_files,
    filename = basename(all_files),
    stringsAsFactors = FALSE
  )
  
  # Extract year from each filename
  file_df$year <- sapply(file_df$filename, extract_year)
  
  # Extract tile coordinates for debugging
  file_df$tile_coords <- str_extract(file_df$filename, "-\\d{10}-\\d{10}")
  
  # Remove files without valid years
  invalid_files <- file_df[is.na(file_df$year), ]
  if (nrow(invalid_files) > 0) {
    cat("\nWARNING: Could not extract year from these files:\n")
    cat(paste("  -", invalid_files$filename, collapse = "\n"), "\n")
  }
  
  file_df <- file_df[!is.na(file_df$year), ]
  
  # Group by year
  files_by_year <- split(file_df$filepath, file_df$year)
  
  # Report tile counts per year
  cat("\nTiles per year:\n")
  for (year in sort(as.numeric(names(files_by_year)))) {
    year_files <- files_by_year[[as.character(year)]]
    cat(sprintf("  %d: %d tiles\n", year, length(year_files)))
  }
  
  return(files_by_year)
}

# ==============================================================================
# VALIDATION FUNCTIONS
# ==============================================================================

# Function to validate CRS consistency
check_crs_consistency <- function(raster_list) {
  crs_list <- lapply(raster_list, terra::crs)
  
  # Check if all CRS are the same
  all_same <- all(sapply(crs_list[-1], function(x) identical(x, crs_list[[1]])))
  
  if (!all_same) {
    warning("CRS mismatch detected among tiles!")
    cat("Different CRS found:\n")
    unique_crs <- unique(crs_list)
    for (i in seq_along(unique_crs)) {
      cat(paste("CRS", i, ":", unique_crs[[i]], "\n"))
    }
  }
  
  return(list(consistent = all_same, crs = crs_list[[1]]))
}

# Function to check band alignment with expected bands
check_band_alignment <- function(raster_obj, year) {
  cat("\n--- Band Alignment Check for Year", year, "---\n")
  
  actual_bands <- names(raster_obj)
  n_actual <- length(actual_bands)
  n_expected <- length(EXPECTED_BANDS)
  
  cat(sprintf("  Actual bands: %d\n", n_actual))
  cat(sprintf("  Expected bands: %d\n", n_expected))
  
  if (n_actual == n_expected) {
    cat("  b Band count matches!\n")
    
    # Check if names match (if bands are named)
    if (!all(actual_bands == paste0("Band_", 1:n_actual))) {
      # Bands have meaningful names, compare
      matching <- sum(actual_bands == EXPECTED_BANDS)
      cat(sprintf("  b %d/%d band names match expected\n", matching, n_expected))
    }
  } else if (n_actual < n_expected) {
    cat(sprintf("  b  WARNING: Missing %d bands!\n", n_expected - n_actual))
    cat("  Missing bands may include:", 
        paste(EXPECTED_BANDS[(n_actual+1):n_expected], collapse = ", "), "\n")
  } else {
    cat(sprintf("  b  WARNING: %d extra bands found!\n", n_actual - n_expected))
  }
  
  return(list(
    n_actual = n_actual,
    n_expected = n_expected,
    aligned = n_actual == n_expected
  ))
}

# Function to check for NoData values and statistics
check_nodata_stats <- function(raster_obj, year, bands_to_check = NULL) {
  cat("\n--- NoData and Statistics Check for Year", year, "---\n")
  
  # Get basic information
  ncells <- terra::ncell(raster_obj)
  band_names <- names(raster_obj)
  
  # If bands_to_check not specified, check key bands or first 10
  if (is.null(bands_to_check)) {
    key_indices <- which(band_names %in% KEY_VALIDATION_BANDS)
    if (length(key_indices) > 0) {
      bands_to_check <- key_indices
    } else {
      bands_to_check <- 1:min(10, terra::nlyr(raster_obj))
    }
  }
  
  # Check each specified band for NoData
  band_stats <- list()
  
  for (i in bands_to_check) {
    band_name <- band_names[i]
    band_data <- raster_obj[[i]]
    
    # Calculate statistics
    stats <- terra::global(band_data, c("min", "max", "mean", "sd"), na.rm = TRUE)
    
    # Count NA values
    na_count <- terra::global(is.na(band_data), "sum")[[1]]
    na_percent <- (na_count / ncells) * 100
    
    band_stats[[band_name]] <- list(
      min = stats$min,
      max = stats$max,
      mean = stats$mean,
      sd = stats$sd,
      na_count = na_count,
      na_percent = na_percent
    )
    
    cat(sprintf("  Band %02d %-18s: Min=%8.4f, Max=%8.4f, Mean=%8.4f, NA%%=%5.2f%%\n",
                i, band_name, stats$min, stats$max, stats$mean, na_percent))
  }
  
  return(band_stats)
}

# Function to validate value ranges for key bands
validate_band_ranges <- function(raster_obj, year) {
  cat("\n--- Value Range Validation for Year", year, "---\n")
  
  band_names <- names(raster_obj)
  
  # Define expected ranges for key bands
  expected_ranges <- list(
    # Normalized indices should be [-1, 1]
    NDVI = c(-1, 1),
    NDWI = c(-1, 1),
    MNDWI = c(-1, 1),
    BSI = c(-1, 1),
    NDBI = c(-1, 1),
    CI = c(-1, 1),
    CLAY_INDEX = c(-1, 1),
    
    # AWEI has extended range
    AWEI_NSH = c(-10, 5),
    AWEI_SH = c(-10, 5),
    
    # Other indices
    BI = c(0, 1),
    WATER_PROB = c(0, 3),
    IRON_OXIDE = c(0, 5),
    
    # Reflectance should be [0, 1]
    BLUE = c(0, 1),
    GREEN = c(0, 1),
    RED = c(0, 1),
    NIR = c(0, 1),
    SWIR1 = c(0, 1),
    SWIR2 = c(0, 1),
    
    # Topographic
    ELEVATION = c(-50, 1000),
    SLOPE = c(0, 90),
    ASPECT = c(0, 360)
  )
  
  issues <- list()
  
  for (band_name in names(expected_ranges)) {
    if (band_name %in% band_names) {
      expected <- expected_ranges[[band_name]]
      actual_stats <- terra::global(raster_obj[[band_name]], c("min", "max"), na.rm = TRUE)
      
      out_of_range <- FALSE
      if (actual_stats$min < expected[1] * 1.1 - 0.1 || actual_stats$max > expected[2] * 1.1 + 0.1) {
        out_of_range <- TRUE
        issues[[band_name]] <- list(
          expected = expected,
          actual = c(actual_stats$min, actual_stats$max)
        )
        cat(sprintf("  b  %s: Expected [%.2f, %.2f], Got [%.4f, %.4f]\n",
                    band_name, expected[1], expected[2], actual_stats$min, actual_stats$max))
      }
    }
  }
  
  if (length(issues) == 0) {
    cat("  b All checked bands within expected ranges\n")
  } else {
    cat(sprintf("  b  %d bands with values outside expected ranges\n", length(issues)))
  }
  
  return(issues)
}

# Function to validate merged raster
validate_merged_raster <- function(merged_raster, original_tiles, year, output_dir) {
  cat("\n========================================\n")
  cat("VALIDATION REPORT FOR YEAR", year, "\n")
  cat("========================================\n")
  
  validation_results <- list()
  
  # 1. Check dimensions and resolution
  cat("\n1. SPATIAL PROPERTIES:\n")
  cat("   Dimensions (rows, cols):", terra::nrow(merged_raster), "x", terra::ncol(merged_raster), "\n")
  cat("   Resolution:", terra::res(merged_raster), "\n")
  cat("   Extent:", as.vector(terra::ext(merged_raster)), "\n")
  
  validation_results$dimensions <- c(terra::nrow(merged_raster), terra::ncol(merged_raster))
  validation_results$resolution <- terra::res(merged_raster)
  validation_results$extent <- as.vector(terra::ext(merged_raster))
  
  # 2. Check CRS
  cat("\n2. COORDINATE REFERENCE SYSTEM:\n")
  crs_info <- terra::crs(merged_raster, describe = TRUE)
  cat("   CRS Name:", crs_info$name, "\n")
  cat("   EPSG Code:", crs_info$code, "\n")
  cat("   Units:", crs_info$units, "\n")
  cat("   Is Geographic:", terra::is.lonlat(merged_raster), "\n")
  
  validation_results$crs <- crs_info
  validation_results$is_geographic <- terra::is.lonlat(merged_raster)
  
  # 3. Check band alignment
  band_alignment <- check_band_alignment(merged_raster, year)
  validation_results$band_alignment <- band_alignment
  
  # 4. Assign band names if count matches
  band_names <- names(merged_raster)
  if (band_alignment$aligned) {
    names(merged_raster) <- EXPECTED_BANDS
    band_names <- EXPECTED_BANDS
    cat("\n3. BAND NAMES:\n")
    cat("   b Assigned", length(EXPECTED_BANDS), "band names from template\n")
  } else {
    cat("\n3. BAND NAMES:\n")
    cat("   b  Could not assign band names (count mismatch)\n")
    cat("   Current names:", paste(head(band_names, 10), collapse = ", "), "...\n")
  }
  
  validation_results$bands <- list(
    total = terra::nlyr(merged_raster),
    names = band_names
  )
  
  # 5. Check for NoData and statistics
  band_stats <- check_nodata_stats(merged_raster, year)
  validation_results$band_stats <- band_stats
  
  # 6. Validate value ranges
  range_issues <- validate_band_ranges(merged_raster, year)
  validation_results$range_issues <- range_issues
  
  # 7. Check data type
  cat("\n6. DATA TYPE:\n")
  data_type <- terra::datatype(merged_raster)
  cat("   Data type:", data_type, "\n")
  
  validation_results$data_type <- data_type
  
  # 8. File size comparison
  cat("\n7. FILE SIZE:\n")
  merged_file <- file.path(output_dir, paste0("Mining_Composite_", year, ".tif"))
  if (file.exists(merged_file)) {
    file_size_mb <- file.size(merged_file) / (1024^2)
    cat("   Merged file size:", round(file_size_mb, 2), "MB\n")
    
    # Calculate total size of original tiles
    original_size_mb <- sum(file.size(original_tiles)) / (1024^2)
    cat("   Total original tiles size:", round(original_size_mb, 2), "MB\n")
    cat("   Compression ratio:", round((file_size_mb / original_size_mb) * 100, 1), "%\n")
    
    validation_results$file_size <- list(
      merged = file_size_mb,
      original = original_size_mb,
      compression_ratio = (file_size_mb / original_size_mb) * 100
    )
  }
  
  # Save validation report
  report_file <- file.path(output_dir, "validation_reports", 
                           paste0("validation_report_", year, ".txt"))
  
  sink(report_file)
  cat("VALIDATION REPORT FOR YEAR", year, "\n")
  cat("Generated:", as.character(Sys.time()), "\n")
  cat("Expected bands:", length(EXPECTED_BANDS), "\n")
  cat("========================================\n\n")
  
  for (name in names(validation_results)) {
    cat(toupper(name), ":\n")
    print(validation_results[[name]])
    cat("\n")
  }
  
  # Print full band list
  cat("\nFULL BAND LIST:\n")
  for (i in 1:length(band_names)) {
    cat(sprintf("  %02d. %s\n", i, band_names[i]))
  }
  
  sink()
  
  cat("\nValidation report saved to:", report_file, "\n")
  
  return(validation_results)
}

# ==============================================================================
# MERGING FUNCTION
# ==============================================================================

merge_tiles_for_year <- function(tile_files, year, output_dir) {
  cat("\n========================================\n")
  cat("Processing Year:", year, "\n")
  cat("========================================\n")
  cat("Number of tiles to merge:", length(tile_files), "\n")
  
  # Print tile filenames for verification
  cat("Tiles:\n")
  for (f in tile_files) {
    cat("  -", basename(f), "\n")
  }
  
  # Output name read by the analysis scripts (00_config.R): Mining_Composite_YYYY.tif
  output_file <- file.path(output_dir, paste0("Mining_Composite_", year, ".tif"))
  
  # Check if output already exists
  if (file.exists(output_file)) {
    cat("Output file already exists. Overwriting:", output_file, "\n")
  }
  
  tryCatch({
    # Read all tiles
    cat("\nReading tiles...\n")
    tiles <- lapply(tile_files, function(f) {
      cat("  Reading:", basename(f), "...")
      r <- terra::rast(f)
      cat(" (", terra::nlyr(r), "bands)\n")
      return(r)
    })
    
    # Check band count consistency across tiles
    band_counts <- sapply(tiles, terra::nlyr)
    if (length(unique(band_counts)) > 1) {
      cat("\nb  WARNING: Inconsistent band counts across tiles!\n")
      cat("  Band counts:", paste(band_counts, collapse = ", "), "\n")
      cat("  This may indicate corrupted or incomplete tiles.\n")
    } else {
      cat("\nb All tiles have", unique(band_counts), "bands\n")
    }
    
    # Check CRS consistency
    cat("\nChecking CRS consistency...\n")
    crs_check <- check_crs_consistency(tiles)
    
    if (!crs_check$consistent) {
      cat("WARNING: CRS mismatch detected. Reprojecting all tiles to first tile's CRS...\n")
      target_crs <- terra::crs(tiles[[1]])
      tiles <- lapply(tiles, function(r) {
        if (!identical(terra::crs(r), target_crs)) {
          terra::project(r, target_crs)
        } else {
          r
        }
      })
    }
    
    # Create a SpatRasterCollection
    cat("Creating raster collection...\n")
    src <- terra::sprc(tiles)
    
    # Merge tiles
    cat("Merging tiles (this may take a while)...\n")
    merged <- terra::merge(src)
    
    # Assign band names if count matches
    if (terra::nlyr(merged) == length(EXPECTED_BANDS)) {
      names(merged) <- EXPECTED_BANDS
      cat("b Assigned", length(EXPECTED_BANDS), "band names\n")
    } else {
      cat("b  Band count mismatch - using generic names\n")
      names(merged) <- paste0("Band_", 1:terra::nlyr(merged))
    }
    
    # Handle NoData values if needed
    if (is.na(terra::NAflag(merged))) {
      terra::NAflag(merged) <- -9999
      cat("Set NoData value to -9999\n")
    }
    
    # Write output with compression
    cat("Writing merged raster to disk...\n")
    terra::writeRaster(merged, 
                       output_file, 
                       overwrite = TRUE,
                       datatype = "FLT4S",
                       gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"),
                       NAflag = -9999)
    
    cat("b Successfully merged and saved:", output_file, "\n")
    
    # Validate the merged raster
    validation_results <- validate_merged_raster(merged, tile_files, year, output_dir)
    
    # Clean up memory
    rm(tiles, src)
    gc()
    
    return(list(success = TRUE, merged = merged, validation = validation_results))
    
  }, error = function(e) {
    cat("ERROR during merging:", e$message, "\n")
    return(list(success = FALSE, error = e$message))
  })
}

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

create_comprehensive_visualization <- function(merged_raster, year, output_dir) {
  cat("\nCreating visualizations for year", year, "...\n")
  
  viz_dir <- file.path(output_dir, "visualizations")
  
  # Create both PDF and PNG outputs
  pdf_file <- file.path(viz_dir, paste0("Mining_Composite_", year, "_visualization.pdf"))
  png_file <- file.path(viz_dir, paste0("Mining_Composite_", year, "_overview.png"))
  
  # Get band names
  band_names <- names(merged_raster)
  
  # Start PDF
  pdf(pdf_file, width = 14, height = 10)
  
  # Page 1: RGB and False Color Composites
  if (all(c("RED", "GREEN", "BLUE", "NIR") %in% band_names)) {
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
    
    # Convert to raster for plotRGB
    r_stack <- raster::stack(merged_raster)
    
    # 1. Natural Color RGB
    raster::plotRGB(r_stack, r = "RED", g = "GREEN", b = "BLUE", 
                    stretch = "lin", main = paste("Natural Color RGB -", year))
    
    # 2. False Color (NIR-RED-GREEN)
    raster::plotRGB(r_stack, r = "NIR", g = "RED", b = "GREEN", 
                    stretch = "lin", main = paste("False Color (NIR-R-G) -", year))
    
    # 3. SWIR composite if available
    if (all(c("SWIR1", "SWIR2") %in% band_names)) {
      raster::plotRGB(r_stack, r = "SWIR1", g = "NIR", b = "RED", 
                      stretch = "lin", main = paste("SWIR Composite -", year))
    }
    
    # 4. Add extent and CRS info
    plot.new()
    ext <- terra::ext(merged_raster)
    crs_info <- terra::crs(merged_raster, describe = TRUE)
    text(0.5, 0.85, paste("Mining Composite", year), cex = 1.5, font = 2)
    text(0.5, 0.70, paste("Extent:", round(ext[1], 2), "-", round(ext[2], 2), 
                          "(E),", round(ext[3], 2), "-", round(ext[4], 2), "(N)"), cex = 0.9)
    text(0.5, 0.60, paste("CRS:", crs_info$name), cex = 0.9)
    text(0.5, 0.50, paste("Dimensions:", terra::nrow(merged_raster), "rows x", 
                          terra::ncol(merged_raster), "cols"), cex = 0.9)
    text(0.5, 0.40, paste("Resolution:", paste(round(terra::res(merged_raster), 6), collapse = " x ")), cex = 0.9)
    text(0.5, 0.30, paste("Bands:", terra::nlyr(merged_raster)), cex = 0.9)
    text(0.5, 0.20, paste("Is Geographic:", terra::is.lonlat(merged_raster)), cex = 0.9)
  }
  
  # Page 2: Mining-specific indices (UPDATED for new bands)
  indices_to_plot <- c("NDVI", "BSI", "MNDWI", "NDBI", "AWEI_NSH", "WATER_PROB")
  available_indices <- indices_to_plot[indices_to_plot %in% band_names]
  
  if (length(available_indices) > 0) {
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 5))
    
    for (idx in available_indices) {
      r_idx <- raster::raster(merged_raster[[idx]])
      
      # Choose appropriate color palette
      cols <- switch(idx,
                     "NDVI" = colorRampPalette(c("brown", "yellow", "darkgreen"))(100),
                     "MNDWI" = colorRampPalette(c("brown", "white", "blue", "darkblue"))(100),
                     "BSI" = colorRampPalette(c("darkgreen", "yellow", "orange", "brown"))(100),
                     "NDBI" = colorRampPalette(c("darkgreen", "yellow", "orange", "red"))(100),
                     "AWEI_NSH" = colorRampPalette(c("brown", "white", "lightblue", "darkblue"))(100),
                     "WATER_PROB" = colorRampPalette(c("tan", "lightblue", "blue", "darkblue"))(4),
                     viridis::viridis(100)
      )
      
      plot(r_idx, main = paste(idx, "-", year), col = cols, axes = FALSE)
      
      # Add basic statistics as text
      stats <- terra::global(merged_raster[[idx]], c("mean", "sd"), na.rm = TRUE)
      legend("topright", 
             legend = c(paste("Mean:", round(stats$mean, 3)),
                        paste("SD:", round(stats$sd, 3))),
             bty = "n", cex = 0.8)
    }
  }
  
  # Page 3: Additional indices (CI, BI, IRON_OXIDE, CLAY_INDEX)
  additional_indices <- c("CI", "BI", "IRON_OXIDE", "CLAY_INDEX", "AWEI_SH", "NDVI_ANOMALY")
  available_additional <- additional_indices[additional_indices %in% band_names]
  
  if (length(available_additional) > 0) {
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 5))
    
    for (idx in available_additional) {
      r_idx <- raster::raster(merged_raster[[idx]])
      
      cols <- switch(idx,
                     "CI" = colorRampPalette(c("blue", "white", "red"))(100),
                     "BI" = colorRampPalette(c("black", "gray", "white"))(100),
                     "IRON_OXIDE" = colorRampPalette(c("gray", "yellow", "orange", "red"))(100),
                     "CLAY_INDEX" = colorRampPalette(c("blue", "white", "brown"))(100),
                     "AWEI_SH" = colorRampPalette(c("brown", "white", "lightblue", "darkblue"))(100),
                     "NDVI_ANOMALY" = colorRampPalette(c("red", "white", "darkgreen"))(100),
                     viridis::viridis(100)
      )
      
      plot(r_idx, main = paste(idx, "-", year), col = cols, axes = FALSE)
      
      stats <- terra::global(merged_raster[[idx]], c("mean", "sd"), na.rm = TRUE)
      legend("topright", 
             legend = c(paste("Mean:", round(stats$mean, 3)),
                        paste("SD:", round(stats$sd, 3))),
             bty = "n", cex = 0.8)
    }
  }
  
  # Page 4: Histogram distributions for key bands
  par(mfrow = c(3, 3), mar = c(4, 4, 3, 2))
  
  # Plot histograms for key bands
  bands_to_hist <- c("NDVI", "BSI", "MNDWI", "NDBI", "AWEI_NSH", "AWEI_SH", 
                     "IRON_OXIDE", "WATER_PROB", "NDVI_ANOMALY")
  bands_to_hist <- bands_to_hist[bands_to_hist %in% band_names]
  
  for (band in bands_to_hist[1:min(9, length(bands_to_hist))]) {
    # Sample data for histogram (to avoid memory issues)
    band_data <- terra::spatSample(merged_raster[[band]], 10000, method = "regular", na.rm = TRUE)
    
    if (nrow(band_data) > 0) {
      hist(band_data[,1], 
           main = paste("Histogram:", band),
           xlab = "Value",
           ylab = "Frequency",
           col = "lightblue",
           breaks = 50)
      
      # Add mean and median lines
      abline(v = mean(band_data[,1], na.rm = TRUE), col = "red", lwd = 2)
      abline(v = median(band_data[,1], na.rm = TRUE), col = "blue", lwd = 2)
      
      legend("topright", 
             legend = c("Mean", "Median"),
             col = c("red", "blue"),
             lwd = 2,
             bty = "n",
             cex = 0.7)
    }
  }
  
  dev.off()
  
  # Create a quick overview PNG
  png(png_file, width = 1600, height = 1200, res = 150)
  
  if (all(c("RED", "GREEN", "BLUE") %in% band_names)) {
    r_stack <- raster::stack(merged_raster)
    raster::plotRGB(r_stack, r = "RED", g = "GREEN", b = "BLUE", 
                    stretch = "lin", main = paste("Mining Composite medoid -", year))
  } else {
    plot(raster::raster(merged_raster[[1]]), 
         main = paste("Mining Composite medoid -", year, "- Band 1"))
  }
  
  dev.off()
  
  cat("Visualizations saved to:\n")
  cat("  PDF:", pdf_file, "\n")
  cat("  PNG:", png_file, "\n")
}

# ==============================================================================
# CROSS-YEAR CONSISTENCY CHECK
# ==============================================================================

check_cross_year_consistency <- function(output_dir, years) {
  cat("\n========================================\n")
  cat("CROSS-YEAR CONSISTENCY CHECK\n")
  cat("========================================\n")
  
  consistency_results <- data.frame()
  
  for (year in years) {
    merged_file <- file.path(output_dir, paste0("Mining_Composite_", year, ".tif"))
    
    if (file.exists(merged_file)) {
      r <- terra::rast(merged_file)
      
      row_data <- data.frame(
        year = year,
        n_bands = terra::nlyr(r),
        n_rows = terra::nrow(r),
        n_cols = terra::ncol(r),
        res_x = terra::res(r)[1],
        res_y = terra::res(r)[2],
        crs = terra::crs(r, describe = TRUE)$code,
        file_size_mb = round(file.size(merged_file) / (1024^2), 2)
      )
      
      consistency_results <- rbind(consistency_results, row_data)
      rm(r)
    }
  }
  
  # Check for inconsistencies
  cat("\nYear-by-year summary:\n")
  print(consistency_results)
  
  # Flag any inconsistencies
  if (length(unique(consistency_results$n_bands)) > 1) {
    cat("\nb  WARNING: Inconsistent band counts across years!\n")
  }
  if (length(unique(consistency_results$res_x)) > 1) {
    cat("\nb  WARNING: Inconsistent resolutions across years!\n")
  }
  if (length(unique(consistency_results$crs)) > 1) {
    cat("\nb  WARNING: Inconsistent CRS across years!\n")
  }
  
  # Save consistency report
  write.csv(consistency_results, 
            file.path(output_dir, "validation_reports", "cross_year_consistency.csv"),
            row.names = FALSE)
  
  return(consistency_results)
}

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================

cat("\n==========================================\n")
cat("GEE MINING COMPOSITE MERGER AND VALIDATOR\n")
cat("==========================================\n")
cat("Start time:", as.character(Sys.time()), "\n")
cat("Input directory:", input_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Years to process:", paste(years_to_process, collapse = ", "), "\n")
cat("Expected bands per composite:", length(EXPECTED_BANDS), "\n")
cat("\nFile pattern: Mining_Composite_medoid_YYYY-XXXXXXXXXX-XXXXXXXXXX.tif\n")

# Group files by year
cat("\nGrouping files by year...\n")
files_by_year <- group_files_by_year(input_dir)

# Display found files
cat("\nFiles found per year:\n")
for (year in names(files_by_year)) {
  cat("  Year", year, ":", length(files_by_year[[year]]), "tiles\n")
}

# Process each year
all_results <- list()

for (year in years_to_process) {
  year_char <- as.character(year)
  
  if (!year_char %in% names(files_by_year)) {
    cat("\nWARNING: No files found for year", year, "- skipping\n")
    next
  }
  
  if (length(files_by_year[[year_char]]) != 6) {
    cat("\nWARNING: Expected 6 tiles for year", year, "but found", 
        length(files_by_year[[year_char]]), "- proceeding anyway\n")
  }
  
  # Merge tiles for this year
  result <- merge_tiles_for_year(files_by_year[[year_char]], year, output_dir)
  
  if (result$success) {
    # Create visualizations
    create_comprehensive_visualization(result$merged, year, output_dir)
    
    all_results[[year_char]] <- result
    
    # Clean up memory
    rm(result)
    gc()
  } else {
    cat("ERROR: Failed to process year", year, "\n")
    all_results[[year_char]] <- result
  }
}

# ==============================================================================
# CROSS-YEAR CONSISTENCY CHECK
# ==============================================================================

processed_years <- as.numeric(names(all_results)[sapply(all_results, function(x) x$success)])
if (length(processed_years) > 1) {
  consistency <- check_cross_year_consistency(output_dir, processed_years)
}

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n==========================================\n")
cat("PROCESSING COMPLETE\n")
cat("==========================================\n")

# Summary statistics
successful <- sum(sapply(all_results, function(x) x$success))
failed <- length(all_results) - successful

cat("\nSummary:\n")
cat("  Successfully processed:", successful, "years\n")
cat("  Failed:", failed, "years\n")

if (successful > 0) {
  cat("\nMerged files created:\n")
  merged_files <- list.files(output_dir, pattern = "^Mining_Composite_\\d{4}\\.tif$", full.names = FALSE)
  for (f in merged_files) {
    fpath <- file.path(output_dir, f)
    size_mb <- file.size(fpath) / (1024^2)
    cat("  -", f, "(", round(size_mb, 2), "MB)\n")
  }
  
  cat("\nVisualization files in:", file.path(output_dir, "visualizations"), "\n")
  cat("Validation reports in:", file.path(output_dir, "validation_reports"), "\n")
}

cat("\nEnd time:", as.character(Sys.time()), "\n")

# Save session info for reproducibility
session_file <- file.path(output_dir, "session_info.txt")
sink(session_file)
cat("Session Information\n")
cat("==================\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("R Version:", R.version.string, "\n\n")
cat("Expected bands:", length(EXPECTED_BANDS), "\n")
cat("Band names:\n")
for (i in 1:length(EXPECTED_BANDS)) {
  cat(sprintf("  %02d. %s\n", i, EXPECTED_BANDS[i]))
}
cat("\n")
print(sessionInfo())
sink()

cat("\nSession info saved to:", session_file, "\n")
cat("\nAll processing complete!\n")

# =============================================================
# 2019 GHANA-WIDE COMPOSITE
#
# Outputs in B_MERGED/:
#   Mining_Composite_2019_full.tif
# =============================================================

if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")
library(terra)

# --- CONFIG -------------------------------------------------------
TILES_DIR  <- file.path(DATA_ROOT, "Ghana_Tifs", "A_UNMERGED", "2019_full")
B_MERGED   <- output_dir
OUT_TIF    <- file.path(B_MERGED, "Mining_Composite_2019_full.tif")
OUT_VRT    <- file.path(B_MERGED, "Mining_Composite_2019_full.vrt")
PERZONE_DIR <- file.path(B_MERGED, "perzone_2019_qa")  # for diagnostics

GDAL_OPTS <- c("COMPRESS=DEFLATE","PREDICTOR=2","TILED=YES",
               "BLOCKXSIZE=512","BLOCKYSIZE=512","BIGTIFF=YES")

ZONES <- c("Northern","ForestTransition","Coastal")

EXPECTED_BANDS <- c(
  "BLUE","GREEN","RED","NIR","SWIR1","SWIR2",
  "NDVI","NDWI","MNDWI","AWEI_NSH","AWEI_SH",
  "BSI","NDBI","IRON_OXIDE","CLAY_INDEX","CI","BI","WATER_PROB",
  "ENTROPY","NIR_STD_DEV","NIR_contrast","NIR_corr","NIR_var",
  "ELEVATION","SLOPE","ASPECT","TRI",
  "NIGHTLIGHT","NIGHTLIGHT_NORM","IS_LIT","NIGHTLIGHT_DENSITY",
  "PIXEL_COUNT","CLOUD_PROB","COMPOSITE_YEAR","DATA_MASK",
  "NIR_STD_TEMPORAL","RED_STD_TEMPORAL","SWIR1_STD_TEMPORAL",
  "VALID_PIXEL_COUNT","NDVI_ANOMALY"
)
# -----------------------------------------------------------------

dir.create(B_MERGED, showWarnings = FALSE, recursive = TRUE)
dir.create(PERZONE_DIR, showWarnings = FALSE, recursive = TRUE)

# --- 1. Inventory tiles ------------------------------------------
all_tiles <- list.files(TILES_DIR,
                        pattern = "^Mining_Composite_2019_Zone_.*\\.tif$", full.names = TRUE)
if (length(all_tiles) == 0) stop("No tiles found in ", TILES_DIR)

cat(sprintf("Total tiles found: %d\n", length(all_tiles)))
for (z in ZONES) {
  n <- sum(grepl(sprintf("Zone_%s-", z), basename(all_tiles)))
  cat(sprintf("  %-18s %d tile(s)\n", z, n))
}

# Reference band names from the first tile (used to restore names on VRT)
ref_bands <- names(rast(all_tiles[1]))
cat(sprintf("\nFirst tile bands: %d\n", length(ref_bands)))
missing_bands <- setdiff(EXPECTED_BANDS, ref_bands)
extra_bands   <- setdiff(ref_bands,       EXPECTED_BANDS)
if (length(missing_bands) > 0)
  warning("Missing expected bands: ", paste(missing_bands, collapse = ", "))
if (length(extra_bands) > 0)
  cat("Extra bands beyond expected: ", paste(extra_bands, collapse=", "), "\n")
if (length(missing_bands) == 0 && length(extra_bands) == 0)
  cat("All 40 expected bands present.\n")

# --- 2. Per-zone VRTs (optional, for zone-specific QA) -----------
for (z in ZONES) {
  z_tiles <- all_tiles[grepl(sprintf("Zone_%s-", z), basename(all_tiles))]
  if (length(z_tiles) == 0) next
  p <- file.path(PERZONE_DIR, sprintf("zone_%s_2019.vrt", z))
  v <- vrt(z_tiles, p, overwrite = TRUE)
  if (length(ref_bands) == nlyr(v)) names(v) <- ref_bands
  cat(sprintf("[%s] VRT: %s  (%d x %d, %d bands)\n",
              z, p, nrow(v), ncol(v), nlyr(v)))
}

# --- 3. Build Ghana-wide VRT then write GeoTIFF ------------------
cat("\nBuilding Ghana-wide VRT...\n")
ghana <- vrt(all_tiles, OUT_VRT, overwrite = TRUE)
if (length(ref_bands) == nlyr(ghana)) names(ghana) <- ref_bands
cat(sprintf("  VRT path: %s\n", OUT_VRT))
cat(sprintf("  Dims:     %d rows x %d cols x %d bands\n",
            nrow(ghana), ncol(ghana), nlyr(ghana)))
cat(sprintf("  CRS:      %s   Resolution: %s m\n",
            crs(ghana, describe = TRUE)$code,
            paste(round(res(ghana),2), collapse=" x ")))

cat("\nWriting Ghana-wide GeoTIFF (this can take ~10-20 minutes)...\n")
writeRaster(ghana, OUT_TIF,
            datatype = "FLT4S", gdal = GDAL_OPTS, overwrite = TRUE)
size_gb <- file.info(OUT_TIF)$size / 1e9
cat(sprintf("  Written: %s  (%.2f GB)\n", OUT_TIF, size_gb))

# Reload from disk to verify
ghana <- rast(OUT_TIF)
cat(sprintf("  Reloaded OK: %d bands, dims %d x %d\n",
            nlyr(ghana), nrow(ghana), ncol(ghana)))

# --- 4. Sanity check: training-point coverage --------------------
TRAIN_CSV <- file.path(DATA_ROOT, "Training_data", "Mining_Training_Data_2019.csv")
if (file.exists(TRAIN_CSV)) {
  cat("\n=== Training-point coverage check ===\n")
  train <- read.csv(TRAIN_CSV)
  pts <- vect(train, geom=c("longitude","latitude"), crs="EPSG:4326") |>
    project(crs(ghana))
  vals <- terra::extract(ghana, pts, ID = FALSE)
  na_per_pt <- rowSums(is.na(vals))
  any_na <- na_per_pt > 0
  all_na <- na_per_pt == nlyr(ghana)
  cat(sprintf("  Points with any NA band:        %d / %d (%.1f%%)\n",
              sum(any_na), nrow(train), 100*sum(any_na)/nrow(train)))
  cat(sprintf("  Points fully outside coverage:  %d / %d (%.1f%%)\n",
              sum(all_na), nrow(train), 100*sum(all_na)/nrow(train)))
  cat("\n  By eco_zone x source:\n")
  train$any_na <- any_na
  print(aggregate(any_na ~ eco_zone + source, data = train,
                  FUN = function(x) sprintf("%d/%d", sum(x), length(x))))
  
  # Compare against buffered 2019 composite to quantify the difference
  BUFF_TIF <- file.path(B_MERGED, "Mining_Composite_2019.tif")
  if (file.exists(BUFF_TIF)) {
    buff <- rast(BUFF_TIF)
    vals_b <- terra::extract(buff, pts, ID = FALSE)
    any_na_b <- rowSums(is.na(vals_b)) > 0
    cat(sprintf("\n  Comparison: %d / %d (%.1f%%) points NA on BUFFERED composite\n",
                sum(any_na_b), nrow(train),
                100*sum(any_na_b)/nrow(train)))
    cat(sprintf("  This is the gap the Ghana-wide composite fixes.\n"))
  }
}

cat("\nDone. Use 'Mining_Composite_2019_full.tif' for training.\n")
cat("Predictions remain on the buffered 'Mining_Composite_<year>.tif' files.\n")