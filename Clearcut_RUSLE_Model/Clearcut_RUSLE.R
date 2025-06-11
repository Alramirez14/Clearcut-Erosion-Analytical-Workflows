###############################################################
######### This script runs the RUSLE model for each year   ####
######### within each setup study site                     ####
######### Created 05/25/2025 by Clearcut Erosion Modeling  ####
###############################################################


library(terra)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(ggpmisc)

# Set folders
r_folder <- "processed_Rwisch_factors"      # Folder with R rasters (one per year)
c_folder <- "processed_C_factors"      # Folder with C rasters (one per year)
output_rusle_folder <- "rusle_outputs" # Folder for RUSLE rasters (optional)
dir.create(output_rusle_folder, showWarnings = FALSE)

# Load static rasters
K <- rast("processed_K_factor/K.tif")
LS <- rast("processed_LS_factor/LS.tif")

# List R and C files
r_files <- list.files(r_folder, pattern = "\\.tif$", full.names = TRUE)
c_files <- list.files(c_folder, pattern = "\\.tif$", full.names = TRUE)

# Extract year keys
r_years <- str_extract(basename(r_files), "\\d{4}")
c_years <- str_extract(basename(c_files), "\\d{4}")

# Intersect years available in both R and C
common_years <- intersect(r_years, c_years)

#NDVI for correlation
ndvi_folder <- "sample_data/NDVI"

prev_ndvi <- NA  # initialize outside the loop

ingsAsFactors = FALSE
)

# Loop through each year

rusle_results <- data.frame(
  Year = character(),
  Mean_A = numeric(),
  SD_A = numeric(),
  Mean_R = numeric(),
  Mean_K = numeric(),
  Mean_LS = numeric(),
  Mean_C = numeric(),
  Mean_NDVI = numeric(),
  Delta_NDVI = numeric(),
  stringsAsFactors = FALSE
)

for (year in sort(common_years)) {
  message("Processing year: ", year)
  
  r_path <- r_files[which(r_years == year)]
  c_path <- c_files[which(c_years == year)]
  
  R <- rast(r_path)
  C <- rast(c_path)
  
  # Locate matching NDVI file
  ndvi_path <- list.files(ndvi_folder, pattern = paste0(year, ".*\\.tif$"), full.names = TRUE)
  if (length(ndvi_path) == 1) {
    ndvi_raster <- rast(ndvi_path)
    ndvi_projected <- project(ndvi_raster, "EPSG:32610", method = "near")
    ndvi_resampled <- resample(mask(crop(ndvi_projected, test_site), test_site), reference)
    mean_ndvi <- mean(values(ndvi_resampled, na.rm = TRUE))
  } else {
    warning("NDVI file for year ", year, " not found or ambiguous.")
    mean_ndvi <- NA
  }
  
  # Compute NDVI change
  if (!is.na(prev_ndvi) && !is.na(mean_ndvi)) {
    delta_ndvi <- mean_ndvi - prev_ndvi
  } else {
    delta_ndvi <- NA
  }
  prev_ndvi <- mean_ndvi  # update for next year
  
  # Compute RUSLE A
  A <- R * K * LS * C
  output_path <- file.path(output_rusle_folder, paste0("RUSLE_A_", year, ".tif"))
  writeRaster(A, output_path, overwrite = TRUE)
  
  # Stats
  mean_A <- mean(values(A, na.rm = TRUE))
  sd_A <- sd(values(A, na.rm = TRUE))
  mean_R <- mean(values(R, na.rm = TRUE))
  mean_K <- mean(values(K, na.rm = TRUE))
  mean_LS <- mean(values(LS, na.rm = TRUE))
  mean_C <- mean(values(C, na.rm = TRUE))
  
  rusle_results <- bind_rows(rusle_results, data.frame(
    Year = year,
    Mean_A = mean_A,
    SD_A = sd_A,
    Mean_R = mean_R,
    Mean_K = mean_K,
    Mean_LS = mean_LS,
    Mean_C = mean_C,
    Mean_NDVI = mean_ndvi,
    Delta_NDVI = delta_ndvi
  ))
}


# Save results to CSV
write.csv(rusle_results, "rusle_summary_by_year.csv", row.names = FALSE)

###Calculate for each pixel for higher granularity

rusle_pixel_data <- data.frame()
prev_ndvi_vals <- NULL

for (year in sort(common_years)) {
  message("Processing year: ", year)
  
  # Get file paths for R, C, NDVI
  r_path <- r_files[which(r_years == year)]
  c_path <- c_files[which(c_years == year)]
  ndvi_path <- list.files(ndvi_folder, pattern = paste0(year, ".*\\.tif$"), full.names = TRUE)
  
  # Load rasters
  R <- rast(r_path)
  C <- rast(c_path)
  
  # Check raster alignment of R, K, LS, C before calculation
  # Assuming K and LS are already loaded and aligned (preloaded outside loop)
  # Align R to K
  if (!compareGeom(R, K, stopOnError = FALSE)) {
    R <- resample(R, K, method = "bilinear")
  }
  # Align C to K
  if (!compareGeom(C, K, stopOnError = FALSE)) {
    C <- resample(C, K, method = "bilinear")
  }
  
  # Calculate A raster
  A <- R * K * LS * C
  
  # Check that A is valid and non-empty
  A_vals <- as.vector(values(A))
  if (length(A_vals) == 0 || all(is.na(A_vals))) {
    warning(paste("No valid data in A raster for year", year, "- skipping this year"))
    next
  }
  
  # Load NDVI raster, project, mask, resample
  if (length(ndvi_path) == 1) {
    ndvi_raster <- rast(ndvi_path)
    # Reproject NDVI to match K
    if (!compareGeom(ndvi_raster, K, stopOnError = FALSE)) {
      ndvi_raster <- project(ndvi_raster, crs(K), method = "near")
    }
    # Crop and mask to study area if needed, here assumed 'test_site' is a polygon (SpatVector)
    ndvi_crop <- crop(ndvi_raster, ext(K))
    ndvi_mask <- mask(ndvi_crop, K)
    ndvi_resampled <- resample(ndvi_mask, K, method = "bilinear")
    
    ndvi_vals <- as.vector(values(ndvi_resampled))
  } else {
    warning("NDVI file for year ", year, " not found or ambiguous.")
    ndvi_vals <- rep(NA, ncell(K))
  }
  
  # Calculate delta NDVI compared to previous year
  if (!is.null(prev_ndvi_vals)) {
    # Make sure lengths match
    if (length(ndvi_vals) == length(prev_ndvi_vals)) {
      delta_ndvi_vals <- ndvi_vals - prev_ndvi_vals
    } else {
      warning(paste("NDVI length mismatch for year", year))
      delta_ndvi_vals <- rep(NA, length(ndvi_vals))
    }
  } else {
    delta_ndvi_vals <- rep(NA, length(ndvi_vals))
  }
  
  # Extract values for other rasters
  R_vals <- as.vector(values(R))
  K_vals <- as.vector(values(K))
  LS_vals <- as.vector(values(LS))
  C_vals <- as.vector(values(C))
  # A_vals already extracted
  
  # Confirm all vectors have the same length
  val_lengths <- c(length(R_vals), length(K_vals), length(LS_vals), length(C_vals), length(A_vals), length(ndvi_vals), length(delta_ndvi_vals))
  if (length(unique(val_lengths)) != 1) {
    warning(paste("Raster value lengths mismatch for year", year, "- skipping"))
    next
  }
  
  # Build data frame for this year
  df <- data.frame(
    Year = year,
    cell = 1:length(A_vals),
    R = R_vals,
    K = K_vals,
    LS = LS_vals,
    C = C_vals,
    NDVI = ndvi_vals,
    Delta_NDVI = delta_ndvi_vals,
    A = A_vals
  )
  
  # Remove rows with NA in A (outside study area etc.)
  df <- df %>% filter(!is.na(A))
  
  # Append to combined data frame
  rusle_pixel_data <- bind_rows(rusle_pixel_data, df)
  
  # Store current NDVI values for next year delta calculation
  prev_ndvi_vals <- ndvi_vals
}

# Save or return rusle_pixel_data as needed
write.csv(rusle_pixel_data, "rusle_summary_by_year_pixelwise.csv", row.names = FALSE)
