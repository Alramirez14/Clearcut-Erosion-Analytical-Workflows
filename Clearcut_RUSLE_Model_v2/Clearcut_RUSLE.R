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


# Initialize results data frame
rusle_results <- data.frame(
  Year = character(),
  Mean_A = numeric(),
  SD_A = numeric(),
  Mean_R = numeric(),
  Mean_K = numeric(),
  Mean_LS = numeric(),
  Mean_C = numeric(),
  stringsAsFactors = FALSE
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


