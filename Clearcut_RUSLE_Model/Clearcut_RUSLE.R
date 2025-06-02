###############################################################
######### This script runs the RUSLE model for each year   ####
######### within each setup study site                     ####
######### Created 05/25/2025 by Clearcut Erosion Modeling  ####
###############################################################


library(terra)
library(dplyr)
library(stringr)

# Set folders
r_folder <- "processed_R_factors"      # Folder with R rasters (one per year)
c_folder <- "processed_C_factors"      # Folder with C rasters (one per year)
output_rusle_folder <- "rusle_outputs" # Folder for RUSLE rasters (optional)
dir.create(output_rusle_folder, showWarnings = FALSE)

# Load static rasters
K <- rast("path/to/static_K_resampled.tif")
LS <- rast("path/to/static_LS.tif")

# List R and C files
r_files <- list.files(r_folder, pattern = "\\.tif$", full.names = TRUE)
c_files <- list.files(c_folder, pattern = "\\.tif$", full.names = TRUE)

# Extract year keys
r_years <- str_extract(basename(r_files), "\\d{4}")
c_years <- str_extract(basename(c_files), "\\d{4}")

# Intersect years available in both R and C
common_years <- intersect(r_years, c_years)

# Initialize results data frame
rusle_results <- data.frame(Year = character(), Mean_A = numeric(), SD_A = numeric(), stringsAsFactors = FALSE)

# Loop through each year
for (year in sort(common_years)) {
  message("Processing year: ", year)
  
  r_path <- r_files[which(r_years == year)]
  c_path <- c_files[which(c_years == year)]
  
  R <- rast(r_path)
  C <- rast(c_path)
  
  # Compute RUSLE A = R * K * LS * C
  A <- R * K * LS * C
  
  # Optionally save the raster
  output_path <- file.path(output_rusle_folder, paste0("RUSLE_A_", year, ".tif"))
  writeRaster(A, output_path, overwrite = TRUE)
  
  # Summarize and append to results
  stats <- values(A, na.rm = TRUE)
  mean_A <- mean(stats)
  sd_A <- sd(stats)
  
  rusle_results <- bind_rows(rusle_results, data.frame(Year = year, Mean_A = mean_A, SD_A = sd_A))
}

# Save results to CSV
write.csv(rusle_results, "rusle_summary_by_year.csv", row.names = FALSE)
