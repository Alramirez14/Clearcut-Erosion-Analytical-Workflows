###############################################################
######### This script accomplishes requisite data cleaning ####
######### to run the RUSLE model                           ####
######### Created 05/23/2025 by Clearcut Erosion Modeling  ####
###############################################################

###Setup Environment

require(raster)
require(terra)
require(dplyr)
require(stringr)

###Load Study Site

test_site_path <- ("sample_data/sample_site/Smith_study_polygons.shp")

test_site <- vect(test_site_path) %>% 
  project("EPSG:32610")

###Establish reference image with one NDVI file

first_NDVI_image_path <- ("sample_data/NDVI/Smith_Landsat8_NDVI_2013.tif")

reference <- rast(first_NDVI_image_path)

###Load DEM

dem <- rast("sample_data/DEM/DEM_UTM10N_clipped.tiff") %>% 
  project("EPSG:32610", method = "near")

### Calculate Cover Factor C

# Define folders
ndvi_folder <- "sample_data/NDVI"
output_c_folder <- "processed_C_factors"
dir.create(output_c_folder, showWarnings = FALSE)

# Load all NDVI files
ndvi_files <- list.files(ndvi_folder, pattern = "\\.tif$", full.names = TRUE)

# Define the C factor calculation function (Knijff et al.)
calculate_c_knijff <- function(ndvi) {
  alpha <- 2
  beta <- 0.5
  c <- exp(-alpha * (ndvi / beta - ndvi))
  c <- ifelse(c > 1, 1, c)
  return(c)
}

# Loop through NDVI files
for (file in ndvi_files) {
  year <- str_extract(basename(file), "\\d{4}")
  if (is.na(year)) next
  
  ndvi <- rast(file)
  ndvi_projected <- project(ndvi, "EPSG:32610", method = "near")
  
  # Calculate C factor
  c_factor <- app(ndvi_projected, calculate_c_knijff)
  
  # Clip and mask to study site
  c_clipped <- mask(crop(c_factor, test_site), test_site)
  
  # Resample to match reference raster
  c_resampled <- resample(c_clipped, reference, method = "near")
  
  # Save output
  writeRaster(c_resampled, file.path(output_c_folder, paste0("C_factor_", year, ".tif")), overwrite = TRUE)
}

###Process Precipitation Data

# Set the folder path where annual precipitation rasters are stored

precip_folder <- "sample_data/precip_raw"  # e.g., "data/precip/"
precip_files <- list.files(precip_folder, pattern = "\\.bil$", full.names = TRUE)

# Define the destination folder for processed precipitation rasters

output_precip_folder <- "processed_precip"
dir.create(output_precip_folder, showWarnings = FALSE)

# Create an empty list to store processed rasters
precip_list <- list()

# Loop through each file, assuming filenames include the year (e.g., "precip_2020.bil")

for (file in precip_files) {
  
  # Extract year from filename using regex
  year <- stringr::str_extract(basename(file), "\\d{4}")
  if (is.na(year)) next  # skip files without a valid year
  
  # Load, project, and resample precipitation
  precip_raw <- rast(file)
  precip_projected <- project(precip_raw, "EPSG:32610", method = "near")
  precip_clipped <- mask(crop(precip_projected, test_site), test_site)
  precip_resampled <- resample(precip_clipped, reference, method = "near")
  
  # Save processed raster
  output_path <- file.path(output_precip_folder, paste0("precip_", year, "_processed.tif"))
  writeRaster(precip_resampled, output_path, overwrite = TRUE)
  
  # Optionally store in list for further calculations
  precip_list[[year]] <- precip_resampled
}


### !!!!!! NDVI Data should be prepped externally using the provided Google Earth 
### Engine Script

###K Factor 

K <- rast("sample_data/soil_data/k_utm_v2.tif") %>% 
  '/'(100) #correcting unit scale for RUSLE

K_projected <- project(K, "EPSG:32610", method = "near")

###Clip and mask static for handling

K_clipped <- mask(crop(K_projected, test_site), test_site)

dem_clipped <- mask(crop(dem, test_site), test_site)

###Calculate slope on clipped DEM

slope_clipped <- dem_clipped %>% 
  project("EPSG:32610", method = "near") %>% 
  terrain(v = "slope", unit = "degrees")

###Resample for homogeneous resolution/grid

#---! NEED TO IMPLEMENT LOOP FOR ANNUAL VALUES !---#

slope_resampled <- resample(slope_clipped, reference, method = "near")

K_resampled <- resample(K_clipped, reference, method = "near")

output_K_folder <- "processed_K_factor"
dir.create(output_K_folder, showWarnings = FALSE)

writeRaster(K_resampled, file.path(output_K_folder, "K.tif"), overwrite = TRUE)

dem_resampled <- resample(dem_clipped, reference, method = "near")

#NDVI is implicitly resampled as it is defined as the reference layer, and precip
#was resampled in the batch processing step

###Calculate R Factor (work in progress)

calculate_r_moore <- function(p) {
  ke <- 11.46 * p - 2226
  r <- 0.029 * ke - 26
  r_si <- 17.02 * r
  return(r_si)
}

calculate_R_wischmeier <- function(p_mm) {
  p_in <- p_mm / 25.4  
  R_in <- 0.7397 * p_in^1.847
  return(R_in)
}

output_rmoore_folder <- "processed_Rmoore_factors"
dir.create(output_rmoore_folder, showWarnings = FALSE)

output_rwisch_folder <- "processed_Rwisch_factors"
dir.create(output_rwisch_folder, showWarnings = FALSE)

processed_files <- list.files(output_precip_folder, pattern = "\\.tif$", full.names = TRUE)

for (file in processed_files) {
  year <- stringr::str_extract(basename(file), "\\d{4}")
  if (is.na(year)) next
  
  precip <- terra::rast(file)
  
  r_moore <- terra::app(precip, calculate_r_moore)
  r_wisch <- terra::app(precip, calculate_R_wischmeier)
  
  writeRaster(r_moore, file.path(output_rmoore_folder, paste0("R_moore_", year, ".tif")), overwrite = TRUE)
  writeRaster(r_wisch, file.path(output_rwisch_folder, paste0("R_wisch_", year, ".tif")), overwrite = TRUE)
}


###Create LS and Flow accumulation rasters

flow_dir <- terrain(dem_resampled, 
                            v = "flowdir", 
                            unit = "radians")

flow_acc <- flowAccumulation(flow_dir)

cell_size <- res(dem_resampled)[1]  # get cell size in meters (assuming square cells)

slope_rad <- slope_resampled * (pi / 180)  # if slope in degrees; skip if already radians

m <- 0.5
n <- 1.3

LS <- ((flow_acc * cell_size) / 22.13)^m * (sin(slope_rad) / 0.0896)^n

output_ls_folder <- "processed_LS_factor"
dir.create(output_ls_folder, showWarnings = FALSE)

writeRaster(LS, file.path(output_ls_folder, "LS.tif"), overwrite = TRUE) #tweak for desired output dir



