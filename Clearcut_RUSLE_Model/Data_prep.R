###############################################################
######### This script accomplishes requisite data cleaning ####
######### to run the RUSLE model                           ####
######### Created 05/23/2025 by Clearcut Erosion Modeling  ####
###############################################################

###Setup Environment

require(raster)
require(terra)
require(dplyr)

###Load Study Site

test_site_path <- ("your/test/site/here.shp")

test_site <- vect(test_site_path) %>% 
  project("EPSG:32610")

###Load DEM

dem <- rast("Your/DEM.tif") %>% 
  project("EPSG:32610", method = "near")

###Process Precipitation Data


#---! NEED TO IMPLEMENT LOOP FOR ANNUAL VALUES !---#

precip_raw <- rast("Your.precip.bil")

precip_projected <- project(precip_raw, "EPSG:32610", method = "near")

writeRaster(precip_utm, "precip_projected.tif", overwrite = TRUE)

###Load Preclipped NDVI Data

NDVI_2023 <- rast("Your/NDVI/image.tif")

###K Factor 

K <- rast("Your/K/factor.tif") %>% 
  '/'(100) #correcting unit scale for RUSLE

K_projected <- project(K, "EPSG:32610", method = "near")

###Assign NDVI as reference layer (smallest size, highest res) 

#---! NEED TO IMPLEMENT LOOP FOR ANNUAL VALUES !---#

reference <- NDVI_2023

###Clip and mask early for handling

dem_clipped <- mask(crop(dem, test_site), test_site)

precip_clipped <- mask(crop(precip_projected, test_site), test_site)

K_clipped <- mask(crop(K_projected, test_site), test_site)

dem_clipped <- mask(crop(dem_resampled, test_site), test_site)

###Calculate slope on clipped DEM

slope_clipped <- dem_clipped %>% 
  project("EPSG:32610", method = "near") %>% 
  terrain(v = "slope", unit = "degrees")

###Resample for homogenous resolution/grid

#---! NEED TO IMPLEMENT LOOP FOR ANNUAL VALUES !---#

slope_resampled <- resample(slope_clipped, reference, method = "near")

precip_resampled <- resample(precip_clipped, reference, method = "near")

K_resampled <- resample(K_clipped, reference, method = "near")

dem_resampled <- resample(dem_clipped, reference, method = "near")

#NDVI is implicitly resampled as it is defined as the reference layer

###Calculate R Factor (work in progress)

#Utilizing Moore et al. method (high values)


#---! NEED TO IMPLEMENT LOOP FOR ANNUAL VALUES !---#

calculate_r_moore <- function (p) {
  ke <- 11.46*p - 2226
  r <- 0.029*ke - 26
  r_si <- 17.02*r
  return(r_si)
}

r_moore <- app(x = precip_resampled, fun = calculate_r_moore)

#Utilizing Wischmeier & Smith method (lower but still inflated)
calculate_R_wischmeier <- function(p_mm) {
  p_in <- p_mm / 25.4  
  R_in <- 0.7397 * p_in^1.847  # Wischmeier & Smith formula
  return(R_in)  # Units: MJ mm ha-1 h-1 yr-1
}

r_wischmeier <- app(x = precip_resampled, fun = calculate_R_wischmeier)

###Create LS and Flow accumulation rasters

flow_dir <- terrain(dem_resampled, 
                            v = "flowdir", 
                            unit = "radians")

flow_acc <- flowAccumulation(flow_dir)

cell_size <- res(dem_resampled)[1]  # get cell size in meters (assuming square cells)

slope_rad <- slope_clipped * (pi / 180)  # if slope in degrees; skip if already radians

m <- 0.5
n <- 1.3

LS <- ((flow_acc_clipped * cell_size) / 22.13)^m * (sin(slope_rad) / 0.0896)^n

### Calculate Cover Factor C


#---! NEED TO IMPLEMENT LOOP FOR ANNUAL VALUES !---#

calculate_c_knijff <- function(ndvi) {
  alpha <- 2
  beta <- 0.5  # smaller than 1 to make the curve responsive
  c <- exp(-alpha * (ndvi / beta - ndvi))
  c <- ifelse(c > 1, 1, c)
  return(c)
}


c_factor <- app(NDVI_2023, calculate_c_knijff)

c_resampled <- resample(c_factor, reference, method = "near")

C <- mask(crop(c_resampled, test_site), test_site)

###Write all of these rasters to a folder directory

