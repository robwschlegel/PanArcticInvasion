# code/02_BioOracle_prep.R
# The code used to download and prepare the SDM data layers from Bio-Oracle v3.0


# Setup -------------------------------------------------------------------

# Load packages+functions
source("code/00_functions.R")

# Bio-Oracle access
library(biooracler)

# Set cores
registerDoParallel(cores = 15)

# Disable scientific notation
options(scipen = 999)

# Data requirements
## Use all mean variables
  ## We do not want to start looking at differences between max, min, or long term max or min
## SST, SSS, and Ice thickness for all types of taxonomic groups. 
## Then add different variables according to each group:
## Fish: chlorophyll a and nitrate
## Phytopbenthos and Phytoplankton: nitrate, iron, and PAR
## Zoobenthos and Zooplankton: Chlorophyll a


# Download and prep -------------------------------------------------------

# Mean surface values for: 
# temperature, salinity, sea ice thickness, chlorophyll a, nitrate, iron, and PAR present and future
var_choices <- c("thetao", "so", "sithick", "chl", "no3", "dfe", "par")


## Metadata ---------------------------------------------------------------

# Explore layers in a dataset
BO_layers <- biooracler::list_layers(simplify = TRUE)

# Info
info_layer("thetao_baseline_2000_2019_depthsurf")

# Future scenario conversions
# https://ar5-syr.ipcc.ch/topic_futurechanges.php



## Download all -----------------------------------------------------------

# NB: Only need to run this once

# Baseline 2020
plyr::l_ply(var_choices, BO_v3_dl, .parallel = TRUE)

# RCP 4.5
plyr::l_ply(var_choices, BO_v3_dl, .parallel = TRUE, scenario = "ssp245", decade = 2040)
plyr::l_ply(var_choices, BO_v3_dl, .parallel = TRUE, scenario = "ssp245", decade = 2090)

# RCP 8.5
plyr::l_ply(var_choices, BO_v3_dl, .parallel = TRUE, scenario = "ssp585", decade = 2040)
plyr::l_ply(var_choices, BO_v3_dl, .parallel = TRUE, scenario = "ssp585", decade = 2090)


## Present data -----------------------------------------------------------

# Baseline 2020 files
dir_baseline_2020 <- dir("data/BO_v3", pattern = "baseline_2020.tiff$", full.names = TRUE)

# Load and stack
BO3_baseline_2020 <- terra::rast(dir_baseline_2020)
names(BO3_baseline_2020)

# Plot
plot(BO3_baseline_2020)

# Save
BO_present <- writeRaster(BO3_baseline_2020, file = "data/BO_present.tiff", overwrite = TRUE)


## Future data ------------------------------------------------------------

# RCP8.5 2050 and 2100 files
dir_ssp585_2050 <- c(dir("data/BO_v3", pattern = "ssp585_2050.tiff$", full.names = TRUE), 
                     dir_baseline_2020[4]) # NB: Ensure this is adding the PAR data
dir_ssp585_2100 <- c(dir("data/BO_v3", pattern = "ssp585_2100.tiff$", full.names = TRUE), 
                     dir_baseline_2020[4]) # NB: Ensure this is adding the PAR data

# Load and stack
BO3_ssp585_2050 <- terra::rast(dir_ssp585_2050); names(BO3_ssp585_2050)
BO3_ssp585_2100 <- terra::rast(dir_ssp585_2100); names(BO3_ssp585_2100)

# Plot
plot(BO3_ssp585_2050)
plot(BO3_ssp585_2100)
plot(c(BO3_baseline_2020[[7]], BO3_ssp585_2050[[6]], BO3_ssp585_2100[[6]]))

# Save
BO_2050 <- writeRaster(BO3_ssp585_2050, file = "data/BO_2050.tiff", overwrite = TRUE)
BO_2100 <- writeRaster(BO3_ssp585_2100, file = "data/BO_2100.tiff", overwrite = TRUE)


# Test visuals ------------------------------------------------------------


