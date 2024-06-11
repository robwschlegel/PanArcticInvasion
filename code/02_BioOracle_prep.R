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

# Present baseline tests
plyr::l_ply(var_choices, BO_v3_dl, .parallel = TRUE)

# RCP 4.5
plyr::l_ply(var_choices, BO_v3_dl, .parallel = TRUE, scenario = "ssp245")

# RCP 8.5
plyr::l_ply(var_choices, BO_v3_dl, .parallel = TRUE, scenario = "ssp585")


## Present data -----------------------------------------------------------


var_choices <- c("thetao", "so", "sithick", "chl", "no3", "dfe", "par")


# Load and stack them
BO3_thetao_baseline <- read_rds("data/BO_V3/thetao_baseline.rds")
BO3_so_baseline <- read_rds("data/BO_V3/so_baseline.rds")
BO3_sithick_baseline <- loadRData("data/BO_V3/sithick_baseline.RData")
BO3_chl_baseline <- loadRData("data/BO_V3/chl_baseline.RData")
BO3_no3_baseline <- loadRData("data/BO_V3/no3_baseline.RData")
BO3_dfe_baseline <- loadRData("data/BO_V3/dfe_baseline.RData")
BO3_par_baseline <- loadRData("data/BO_V3/par_baseline.RData")
BO3_baseline <- raster::stack(BO3_thetao_baseline, BO3_so_baseline)
names(BO3_baseline)

# Convert to dataframe
BO_present <- as.data.frame(BO_layers_dl, xy = T) %>% 
  dplyr::rename(lon = x, lat = y) %>% 
  filter(BO22_icethickmean_ss >= 0, lat >= 30)
save(BO_present, file = "data/BO_present.RData")
rm(BO_layers_dl); gc()

# Visualise
ggplot(BO_present, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = BO22_parmean)) +
  borders(fill = "grey70", colour = "black") +
  scale_fill_viridis_c(option = "D") +
  coord_quickmap(xlim = c(-180, 180), ylim = c(25, 90), expand = F) +
  theme(legend.position = "bottom")


## Future data ------------------------------------------------------------

# Load present data for easier joining
if(!exists("BO_present")) load("data/BO_present.RData")

# Download as similar of layers as possible to present data
                                           # SST, SSS, Ice thickness
BO_layers_future <- get_future_layers(c("BO22_tempmean_ss", "BO22_salinitymean_ss", "BO22_icethickmean_ss",
                                        # chlorophyll a 
                                        # NB: no futre layers for nitrates, iron, PAR
                                        "BO22_chlomean_ss"),
                                        # NB: at the moment just getting surface values, but may want bottom, too
                                        # "BO22_chlomean_bdmax"), 
                                      scenario = "RCP85", year = c(2050, 2100))

# NB: This will download the data, but throws an error when stacking them due to different chl a extent
BO_layers_future_dl <- load_layers(BO_layers_future$layer_code)
BO_future_chlo <- as.data.frame(load_layers(BO_layers_future$layer_code[c(1,5)]), xy = T)
BO_future_other <- as.data.frame(load_layers(BO_layers_future$layer_code[c(2,3,4,6,7,8)]), xy = T)

# Convert to data.frame
BO_future <- left_join(BO_future_other, BO_future_chlo, by = c("x", "y")) |> 
  dplyr::rename(lon = x, lat = y) |> 
  filter(BO22_RCP85_2100_icethickmean_ss >= 0, lat >= 30)
save(BO_future, file = "data/BO_future.RData")
rm(BO_layers_future, BO_layers_future_dl, BO_future_chlo, BO_future_other); gc()

# Create 2050 data.frame
BO_2050 <- BO_present %>% 
  dplyr::select(lon, lat, BO22_nitratemean_ss:BO22_parmean) %>% 
  left_join(dplyr::select(BO_future, lon, lat, BO22_RCP85_2050_icethickmean_ss:BO22_RCP85_2050_tempmean_ss,
                          BO22_RCP85_2050_chlomean_ss), by = c("lon", "lat"))
colnames(BO_2050) <- sub("RCP85_2050_", "", colnames(BO_2050))
save(BO_2050, file = "data/BO_2050.RData")

# Create 2100 data.frame
BO_2100 <- BO_present %>% 
  dplyr::select(lon, lat, BO22_nitratemean_ss:BO22_parmean) %>% 
  left_join(dplyr::select(BO_future, lon, lat, BO22_RCP85_2100_icethickmean_ss:BO22_RCP85_2100_tempmean_ss,
                          BO22_RCP85_2100_chlomean_ss), by = c("lon", "lat"))
colnames(BO_2100) <- sub("RCP85_2100_", "", colnames(BO_2100))
save(BO_2100, file = "data/BO_2100.RData")

# Visualise
ggplot(BO_2100, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = BO22_tempmean_ss)) +
  borders(fill = "grey70", colour = "black") +
  scale_fill_viridis_c(option = "D") +
  coord_quickmap(xlim = c(-180, 180), ylim = c(25, 90), expand = F) +
  theme(legend.position = "bottom")


# Test visuals ------------------------------------------------------------


