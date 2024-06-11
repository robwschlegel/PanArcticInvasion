# tests/test_BO.R
# The purpose of this script is to test the v3.0 Bio-Oracle data
# More tests available here:
# https://github.com/robwschlegel/ArcticKelp/blob/master/tests/test_BO.R


# Setup -------------------------------------------------------------------

# Load tidyverse
library(tidyverse)

# Bio-Oracle access
library(sdmpredictors)
library(biooracler) # New version

# Work with raster layers
library(raster)

# Set cores
library(doParallel)
registerDoParallel(cores = 15)

# Disable scientific notation
options(scipen = 999)

# Check layers
BO_layers <- biooracler::list_layers(simplify = TRUE)

# Set temporary directory for testing
dir <- tempdir()

# Info
info_layer("thetao_baseline_2000_2019_depthsurf")


# Functions ---------------------------------------------------------------

# var_choice <- "thetao"; scenario <- "baseline"
BO_v3_test <- function(var_choice, scenario = "baseline"){
  
  # Get info if necessary
  if(!exists("BO_layers")) BO_layers <- biooracler::list_layers(simplify = TRUE)
  
  # Extract dataset id
  # NB: This more complicated method is used because not all baseline periods have
  # the same year range. E.g. 2000-2018 OR 2000-2019
  dataset_choice <- BO_layers |> 
    filter(grepl(paste0(var_choice,"_"), dataset_id)) |> 
    filter(grepl(paste0(scenario,"_"), dataset_id)) |> 
    filter(grepl("depthsurf", dataset_id)) |> 
    filter(!grepl("kdpar", dataset_id))
  
  # Exit if no values
  if(nrow(dataset_choice) == 0) return()
  
  # Set constraints
  if(scenario == "baseline") {
    time_con <-  c("2001-01-01T00:00:00Z", "2010-01-01T00:00:00Z")
  } else {
    time_con <-  c("2020-01-01T00:00:00Z", "2090-01-01T00:00:00Z")
  }
  lat_con <- c(-20, 89.975); lon_con <- c(-90, 90)
  set_con <- list(time_con, lat_con, lon_con); names(set_con) = c("time", "latitude", "longitude")
  
  # Create variable names
  if(var_choice == "par"){
    var_names <- c("PAR_mean_ltmax", "PAR_mean_ltmin")
    var_name <- "PAR_mean"
  } else {
    var_names <- c(paste0(var_choice,"_ltmax"), paste0(var_choice,"_ltmin"))
    var_name <- var_choice
  }
  
  # Test bottom currents
  choice_layers <- download_layers(dataset_id = dataset_choice$dataset_id,
                                   variables = var_names, constraints = set_con)
  choice_df <- as.data.frame(choice_layers, xy = T) |> 
    rename_at(.vars = vars(starts_with(var_name)),
              .funs = list(~ sub(paste0(var_name,"_"), "", .))) |> 
    dplyr::rename(lon = x, lat = y) |> na.omit() |> 
    mutate(lon = round(lon, 4), lat = round(lat, 4)) |>
    mutate(max_min_1 = ifelse(ltmax_1 >= ltmin_1, TRUE, FALSE),
           max_min_2 = ifelse(ltmax_2 >= ltmin_2, TRUE, FALSE)) |> 
    dplyr::select(-ltmax_1, -ltmax_2, -ltmin_1, -ltmin_2) |> 
    pivot_longer(max_min_1:max_min_2)
  
  # Visualise pixels where the max and min values are not as expected
  plot_choice <- ggplot(data = choice_df, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = value)) + coord_quickmap(expand = F) +
    labs(fill = "Max greater than min", x = NULL, y = NULL,
         title = paste0("Bio-Oracle v3.0: ", dataset_choice$dataset_id)) +
    facet_wrap("name") +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          legend.position = "bottom")
  ggsave(plot = plot_choice, filename = paste0("tests/",dataset_choice$dataset_id,".png"), height = 5, width = 9)
  # rm(var_choice, scenario, dataset_choice, var_names, var_name, time_con, lat_con, lon_con, set_con, choice_layers, choice_df, plot_choice); gc()
}


# Test key variables ------------------------------------------------------

# NB: In an earlier version of this script (see history via GitHub) all of the different data
# types were also tested, e.g. raster, NetCDF, and CSV
# There were some differences between the different file types, 
# but we only test the raster data type in this version of the script

# Mean surface values for: 
# temperature, salinity, sea ice thickness, chlorophyll a, nitrate, iron, and PAR present and future
var_choices <- c("thetao", "so", "sithick", "chl", "no3", "dfe", "par")

# Present baseline tests
plyr::l_ply(var_choices, BO_v3_test, .parallel = TRUE)

# RCP 4.5
plyr::l_ply(var_choices, BO_v3_test, .parallel = TRUE, scenario = "ssp245")

# RCP 8.5
plyr::l_ply(var_choices, BO_v3_test, .parallel = TRUE, scenario = "ssp585")

