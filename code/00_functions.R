# code/00_functions.R
# Functions used throughout the project


# Setup -------------------------------------------------------------------

library(tidyverse)
library(geosphere)
library(raster)
library(terra)
library(arrow)
library(doParallel); registerDoParallel(cores = detectCores()-1)


# BO v3 -------------------------------------------------------------------

# var_choice <- "thetao"; scenario <- "ssp245"; decade <- 2040
BO_v3_dl <- function(var_choice, scenario = "baseline", decade = 2010){
  
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
  
  # Create variable names
  if(var_choice == "par"){
    var_name <- "PAR_mean_mean"
  } else {
    var_name <- paste0(var_choice,"_mean")
  }
  
  # Set decade
  # NB: this is necessary for below
  decade_chr <- as.character(decade+10)
  file_name <- paste(var_name, scenario, decade_chr, sep = "_")
  
  # Set constraints
  # NB: Time must be two values, even if the same
  time_con <-  rep(paste0(decade,"-01-01T00:00:00Z"), 2)
  lat_con <- c(-89.975, 89.975); lon_con <- c(-179.975, 179.975) 
  set_con <- list(time_con, lat_con, lon_con); names(set_con) = c("time", "latitude", "longitude")
  
  # Get data and exit
  # NB: Future projections are too beefy to download 2050 and 2100 at the same time
  choice_layers <- download_layers(dataset_id = dataset_choice$dataset_id,
                                   variables = var_name, constraints = set_con)
  writeRaster(choice_layers, file = paste0("data/BO_v3/",file_name,".tiff"), overwrite = TRUE)
  # rm(var_choice, scenario, dataset_choice, layer_names, time_con, lat_con, lon_con, set_con, set_con_1, set_con_2, var_name,
     # choice_layers, choice_layers_1, choice_layers_2, decade, decade_chr, file_name)
}


# AIS data ----------------------------------------------------------------

# Function for loading a file to extract unique port info
extract_ports <- function(file_name){
  port_df <- arrow::read_delim_arrow(file_name, delim = "\t") |> 
    dplyr::select(Country:AnchorageParentName) |> 
    distinct()
  return(port_df)
}

# Receive a data.frame and row index to process last two ports
previous_ports_wide <- function(row_idx, df){
  
  # process data
  ship_wide <- df[c(row_idx, row_idx+1, row_idx+2),] |> 
    mutate(port_idx = 1:n()) |> 
    pivot_longer(c(-ShipName, -port_idx), names_to = "key", values_to = "val",
                 values_transform = list(val = as.character)) |> 
    pivot_wider(values_from = val, names_from = c(key, port_idx)) |> 
    filter(!is.na(ShipName)) |> # Catch edge cases when last port is at end of time series
    mutate(across(contains(c("CallID", "lon", "lat", "Draught", "Hours", "Arctic")), as.numeric)) |> 
    mutate(across(contains("Date"), as.POSIXct)) |> 
    mutate(Arctic_count = sum(c(Arctic_1, Arctic_2, Arctic_3), na.rm = T),
           dist_2 = case_when(!is.na(lon_2) ~ distm(c(lon_1, lat_1), c(lon_2, lat_2), 
                                                    fun = distHaversine)[1]/1000, TRUE ~ as.numeric(NA)),
           dist_3 = case_when(!is.na(lon_3) ~ distm(c(lon_2, lat_2), c(lon_3, lat_3), 
                                                    fun = distHaversine)[1]/1000, TRUE ~ as.numeric(NA)),
           time_2 = ArrivalDateFull_1 - SailDateFull_2, time_3 = ArrivalDateFull_2 - SailDateFull_3) |> 
    dplyr::select(ShipName, Arctic_count, 
                  CallID_1, ZoneName_1, Arctic_1, ArrivalDateFull_1, HoursinPort_1, 
                  SailDateFull_1, ArrivalDraught_1, DepartureDraught_1,
                  CallID_2, ZoneName_2, Arctic_2, dist_2, time_2, ArrivalDateFull_2, HoursinPort_2, 
                  SailDateFull_2, ArrivalDraught_2, DepartureDraught_2,
                  CallID_3, ZoneName_3, Arctic_3, dist_3, time_3, ArrivalDateFull_3, HoursinPort_3, 
                  SailDateFull_3, ArrivalDraught_3, DepartureDraught_3)
  return(ship_wide)
  # rm(row_idx, df, ship_sub, ship_wide)
}

# Determine ships that arrived in a port and the previous two ports
previous_ports <- function(ship_name){
  
  # Filter one ship of data
  ship_df <- all_AIS |> 
    filter(ShipName == ship_name) |> arrange(desc(ArrivalDateFull)) |> 
    left_join(arctic_ports, by = join_by(Country, PortID, PortName, PortGeoID, ZoneName)) |>
    mutate(Arctic = case_when(REALM == "Arctic" ~ 1, TRUE ~ 0)) |> 
    dplyr::select(CallID, ShipName, ZoneName, PortLatitudeDecimal, PortLongitudeDecimal, 
                  ArrivalDateFull, SailDateFull, ArrivalDraught, DepartureDraught, HoursinPort, Arctic) |> 
    dplyr::rename(lon = PortLongitudeDecimal, lat = PortLatitudeDecimal)
  
  # Get row index of when the ship docked at an Arctic port
  ship_idx <- which(ship_df$Arctic == TRUE)
  
  # Extract and process data for last two ports
  # system.time(
  ship_res <- plyr::ldply(ship_idx, previous_ports_wide, .parallel = F, df = ship_df)
  # ) # 4 seconds single core, 2 seconds multi-core
  
  # Exit
  return(ship_res)
  # rm(ship_name, ship_df, ship_idx, ship_res)
}


# Various -----------------------------------------------------------------

# Function for re-loading .RData files as necessary
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Convenience function for loading .asc files
load_asc <- function(file_name, col_name){
  df <- as.data.frame(raster(file_name), xy = T) |> 
    `colnames<-`(c("lon", "lat", col_name))
}

# Quickly convert biomod2 raster to a data.frame
rast_df <- function(rast){
  df_out <- as.data.frame(rast[[1]], xy = T) |> 
    `colnames<-`(c("lon", "lat", "presence")) |>  
    # mutate(lon = round(lon, 4), lat = round(lat, 4)) %>% 
    # left_join(bathy_layer, by = c("lon", "lat")) %>% 
    na.omit() 
}

