# code/00_functions.R
# Functions used throughout the project


# Setup -------------------------------------------------------------------

library(tidyverse)
library(geosphere)
library(arrow)
library(doParallel); registerDoParallel(cores = 15)


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

