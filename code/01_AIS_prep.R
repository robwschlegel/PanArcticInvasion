# code/01_AIS_prep.R
# This script is designed to prepare/process ship AIS data from a local drive 
# into a tidy data.frame that can be used for further investigations
# of the movement of ships into Arctic ports


# Setup -------------------------------------------------------------------

# Load packages+functions
source("code/00_functions.R")
library(sf)

# The movement files
movement_files <- dir("data", pattern = "Movement*.*txt", full.names = T)

# Arctic MEOW
arctic_MEOW <- read_sf("metadata/MEOW/meow_ecos.shp") |> filter(REALM == "Arctic")


# Arctic ports ------------------------------------------------------------

# Extract full list of ports from data
# all_ports <- plyr::ldply(movement_files, extract_ports, .parallel = T) |> distinct()
# write_csv(all_ports, file = "metadata/all_ports.csv")
all_ports <- read_csv("metadata/all_ports.csv")

# Determine which ports are in an Arctic MEOW
# arctic_ports <- all_ports |>
#   st_as_sf(coords = c("PortLongitudeDecimal", "PortLatitudeDecimal"), crs = 4326) |>
#   st_join(arctic_MEOW) |> filter(!is.na(REALM)) |>
#   extract(geometry, into = c('lon', 'lat'), '\\((.*),(.*)\\)', conv = T) |>
#   dplyr::select(lon, lat, Country:ZoneName, ECOREGION, PROVINCE, REALM)
# write_csv(arctic_ports, file = "metadata/arctic_ports.csv")
arctic_ports <- read_csv("metadata/arctic_ports.csv")


# AIS data ----------------------------------------------------------------

# Load all AIS data
# all_AIS <- plyr::ldply(movement_files, read_delim_arrow, delim = "\t", .parallel = T)
# write_csv_arrow(all_AIS, file = "data/all_AIS.csv")
all_AIS <- read_csv_arrow("data/all_AIS.csv")

# Get ships that have been in the Arctic
# arctic_ships <- all_AIS |>
#   filter(ZoneName %in% arctic_ports$ZoneName) |>
#   dplyr::select(LRNOIMOShipNo:ShipType) |> distinct()
# write_csv_arrow(arctic_ships, file = "metadata/arctic_ships.csv")
arctic_ships <- read_csv_arrow("metadata/arctic_ships.csv")

# Previous two ports of call
ship_1 <- all_AIS |> 
  filter(ShipName == "ANTIGUA") |> arrange(desc(ArrivalDateFull)) |> 
  left_join(arctic_ports, by = join_by(Country, PortID, PortName, PortGeoID, ZoneName)) |>
  mutate(Arctic = case_when(REALM == "Arctic" ~ 1, TRUE ~ 0)) |> 
  dplyr::select(CallID, ShipName, ZoneName, PortLatitudeDecimal, PortLongitudeDecimal, 
                ArrivalDateFull, SailDateFull, ArrivalDraught, DepartureDraught, HoursinPort, Arctic) |> 
  dplyr::rename(lon = PortLongitudeDecimal, lat = PortLatitudeDecimal)
ship_2 <- which(ship_1$Arctic == TRUE)
ship_3 <- ship_1[c(ship_2[1], ship_2[1]+1, ship_2[1]+2),] |> 
  filter(!is.na(ShipName)) |> # Catch edge cases when last port is at end of time series
  mutate(port_idx = 1:n())
ship_4 <- ship_3 |> 
  pivot_longer(c(-ShipName, -port_idx), names_to = "key", values_to = "val", 
               values_transform = list(val = as.character)) |> 
  pivot_wider(values_from = val, names_from = c(key, port_idx)) |> 
  mutate(across(contains(c("CallID", "lon", "lat", "Draught", "Hours", "Arctic")), as.numeric)) |> 
  mutate(across(contains("Date"), as.POSIXct)) |> 
  mutate(Arctic_sum = Arctic_1+ Arctic_2 + Arctic_3,
         dist_1 = 0,
         dist_2 = distm(c(lon_1, lat_1), c(lon_2, lat_2), fun = distHaversine)[1]/1000,
         dist_3 = distm(c(lon_2, lat_2), c(lon_3, lat_3), fun = distHaversine)[1]/1000,
         time_1 = 0,
         time_2 = ArrivalDateFull_1 - SailDateFull_2,
         time_3 = ArrivalDateFull_2 - SailDateFull_3)

