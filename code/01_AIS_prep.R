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
system.time(
arctic_previous_ports <- plyr::ldply(unique(arctic_ships$ShipName)[1:20], previous_ports, .parallel = T)
) # 38 seconds for first 20 with multi-core, 54 without
write_csv_arrow(arctic_previous_ports, file = "data/arctic_previous_ports.csv")

