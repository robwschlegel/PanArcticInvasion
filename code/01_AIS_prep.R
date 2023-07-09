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
# NB: It appears to be the same if one runs this on all files, or just one
# i.e. any given file appears to contain records for all ports
# all_ports <- plyr::ldply(movement_files, extract_ports, .parallel = T) |> distinct()
# write_csv(arctic_ports, file = "metadata/all_ports.csv")
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
arctic_ships <- all_AIS |> 
  filter(ZoneName %in% arctic_ports$ZoneName) |> 
  dplyr::select(LRNOIMOShipNo:ShipType) |> distinct()
write_csv_arrow(arctic_ships, file = "metadata/arctic_ships.csv")

# Previous two ports of call
AIS_1 <- AIS_all |> 
  filter(ZoneName == arctic_ports$ZoneName[1]) |> 
  arrange(-ArrivalDateFull)
AIS_2 <- AIS_all |> 
  filter(ShipName == "ANTIGUA")
AIS_3 <- AIS_2 |> 
  left_join(arctic_ports, by = join_by(Country, PortID, PortName, PortGeoID, ZoneName))
AIS_4 <- which(AIS_3$REALM == "Arctic")
AIS_5 <- AIS_3[c(AIS_4[1], AIS_4[1]+1, AIS_4[1]+2),] |> 
  pivot_longer()

