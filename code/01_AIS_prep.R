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
#   extract(geometry, into = c('lon', 'lat'), '\\((.*),(.*)\\)', conv = T)
# write_csv(arctic_ports, file = "metadata/arctic_ports.csv")
arctic_ports <- read_csv("metadata/arctic_ports.csv")


# Data --------------------------------------------------------------------



# Open first few lines to inspect file
AIS_2015 <- read_delim("data/MovementData_2015H1.txt", n_max = 20, delim = "\t")
AIS_2017 <- read_delim("data/CombinedHistory_2017H1.txt", n_max = 20, delim = "\t")

# Load all data for on file
AIS_2015 <- read_delim("data/MovementData_2015H1.txt", delim = "\t")

# Look at list of unique zone names (ports)
unique(AIS_2015$ZoneName)
