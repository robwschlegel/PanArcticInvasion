# code/01_AIS_prep.R
# This script is designed to prepare/process ship AIS data from a local drive 
# into a tidy data.frame that can be used for further investigations
# of the movement of ships into Arctic ports


# Setup -------------------------------------------------------------------

# Load packages+functions
source("code/00_functions.R")

# The movement files
movement_files <- dir("data", pattern = "Movement*.*txt", full.names = T)


# Data --------------------------------------------------------------------

# Open first few lines to inspect file
AIS_2015 <- read_delim("data/MovementData_2015H1.txt", n_max = 20, delim = "\t")
AIS_2017 <- read_delim("data/CombinedHistory_2017H1.txt", n_max = 20, delim = "\t")

# Load all data for on file
AIS_2015 <- read_delim("data/MovementData_2015H1.txt", delim = "\t")

# Look at list of unique zone names (ports)
unique(AIS_2015$ZoneName)
