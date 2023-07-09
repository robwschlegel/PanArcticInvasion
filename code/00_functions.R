# code/00_functions.R
# Functions used throughout the project


# Setup -------------------------------------------------------------------

library(tidyverse)
library(arrow)
library(doParallel); registerDoParallel(cores = 15)


# AIS data ----------------------------------------------------------------

# Function for loading a file to extract unique port info
extract_ports <- function(file_name){
  port_df <- arrow::read_delim_arrow("data/MovementData_2015H1.txt", delim = "\t") |> 
    dplyr::select(Country:AnchorageParentName) |> 
    distinct()
  return(port_df)
}

# Determine ships that arrived in a port and the previous two ports
