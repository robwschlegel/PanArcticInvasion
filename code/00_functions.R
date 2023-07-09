# code/00_functions.R
# Functions used throughout the project


# Setup -------------------------------------------------------------------

library(tidyverse)
library(arrow)


# AIS data ----------------------------------------------------------------

# Function for loading a file to extract unique port info
extract_ports <- function(file_name){
  port_df <- arrow::read_delim_arrow("data/MovementData_2015H1.txt", delim = "\t") |> 
    dplyr::select(Country:AnchorageParentName)
}
