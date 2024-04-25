# tests/test_BO.R
# The purpose of this script is to test the Bio-Oracle data
# More tests available here:
# https://github.com/robwschlegel/ArcticKelp/blob/master/tests/test_BO.R

# Load tidyverse
library(tidyverse)

# Bio-Oracle access
library(sdmpredictors)
library(biooracler) # New version

# Work with raster layers
library(raster)

# Set cores
registerDoParallel(cores = 15)

# Disable scientific notation
options(scipen = 999)

# Check layers
BO_layers <- biooracler::list_layers(simplify = FALSE)

# Set global constraints in advance
time_con <-  c("2001-01-01T00:00:00Z", "2010-01-01T00:00:00Z")
lat_con <- c(-89.975, 89.975)
lon_con <- c(-179.975, 179.975)
set_con <- list(time_con, lat_con, lon_con)
names(set_con) = c("time", "latitude", "longitude")

# Set temporary directory for testing
dir <- tempdir()


# Tests of bottom current layers ------------------------------------------

# Info
info_layer("sws_baseline_2000_2019_depthmax")

# Test bottom currents
current_bdmax_layers <- download_layers(dataset_id = "sws_baseline_2000_2019_depthmax", 
                                        variables = c("sws_ltmax", "sws_ltmin"), 
                                        constraints = set_con)
current_bdmax_test <- as.data.frame(current_bdmax_layers, xy = T) |> 
  dplyr::rename(lon = x, lat = y) |> na.omit() |> 
  mutate(lon = round(lon, 4), lat = round(lat, 4)) |>
  mutate(max_min_rast_2000 = ifelse(sws_ltmax_1 > sws_ltmin_1, TRUE, FALSE),
         max_min_rast_2010 = ifelse(sws_ltmax_2 > sws_ltmin_2, TRUE, FALSE)) |> 
  dplyr::select(-sws_ltmax_1, -sws_ltmax_2, -sws_ltmin_1, -sws_ltmin_2)

# Visualise pixels where the max and min values are not as expected
plot_sws_depthmax_2000 <- ggplot(data = current_bdmax_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min_rast_2000)) + coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Bottom currents v3.0; rast_2000: max vs min") +
  theme(legend.position = "bottom")
ggsave(plot = plot_sws_depthmax_2000, filename = "tests/sws_depthmax_rast_2000.png", height = 5, width = 8)
plot_sws_depthmax_2010 <- ggplot(data = current_bdmax_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min_rast_2010)) + coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Bottom currents v3.0; rast_2010: max vs min") +
  theme(legend.position = "bottom")
ggsave(plot = plot_sws_depthmax_2010, filename = "tests/sws_depthmax_rast_2010.png", height = 5, width = 8)

# Same test but for csv format download
# NB: Slow, but does run. Takes a few minutes.
current_bdmax_csv <- download_layers(dataset_id = "sws_baseline_2000_2019_depthmax", 
                                     variables = c("sws_ltmax", "sws_ltmin"), 
                                     constraints = set_con, fmt = "csv", directory = dir)
current_bdmax_csv_test <- current_bdmax_csv |> na.omit() |> 
  dplyr::rename(lon = longitude, lat = latitude) |> 
  mutate(lon = round(lon, 4), lat = round(lat, 4)) |>
  mutate(max_min_csv = ifelse(sws_ltmax > sws_ltmin, TRUE, FALSE)) |> 
  mutate(variable = case_when(time == "2000-01-01T00:00:00Z" ~ "max_min_csv_2000",
                              time == "2010-01-01T00:00:00Z" ~ "max_min_csv_2010", 
                              TRUE ~ as.character(NA))) |> 
  dplyr::select(-time, -sws_ltmax, -sws_ltmin) |> 
  pivot_wider(names_from = variable, values_from = max_min_csv)
plot_sws_depthmax_csv_2000 <- ggplot(data = current_bdmax_csv_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min_csv_2000)) + coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Bottom currents v3.0; csv_2000: max vs min") +
  theme(legend.position = "bottom")
ggsave(plot = plot_sws_depthmax_csv_2000, filename = "tests/sws_depthmax_csv_2000.png", height = 5, width = 8)
plot_sws_depthmax_csv_2010 <- ggplot(data = current_bdmax_csv_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min_csv_2010)) + coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Bottom currents v3.0; csv_2010: max vs min") +
  theme(legend.position = "bottom")
ggsave(plot = plot_sws_depthmax_csv_2010, filename = "tests/sws_depthmax_csv_2010.png", height = 5, width = 8)

# Same test but for NetCDF format download
current_bdmax_nc <- download_layers(dataset_id = "sws_baseline_2000_2019_depthmax", 
                                     variables = c("sws_ltmax", "sws_ltmin"), 
                                     constraints = set_con, fmt = "nc", directory = dir)
current_bdmax_nc_test <- current_bdmax_nc$data |> na.omit() |> 
  dplyr::rename(lon = longitude, lat = latitude) |> 
  mutate(lon = round(lon, 4), lat = round(lat, 4)) |>
  mutate(max_min_nc = ifelse(sws_ltmax > sws_ltmin, TRUE, FALSE)) |> 
  mutate(variable = case_when(time == "2000-01-01T00:00:00Z" ~ "max_min_nc_2000",
                              time == "2010-01-01T00:00:00Z" ~ "max_min_nc_2010", 
                              TRUE ~ as.character(NA))) |> 
  dplyr::select(-time, -sws_ltmax, -sws_ltmin) |> 
  pivot_wider(names_from = variable, values_from = max_min_nc)
plot_sws_depthmax_nc_2000 <- ggplot(data = current_bdmax_nc_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min_nc_2000)) + coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Bottom currents v3.0; nc_2000: max vs min") +
  theme(legend.position = "bottom")
ggsave(plot = plot_sws_depthmax_nc_2000, filename = "tests/sws_depthmax_nc_2000.png", height = 5, width = 8)
plot_sws_depthmax_nc_2010 <- ggplot(data = current_bdmax_nc_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min_nc_2010)) + coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Bottom currents v3.0; nc_2010: max vs min") +
  theme(legend.position = "bottom")
ggsave(plot = plot_sws_depthmax_nc_2010, filename = "tests/sws_depthmax_nc_2010.png", height = 5, width = 8)

# Join all tests for further plotting
current_bdmax_all_test <- left_join(current_bdmax_test, current_bdmax_csv_test) |> 
  left_join(current_bdmax_nc_test)

# Compare rast_2000 to rast_2010
plot_sws_depthmax_rast_2000_v_rast_2010 <- current_bdmax_all_test |> 
  mutate(max_min_rast_2000_v_max_min_rast_2010 = if_else(max_min_rast_2000 == max_min_rast_2010, TRUE, FALSE)) |> 
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min_rast_2000_v_max_min_rast_2010)) + coord_quickmap(expand = F) +
  labs(fill = "Mismatches are the same", x = NULL, y = NULL,
       title = "Bottom currents v3.0; rast_2000 vs rast_2010") +
  theme(legend.position = "bottom")
ggsave(plot = plot_sws_depthmax_rast_2000_v_rast_2010, filename = "tests/sws_depthmax_rast_2000_v_rast_2010.png", height = 5, width = 8)

# Compare rast_2000 to csv_2000
plot_sws_depthmax_rast_2000_v_csv_2000 <- current_bdmax_all_test |> 
  mutate(max_min_rast_2000_v_max_min_csv_2000 = if_else(max_min_rast_2000 == max_min_csv_2000, TRUE, FALSE)) |> 
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min_rast_2000_v_max_min_csv_2000)) + coord_quickmap(expand = F) +
  labs(fill = "Mismatches are the same", x = NULL, y = NULL,
       title = "Bottom currents v3.0; rast_2000 vs csv_2000") +
  theme(legend.position = "bottom")
ggsave(plot = plot_sws_depthmax_rast_2000_v_csv_2000, filename = "tests/sws_depthmax_rast_2000_v_csv_2000.png", height = 5, width = 8)

# Compare rast_2000 to csv_2000
plot_sws_depthmax_rast_2010_v_nc_2010 <- current_bdmax_all_test |> 
  mutate(max_min_rast_2010_v_max_min_nc_2010 = if_else(max_min_rast_2010 == max_min_nc_2010, TRUE, FALSE)) |> 
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min_rast_2010_v_max_min_nc_2010)) + coord_quickmap(expand = F) +
  labs(fill = "Mismatches are the same", x = NULL, y = NULL,
       title = "Bottom currents v3.0; rast_2010 vs nc_2010") +
  theme(legend.position = "bottom")
ggsave(plot = plot_sws_depthmax_rast_2010_v_nc_2010, filename = "tests/sws_depthmax_rast_2010_v_nc_2010.png", height = 5, width = 8)


# Test surface current layer ----------------------------------------------

# Info
info_layer("sws_baseline_2000_2019_depthsurf")

# Test surface currents
current_surface_layers <- download_layers(dataset_id = "sws_baseline_2000_2019_depthsurf", 
                                          variables = c("sws_ltmax", "sws_ltmin"), 
                                          constraints = set_con)
curvel_ss_test <- as.data.frame(current_surface_layers, xy = T) |> 
  dplyr::rename(lon = x, lat = y) |> na.omit() |> 
  mutate(lon = round(lon, 4), lat = round(lat, 4)) |> 
  mutate(max_min_rast_2000 = ifelse(sws_ltmax_1 > sws_ltmin_1, TRUE, FALSE),
         max_min_rast_2010 = ifelse(sws_ltmax_2 > sws_ltmin_2, TRUE, FALSE)) |> 
  dplyr::select(-sws_ltmax_1, -sws_ltmax_2, -sws_ltmin_1, -sws_ltmin_2)

# Visualise pixels where the max and min values are not as expected
plot_sws_depthsurf_1 <- ggplot(data = curvel_ss_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min_rast_2000)) + coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Surface currents v3.0; rast_2000: max vs min") +
  theme(legend.position = "bottom")
ggsave(plot = plot_sws_depthsurf_1, filename = "tests/sws_depthsurf_rast_2000.png", height = 5, width = 8)

# Visualise pixels where the max and min values are not as expected
plot_sws_depthsurf_2 <- ggplot(data = curvel_ss_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min_rast_2010)) + coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Surface currents v3.0; rast_2010: max vs min") +
  theme(legend.position = "bottom")
ggsave(plot = plot_sws_depthsurf_2, filename = "tests/sws_depthsurf_ltmax_2.png", height = 5, width = 8)

# Same test but for csv format download
# NB: Slow, but does run. Takes a few minutes.
current_surf_csv <- download_layers(dataset_id = "sws_baseline_2000_2019_depthsurf", 
                                     variables = c("sws_ltmax", "sws_ltmin"), 
                                     constraints = set_con, fmt = "csv", directory = dir)
current_surf_csv_test <- current_surf_csv |> na.omit() |> 
  dplyr::rename(lon = longitude, lat = latitude) |> 
  mutate(lon = round(lon, 4), lat = round(lat, 4)) |>
  mutate(max_min = ifelse(sws_ltmax > sws_ltmin, TRUE, FALSE))
plot_sws_surf_ltmax_csv <- ggplot(data = current_surf_csv_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min)) + coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "BottomSurface currents v3.0; sws_ltmax_csv vs sws_ltmin_csv") +
  theme(legend.position = "bottom")
ggsave(plot = plot_sws_surf_ltmax_csv, filename = "tests/sws_surf_ltmax_csv.png", height = 5, width = 8)

# Same test but for NetCDF format download
current_bdmax_nc <- download_layers(dataset_id = "sws_baseline_2000_2019_depthmax", 
                                    variables = c("sws_ltmax", "sws_ltmin"), 
                                    constraints = set_con, fmt = "nc", directory = dir)
current_bdmax_nc_test <- current_bdmax_nc$data |> na.omit() |> 
  dplyr::rename(lon = longitude, lat = latitude) |> 
  mutate(lon = round(lon, 4), lat = round(lat, 4)) |>
  mutate(max_min = ifelse(sws_ltmax > sws_ltmin, TRUE, FALSE))
plot_sws_depthmax_ltmax_nc <- ggplot(data = current_bdmax_nc_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min)) + coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Bottom currents v3.0; sws_ltmax_nc vs sws_ltmin_nc") +
  theme(legend.position = "bottom")
ggsave(plot = plot_sws_depthmax_ltmax_nc, filename = "tests/sws_depthmax_ltmax_nc.png", height = 5, width = 8)


# Test SST layer ----------------------------------------------------------

# Info
info_layer("thetao_baseline_2000_2019_depthsurf")

# Test surface currents
thetao_surface_layers <- download_layers(dataset_id = "thetao_baseline_2000_2019_depthsurf", 
                                          variables = c("thetao_ltmax", "thetao_ltmin"), 
                                          constraints = set_con)
thetao_ss_test <- as.data.frame(thetao_surface_layers, xy = T) |> 
  dplyr::rename(lon = x, lat = y) |> na.omit() |> 
  mutate(lon = round(lon, 4), lat = round(lat, 4)) |> 
  mutate(max_min_1 = ifelse(thetao_ltmax_1 > thetao_ltmin_1, TRUE, FALSE),
         max_min_2 = ifelse(thetao_ltmax_2 > thetao_ltmin_2, TRUE, FALSE))

# Visualise pixels where the max and min values are not as expected
plot_thetao_depthsurf_1 <- ggplot(data = thetao_ss_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min_1)) + coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Surface temperature v3.0; thetao_ltmax_1 vs thetao_ltmin_2") +
  theme(legend.position = "bottom")
ggsave(plot = plot_thetao_depthsurf_1, filename = "tests/thetao_depthsurf_ltmax_1.png", height = 5, width = 8)

# Visualise pixels where the max and min values are not as expected
plot_thetao_depthsurf_2 <- ggplot(data = curvel_ss_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min_2)) + coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Surface temperature v3.0; thetao_ltmax_2 vs thetao_ltmin_2") +
  theme(legend.position = "bottom")
ggsave(plot = plot_thetao_depthsurf_2, filename = "tests/thetao_depthsurf_ltmax_2.png", height = 5, width = 8)

