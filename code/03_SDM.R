# code/03_SDM.R
# Herein the code for running the SDM
# The order of operations is:
# 1: Setup the environment
# 2: Load data
# 3: Prep data for modelling
# 4: Run model ensemble
# 5: Present projections
# 6: Future projections
# 7: Run the pipeline
# 8: Analyse model output
# 9: Model visualisations
# 10: Copy files to a folder


# 1: Setup ----------------------------------------------------------------

# Load study sites and base packages
source("code/00_functions.R")

# Load additional packages
library(biomod2)
library(FNN)
# library(doParallel) # NB: DO NOT use this library. It has strange interacitons with biomod2 multi-core behaviour.

# The species occurrence data
sps_list <- read_csv("~/PanArcticInvasion/metadata/Species distribution key_functional groups_Mar6_2025.csv") |> arrange(Acronym)
sps_files <- dir("~/PanArcticInvasion/data/spp_presence", full.names = T, pattern = "final")

# The acronyms from the presence files
sps_accs <- str_remove(dir("~/PanArcticInvasion/data/spp_presence", 
                           full.names = F, pattern = "final"), pattern = "_final.csv")

# Ensure files are present for all species in the list
# sps_check <- sps_list |> dplyr::filter(!Acronym %in% sps_accs)
# sps_files_check <- data.frame(sps_acronym = sps_accs) |> dplyr::filter(!sps_acronym %in% sps_list$Acronym)
# length(unique(sps_list$Acronym)); length(sps_files_check$sps_acronym)

# Test .csv files for 'Sps' 'Long' and 'Lat' columns, print file names if not
# sps_test <- lapply(sps_files, function(x) read_csv(x) |> dplyr::select(Sps, Long, Lat))
# sps_test <- read_csv(sps_files[175])# |> dplyr::select(Sps, Long, Lat)


# 2: Load data ------------------------------------------------------------

# Load TIFF files for present, 2050, 2100
if(!exists("stack_present")) stack_present <- terra::rast("data/BO_v3/BO_present.tiff")
if(!exists("stack_2050")) stack_2050 <- terra::rast("data/BO_v3/BO_2050.tiff")
if(!exists("stack_2100")) stack_2100 <- terra::rast("data/BO_v3/BO_2100.tiff")
# plot(stack_present) # Visualise raster stack

# NB: These stack files contain the full global data
# If it is decided to stay with the lat 30 to 90 range they must be cropped

# Coordinates only
global_coords <- as.data.frame(stack_present[[1]], xy = T) |> 
  dplyr::select(x,y) |> dplyr::rename(lon = x, lat = y) |> 
  distinct() |> mutate(env_index = 1:n())

# Choose a species for testing the code
# sps_choice <- sps_files[1]

# The full pipeline wrapped into a function
biomod_pipeline <- function(sps_choice, force_run = FALSE, test_run = FALSE, n_cores = 4){
  
  # The species abbreviation
  # NB: This will change based on the users folder structure
  sps_acc <- str_remove(sapply(strsplit(sps_choice, "/"), "[[", 7), "_final.csv")
  
  # Skip species if a folder has already been created for it and force_run is FALSE
  if(file.exists(paste0("~/pCloud Drive/PanArctic/data/spp_projection/",sps_acc)) & !force_run){
    return("Species has already been modelled and force_run is FALSE")
  }
  
  # Message for current run
  print(paste0("Began run on ",sps_acc))
  
  # Load the species
  print(paste0("Prepping species data at ", Sys.time()))
  sps <- read_csv(sps_choice) |>
    filter(!is.na(Long), !is.na(Lat)) |> 
    dplyr::select(Sps, Long, Lat) |> 
    dplyr::mutate(Long = as.numeric(Long), Lat = as.numeric(Lat))
  sps$env_index = as.vector(knnx.index(as.matrix(global_coords[,c("lon", "lat")]),
                                       as.matrix(sps[,2:3]), k = 1))
  sps <- left_join(sps, global_coords, by = "env_index") |> 
    dplyr::select(Sps, lon, lat) |> distinct()
  
  # Set temp folder save locations
  # See: http://www.r-forge.r-project.org/forum/forum.php?thread_id=30946&forum_id=995&group_id=302
  sps_path <- paste0("~/pCloud Drive/PanArctic/data/spp_projection/",sps_acc)
  dir.create(file.path(sps_path), showWarnings = FALSE)
  dir.create(file.path(paste0(sps_path,"/Temp")), showWarnings = FALSE)
  rasterOptions(tmpdir = paste0(sps_path,"/Temp"))
  
  
  # 3: Prep data ------------------------------------------------------------
  
  # Filter environmental data based on species classification
  # NB: This is based on the Feb 2025 list
  sps_list_sub <- filter(sps_list,  Acronym == sps_acc)
  if(sps_list_sub$`Functional Group`[1] == "Fish"){
    stack_layers <- c("thetao_mean", "so_mean", "sithick_mean", "chl_mean", "no3_mean")
  } else if(sps_list_sub$`Functional Group`[1] %in% c("Phytobenthos", "Phytoplankton")){
    stack_layers <- c("thetao_mean", "so_mean", "sithick_mean", "no3_mean", "dfe_mean", "par_mean_mean")
  } else if(sps_list_sub$`Functional Group`[1] %in% c("Zoobenthos", "Zooplankton")){
    stack_layers <- c("thetao_mean", "so_mean", "sithick_mean", "chl_mean")
  } else {
    return("Species not classified")
  }
  
  # Filter raster stacks by layers
  stack_present_sub <- stack_present[[stack_layers]]
  stack_2050_sub <- stack_2050[[stack_layers]]
  stack_2100_sub <- stack_2100[[stack_layers]]
  
  if(test_run){
    # Very small area for testing
    stack_present_sub <- crop(stack_present_sub, ext(-2, 2, 50, 52))
    stack_2050_sub <- crop(stack_2050_sub, ext(-2, 2, 50, 52))
    stack_2100_sub <- crop(stack_2100_sub, ext(-2, 2, 50, 52))
    # plot(stack_present_baby)
  }

  # Split up data.frame for testing
  ## NB: BIOMOD_FormatingData throws errors for any eval datasets that don't contain absences
  ## So this is not used here and rather the built-in cross-validation is performed in the following step
  # set.seed(13)
  # sps_resp <- sps[sample(nrow(sps)*0.7),]
  # sps_eval <- sps[sample(nrow(sps)*0.3),]
  
  # Prep data for modelling
  print(paste0("Prepping model data at ", Sys.time()))
  if(test_run){
    abs_count <- 100
  } else {
    abs_count <- 5000
  }
  # system.time(
  biomod_data <- BIOMOD_FormatingData(
    resp.name = sps_acc,
    resp.var = rep(1, nrow(sps)),
    expl.var = stack_present_sub,
    resp.xy = as.matrix(sps[,2:3]),
    # eval.resp.var = rep(1, nrow(sps_eval)),
    # eval.expl.var = stack_present_sub,
    # eval.resp.xy = as.matrix(sps_eval[,2:3]),
    filter.raster = TRUE, # Removes duplicate points
    PA.strategy = "random",
    PA.nb.rep = 1, # Intentionally only one set of 5,000 pseudo absences is generated
    PA.nb.absences = abs_count
    )
  # ) # 0 seconds on test, 188 seconds on full run
  
  # Visualise
  # biomod_data # object summary
  # plot(biomod_data) # plot selected pseudo-absences
  
  # Save the pre-model data for possible later use
  saveRDS(biomod_data, file = paste0(sps_path,"/",sps_acc,".base.Rds"))
  # biomod_data <- readRDS(paste0(sps_name,"/",sps_name,".base.Rds"))
  
  
  # 4: Model ----------------------------------------------------------------
  
  # Run the model
  print(paste0("Modelling species distribution at ", Sys.time()))
  if(test_run){
    nb_rep <- 1; var_imp <- 1
  } else {
    nb_rep <- 5; var_imp <- 3
  }
  # system.time(
  biomod_model <- BIOMOD_Modeling(
    bm.format = biomod_data,
    modeling.id = sps_acc,
    models = c("RF", "GLM", "GAM", "ANN", "MARS"),
    CV.strategy = "random",
    CV.nb.rep = nb_rep,
    CV.perc = 0.7,
    metric.eval = c("TSS", "ROC"),#, "KAPPA", "ACCURACY"),
    var.import = var_imp,
    scale.models = TRUE,
    CV.do.full.models = FALSE,
    nb.cpu = n_cores,
    do.progress = TRUE)
  # ) # 3 seconds for test run, 20 seconds for full run
  
  # Build the ensemble models
  print(paste0("Building ensemble model at ", Sys.time()))
  # system.time(
  biomod_ensemble <- BIOMOD_EnsembleModeling(
    bm.mod = biomod_model,
    models.chosen = "all",
    em.by = "all",
    em.algo = c("EMmean", "EMcv", "EMci"),
    # There are potential issues here after biomod2 v4.2+
    # https://github.com/biomodhub/biomod2/issues/166
    metric.select = "TSS",
    metric.select.thresh = 0.7, # -2, # Allows everything
    #
    metric.eval = c("TSS", "ROC"),#, "KAPPA", "ACCURACY"),
    var.import = 5, # 5 final # 0 testing
    EMci.alpha = 0.05,
    nb.cpu = n_cores)
  # ) # 17 seconds for full ensemble run

  # If all models fail the TSS > 0.7 cutoff exit the pipeline
  if(!file.exists(paste0(sps_path,"/",sps_acc,".",sps_acc,".ensemble.models.out"))){
    return("All models failed")
  } 
  
  
  # 5. Present projections --------------------------------------------------
  
  # Create projections
  print(paste0("Creating present-day projections at ", Sys.time()))
  # system.time(
  biomod_projection <- BIOMOD_Projection(
    bm.mod = biomod_model,
    proj.name = "present",
    new.env = stack_present_sub,
    models.chosen = "all",
    metric.binary = "TSS",
    compress = "gzip",
    build.clamping.mask = FALSE,
    nb.cpu = n_cores,
    output.format = ".tif")
  # ) # ~8 minutes for full model run
  # plot(biomod_projection)
  
  # Create ensemble projections
  print(paste0("Creating present-day ensemble at ", Sys.time()))
  # system.time(
  biomod_ensemble_projection <- BIOMOD_EnsembleForecasting(
    bm.em = biomod_ensemble,
    bm.proj = biomod_projection,
    metric.binary = "TSS",
    compress = "gzip",
    nb.cpu = n_cores,
    output.format = ".tif")
  # ) # ~4 minutes for full model run
  
  # Visualise
  # plot(biomod_ensemble_projection)
  
  # Clean out some space
  rm(biomod_projection, biomod_ensemble_projection); gc()
  
  # Flush local tmp drive. Better not to do this if running on multiple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  # dir(tempdir())
  
  
  # 6: Future projections ---------------------------------------------------
  
  # Run 2050 projections
  print(paste0("Creating 2050 projections at ", Sys.time()))
  # system.time(
  biomod_projection_2050 <- BIOMOD_Projection(
    bm.mod = biomod_model,
    proj.name = "2050",
    new.env = stack_2050_sub,
    models.chosen = "all",
    metric.binary = "TSS",
    compress = "gzip",
    build.clamping.mask = FALSE,
    nb.cpu = n_cores,
    output.format = ".tif")
  # ) # ~8 minutes for full model run
  # plot(biomod_projection_2050)
  
  # Create 2050 ensemble projections
  print(paste0("Creating 2050 ensemble at ", Sys.time()))
  # system.time(
  biomod_ensemble_projection_2050 <- BIOMOD_EnsembleForecasting(
    bm.em = biomod_ensemble,
    bm.proj = biomod_projection_2050,
    metric.binary = "TSS",
    compress = "gzip",
    nb.cpu = n_cores,
    output.format = ".tif")
  # ) # ~4 minutes for full model run
  # plot(biomod_ensemble_projection_2050)
  
  # Clean out 2050
  rm(biomod_projection_2050, biomod_ensemble_projection_2050); gc()
  
  # Flush local tmp drive. Better not to do this if running on multiple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  
  # Run 2100 projections
  print(paste0("Creating 2100 projections at ", Sys.time()))
  # system.time(
  biomod_projection_2100 <- BIOMOD_Projection(
    bm.mod = biomod_model,
    proj.name = "2100",
    new.env = stack_2100_sub,
    models.chosen = "all",
    metric.binary = "TSS",
    compress = "gzip",
    build.clamping.mask = FALSE,
    nb.cpu = n_cores,
    output.format = ".tif")
  # ) # ~8 minutes for full model run
  # plot(biomod_projection_2100)
  
  # Create 2100 ensemble projections
  print(paste0("Creating 2100 ensemble at ", Sys.time()))
  # system.time(
  biomod_ensemble_projection_2100 <- BIOMOD_EnsembleForecasting(
    bm.em = biomod_ensemble,
    bm.proj = biomod_projection_2100,
    metric.binary = "TSS",
    compress = "gzip",
    nb.cpu = n_cores,
    output.format = ".tif")
  # ) # ~4 minutes for full model run
  # plot(biomod_ensemble_projection_2100)
  
  # Clean out 2100
  rm(biomod_projection_2100, biomod_ensemble_projection_2100); gc()
  
  # Flush local tmp drive. Better not to do this if running on mulitple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  
  # Unlink Temp folder - this also shows that the pipeline finished successfully
  unlink(paste0(sps_path,"/Temp"), recursive = TRUE)
}


# 7: Run the pipeline -----------------------------------------------------


# Change the working directory so that biomod2 saves the results in a convenient folder
setwd("~/pCloud Drive/PanArctic/data/spp_projection/")

# Run one test
# system.time(biomod_pipeline(sps_files[73], test_run = TRUE))
# 19 seconds for 1 species on minimum reps and baby range
# system.time(biomod_pipeline(sps_files[1], test_run = FALSE))
# ~35 minutes for a full run

# Run them all
# NB: Do not run in parallel
plyr::l_ply(sps_files[1:105], biomod_pipeline, .parallel = FALSE,
            force_run = FALSE, test_run = FALSE, n_cores = 4)

## Error log

## Species that did not run
### By file count
# sps_file_count <- data.frame(files = list.files("~/PanArcticInvasion/data/spp_projection", full.names = T, recursive = T)) |> 
#   mutate(dir = sapply(strsplit(files, "/"), "[[", 7)) |> count(dir) |> filter(n == 24)
# sps_file_rerun <- sps_files[which(!sps_names %in% sps_file_count$dir)]
# plyr::l_ply(sps_file_rerun, biomod_pipeline, .parallel = FALSE)

### By date
# sps_date <- file.info(dir("~/PanArcticInvasion/data/spp_projection", full.names = T)) |> 
#   rownames_to_column(var = "folder_name") |> 
#   mutate(folder_name = sapply(strsplit(folder_name, "/"), "[[", 7)) |> 
#   filter(ctime < Sys.Date()) # Number of days in the past, may need to be adjusted
# sps_date_rerun <- sps_files[which(sps_folders$folder_name %in% sps_names)]
# plyr::l_ply(sps_date_rerun, biomod_pipeline, .parallel = TRUE)


# Set working directory back to project base
setwd("~/PanArcticInvasion/")


# 8: Analyse model output -------------------------------------------------

# Choose a species for the following code
# sps_choice <- sps_names[73]

# Load chosen biomod_model and print evaluation scores
# biomod_model <- loadRData(paste0("~/pCloud Drive/PanArctic/data/spp_projection/",sps_choice,"/",sps_choice,".",sps_choice,".models.out"))
# biomod_model
# get_evaluations(biomod_model)

# Model evaluation by 
# bm_PlotEvalMean(biomod_model, dataset = "calibration", metric.eval = c('ROC','TSS'))

# Model evaluation by 
# bm_PlotEvalMean(biomod_model, dataset = "validation", metric.eval = c('ROC','TSS'))

# Model evaluation by 
# bm_PlotEvalMean(biomod_model, dataset = "evaluation", metric.eval = c('ROC','TSS'))

# NB: Consider using evaluation code from '~/HBCproject/analyses/2.EvaluationScores.R'


# 9: Visualise ensemble models --------------------------------------------

# Choose a species
# sps_choice <- sps_names[73]

# Function that outputs BIOMOD projection comparison figures
# TODO: Run this step by step to ensure it stills behaves as expected
plot_biomod <- function(sps_choice){
 
  # File name
  sps_file <- sps_files[str_which(sps_files, sps_choice)]
  sps_dot_name <- str_replace(sps_choice, "_", ".")
  
   # Load the species points
  # Load the species
  sps_points <- read_csv(sps_file) |> filter(!is.na(Long), !is.na(Lat))
  sps_points$env_index = as.vector(knnx.index(as.matrix(global_coords[,c("lon", "lat")]),
                                       as.matrix(sps_points[,3:4]), k = 1))
  sps_points <- left_join(sps_points, global_coords, by = "env_index") |> 
    dplyr::select(Sps, lon, lat) |> distinct()
  
  # Load the ensemble projections
  biomod_project_present <- raster(paste0("~/pCloud Drive/PanArctic/data/spp_projection/",sps_dot_name,
                                          "/proj_present/proj_present_",sps_dot_name,"_ensemble_TSSbin.tif"))
  biomod_project_2050 <- raster(paste0("~/pCloud Drive/PanArctic/data/spp_projection/",sps_dot_name,
                                       "/proj_2050/proj_2050_",sps_dot_name,"_ensemble_TSSbin.tif"))
  biomod_project_2100 <- raster(paste0("~/pCloud Drive/PanArctic/data/spp_projection/",sps_dot_name,
                                       "/proj_2100/proj_2100_",sps_dot_name,"_ensemble_TSSbin.tif"))
  
  # Convert to data.frames
  df_project_present <- rast_df(biomod_project_present[[1]])
  df_project_2050 <- rast_df(biomod_project_2050[[1]])
  df_project_2100 <- rast_df(biomod_project_2100[[1]])
  
  # Visualise present data
  plot_present <- df_project_present %>% 
    filter(presence == 1) %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = as.factor(presence))) +
    borders(fill = "grey90", colour = "black") +
    geom_point(data = sps_points, colour = "yellow", size = 0.5) +
    scale_y_continuous(breaks = c(40, 60, 80), labels = c("40°N", "60°N", "80°N")) +
    scale_x_continuous(breaks = c(-160, -80, 0, 80, 160), labels = c("160°W", "80°W", "0°E", "60°E", "160°E")) +
    coord_quickmap(xlim = c(-180, 180), ylim = c(30, 90), expand = F) +
    scale_fill_manual(values = c("forestgreen")) +
    labs(x = NULL, y = NULL, fill = "presence", 
         title = paste0(sps_points$Sps)) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # Function for visualising changes over time
  plot_diff <- function(df_future, year_label){
    plot_out <- left_join(df_project_present, df_future, 
                          by = c("lon", "lat")) |> 
      mutate(change = presence.x - presence.y,
             change = ifelse(presence.x == F & presence.y == F, NA, change),
             change =  factor(change, 
                              levels = c("-1", "0", "1"),
                              labels = c("gain", "no change", "loss"))) |>  
      na.omit() |> 
      ggplot(aes(x = lon, y = lat)) +
      geom_tile(aes(fill = change)) +
      borders(fill = "grey90", colour = "black") +
      scale_y_continuous(breaks = c(40, 60, 80), labels = c("40°N", "60°N", "80°N")) +
      scale_x_continuous(breaks = c(-160, -80, 0, 80, 160), labels = c("160°W", "80°W", "0°E", "60°E", "160°E")) +
      coord_quickmap(xlim = c(-180, 180), ylim = c(30, 90), expand = F) +
      scale_fill_brewer(palette = "Set1", direction = -1) +
      labs(x = NULL, y = NULL, title = paste0("Present - ",year_label)) +
      theme_bw() +
      theme(legend.position = "bottom")#,
            # axis.text.x = element_blank(),
            # axis.ticks.x = element_blank())
  }
  
  # Visualise present - 2050
  plot_2050 <- plot_diff(df_project_2050, "2050")
  
  # Visualise present - 2100
  plot_2100 <- plot_diff(df_project_2100, "2100")
  
  # Combine and save
  plot_ALL <- cowplot::plot_grid(
    cowplot::plot_grid(
      plot_present + theme(legend.position = "none"),
      plot_2050 + theme(legend.position = "none"),
      plot_2100 + theme(legend.position = "none"),
      nrow = 3,
      align = "h"),
    cowplot::plot_grid(
      cowplot::get_legend(plot_2100)),
    nrow = 2, rel_heights = c(10,1)
  )
  cowplot::save_plot(paste0("figures/biomod_diff_",sps_choice,".png"), plot_ALL, base_width = 10, base_height = 12)
  # ggsave(paste0("figures/biomod_diff_",sps_choice,".png"), plot_ALL, width = 12, height = 12)
}

# Create all visuals
# registerDoParallel(cores = detectCores()-1)
# plyr::l_ply(sps_names, plot_biomod, .parallel = TRUE)

# Check that all figures run
# sps_fig_check <- str_remove(str_remove(dir("~/PanArcticInvasion/figures", 
#                                            full.names = FALSE), "biomod_diff_"), ".png")
# sps_fig_rerun <- sps_fig_check[which(!sps_names %in% sps_fig_check)]
# plyr::l_ply(sps_fig_rerun, plot_biomod, .parallel = FALSE)


# 10: Copy files to folder for sharing ------------------------------------

# Identify the folders
# spp_folder <- "data/spp_projection"
# new_folder <- "H:/Where I want my files to be copied to"

# Find the files that you want
# proj_files <- list.files(path = spp_folder, recursive = TRUE,
#                          pattern = "ensemble_TSSbin.tif", full.names = TRUE)

# Copy the files to the new folder
# file.copy(proj_files, "data/spp_projection_ensemble_TSSbin")

