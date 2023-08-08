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
# 10: Convert .RData to .grd


# 1: Setup ----------------------------------------------------------------

# Load study sites and base packages
source("code/00_functions.R")

# Load additional packages
library(biomod2)
library(raster)
library(FNN)
library(doParallel)
# library(usdm)
# library(corrplot)

# The species occurrence data
sps_files <- dir("data/spp_presence", full.names = T, pattern = "final")
sps_names <- str_remove(dir("data/spp_presence", full.names = F, pattern = "final"), pattern = "_final.csv")


# 2: Load data ------------------------------------------------------------

# Load RData files for present, 2050, 21000
load("data/BO_present.RData")
load("data/BO_2050.RData")
load("data/BO_2100.RData")

# Coordinates only
global_coords <- dplyr::select(BO_present, lon, lat) |> 
  mutate(env_index = 1:nrow(BO_present))

# Convert to raster stacks
stack_present <- stack(rasterFromXYZ(BO_present, crs = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
stack_2050 <- stack(rasterFromXYZ(BO_2050, crs = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
stack_2100 <- stack(rasterFromXYZ(BO_2100, crs = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
plot(stack_present) # Visualise raster stack

# Very small area for testing
stack_present_baby <- stack(rasterFromXYZ(filter(BO_present, lon >= -2, lon <= 2, lat >= 50, lat <= 52), 
                                          crs = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
plot(stack_present_baby)

# Choose a species for testing the code
# sps_choice <- sps_files[1]

# The full pipeline wrapped into a function
biomod_pipeline <- function(sps_choice){
  
  # The species abbreviation
  sps_name <- str_remove(sapply(strsplit(sps_choice, "/"), "[[", 3), "_final.csv")
  
  # Message for current run
  print(paste0("Began run on ",sps_name))
  
  # Load the species
  sps <- read_csv(sps_choice) %>% 
    mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("lon", "lat")]),
                                            as.matrix(.[,3:4]), k = 1))) %>%
    left_join(global_coords, by = "env_index") %>% 
    dplyr::select(Sps, lon, lat) |> 
    distinct()
  
  # Set temp folder save locations
  # See: http://www.r-forge.r-project.org/forum/forum.php?thread_id=30946&forum_id=995&group_id=302
  # sps_path <- paste0("data/spp_projection/",sps_name)
  dir.create(file.path(sps_name), showWarnings = FALSE)
  dir.create(file.path(paste0(sps_name,"/Temp")), showWarnings = FALSE)
  rasterOptions(tmpdir = paste0(sps_name,"/Temp"))
  
  
  # 3: Prep data ------------------------------------------------------------
  
  # Prep data for modeling
  # system.time(
  biomod_data <- BIOMOD_FormatingData(
    resp.var = rep(1, nrow(sps)),
    resp.xy = as.matrix(sps[,2:3]),
    resp.name = sps_name,
    expl.var = stack_present, #stack_present_baby
    PA.strategy = 'random', 
    PA.nb.rep = 3, #3, 
    PA.nb.absences = 5000)
  # ) # 18 seconds
  
  # biomod_data # object summary
  # plot(biomod_data) # plot selected pseudo-absences
  
  # Save the pre-model data for possible later use
  saveRDS(biomod_data, file = paste0(sps_name,"/",sps_name,".base.Rds"))
  # biomod_data <- readRDS(paste0(sps_name,"/",sps_name,".base.Rds"))
  
  
  # 4: Model ----------------------------------------------------------------
  
  # Model options
  biomod_option <- BIOMOD_ModelingOptions()
  
  # Load if the model has already been run
  # if(file.exists())
  # biomod_model <- loadRData(paste0(sps_path,"/",sps_name,".",sps_name,".models.out"))
  
  # Run the model
  # system.time(
  biomod_model <- BIOMOD_Modeling(
    bm.format = biomod_data,
    modeling.id = sps_name,
    models = c("RF", "GLM", "GAM", "ANN", "MARS"),
    bm.options = biomod_option,
    CV.strategy = "random",
    CV.nb.rep = 5, #5,
    CV.perc = 0.7,
    metric.eval = c("TSS", "ROC", "KAPPA", "ACCURACY"),
    var.import = 3, # Number of permutations to estimate variable importance
    scale.models = TRUE,
    CV.do.full.models = FALSE,
    nb.cpu = 15,
    do.progress = TRUE)
  # ) # 379 seconds for full data 75 model run
  
  # Load if the model has already been run
  # if(file.exists())
  # biomod_ensemble <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,"ensemble.models.out"))
  
  # Build the ensemble models
  # system.time(
  biomod_ensemble <- BIOMOD_EnsembleModeling(
    bm.mod = biomod_model,
    models.chosen = "all",  # defines models kept (useful for removing non-preferred models)
    em.by = "all",
    em.algo = c("EMmean", "EMcv", "EMci"), 
    metric.select = c("TSS"),
    metric.select.thresh = c(0.7), # Turn this off during testing if the ensemble won't run...
    metric.eval = c("TSS", "ROC"), #, "KAPPA, "ACCURACY"),
    EMci.alpha = 0.05,
    var.import = 3, #0 # This takes a massive amount of time, the default is 0
    nb.cpu = 15)
  # ) # ~40 minutes for full data 75 model run
  
  
  # 5. Present projections --------------------------------------------------
  
  # Create projections
  # system.time(
  biomod_projection <- BIOMOD_Projection(
    bm.mod = biomod_model,
    proj.name = "present",
    new.env = stack_present, #stack_present_baby,
    models.chosen = "all",
    metric.binary = "TSS",
    compress = "gzip",
    build.clamping.mask = FALSE,
    nb.cpu = 15,
    output.format = ".tif")
  # ) # ~22 minutes for full data 75 model run
  # plot(biomod_projection)
  
  # Create ensemble projections
  system.time(
  biomod_ensemble_projection <- BIOMOD_EnsembleForecasting(
    bm.em = biomod_ensemble,
    bm.proj = biomod_projection,
    metric.binary = "TSS",
    compress = "gzip",
    nb.cpu = 15,
    output.format = ".tif")
  ) # xxx seconds for full data 75 model run
  
  # Visualise
  # plot(biomod_ensemble_projection)
  
  # Clean out some space
  rm(biomod_projection, biomod_ensemble_projection); gc()
  
  # Flush local tmp drive. Better not to do this if running on multiple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  # dir(tempdir())
  
  
  # 6: Future projections ---------------------------------------------------
  
  # Run 2050 projections
  # system.time(
  system.time(
    biomod_projection_2050 <- BIOMOD_Projection(
      bm.mod = biomod_model,
      proj.name = "2050",
      new.env = stack_2050, #stack_present_baby,
      models.chosen = "all",
      metric.binary = "TSS",
      compress = "gzip",
      build.clamping.mask = FALSE,
      nb.cpu = 15,
      output.format = ".tif")
  ) # xxx seconds for full data 75 model run
  # plot(biomod_projection_2050)
  
  # Create 2050 ensemble projections
  system.time(
    biomod_ensemble_projection_2050 <- BIOMOD_EnsembleForecasting(
      bm.em = biomod_ensemble,
      bm.proj = biomod_projection_2050,
      metric.binary = "TSS",
      compress = "gzip",
      nb.cpu = 15,
      output.format = ".tif")
  ) # xxx seconds for full data 75 model run
  # plot(biomod_ensemble_projection_2050)
  
  # Clean out 2050
  rm(biomod_projection_2050, biomod_ensemble_projection_2050); gc()
  
  # Flush local tmp drive. Better not to do this if running on multiple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  
  # Run 2100 projections
  # system.time(
  system.time(
    biomod_projection_2100 <- BIOMOD_Projection(
      bm.mod = biomod_model,
      proj.name = "2100",
      new.env = stack_2100, #stack_present_baby,
      models.chosen = "all",
      metric.binary = "TSS",
      compress = "gzip",
      build.clamping.mask = FALSE,
      nb.cpu = 15,
      output.format = ".tif")
  ) # xx seconds for full data 75 model run
  # plot(biomod_projection_2100)
  
  # Create 2050 ensemble projections
  system.time(
    biomod_ensemble_projection_2100 <- BIOMOD_EnsembleForecasting(
      bm.em = biomod_ensemble,
      bm.proj = biomod_projection_2100,
      metric.binary = "TSS",
      compress = "gzip",
      nb.cpu = 15,
      output.format = ".tif")
  ) # xxx seconds for full data 75 model run
  # plot(biomod_ensemble_projection_2100)
  
  # Clean out 2100
  rm(biomod_projection_2100, biomod_ensemble_projection_2100); gc()
  
  # Flush local tmp drive. Better not to do this if running on mulitple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  
  # Unlink Temp folder
  unlink(paste0(sps_name,"/Temp"), recursive = TRUE)
}


# 7: Run the pipeline -----------------------------------------------------

# Detect available cores automagically and set accordingly
# registerDoParallel(cores = detectCores()-1)

# Run one
# registerDoParallel(cores = 1)
# biomod_pipeline(sps_files[1])

# Run them all
registerDoParallel(cores = 15)
plyr::l_ply(sps_files, biomod_pipeline, .parallel = TRUE)


# 8: Analyse model output -------------------------------------------------

# Choose a species for the following code
sps_choice <- sps_names[1]

# Load chosen biomod_model and print evaluation scores
biomod_model <- loadRData(paste0(sps_choice,"/",sps_choice,".",sps_choice,".models.out"))
(Model_scores <- get_evaluations(biomod_model))
apply(Model_scores, c(1,2,3), mean, na.rm = T)
# dim(Model_scores)
# dimnames(Model_scores)

# Model evaluation by algorithm
models_scores_graph(biomod_model, by = "models", metrics = c('ROC','TSS'),
                    xlim = c(0.5,1), ylim = c(0.5,1)) + ggtitle("Algorithm") +
  geom_hline(aes(yintercept = 0.7), colour = "red", size = 2)

# Model evaluation by cross-validation
models_scores_graph(biomod_model, by = "cv_run", metrics = c('ROC','TSS'),
                    xlim = c(0.5,1), ylim = c(0.5,1)) + ggtitle("Run") +
  geom_hline(aes(yintercept = 0.7), colour = "red", size = 2)

# Model evaluation by dataset
models_scores_graph(biomod_model, by = "data_set", metrics = c('ROC','TSS'),
                    xlim = c(0.5,1), ylim = c(0.5,1)) + ggtitle("PA") +
  geom_hline(aes(yintercept = 0.7), colour = "red", size = 2)

## Calculate mean of variable importance by algorithm
# JG: I have read that scores reported are raw in the table (to be easier to interpret, 
# it should be normalized on our own - sum to 1 across algorithms)
(models_var_import <- get_variables_importance(biomod_model))
apply(models_var_import, c(1,2), mean, na.rm = T)
apply(apply(models_var_import, c(1,2), mean, na.rm = T), 1, mean) # Overall mean per variable

# Visualize species' modeled response to the given variables
# These necessary files are not on GitHub as they are too large
sp_name_Maxent <- BIOMOD_LoadModels(biomod_model, models = 'MAXENT.Phillips')
sp_name_GLM <- BIOMOD_LoadModels(biomod_model, models = 'GLM')
sp_name_ANN <- BIOMOD_LoadModels(biomod_model, models = 'ANN')
sp_name_RF <- BIOMOD_LoadModels(biomod_model, models = 'RF')
sp_name_GAM <- BIOMOD_LoadModels(biomod_model, models = 'GAM')
sp_name_ALL <- BIOMOD_LoadModels(biomod_model, models = c('MAXENT.Phillips', 'GLM', 'ANN', 'RF', 'GAM'))

# Evaluate individual models
Maxent_eval_strip <- biomod2::response.plot2(
  models = sp_name_Maxent,
  Data = get_formal_data(biomod_model, 'expl.var'),
  show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  do.bivariate = F,
  fixed.var.metric = 'mean',
  legend = F,
  display_title = F,
  data_species = get_formal_data(biomod_model, 'resp.var')
)
ggsave("", Maxent_eval_strip)
GLM_eval_strip <- biomod2::response.plot2(
  models = sp_name_GLM,
  Data = get_formal_data(biomod_model, 'expl.var'),
  show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  do.bivariate = F,
  fixed.var.metric = 'mean',
  legend = F,
  display_title = F,
  data_species = get_formal_data(biomod_model, 'resp.var')
)
ANN_eval_strip <- biomod2::response.plot2(
  models = sp_name_ANN,
  Data = get_formal_data(biomod_model, 'expl.var'),
  show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  do.bivariate = F,
  fixed.var.metric = 'mean',
  legend = F,
  display_title = F,
  data_species = get_formal_data(biomod_model, 'resp.var')
)
RF_eval_strip <- biomod2::response.plot2(
  models = sp_name_RF,
  Data = get_formal_data(biomod_model, 'expl.var'),
  show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  do.bivariate = F,
  fixed.var.metric = 'mean',
  legend = F,
  display_title = F,
  data_species = get_formal_data(biomod_model, 'resp.var')
)
GAM_eval_strip <- biomod2::response.plot2(
  models = sp_name_GAM,
  Data = get_formal_data(biomod_model, 'expl.var'),
  show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  do.bivariate = F,
  fixed.var.metric = 'mean',
  legend = F,
  display_title = F,
  data_species = get_formal_data(biomod_model, 'resp.var')
)

# Load and print ensemble model results
biomod_ensemble <- loadRData(paste0(sps_choice,"/",sps_choice,".",sps_choice,"ensemble.models.out"))
(models_scores_biomod_ensemble <- get_evaluations(biomod_ensemble))
(models_var_import <- get_variables_importance(biomod_ensemble))
apply(models_var_import, c(1,2), mean, na.rm = T)
apply(apply(models_var_import, c(1,2), mean, na.rm = T), 1, mean) # Overall mean per variable


# 9: Visualise ensemble models --------------------------------------------

# Load data used for maps etc.
load("data/Arctic_AM.RData")
colnames(Arctic_AM)[4] <- "depth"
Arctic_AM <- Arctic_AM %>%
  mutate(lon = round(lon, 4), lat = round(lat, 4))

# Choose a species
# sps_choice <- sps_names[1]

# Function that outputs BIOMOD projection comparison figures
plot_biomod <- function(sps_choice){
  # Load the species points
  sps_points <- read_csv(sps_files[str_which(sps_files,sps_choice)]) %>% 
    mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("lon", "lat")]),
                                            as.matrix(.[,2:3]), k = 1))) %>%
    left_join(global_coords, by = "env_index") %>% 
    dplyr::select(Sp, lon.y, lat.y) %>%
    dplyr::rename(lon = lon.y, lat = lat.y)
  
  # Load the ensemble projections
  biomod_project_present <- loadRData(paste0(sps_choice,"/proj_present/proj_present_",sps_choice,"_ensemble_TSSbin.RData"))
  biomod_project_2050 <- loadRData(paste0(sps_choice,"/proj_2050/proj_2050_",sps_choice,"_ensemble_TSSbin.RData"))
  biomod_project_2100 <- loadRData(paste0(sps_choice,"/proj_2100/proj_2100_",sps_choice,"_ensemble_TSSbin.RData"))
  
  # Convert to data.frames
  rast_df <- function(rast){
    df_out <- as.data.frame(rast[[1]], xy = T) %>% 
      `colnames<-`(c("lon", "lat", "presence")) %>% 
      mutate(lon = round(lon, 4), lat = round(lat, 4)) %>% 
      left_join(Arctic_AM, by = c("lon", "lat")) %>% 
      na.omit() 
  }
  df_project_present <- rast_df(biomod_project_present[[1]])
  df_project_2050 <- rast_df(biomod_project_2050[[1]])
  df_project_2100 <- rast_df(biomod_project_2100[[1]])
  
  # Visualise present data
  plot_present <- df_project_present %>% 
    filter(land_distance <= 50 | depth <= 100) %>% 
    filter(presence == TRUE) %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = presence)) +
    borders(fill = "grey90", colour = "black") +
    geom_point(data = sps_points, colour = "yellow", size = 0.5) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    scale_fill_manual(values = c("forestgreen")) +
    labs(x = NULL, y = NULL, title = paste0(sps_choice,": Present")) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # Function for visualising changes over time
  plot_diff <- function(df_future, year_label){
    plot_out <- left_join(df_project_present, df_future, 
                          by = c("lon", "lat", "land_distance", "depth")) %>% 
      mutate(change = presence.x - presence.y,
             change = ifelse(presence.x == F & presence.y == F, NA, change),
             change =  factor(change, 
                              levels = c("-1", "0", "1"),
                              labels = c("gain", "no change", "loss"))) %>% 
      na.omit() %>% 
      filter(land_distance <= 100 | depth <= 100) %>% 
      ggplot(aes(x = lon, y = lat)) +
      geom_tile(aes(fill = change)) +
      borders(fill = "grey90", colour = "black") +
      scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
      scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
      coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                     ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
      scale_fill_brewer(palette = "Set1", direction = -1) +
      labs(x = NULL, y = NULL, title = paste0(sps_choice,": Present - ",year_label)) +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
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
      ncol = 3,
      align = "v"),
    cowplot::plot_grid(
      cowplot::get_legend(plot_present),
      ggplot() + theme_void(),
      cowplot::get_legend(plot_2100),
      ncol = 3, rel_widths = c(1, 0, 1.5)),
    nrow = 2, rel_heights = c(10,1)
  )
  ggsave(paste0("graph/biomod_diff_",sps_choice,".png"), plot_ALL, width = 8, height = 4.7)
}

# Create all visuals
registerDoParallel(cores = 5)
plyr::l_ply(sps_names, plot_biomod, .parallel = T)


# 10: Save ensemble models as .grd files ----------------------------------

# Write function to load results saved as .RData files and save them as .grd files

# RWS: NB: I can't guarantee this code won't break when the required packages are updated by the authors.

# Save as a raster file
RData_to_grd <- function(df){
  sps_choice <- df$sps_choice; proj_choice <- df$proj_choice
  proj_data <- loadRData(paste0(sps_choice,"/proj_",proj_choice,"/proj_",proj_choice,
                                "_",sps_choice,"_ensemble_TSSbin.RData"))
  proj_names <- sapply(strsplit(names(proj_data), "_"), "[[", 2)
  writeRaster(x = proj_data, bylayer = TRUE, suffix = proj_names, overwrite = TRUE,
              filename = paste0("data/ascii_results/proj_",proj_choice,"_",sps_choice,"_ensemble_TSSbin.asc"))
}

# Run them all
quick_grid <- expand.grid(sps_choice = sps_names, 
                          proj_choice = c("present", "2050", "2100"), stringsAsFactors = F) %>% 
  mutate(plyr_id = 1:n()) %>% 
  data.frame()
plyr::d_ply(quick_grid, c("plyr_id"), RData_to_grd, .parallel = T)

# Check the output
test_rast <- raster("data/ascii_results/proj_2050_Acla_ensemble_TSSbin_EMmeanByTSS.asc")
plot(test_rast)

