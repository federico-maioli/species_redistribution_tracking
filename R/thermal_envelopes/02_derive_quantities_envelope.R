# load libraries ----------------------------------------------------------
library(tidyverse)
library(here)
library(sdmTMB)
library(furrr)

# define paths ------------------------------------------------------------
models_folder <- here("R/thermal_envelopes/fitted")
model_files <- list.files(models_folder, pattern = "\\.rds$", full.names = TRUE)

# define helper function --------------------------------------------------
process_model_file <- function(model_file) {
  
  # load spatial grid -----------------------------------------------------
  grid <- readRDS(here("R/data/processed/prediction_grid.rds")) %>%
    droplevels()
  
  # parse model metadata --------------------------------------------------
  file_name <- basename(model_file)
  parts <- strsplit(file_name, "_")[[1]]
  this_region <- parts[1]
  this_species <- parts[2]
  
  # load and prepare model ------------------------------------------------
  model <- readRDS(model_file)
  model <- sdmTMB:::reload_model(model)
  years <- unique(model$data$year)
  
  # prepare region-specific grid ------------------------------------------
  main_survey <- model$data %>%
    count(survey, sort = TRUE) %>%
    slice(1) %>%
    pull(survey) %>%
    as.factor()
  
  # this is a formula wgt_cpua ~ 0 + as.factor(year) + s(mean_temp_std, k = 3) + s(warmest_temp_std, k = 3)
  
  region_grid <- grid %>%
    filter(region_short == this_region, year %in% years) %>%
    mutate( # standardize variables
      mean_temp_std = (mean_temp - mean(model$data$mean_temp, na.rm = TRUE)) / 
        sd(model$data$mean_temp, na.rm = TRUE),
      warmest_temp_std = (warmest_temp - mean(model$data$warmest_temp, na.rm = TRUE)) /
        sd(model$data$warmest_temp, na.rm = TRUE),
      survey = main_survey,
      area = 16
    ) %>%
    na.omit() %>%
    droplevels()
  
  if (nrow(region_grid) == 0) {
    message("! Skipping model for region: ", this_region,
            " — no grid data for years: ", paste(years, collapse = ", "))
    return(NULL)
  }
  
  gc()
  
  # predict from model ----------------------------------------------------
  pred <- predict(model, newdata = region_grid, return_tmb_object = TRUE)
  gc()
  
  # derive quantities -----------------------------------------------------
  index <- get_index(pred, area = 16, bias_correct = TRUE) %>%
    transmute(year, index = log_est, index_se = se)
  
  cog <- get_cog(pred, area = 16, format = "wide") %>%
    transmute(
      year,
      cog_y = est_y,
      cog_x = est_x,
      cog_y_se = se_y,
      cog_x_se = se_x
    )
  
  depth_niche <- get_weighted_average(pred, vector = pred$data$depth, area = 16) %>%
    transmute(year, depth_niche = est, depth_niche_se = se)
  
  thermal_niche <- get_weighted_average(pred, vector = pred$data$mean_temp, area = 16) %>%
    transmute(year, thermal_niche = est, thermal_niche_se = se)
  
  # combine results -------------------------------------------------------
  out <- index %>%
    inner_join(cog, by = "year") %>% inner_join(depth_niche, by = "year") %>%
    inner_join(thermal_niche, by = "year")
  
  out <- out |> mutate(species = this_species, region = this_region)
  
  
  return(out)
}

# process models in parallel ----------------------------------------------
# Set up parallel processing with furrr
plan(multisession, workers = 2)

# Ensure reproducibility
options(future.seed = TRUE)

set.seed(123)

# Process models in parallel
data <- future_map(model_files, process_model_file, .progress = TRUE,.options = furrr_options(seed = TRUE))

# Combine all results into a single data frame
data <- bind_rows(data)

saveRDS(data,here('R/data/processed/derived_quantities_te.rds')) 
