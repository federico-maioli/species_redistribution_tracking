# load libraries ----------------------------------------------------------
library(tidyverse)
library(here)
library(sdmTMB)
library(furrr)


# define paths ------------------------------------------------------------
models_folder <- here("climate_velocity/fitted")
model_files <- list.files(models_folder, pattern = "\\.rds$", full.names = TRUE)

# define helper function --------------------------------------------------
process_model_file <- function(model_file) {
  
  # load spatial grid -----------------------------------------------------
  grid <- readRDS(here("data/processed/prediction_grid.rds")) %>%
    mutate(month_f = as.factor(min_month)) %>%
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
    mutate( # standatdize variables
      #logdepth = (log(depth) - mean(log(model$data$depth), na.rm = TRUE)) / 
      #  sd(log(model$data$depth), na.rm = TRUE), 
      #logdepth2 = logdepth^2,
      mean_temp_std = (mean_temp - mean(model$data$mean_temp, na.rm = TRUE)) / 
        sd(model$data$mean_temp, na.rm = TRUE),
      warmest_temp_std = (warmest_temp - mean(model$data$warmest_temp, na.rm = TRUE)) /
        sd(model$data$warmest_temp, na.rm = TRUE),
      survey = main_survey,
      area = 16
    ) %>%
    filter(!is.na(month_f)) %>%
    na.omit() %>%
    droplevels()
  
  if (nrow(region_grid) == 0) {
    message("! Skipping model for region: ", this_region,
            " — no grid data for years: ", paste(years, collapse = ", "))
    return(NULL)
  }
  
  # set model month -------------------------------------------------------
  min_model_month <- min(as.numeric(as.character(model$data$month_f)), na.rm = TRUE)
  region_grid$month_f <- factor(min_model_month)
  
  gc()
  
  # predict from model ----------------------------------------------------
  pred <- predict(model, newdata = region_grid, return_tmb_object = TRUE)
  gc()
  
  # derive quantities -----------------------------------------------------
  index <- get_index(pred, area = 16, bias_correct = TRUE) %>%
    transmute(year, index = log_est, index_se = se)
  
  eao <- get_eao(pred, area = 16) %>%
    transmute(year, eao = log_est, eao_se = se)
  
  cog <- get_cog(pred, area = 16, format = "wide") %>%
    transmute(
      year,
      cog_y = est_y,
      cog_x = est_x,
      cog_y_se = se_y,
      cog_x_se = se_x
    )
  
  # combine results -------------------------------------------------------
  out <- index %>%
    inner_join(cog, by = "year") %>%
    inner_join(eao, by = "year")
  
  niche <- pred$data %>% mutate(biomass = exp(est + log(pred$data$area))) |> # area expansion
    group_by(year) %>%
    summarise(
      depth_niche = sum(depth * biomass) / sum(biomass),
      thermal_niche = sum(mean_temp * biomass) / sum(biomass),
      #thermal_niche_min = sum(coldest_temp * biomass) / sum(biomass),
      #thermal_niche_max = sum(warmest_temp * biomass) / sum(biomass),
      #upr90_thermal = Hmisc::wtd.quantile(mean_temp, weights = biomass, probs = 0.9), # remember to reduce dependencies
      #lwr10_thermal = Hmisc::wtd.quantile(mean_temp, weights = biomass, probs = 0.1),
      #upr90_depth = Hmisc::wtd.quantile(depth, weights = biomass, probs = 0.9),
      #lwr10_depth = Hmisc::wtd.quantile(depth, weights = biomass, probs = 0.1),
    ) 
  
  out <- out |> left_join(niche) |> mutate(species = this_species, region = this_region)
  
  
  return(out)
}

# process models in parallel ----------------------------------------------
# Set up parallel processing with furrr
plan(multisession, workers = 2)

# Ensure reproducibility
options(future.seed = TRUE)

set.seed(123)

# Process models in parallel
data <- future_map(model_files, process_model_file, .progress = TRUE)

# Combine all results into a single data frame
data <- bind_rows(data)

saveRDS(data,here('data/processed/derived_quantities_envelope.rds')) 
