# load libraries ----------------------------------------------------------
library(tidyverse)
library(here)
library(sdmTMB)
library(furrr)


# define paths ------------------------------------------------------------
models_folder <- here("R/02_sdm_modeling/fitted")
model_files <- list.files(models_folder, pattern = "\\.rds$", full.names = TRUE)

# define function --------------------------------------------------
process_model_file <- function(model_file) {
  
  # load spatial grid -----------------------------------------------------
  grid <- readRDS(here("R/data/processed/prediction_grid.rds")) 
  # %>% mutate(month_f = as.factor(min_month)) %>% droplevels()
  
  # parse model metadata --------------------------------------------------
  file_name <- basename(model_file)
  parts <- strsplit(file_name, "_")[[1]]
  this_region <- parts[1]
  this_species <- parts[2]
  
  # load and prepare model ------------------------------------------------
  model <- readRDS(model_file)
  model <- sdmTMB:::reload_model(model)
  # Determine temporal domain for prediction
  if (this_region %in% c("BC", "GOA")) {
    # Continuous RW — predict full span
    years <- seq(min(model$data$year), max(model$data$year))
  } else {
    # Only predict observed years
    years <- sort(unique(model$data$year))
  }
  
  # prepare region-specific grid ------------------------------------------
  
  # standardize depth
  region_grid <- grid %>%
    filter(region_short == this_region, year %in% years) %>%
    mutate(
      logdepth = (log(depth) - mean(log(model$data$depth), na.rm = TRUE)) / 
        sd(log(model$data$depth), na.rm = TRUE), # scale by data
      logdepth2 = logdepth^2,
      area = 16,
      cell_id = paste(X,Y,sep='_')
    ) %>%
    #filter(!is.na(month_f)) %>%
    na.omit() %>%
    droplevels()
  
  # if something goes wrong here, skip
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
  # index <- get_index(pred, area = 16, bias_correct = TRUE) %>%
  #   transmute(year, index = log_est, index_se = se)
  
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
  

# thermal niche with constant density -------------------------------------

  # Extract prediction data
  pdat <- pred$data %>%
    mutate(
      biomass = exp(est) * area
    )
  
  # Compute average biomass per cell across all years
  pdat <- pdat %>%
    group_by(cell_id) %>%
    mutate(avg_density = mean(biomass, na.rm = TRUE)) %>%
    ungroup()
  
  # Compute yearly thermal niche using static density
  thermal_niche_constant <- pdat %>%
    group_by(year) %>%
    summarise(
      thermal_niche_constant_density =
        sum(mean_temp * avg_density, na.rm = TRUE) /
        sum(avg_density, na.rm = TRUE),
      .groups = "drop"
    )
  
  # combine results -------------------------------------------------------
  out <- cog %>% inner_join(depth_niche, by = "year") %>%
    inner_join(thermal_niche, by = "year") %>% inner_join(thermal_niche_constant, by = "year")
  
  # no longer needed
  # niche <- pred$data %>% mutate(biomass = exp(est + log(pred$data$area))) |> # area expansion
  #   group_by(year) %>%
  #   summarise(
  #     depth_niche = sum(depth * biomass) / sum(biomass),
  #     thermal_niche = sum(mean_temp * biomass) / sum(biomass),
  #     thermal_niche_min = sum(coldest_temp * biomass) / sum(biomass),
  #     thermal_niche_max = sum(warmest_temp * biomass) / sum(biomass),
  #     upr90_thermal = Hmisc::wtd.quantile(mean_temp, weights = biomass, probs = 0.9), # remember to reduce dependencies
  #     lwr10_thermal = Hmisc::wtd.quantile(mean_temp, weights = biomass, probs = 0.1),
  #     upr90_depth = Hmisc::wtd.quantile(depth, weights = biomass, probs = 0.9),
  #     lwr10_depth = Hmisc::wtd.quantile(depth, weights = biomass, probs = 0.1),
  #   ) 
  
  out <- out %>% mutate(species = this_species, region = this_region)
  
  
  return(out)
}

# process models in parallel ----------------------------------------------
# Set up parallel processing with furrr
plan(multisession, workers = 2)

# Ensure reproducibility
#set.seed(42)

# Process models in parallel
data <- future_map(model_files, process_model_file, .progress = TRUE, .options = furrr_options(seed = 42))

# Combine all results into a single data frame
data <- bind_rows(data)

saveRDS(data,here('R/data/processed/derived_quantities_sdm.rds')) 
