# load libraries ----------------------------------------------------------
library(tidyverse)
library(here)
library(sdmTMB)
library(furrr)


# define paths ------------------------------------------------------------
models_folder <- here("sdm_modeling/fitted")
model_files <- list.files(models_folder, pattern = "\\.rds$", full.names = TRUE)

# define function --------------------------------------------------
process_model_file <- function(model_file) {
  
  # load spatial grid -----------------------------------------------------
  grid <- readRDS(here("data/processed/prediction_grid.rds")) 
  # %>% mutate(month_f = as.factor(min_month)) %>% droplevels()
  
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
  
  # get the main survey
  main_survey <- model$data %>%
    count(survey, sort = TRUE) %>%
    slice(1) %>%
    pull(survey) %>%
    as.factor()
  
  # quarter order
  quarter_levels <- c("Q1", "Q2", "Q2_Q3", "Q3", "Q3_Q4", "Q4")
  
  # get minimum quarter for this region
  min_quarter <- model$data %>%
    filter(region_short == this_region) %>%
    mutate(quarter = factor(quarter, levels = quarter_levels, ordered = TRUE)) %>%
    pull(quarter) %>%
    na.omit() %>%
    min()
  
  # standardize depth
  region_grid <- grid %>%
    filter(region_short == this_region, year %in% years) %>%
    mutate(
      logdepth = (log(depth) - mean(log(model$data$depth), na.rm = TRUE)) / 
        sd(log(model$data$depth), na.rm = TRUE), # scale by data
      logdepth2 = logdepth^2,
      survey = main_survey,
      area = 16,
      quarter = min_quarter
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
  
  depth_niche <- get_weighted_average(pred, vector = pred$data$depth, area = 16) %>%
    transmute(year, depth_niche = est, depth_niche_se = se)
  
  thermal_niche <- get_weighted_average(pred, vector = pred$data$mean_temp, area = 16) %>%
    transmute(year, thermal_niche = est, thermal_niche_se = se)
  
  # combine results -------------------------------------------------------
  out <- index %>%
    inner_join(cog, by = "year") %>%
    inner_join(eao, by = "year") %>% inner_join(depth_niche, by = "year") %>%
    inner_join(thermal_niche, by = "year")
  
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
options(future.seed = TRUE)

set.seed(42)

# Process models in parallel
data <- future_map(model_files, process_model_file, .progress = TRUE)

# Combine all results into a single data frame
data <- bind_rows(data)

saveRDS(data,here('data/processed/derived_quantities.rds')) 


# pred$data$biomass_threshold <- quantile(exp(pred$data$est), probs = 0.05, na.rm = TRUE)
# 
# pred$data <- pred$data %>%
#   mutate(biomass = exp(est),
#          occupied = biomass > biomass_threshold)
# 
# # 3. Calculate weighted means (already in your code, so keep)
# pred$data <- pred$data %>%
#   group_by(year) %>%
#   mutate(
#     weighted_depth = sum(depth * biomass) / sum(biomass),
#     weighted_temp_mean = sum(mean_bottomT * biomass) / sum(biomass),
#     weighted_temp_min = sum(coldest_monthT * biomass) / sum(biomass),
#     weighted_temp_max = sum(warmest_monthT * biomass) / sum(biomass),
#     weighted_q90_temp = Hmisc::wtd.quantile(warmest_monthT, weights = biomass, probs = 0.9),
#     weighted_q90_depth = Hmisc::wtd.quantile(depth, weights = biomass, probs = 0.9),
#     weighted_q10_depth = Hmisc::wtd.quantile(depth, weights = biomass, probs = 0.1),# need for interquartile range
#     weighted_q10_temp = Hmisc::wtd.quantile(coldest_monthT, weights = biomass, probs = 0.1)
#   ) 
# 
# # 4. Calculate extent of occurrence
# occupied_area <- pred$data %>%
#   group_by(year) %>%
#   summarise(
#     occupied_area = sum(area[occupied], na.rm = TRUE),  # area is 25 km² per cell
#     log_extent = log(occupied_area)
#   ) |> select(-occupied_area)
# 
# # 5. range edges
# edges <- pred$data %>%
#   group_by(year) %>%
#   # Sort by Y before cumulative biomass calc
#   arrange(Y) %>%
#   mutate(
#     cum_biomass = cumsum(biomass) / sum(biomass)
#   ) %>%
#   summarise(
#     trailing_edge_Y = approx(cum_biomass, Y, xout = 0.01, ties = mean)$y,
#     leading_edge_Y = approx(cum_biomass, Y, xout = 0.99, ties = mean)$y,
#     .groups = "drop"
#   )
# 
# minY <- min(pred$data$Y, na.rm = TRUE)
# maxY <- max(pred$data$Y, na.rm = TRUE)
# buffer <- 0.10 * (maxY - minY)
# 
# lower_bound <- minY + buffer
# upper_bound <- maxY - buffer
# 
# # Compute mean edge positions
# mean_trailing <- mean(edges$trailing_edge_Y, na.rm = TRUE)
# mean_leading <- mean(edges$leading_edge_Y, na.rm = TRUE)
# 
# # Check if edges are too close to boundaries
# edges <- edges %>%
#   mutate(
#     trailing_edge_Y = if (mean_trailing < lower_bound) NA else trailing_edge_Y,
#     leading_edge_Y = if (mean_leading > upper_bound) NA else leading_edge_Y
#   )
# 
# 
# 
# derived_quantities = derived_quantities |>left_join(occupied_area) |> left_join(distinct(pred$data,year,weighted_depth,weighted_temp_mean,weighted_temp_min,weighted_temp_max,weighted_q90_temp, weighted_q90_depth, weighted_q10_temp,weighted_q10_depth )) |> left_join(edges)
# 
# # now get the core area temperature and so on
# 
# pred$data$id_cell = paste(pred$data$X,pred$data$Y,sep='_'
# )
# 
# core_area_cells <- pred$data %>%
#   group_by(year) %>%
#   arrange(desc(biomass), .by_group = TRUE) %>%
#   mutate(
#     total_biomass_year = sum(biomass, na.rm = TRUE),
#     cumulative_biomass_year = cumsum(biomass) / total_biomass_year
#   ) %>%
#   filter(cumulative_biomass_year <= 0.95) %>%
#   ungroup() %>%
#   group_by(id_cell) %>%
#   summarise(core_presence = n(), total_years = n_distinct(year), .groups = "drop") %>%
#   mutate(weighted_core_score = core_presence / total_years) %>%
#   filter(weighted_core_score > 0.5) %>%
#   pull(id_cell)
# 
# out <- derived_quantities |> left_join(
#   pred$data %>%
#     filter(id_cell %in% core_area_cells) %>%
#     group_by(year) %>%
#     summarise(
#       mean_bottomT = mean(mean_bottomT, na.rm = TRUE),
#       coldest_monthT = mean(coldest_monthT, na.rm = TRUE),
#       warmest_monthT = mean(warmest_monthT, na.rm = TRUE)
#     ), by = "year"
# ) |> mutate(species = this_species, region = this_region)
# 
# gc()  # Final garbage collection before returning

