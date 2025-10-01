library(here)
library(tidyverse)
library(sdmTMB)
library(lubridate)
library(future)
library(furrr)
library(purrr)

# load data and rename column
data <- read_rds(here('data/processed/fishglob_clean.rds')) |> 
  rename(species = accepted_name)

# add quarter, month_f (factor), and survey (factor) columns
data <- data %>% 
  mutate(
    month = as.numeric(month),
    quarter = as.factor(case_when(
      month %in% 1:3 ~ 1,
      month %in% 4:6 ~ 2,
      month %in% 7:9 ~ 3,
      month %in% 10:12 ~ 4
    )),
    month_f = as.factor(month),
    survey = as.factor(survey)
  )


# set up parallel processing with 2 workers
plan(multisession, workers = 2)

# create a table of unique species-region combinations
species_region_table <- data |> distinct(species, region_short)
# Function to process each species-region subset
process_species_region <- function(i) {
  tryCatch({
    this_species <- species_region_table$species[i]
    this_region <- species_region_table$region_short[i]
    
    message("Processing: ", this_species, " - ", this_region)
    
    # Filter data for current species and region, drop unused factor levels
    sub <- data %>%
      filter(species == this_species, region_short == this_region, !is.na(depth), !is.na(mean_temp)) %>%
      droplevels()
    
    # Skip if no data available
    if (nrow(sub) == 0) return(NULL)
    
    # Standarize variables
    sub <- sub %>%
      mutate(
        log_depth_std = scale(log(depth))[, 1],
        mean_temp_std = scale(mean_temp)[, 1],
        warmest_temp_std = scale(warmest_temp)[, 1],
        coldest_temp_std = scale(coldest_temp)[, 1]
      )
    
    # Check if multiple surveys and sufficient months are available
    has_multiple_surveys <- n_distinct(sub$survey) > 1
    has_enough_months <- n_distinct(sub$month_f) >= 3
    
    # Define base formula
    formula <- wgt_cpua ~ 0 + as.factor(year) + s(mean_temp_std, k = 3) + s(warmest_temp_std, k = 3)
    
    # Add random effects to formula based on data structure
    if (has_multiple_surveys && has_enough_months) {
      formula <- update(formula, . ~ . + as.factor(survey) + (1 | month_f))
    } else if (has_multiple_surveys) {
      formula <- update(formula, . ~ . + as.factor(survey))
    } else if (has_enough_months) {
      formula <- update(formula, . ~ . + (1 | month_f))
    }
    
    # Fit the spatial-temporal model using sdmTMB
    fit <- sdmTMB(
      formula = formula,
      family = tweedie(link = "log"),
      data = sub,
      time = 'year',
      spatiotemporal = 'off',
      spatial = "off"
    )
    
    # Run extra optimization if gradient is large
    if (max(fit$gradients) > 0.001) {
      fit <- sdmTMB::run_extra_optimization(fit, nlminb_loops = 1L, newton_loops = 1L)
    }
    
    # Perform sanity check on model fit
    sanity_check <- sanity(fit, silent = TRUE)
    if (!sanity_check$hessian_ok) {
      message("Skipping: ", this_species, " - ", this_region, " (Hessian check failed)")
      return(NULL)
    }
    
    
    # Save the fitted model object
    saveRDS(fit, file = here("climate_velocity/fitted", paste0(this_region, "_", this_species, "_thermal_envelope.rds")))
    
    message("Completed: ", this_species, " - ", this_region)
    return(NULL)
    
  }, error = function(e) {
    message("Error in: ", this_species, " - ", this_region, " -> ", e$message)
    return(NULL)
  })
}

# Run the process in parallel over all species-region combos with progress bar
furrr::future_map(1:nrow(species_region_table), process_species_region, .progress = TRUE)
