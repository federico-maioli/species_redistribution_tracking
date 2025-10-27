# libraries --------------------------------------------------------------
library(here)
library(tidyverse)
library(sdmTMB)
library(lubridate)
library(future)
library(furrr)
library(purrr)

# load data and rename column
data <- read_rds(here('R/data/processed/fishglob_clean.rds'))

options(future.globals.maxSize = 1024 * 1024^2)  # 1 GB

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
      filter(species == this_species, region_short == this_region, !is.na(depth)) %>%
      droplevels()
    
    # Skip if no data available
    if (nrow(sub) == 0) return(NULL)
    
    # scale
    sub <- sub %>%
     mutate(
       warmest_temp_std = scale(warmest_temp)[, 1],
       mean_temp_std = scale(mean_temp)[,1]
     )
    
    # Add UTM coordinates (assuming add_utm_columns() is defined elsewhere)
    sub <- add_utm_columns(sub, ll_names = c("longitude", "latitude"), units = "km")
    
    # Create non-convex hull boundary for mesh construction
    #bnd <- INLA::inla.nonconvex.hull(cbind(sub$X, sub$Y), convex = -0.05)
    bnd <- fmesher::fm_nonconvex_hull(cbind(sub$X, sub$Y), convex = -0.05,format = "fm")
    
    # Build 2D mesh with custom parameters
    inla_mesh <- fmesher::fm_mesh_2d_inla(
      boundary = bnd,
      loc = cbind(sub$X, sub$Y),
      max.edge = c(80, 200),
      offset = c(-0.1, 50),
      cutoff = 20, # maximum 20 km between vertices
      min.angle = 5 
    )
    
    # Create SPDE mesh for spatial modeling
    spde <- make_mesh(sub, c("X", "Y"), mesh = inla_mesh)
    
    # Check if multiple surveys and sufficient months are available
    multi_survey_regions <- c("CBS",'NEUS-SS')
    multi_season_regions <- c("BAL", "NS", "NEUS-SS", "CBS", "GOM") # Identify whether the current region should use quarter
    use_quarter <- this_region %in% multi_season_regions
    use_survey <- this_region %in% multi_survey_regions
    
    # Define base formula
    formula <- kg_km2 ~ 0 + as.factor(year) + s(mean_temp_std, k=3) + s(warmest_temp_std, k=3)
    
    # Add random effects to formula based on data structure
    if (use_survey && use_quarter) {
      formula <- update(formula, . ~ . + as.factor(survey) + as.factor(quarter))
    } else if (use_survey) {
      formula <- update(formula, . ~ . + as.factor(survey))
    } else if (use_quarter) {
      formula <- update(formula, . ~ . + as.factor(quarter))
    }
    
    
    # Fit the spatial-temporal model using sdmTMB
    fit <- sdmTMB(
      formula = formula,
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      spatial = "on",
      spatiotemporal = 'off'
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
    saveRDS(fit, file = here("R/thermal_envelopes/fitted", paste0(this_region, "_", this_species, "_te.rds")))
    
    message("Completed: ", this_species, " - ", this_region)
    return(NULL)
    
  }, error = function(e) {
    message("Error in: ", this_species, " - ", this_region, " -> ", e$message)
    return(NULL)
  })
}

# Run the process in parallel over all species-region combos with progress bar
future_map(1:nrow(species_region_table), process_species_region, .progress = TRUE)

#process_species_region(2)
