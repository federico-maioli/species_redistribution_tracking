# libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(sdmTMB)
library(lubridate)
library(future)
library(furrr)
library(purrr)

# load data and rename column
data <- read_rds(here('data/trawl_surveys/fishglob_clean.rds')) |> 
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
      filter(species == this_species, region_short == this_region, !is.na(depth)) %>%
      droplevels()
    
    # Skip if no data available
    if (nrow(sub) == 0) return(NULL)
    
    # Create transformed depth variables for modeling
    sub <- sub %>%
      mutate(
        logdepth = scale(log(depth))[, 1],
        logdepth2 = logdepth^2
      )
    
    # Add UTM coordinates (assuming add_utm_columns() is defined elsewhere)
    sub <- add_utm_columns(sub, ll_names = c("longitude", "latitude"), units = "km")
    
    # Create non-convex hull boundary for mesh construction
    bnd <- INLA::inla.nonconvex.hull(cbind(sub$X, sub$Y), convex = -0.05)
    
    # Build 2D mesh with custom parameters
    inla_mesh <- fmesher::fm_mesh_2d_inla(
      boundary = bnd,
      loc = cbind(sub$X, sub$Y),
      max.edge = c(80, 200),
      offset = c(-0.1, 50),
      cutoff = 20,
      min.angle = 5 # from Eric's paper
    )
    
    # Create SPDE mesh for spatial modeling
    spde <- make_mesh(sub, c("X", "Y"), mesh = inla_mesh)
    
    # Choose spatiotemporal model based on region
    st_type <- ifelse(this_region == "BC", "rw", "iid")
    
    # Check if multiple surveys and sufficient months are available
    has_multiple_surveys <- n_distinct(sub$survey) > 1
    has_enough_months <- n_distinct(sub$month_f) >= 3
    
    # Define base formula
    formula <- wgt_cpua ~ 0 + as.factor(year) + logdepth + logdepth2
    
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
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = st_type
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
    saveRDS(fit, file = here("SDM/fitted", paste0(this_region, "_", this_species, "_sdm.rds")))
    
    message("Completed: ", this_species, " - ", this_region)
    return(NULL)
    
  }, error = function(e) {
    message("Error in: ", this_species, " - ", this_region, " -> ", e$message)
    return(NULL)
  })
}

# Run the process in parallel over all species-region combos with progress bar
future_map(1:nrow(species_region_table), process_species_region, .progress = TRUE)
