# libraries ---------------------------------------------------------------

library(brms)
library(here)
library(tidyverse)

# Ensure cmdstanr is set up
cmdstanr::set_cmdstan_path()

# load data ---------------------------------------------------------------

data <- read_rds(here('data/meta_regression/processed/shift_data_processed.rds'))


# modeling ----------------------------------------------------------------
set.seed(123)

# Shared model settings, increase this for final run
model_control <- list(adapt_delta = 0.97, max_treedepth = 15)
n_chains <- 2
n_iter <- 2000
n_warmup <- 1000
n_cores <- 2

### We test 4 models, species nested in region and sp_region only, gaussian and student

m_gaussian <- brm(
  formula = bf(
    mvbind(cog_y_c, cog_x_c, depth_niche_c, thermal_niche_c) ~
      0 + year_c + (0 + year_c | p | sp_region)
  ) + set_rescor(FALSE),
  data = data,
  family = gaussian(),
  chains = n_chains,
  iter = n_iter,
  warmup = n_warmup,
  cores = n_cores,
  control = model_control,
  backend = "cmdstanr",
  file = here("meta_regression/model_fit/m_normal.rds")
)

m_student <- brm(
  formula = bf(
    mvbind(cog_y_c, cog_x_c, depth_niche_c, thermal_niche_c) ~
      0 + year_c + (0 + year_c | p | sp_region)
  ) + set_rescor(FALSE),
  data = data,
  family = student(),
  chains = n_chains,
  iter = n_iter,
  warmup = n_warmup,
  cores = n_cores,
  control = model_control,
  backend = "cmdstanr",
  file = here("meta_regression/model_fit/m_student.rds")
)
