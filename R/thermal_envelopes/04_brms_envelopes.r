# libraries ---------------------------------------------------------------
library(brms)
library(here)
library(tidyverse)
library(broom)
library(ggdist)

# load raw data -----------------------------------------------------------

data <- readRDS(here("R/data/processed/derived_quantities_te.rds"))

# cleaning -----------------------------------------------------------

data <- data |> mutate(species = str_replace_all(species, " ", "_"), # clean species names 
                       sp_region = paste(species, region, sep = "_") 
)

# standardize and calculate anomalies -------------------------------------------------------------

data <- data |> group_by(region) |> mutate(
  year_c = (year - mean(year, na.rm = TRUE))/10 # center year
)  |> ungroup() |> 
  group_by(sp_region) |> # std index of abundance, and calculate species-specific anomalies
  mutate(
    index_std = scale(index)[, 1]) |> 
  mutate(across(
    c(cog_y, cog_x, depth_niche, thermal_niche),
    ~ .x - mean(.x, na.rm = TRUE),
    .names = "{.col}_c"
  )) |>  # calculated anomalies
  ungroup()

# Ensure cmdstanr is set up
cmdstanr::set_cmdstan_path()


# modeling ----------------------------------------------------------------

set.seed(123)

# Shared model settings, increase this for final run
model_control <- list(adapt_delta = 0.99, max_treedepth = 20)
n_chains <- 4
n_iter <- 4000
n_warmup <- 2000
n_cores <- 4


## priors ------------------------------------------------------------------


priors_custom <- c(
  # fixed slopes (global)
  prior(normal(0, 50), class = "b", coef = "year_c", resp = "cogyc"),
  prior(normal(0, 50), class = "b", coef = "year_c", resp = "cogxc"),
  
  prior(normal(0, 5), class = "b", coef = "year_c", resp = "depthnichec"),
  prior(normal(0, 0.5), class = "b", coef = "year_c", resp = "thermalnichec"),
  
  # Region-level priors for slopes
  prior(student_t(3, 0, 30), class = "sd", group = "region", coef = "year_c", resp = "cogyc"),
  prior(student_t(3, 0, 30), class = "sd", group = "region", coef = "year_c", resp = "cogxc"),
  prior(student_t(3, 0, 10), class = "sd", group = "region", coef = "year_c", resp = "depthnichec"),
  prior(student_t(3, 0, 0.2), class = "sd", group = "region", coef = "year_c", resp = "thermalnichec"),
  
  # Species nested in region priors for slopes
  prior(student_t(3, 0, 40), class = "sd", group = "region:species", coef = "year_c", resp = "cogyc"),
  prior(student_t(3, 0, 40), class = "sd", group = "region:species", coef = "year_c", resp = "cogxc"),
  prior(student_t(3, 0, 20), class = "sd", group = "region:species", coef = "year_c", resp = "depthnichec"),
  prior(student_t(3, 0, 0.4), class = "sd", group = "region:species", coef = "year_c", resp = "thermalnichec")
)



stud <- brm(
  formula = bf(
    cog_y_c |
      se(cog_y_se, sigma = TRUE) ~ 0 + year_c + (0 + year_c |
                                                   region) + (0 + year_c | ID | gr(region:species, by = region))
  ) +
    bf(
      cog_x_c |
        se(cog_x_se, sigma = TRUE) ~ 0 + year_c +  (0 + year_c  |
                                                      region) + (0 + year_c | ID | gr(region:species, by = region))
    ) +
    bf(
      depth_niche_c |
        se(depth_niche_se, sigma = TRUE) ~ 0 + year_c + (0 + year_c  |
                                                           region) + (0 + year_c | ID | gr(region:species, by = region))
    ) +
    bf(
      thermal_niche_c |
        se(thermal_niche_se, sigma = TRUE) ~ 0 + year_c + (0 + year_c |
                                                             region) + (0 + year_c | ID | gr(region:species, by = region))
    ) +
    set_rescor(FALSE), # remove residual correlation
  data = data,
  family = student(),
  chains = n_chains,
  iter = n_iter,
  seed = 123,
  warmup = n_warmup,
  cores = n_cores,
  control = model_control,
  prior = priors_custom,
  backend = "cmdstanr"
)

summary(stud)
coef(stud)
