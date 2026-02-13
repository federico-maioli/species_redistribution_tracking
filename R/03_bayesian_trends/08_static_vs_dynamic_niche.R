library(here)
library(tidyverse)
library(broom)
library(ggdist)
library(purrr)

# load raw data -----------------------------------------------------------

data <- readRDS(here("R/data/processed/derived_quantities_sdm.rds"))

# cleaning -----------------------------------------------------------

data <- data |> mutate(species = str_replace_all(species, " ", "_"), # clean species names 
                       sp_region = paste(species, region, sep = "_") 
)

# standardize and calculate anomalies -------------------------------------------------------------

data <- data |> group_by(region) |> mutate(
  year_c = (year - mean(year, na.rm = TRUE))/10 # center year
)  |> ungroup() |> 
  group_by(sp_region) |> # std index of abundance, and calculate species-specific anomalies
  # mutate(
  #   index_std = scale(index)[, 1]) |> 
  mutate(across(
    c(cog_y, cog_x, depth_niche, thermal_niche,thermal_niche_constant_density),
    ~ .x - mean(.x, na.rm = TRUE),
    .names = "{.col}_c"
  )) |>  # calculated anomalies
  ungroup()

slopes <- data %>%
  group_by(sp_region,region) %>%
  nest() %>%
  mutate(
    mean_realized = map(data, ~ lm(thermal_niche_c ~ 0 + year_c, data = .x)),
    mean_static   = map(data, ~ lm(thermal_niche_constant_density_c ~ 0 + year_c, data = .x)),
    beta_realized = map_dbl(mean_realized, ~ coef(.x)["year_c"]),
    beta_static   = map_dbl(mean_static,   ~ coef(.x)["year_c"]),
    delta_beta    = beta_realized - beta_static
  )

region_summary <- slopes %>%
  group_by(region) %>%
  summarise(
    n_species      = n(),
    mean_realized  = mean(beta_realized, na.rm = TRUE),
    mean_static    = mean(beta_static, na.rm = TRUE),
    mean_delta     = mean(beta_realized - beta_static, na.rm = TRUE),
    prop_buffering = mean((beta_realized - beta_static) < 0, na.rm = TRUE),
    cor_slopes     = cor(beta_realized, beta_static, 
                         use = "complete.obs"),
    .groups = "drop"
  )

ggplot(region_summary, 
       aes(x = mean_static, y = mean_realized, label = region)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 3) +
  geom_text(nudge_y = 0.03, size = 3) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6)
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6)) +
  labs(
    x = expression("Static thermal-niche trend (" * degree * "C decade"^{-1} * ")"),
    y = expression("Realized thermal-niche trend (" * degree * "C decade"^{-1} * ")"),
  ) +
  theme_light()

ggsave(
  here('output/figures/supp/static_vs_dynamic_niche.png'),
  #plot = pp_plot,  
  width = 180,  
  height = 110,
  dpi = 600,
  units = "mm"
)
