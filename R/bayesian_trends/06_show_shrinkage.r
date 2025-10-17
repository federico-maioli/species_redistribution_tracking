# plot shrinkage ---------------------------------------------------------
library(here)
library(broom)
library(purrr)
library(tidyverse)

data <- readRDS(here("R/data/processed/derived_quantities.rds"))


data <- data |> mutate(species = str_replace_all(species, " ", "_"), # clean species names 
                       sp_region = paste(species, region, sep = "_") 
)

# center data

data <- data |> group_by(region) |> mutate(
  year_c = (year - mean(year, na.rm = TRUE))/10 # center year
)  |> ungroup() |> 
  group_by(sp_region) |> # std index of abundance, and calculate species-specific anomalies
  mutate(
    index_std = scale(index)[, 1],
    eao_std = scale(eao)[, 1]) |> 
  mutate(across(
    c(cog_y, cog_x, depth_niche, thermal_niche),
    ~ .x - mean(.x, na.rm = TRUE),
    .names = "{.col}_c"
  )) |>  # calculated anomalies
  ungroup()


# outcomes to loop over
outcomes <- c("cog_y_c", "cog_x_c", "depth_niche_c", "thermal_niche_c")


# Function to extract LM slopes
get_lm_slopes <- function(outcome) {
  data %>%
    group_by(sp_region) %>%
    do(broom::tidy(lm(as.formula(paste0(outcome, " ~ 0 + year_c")), data = .))) %>%
    filter(term == "year_c") %>%
    rename(mean_slope = estimate) %>%
    mutate(
      lower_95 = mean_slope - 1.96 * std.error,
      upper_95 = mean_slope + 1.96 * std.error,
      model = "lm",
      outcome = outcome
    ) %>%
    select(sp_region, mean_slope, lower_95, upper_95, model, outcome)
}

# load Bayesian slopes
bayes_slopes <- readRDS(here('R/data/processed/bayesian_species_trends.rds')) %>% mutate(model = 'Bayesian GLMM') %>% mutate(sp_region = paste(species, region, sep = "_"))

# Combine all LM and Bayesian slopes
lm_slopes <- map_dfr(outcomes, function(outcome) {
  lm_part <- get_lm_slopes(outcome)
  #bayes_part <- get_bayes_slopes(outcome, bayes_resp_map[outcome])
  #bind_rows(lm_part, bayes_part)
})

all_slopes <- bind_rows(bayes_slopes %>% select(sp_region,outcome,mean_slope,lower_95,upper_95,model), lm_slopes %>% select(sp_region,outcome,mean_slope,lower_95,upper_95,model))

# Rename outcomes for prettier plot labels
all_slopes$outcome <- recode(
  all_slopes$outcome,
  lat_shift = "Latitudinal centroid",
  lon_shift = "Longitudinal centroid",
  depth_shift = "Depth niche",
  thermal_shift = "Thermal niche",
  cog_y_c = "Latitudinal centroid",
  cog_x_c = "Longitudinal centroid",
  depth_niche_c = "Depth niche",
  thermal_niche_c = "Thermal niche"
)

# force order
all_slopes$outcome <- factor(
  all_slopes$outcome,
  levels = c("Latitudinal centroid",
             "Longitudinal centroid",
             "Depth niche",
             "Thermal niche")
)

model_colors <- c("lm" = "#2c7fb8", "Bayesian GLMM" = "#f03b20")  # blue & orange

set.seed(1)  # for reproducibility
selected_species <- all_slopes %>%
  distinct(sp_region) %>%        # get unique species
  sample_n(50) %>%             # randomly pick 50 species
  pull(sp_region)                # extract as vector

# Filter the dataset for these species
subsampled_slopes <- all_slopes %>%
  filter(sp_region %in% selected_species)

# Plot with position dodge to separate the points slightly
ggplot(subsampled_slopes, aes(
  x = mean_slope,
  y = reorder(sp_region, mean_slope),
  color = model,
  group = sp_region
)) +
  geom_pointrange(
    aes(xmin = lower_95, xmax = upper_95),
    size = 0.1,
    position = position_dodge2(width = 0.5, padding = 0.2),
    alpha = 0.9
  ) +
  facet_wrap(~outcome, scales = "free", ncol = 1) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = model_colors) +
  ggstats::geom_stripped_rows(
   data = subsampled_slopes,
   aes(y = sp_region),
   odd = "white", even = "grey90", alpha = .4, inherit.aes = FALSE
  ) +
  labs(
    y = 'Species',
    x = "Estimate of decadal change",
    color = "Model"
  ) + geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  coord_flip()

# save plot 
ggsave(
  here('output/figures/supp/shrinkage_supp.png'),
  width = 180,
  height = 190,
  dpi = 600,
  units = "mm",
  bg = "white"
)

