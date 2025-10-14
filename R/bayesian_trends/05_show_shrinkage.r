# plot shrinkage ---------------------------------------------------------

library(broom)
library(purrr)


data <- readRDS(here("data/processed/derived_quantities.rds"))


data <- data |> mutate(thermal_niche_breadth = upr90_thermal - lwr10_thermal, depth_niche_breadth = upr90_depth - lwr10_depth, 
                       species = str_replace_all(species, " ", "_"), # clean species names 
                       region_species = paste(region,species, sep = "_") 
)


data <- data |> group_by(region) |> mutate(
  year_c = (year - mean(year, na.rm = TRUE))/10 # center year
)  |> ungroup() |> 
  group_by(region_species) |> # std index of abundance, and calculate species-specific anomalies
  mutate(
    index_std = scale(index)[, 1],
    eao_std = scale(eao)[, 1]) |> 
  mutate(across(
    c(cog_y, cog_x, depth_niche, thermal_niche,
      thermal_niche_min, thermal_niche_max,
      thermal_niche_breadth, depth_niche_breadth),
    ~ .x - mean(.x, na.rm = TRUE),
    .names = "{.col}_c"
  )) |>  # calculated anomalies
  ungroup()

# Outcomes to loop over
outcomes <- c("cog_y_c", "cog_x_c", "depth_niche_c", "thermal_niche_c")

# bayes_resp_map <- c(
#   cog_y_c = "cogyc_year_c",
#   cog_x_c = "cogxc_year_c",
#   depth_niche_c = "depthnichec_year_c",
#   thermal_niche_c = "thermalnichec_year_c"
# )

# Function to extract LM slopes
get_lm_slopes <- function(outcome) {
  data %>%
    group_by(region_species) %>%
    do(tidy(lm(as.formula(paste0(outcome, " ~ year_c")), data = .))) %>%
    filter(term == "year_c") %>%
    rename(mean_slope = estimate) %>%
    mutate(
      lower_95 = mean_slope - 1.96 * std.error,
      upper_95 = mean_slope + 1.96 * std.error,
      model = "lm",
      outcome = outcome
    ) %>%
    select(region_species, mean_slope, lower_95, upper_95, model, outcome)
}

# load Bayesian slopes
bayes_slopes <- readRDS(here("data/processed/bayesian_species_trends.rds")) %>% mutate(model = 'Bayesian GLMM')


# Combine all LM and Bayesian slopes
lm_slopes <- map_dfr(outcomes, function(outcome) {
  lm_part <- get_lm_slopes(outcome)
  #bayes_part <- get_bayes_slopes(outcome, bayes_resp_map[outcome])
  #bind_rows(lm_part, bayes_part)
})

all_slopes <- bind_rows(bayes_slopes %>% select(region_species,outcome,mean_slope,lower_95,upper_95,model), lm_slopes %>% select(region_species,outcome,mean_slope,lower_95,upper_95,model))

# Rename outcomes for prettier plot labels
all_slopes$outcome <- recode(
  all_slopes$outcome,
  cog_y_c = "Latitudinal centroid",
  cog_x_c = "Longitudinal centroid",
  depth_niche_c = "Depth niche",
  thermal_niche_c = "Thermal niche",
  cogyc = "Latitudinal centroid",
  cogxc = "Longitudinal centroid",
  depthnichec = "Depth niche",
  thermalnichec = "Thermal niche"
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

# Plot with position dodge to separate the points slightly
ggplot(all_slopes, aes(
  x = mean_slope,
  y = reorder(region_species, mean_slope),
  color = model
)) +
  geom_pointrange(
    aes(xmin = lower_95, xmax = upper_95),
    size = 0.2,
    position = position_dodge(width = 0.5),alpha = .8
  ) +
  facet_wrap(~outcome, scales = "free", ncol = 1) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = model_colors) +
  ggstats::geom_stripped_rows(
    data = all_slopes,
    aes(y = region_species),
    odd = "white", even = "grey90", alpha = .4, inherit.aes = FALSE
  ) +
  labs(
    y = 'Species',
    x = "Estimate of decade effect",
    color = "Model"
  ) +
  theme(
    axis.text.x = element_blank(),  # removes y-axis text (species names)
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank()
  ) +
  coord_flip()

# save plot 
ggsave(
  here("bayesian_trends/figures/supp/shrinkage_slopes.png"),
  width = 8, height = 10, dpi = 450
)