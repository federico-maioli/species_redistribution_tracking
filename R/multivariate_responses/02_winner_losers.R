# libraries
library(broom)
library(brms)
library(here)
library(patchwork)
library(tidyverse)
library(cowplot)
library(tidybayes)

# load slopes -----------------------------------------------------------

slopes <- read_rds(here('R/data/processed/bayesian_species_trends.rds'))

# load pca scores ---------------------------------------------------------

pc_scores <- read_rds(here('R/data/processed/pca_scores.rds'))

# read in data ------------------------------------------------------------
data = read_rds(here('R/data/processed/derived_quantities.rds'))

data <- data |> mutate(species = str_replace_all(species, " ", "_"), # clean species names 
                       sp_region = paste(species, region, sep = "_") 
)

data <- data |> group_by(region) |> mutate(
  year_c = (year - mean(year, na.rm = TRUE))/10 # center year
)  |> ungroup() |> 
  group_by(sp_region) |> # std index of abundance, and calculate species-specific anomalies
  mutate(
    index_c = index - mean(index)) |> 
  mutate(across(
    c(cog_y, cog_x, depth_niche, thermal_niche),
    ~ .x - mean(.x, na.rm = TRUE),
    .names = "{.col}_c"
  )) |>  # calculated anomalies
  ungroup()

ggplot(data %>% filter(sp_region=='Centropristis_philadelphica_GMX'),aes(year,index))+geom_point()

# compute abundance trends ------------------------------------------------

abundance_trends <- data %>%
  group_by(species,region) %>%
  do(model = lm(index ~ 1 + year_c, data = .)) %>%
  mutate(slope = coef(model)[["year_c"]])

# Ensure cmdstanr is set up
cmdstanr::set_cmdstan_path()

fit <- brm(
  index_c ~ 0 + year_c + (0 + year_c | region:species),
  data = data,
  prior = c(
    prior(normal(0, 0.3), class = "b", coef = "year_c"),  # strong prior on slope
    prior(exponential(1), class = "sd")                    # reasonable random-effect scale
  ),
  cores = 4, iter = 4000, control = list(adapt_delta = 0.95),
  backend = "cmdstanr"
)

posterior_slopes <- fit %>%
  spread_draws(
    b_year_c,                       # global fixed slope
    `r_region:species`[species_region, year_c]  # species × region random slopes
  )


species_slopes <- posterior_slopes %>% mutate(
    slope_species = b_year_c + `r_region:species`
  )

species_summary <- species_slopes %>%
  group_by(species_region) %>%
  summarise(
    slope_median = median(slope_species),
    slope_mean   = mean(slope_species),
    slope_lower  = quantile(slope_species, 0.025),
    slope_upper  = quantile(slope_species, 0.975),
    .groups = "drop"
  )

# merge datasets ----------------------------------------------------------

win_lose <- slopes %>% select(species, region, outcome, mean_slope) %>% left_join(abundance_trends %>% mutate(index = slope) %>% select(species,region,index)) %>% left_join(pc_scores %>% select(Dim.1,species,region))


# winner and losers -----------------------------------------------------


# regression stats for annotation
# Correct label expressions (with balanced parentheses)
outcomes_labels <- list(
  "lat_shift"     = expression("Latitudinal shift (km decade"^-1*")"),
  "lon_shift"     = expression("Longitudinal shift (km decade"^-1*")"),
  "depth_shift"   = expression("Depth shift (m decade"^-1*")"),
  "thermal_shift" = expression("Thermal niche shift ("*degree*C~decade^-1*")")
)

# Colors
predictor_colors <- c(
  "lat_shift" = "#762A83",
  "lon_shift" = "#543005",
  "depth_shift" = "#8C510A",
  "thermal_shift" = "#2166AC",
  "Dim.1" = "#1B9E77"
)

# Compute lm_stats safely
lm_stats <- win_lose %>%
  group_by(outcome) %>%
  summarise(model = list(lm(index ~ mean_slope)), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    slope = coef(model)[[2]],
    r2 = summary(model)$r.squared,
    p_value = summary(model)$coefficients[2, 4]
  ) %>%
  ungroup() %>%
  mutate(
    label = paste0(
      "Slope = ", round(slope, 3),
      "\nR² = ", round(r2, 2),
      "\np = ", signif(p_value, 2)
    )
  )

# The plotting function
make_plot <- function(var_name) {
  df <- win_lose %>% filter(outcome == var_name)
  stats <- lm_stats %>% filter(outcome == var_name)
  
  ggplot(df, aes(x = mean_slope, y = index)) +
    geom_point(alpha = 0.3, color = predictor_colors[var_name], size = 2) +
    geom_smooth(
      method = "lm", se = TRUE,
      color = predictor_colors[var_name],
      fill = predictor_colors[var_name],
      alpha = 0.15
    ) +
    geom_text(
      data = stats,
      aes(x = Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = 1.1, vjust = 1.1,
      size = 3, fontface = "bold"
    ) +
    labs(
      x = outcomes_labels[[var_name]],
      y =  expression("Change in log(abundance) per decade")
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 11, face = "bold"),
      axis.title.y = element_text(size = 11, face = "bold")
    )
}

# Make all four plots
p1 <- make_plot("lat_shift")
p2 <- make_plot("lon_shift")
p3 <- make_plot("depth_shift")
p4 <- make_plot("thermal_shift")
p5 <- ggplot() + theme_void()

p_top <- p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 5, axes = "collect") + theme(axis.title.y = element_text(size = 9), axis.title.x = element_text(size = 9))


# now by region -----------------------------------------------------

# add dim 1 

win_lose <- win_lose %>%
# create a new row for Dim.1
bind_rows(
  win_lose %>%
    select(species, region, index, Dim.1) %>%
    mutate(
      outcome = "Dim.1",
      mean_slope = Dim.1
    ) %>%
    select(names(win_lose))  # make column order consistent
)

models_by_region <- win_lose %>%
  group_by(region, outcome) %>%      # real groups that exist in data
  filter(!is.na(index), !is.na(mean_slope)) %>%
  tidyr::nest() %>%
  mutate(
    model  = map(data, ~ lm(index ~ mean_slope, data = .x)),
    tidy   = map(model, broom::tidy, conf.int = TRUE, conf.level = 0.95),
    glance = map(model, broom::glance),
    n      = map_int(data, nrow)
  )

coef_df <- models_by_region %>%
  tidyr::unnest(cols = tidy) %>%
  filter(term == "mean_slope") %>%
  select(region, outcome, estimate, std.error, conf.low, conf.high, statistic, p.value, n) %>%
  left_join(
    models_by_region %>%
      tidyr::unnest(cols = glance) %>%
      select(region, outcome, r.squared),
    by = c("region","outcome")
  ) %>%
  mutate(
    sig = if_else(p.value < .05, "sig", "ns")
  ) %>%
  mutate(
    est_pct = (exp(estimate) - 1) * 100,
    conf.low_pct = (exp(conf.low) - 1) * 100,
    conf.high_pct = (exp(conf.high) - 1) * 100
  )


# bottom row plot slopes --------------------------------------------------

# nice math-friendly x-axis labels for slope plot
predictor_axis_labels <- list(
  lat_shift = expression(atop(Delta * " Abundance decade"^-1,
                          "with a km decade"^-1 * " latitudinal shift")),
  lon_shift = expression(atop(Delta * " Abundance decade"^-1,
                          "with a km decade"^-1 * " longitudinal shift")),
  depth_shift = expression(atop(Delta * " Abundance decade"^-1,
                                "with a m decade"^-1 * " depth shift")),
  thermal_shift = expression(atop(Delta * " Abundance decade"^-1,
                                  "with a " * degree * "C decade"^-1 * " thermal niche shift")),
  Dim.1 = expression(atop(Delta * " Abundance decade"^-1,
                          "with a unit PC1 change"))
)
  
# plot coefs
make_coef_plot <- function(outcome_name) {
  df <- coef_df %>% filter(outcome ==   outcome_name)

  ggplot(df, aes(x = est_pct, y = region)) +
    ggstats::geom_stripped_rows(
      data = df,
      aes(y = region),
      odd = "white", even = "grey90", alpha = .2, inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_pointrange(
      aes(xmin = conf.low_pct,
          xmax = conf.high_pct,
          alpha = sig),
      color = predictor_colors[[outcome_name]],  # outcome-based colors
      size = 0.8
    ) +
    scale_alpha_manual(values = c("sig" = 1, "ns" = 0.4)) +
    labs(x = predictor_axis_labels[[outcome_name]],y=NULL) +
    theme_bw(base_size = 10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 11),
          legend.position = "none")
}

# make coefficient plots for all predictors
c1 <- make_coef_plot("lat_shift")
c2 <- make_coef_plot("lon_shift")
c3 <- make_coef_plot("depth_shift")
c4 <- make_coef_plot("thermal_shift")
c5 <- make_coef_plot("Dim.1")

p_bottom <- c1 + c2 + c3 + c4 + c5 + plot_layout(ncol = 5, axes = 'collect')

p_top / p_bottom + plot_layout(heights = c(1, 1)) 

cowplot::plot_grid(p_top, p_bottom, ncol = 1, align = "v", rel_heights = c(1,1 ), axis = "ll", labels = c("a)", "b)"),
                   label_size = 12)

ggsave(here('multivariate_responses/figures/win_lose.png'),   width = 14,    # inches, matches \textwidth in Overleaf
       height = 10,   # adjust for aspect ratio
       dpi = 600, units = 'in')

