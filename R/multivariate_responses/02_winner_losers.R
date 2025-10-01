# libraries
library(broom)
library(brms)
library(here)
library(patchwork)
library(tidyverse)
library(cowplot)

# load slopes -----------------------------------------------------------

slopes <- read_rds(here('data/processed/bayesian_species_trends.rds'))


# load oca scores ---------------------------------------------------------

pc_scores <- read_rds(here('data/processed/pca_scores.rds'))

# read in data ------------------------------------------------------------
data = read_rds(here('data/processed/derived_quantities.rds'))

data <- data |> mutate(thermal_niche_breadth = upr90_thermal - lwr10_thermal, depth_niche_breadth = upr90_depth - lwr10_depth, 
                       species = str_replace_all(species, " ", "_"), # clean species names 
                       sp_region = paste(species, region, sep = "_") 
)

data <- data |> group_by(region) |> mutate(
  year_c = (year - mean(year, na.rm = TRUE))/10 # center year
)  |> ungroup() |> 
  group_by(sp_region) |> # std index of abundance, and calculate species-specific anomalies
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


# compute abundance trends ------------------------------------------------

abundance_trends <- data %>%
  group_by(species,region) %>%
  do(model = lm(index_std ~ 1 + year_c, data = .)) %>%
  mutate(slope = coef(model)[["year_c"]])



# merge datasets ----------------------------------------------------------

win_lose <- slopes %>% select(species, region, outcome, mean_slope) %>% left_join(abundance_trends %>% mutate(index = slope) %>% select(species,region,index)) %>% left_join(pc_scores %>% select(Dim.1,species,region))


# winner and losers -----------------------------------------------------

# top row plot 

# 
# outcomes_pretty <- c(
#   "cogyc" = "Latitudinal shift",
#   "cogxc" = "Longitudinal shift",
#   "depthnichec" = "Depth shift",
#   "thermalnichec" = "Thermal niche shift"
# )

#  "cogyc"        = expression("km decade"^-1),
#  "cogxc"        = expression("km decade"^-1),
#  "depthnichec"  = expression("m decade"^-1),
#  "thermalnichec"= expression(degree*C ~ decade^-1)
#)

outcomes_labels <- list(
  "cogyc"        = expression(atop("Latitudinal shift", "(km decade"^-1)),
  "cogxc"        = expression(atop("Longitudinal shift", "(km decade"^-1)),
  "depthnichec"  = expression(atop("Depth shift", "(m decade"^-1)),
  "thermalnichec"= expression(atop("Thermal niche shift", degree*C~decade^-1))
)

win_lose <- win_lose %>%
  mutate(pretty_outcome = outcomes_pretty[outcome])

# regression stats for annotation

lm_stats <- win_lose %>%
  group_by(outcome) %>%
  summarise(model = list(lm(index ~ 1 + mean_slope)), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    slope = coef(model)[[2]],
    r2 = summary(model)$r.squared,
    p_value = summary(model)$coefficients[2, 4]
  ) %>%
  ungroup() %>%
  mutate(
    pretty_outcome = outcomes_pretty[outcome],
    label = paste0("Slope = ", round(slope, 3),
                   "\n R² = ", round(r2, 2),
                   "\n p = ", signif(p_value, 2))
  )

# color and labels

predictor_colors <- c(
  "cogyc" = "#762A83",
  "cogxc" = "#543005",
  "depthnichec" = "#8C510A",
  "thermalnichec" = "#2166AC",
  "Dim.1" = "#1B9E77"   # new one for Dim.1
)

#predictor_labels <- c(
#  #"Dim.1" = "Community position (Dim.1)",
#  "cogyc" = "Latitudinal shift (km/yr)",
#  "cogxc" = "Longitudinal shift (km/yr)",
#  "depthnichec" = "Depth shift (m/yr)",
#  "thermalnichec" = "Thermal niche shift (°C/yr)"
#)


make_plot <- function(var_name) {
  df <- win_lose %>% filter(outcome == var_name)
  stats <- lm_stats %>% filter(outcome == var_name)
  
  ggplot(df, aes(x = mean_slope, y = index)) +
    geom_point(alpha = 0.3, color = predictor_colors[var_name], size = 2) +
    geom_smooth(method = "lm", se = TRUE,
                color = predictor_colors[var_name],
                fill = predictor_colors[var_name],
                alpha = 0.15) +
    geom_text(
      data = stats,
      aes(x = Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = 1.1, vjust = 1.1,
      size = 3, fontface = "bold"
    ) +
    labs(
      x = outcomes_labels[[var_name]],
      y = expression(atop("Change in standardized abundance index",
                          "(decade"^-1*")"))
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

# make scatterplots for all predictors
p1 <- make_plot("cogyc")
p2 <- make_plot("cogxc")
p3 <- make_plot("depthnichec")
p4 <- make_plot("thermalnichec")
p5 <- ggplot() + theme_void()

p_top <- p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 5, axes = "collect")  & theme(axis.title.y = element_text(size = 9), axis.title.x = element_text(size = 9))


# coef plot by region -----------------------------------------------------

# add dim 1 

win_lose <- win_lose %>%
# create a new row for Dim.1
bind_rows(
  win_lose %>%
    select(species, region, index, Dim.1) %>%
    mutate(
      outcome = "Dim.1",
      mean_slope = Dim.1,
      pretty_outcome = "Dim.1"
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
  )


# bottom row plot slopes --------------------------------------------------

# nice math-friendly x-axis labels for slope plot
predictor_axis_labels <- list(
  cogyc = expression(atop(Delta * " Abundance decade"^-1,
                          "with a km decade"^-1 * " latitudinal shift")),
  cogxc = expression(atop(Delta * " Abundance decade"^-1,
                          "with a km decade"^-1 * " longitudinal shift")),
  depthnichec = expression(atop(Delta * " Abundance decade"^-1,
                                "with a m decade"^-1 * " depth shift")),
  thermalnichec = expression(atop(Delta * " Abundance decade"^-1,
                                  "with a " * degree * "C decade"^-1 * " thermal niche shift")),
  Dim.1 = expression(atop(Delta * " Abundance decade"^-1,
                          "with a unit PC1 change"))
)
  
# plot coefs
make_coef_plot <- function(outcome_name) {
  df <- coef_df %>% filter(outcome ==   outcome_name)

  ggplot(df, aes(x = estimate, y = region)) +
    ggstats::geom_stripped_rows(
      data = df,
      aes(y = region),
      odd = "white", even = "grey90", alpha = .2, inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_pointrange(
      aes(xmin = conf.low,
          xmax = conf.high,
          alpha = sig),
      color = predictor_colors[[outcome_name]],  # outcome-based colors
      size = 0.8
    ) +
    scale_alpha_manual(values = c("sig" = 1, "ns" = 0.4)) +
    labs(x = predictor_axis_labels[[outcome_name]],y=NULL) +
    theme_bw(base_size = 14)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 11),
          legend.position = "none")
}

# make coefficient plots for all predictors
c1 <- make_coef_plot("cogyc")
c2 <- make_coef_plot("cogxc")
c3 <- make_coef_plot("depthnichec")
c4 <- make_coef_plot("thermalnichec")
c5 <- make_coef_plot("Dim.1")

p_bottom <- c1 + c2 + c3 + c4 + c5 + plot_layout(ncol = 5, axes = 'collect')

p_top / p_bottom + plot_layout(heights = c(1, 1))

cowplot::plot_grid(p_top, p_bottom, ncol = 1, align = "v", rel_heights = c(1,1 ), axis = "ll", labels = c("a)", "b)"),
                   label_size = 12)

ggsave(here('multivariate_responses/figures/win_lose.png'),   width = 14,    # inches, matches \textwidth in Overleaf
       height = 10,   # adjust for aspect ratio
       dpi = 600, units = 'in')

# coef_df %>% filter(sig == 'sig')
# # ────────────────────────────────
# # 📦 6. Combine with patchwork
# # ────────────────────────────────
# 
# # Add an empty plot as a placeholder for alignment
# empty_plot <- ggplot() + theme_void()
# 
# # Pad the top row (4 panels) with an empty plot so it has 5 columns
# top_row <- plot_grid(
#   p_global , empty_plot,  # 5 panels total
#   ncol = 5,
#   align = "v",
#   axis = "tb"
# )
# 
# # Bottom row (already 5 panels)
# bottom_row <- plot_grid(
#   c1, c2, c3, c4, c5,
#   ncol = 5,
#   align = "v",
#   axis = "tb"
# )
# 
# # Combine top and bottom rows vertically
# final_plot <- plot_grid(
#   top_row,
#   bottom_row,
#   ncol = 1,
#   align = "v",
#   axis = "l"
# )
# 
# # Adjust common theme elements if needed
# final_plot
# 
# cowplot::plot_grid(p_global, p_coef, ncol = 1, align = "v", rel_heights = c(1,1 ), axis = "ll", labels = c("a)", "b)"),
#                    label_size = 12)
# 
# 
# ggsave(here('plots/main/winner_loser.png'),   width = 14,    # inches, matches \textwidth in Overleaf
#        height = 10,   # adjust for aspect ratio
#        dpi = 600, units = 'in')
# 
# 
# write_csv(coef_df,here('meta_regression/analysis/coef.csv'))
# 
# 
# plot_data |> filter(region == 'NIC', variable == 'thermalnichec') |> ggplot(aes(value, abu_trend))+geom_point()+geom_smooth(method = 'lm')
