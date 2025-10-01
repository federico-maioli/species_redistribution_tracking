library(tidyverse)
library(brms)
library(here)
library(broom)
library(ggdist)
library(patchwork)
library(cowplot)

# read in model  ----------------------------------------------------------

fit <- read_rds(here('bayesian_trends/fitted/m_final.rds'))

get_prior(fit)

# global slopes -----------------------------------------------------------

# Define priors
priors <- c(
  prior(normal(0, 50), class = b, coef = year_c, resp = cogxc),
  prior(normal(0, 50), class = b, coef = year_c, resp = cogyc),
  prior(normal(0, 5), class = b, coef = year_c, resp = depthnichec),
  prior(normal(0, 0.5), class = b, coef = year_c, resp = thermalnichec)
)

# Metadata for plotting
priors_meta <- tibble::tibble(
  coef  = "year_c",
  resp  = c("cogxc", "cogyc", "depthnichec", "thermalnichec"),
  vline = c(30.6, 30.6, 3.6, 0.37),
  paper = c(
    "Poloczanska et al. (2013)",
    "Poloczanska et al. (2013)",
    "Dulvy et al. (2008)",
    "Chen et al. (2020)"
  ),
  resp_label = c(
    "beta^{latitude~centroid}",
    "beta^{longitude~centroid}",
    "beta^{depth~niche}",
    "beta^{thermal~niche}"
  )
)

priors_df <- priors %>%
  parse_dist(prior) %>%
  left_join(priors_meta, by = c("resp", "coef")) |> mutate( prior = case_when(
    prior == "normal(0, 50)"   ~ "Normal(0, 50)",
    prior == "normal(0, 5)"   ~ "Normal(0, 5)",
    prior == "normal(0, 0.5)"   ~ "Normal(0, 0.5)",
    TRUE ~ prior))

p1 = ggplot(priors_df, aes(
  y = resp_label,
  xdist = .dist_obj
)) +
  stat_halfeye(normalize = "xy", color = "steelblue4", fill = "steelblue1") +
  geom_vline(aes(xintercept = vline), colour = "black", linetype = "dashed") +
  geom_text(aes(x = vline, y = 0.2, label = paper),
            color = "black", angle = 90, vjust = -0.6, hjust = 0, size = 3) +
  facet_wrap(~prior, scales = "free") + scale_y_discrete(labels = function(x) parse(text = x)) +
  labs(
    x = NULL,
    y = NULL
  ) + theme_minimal(base_size = 12)



#  for region -------------------------------------------------------------

# get data

data <- readRDS(here("data/processed/derived_quantities.rds"))


data <- data |> mutate(
  species = str_replace_all(species, " ", "_"),
  # clean species names
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
    c(cog_y, cog_x, depth_niche, thermal_niche),
    ~ .x - mean(.x, na.rm = TRUE),
    .names = "{.col}_c"
  )) |>  # calculated anomalies
  ungroup()

sd_ols_region <- data %>%
  pivot_longer(
    cols = c(cog_y_c, cog_x_c, depth_niche_c, thermal_niche_c),
    names_to = "var", values_to = "value"
  ) %>%
  group_by(region, var) %>%
  do(tidy(lm(value ~ 0 + year_c, data = .))) %>%
  filter(term == "year_c") %>%
  group_by(var) %>%
  summarise(sd = sd(estimate, na.rm = TRUE), .groups = "drop")

sd_ols_region


priors <- c(
  prior(
    student_t(3, 0, 30),
    class = "sd",
    group = "region",
    coef = "year_c",
    resp = "cogyc", lb =0 
  ),
  prior(
    student_t(3, 0, 30),
    class = "sd",
    group = "region",
    coef = "year_c",
    resp = "cogxc", lb =0 
  ),
  prior(
    student_t(3, 0, 10),
    class = "sd",
    group = "region",
    coef = "year_c",
    resp = "depthnichec", lb =0 
  ),
  prior(
    student_t(3, 0, 0.2),
    class = "sd",
    group = "region",
    coef = "year_c",
    resp = "thermalnichec", lb = 0
  )
)


# Metadata for plotting
priors_meta <- tibble::tibble(
  coef  = "year_c",
  resp  = c("cogxc", "cogyc", "depthnichec", "thermalnichec"),
  vline = c(9.24, 13.9, 1.68, 0.178),
  resp_label = c(
    "beta['region']^{latitude~centroid}",
    "beta['region']^{longitude~centroid}",
    "beta['region']^{depth~niche}",
    "beta['region']^{thermal~niche}"
  )
)

priors_df <- priors %>%
  parse_dist(prior) %>%
  left_join(priors_meta, by = c("resp", "coef")) %>%   mutate(
    # Give each row a unique facet ID based on prior + resp
    facet_id = paste0(prior, "_", resp),
    
    # Assign the label for the facet strip
    # Facet label now is the prior distribution
    facet_label = case_when(
      prior == "student_t(3, 0, 3)"   ~ "Student t(3, 0, 3)",
      prior == "student_t(3, 0, 1)"   ~ "Student t(3, 0, 1)",
      prior == "student_t(3, 0, 0.02)" ~ "Student t(3, 0, 0.02)",
      TRUE ~ prior
    ))

# plot 
p2 = ggplot(priors_df, aes(
  y = resp_label,
  xdist = .dist_obj
)) +
  stat_halfeye(normalize = "xy", color = "steelblue4", fill = "steelblue1", p_limits = c(0, 0.99)) +
  geom_vline(aes(xintercept = vline), colour = "orange", linetype = "dashed",size = 1) +
  facet_wrap(
    ~facet_id,
    scales = "free",
    labeller = labeller(facet_id = setNames(priors_df$facet_label, priors_df$facet_id))
  ) + scale_y_discrete(labels = function(x) parse(text = x)) +
  labs(
    x = NULL,
    y = NULL
  ) + theme_minimal(base_size = 12)

# species in region -------------------------------------------------------

# maximum likelihood estimates of SD
sd_ols <- data %>%
  tidyr::pivot_longer(
    cols = c(cog_y_c, cog_x_c, depth_niche_c, thermal_niche_c),
    names_to = "var", values_to = "value"
  ) %>%
  group_by(sp_region, var) %>%
  do(tidy(lm(value ~ 0 + year_c, data = .))) %>%
  filter(term == "year_c") %>%
  group_by(var) %>%
  summarise(sd = sd(estimate, na.rm = TRUE), .groups = "drop")

# priors
# Species nested in region priors for slopes
priors <- c(
  prior(
    student_t(3, 0, 40),
    class = "sd",
    group = "region:species",
    coef = "year_c",
    resp = "cogyc", lb =0 
  ),
  prior(
    student_t(3, 0, 40),
    class = "sd",
    group = "region:species",
    coef = "year_c",
    resp = "cogxc", lb =0 
  ),
  prior(
    student_t(3, 0, 20),
    class = "sd",
    group = "region:species",
    coef = "year_c",
    resp = "depthnichec", lb =0 
  ),
  prior(
    student_t(3, 0, 0.4),
    class = "sd",
    group = "region:species",
    coef = "year_c",
    resp = "thermalnichec", lb = 0
  )
)

# Metadata for plotting
priors_meta <- tibble::tibble(
  coef  = "year_c",
  resp  = c("cogxc", "cogyc", "depthnichec", "thermalnichec"),
  vline = c(31.7, 42.4, 8.06, 0.253),
  resp_label = c(
    "beta['region:species']^{latitude~centroid}",
    "beta['region:species']^{longitude~centroid}",
    "beta['region:species']^{depth~niche}",
    "beta['region:species']^{thermal~niche}"
  )
)

priors_df <- priors %>%
  parse_dist(prior) %>%
  left_join(priors_meta, by = c("resp", "coef")) %>%   mutate(
    # Give each row a unique facet ID based on prior + resp
    facet_id = paste0(prior, "_", resp),
    
    # Assign the label for the facet strip
    # Facet label now is the prior distribution
    facet_label = case_when(
      prior == "student_t(3, 0, 40)"   ~ "Student t(3, 0, 40)",
      prior == "student_t(3, 0, 20)"   ~ "Student t(3, 0, 20)",
      prior == "student_t(3, 0, 0.4)" ~ "Student t(3, 0, 0.4)",
      TRUE ~ prior
    ))

# plot 
p3 = ggplot(priors_df, aes(
  y = resp_label,
  xdist = .dist_obj
)) +
  stat_halfeye(normalize = "xy", color = "steelblue4", fill = "steelblue1", p_limits = c(0, 0.99)) +
  geom_vline(aes(xintercept = vline), colour = "orange", linetype = "dashed",size = 1) +
  facet_wrap(
    ~facet_id,
    scales = "free",
    labeller = labeller(facet_id = setNames(priors_df$facet_label, priors_df$facet_id))
  ) + scale_y_discrete(labels = function(x) parse(text = x)) +
  labs(
    x = NULL,
    y = NULL
  ) + theme_minimal(base_size = 12)


#cowplot::plot_grid(p1, p2, p3, ncol = 1, align = "v", axis = "lr", labels = c("a)", "b)", "c)"),label_size = 12)

combined <- plot_grid(
  p1, p2, p3,
  ncol = 1, align = "v", axis = "lr",
  labels = c("a)", "b)", "c)"), label_size = 10
)

# Overlay horizontal divider lines
ggdraw(combined) +
  draw_line(x = c(0, 1), y = c(0.67, 0.67), color = "gray90", size = 0.5) +
  draw_line(x = c(0, 1), y = c(0.34, 0.34), color = "gray90", size = 0.5)

ggsave(
  here("bayesian_trends/figures/supp/prior_checks.png"),
  width = 8, height = 9, dpi = 600, bg = "white"
)
