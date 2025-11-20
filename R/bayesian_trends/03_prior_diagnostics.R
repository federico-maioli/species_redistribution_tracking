library(tidyverse)
library(brms)
library(here)
library(broom)
library(ggdist)
library(patchwork)
#library(cowplot)
library(knitr)

# print priors  ----------------------------------------------------------

fit <- read_rds(here('R/bayesian_trends/fitted/m_stud.rds'))

clean_priors <- prior_summary(fit) %>% as.data.frame() %>%
  filter(prior != "") %>%          # drop rows with no prior
  select(prior, class, coef, group, resp, source) %>%  # keep key columns
  arrange(class, group, coef) %>%    # keep only rows with an explicit prior
  mutate(
    # convert distribution names to nicer LaTeX-friendly text
    prior = str_replace_all(prior, fixed("student_t"), "Student-t"),
    prior = str_replace_all(prior, fixed("normal"), "Normal"),
    prior = str_replace_all(prior, fixed("gamma"), "Gamma"),
    prior = str_replace_all(prior, fixed("lkj_corr_cholesky"), "LKJ"),
    
    # replace with Greek symbols
    class = case_when(
      class == "b" ~ "$\\beta$",
      class == "sd" ~ "$\\sigma$",
      class == "sigma" ~ "$\\sigma$",
      class == "nu" ~ "$\\nu$",
      class == "L" ~ "$\\L$",
      TRUE ~ class
    ),
    # replace coef variable names
    coef = case_when(coef == "year_c" ~ "decade", TRUE ~ coef),
    # replace resp
    resp = case_when(
      resp == "cogxc" ~ "lat centroid",
      resp == "cogyc" ~ "lon centroid",
      resp == "depthnichec" ~ "depth niche",
      resp == "thermalnichec" ~ "thermal niche",
      TRUE ~ resp
    )
  ) %>%
  rename(
    Prior = prior,
    Parameter = class,
    Coefficient = coef,
    Group  = group,
    Response = resp,
    Source = source
  )

# Export LaTeX table
table_priors <- kable(
  clean_priors, 
  format = "latex", 
  label = "priors",
  booktabs = TRUE,
  escape = FALSE,
  # caption = "Prior distributions for model parameters. Fixed effects ($\\beta$), group-level and residual standard deviations ($\\sigma$), degrees of freedom ($\\nu$) for the Student-t distribution, and the correlation matrix of group-level effects parameterized via its Cholesky factor ($L$). The source column indicates whether the prior was changed from the default.",
  linesep = "" 
)

writeLines(table_priors, here('output/tables/supp/table_priors_supp.tex'))

# prior validation -----------------------------------------------------------

## global slopes ------------------------------------------------------------

# define priors
priors <- c(
  prior(normal(0, 50), class = b, coef = year_c, resp = cogxc),
  prior(normal(0, 50), class = b, coef = year_c, resp = cogyc),
  prior(normal(0, 5), class = b, coef = year_c, resp = depthnichec),
  prior(normal(0, 0.5), class = b, coef = year_c, resp = thermalnichec)
)

# metadata for plotting
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
    "beta^{lat~centroid}",
    "beta^{lon~centroid}",
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
  stat_halfeye(normalize = "xy", color = "steelblue4", fill = "steelblue1",alpha=.7) +
  geom_vline(aes(xintercept = vline), colour = "black", linetype = "dashed", size = .6) +
  geom_text(aes(x = vline, y = 0.25, label = paper),
            color = "black", angle = 90, vjust = -0.6, hjust = 0, size = 3) +
  facet_wrap(~prior, scales = "free") + scale_y_discrete(labels = function(x) parse(text = x)) +
  labs(
    x = NULL,
    y = NULL
  ) + theme_bw() + theme(
    strip.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  ) 

p1

#  for region -------------------------------------------------------------

# get data

data <- readRDS(here("R/data/processed/derived_quantities_sdm.rds"))


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
    index_std = scale(index)[, 1]) |> 
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
    "beta['region']^{lat~centroid}",
    "beta['region']^{lon~centroid}",
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
      prior == "student_t(3, 0, 30)"   ~ "Student-t(3, 0, 3)",
      prior == "student_t(3, 0, 10)"   ~ "Student-t(3, 0, 1)",
      prior == "student_t(3, 0, 0.2)" ~ "Student-t(3, 0, 0.2)",
      TRUE ~ prior
    ))

facet_labels <- priors_df %>%
  distinct(facet_id, facet_label) %>%
  deframe

# plot 
p2 = ggplot(priors_df, aes(
  y = resp_label,
  xdist = .dist_obj
)) +
  stat_halfeye(normalize = "xy", color = "steelblue4", fill = "steelblue1", p_limits = c(0, 0.99),alpha = .7) +
  geom_vline(aes(xintercept = vline), colour = "orange", linetype = "dashed",size = .6) +
  facet_wrap(
    ~facet_id,
    labeller = labeller(facet_id = facet_labels),
    scales = "free"
  ) + scale_y_discrete(labels = function(x) parse(text = x)) + scale_x_continuous(expand = c(0,0)) +
  labs(
    x = NULL,
    y = NULL
  )  + theme_bw() + theme(
    strip.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  ) 

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
    "beta['region:species']^{lat~centroid}",
    "beta['region:species']^{lon~centroid}",
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
      prior == "student_t(3, 0, 40)"   ~ "Student-t(3, 0, 40)",
      prior == "student_t(3, 0, 20)"   ~ "Student-t(3, 0, 20)",
      prior == "student_t(3, 0, 0.4)" ~ "Student-t(3, 0, 0.4)",
      TRUE ~ prior
    ))

# plot 
p3 = ggplot(priors_df, aes(
  y = resp_label,
  xdist = .dist_obj
)) +
  stat_halfeye(normalize = "xy", color = "steelblue4", fill = "steelblue1", p_limits = c(0, 0.99),alpha =.7) +
  geom_vline(aes(xintercept = vline), colour = "orange", linetype = "dashed",size = .6) +
  facet_wrap(
    ~facet_id,
    scales = "free",
    labeller = labeller(facet_id = setNames(priors_df$facet_label, priors_df$facet_id))
  ) + scale_y_discrete(labels = function(x) parse(text = x)) + scale_x_continuous(expand = c(0,0)) +
  labs(
    x = NULL,
    y = NULL
  ) + theme_bw() + theme(
    strip.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  ) 

combined <- (p1 / p2 /p3) & plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 12, face='bold'))

#ggdraw(combined) +
#  draw_line(x = c(0, 1), y = c(0.67, 0.67), color = "gray90", size = 0.5) +
#  draw_line(x = c(0, 1), y = c(0.34, 0.34), color = "gray90", size = 0.5)

combined

ggsave(
  here('output/figures/supp/priors_supp.png'),
  #plot = pp_plot,  
  width = 180,  
  height = 170,
  dpi = 600,
  units = "mm"
)

