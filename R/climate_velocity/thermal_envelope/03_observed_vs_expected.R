# libraries  --------------------------------------------------------------
library(tidyverse)
library(here)
library(knitr)

#  read in data -----------------------------------------------------------

observed <-  read_rds(here('data/processed/derived_quantities.rds')) %>% select(year,cog_y,cog_x, depth_niche, species,region)
expected <- read_rds(here('data/processed/derived_quantities_envelope.rds'))%>% select(year,cog_y,cog_x, depth_niche, species,region)

data <- data |> mutate(thermal_niche_breadth = upr90_thermal - lwr10_thermal, depth_niche_breadth = upr90_depth - lwr10_depth, 
                       species = str_replace_all(species, " ", "_"), # clean species names 
                       sp_region = paste(species, region, sep = "_") 
)


# standardize -------------------------------------------------------------

# helper function for standardization and anomalies
process_data <- function(df) {
  df %>%
    mutate(
      species = str_replace_all(species, " ", "_"),
      sp_region = paste(species, region, sep = "_")
    ) %>%
    group_by(region) %>%
    mutate(year_c = (year - mean(year, na.rm = TRUE)) / 10) %>%
    ungroup() %>%
    group_by(sp_region) %>%
    mutate(across(
      c(cog_y, cog_x, depth_niche),
      ~ .x - mean(.x, na.rm = TRUE),
      .names = "{.col}_c"
    )) %>%
    ungroup()
}

# apply to both observed and expected
observed_std <- process_data(observed)
expected_std <- process_data(expected)


# get slopes from lms -----------------------------------------------------

# function to fit slopes for given dataset
get_slopes <- function(df, type_label) {
  df %>%
    group_by(species, region, sp_region) %>%
    # pivot variables into long form to loop across them
    pivot_longer(cols = c(cog_y_c, cog_x_c, depth_niche_c),
                 names_to = "variable", values_to = "value") %>%
    group_by(species, region, sp_region, variable) %>%
    do(tidy(lm(value ~ year_c, data = .))) %>%
    filter(term == "year_c") %>%   # keep only slope
    select(species, region, sp_region, variable, slope = estimate) %>%
    mutate(type = type_label) %>%
    ungroup()
}

# apply to observed and expected
observed_slopes <- get_slopes(observed_std, "observed")
expected_slopes <- get_slopes(expected_std, "expected")

# combine
slopes_all <- bind_rows(observed_slopes, expected_slopes)

# reshape to wide format: observed + expected side by side
slopes_wide <- slopes_all %>%
  pivot_wider(
    names_from = type,
    values_from = slope
  )  %>% mutate(
    sign_obs = sign(observed),
    sign_exp = sign(expected),
    sign_agree = ifelse(sign_obs == sign_exp, 1, 0)
  )

# overall proportion of agreement per variable
sign_agree_overall <- slopes_wide %>%
  group_by(variable) %>%
  summarise(
    prop_agree = round(mean(sign_agree, na.rm = TRUE),2),
    n = sum(!is.na(sign_agree)),
    .groups = "drop"
  )

cor_results <- slopes_wide %>%
  group_by(variable) %>%
  summarise(
    test = list(cor.test(observed, expected, use = "complete.obs")),
    n = sum(!is.na(observed) & !is.na(expected)),
    .groups = "drop"
  ) %>%
  mutate(test = map(test, tidy)) %>%   # extract estimate & p-value
  unnest(test) %>%
  select(
    variable,
    cor_obs_exp = estimate,
    p.value,
    n
  ) %>%
  mutate(cor_obs_exp = round(cor_obs_exp, 2),
         p.value = signif(p.value, 2))


overall_tracking <- sign_agree_overall %>% left_join(cor_results) 

# proportion of agreement per variable and region
sign_agree_region <- slopes_wide %>%
  group_by(region, variable) %>%
  summarise(
    prop_agree = round(mean(sign_agree, na.rm = TRUE), 2),
    n = sum(!is.na(sign_agree)),
    .groups = "drop"
  )

# by region
cor_results_region <- slopes_wide %>%
  group_by(variable,region) %>%
  summarise(
    test = list(cor.test(observed, expected, use = "complete.obs")),
    n = sum(!is.na(observed) & !is.na(expected)),
    .groups = "drop"
  ) %>%
  mutate(test = map(test, tidy)) %>%   # extract estimate & p-value
  unnest(test) %>%
  select(
    region,
    variable,
    cor_obs_exp = estimate,
    p.value,
    n
  ) %>%
  mutate(cor_obs_exp = round(cor_obs_exp, 2),
         p.value = signif(p.value, 2))


regional_tracking <- sign_agree_region %>% left_join(cor_results_region) 

knitr::kable(regional_tracking, digits = 2, caption = "Observed vs Expected Slopes by Region and Variable")

ggplot(slopes_wide, aes(x = observed, y = expected)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  facet_wrap(~ variable, scales = "free") +
  labs(
    x = "Observed Slopes",
    y = "Expected Slopes",
    title = "Observed vs Expected Slopes by Variable"
  ) +
  theme_minimal()


