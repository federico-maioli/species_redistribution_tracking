# libraries  --------------------------------------------------------------
library(tidyverse)
library(here)
library(knitr)
library(broom)

#  read in data -----------------------------------------------------------

observed <-  read_rds(here('R/data/processed/derived_quantities_sdm.rds')) %>% select(year,cog_y,cog_x, depth_niche, species,region)
expected <- read_rds(here('R/data/processed/derived_quantities_te.rds'))%>% select(year,cog_y,cog_x, depth_niche, species,region)


# standardize -------------------------------------------------------------

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
    do(broom::tidy(lm(value ~ year_c, data = .))) %>%
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
  ) %>% na.omit() # remove species if not in obs and exp

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
         p.value = round(signif(p.value, 2),3))


overall_tracking <- sign_agree_overall %>% left_join(cor_results) 

clean_tracking <- overall_tracking %>%
  mutate(
    # Replace variable names with readable "Shift type" labels
    variable = case_when(
      variable == "cog_y_c" ~ "Latitudinal",
      variable == "cog_x_c" ~ "Longitudinal",
      variable == "depth_niche_c" ~ "Depth",
      TRUE ~ variable
    ),
    # Set factor levels to control order
    variable = factor(variable, levels = c("Latitudinal", "Longitudinal", "Depth"))
  ) %>%
  rename(
    "Shift type" = variable,
    "Proportion of aligned responses" = prop_agree,
    "N" = n,
    "$\\rho$" = cor_obs_exp,
    "$p$-value" = p.value
  )

# Export LaTeX table
table_exp_obs <- kable(
  clean_tracking, 
  format = "latex", 
  label = "exp_obs",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Correlation ($\rho$) and directional agreement between observed species shifts and those expected under thermal envelope tracking. The proportion of aligned responses indicates the fraction of species whose observed and expected shifts occurred in the same direction (i.e., shared the same sign). N denotes the number of species compared for each shift type."
)

writeLines(table_exp_obs, here('output/tables/main/table_exp_obs.tex'))


# now by region


# proportion of agreement per variable and region
sign_agree_region <- slopes_wide %>%
  group_by(region, variable) %>%
  summarise(
    prop_agree = round(mean(sign_agree, na.rm = TRUE), 2),
    n = sum(!is.na(sign_agree)),
    .groups = "drop"
  )


# By region ---------------------------------------------------------------


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
         p.value = round(signif(p.value, 2),3))


regional_tracking <- sign_agree_region %>% left_join(cor_results_region) 

clean_regional_tracking <- regional_tracking %>%
  mutate(
    # Replace variable codes with readable shift types
    variable = case_when(
      variable == "cog_y_c" ~ "Latitudinal",
      variable == "cog_x_c" ~ "Longitudinal",
      variable == "depth_niche_c" ~ "Depth",
      TRUE ~ variable
    ),
    region = case_when(
      region == "GOA"      ~ "Gulf of Alaska (GOA)",
      region == "BC"       ~ "British Columbia (BC)",
      region == "USWC"     ~ "U.S. West Coast (USWC)",
      region == "BAL"      ~ "Baltic Sea (BAL)",
      region == "NS"       ~ "North Sea (NS)",
      region == "BS"       ~ "Barents Sea (BS)",
      region == "EBS"      ~ "East Bering Sea (EBS)",
      region == "GOM"      ~ "Gulf of Mexico (GOM)",
      region == "CBS"      ~ "Celtic–Biscay Shelf (CBS)",
      region == "NIC"      ~ "North Iberian Coast (NIC)",
      region == "NEUS-SS"  ~ "NE U.S. & Scotian Shelf (NEUS–SS)",
      TRUE ~ region),
    # Set factor levels to ensure desired order
    variable = factor(variable, levels = c("Latitudinal", "Longitudinal", "Depth"))
  ) %>%
  rename(
    "Region" = region,
    "Shift type" = variable,
    "Proportion of aligned responses" = prop_agree,
    "N" = n,
    "$\\rho$" = cor_obs_exp,
    "$p$-value" = p.value
  ) %>%
  arrange(-`$\\rho$`)

# Export LaTeX table
table_exp_obs_region <- kable(
  clean_regional_tracking, 
  format = "latex", 
  label = "exp_obs_region",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Correlation ($\rho$) and directional agreement between observed species shifts and those expected under thermal envelope tracking. The proportion of aligned responses indicates the fraction of species whose observed and expected shifts occurred in the same direction (i.e., shared the same sign). N denotes the number of species compared for each shift type within each region."
)

writeLines(table_exp_obs_region, here('output/tables/supp/table_exp_obs_region.tex'))

# no write macro for all the values
mround <- function(x, digits) sprintf(paste0("%.", digits, "f"), round(x, digits))

# LaTeX macro writer
write_tex <- function(x, macro, append = TRUE) {
  paste0("\\newcommand{\\", macro, "}{", x, "}") |>
    readr::write_lines("output/values/exp_obs_region.tex", append = append)
}

unlink("output/values/exp_obs_region.tex")

# Helper: round correlation and format p-values nicely
format_rho_p <- function(rho, n, p) {
  rho_str <- round(rho, 2)
  
  # Format p-values according to Nature style
  p_str <- ifelse(p < 0.001, "< 0.001", paste0("= ", signif(p, 2)))
  
  paste0("($\\rho$ = ", rho_str, ", n = ", n, ", \\textit{p} ", p_str, ")")
}

# Loop over each row of regional_tracking (or overall_tracking)
for (i in seq_len(nrow(regional_tracking))) {
  variable <- regional_tracking$variable[i]
  region <- regional_tracking$region[i]
  rho <- regional_tracking$cor_obs_exp[i]
  n <- regional_tracking$n[i]
  p <- regional_tracking$p.value[i]
  
  # Convert variable names to a compact version for macro name
  var_macro <- gsub("_c", "", variable)  # remove trailing _c
  var_macro <- gsub("_", "", var_macro)  # remove underscores
  
  # Construct macro name
  macro_name <- paste0("rho_", var_macro, "_", region)
  
  # Construct LaTeX content
  macro_value <- format_rho_p(rho, n, p)
  
  # Write to .tex file
  write_tex(macro_value, macro_name)
}

