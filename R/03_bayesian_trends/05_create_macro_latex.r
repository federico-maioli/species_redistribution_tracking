library(here)
library(tidyverse)
library(readr)

# Read in the trend data
global <- readRDS(here('R/data/processed/bayesian_global_trends.rds')) %>%  mutate(
  outcome = case_when(
    outcome == "cogyc" ~ "lat",
    outcome == "cogxc" ~ "lon",
    outcome == "depthnichec" ~ "depth",
    outcome == "thermalnichec" ~ "thermal",
    TRUE ~ as.character(outcome)
  )
)
region <- readRDS(here('R/data/processed/bayesian_region_trends.rds')) %>%  mutate(
  outcome = case_when(
    outcome == "cogyc" ~ "lat",
    outcome == "cogxc" ~ "lon",
    outcome == "depthnichec" ~ "depth",
    outcome == "thermalnichec" ~ "thermal",
    TRUE ~ as.character(outcome)
  ), region = case_when(
    region == "NEUS-SS" ~ "NEUS",
    TRUE ~ as.character(region)
  )
)
species <- readRDS(here('R/data/processed/bayesian_species_trends.rds'))


# helper functions form Sean https://github.com/pbs-assess/dogfish-trends/blob/main/analysis/080-values.R ---------------------------------------------

# rounding helper
mround <- function(x, digits) sprintf(paste0("%.", digits, "f"), round(x, digits))


# global slope ------------------------------------------------------------

# LaTeX macro writer
write_tex <- function(x, macro, append = TRUE) {
  paste0("\\newcommand{\\", macro, "}{", x, "}") |>
    readr::write_lines("output/values/bayesian_trend_analysis/global_slopes.tex", append = append)
}

unlink("output/values/bayesian_trend_analysis/global_slopes.tex")

for (i in seq_len(nrow(global))) {
  oc <- tools::toTitleCase(tolower(global$outcome[i])) 
  
  # write the median slope macro
  write_tex(
    mround(global$median_slope[i], 2),
    paste0(oc, "Median")
  )
  
  # write the CI-only macro in parentheses
  write_tex(
    paste0(
      "[95\\% CI: ",
      mround(global$lower_95[i], 2), ", ",
      mround(global$upper_95[i], 2), "]"
    ),
    paste0(oc, "CI")
  )
}


# region slope ------------------------------------------------------------

# LaTeX macro writer
write_tex <- function(x, macro, append = TRUE) {
  paste0("\\newcommand{\\", macro, "}{", x, "}") |>
    readr::write_lines("output/values/bayesian_trend_analysis/region_slopes.tex", append = append)
}

# remove existing file
unlink("output/values/bayesian_trend_analysis/region_slopes.tex")

# loop through rows
for (i in seq_len(nrow(region))) {
  oc <- tools::toTitleCase(tolower(region$outcome[i])) 
  re <- region$region[i]
  
  # combine without underscore
  macro_base <- paste0(oc, re)
  
  # median macro 
  write_tex(
    mround(region$median_slope[i], 2),
    paste0(macro_base, "Median")
  )
  
  # CI-only macro in parentheses 
  write_tex(
    paste0(
      "[95\\% CI: ",
      mround(region$lower_95[i], 2), ", ", 
      mround(region$upper_95[i], 2), "]"
    ),
    paste0(macro_base, "CI")
  )
}

# species slopes, just min and max ----------------------------------------------------------

# prepare species names in LaTeX italics
# get one-row tables for max and min per outcome
max_df <- species %>%
  group_by(outcome) %>%
  slice_max(order_by = median_slope, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    outcome,
    max_region = region,
    max_species = species,
    max_species_tex = paste0("\\textit{", gsub("_", " ", species), "}"),
    max_median = median_slope,
    max_lower_95 = lower_95,
    max_upper_95 = upper_95
  )

min_df <- species %>%
  group_by(outcome) %>%
  slice_min(order_by = median_slope, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    outcome,
    min_region = region,
    min_species = species,
    min_species_tex = paste0("\\textit{", gsub("_", " ", species), "}"),
    min_median = median_slope,
    min_lower_95 = lower_95,
    min_upper_95 = upper_95
  )

# join by outcome (now tidy: one row per outcome, clear column names)
extremes <- left_join(max_df, min_df, by = "outcome")

extremes <- extremes %>% mutate(
  outcome = case_when(
    outcome == "lat_shift" ~ "lat",
    outcome == "lon_shift" ~ "lon",
    outcome == "depth_shift" ~ "depth",
    outcome == "thermal_shift" ~ "thermal",
    TRUE ~ as.character(outcome)
  )
)

# remove existing file
unlink("output/values/bayesian_trend_analysis/species_slopes.tex")
# LaTeX macro writer
write_tex <- function(x, macro, append = TRUE) {
  paste0("\\newcommand{\\", macro, "}{", x, "}") |>
    readr::write_lines("output/values/bayesian_trend_analysis/species_slopes.tex", append = append)
}

# loop through each outcome
for (i in seq_len(nrow(extremes))) {
  oc <- tools::toTitleCase(tolower(extremes$outcome[i]))  # Capitalized outcome
  
  ## ----- MAX macros -----
  
  # species and region
  write_tex(extremes$max_species_tex[i], paste0("Max", oc, "Species"))
  write_tex(extremes$max_region[i], paste0("Max", oc, "Region"))
  
  # median
  write_tex(
    mround(extremes$max_median[i], 0),
    paste0("Max", oc, "Median")
  )
  
  # CI only
  write_tex(
    paste0(
      "[95\\% CI: ",
      mround(extremes$max_lower_95[i], 0), ", ",
      mround(extremes$max_upper_95[i], 0), "]"
    ),
    paste0("Max", oc, "CI")
  )
  
  ## ----- MIN macros -----
  
  # species and region
  write_tex(extremes$min_species_tex[i], paste0("Min", oc, "Species"))
  write_tex(extremes$min_region[i], paste0("Min", oc, "Region"))
  
  # median
  write_tex(
    mround(extremes$min_median[i], 0),
    paste0("Min", oc, "Median")
  )
  
  # CI only
  write_tex(
    paste0(
      "[95\\% CI: ",
      mround(extremes$min_lower_95[i], 0), ", ",
      mround(extremes$min_upper_95[i], 0), "]"
    ),
    paste0("Min", oc, "CI")
  )
}


# proportion ----------------------------------------------------

prop_by_outcome <- species %>%
  group_by(outcome) %>%
  summarise(
    n = n(),
    n_sig = sum(significant == "yes", na.rm = TRUE),
    prop_sig = (n_sig / n) * 100
  )  %>% mutate(outcome = case_when(
    outcome == "lat_shift" ~ "lat",
    outcome == "lon_shift" ~ "lon",
    outcome == "depth_shift" ~ "depth",
    outcome == "thermal_shift" ~ "thermal",
    TRUE ~ as.character(outcome)
  ))

# Define LaTeX macro writer
write_tex <- function(x, macro, append = TRUE) {
  paste0("\\newcommand{\\", macro, "}{", x, "}") |>
    write_lines("output/values/bayesian_trend_analysis/prop_signif.tex", append = append)
}

# Remove old output file
unlink("output/values/bayesian_trend_analysis/prop_signif.tex")


for (i in seq_len(nrow(prop_by_outcome))) {
  oc <- tools::toTitleCase(tolower(prop_by_outcome$outcome[i]))
  write_tex(
    paste0(mround(prop_by_outcome$prop_sig[i],0), "\\%"),
    paste0(oc, "Perc")
  )
}


# proportion by direction -------------------------------------------------

prop_by_outcome_dir <- species %>%
  filter(direction_label != "Not significant") %>%  # keep only meaningful directions
  group_by(outcome, direction_label) %>%
  summarise(
    n_sig = sum(significant == "yes", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # get total rows per outcome
  left_join(
    species %>% group_by(outcome) %>% summarise(n_total = n(), .groups = "drop"),
    by = "outcome"
  ) %>%
  mutate(
    prop_sig = 100 * n_sig / n_total,
    outcome_short = case_when(
      outcome == "lat_shift" ~ "lat",
      outcome == "lon_shift" ~ "lon",
      outcome == "depth_shift" ~ "depth",
      outcome == "thermal_shift" ~ "thermal",
      TRUE ~ as.character(outcome)
    ),
    direction_clean = gsub(" ", "", direction_label)
  )
# View
prop_by_outcome_dir

out_file <- "output/values/bayesian_trend_analysis/prop_signif_by_dir.tex"
unlink(out_file)

write_tex <- function(x, macro, append = TRUE) {
  paste0("\\newcommand{\\", macro, "}{", x, "}") |>
    write_lines(out_file, append = append)
}

for (i in seq_len(nrow(prop_by_outcome_dir))) {
  row <- prop_by_outcome_dir[i, ]
  
  # e.g., \LatNorthingPerc
  macro_name <- paste0(
    tools::toTitleCase(row$outcome_short),
    tools::toTitleCase(row$direction_clean),
    "Perc"
  )
  
  write_tex(paste0(mround(row$prop_sig,0), "\\%"), macro_name)
}

# thermal niche warming in pace with local warming ------------------------

# get regional slopes

grid <- read_rds(here('R/data/processed/prediction_grid.rds'))

# summarize mean temperature per year-region
temp_summary <- grid %>%
  group_by(year, region_short) %>%
  summarise(mean_temp = mean(mean_year_temp, na.rm = TRUE), .groups = "drop") 

# compute slopes per region
slopes <- temp_summary %>%
  group_by(region_short) %>%
  summarise(
    slope_decade = coef(lm(mean_temp ~ year))[2] * 10,
    .groups = "drop"
  ) 


region_termal_niche_slope <- region %>% filter(outcome == 'thermal') %>% select(region,median_slope)

comb <- slopes %>% left_join(region_termal_niche_slope, by = c('region_short'='region'))

# calculate corr
test_result <- cor.test(comb$slope_decade, comb$median_slope, use = "complete.obs", method = "pearson")

