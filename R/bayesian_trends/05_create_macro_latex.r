library(here)
library(tidyverse)
library(readr)

# Read in the trend data
data <- readRDS(here('R/data/processed/bayesian_species_trends.rds'))


# helper functions form Sean https://github.com/pbs-assess/dogfish-trends/blob/main/analysis/080-values.R ---------------------------------------------

# rounding helper
mround <- function(x, digits) sprintf(paste0("%.", digits, "f"), round(x, digits))

# LaTeX macro writer
write_tex <- function(x, macro, append = TRUE) {
  paste0("\\newcommand{\\", macro, "}{", x, "}") |>
    readr::write_lines("output/values/species_values.tex", append = append)
}

# prepare species names in LaTeX italics
# get one-row tables for max and min per outcome
max_df <- data %>%
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

min_df <- data %>%
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

unlink("output/values/species_values.tex")


# write macros (one set per outcome)
for (i in seq_len(nrow(extremes))) {
  oc <- extremes$outcome[i]
  
  # max macros
  write_tex(extremes$max_species_tex[i], paste0("max", oc, "Species"))
  write_tex(extremes$max_region[i], paste0("max", oc, "Region"))
  write_tex(
    paste0(
      mround(extremes$max_median[i], 2),
      " (95\\% CI: ",
      mround(extremes$max_lower_95[i], 2), "--", mround(extremes$max_upper_95[i], 2), ")"
    ),
    paste0("max", oc, "Value")
  )
  
  # min macros
  write_tex(extremes$min_species_tex[i], paste0("min", oc, "Species"))
  write_tex(extremes$min_region[i], paste0("min", oc, "Region"))
  write_tex(
    paste0(
      mround(extremes$min_median[i], 2),
      " (95\\% CI: ",
      mround(extremes$min_lower_95[i], 2), "--", mround(extremes$min_upper_95[i], 2), ")"
    ),
    paste0("min", oc, "Value")
  )
}
