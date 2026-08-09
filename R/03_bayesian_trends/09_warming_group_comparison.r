# libraries ---------------------------------------------------------------

library(tidyverse)
library(here)
library(brms)
library(tidybayes)
library(posterior)
library(ggtext)
library(patchwork)
library(ggrepel)


# load model --------------------------------------------------------------

fit <- read_rds(here('R/03_bayesian_trends/fitted/m_stud.rds'))

# confirm parameter names before extraction
#model_vars <- variables(fit)
#print(grep("^b_.*year_c$", model_vars, value = TRUE))
#print(head(grep("^r_region__.*\\[.*year_c\\]$", model_vars, value = TRUE)))


# warming vs non-warming region grouping ----------------------------------

# regions with a significant bottom-temperature trend (p < 0.05)
warming_regions <- c("BC", "USWC", "NEUS-SS", "NS", "BAL")

# regions with no significant temperature trend
non_warming_regions <- c("EBS", "GOA", "GOM", "BS", "CBS", "NIC")


# extract posterior draws -------------------------------------------------

# note: year_c was scaled to decades before fitting (see 03_prior_diagnostics.R
# and 04_plot_model_slopes.r), so the slopes below are already in per-decade
# units and require no further conversion.
region_draws <- fit %>%
  spread_draws(
    b_cogyc_year_c, b_cogxc_year_c,
    b_depthnichec_year_c, b_thermalnichec_year_c,
    r_region__cogyc[region, year_c],
    r_region__cogxc[region, year_c],
    r_region__depthnichec[region, year_c],
    r_region__thermalnichec[region, year_c]
  )


# reconstruct region-level slopes -----------------------------------------

region_slopes <- region_draws %>%
  mutate(
    lat     = b_cogyc_year_c        + r_region__cogyc,
    lon     = b_cogxc_year_c        + r_region__cogxc,
    depth   = b_depthnichec_year_c  + r_region__depthnichec,
    thermal = b_thermalnichec_year_c + r_region__thermalnichec
  ) %>%
  select(.draw, region, lat, lon, depth, thermal) %>%
  pivot_longer(cols = c(lat, lon, depth, thermal),
               names_to = "outcome", values_to = "slope") %>%
  mutate(
    warming_group = case_when(
      region %in% warming_regions     ~ "warming",
      region %in% non_warming_regions ~ "non_warming",
      TRUE                            ~ NA_character_
    )
  )


# per-draw group means and differences ------------------------------------

# each region weighted equally within its group
per_draw <- region_slopes %>%
  group_by(.draw, outcome, warming_group) %>%
  summarise(group_mean = mean(slope), .groups = "drop") %>%
  pivot_wider(names_from = warming_group, values_from = group_mean) %>%
  mutate(diff = warming - non_warming)


# summarise across draws --------------------------------------------------

summary_tbl <- per_draw %>%
  pivot_longer(cols = c(warming, non_warming, diff),
               names_to = "quantity", values_to = "value") %>%
  group_by(outcome, quantity) %>%
  summarise(
    median   = median(value),
    lower_95 = quantile(value, 0.025),
    upper_95 = quantile(value, 0.975),
    .groups  = "drop"
  ) %>%
  mutate(
    outcome  = factor(outcome,  levels = c("lat", "lon", "depth", "thermal")),
    quantity = factor(quantity, levels = c("warming", "non_warming", "diff"))
  ) %>%
  arrange(outcome, quantity)

print(summary_tbl, n = Inf)


# write LaTeX macros ------------------------------------------------------

mround <- function(x, digits) sprintf(paste0("%.", digits, "f"), round(x, digits))

out_file <- here("output/values/bayesian_trend_analysis/warming_group_slopes.tex")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
unlink(out_file)

write_tex <- function(x, macro, append = TRUE) {
  paste0("\\newcommand{\\", macro, "}{", x, "}") |>
    write_lines(out_file, append = append)
}

quantity_tag <- c(warming = "Warm", non_warming = "NonWarm", diff = "Diff")

for (i in seq_len(nrow(summary_tbl))) {
  oc <- tools::toTitleCase(as.character(summary_tbl$outcome[i]))
  qt <- quantity_tag[as.character(summary_tbl$quantity[i])]

  base <- paste0("WarmGrp", qt, oc)

  write_tex(mround(summary_tbl$median[i], 2), paste0(base, "Median"))
  write_tex(
    paste0("[95\\% CI: ",
           mround(summary_tbl$lower_95[i], 2), ", ",
           mround(summary_tbl$upper_95[i], 2), "]"),
    paste0(base, "CI")
  )
}


# regional temperature-trend correlations --------------------------------

# posterior-median regional slopes for each outcome (in per-decade units)
region_median_slopes <- region_slopes %>%
  group_by(outcome, region) %>%
  summarise(median_slope = median(slope), .groups = "drop")

# regional decadal temperature slopes (from prediction grid)
grid <- read_rds(here('R/data/processed/prediction_grid.rds'))

temp_summary <- grid %>%
  group_by(year, region_short) %>%
  summarise(mean_temp = mean(mean_year_temp, na.rm = TRUE), .groups = "drop") %>%
  filter(
    (region_short == "NIC"     & year >= 1994 & year <= 2019) |
    (region_short == "CBS"     & year >= 1994 & year <= 2020) |
    (region_short == "BAL"     & year >= 1999 & year <= 2020) |
    (region_short == "NS"      & year >= 1994 & year <= 2020) |
    (region_short == "BS"      & year >= 2004 & year <= 2021) |
    (region_short == "GOM"     & year >= 1994 & year <= 2020) |
    (region_short == "NEUS-SS" & year >= 1994 & year <= 2020) |
    (region_short == "USWC"    & year >= 2003 & year <= 2018) |
    (region_short == "BC"      & year >= 2003 & year <= 2019) |
    (region_short == "GOA"     & year >= 1996 & year <= 2019) |
    (region_short == "EBS"     & year >= 1994 & year <= 2019)
  ) %>%
  group_by(region_short) %>%
  summarise(temp_slope = coef(lm(mean_temp ~ year))[2] * 10, .groups = "drop")

# join per-outcome slopes with temperature slopes
cor_data <- region_median_slopes %>%
  left_join(temp_summary, by = c("region" = "region_short")) %>%
  mutate(outcome = factor(outcome, levels = c("lat", "lon", "depth", "thermal")))

cor_results <- cor_data %>%
  group_by(outcome) %>%
  summarise(
    r       = cor.test(temp_slope, median_slope)$estimate |> unname(),
    p_value = cor.test(temp_slope, median_slope)$p.value,
    ci_low  = cor.test(temp_slope, median_slope)$conf.int[1],
    ci_high = cor.test(temp_slope, median_slope)$conf.int[2],
    .groups = "drop"
  )

print(cor_results)

# write correlation macros
cor_file <- here("output/values/bayesian_trend_analysis/temp_correlations.tex")
dir.create(dirname(cor_file), recursive = TRUE, showWarnings = FALSE)
unlink(cor_file)

write_tex_cor <- function(x, macro, append = TRUE) {
  paste0("\\newcommand{\\", macro, "}{", x, "}") |>
    write_lines(cor_file, append = append)
}

for (i in seq_len(nrow(cor_results))) {
  oc      <- tools::toTitleCase(as.character(cor_results$outcome[i]))
  r_str   <- mround(cor_results$r[i], 2)
  is_tiny <- cor_results$p_value[i] < 0.001
  p_str   <- if (is_tiny) "$<$0.001" else mround(cor_results$p_value[i], 2)
  p_math  <- if (is_tiny) "< 0.001" else paste0("= ", mround(cor_results$p_value[i], 2))

  write_tex_cor(r_str, paste0("Cor", oc, "R"))
  write_tex_cor(p_str, paste0("Cor", oc, "P"))
  write_tex_cor(paste0("$r = ", r_str, "$, $p ", p_math, "$"), paste0("Cor", oc))
  write_tex_cor(
    paste0("[95\\% CI: ",
           mround(cor_results$ci_low[i], 2), ", ",
           mround(cor_results$ci_high[i], 2), "]"),
    paste0("Cor", oc, "CI")
  )
}


# plot --------------------------------------------------------------------

outcome_titles <- c(
  lat     = "**Latitudinal shift<br>(km decade<sup>-1</sup>)**",
  lon     = "**Longitudinal shift<br>(km decade<sup>-1</sup>)**",
  depth   = "**Depth shift<br>(m decade<sup>-1</sup>)**",
  thermal = "**Occupied temperature<br>change (&deg;C decade<sup>-1</sup>)**"
)

group_colors <- c(
  "Warming (n = 5)"     = "#B2182B",
  "Non-warming (n = 6)" = "#2166AC",
  "Difference"          = "grey30"
)

plot_data <- per_draw %>%
  pivot_longer(cols = c(warming, non_warming, diff),
               names_to = "quantity", values_to = "value") %>%
  mutate(
    outcome = factor(outcome, levels = c("lat", "lon", "depth", "thermal")),
    quantity_label = recode(quantity,
                            warming     = "Warming (n = 5)",
                            non_warming = "Non-warming (n = 6)",
                            diff        = "Difference"),
    quantity_label = factor(quantity_label,
                            levels = rev(c("Warming (n = 5)",
                                           "Non-warming (n = 6)",
                                           "Difference")))
  )

make_dist_panel <- function(oc, show_y_axis = FALSE) {
  d <- plot_data %>% filter(outcome == oc)

  p <- ggplot(d, aes(x = value, y = quantity_label,
                     fill = quantity_label, color = quantity_label)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    tidybayes::stat_halfeye(
      .width        = 0.95,
      alpha         = 0.75,
      point_size    = 1.3,
      interval_size = 1.0,
      normalize     = "xy",
      show.legend   = FALSE
    ) +
    scale_fill_manual(values  = group_colors) +
    scale_color_manual(values = group_colors) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    labs(title = outcome_titles[[oc]], x = NULL, y = NULL) +
    theme_bw(base_size = 9, base_family = "sans") +
    theme(
      plot.title        = ggtext::element_markdown(hjust = 0.5, size = 9,
                                                   family = "sans"),
      panel.grid        = element_blank(),
      axis.ticks        = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.title        = element_text(size = 9),
      axis.text         = element_text(size = 8),
      axis.text.y       = element_text(face = "bold", size = 8),
      plot.margin       = margin(2, 6, 2, 6),
      panel.border      = element_blank(),
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.background   = element_rect(fill = "white", colour = NA)
    )

  if (!show_y_axis) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }
  p
}

# correlation scatter panels ---------------------------------------------

warming_lookup <- tibble(
  region = c(warming_regions, non_warming_regions),
  warming_group = c(rep("warming", length(warming_regions)),
                    rep("non_warming", length(non_warming_regions)))
)

cor_plot_data <- cor_data %>%
  left_join(warming_lookup, by = "region")

cor_ann <- cor_results %>%
  mutate(
    label = ifelse(p_value < 0.001,
                   paste0("italic(r) == ", mround(r, 2), " * ',' ~ italic(p) < 0.001"),
                   paste0("italic(r) == ", mround(r, 2), " * ',' ~ italic(p) == ", mround(p_value, 2)))
  )

cor_y_titles <- c(
  lat     = "Latitudinal shift<br>(km decade<sup>-1</sup>)",
  lon     = "Longitudinal shift<br>(km decade<sup>-1</sup>)",
  depth   = "Depth shift<br>(m decade<sup>-1</sup>)",
  thermal = "Occupied temperature<br>change (&deg;C decade<sup>-1</sup>)"
)

make_cor_panel <- function(oc) {
  d <- cor_plot_data %>% filter(outcome == oc)
  ann_row <- cor_ann %>% filter(outcome == oc)

  ggplot(d, aes(x = temp_slope, y = median_slope)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_smooth(method = "lm", se = FALSE,
                color = "grey40", linewidth = 0.4) +
    geom_point(aes(fill = warming_group),
               size = 1.8, shape = 21, colour = "black", stroke = 0.3) +
    ggrepel::geom_text_repel(aes(label = region),
                             size = 7 / .pt, family = "sans",
                             color = "grey20",
                             segment.size = 0.2, segment.color = "grey60",
                             min.segment.length = 0.1, box.padding = 0.25,
                             seed = 1, max.overlaps = Inf,
                             show.legend = FALSE) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.4,
             label = ann_row$label, parse = TRUE,
             size = 7 / .pt, family = "sans", color = "grey20") +
    scale_fill_manual(values = c(warming = "#B2182B", non_warming = "#2166AC")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
    labs(x = "Temperature trend<br>(&deg;C decade<sup>-1</sup>)",
         y = cor_y_titles[[oc]]) +
    theme_bw(base_size = 9, base_family = "sans") +
    theme(
      panel.grid        = element_blank(),
      axis.ticks        = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.title.x      = ggtext::element_markdown(size = 9, hjust = 0.5),
      axis.title.y      = ggtext::element_markdown(size = 9, hjust = 0.5),
      axis.text         = element_text(size = 8),
      plot.margin       = margin(2, 6, 2, 6),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.3),
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.background   = element_rect(fill = "white", colour = NA),
      legend.position   = "none"
    )
}


# combine paired panels row-by-row ---------------------------------------

# 4 outcome rows, each with the group-density plot on the LEFT and the
# temperature-correlation scatter on the RIGHT.
outcome_names <- names(outcome_titles)

row_pairs <- purrr::map(outcome_names, function(oc) {
  left  <- make_dist_panel(oc, show_y_axis = TRUE)
  right <- make_cor_panel(oc)
  cowplot::plot_grid(left, right, nrow = 1, align = "h", axis = "tb",
                     rel_widths = c(1.1, 1))
})

final_plot <- cowplot::plot_grid(
  plotlist = row_pairs,
  ncol = 1,
  labels = c("a", "b", "c", "d"),
  label_size = 10, label_fontfamily = "sans",
  hjust = -0.4, vjust = 1.3
)

fig_path <- here("output/figures/supp/warming_group_slopes.png")

ggsave(fig_path, plot = final_plot,
       width = 178, height = 220,
       dpi = 600, units = "mm", bg = "white")
