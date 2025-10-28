# libraries ----------------------------------------------------------
library(tidyverse)
library(here)
library(brms)
library(tidybayes)
library(patchwork)
library(ggnewscale)
library(paletteer)
library(ggtext)
library(ggstats)
library(cowplot)

# read model ---------------------------------------------------------
m <- read_rds(here('R/bayesian_trends/fitted/m_stud.rds'))

# outcome names and pretty labels -----------------------------------
outcomes <- c("cogyc", "cogxc", "depthnichec", "thermalnichec")
outcomes_pretty <- c(
  "**Latitudinal shift<br>(km decade<sup>-1</sup>)**",
  "**Longitudinal shift<br>(km decade<sup>-1</sup>)**",
  "**Depth shift<br>(m decade<sup>-1</sup>)**",
  "**Thermal niche shift<br>(°C decade<sup>-1</sup>)**"
)

# color scales ------------------------------------------------------
fill_gradients <- list(
  cogyc = scale_fill_gradient2(low = "#762A83", mid = "#E0E0E0", high = "#1B7837", midpoint = 0),
  cogxc = scale_fill_gradient2(low = "#543005", mid = "#CFCFCF", high = "#003C30", midpoint = 0),
  depthnichec = scale_fill_gradient2(low = "#8C510A", mid = "#CFCFCF", high = "#01665E", midpoint = 0),
  thermalnichec = scale_fill_gradient2(low = "#2166AC", mid = "#DDDDDD", high = "#B2182B", midpoint = 0)
)
col_gradients <- list(
  cogyc = scale_color_gradient2(low = "#762A83", mid = "#E0E0E0", high = "#1B7837", midpoint = 0),
  cogxc = scale_color_gradient2(low = "#543005", mid = "#CFCFCF", high = "#003C30", midpoint = 0),
  depthnichec = scale_color_gradient2(low = "#8C510A", mid = "#CFCFCF", high = "#01665E", midpoint = 0),
  thermalnichec = scale_color_gradient2(low = "#2166AC", mid = "#DDDDDD", high = "#B2182B", midpoint = 0)
)

# global slopes -----------------------------

global_slopes <- m %>%
  spread_draws(b_cogyc_year_c, b_cogxc_year_c, b_depthnichec_year_c, b_thermalnichec_year_c)

global_slopes_long <- global_slopes %>%
  pivot_longer(
    cols = starts_with("b_"),
    names_to = "parameter",
    values_to = "slope"
  ) %>%
  mutate(
    outcome = str_remove(parameter, "b_"),
    outcome = str_remove(outcome, "_year_c")
  )

xlims_df <- global_slopes_long %>%
  group_by(outcome) %>%
  summarise(max_abs = max(abs(slope)), .groups = "drop") %>%
  mutate(max_abs = ifelse(max_abs == 0, 0.01, max_abs),
         xlim_low = -max_abs, xlim_high = max_abs)

global_plot_list <- list()
for (i in seq_along(outcomes)) {
  outcome_name <- outcomes[i]
  outcome_title <- outcomes_pretty[i]
  lims <- xlims_df %>% filter(outcome == outcome_name)
  
  p <- global_slopes_long %>%
    filter(outcome == outcome_name) %>%
    ggplot(aes(x = slope, y = 1)) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.4, color ='grey40') +
    tidybayes::stat_halfeye(
      aes(fill = after_stat(x)),
      .width = 0.95,
      normalize = 'all',
      alpha = 0.9,
      fill_type = "gradient",
      show.legend = FALSE,
      point_size = 1.3,
      interval_size = 1.2
    ) +
    fill_gradients[[outcome_name]] +
    scale_x_continuous(limits = c(lims$xlim_low, lims$xlim_high), breaks = scales::pretty_breaks(n = 5)) +
    scale_y_discrete(expand = c(0, 0.05)) +
    labs(title = outcome_title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = ggtext::element_markdown(hjust = 0.5, size = 10),
      panel.grid = element_blank(),
      axis.ticks.x = element_line(color = "black", size = 0.3),
      axis.ticks.length.x = unit(2, "pt"),
      plot.margin = margin(2, 5, 2, 5)
    )
  global_plot_list[[i]] <- p
}
p_global <- wrap_plots(global_plot_list, nrow = 1) +
  plot_layout(axis = "collect") &
  theme(panel.border = element_blank(), panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))

# region slopes -----------------------------

region_slopes <- m %>%
  spread_draws(
    b_cogyc_year_c, b_cogxc_year_c, b_depthnichec_year_c, b_thermalnichec_year_c,
    r_region__cogyc[region, year_c],
    r_region__cogxc[region, year_c],
    r_region__depthnichec[region, year_c],
    r_region__thermalnichec[region, year_c]
  )

region_slopes_long <- region_slopes %>%
  mutate(
    slope_cogyc = b_cogyc_year_c + r_region__cogyc,
    slope_cogxc = b_cogxc_year_c + r_region__cogxc,
    slope_depthnichec = b_depthnichec_year_c + r_region__depthnichec,
    slope_thermalnichec = b_thermalnichec_year_c + r_region__thermalnichec
  ) %>%
  select(.draw, region, starts_with("slope_")) %>%
  pivot_longer(cols = starts_with("slope_"), names_to = "outcome", names_prefix = "slope_", values_to = "slope")

xlims_region <- region_slopes_long %>%
  group_by(outcome) %>%
  summarise(max_abs = max(quantile(abs(slope), probs = .999), na.rm = TRUE), .groups = "drop") %>%
  mutate(max_abs = ifelse(max_abs == 0, 0.01, max_abs),
         xlim_low = -max_abs, xlim_high = max_abs)

region_plot_list <- list()
for (i in seq_along(outcomes)) {
  outcome_name <- outcomes[i]
  lims <- xlims_region %>% filter(outcome == outcome_name)
  
  p <- region_slopes_long %>%
    filter(outcome == outcome_name) %>%
    group_by(region) %>% mutate(median_slope = median(slope)) %>% ungroup() %>%
    mutate(region = factor(region, levels = rev(c('EBS','GOA','BC','COW','NEUS','GMX','BS','NS','CBS', 'BAL', 'NIC')))) %>%
    ggplot(aes(x = slope, y = region, fill = slope)) +
    ggstats::geom_stripped_rows(aes(y = region), odd = "white", even = "grey90", alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
    stat_pointinterval(aes(color = median_slope), .width = 0.95, point_size = 0.7, interval_size = 0.2, show.legend = FALSE) +
    col_gradients[[outcome_name]] +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    coord_cartesian(xlim = c(lims$xlim_low, lims$xlim_high)) +
    labs(x = NULL, y = "Region") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), legend.position = "none",
          axis.ticks = element_line(color = "black", size = 0.3),
          axis.ticks.length = unit(2, "pt"),
          plot.margin = margin(2, 5, 2, 5))
  region_plot_list[[outcome_name]] <- p
}

p_region <- wrap_plots(region_plot_list, nrow = 1) +
  plot_layout(axis_titles = "collect", axes = "collect_y") &
  theme(panel.border = element_blank(), panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))


# species slopes ----------------------------------------------------------

species_slopes <- m %>%
  spread_draws(
    `r_region:species__cogyc`[region_species, year_c],
    `r_region:species__cogxc`[region_species, year_c],
    `r_region:species__depthnichec`[region_species, year_c],
    `r_region:species__thermalnichec`[region_species, year_c]
  ) %>%
  mutate(region = sub("_.*", "", region_species))

species_slopes_long <- species_slopes %>%
  left_join(region_slopes, by = c(".draw", "region", "year_c")) %>%
  mutate(
    slope_cogyc   = b_cogyc_year_c + r_region__cogyc + `r_region:species__cogyc`,
    slope_cogxc   = b_cogxc_year_c + r_region__cogxc + `r_region:species__cogxc`,
    slope_depthnichec   = b_depthnichec_year_c + r_region__depthnichec + `r_region:species__depthnichec`,
    slope_thermalnichec = b_thermalnichec_year_c + r_region__thermalnichec + `r_region:species__thermalnichec`
  ) %>%
  pivot_longer(cols = starts_with("slope_"), names_to = "outcome", names_prefix = "slope_", values_to = "full_slope") %>%
  mutate(region = factor(region, levels = c('EBS','GOA','BC','COW','NEUS','GMX','BS','NS','CBS', 'BAL', 'NIC'))) %>%
  group_by(region, region_species) %>% mutate(sp_id = cur_group_id()) %>% ungroup() %>%
  mutate(sp_id = factor(sp_id, levels = sort(unique(sp_id), decreasing = TRUE)),
         outcome = factor(outcome, levels = outcomes))

xlims_species <- species_slopes_long %>%
  group_by(outcome) %>%
  summarise(max_abs = max(quantile(abs(full_slope), probs = .999), na.rm = TRUE), .groups = "drop") %>%
  mutate(max_abs = ifelse(max_abs == 0, 0.01, max_abs),
         xlim_low = -max_abs, xlim_high = max_abs)

every5 <- levels(species_slopes_long$sp_id)[as.integer(levels(species_slopes_long$sp_id)) %% 5 == 0]

species_plot_list <- list()
for (i in seq_along(outcomes)) {
  outcome_name <- outcomes[i]
  lims <- xlims_species %>% filter(outcome == outcome_name)
  species_data <- species_slopes_long %>% filter(outcome == outcome_name) %>% group_by(region, region_species) %>% mutate(median_slope = median(full_slope)) %>% ungroup()
  
  p <- species_data %>%
    ggplot(aes(x = full_slope, y = sp_id)) +
    ggstats::geom_stripped_rows(aes(y = sp_id), odd = "white", even = "grey90", alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
    stat_pointinterval(aes(color = median_slope), .width = 0.95, point_size = 0.1, interval_size = 0.1, show.legend = FALSE) +
    col_gradients[[outcome_name]] +
    scale_y_discrete(breaks = every5, labels = every5) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    coord_cartesian(xlim = c(lims$xlim_low, lims$xlim_high)) +
    facet_grid(region ~ ., scales = "free", space = "free", switch = "y") +
    labs(x = "Posterior estimate", y = "Species ID") +
    theme_minimal(base_size = 12) +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, size = 8),
          axis.text.y = element_text(size = 8),
          panel.spacing.y = unit(0.1, "cm"),
          strip.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          panel.grid = element_blank(),
          axis.ticks = element_line(color = "black", size = 0.3),
          axis.ticks.length = unit(2, "pt"),
          plot.margin = margin(2, 5, 2, 5))
  
  if (i > 1) p <- p + theme(strip.text.y.left = element_blank(), strip.background = element_blank())
  species_plot_list[[outcome_name]] <- p
}

p_species <- wrap_plots(species_plot_list, nrow = 1) +
  plot_layout(axis_titles = "collect", axes = "collect_y") &
  theme(panel.border = element_blank(), panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))


# final plot --------------------------------------------------------------

cowplot::plot_grid(
  p_global, p_region, p_species,
  align = "v", ncol = 1, axis = 'l',
  rel_heights = c(2, 2.5, 10),
  labels = c("a", "b", "c"), label_size = 12
)

# save
ggsave(
  here('output/figures/main/posterior_slopes.png'),
  width = 180,
  height = 250,
  dpi = 600,
  units = "mm",
  bg = "white"
)


# redo panel C for supp info ---------------------------------------------------

species_plot_list_supp <- list()

for (i in seq_along(outcomes)) {
  outcome_name <- outcomes[i]
  outcome_title <- outcomes_pretty[i]
  lims <- xlims_species %>% filter(outcome == outcome_name)
  
  species_data <- species_slopes_long %>%
    filter(outcome == outcome_name) %>%
    mutate(
      # extract species name from region_species
      species = region_species %>%
        str_remove("^[^_]+_") %>%
        str_replace_all("_", " "),
      # create parsed label: number + italic species
      species_label = paste0(sp_id, "~italic('", species, "')")
    ) %>%
    group_by(region, region_species) %>%
    mutate(median_slope = median(full_slope)) %>%
    ungroup()
  
  p <- species_data %>%
    ggplot(aes(x = full_slope, y = sp_id)) +  # keep order by sp_id
    ggstats::geom_stripped_rows(aes(y = sp_id),
                                odd = "white", even = "grey90", alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    stat_pointinterval(aes(color = median_slope),
                       .width = 0.95, point_size = 0.1, interval_size = 0.1,
                       show.legend = FALSE) +
    col_gradients[[outcome_name]] +
    # map labels manually while keeping order
    scale_y_discrete(
      labels = function(x) parse(text = species_data$species_label[match(x, species_data$sp_id)])
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    coord_cartesian(xlim = c(lims$xlim_low, lims$xlim_high)) +
    facet_grid(region ~ ., scales = "free", space = "free", switch = "y") +
    labs(title = outcome_title, x = "Posterior estimate", y = "Species") +
    theme_minimal(base_size = 9) +
    theme(
      strip.placement = "outside",
      plot.title = ggtext::element_markdown(hjust = 0.5, size = 8),
      strip.text.y.left = element_text(angle = 0, size = 5),
      axis.text.y = element_text(size = 3),
      panel.spacing.y = unit(0.1, "cm"),
      strip.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(2, "pt"),
      plot.margin = margin(2, 5, 2, 5)
    )
  
  if (i > 1)
    p <- p + theme(strip.text.y.left = element_blank(),
                   strip.background = element_blank())
  
  species_plot_list_supp[[outcome_name]] <- p
}

p_species_supp <- wrap_plots(species_plot_list_supp, nrow = 1) +
  plot_layout(axis_titles = "collect", axes = "collect_y") &
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))

p_species_supp

# save
ggsave(
  here('output/figures/supp/posterior_slopes_supp.png'),
  width = 180,
  height = 220,
  dpi = 600,
  units = "mm",
  bg = "white"
)


# prop of significant shifts per region -----------------------------------

species_summary <- species_slopes_long %>%
  group_by(outcome, region_species, sp_id) %>% # Group by outcome and region_species
  summarise(
    median_slope = median(full_slope, na.rm = TRUE), # median
    mean_slope = mean(full_slope, na.rm = TRUE), # mean
    lower_95 = quantile(full_slope, 0.025, na.rm = TRUE), # 2.5% quantile
    upper_95 = quantile(full_slope, 0.975, na.rm = TRUE), # 97.5% quantile
    .groups = "drop"
  ) %>%
  mutate(
    significant = case_when(
      lower_95 > 0 | upper_95 < 0 ~ "yes",   # credible interval does not include zero
      TRUE ~ "no"
    )
  ) %>%  mutate(
    region = sub("_.*", "", region_species),
    species = sub(".*?_", "", region_species),
    direction_label = case_when(
      significant == "no" ~ "Not significant",
      
      outcome == "cogyc" & median_slope > 0 ~ "Northing",
      outcome == "cogyc" & median_slope < 0 ~ "Southing",
      outcome == "cogxc" & median_slope > 0 ~ "Easting",
      outcome == "cogxc" & median_slope < 0 ~ "Westing",
      
      outcome == "depthnichec" & median_slope > 0 ~ "Deepening",
      outcome == "depthnichec" & median_slope < 0 ~ "Shallowing",
      outcome == "thermalnichec" & median_slope > 0 ~ "Warming",
      outcome == "thermalnichec" & median_slope < 0 ~ "Cooling",
      TRUE ~ "Other"
    )
  )

all_directions <- c(
  "Northing",
  "Southing",
  "Easting",
  "Westing",
  "Deepening",
  'Shallowing',
  "Warming",
  "Cooling",
  "Not significant"
)

# total number of species-region combinations per region
total_counts <- species_summary %>%
  distinct(region, species) %>%
  count(region, name = "total_species")

# count by region, outcome, direction label
summary_data <- species_summary %>%
  mutate(direction_label = factor(direction_label, levels = all_directions)) %>%
  group_by(region, outcome, direction_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(region, outcome, direction_label = all_directions, fill = list(n = 0)) %>%
  left_join(total_counts, by = "region") %>%
  mutate(proportion = n / total_species)

# label outcomes for legend titles
summary_data <- summary_data %>%
  mutate(outcome_label = recode(outcome,
                                "cogyc" = "Latitudinal shift",
                                "cogxc" = "Longitudinal shift",
                                "depthnichec" = 'Depth shift',
                                "thermalnichec" = 'Thermal niche shift'
  ))

# define consistent outcome and label mapping
outcome_order <- c(
  "Latitudinal shift", "Longitudinal shift", "Depth shift", "Thermal niche shift"
)


direction_map <- list(
  "Latitudinal shift"      = c("Northing", "Not significant", "Southing"),
  "Longitudinal shift"      = c("Easting", "Not significant", "Westing"),
  "Depth shift"            = c("Deepening", "Not significant", "Shallowing"),
  "Thermal niche shift"      = c("Warming", "Not significant", "Cooling")
)

summary_data = summary_data |> mutate(
  region = recode(
    region,
    "BAL" = "Baltic Sea",
    "BC" = 'British Columbia',
    "BS" = 'Barents Sea',
    "CBS" = 'Celtic-Biscay Shelf',
    "COW" = "U.S. West Coast", 
    "EBS" = 'Eastern Bering Sea',
    "GMX" = 'Gulf of Mexico',
    "GOA" = 'Gulf of Alaska',
    "NEUS" = 'NE US & Scotian Shelf',
    "NIC" = 'Northern Iberian Coast',
    "NS" = 'North Sea'
  )
)

summary_data <- summary_data %>%
  mutate(
    outcome_label = factor(outcome_label, levels = outcome_order)
  )


color_palette <- c(
  "Northing" = "#1B7837", "Southing" = "#762A83",
  "Easting" = "#003C30", "Westing" = "#543005",
  "Deepening" = "#01665E", "Shallowing" = "#8C510A",
  "Warming" = "#B2182B", "Cooling" = "#2166AC",
  "Not significant" = "grey50"
)

# base plot
p <- ggplot()

# loop over each outcome and layer with its own fill scale
for (i in seq_along(outcome_order)) {
  oc <- outcome_order[i]
  
  df_sub <- summary_data %>%
    filter(outcome_label == oc) %>%
    filter(direction_label %in% direction_map[[oc]]) %>%
    mutate(direction_label = factor(direction_label, levels = direction_map[[oc]]))
  
  valid_colors <- color_palette[names(color_palette) %in% direction_map[[oc]]]
  
  p <- p +
    geom_bar(
      data = df_sub,
      aes(x = outcome_label, y = proportion, fill = direction_label),
      stat = "identity",
      position = "stack", alpha = 0.5
    ) + geom_text(
      data = df_sub,
      aes(
        x = outcome_label,
        y = proportion,
        label = ifelse(proportion > 0, scales::percent(proportion, accuracy = 1), NA),
        group = direction_label
      ),
      stat = "identity",
      position = position_stack(vjust = .5),
      # centers labels in each segment
      size = 2.5,
      color = "black"
    ) +
    scale_fill_manual(
      values = valid_colors[direction_map[[oc]]],
      name = oc,
      drop = FALSE,
      guide = guide_legend(order = i)  # force the order here
    ) +
    scale_y_continuous(labels = NULL) +
    new_scale_fill() + coord_cartesian(expand = FALSE) + coord_flip() + scale_x_discrete(limits = rev(outcome_order))
}
# add theme and facet
p +
  facet_wrap(~ region,ncol=3) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    # legend.position = "bottom",
    # legend.box = "vertical",    
    strip.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),  
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(.1, "lines"),
    legend.title = element_text(size = 7) 
  ) +
  labs(x = "", y = "Species (%) per directional shift") 

# save
ggsave(
  here('output/figures/supp/prop_significant_supp.png'),
  width = 180,
  height = 150,
  dpi = 600,
  units = "mm",
  bg = "white"
)


# save slopes -------------------------------------------------------------

# a dataset with species-specific trend slopes

species_summary <- species_summary %>%
  mutate(
    outcome = recode(
      outcome,
      "cogyc" = "lat_shift",
      "cogxc" = "lon_shift",
      "thermalnichec" = "thermal_shift",
      "depthnichec" = "depth_shift"
    )
  ) %>%
  select(-region_species) %>%   # remove redundant column
  relocate(outcome, region, species, sp_id, .before = median_slope) %>% arrange(desc(sp_id))

# save as rds and csv
write_rds(species_summary, here('R/data/processed/bayesian_species_trends.rds'))
write_csv(species_summary, here('output/supp_data/bayesian_species_trends.csv'))









# from here ---------------------------------------------------------------






species_slopes_summary <- species_slopes_long %>%
  group_by(region, region_species, sp_id, outcome) %>%
  summarise(
    mean_slope = mean(full_slope),
    median_slope = median(full_slope),
    lower_95 = quantile(full_slope, 0.025),
    upper_95 = quantile(full_slope, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    significant = if_else(lower_95 > 0 | upper_95 < 0, "Yes", "No"),
    direction = case_when(
      significant == "Yes" & mean_slope > 0 ~ "Positive",
      significant == "Yes" & mean_slope < 0 ~ "Negative",
      TRUE ~ "Not significant"
    )
  ) %>%
  group_by(outcome) %>%
  summarize(
    # Species with min and max slope
    most_negative_species = region_species[which.min(mean_slope)],
    most_negative_value   = min(mean_slope),
    most_positive_species = region_species[which.max(mean_slope)],
    most_positive_value   = max(mean_slope),
    
    # Proportions
    n_species = n(),
    n_positive = sum(direction == "Positive"),
    n_negative = sum(direction == "Negative"),
    prop_positive = n_positive / n_species,
    prop_negative = n_negative / n_species,
    total_prop_significant = (n_positive + n_negative) / n_species,
    
    .groups = "drop"
  )













# summarise slopes for results --------------------------------------------

#global
global_slopes_summary <- global_slopes_long %>%
  group_by(outcome, parameter) %>%
  summarise(
    mean_slope = mean(slope),
    lower_95 = quantile(slope, 0.025),
    upper_95 = quantile(slope, 0.975),
    .groups = "drop"
  )

#region
region_slopes_summary <- region_slopes_long %>%
  group_by(region, outcome) %>%
  summarise(
    mean_slope = mean(slope),
    lower_95 = quantile(slope, 0.025),
    upper_95 = quantile(slope, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    significant = if_else(lower_95 > 0 | upper_95 < 0, "Yes", "No")
    )


prop_significant_region <- species_slopes_long %>%
  group_by(region, region_species, outcome) %>%
  summarise(
    mean_slope = mean(full_slope),
    lower_95 = quantile(full_slope, 0.025),
    upper_95 = quantile(full_slope, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    significant = if_else(lower_95 > 0 | upper_95 < 0, "Yes", "No"),
    direction = case_when(
      significant == "Yes" & mean_slope > 0 ~ "Positive",
      significant == "Yes" & mean_slope < 0 ~ "Negative",
      TRUE ~ "Not significant"
    )
  ) %>%
  group_by(outcome,region) %>%
  summarize(
    # Proportions
    n_species = n(),
    n_positive = sum(direction == "Positive"),
    n_negative = sum(direction == "Negative"),
    prop_positive = n_positive / n_species,
    prop_negative = n_negative / n_species,
    total_prop_significant = (n_positive + n_negative) / n_species,
    .groups = "drop"
  )

write_csv(prop_significant_region, here('bayesian_trends/tables/prop_by_region.csv'))

species_slopes_summary <- species_slopes_long %>%
  group_by(region, region_species, outcome) %>%
  summarise(
    mean_slope = mean(full_slope),
    lower_95 = quantile(full_slope, 0.025),
    upper_95 = quantile(full_slope, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    significant = if_else(lower_95 > 0 | upper_95 < 0, "Yes", "No"),
    direction = case_when(
      significant == "Yes" & mean_slope > 0 ~ "Positive",
      significant == "Yes" & mean_slope < 0 ~ "Negative",
      TRUE ~ "Not significant"
    )
  ) %>%
  group_by(outcome) %>%
  summarize(
    # Species with min and max slope
    most_negative_species = region_species[which.min(mean_slope)],
    most_negative_value   = min(mean_slope),
    most_positive_species = region_species[which.max(mean_slope)],
    most_positive_value   = max(mean_slope),
    
    # Proportions
    n_species = n(),
    n_positive = sum(direction == "Positive"),
    n_negative = sum(direction == "Negative"),
    prop_positive = n_positive / n_species,
    prop_negative = n_negative / n_species,
    total_prop_significant = (n_positive + n_negative) / n_species,
    
    .groups = "drop"
  )


write_csv(species_slopes_summary, here('bayesian_trends/tables/species_slopes_summary.csv'))


# export slopes for processing --------------------------------------------

bayesian_species_trends <- species_slopes_long %>%
  group_by(region, region_species, outcome) %>%
  summarise(
    mean_slope = mean(full_slope),
    lower_95 = quantile(full_slope, 0.025),
    upper_95 = quantile(full_slope, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    significant = if_else(lower_95 > 0 | upper_95 < 0, "Yes", "No"),
    direction = case_when(
      significant == "Yes" & mean_slope > 0 ~ "Positive",
      significant == "Yes" & mean_slope < 0 ~ "Negative",
      TRUE ~ "Not significant"
    )
  ) %>%
  separate(
    region_species,
    into = c("region", "species"),
    sep = "_",
    extra = "merge",
    remove = FALSE
  )

write_rds(bayesian_species_trends, here('data/processed/bayesian_species_trends.rds'))


# export species-region slopes ---------------------------------------------
bayesian_species_trends <- species_slopes_long %>%
  group_by(region, region_species, outcome) %>%
  summarise(
    mean_slope = mean(full_slope),
    lower_95 = quantile(full_slope, 0.025),
    upper_95 = quantile(full_slope, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    significant = if_else(lower_95 > 0 | upper_95 < 0, "Yes", "No"),
    direction = case_when(
      significant == "Yes" & mean_slope > 0 ~ "Positive",
      significant == "Yes" & mean_slope < 0 ~ "Negative",
      TRUE ~ "Not significant"
    )
  ) %>%
  separate(
    region_species,
    into = c("region", "species"),
    sep = "_",
    extra = "merge"
  ) %>%
  mutate(
    # Clean species names and format for LaTeX
    species = str_replace_all(species, "_", " "),
    species = str_to_title(species),
    Species = paste0("\\textit{", species, "}"),
    
    # Replace outcome codes with descriptive labels
    outcome = recode(
      outcome,
      cogyc = "Latitudinal shift",
      cogxc = "Longitudinal shift",
      depthnichec = "Depth shift",
      thermalnichec = "Thermal niche shift"
    ),
    
    # Round numeric columns to 2 decimals
    mean_slope = round(mean_slope, 2),
    lower_95 = round(lower_95, 2),
    upper_95 = round(upper_95, 2)
  ) %>%
  # Rename columns for LaTeX-friendly table
  rename(
    Region = region,
    `Type of shift` = outcome,
    Estimate = mean_slope,
    `Lower 95 % CI` = lower_95,
    `Upper 95 % CI` = upper_95,
    Significant = significant,
    Direction = direction
  ) %>%
  select(Region, Species, `Type of shift`, Estimate, `Lower 95 % CI`, `Upper 95 % CI`, Significant, Direction)

write_csv(bayesian_species_trends, here('data/processed/bayesian_species_trends.rds'))

# Export LaTeX table
kable(bayesian_species_trends,
      format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      caption = "Posterior estimates of species-level shifts. For each species, the table shows the estimated rate of change (Estimate) with 95% credible intervals (CI), whether the shift was statistically significant, and the direction of change. Outcomes include latitudinal, longitudinal, depth, and thermal niche shifts."
) %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  cat(file = here('bayesian_trends/tables/bayesian_species_trends.tex'))






# export species-specific slopes ------------------------------------------
library(knitr)
library(kableExtra)

species_slopes_summary <- species_slopes_long %>%
  group_by(region, region_species, outcome, sp_id) %>%
  summarise(
    mean_slope = mean(full_slope),
    lower_95 = quantile(full_slope, 0.025),
    upper_95 = quantile(full_slope, 0.975),
    .groups = "drop"
  ) %>%   mutate(outcome = recode(outcome,
                          cogyc = "Latitudinal shift",
                          cogxc = "Longitudinal shift",
                          depthnichec = "Depth shift",
                          thermalnichec = "Thermal niche shift"
  )) %>%
  # format mean ± CI
  mutate(value = paste0(round(mean_slope, 2), 
                        " [", round(lower_95, 2), "–", round(upper_95, 2), "]")) %>%
  # keep only species name (after first "_")
  mutate(species = str_remove(region_species, "^[^_]+_")) %>%
  # pivot outcomes to wide format
  select(region, species, sp_id, outcome, value) %>%
  pivot_wider(names_from = outcome, values_from = value) %>% rename(
    Region = region,
    Species = species,
    `Species ID` = sp_id
  ) %>%
  # replace underscores with spaces in species names
  mutate(Species = gsub("_", " ", Species))

# Create LaTeX table
# Create LaTeX table for multi-page
latex_table <- species_slopes_summary %>%
  kable(format = "latex",
        booktabs = TRUE,
        longtable = TRUE,          # enables multi-page
        caption = "Species Slopes Summary",
        escape = FALSE) %>%
  kable_styling(latex_options = c("repeat_header", "scale_down"))

# Export to .tex
cat(latex_table, file = here("bayesian_trends/tables/species_slopes_summary.tex"))
