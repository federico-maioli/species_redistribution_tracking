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
library(geomtextpath)
library(knitr)
library(kableExtra)

## definitely not the most efficient code you ever saw for plotting
#https://commons.wikimedia.org/wiki/Category:It_ain%27t_much,_but_it%27s_honest_work#/media/File:Farmer_meme_with_apostrophe.jpg

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
  depthnichec = scale_fill_gradient2(low = "#A6CEE3", mid = "#CFCFCF", high = "#003C8F", midpoint = 0),
  thermalnichec = scale_fill_gradient2(low = "#2166AC", mid = "#DDDDDD", high = "#B2182B", midpoint = 0)
)
col_gradients <- list(
  cogyc = scale_color_gradient2(low = "#762A83", mid = "#E0E0E0", high = "#1B7837", midpoint = 0),
  cogxc = scale_color_gradient2(low = "#543005", mid = "#CFCFCF", high = "#003C30", midpoint = 0),
  depthnichec = scale_color_gradient2(low = "#A6CEE3", mid = "#CFCFCF", high =  "#003C8F", midpoint = 0),
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
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.3, alpha = .3) +
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
    scale_x_continuous(limits = c(lims$xlim_low, lims$xlim_high), breaks = scales::pretty_breaks(n = 3)) +
    scale_y_discrete(expand = c(0, 0.05)) +
    labs(title = outcome_title, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = ggtext::element_markdown(hjust = 0.5, size = 10),
      panel.grid = element_blank(),
      axis.ticks.x = element_line(color = "black", size = 0.3),
      axis.ticks.length.x = unit(2, "pt"),
      axis.text.y = element_text(face = "bold"),#  axis.text.x = element_text(angle = 45, hjust = 1),
     plot.margin = margin(5, 9, 5, 9)
    )
  global_plot_list[[i]] <- p
}
p_global <- wrap_plots(global_plot_list, nrow = 1) +
  plot_layout(axis = "collect") &
  theme(panel.border = element_blank(), panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))

summary_global_slopes <- global_slopes_long  %>%  group_by(outcome) %>% # Group by outcome and region_species
  summarise(
    median_slope = median(slope, na.rm = TRUE), # median
    mean_slope = mean(slope, na.rm = TRUE), # mean
    lower_95 = quantile(slope, 0.025, na.rm = TRUE), # 2.5% quantile
    upper_95 = quantile(slope, 0.975, na.rm = TRUE), # 97.5% quantile
    .groups = "drop"
  )

write_rds(summary_global_slopes, here('R/data/processed/bayesian_global_trends.rds'))

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
    mutate(region = factor(region, levels = rev(c('EBS','GOA','BC','USWC','NEUS-SS','GOM','BS','NS','CBS', 'BAL', 'NIC')))) %>%
    ggplot(aes(x = slope, y = region, fill = slope)) +
    ggstats::geom_stripped_rows(aes(y = region), odd = "grey90", even = 'white' , alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.3, alpha = .3) +
    stat_pointinterval(aes(color = median_slope), .width = 0.95, point_size = 1.6, interval_size = 0.9, show.legend = FALSE) +
    col_gradients[[outcome_name]] +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    coord_cartesian(xlim = c(lims$xlim_low, lims$xlim_high)) +
    labs(x = NULL, y = "Region") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "none",
          axis.ticks = element_line(color = "black", size = 0.3),
          axis.ticks.length = unit(2, "pt"),
          #panel.spacing = unit(200, "lines"),
          plot.margin = margin(5, 9, 5, 9),
          #axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(face = "bold") )
  region_plot_list[[outcome_name]] <- p
}

p_region <- wrap_plots(region_plot_list, nrow = 1) +
  plot_layout(axis_titles = "collect", axes = "collect_y") &
  theme(panel.border = element_blank(), panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))

# export slopes 
summary_region_slopes <- region_slopes_long  %>%  group_by(outcome,region) %>% # Group by outcome and region_species
  summarise(
    median_slope = median(slope, na.rm = TRUE), # median
    mean_slope = mean(slope, na.rm = TRUE), # mean
    lower_95 = quantile(slope, 0.025, na.rm = TRUE), # 2.5% quantile
    upper_95 = quantile(slope, 0.975, na.rm = TRUE), # 97.5% quantile
    .groups = "drop"
  )

write_rds(summary_region_slopes, here('R/data/processed/bayesian_region_trends.rds'))


# combine -----------------------------------------------------------------

cowplot::plot_grid(
  p_global, p_region,
  align = "v", ncol = 1, axis = 'l',
  rel_heights = c(1,2.3),
  labels = c("a", "b"), label_size = 14
)

ggsave(
  here('output/figures/main/posterior_slopes.png'),
  width = 180,
  height = 140,
  dpi = 600,
  units = "mm",
  bg = "white"
)

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
  mutate(region = factor(region, levels = c('EBS','GOA','BC','USWC','NEUS-SS','GOM','BS','NS','CBS', 'BAL', 'NIC'))) %>%
  group_by(region, region_species) %>% mutate(sp_id = cur_group_id()) %>% ungroup() %>%
  mutate(sp_id = factor(sp_id, levels = sort(unique(sp_id), decreasing = TRUE)),
         outcome = factor(outcome, levels = outcomes))

# xlims_species <- species_slopes_long %>%
#   group_by(outcome) %>%
#   summarise(max_abs = max(quantile(abs(full_slope), probs = .999), na.rm = TRUE), .groups = "drop") %>%
#   mutate(max_abs = ifelse(max_abs == 0, 0.01, max_abs),
#          xlim_low = -max_abs, xlim_high = max_abs)
# 
# every5 <- levels(species_slopes_long$sp_id)[as.integer(levels(species_slopes_long$sp_id)) %% 5 == 0]
# 
# species_plot_list <- list()
# for (i in seq_along(outcomes)) {
#   outcome_name <- outcomes[i]
#   lims <- xlims_species %>% filter(outcome == outcome_name)
#   species_data <- species_slopes_long %>% filter(outcome == outcome_name) %>% group_by(region, region_species) %>% mutate(median_slope = median(full_slope)) %>% ungroup()
#   
#   p <- species_data %>%
#     ggplot(aes(x = full_slope, y = sp_id)) +
#     ggstats::geom_stripped_rows(aes(y = sp_id), odd = "white", even = "grey90", alpha = 0.2) +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
#     stat_pointinterval(aes(color = median_slope), .width = 0.95, point_size = 0.1, interval_size = 0.1, show.legend = FALSE) +
#     col_gradients[[outcome_name]] +
#     scale_y_discrete(breaks = every5, labels = every5) +
#     scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
#     coord_cartesian(xlim = c(lims$xlim_low, lims$xlim_high)) +
#     facet_grid(region ~ ., scales = "free", space = "free", switch = "y") +
#     labs(x = "Posterior estimate", y = "Species ID") +
#     theme_minimal(base_size = 12) +
#     theme(strip.placement = "outside",
#           strip.text.y.left = element_text(angle = 0, size = 8),
#           axis.text.y = element_text(size = 8),
#           panel.spacing.y = unit(0.1, "cm"),
#           strip.background = element_rect(fill = "white", color = NA),
#           plot.background = element_rect(fill = "white", color = NA),
#           panel.background = element_rect(fill = "white", color = NA),
#           panel.grid = element_blank(),
#           axis.ticks = element_line(color = "black", size = 0.3),
#           axis.ticks.length = unit(2, "pt"),
#           plot.margin = margin(2, 5, 2, 5))
#   
#   if (i > 1) p <- p + theme(strip.text.y.left = element_blank(), strip.background = element_blank())
#   species_plot_list[[outcome_name]] <- p
# }
# 
# p_species <- wrap_plots(species_plot_list, nrow = 1) +
#   plot_layout(axis_titles = "collect", axes = "collect_y") &
#   theme(panel.border = element_blank(), panel.background = element_rect(fill = "white", colour = NA),
#         plot.background  = element_rect(fill = "white", colour = NA))



# species-specfic slopes 4 supp info ---------------------------------------------------

xlims_species <- species_slopes_long %>%
  group_by(outcome) %>%
  summarise(max_abs = max(quantile(abs(full_slope), probs = .999), na.rm = TRUE), .groups = "drop") %>%
  mutate(max_abs = ifelse(max_abs == 0, 0.01, max_abs),
         xlim_low = -max_abs, xlim_high = max_abs)

species_plot_list_supp <- list()

for (i in seq_along(outcomes)) {
  outcome_name <- outcomes[i]
  outcome_title <- outcomes_pretty[i]
  lims <- xlims_species %>% filter(outcome == outcome_name)
  
  species_data <- species_slopes_long %>%
    filter(outcome == outcome_name) %>% #%in% c('NS')
    mutate(
      # extract species name from region_species
      species = region_species %>%
        str_remove("^[^_]+_") %>%
        str_replace_all("_", " "),
      # create parsed label: number + italic species
      species_label = paste0("~italic('", species, "')")
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
  height = 240,
  dpi = 600,
  units = "mm",
  bg = "white"
)


# one example for presentation --------------------------------------------

# xlims_species <- species_slopes_long %>%
#   group_by(outcome) %>%
#   summarise(max_abs = max(quantile(abs(full_slope), probs = .999), na.rm = TRUE), .groups = "drop") %>%
#   mutate(max_abs = ifelse(max_abs == 0, 0.01, max_abs),
#          xlim_low = -max_abs, xlim_high = max_abs)
# 
# species_plot_list_supp <- list()
# 
# for (i in seq_along(outcomes)) {
#   outcome_name <- outcomes[i]
#   outcome_title <- outcomes_pretty[i]
#   lims <- xlims_species %>% filter(outcome == outcome_name)
#   
#   species_data <- species_slopes_long %>%
#     filter(outcome == outcome_name,region %in% c('NS')) %>%
#     mutate(
#       # extract species name from region_species
#       species = region_species %>%
#         str_remove("^[^_]+_") %>%
#         str_replace_all("_", " "),
#       # create parsed label: number + italic species
#       species_label = paste0("~italic('", species, "')")
#     ) %>%
#     group_by(region, region_species) %>%
#     mutate(median_slope = median(full_slope)) %>%
#     ungroup()
#   
#   p <- species_data %>%
#     ggplot(aes(x = full_slope, y = sp_id)) +  # keep order by sp_id
#     ggstats::geom_stripped_rows(aes(y = sp_id),
#                                 odd = "white", even = "grey90", alpha = 0.2) +
#     geom_vline(xintercept = 0, linetype = "dashed",
#                color = "grey90", linewidth = 0.3) +
#     stat_pointinterval(aes(color = median_slope),
#                        .width = .95,
#                        point_size = .5, 
#                        interval_size = .5,
#                        show.legend = FALSE) +
#     col_gradients[[outcome_name]] +
#     # map labels manually while keeping order
#     scale_y_discrete(
#       labels = function(x) parse(text = species_data$species_label[match(x, species_data$sp_id)])
#     ) +
#     scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
#     #coord_cartesian(xlim = c(lims$xlim_low, lims$xlim_high)) +
#     facet_grid(region ~ ., scales = "free", space = "free", switch = "y") +
#     labs(title = outcome_title, x = "Posterior estimate", y = "Species") +
#     theme_minimal(base_size = 9) +
#     theme(
#       strip.placement = "outside",
#       plot.title = ggtext::element_markdown(hjust = 0.5, size = 8),
#       strip.text.y.left = element_text(angle = 0, size = 8),
#       axis.text.y = element_text(size = 7),
#       panel.spacing.y = unit(0.1, "cm"),
#       strip.background = element_rect(fill = "white", color = NA),
#       plot.background = element_rect(fill = "white", color = NA),
#       panel.background = element_rect(fill = "white", color = NA),
#       panel.grid = element_blank(),
#       axis.ticks = element_line(color = "black", size = 0.3),
#       axis.ticks.length = unit(2, "pt"),
#       plot.margin = margin(2, 5, 2, 5)
#     )
#   
#   if (i > 1)
#     p <- p + theme(strip.text.y.left = element_blank(),
#                    strip.background = element_blank())
#   
#   species_plot_list_supp[[outcome_name]] <- p
# }
# 
# p_species_supp <- wrap_plots(species_plot_list_supp, nrow = 1) +
#   plot_layout(axis_titles = "collect", axes = "collect_y") &
#   theme(panel.border = element_blank(),
#         panel.background = element_rect(fill = "white", colour = NA),
#         plot.background = element_rect(fill = "white", colour = NA))
# 
# p_species_supp
# 
# # save
# ggsave(
#   here('output/figures/supp/posterior_slopes_supp.png'),
#   width = 180,
#   height = 150,
#   dpi = 600,
#   units = "mm",
#   bg = "white"
# )



#  proportion of significant shifts plot ---------------------------------------

# species summary 
# compute median, mean, 95% CI per species x region x outcome
species_summary <- species_slopes_long %>%
  group_by(outcome, region_species, sp_id) %>%
  summarise(
    median_slope = median(full_slope, na.rm = TRUE),
    mean_slope   = mean(full_slope, na.rm = TRUE),
    lower_95     = quantile(full_slope, 0.025, na.rm = TRUE),
    upper_95     = quantile(full_slope, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # assign significance
    significant = ifelse(lower_95 > 0 | upper_95 < 0, "yes", "no"),
    # split region and species
    region  = sub("_.*", "", region_species),
    species = sub(".*?_", "", region_species),
    # assign direction labels
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

# all possible directions
all_directions <- c(
  "Northing","Southing","Easting","Westing",
  "Deepening","Shallowing","Warming","Cooling",
  "Not significant"
)

# total species per region
total_counts <- species_summary %>%
  distinct(region, species) %>%
  count(region, name = "total_species")

# overall counts 
total_counts_overall <- tibble(
  region = "Overall",
  total_species = sum(total_counts$total_species)
)

# summary data by region
summary_data <- species_summary %>%
  mutate(direction_label = factor(direction_label, levels = all_directions)) %>%
  group_by(region, outcome, direction_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(region, outcome, direction_label = all_directions, fill = list(n = 0)) %>%
  left_join(total_counts, by = "region") %>%
  mutate(proportion = n / total_species)

# summary data overall
summary_data_overall <- species_summary %>%
  mutate(direction_label = factor(direction_label, levels = all_directions)) %>%
  group_by(outcome, direction_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(outcome, direction_label = all_directions, fill = list(n = 0)) %>%
  mutate(region = "Overall") %>%
  left_join(total_counts_overall, by = "region") %>%
  mutate(proportion = n / total_species)

# outcome labeling 
summary_data <- summary_data %>%
  mutate(outcome_label = recode(outcome,
                                "cogyc" = "Latitude",
                                "cogxc" = "Longitude",
                                "depthnichec" = "Depth",
                                "thermalnichec" = "Thermal niche"))

summary_data_overall <- summary_data_overall %>%
  mutate(outcome_label = recode(outcome,
                                "cogyc" = "Latitude",
                                "cogxc" = "Longitude",
                                "depthnichec" = "Depth",
                                "thermalnichec" = "Thermal niche"))

# bind region + overall
summary_data <- bind_rows(summary_data, summary_data_overall)

# region labels and ordering
region_order <- c('EBS','GOA','BC','USWC','NEUS-SS','GOM','BS','NS','CBS','BAL','NIC')

summary_data <- summary_data %>%
  mutate(
    region_label = paste0(region, " (n = ", total_species, ")"),
    region_code  = ifelse(region == "Overall", "0", region),
    region_code  = factor(region_code, levels = c("0", region_order))
  ) %>%
  arrange(region_code) %>%
  mutate(region_label = factor(region_label, levels = unique(region_label[order(region_code)])))

# direction maps 
direction_map <- list(
  "Latitude"      = c("Northing", "Southing", "Not significant"),
  "Longitude"     = c("Easting", "Westing", "Not significant"),
  "Depth"         = c("Deepening", "Shallowing", "Not significant"),
  "Thermal niche" = c("Warming", "Cooling", "Not significant")
)

direction_map_tidy <- tibble(
  outcome_label   = rep(names(direction_map), lengths(direction_map)),
  direction_label = unlist(direction_map)
)

# filter relevant rows
summary_data_filtered <- summary_data %>%
  inner_join(direction_map_tidy, by = c("outcome_label", "direction_label"))

# color palette 
color_palette <- c(
  "Northing" = "#1B7837", "Southing" = "#762A83",
  "Easting" = "#003C30", "Westing" = "#543005",
  "Deepening" = "#003C8F", "Shallowing" = "#A6CEE3",
  "Warming"  = "#B2182B", "Cooling" = "#2166AC",
  "NS_Latitude" = "grey80", "NS_Longitude" = "grey80",
  "NS_Depth" = "grey80", "NS_Thermal niche" = "grey80"
)

# create NS-unique category + ordering
df <- summary_data_filtered %>%
  mutate(
    direction_label_unique = ifelse(
      direction_label == "Not significant",
      paste0("NS_", outcome_label),
      direction_label
    ),
    direction_label_unique = case_when(
      outcome_label == "Depth"         ~ factor(direction_label_unique,
                                                c("Shallowing", "NS_Depth", "Deepening")),
      outcome_label == "Latitude"      ~ factor(direction_label_unique,
                                                c("Southing", "NS_Latitude", "Northing")),
      outcome_label == "Longitude"     ~ factor(direction_label_unique,
                                                c("Westing", "NS_Longitude", "Easting")),
      outcome_label == "Thermal niche" ~ factor(direction_label_unique,
                                                c("Cooling", "NS_Thermal niche", "Warming"))
    )
  )

# add spacing rows manually
n_empty <- 2

df_spaced <- df %>%
  mutate(
    outcome_label = factor(outcome_label,
                           levels = c("Latitude", "Longitude", "Depth", "Thermal niche"))
  ) %>%
  arrange(region, outcome_label, direction_label_unique) %>%
  group_by(region, outcome_label) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup() %>%
  group_by(region, group_id) %>%
  group_modify(~ bind_rows(
    .x,
    tibble(
      outcome_label = NA,
      direction_label = NA,
      direction_label_unique = NA,
      proportion = 0,
      id = NA
    )[rep(1, n_empty), ]
  )) %>%
  ungroup() %>%
  group_by(region) %>%
  mutate(id = row_number()) %>%
  ungroup()

# angles and label justification
df_spaced <- df_spaced %>%
  group_by(region) %>%
  mutate(
    n_bar = n(),
    angle = 90 - 360 * (id - 0.5) / n_bar,
    hjust = ifelse(angle < -90, 1, 0),
    angle = ifelse(angle < -90, angle + 180, angle)
  ) %>%
  ungroup()

# base data for outcome label arcs
base_data_main <- df_spaced %>%
  filter(!is.na(outcome_label)) %>%
  group_by(region, outcome_label) %>%
  summarise(start = min(id), end = max(id), .groups = "drop") %>%
  mutate(title = (start + end) / 2)

# circular region plot function
plot_region_circular <- function(region_name) {
  df_region   <- df_spaced       %>% filter(region == region_name)
  base_region <- base_data_main  %>% filter(region == region_name)
  
  ggplot(df_region, aes(factor(id), proportion, fill = direction_label_unique)) +
    geom_bar(stat = "identity", width = 1, alpha = 0.7) +
    scale_fill_manual(values = color_palette) +
    geom_text(
      aes(
        label = ifelse(
          proportion > 0 & !is.na(proportion),
          ifelse(proportion < 0.01, "<1%", scales::percent(proportion, accuracy = 1)),
          ""
        ),
        angle = angle, hjust = hjust
      ),
      vjust = .5, size = 2
    ) +
    coord_polar() +
    ylim(-1.2, max(df_region$proportion, na.rm = TRUE) + 0.1) +
    theme_minimal(base_size = 7) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      plot.title = element_text(hjust = 0.5, vjust = -5, face = "bold")
    ) +
    geom_segment(
      data = base_region,
      aes(x = start-.5, y = 0, xend = end+.5, yend = 0),
      inherit.aes = FALSE, color = "black", linewidth = 0.3
    ) +
    geom_textpath(
      data = base_region,
      aes(x = title, y = -0.15, label = outcome_label),
      inherit.aes = FALSE,
      size = 2.1
    ) +
    ggtitle(region_name)
}

# build plot list
regions   <- unique(df_spaced$region)
plot_list <- setNames(lapply(regions, plot_region_circular), regions)

# custom legend 
legend_data <- tibble(
  outcome_label = rep(c("Latitude", "Longitude", "Depth", "Thermal niche"), each = 3),
  direction_label = c("Southing", "n.s.", "Northing",
                      "Westing", "n.s.", "Easting",
                      "Shallowing", "n.s.", "Deepening",
                      "Cooling", "n.s.", "Warming")
) %>%
  mutate(outcome_label = factor(
    outcome_label,
    levels = c("Latitude", "Longitude", "Depth", "Thermal niche")
  ))

legend_spaced <- legend_data %>%
  group_split(outcome_label) %>%
  map_df(~ bind_rows(.x,
                     tibble(outcome_label = unique(.x$outcome_label),
                            direction_label = NA)[rep(1, n_empty),])) %>%
  mutate(
    id    = row_number(),
    value = .5,
    angle = 90 - 360 * (id - 0.5) / n(),
    hjust = ifelse(angle < -90, 1, 0),
    angle = ifelse(angle < -90, angle + 180, angle)
  )

base_data_leg <- legend_spaced %>%
  filter(!is.na(direction_label)) %>%
  group_by(outcome_label) %>%
  summarise(start = min(id), end = max(id), .groups = "drop") %>%
  mutate(title = (start + end) / 2)

legend_colors <- c(
  "Southing" = "#762A83", "Northing" = "#1B7837",
  "Westing"  = "#543005", "Easting" = "#003C30",
  "Shallowing" = "#A6CEE3", "Deepening" = "#003C8F",
  "Cooling" = "#2166AC", "Warming" = "#B2182B",
  "n.s." = "grey80"
)

legend_plot <- ggplot(legend_spaced, aes(factor(id), value, fill = direction_label)) +
  geom_bar(stat = "identity", width = 1, alpha=.7) +
  scale_fill_manual(values = legend_colors, na.value = NA) +
  coord_polar() +
  ylim(-0.5, .8) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0), legend.position = 'none') +
  geom_text(
    data = legend_spaced %>% filter(!is.na(direction_label)),
    aes(label = direction_label, angle = angle, hjust = hjust),
    vjust = .5, size = 2.4, fontface = "bold"
  ) +
  geom_segment(
    data = base_data_leg,
    aes(x = start-.5, y = 0, xend = end+.5, yend = 0),
    inherit.aes = FALSE, color = "black", linewidth = .3
  ) +
  geom_textpath(
    data = base_data_leg,
    aes(x = title, y = -0.08, label = outcome_label),
    inherit.aes = FALSE,
    size = 2.7,
    fontface = "bold"
  )

# combine main figure + legend
p <- (plot_list[['Overall']] + plot_list[['EBS']] + plot_list[["GOA"]] + plot_list[['BC']] + plot_list[['USWC']] +
        plot_list[['NEUS-SS']] + plot_list[['GOM']] + plot_list[['BS']] + free(legend_plot) + plot_list[['NS']] +
        plot_list[['CBS']] + plot_list[['BAL']] + plot_list[['NIC']] ) + 
  plot_layout(
    design = "
    ABCD
    EFGH
    IIKL
    IIMN
    "
  )

p

# save figure 
ggsave(
  here("output/figures/main/prop_significant.png"),
  width = 180, height = 190,
  dpi = 600, units = "mm", bg = "white"
)



# save slopes for supp data-------------------------------------------------------------
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


supp_table <- species_summary %>%
  mutate(
    estimate = sprintf(
      "%.2f [%.2f, %.2f]",
      median_slope, lower_95, upper_95
    )
  ) %>%
  select(
    Region = region,
    Species = species,
    Outcome = outcome,
    Estimate = estimate
  )

supp_table_wide <- supp_table %>%
  pivot_wider(
    names_from = Outcome,
    values_from = Estimate
  )

supp_table_wide <- supp_table_wide %>%
  rename(
    `Latitudinal shift` = lat_shift,
    `Longitudinal shift` = lon_shift,
    `Depth shift` = depth_shift,
    `Thermal-niche shift` = thermal_shift
  ) %>% mutate(
    Species = gsub("_", " ", Species),
    Species = paste0("\\textit{", Species, "}")
  )


kbl(supp_table_wide, booktabs = T, "latex",escape=FALSE,longtable = T,
    linesep = "") %>% kable_styling(latex_options = c("repeat_header","striped"), font_size = 7) %>% writeLines( here('output/tables/supp/species_slopes.tex'))

# save as rds and csv
write_rds(species_summary, here('R/data/processed/bayesian_species_trends.rds'))
write_csv(species_summary, here('output/supp_data/bayesian_species_trends.csv'))


