# libraries ---------------------------------------------------------------

library(tidyverse)
library(here)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggrepel)
library(paletteer)
library(ggpubr)
library(patchwork)
library(geomtextpath)


# load and prepare data ---------------------------------------------------

# read processed fishglob data
data <- read_rds(here('R/data/processed/fishglob_clean.rds'))

# extract unique hauls, filter out invalid SCS points
unique_hauls <- data %>%
  select(survey, region_full, region_short, year, latitude, longitude, haul_id) %>%
  distinct() %>%
  filter(!(survey == "SCS" & latitude < 38))

# convert to sf object
unique_hauls_sf <- st_as_sf(unique_hauls, coords = c("longitude", "latitude"), crs = 4326)


# create survey hulls -----------------------------------------------------

survey_hulls <- unique_hauls_sf %>%
  st_cast("POINT") %>%
  group_by(survey, region_full, region_short) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  mutate(geometry = st_concave_hull(geometry, ratio = 0.1)) %>%
  st_make_valid()

# combine hulls by region
survey_regions <- survey_hulls %>%
  group_by(region_short, region_full) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

# define Lambert Conformal Conic projection
lcc_crs <- "+proj=lcc +lat_1=30 +lat_2=60 +lat_0=40 +lon_0=-60 +datum=WGS84 +units=km"

# assign colors to regions
base_colors <- paletteer_d("tidyquant::tq_dark", n = 12) |> as.character()
region_colors <- c(
  "GOA"     = base_colors[1],
  "BC"      = base_colors[2],
  "USWC"    = base_colors[3],
  "BAL"     = base_colors[4],
  "NS"      = base_colors[5],
  "BS"      = base_colors[6],
  "EBS"     = base_colors[7],
  "GOM"     = base_colors[8],
  "CBS"     = base_colors[9],
  "NIC"     = base_colors[10],
  "NEUS-SS" = base_colors[11]
)


# create big map ----------------------------------------------------------

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_union() %>%              # merge all countries into one geometry
  st_sf(geometry = .)  

# manually set labels 
axis_labels <- rbind(
  data.frame(
    lon   = rep(-60, 3),
    lat   = c(40, 60, 70),
    label = paste0(c(40, 60, 70), "\u00B0N")
  ),
  data.frame(
    lon = c(-120, -60, 0),
    lat = rep(85, 3),
    label = c(
      paste0(abs(-120), "\u00B0W"),
      paste0(abs(-60), "\u00B0W"),
      "0\u00B0"  # no E/W for 0
    )
  )
)

# Convert to sf and project
axis_labels_sf <- st_as_sf(axis_labels, coords = c("lon","lat"), crs = 4326) |>
  st_transform(lcc_crs)

map <- ggplot() +
  geom_sf(data = survey_regions, aes(fill = region_short), alpha = 0.7) +
  geom_sf(data = world, color = 'grey40', fill = 'grey90', size = 0.1) +
  
  # Region labels
  geom_text_repel(
    data = survey_regions,
    aes(label = region_short, geometry = geometry),
    size = 3.4,
    stat = "sf_coordinates",
    min.segment.length = 0.5,
    seed = 1,
    box.padding = 0.3,
    color = "black",
    bg.color = "white",
    bg.r = 0.1
  ) +
  geom_sf_text(
    data = axis_labels_sf,
    aes(label = label),
    size = 1.5,
    color = "grey50"
  ) +
  
  scale_fill_manual(values = region_colors, name = "Region") +
  coord_sf(crs = lcc_crs, xlim = c(-5100, 5100), ylim = c(-600, 6100)) +
  theme_light(base_size = 10) +
  theme(
    legend.position = 'none',
    plot.margin = margin(0, 0, 0, 0),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid.major = element_line(linewidth = 0.05, color = "grey85"),
    panel.grid.minor = element_line(linewidth = 0.05, color = "grey85"),
  ) +
  labs(x = NULL, y = NULL)

map
# prepare temperature trend data -----------------------------------------

grid <- read_rds(here('R/data/processed/prediction_grid.rds'))

grid <- grid %>%
  mutate(region_full = case_when(
    region_short == "GOA"      ~ "Gulf of Alaska (GOA)",
    region_short == "BC"       ~ "British Columbia (BC)",
    region_short == "USWC"     ~ "U.S. West Coast (USWC)",
    region_short == "BAL"      ~ "Baltic Sea (BAL)",
    region_short == "NS"       ~ "North Sea (NS)",
    region_short == "BS"       ~ "Barents Sea (BS)",
    region_short == "EBS"      ~ "East Bering Sea (EBS)",
    region_short == "GOM"      ~ "Gulf of Mexico (GOM)",
    region_short == "CBS"      ~ "Celtic-Biscay Shelf (CBS)",
    region_short == "NIC"      ~ "North Iberian Coast (NIC)",
    region_short == "NEUS-SS"  ~ "Northeast U.S. and\nScotian Shelf (NEUS-SS)",
    TRUE ~ NA_character_
  ))

# summarize mean temperature per year-region
temp_summary <- grid %>%
  group_by(year, region_full, region_short) %>%
  summarise(mean_temp = mean(mean_year_temp, na.rm = TRUE), .groups = "drop")

temp_summary <- temp_summary %>%
  filter(
    (region_short == "NIC"      & year >= 1994 & year <= 2019) |
      (region_short == "CBS"      & year >= 1994 & year <= 2020) |
      (region_short == "BAL"      & year >= 1999 & year <= 2020) |
      (region_short == "NS"       & year >= 1994 & year <= 2020) |
      (region_short == "BS"       & year >= 2004 & year <= 2021) |
      (region_short == "GOM"      & year >= 1994 & year <= 2020) |
      (region_short == "NEUS-SS"  & year >= 1994 & year <= 2020) |
      (region_short == "USWC"     & year >= 2003 & year <= 2018) |
      (region_short == "BC"       & year >= 2003 & year <= 2019) |
      (region_short == "GOA"      & year >= 1996 & year <= 2019) |
      (region_short == "EBS"      & year >= 1994 & year <= 2019)
  )

# compute slopes per region
slopes <- temp_summary %>%
  group_by(region_full, region_short) %>%
  summarise(
    slope_decade = coef(lm(mean_temp ~ year))[2] * 10,
    .groups = "drop"
  ) %>%
  mutate(
    slope_label = sprintf("%+.2f", slope_decade),
    slope_expr = paste0("bold('", slope_label, "'*degree*C~decade^'-1')")
  )

# join slope info back
temp_summary <- temp_summary %>%
  left_join(slopes %>% select(region_full, region_short, slope_decade, slope_label, slope_expr),
            by = c("region_full", 'region_short'))


# create regional temperature plots --------------------------------------

regions <- unique(temp_summary$region_short)

region_plots <- map(regions, function(r) {
  data_r <- temp_summary %>% filter(region_short == r)
  y_breaks <- scales::breaks_extended(n = 5)(range(data_r$mean_temp, na.rm = TRUE))
  
  ggplot(data_r, aes(x = year, y = mean_temp)) +
    geom_line(color = region_colors[r], alpha = .9, linewidth = .4) +
    annotate(
      "text",
      x = Inf, y = -Inf,
      label = unique(data_r$slope_expr),
      hjust = 1.05, vjust = -0.4,
      size = 2.2,
      colour = "grey30",
      parse = TRUE
    ) +
    # geom_textsmooth(
    #   aes(label = slope_expr),
    #   method = "lm",
    #   se = FALSE,
    #   size = 2.2,
    #   color = "black",
    #   parse = TRUE,
    #   vjust = -0.3,
    #   linetype = 3,
    #   hjust = 1,
    #   linewidth = 0.3
    # ) +
    scale_y_continuous(
      breaks = y_breaks,
        expand = expansion(mult = c(0.15, 0.02)),
      labels = sapply(y_breaks, function(y) as.expression(bquote(.(y)*degree*C)))
    ) +
    theme_bw(base_size = 7) +
    labs(
      title = unique(data_r$region_full),
      x = NULL,
      y = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 6.5),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(2, "pt"),
      axis.text = element_text(size = 6.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x  = element_text(angle = 45, hjust = 1),
     # theme(aspect.ratio = 0.4) 
     plot.margin = margin(2,2,5,2)
    )
})

names(region_plots) <- regions


# create haul counts heatmap ---------------------------------------------

haul_counts <- unique_hauls %>%
  group_by(region_short, year) %>%
  summarise(n_hauls = n(), .groups = "drop") %>%
  mutate(region_short = factor(
    region_short,
    levels = rev(c("EBS","GOA","BC","USWC","NEUS-SS","GOM","BS","NS","BAL","CBS","NIC"))
  ))

haul_plot <- ggplot(haul_counts, aes(x = year, y = region_short, fill = n_hauls)) +
  geom_tile(color = "black", alpha = 0.8) +
  scale_fill_gradient(
    low  = "#f1f5f9",  # very light blue-grey
    high = "#5b84b1",  # muted steel blue
    name = "Number of hauls",
    guide = guide_colorbar(
      title.position = "left",
      title.vjust = 1,
      barwidth = 15,
      barheight = 0.3
    )
  ) +
  # scale_fill_viridis_c(option = "mako", direction = -1, name = "Number of hauls",
  #                      guide = guide_colorbar(
  #                        title.position = "left",
  #                        title.vjust = 1,
  #                        barwidth = 15,
  #                        barheight = 0.3
  #                      )) +
  labs(x = NULL, y = NULL,title = "Temporal coverage and sampling intensity") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_minimal(base_size = 8) +
  theme(
    panel.spacing = unit(1.2, "lines"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.box.margin = margin(-5, 0, 0, 0),
    legend.margin = margin(1, 0, 0, 0),
    plot.margin = margin(0, 5, 0, 5),
    plot.title = element_text(
      size = 7.5,
      hjust = 0.5,
      face = "bold",
      colour = "black",
      margin = margin(b = 4)
    ),
    plot.title.position = "plot"
)


# assemble final plot ----------------------------------------------------

layout <- "
AABBCCDDEE
FFGGGGGGHH
IIGGGGGGLL
MMNNNNNNOO
"


final_plot <- (
  free(region_plots[['EBS']]) + 
    free(region_plots[['NEUS-SS']]) + 
    free(region_plots[['GOM']]) + 
    free(region_plots[['BS']]) + 
    free(region_plots[['NS']]) +
    free(region_plots[['GOA']]) + 
    free(map) + 
    free(region_plots[['BAL']]) + 
    free(region_plots[['BC']]) + 
    free(region_plots[['CBS']]) +
    free(region_plots[['USWC']]) +
    free(haul_plot) + 
    free(region_plots[['NIC']])
) + plot_layout(design = layout)

# 
# final_plot <- (
#   region_plots[['EBS']] + region_plots[['NEUS-SS']] + region_plots[['GOM']] + region_plots[['BS']] + region_plots[['NS']] +
#     region_plots[['GOA']] + map + region_plots[['BAL']] + region_plots[['BC']] + region_plots[['CBS']] +
#     region_plots[['USWC']] +
#     haul_plot + region_plots[['NIC']]
# ) + plot_layout(design = layout) 

final_plot

# save final plot
ggsave(
  filename = "output/figures/main/map.png",
  width = 180,
  height = 150,
  dpi = 600,
  units = "mm"
)

