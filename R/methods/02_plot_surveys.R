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
  st_make_valid() %>% filter(region_short %in% c('BC','CBS','NEUS-SS'))

# combine hulls by region
# survey_regions <- survey_hulls %>%
#   group_by(region_short, region_full) %>%
#   summarise(geometry = st_union(geometry), .groups = "drop") %>%
#   st_make_valid()

# define Lambert Conformal Conic projection
lcc_crs <- "+proj=lcc +lat_1=30 +lat_2=60 +lat_0=40 +lon_0=-60 +datum=WGS84 +units=km"

# assign colors to regions
base_colors <- paletteer_d("tidyquant::tq_dark", n = 12) |> as.character()
region_colors <- c(
  #"GOA"     = base_colors[1],
  "British Columbia"      = base_colors[2],
  #"USWC"    = base_colors[3],
  #"BAL"     = base_colors[4],
  #"NS"      = base_colors[5],
  #"BS"      = base_colors[6],
  #"EBS"     = base_colors[7],
  #"GOM"     = base_colors[8],
  "Celtic-Biscay Shelf"     = base_colors[9],
  #"NIC"     = base_colors[10],
  "Northeast U.S. & Scotian Shelf" = base_colors[11]
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
  geom_sf(data = survey_hulls, aes(fill = region_full), alpha = 0.3) +
  geom_sf(data = world, color = 'black', fill = 'grey80', size = 0.1) +
  
  # Region labels
  geom_text_repel(
    data = survey_hulls,
    aes(label = survey, geometry = geometry),
    size = 2.8,
    stat = "sf_coordinates",
    min.segment.length = 0.5,
    seed = 1,
    box.padding = 0.9,
    color = "black",
    bg.color = "white",
    bg.r = 0.2
  ) +
  geom_sf_text(
    data = axis_labels_sf,
    aes(label = label),
    size = 1.5,
    color = "grey60"
  ) +
  
  scale_fill_manual(values = region_colors, name = "Region") +
  coord_sf(crs = lcc_crs, xlim = c(-5100, 5100), ylim = c(-600, 6100)) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = 'bottom',
   # plot.margin = margin(0, 0, 0, 0),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
  ) +
  labs(x = NULL, y = NULL)

map


# create haul counts heatmap ---------------------------------------------

haul_counts <- data %>%
  select(survey, region_full, region_short, year, month, latitude, longitude, haul_id) %>%
  distinct() %>%
  filter(!(survey == "SCS" & latitude < 38)) %>%
  group_by(region_short,region_full,survey, year, month) %>%
  summarise(n_hauls = n(), .groups = "drop") %>%
  mutate(region_short = factor(
    region_short,
    levels = rev(c("EBS","GOA","BC","USWC","NEUS-SS","GOM","BS","NS","BAL","CBS","NIC"))
  ))

haul_counts2 <- haul_counts |> filter(region_short %in% c('BC','CBS','NEUS-SS')) %>% 
  dplyr::group_by(region_full) |>
  dplyr::mutate(survey = droplevels(survey)) |>
  dplyr::ungroup()


haul_plot <- ggplot(haul_counts2,
       aes(x = month,
           y = survey,
           fill = n_hauls)) +
  geom_tile(color = "black", alpha = 0.8) +
  scale_x_continuous(breaks = 1:12,expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Number of hauls",
                       guide = guide_colorbar(
                         title.position = "left",
                         title.vjust = 1,
                         barwidth = 15,
                         barheight = 0.3)) +
  labs(x = "Month", y = "Survey") +
  facet_wrap(~ region_full, scales = 'free_y',space = "free_y") +
  theme_bw(base_size = 8) +
  theme(
    strip.background = element_blank(),  # remove rectangle behind facet label
    strip.text = element_text(face = "bold"), # keep facet label text bold
    panel.spacing.y = unit(1.2, "mm"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.box.margin = margin(-5, 0, 0, 0),
    legend.margin = margin(2, 0, 0, 0),
   # plot.margin = margin(2, 8, 2, 5)
  )

map / haul_plot & plot_layout(heights = c(3,1)) &
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 12,face='bold'))

#

ggsave(
  filename = "output/figures/supp/surveys_supp.png",
  width = 180,
  height = 140,
  dpi = 600,
  units = "mm"
)


# assemble final plot ----------------------------------------------------

haul_counts_complete <- haul_counts2 %>%
  group_by(region_full, survey) %>%
  complete(month = 1:12, fill = list(n_hauls = NA)) %>%
  ungroup()

# plot
haul_plot <- ggplot(haul_counts_complete,
                    aes(x = month,
                        y = survey,
                        fill = n_hauls)) +
  geom_tile(color = "black", alpha = 0.8) +
  scale_x_continuous(breaks = 1:12, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Number of hauls",
                       guide = guide_colorbar(
                         title.position = "left",
                         title.vjust = 1,
                         barwidth = 15,
                         barheight = 0.3)) +
  labs(x = "Month", y = "Survey") +
  facet_wrap(~ region_full, scales = 'free_y', space = "free_y") +
  theme_bw(base_size = 8) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.4, "lines"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.box.margin = margin(-5, 0, 0, 0),
    legend.margin = margin(2, 0, 0, 0),
    plot.margin = margin(2, 5, 2, 5)
  )

