
# libraries ---------------------------------------------------------------

library(tidyverse)
library(here)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggrepel)
library(paletteer)
library(ggpubr)

# Get coastline data from Natural Earth
world <- ne_countries(scale = "medium", returnclass = "sf")

# load data ---------------------------------------------------------------

data <- read_rds(here('data/processed/fishglob_clean.rds'))

# no species region
# no_species <- data %>% group_by(region) %>%
#   summarise(`Number of species` = n_distinct(accepted_name), .groups = "drop") %>% rename(Region = region) 

# get the table
#table <-  ggtexttable(no_species, cols = c('Region', 'Number of species'), rows = NULL, theme = ttheme("mBlue"))

data <- data %>%
  left_join(dfo_surveys, by = "haul_id", suffix = c("", "_new")) %>%
  mutate(
    survey = ifelse(!is.na(survey_new), survey_new, survey)
  ) %>%
  select(-survey_new)  # remove the temporary column


unique_hauls <-  data %>%
  select(survey,region,region_short,year,latitude,longitude,haul_id) %>%
  distinct() %>% filter(!(survey == "SCS" & latitude < 38))

unique_hauls_sf <- st_as_sf(unique_hauls, coords = c("longitude", "latitude"), crs = 4326)

region_colors <- paletteer_d("tidyquant::tq_dark", n = 12) |> as.character()

# create survey hull
survey_hulls <- unique_hauls_sf |>
  st_cast("POINT") |>   # ensure geometries are POINTS
  group_by(survey,region) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  mutate(geometry = st_concave_hull(geometry, ratio = 0.1)) |>  # ratio ∈ [0,1], lower = more concave
  st_make_valid()

# Example Lambert Conformal Conic projection
lcc_crs <- "+proj=lcc +lat_1=30 +lat_2=60 +lat_0=40 +lon_0=-60 +datum=WGS84 +units=km"


map_plot <-  ggplot() +
  geom_sf(data = survey_hulls, aes(fill = region), alpha = 0.5) +
  geom_sf(data = world, color = 'grey20', fill = 'antiquewhite', size = .1) +
  geom_text_repel(data = survey_hulls, 
                  aes(label = survey, geometry = geometry),
                  stat = "sf_coordinates", min.segment.length = 0, seed = 1, box.padding = 0.9) + 
  scale_fill_manual(values = region_colors, name = "Region") + coord_sf(xlim = c(-180,60),ylim = c(20, 85)) +
  coord_sf(crs = lcc_crs, xlim = c(-5000, 5000), ylim = c(-500, 6000)) + theme_light()  +
  labs(x = "Longitude", y = "Latitude")

haul_counts <- unique_hauls %>%
  group_by(region, year, survey) %>%
  summarise(n_hauls = n(), .groups = "drop")


heatmap_plot <- ggplot(haul_counts, aes(x = year, y = survey, fill = n_hauls)) +
  geom_tile(color = "black") +
  scale_fill_viridis_c(option = "mako", name = "Number of hauls" , guide = guide_colorbar(
    title.position = "top",   # put title above the keys
    title.hjust = 0.5,       # center the title
    barwidth = 15,            # width of the legend keys
    barheight = 0.3           # height of the legend keys
  )) +
  labs(x = "Year", y = NULL) +
  scale_x_continuous(expand = c(0,0)) + 
  facet_grid(region ~ ., scales = "free_y", space = "free", switch = "x") +
  theme_minimal(base_size = 7) +
  theme(
    strip.placement = "outside",                   # put strip outside axis
    strip.text.y.right = element_text(angle = 0,hjust = 0),   # readable horizontal labels
    panel.spacing = unit(0.4, "lines"),             # tighter panels
    legend.position = "bottom",                # legend at the bottom
#    legend.key.width = unit(2, "cm"),   # stretch horizontally
#    legend.key.height = unit(0.15, "cm"),  # make thinner
    legend.direction = "horizontal" ,
    panel.grid.major = element_blank(),  # remove major grid lines
    panel.grid.minor = element_blank()   # remove minor grid lines
  )


 

plot_grid(map_plot, heatmap_plot, rel_heights = c(1.8, 1), # map taller, heatmap shorter
          align = "v", 
          axis = "r",
          ncol = 1, labels = c('a)', 'b)'), label_size = 12)

# assemble plot
#top <- plot_grid(map_plot, table, rel_widths = c(1,.2), rel_heights = c(1, .3))



ggsave(here('results/figures/map.png'), bg = 'white', width = 10, height = 8, dpi = 600)



# ggdraw(map_plot) +
#   draw_plot(heatmap_plot, x = 0.57, y = .05, width = 0.40, height = 0.32)

#+
  # draw_plot_label(
  #   c("A", "B"),
  #   c(0, 0.45),
  #   c(1, 0.95),
  #   size = 12
  # )



