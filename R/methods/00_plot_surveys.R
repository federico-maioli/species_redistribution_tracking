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
library(broom)
# This works only after the prediction grid has been created

# Get coastline data from Natural Earth
world <- ne_countries(scale = "medium", returnclass = "sf")

# load data ---------------------------------------------------------------

data <- read_rds(here('R/data/processed/fishglob_clean.rds'))


# get unique hauls
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
                  stat = "sf_coordinates", min.segment.length = 0, seed = 1, box.padding = 0.9, color = "black", bg.color = "white", # shadow color
                  bg.r = 0.15 ) +          # shadow radius) + 
  scale_fill_manual(values = region_colors, name = "Region") + coord_sf(xlim = c(-180,60),ylim = c(20, 85)) + theme_light(base_size=10) +
  coord_sf(crs = lcc_crs, xlim = c(-5000, 5000), ylim = c(-500, 6000)) + theme(legend.position = 'none' , plot.margin = margin(0, 0, 0, 0) )  +
  labs(x = NULL, y = NULL)

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

# get grid data
grid <- read_rds(here('R/data/processed/prediction_grid.rds'))

ggplot(grid %>% group_by(year,region_short) %>% summarise(mean_temp = mean(mean_temp)), aes(x = year, y = mean_temp, color = region_short)) +
  geom_line(alpha = 0.6) +                          # temperature line
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +  # linear trend
  facet_wrap(~ region_short, scales = "free_y") +   # one panel per region
  theme_minimal(base_size = 13) +
  labs(
    x = "Date",
    y = "Mean Temperature (°C)",
    title = "Temperature Trends by Region",
    subtitle = "With linear trend lines per region"
  )


temp_summary <- grid %>%
  group_by(year, region_short) %>%
  summarise(mean_temp = mean(mean_temp, na.rm = TRUE), .groups = "drop")

# ensure region order matches your palette
regions <- unique(temp_summary$region_short)

slopes <- temp_summary %>%
  group_by(region_short) %>%
  summarise(
    slope_decade = coef(lm(mean_temp ~ year))[2] * 10,
    .groups = "drop"
  ) %>%
  mutate(
    slope_label = sprintf("%+.2f", slope_decade),
    # This builds a character string like "+0.43*degree*C~decade^-1"
    slope_expr = paste0("bold('", slope_label, "'*degree*C~decade^'-1')")
  )

temp_summary <- temp_summary %>%
  left_join(slopes %>% select(region_short, slope_decade, slope_label,slope_expr),
            by = "region_short")

region_plots <- map(regions, function(r) {
  ggplot(temp_summary %>% filter(region_short == r),
         aes(x = year, y = mean_temp)) +
    geom_line(color = region_colors[match(r, regions)], alpha = .3, linewidth = .3) +
    #geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
    geom_textsmooth(
      aes(label = slope_expr),
      method = "lm",
      se = FALSE,
      size = 2.2,
      color = "black",
      parse = TRUE,
      vjust = - .3,
      linetype = 3,
      hjust = 1,
      linewidth = 0.3
    ) + scale_y_continuous(breaks = scales::breaks_pretty(n = 4),
      labels = function(x) parse(text = paste0(x, "*degree*C"))
    ) +
    theme_light(base_size = 5) +
    labs(
      title = r,
      x = "Year",
      y = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 6),
      axis.ticks = element_line(color = "black"),   # show ticks
      axis.ticks.length = unit(2, "pt"),
      axis.text = element_text(size = 6),
      #axis.title.y = element_text(margin = margin(r = 5)),
      #plot.margin = margin(5, 5, 5, 5),
      panel.grid.major = element_blank(),  # remove major grid lines
      panel.grid.minor = element_blank()
    )
})

region_plots[[2]]

layout <- "
AABBCCDDJJ
EEFFFFFFGG
HHFFFFFFII
LLFFFFFFOO
"
final_plot <- (
  region_plots[[1]] + region_plots[[2]] + region_plots[[3]] +  # top-left group
    region_plots[[4]]  + region_plots[[5]]+ free(map_plot) + region_plots[[6]] +  # right side
    region_plots[[7]] + region_plots[[8]] + region_plots[[9]] +  # left side
    region_plots[[10]] + region_plots[[11]]  + plot_spacer()                  # bottom
                       # your central map
) + plot_layout(design = layout, 
                guides = "collect", 
                axes = "collect"
                )

final_plot

ggsave(
  filename = "final_patchwork.png",
  width = 180,
  height = 120,
  dpi = 600,
  units = "mm"
)

free(final_plot) / heatmap_plot

ggsave(
  filename = "final_patchwork.png",
  width = 180,
  height = 190,
  dpi = 600,
  units = "mm"
)

# ggdraw(map_plot) +
#   draw_plot(heatmap_plot, x = 0.57, y = .05, width = 0.40, height = 0.32)

#+
  # draw_plot_label(
  #   c("A", "B"),
  #   c(0, 0.45),
  #   c(1, 0.95),
  #   size = 12
  # )
layout <- "
AAAABBBB
AAAABBBB
LLLFFFFR
LLLFFFFR
LLLFFFFR
"


