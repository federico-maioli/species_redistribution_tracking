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

# first make the big map --------------------------------------------------
# Get coastline data from Natural Earth
world <- ne_countries(scale = "medium", returnclass = "sf")

data <- read_rds(here('R/data/processed/fishglob_clean.rds'))


# get unique hauls
unique_hauls <-  data %>%
  select(survey,region,region_short,year,latitude,longitude,haul_id) %>%
  distinct() %>% filter(!(survey == "SCS" & latitude < 38))

unique_hauls_sf <- st_as_sf(unique_hauls, coords = c("longitude", "latitude"), crs = 4326)

# create survey hull
survey_hulls <- unique_hauls_sf |>
  st_cast("POINT") |>   # ensure geometries are POINTS
  group_by(survey,region) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  mutate(geometry = st_concave_hull(geometry, ratio = 0.1)) |>  # ratio ∈ [0,1], lower = more concave
  st_make_valid()

survey_hulls <- survey_hulls %>% mutate(region_short = case_when(
  region == "Gulf of Alaska"        ~ "GOA",
  region == "British Columbia"       ~ "BC",
  region == "U.S. West Coast"       ~ "USWC",
  region == "Baltic Sea"            ~ "BAL",
  region == "North Sea"             ~ "NS",
  region == "Barents Sea"           ~ "BS",
  region == "East Bering Sea"       ~ "EBS",
  region == "Gulf of Mexico"        ~ "GOM",
  region == "Celtic-Biscay Shelf"   ~ "CBS",
  region == "North Iberian Coast"   ~ "NIC",
  region == "Northeast U.S. and Scotian Shelf" ~ "NEUS-SS",
  TRUE ~ NA_character_
))

survey_regions <- survey_hulls %>%
  group_by(region_short, region) %>%   # keep region name for labeling
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()  # ensure valid geometries

# Example Lambert Conformal Conic projection
lcc_crs <- "+proj=lcc +lat_1=30 +lat_2=60 +lat_0=40 +lon_0=-60 +datum=WGS84 +units=km"

base_colors <- paletteer_d("tidyquant::tq_dark", n = 12) |> as.character()

# assign them explicitly to regions
region_colors <- c(
  "GOA"     = base_colors[1],  # Gulf of Alaska
  "BC"      = base_colors[2],  # British Columbia
  "USWC"    = base_colors[3],  # U.S. West Coast
  "BAL"     = base_colors[4],  # Baltic Sea
  "NS"      = base_colors[5],  # North Sea
  "BS"      = base_colors[6],  # Barents Sea
  "EBS"     = base_colors[7],  # East Bering Sea
  "GOM"     = base_colors[8],  # Gulf of Mexico
  "CBS"     = base_colors[9],  # Celtic-Biscay Shelf
  "NIC"     = base_colors[10], # North Iberian Coast
  "NEUS-SS" = base_colors[11]  # Northeast U.S. & Scotian Shelf
)

map <-  ggplot() +
  geom_sf(data = survey_regions, aes(fill = region_short), alpha = 0.5) +
  geom_sf(data = world, color = 'grey20', fill = 'antiquewhite', size = .1) +
  geom_text_repel(data = survey_regions, 
                  aes(label = region_short, geometry = geometry),
                  size = 2,
                  stat = "sf_coordinates", min.segment.length = 0.5, seed = 1, box.padding = 0.3, color = "black", bg.color = "white", # shadow color
                  bg.r = 0.15 ) +          # shadow radius) + 
  scale_fill_manual(values = region_colors, name = "Region") + coord_sf(xlim = c(-180,60),ylim = c(20, 85)) + theme_light(base_size=10) +
  coord_sf(crs = lcc_crs, xlim = c(-5000, 5000), ylim = c(-500, 6000)) + theme(legend.position = 'none' , plot.margin = margin(0, 0, 0, 0) )  +
  labs(x = NULL, y = NULL)


# now make the regional temp trend ----------------------------------------

grid <- read_rds(here('R/data/processed/prediction_grid.rds'))


grid <- grid %>% mutate(region_full = case_when(
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

temp_summary <- grid %>%
  group_by(year, region_full, region_short) %>%
  summarise(mean_temp = mean(mean_temp, na.rm = TRUE), .groups = "drop")

slopes <- temp_summary %>%
  group_by(region_full,region_short) %>%
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
  left_join(slopes %>% select(region_full, region_short, slope_decade, slope_label,slope_expr),
            by = c("region_full", 'region_short'))

regions <- unique(temp_summary$region_short)

region_plots <- map(regions, function(r) {
  ggplot(temp_summary %>% filter(region_short == r),
         aes(x = year, y = mean_temp)) +
    geom_line(color = region_colors[r], alpha = .4, linewidth = .3) +
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
    ) + scale_y_continuous(breaks = scales::breaks_extended(n = 4),
                           labels = function(x) parse(text = paste0(x, "*degree*C"))
    ) +
    theme_light(base_size = 5) +
    labs(
      title = unique(temp_summary$region_full[temp_summary$region_short == r]),
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

names(region_plots) <- regions

wrap_plots(region_plots)

# now haul plots
haul_counts <- unique_hauls %>%
  group_by(region_short, year) %>%
  summarise(n_hauls = n(), .groups = "drop") %>%  mutate(region_short = factor(
    region_short,
    levels = rev(c(
      "EBS",
      "GOA",
      "BC",
      "USWC",
      "NEUS-SS",
      "GOM",
      "BS",
      "NS",
      "BAL",
      "CBS",
      "NIC"
    ))
  ))

haul_plot <- ggplot(haul_counts, aes(x = year, y = region_short, fill = n_hauls)) +
  geom_tile(color = "black") +
  scale_fill_viridis_c(option = "mako", name = "Number of hauls" , guide = guide_colorbar(
    title.position = "left",   # put title above the keys
    title.vjust = 1,  
    barwidth = 15,            # width of the legend keys
    barheight = 0.3           # height of the legend keys
  )) +
  labs(x = "Year", y = NULL) +
  scale_x_continuous(expand = c(0,0)) + 
  #facet_grid(region ~ ., scales = "free_y", space = "free", switch = "x") +
  theme_minimal(base_size = 7) +
  theme(
    #strip.placement = "outside",                   # put strip outside axis
    #strip.text.y.right = element_text(angle = 0,hjust = 0),   # readable horizontal labels
    panel.spacing = unit(0.4, "lines"),             # tighter panels
    legend.position = "bottom",                # legend at the bottom
    #    legend.key.width = unit(2, "cm"),   # stretch horizontally
    #    legend.key.height = unit(0.15, "cm"),  # make thinner
    legend.direction = "horizontal" ,
    panel.grid.major = element_blank(),  # remove major grid lines
    panel.grid.minor = element_blank(),   # remove minor grid lines
    legend.box.margin = margin(-5, 0, 0, 0),   # pull legend closer to plot (top, right, bottom, left)
    legend.margin = margin(0, 0, 0, 0),        # remove inner legend padding
    plot.margin = margin(0, 5, 0, 5)    
  )

layout <- "
AABBCCDDEE
FFGGGGGGHH
IIGGGGGGLL
MMNNNNNNOO
"
final_plot <- (
  region_plots[['EBS']] + region_plots[['NEUS-SS']] + region_plots[['GOM']] + region_plots[['BS']]  + region_plots[['NS']] # top
  + region_plots[['GOA']] + free(map) +
    region_plots[['BAL']] + region_plots[['BC']] + region_plots[['CBS']] + 
    region_plots[["USWC"]] + 
    haul_plot + region_plots[["NIC"]]  #+ plot_spacer()                  # bottom
  # your central map
) + plot_layout(
  design = layout,
  #guides = "collect",
  #axes = "collect"
)

final_plot

ggsave(
  filename = "final_patchwork.png",
  width = 180,
  height = 150,
  dpi = 600,
  units = "mm"
)
