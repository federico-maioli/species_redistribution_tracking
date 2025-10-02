# libraries ---------------------------------------------------------------
library(tidyverse)
library(here)
library(sdmTMB)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(terra)
library(ncdf4)

# load function and data -------------------------------------------------

source(here('R/utils/exclude_too_far.R'))
data <- readRDS(here('R/data/processed/fishglob_clean.rds'))

# prepare survey spatial grid ----------------------------------------------

unique_locations <- data |> distinct(region_short, depth, latitude, longitude)
#unique_locations <- data |> distinct(region_short, depth, latitude, longitude, year)
split_data <- split(unique_locations, unique_locations$region_short)

# grid settings
grids <- list()
cell_width <- 4  # grid resolution (4 km)

round_down_even <- function(x, base = 1000) base * floor(x / base)

# load coast and depth data
coastline <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf") |> 
  st_transform(crs = 4326) |> 
  st_make_valid()

depth_raster <- terra::rast(here('R/data/environmental/depth/GEBCO_2023_sub_ice_topo.nc'))

depth_raster[depth_raster > 0] <- NA # remove positive values

depth_raster <- terra::focal(
  depth_raster,
  w = 5,                 # 5x5 neighborhood
  fun = mean,             # average nearby cell values
  na.rm = TRUE,
  na.policy = "only"      # ONLY fill NA cells
)

names(depth_raster) <- 'elevation' # focal changes name

# process regions ----------------------------------------------------------

grids = list()

for (this_region in names(split_data)) {
  
  sub <- split_data[[this_region]]
  
  sub <- sdmTMB::add_utm_columns(sub, ll_names = c("longitude", "latitude"), units = 'm')
  crs <- sdmTMB::get_crs(sub, c("longitude", "latitude"))
  sub <- sub |> distinct(X, Y, .keep_all = TRUE)
  
  # create concave hull
  sf_use_s2(FALSE)
  sf_points <- st_as_sf(sub, coords = c("X", "Y"), crs = crs)
  boundary <- st_buffer(st_concave_hull(st_union(sf_points), ratio = 0.05), dist = 10000)
  
  # generate regular grid
  pred_grid <- expand.grid(
    X = seq(round_down_even(min(sub$X)), max(sub$X), cell_width * 1000),
    Y = seq(round_down_even(min(sub$Y)), max(sub$Y), cell_width * 1000)
  )
  
  # keep only points inside boundary
  sf_object <- st_as_sf(pred_grid, coords = c("X", "Y"), crs = crs)
  inside_boundary <- st_intersects(sf_object, boundary, sparse = FALSE)
  pred_grid <- pred_grid[apply(inside_boundary, 1, any), ]
  
  # transform to lat/lon and remove land points
  sf_use_s2(TRUE)
  sf_wgs84 <- st_transform(st_as_sf(pred_grid, coords = c("X", "Y"), crs = crs), crs = 4326)
  on_land <- st_intersects(sf_wgs84, coastline, sparse = FALSE)
  pred_grid$longitude <- st_coordinates(sf_wgs84)[, 1]
  pred_grid$latitude  <- st_coordinates(sf_wgs84)[, 2]
  pred_grid <- pred_grid[!apply(on_land, 1, any), ]
  
  # apply depth filter
  pred_grid$depth <- terra::extract(depth_raster$elevation, 
                                    pred_grid |> select(longitude, latitude),method = 'simple')$elevation
  pred_grid <- pred_grid |> filter(depth <= 0)
  pred_grid$depth <- abs(pred_grid$depth)
  
  depth_percentiles <- quantile(sub$depth, probs = c(0.01, 0.99), na.rm = TRUE)
  pred_grid <- pred_grid |>
    filter(depth >= depth_percentiles[1], depth <= depth_percentiles[2])
  
  # exclude points too far from observations
  dist_limit <- if (this_region == 'BS') 40000 else 30000
  too_far <- exclude_too_far(pred_grid$X, pred_grid$Y, sub$X, sub$Y, dist = dist_limit)
  pred_grid <- pred_grid[!too_far, ]
  
  # store grid (convert X/Y to km)
  grids[[this_region]] <- pred_grid |> mutate(X = X / 1000, Y = Y / 1000)
}

grid <-  grids %>% bind_rows(.id = "region_short") %>%
  filter(region_short %in% unique(data$region_short)) %>%
  left_join(distinct(data, region_short, min_month), by = "region_short")

grid <- replicate_df(grid, time_name = "year", time_values = 1994:2021) %>%
  mutate(
    month_year_date = as.Date(paste0(year, "-", sprintf("%02d", min_month), "-01"))
  )


# match temperature data --------------------------------------------------

nc_path <- here('R/data/environmental/temperature/cmems_mod_glo_phy_my_0.083deg_P1M-m_bottomT_180.00W-179.92E_20.00N-87.00N_1993-01-01-2021-06-01.nc')

# source time to date function for .nc files
source(here('R/utils/nc_time_to_date.R'))

# Extract time dimension
nc_temp <- nc_open(nc_path)
dates <- nc_time_to_date(nc_temp$dim$time$vals, nc_temp$dim$time$units)
nc_close(nc_temp)

# Load rasters
temp_rast <- rast(nc_path)

names(temp_rast) <- as.character(dates)

temp_rast <- terra::focal(
  temp_rast,
  w = 7,                 # 5x5 neighborhood
  fun = mean,             # average nearby cell values
  na.rm = TRUE,
  na.policy = "only"      # ONLY fill NA cells
)

temp_list <- list()

unique_dates <- sort(unique(grid$month_year_date))

for (i in seq_along(unique_dates)) {
  
  this_date <- unique_dates[i]
  
  data_date <- filter(grid, month_year_date == this_date)
  
  past_12_months <- seq(from = this_date %m-% months(12), to = this_date %m-% months(1), by = "1 month")
  past_12_months_str <- as.character(past_12_months)
  
  selected_rasters <- temp_rast[[which(names(temp_rast) %in% past_12_months_str)]]
  
  if (is.null(selected_rasters) || nlyr(selected_rasters) == 0) {
    warning(paste("No raster data available for region", this_region, "and date", this_date))
    next
  }
  
  extracted_values <- terra::extract(
    selected_rasters,
    data_date %>% select(longitude, latitude), method = 'simple'
  )
  
  if (nrow(extracted_values) == 0 || ncol(extracted_values) <= 1) {
    warning(paste("No extracted data for region", this_region, "and date", this_date))
    next
  }
  
  extracted_matrix <- as.matrix(extracted_values[, -1])
  mean_temp <- colMeans(extracted_matrix, na.rm = TRUE) # get the mean per month slice in the previous 12 months
  
  coldest_month_idx <- which.min(mean_temp) # get the month with the coldest temp in the previous 12 months
  warmest_month_idx <- which.max(mean_temp) # get the month with the warmest temp in the previous 12 months
  
  coldest_monthT <- extracted_matrix[, coldest_month_idx]
  warmest_monthT <- extracted_matrix[, warmest_month_idx]
  
  data_date <- data_date %>%
    mutate(
      mean_temp = rowMeans(extracted_matrix, na.rm = TRUE),
      #sd_temp = apply(extracted_matrix, 1, sd, na.rm = TRUE),
      coldest_temp = coldest_monthT,
      warmest_temp = warmest_monthT
    )
  
  temp_list[[this_date]] <- data_date
}

grid <- bind_rows(temp_list) %>% filter(!is.na(mean_temp)) # remove if NA values, just to be safe

saveRDS(grid, here("R/data/processed/prediction_grid.rds"))


# ### Plot and save survey grids ----------------------------------------------
# 
# for (this_region in names(grids)) {
#   this_region <- "GOA"
#   df <- grids[[this_region]]  # Get survey grid data
#   
#   # Get the CRS from the data
#   crs_df <- sdmTMB::get_crs(df, c("longitude", "latitude"))
#   
#   # Load and transform coastline
#   coastline <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")
#   coastline_transformed <- st_transform(coastline, crs = crs_df)
#   coastline_transformed <- st_make_valid(coastline_transformed)
#   # Convert X and Y to meters for plotting
#   df <- df |> mutate(X_meters = X * 1000, Y_meters = Y * 1000)
#   
#   # Create and save plot
#   plot <- ggplot(data = coastline_transformed) +
#     geom_sf() +
#     geom_raster(data = df, aes(x = X_meters, y = Y_meters, fill = depth), inherit.aes = FALSE) +
#     xlim(min(df$X_meters), max(df$X_meters)) +
#     ylim(min(df$Y_meters), max(df$Y_meters)) +
#     ggtitle(paste0(this_region)) +
#     scale_fill_viridis_c()
#   
#   plot
#   
#   ggsave(plot = plot, filename = paste0(here('output/region_domain/'), this_region, '.png'))
# }

