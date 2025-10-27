# libraries ---------------------------------------------------------------
library(tidyverse)
library(here)
library(sdmTMB)
library(sf)
library(terra)
library(ncdf4)
library(utils)

# load function and data -------------------------------------------------

source(here('R/utils/exclude_too_far.R'))
data <- readRDS(here('R/data/processed/fishglob_clean.rds'))

# prepare survey spatial grid ----------------------------------------------

unique_locations <- data %>%  distinct(region_short, depth, latitude, longitude)
#unique_locations <- data %>%  distinct(region_short, depth, latitude, longitude, year)
split_data <- split(unique_locations, unique_locations$region_short)

# grid settings
grids <- list()
cell_width <- 4  # grid resolution (4 km)

round_down_even <- function(x, base = 1000) base * floor(x / base)


## coastline ---------------------------------------------------------------
# not on cran so need manual download :(
# URL for the Natural Earth 1:10m countries shapefile
url <- "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip"
# Temporary paths
temp <- tempfile(fileext = ".zip")
unzip_dir <- tempdir()
# Download
download.file(url, temp, mode = "wb")
# Unzip *all* files, not just the .shp
unzip(temp, exdir = unzip_dir)
# Find the shapefile (the .shp file inside the extracted folder)
shp_path <- list.files(unzip_dir, pattern = "\\.shp$", full.names = TRUE)

# Read and transform
coastline <- st_read(shp_path) |>
  st_transform(crs = 4326) |>
  st_make_valid()
 

## depth -------------------------------------------------------------------
# get it from GEBCO https://www.gebco.net/data-products/gridded-bathymetry-data
depth_raster <- terra::rast(here('R/data/environmental/depth/GEBCO_2023_sub_ice_topo.nc'))

depth_raster[depth_raster > 0] <- NA # remove positive values

depth_raster <- terra::focal(
  depth_raster,
  w = 5,                 # 5x5 neighborhood
  fun = mean,             # average nearby cell values
  na.rm = TRUE,
  na.policy = "only"      # ONLY fill NA cells
) # takes quite a while!!

names(depth_raster) <- 'elevation' # focal changes name for some reasons

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

nc_path <- here('R/data/environmental/temperature/cmems_mod_glo_phy_my_0.083deg_P1M-m_bottomT_180.00W-179.92E_20.00N-87.00N_1993-01-01-2024-01-01.nc')

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
  w = 7,                 # 7x7 neighborhood
  fun = mean,             # average nearby cell values
  na.rm = TRUE,
  na.policy = "only"      # ONLY fill NA cells
)

temp_list <- list()

unique_dates <- sort(unique(grid$month_year_date))

for (i in seq_along(unique_dates)) {
  
  this_date <- unique_dates[i]
  this_year <- year(this_date)
  
  data_date <- filter(grid, month_year_date == this_date)
  
  # past 12 months
  past_12_months <- seq(from = this_date %m-% months(12),
                        to   = this_date %m-% months(1),
                        by   = "1 month")
  past_12_months_str <- as.character(past_12_months)
  
  selected_rasters_past <- temp_rast[[which(names(temp_rast) %in% past_12_months_str)]]
  
  if (is.null(selected_rasters_past) || nlyr(selected_rasters_past) == 0) {
    warning(paste("No raster data available for region", this_region, "and date", this_date))
    next
  }
  
  extracted_values_past <- terra::extract(
    selected_rasters_past,
    data_date %>% select(longitude, latitude),
    method = 'simple'
  )
  
  extracted_matrix_past <- as.matrix(extracted_values_past[, -1])
  mean_temp <- rowMeans(extracted_matrix_past, na.rm = TRUE)
  
  monthly_means_past <- colMeans(extracted_matrix_past, na.rm = TRUE)
  coldest_month_idx <- which.min(monthly_means_past)
  warmest_month_idx <- which.max(monthly_means_past)
  
  coldest_monthT <- extracted_matrix_past[, coldest_month_idx]
  warmest_monthT <- extracted_matrix_past[, warmest_month_idx]
  
  # calendar year months
  this_year_months <- seq(from = as.Date(paste0(this_year, "-01-01")),
                          to   = as.Date(paste0(this_year, "-12-01")),
                          by   = "1 month")
  this_year_months_str <- as.character(this_year_months)
  
  selected_rasters_year <- temp_rast[[which(names(temp_rast) %in% this_year_months_str)]]
  
  mean_year_temp <- NA  # default
  
  extracted_values_year <- terra::extract(
    selected_rasters_year,
    data_date %>% select(longitude, latitude),
    method = 'simple'
  )
  extracted_matrix_year <- as.matrix(extracted_values_year[, -1])
  
  # compute the mean across all 12 months of that year
  mean_year_temp <- rowMeans(extracted_matrix_year, na.rm = TRUE)
  
  # combine everything
  data_date <- data_date %>%
    mutate(
      mean_temp       = mean_temp,        # past 12-month mean
      coldest_temp    = coldest_monthT,
      warmest_temp    = warmest_monthT,
      mean_year_temp  = mean_year_temp    # mean of 12 months in same calendar year
    )
  
  temp_list[[this_date]] <- data_date
}

grid <- bind_rows(temp_list) %>% filter(!is.na(mean_temp))

# add quarter and survey info
survey_map <- data %>%
  group_by(region_short, survey) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(region_short) %>%
  slice_max(n, with_ties = FALSE) %>%
  select(region_short, survey)

quarter_map <- data %>%
  filter(month == min_month) %>%
  distinct(region_short, min_month = month, quarter) %>% select(region_short,quarter)

# join back to the grid
grid <- grid %>%
  left_join(survey_map, by = "region_short") %>%
  # Join on both region and month for quarter info
  left_join(quarter_map, by = "region_short")

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

