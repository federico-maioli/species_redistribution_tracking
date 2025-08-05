# load libraries ----------------------------------------------------------
library(tidyverse)
library(here)
library(sdmTMB)
library(sf)
library(terra)
library(lubridate)

# load data ---------------------------------------------------------------
data <- readRDS(here("data/trawl_surveys/fishglob_clean.rds"))

data <- data %>%
  mutate(
    month = as.numeric(month),
    quarter = quarter(month)
  ) %>%
  group_by(region) %>%
  mutate(min_month = min(month)) %>%
  ungroup()

# load grid ---------------------------------------------------------------
grid <- readRDS(here("data/spatial_grid/grid.rds")) %>%
  bind_rows(.id = "region_short") %>%
  filter(region_short %in% unique(data$region_short)) %>%
  left_join(distinct(data, region_short, min_month), by = "region_short")

grid <- replicate_df(grid, time_name = "year", time_values = 1994:2021) %>%
  mutate(
    month_year_date = as.Date(paste0(year, "-", sprintf("%02d", min_month), "-01"))
  )

# load temperature rasters ------------------------------------------------
temp_files <- list.files(here("data/environmental/temperature"), full.names = TRUE)
raster_list <- list()

for (file in temp_files) {
  file_date <- gsub(".*bottomT_|\\.csv", "", basename(file))
  file_date <- as.Date(paste0("01-", file_date), format = "%d-%m-%Y")
  
  temp_tibble <- read.csv(file) %>%
    select(longitude, latitude, bottomT)
  
  temp_points <- vect(temp_tibble, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  template_raster <- rast(temp_points, resolution = 0.083)
  raster_layer <- rasterize(temp_points, template_raster, field = "bottomT", fun = mean)
  
  raster_list[[as.character(file_date)]] <- raster_layer
}

temp_raster_stack <- rast(raster_list)

plot(temp_raster_stack)

# process regions ---------------------------------------------------------
temp_list <- list()

for (this_region in unique(grid$region_short)) {
  grid_region <- filter(grid, region_short == this_region)
  
  unique_dates <- sort(unique(grid_region$month_year_date))
  
  for (i in seq_along(unique_dates)) {
    
    this_date <- unique_dates[i]
    
    grid_region_date <- filter(grid_region, month_year_date == this_date)
    
    past_12_months <- seq(from = this_date %m-% months(12), to = this_date %m-% months(1), by = "1 month")
    past_12_months_str <- as.character(past_12_months)
    
    selected_rasters <- temp_raster_stack[[which(names(temp_raster_stack) %in% past_12_months_str)]]
    
    if (is.null(selected_rasters) || nlyr(selected_rasters) == 0) {
      warning(paste("No raster data available for region", this_region, "and date", this_date))
      next
    }
    
    extracted_values <- terra::extract(
      selected_rasters,
      grid_region_date %>% select(longitude, latitude)
    )
    
    if (nrow(extracted_values) == 0 || ncol(extracted_values) <= 1) {
      warning(paste("No extracted data for region", this_region, "and date", this_date))
      next
    }
    
    extracted_matrix <- as.matrix(extracted_values[, -1])
    region_mean_temp <- colMeans(extracted_matrix, na.rm = TRUE)
    
    coldest_month_idx <- which.min(region_mean_temp)
    warmest_month_idx <- which.max(region_mean_temp)
    
    coldest_monthT <- extracted_matrix[, coldest_month_idx]
    warmest_monthT <- extracted_matrix[, warmest_month_idx]
    
    grid_region_date <- grid_region_date %>%
      mutate(
        mean_temp = rowMeans(extracted_matrix, na.rm = TRUE),
        sd_temp = apply(extracted_matrix, 1, sd, na.rm = TRUE),
        coldest_temp = coldest_monthT,
        warmest_temp = warmest_monthT
      )
    
    temp_list[[paste(this_region, this_date, sep = "_")]] <- grid_region_date
  }
}

# combine and save results ------------------------------------------------
grid <- bind_rows(temp_list)

saveRDS(grid, here("data/spatial_grid/grid_with_temp.rds"))

