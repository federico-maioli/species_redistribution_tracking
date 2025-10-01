library(tidyverse)
library(here)
library(terra)

# Get the prediction grid

grid <- readRDS(here("data/processed/prediction_grid.rds")) %>%
  mutate(month_f = as.factor(min_month)) %>%
  droplevels()


#  try with the baltic sea ------------------------------------------------


# Filter only BAL region
df_bal <- grid %>% filter(region_short == "BAL")


years_bal <- sort(unique(df_bal$year))

r_list <- list()

for (yr in years_bal) {
  df_yr <- df_bal %>% filter(year == yr)
  
  # Create raster using X and Y directly
  r <- rast(
    extent = c(min(df_yr$X), max(df_yr$X),
               min(df_yr$Y), max(df_yr$Y)),
    ncol = length(unique(df_bal$X)),
    nrow = length(unique(df_bal$Y)),
    crs = NA # no CRS since it's in a local metric grid
  )
  
  # Assign temperature to correct cells
  cell_idx <- cellFromXY(r, cbind(df_yr$X, df_yr$Y))
  values(r)[cell_idx] <- df_yr$mean_temp
  
  r_list[[as.character(yr)]] <- r
}

# Combine into raster stack
r_stack <- rast(r_list)
names(r_stack) <- paste0("year_", years_bal)

# Check resolution
print(r_stack)
plot(r_stack[[2]], main = paste("BAL Temperature", years_bal[2]))

# Custom function to fit linear trend
slope_fun <- function(v) {
  if (all(is.na(v))) return(NA_real_)
  ok <- !is.na(v)
  if (sum(ok) < 2) return(NA_real_)
  m <- lm(v[ok] ~ years_bal[ok])
  coef(m)[2]  # slope in °C/year
}

# Apply to all pixels
slope_per_year <- app(r_stack, fun = slope_fun)

# Convert to °C/decade
slope_per_decade <- slope_per_year * 10

plot(slope_per_decade, main = "BAL Temperature Trend (°C/decade)")

# Select last two years for spatial gradient
last_two <- 2019:2020 # get it back as before

cat("Last two years:", last_two, "\n")

mean_last_two <- mean(r_stack[[paste0("year_", last_two)]])

# dx and dy are grid spacing in km
dx <- xres(mean_last_two)
dy <- yres(mean_last_two)

grad_vector_fun <- function(v) {
  if (all(is.na(v))) return(c(NA_real_, NA_real_))
  
  east <- v[6]; west <- v[4]
  north <- v[2]; south <- v[8]
  
  if (any(is.na(c(east, west, north, south)))) return(c(NA_real_, NA_real_))
  
  dx_val <- (east - west) / (2 * dx)  # °C/km
  dy_val <- (north - south) / (2 * dy) # °C/km
  
  c(dx_val, dy_val)
}

# Apply focal with multiple outputs
grad_stack <- focal(mean_last_two, w = matrix(1, 3,3), fun = grad_vector_fun, na.policy="omit")

grad_eps <- 1e-9  # avoid dividing by zero

dx_layer <- grad_stack[[1]]
dy_layer <- grad_stack[[2]]

Vx_km_decade <- slope_per_decade / dx_layer
Vy_km_decade <- slope_per_decade / dy_layer

# Set small gradients to NA
Vx_km_decade[!is.na(dx_layer) & abs(dx_layer) < grad_eps] <- NA
Vy_km_decade[!is.na(dy_layer) & abs(dy_layer) < grad_eps] <- NA

# Convert to meters/decade
Vx_m_decade <- Vx_km_decade * 1000
Vy_m_decade <- Vy_km_decade * 1000
velocity_m_decade <- sqrt(Vx_m_decade^2 + Vy_m_decade^2)

angle_rad <- atan2(Vy_km_decade, Vx_km_decade)  # radians
angle_deg <- angle_rad * 180 / pi

# 0° = east, 90° = north, 180° = west, -90° = south
plot(angle_deg, main="Direction of Climate Velocity (degrees)")

# Assuming Vy_m_decade is your northward velocity raster
# Extract values
vals_north <- values(Vy_km_decade)
vals_non_na <- vals_north[!is.na(vals_north)]

# Calculate 0.5% and 99.5% quantiles
lo <- quantile(vals_non_na, 0.005, na.rm = TRUE)
hi <- quantile(vals_non_na, 0.995, na.rm = TRUE)

# Clamp northward velocity
Vy_clamped <- clamp(Vy_m_decade, lower = lo, upper = hi, values = TRUE)
# Diverging color palette: blue = south, white = no shift, red = north
cols <- colorRampPalette(c("blue", "white", "red"))(50)

plot(Vy_clamped,
     main = "BAL Northward Climate Velocity (m/decade, clamped)",
     col = cols)
