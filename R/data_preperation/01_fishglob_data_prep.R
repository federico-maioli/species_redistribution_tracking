# libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(tidylog)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyterra)
library(terra)

# load survey data --------------------------------------------------------
fishglob_url <- "https://github.com/AquaAuma/FishGlob_data/raw/d71dfa03c2912b4e9d9cd10412ae2af52ba56ae5/outputs/Compiled_data/FishGlob_public_std_clean.RData" # stick with this version
options(timeout = 300) # increase timeout to 5 minutes
load(url(fishglob_url))

# exclude surveys -------------------------------------------------
excl_surveys <- c('GSL-N','GSL-S','ROCKALL','AI','DFO-SOG','WCTRI','SEUS','FR-CGFS','SP-PORC') # sampled inconsistently, rarely or started recently

data <- data %>%
  filter(!survey %in% excl_surveys) # exclude the surveys
#   mutate(survey = case_when(
#     survey %in% c("DFO-HS", "DFO-QCS", "DFO-WCHG", "DFO-WCVI") ~ "BC",
#     TRUE ~ survey
#   )) # British Columbia survey go together

# remove bad hauls --------------------------------------------------------
# most of the following filters are based on https://www.nature.com/articles/s41586-023-06449-y https://github.com/afredston/marine_heatwaves_trawl
# some are haul removal is redundant, but better safe than sorry
goa_del <- data %>% filter(survey == 'GOA' & year < 1994) %>% pull(haul_id) |> unique() # we don't have temperature data before 1993

gmex_del <- data %>% filter(survey == 'GMEX' & longitude > -87.9) %>% pull(haul_id) |> unique() # remove hauls on the east part of GMEX, they started in 2008

neus_del <- data %>%
  filter((survey == "NEUS" & year < 2009 & (haul_dur < 0.42 | haul_dur > 0.58)) |
           (survey == "NEUS" & year >= 2009 & (haul_dur < 0.25 | haul_dur > 0.42)) & (year < 1968 | year > 2019)) %>%
  pull(haul_id) |> unique() # based on this https://github.com/fishglob/FishGlob_data/issues/55

bits_del <- data %>% filter(survey == 'BITS' & year < 1999) %>% pull(haul_id) |> unique() # consistent spatial coverage after 1999

nor_bts_del <- c(
  data %>% filter(survey == 'Nor-BTS' & year < 2004) %>% pull(haul_id), # starts in 2004
  data %>% filter(survey == 'Nor-BTS' & latitude < 68) %>% pull(haul_id) # only north of 68N
) |> unique()

evhoe_del <- data %>% filter(survey == 'EVHOE' & year == 2017) %>% pull(haul_id) |> unique() # 2017 has very few hauls

nigfs_del <- data %>% filter(survey == 'NIGFS' & year < 2008) %>% pull(haul_id) |> unique() # starts in 2008

ns_ibts_del <- data %>% filter(survey == 'NS-IBTS' & year < 1980) %>% pull(haul_id) |> unique() # consistent spatial coverage after 1980

swc_ibts_del <- data %>% filter(survey == 'SWC-IBTS' & year == 1995 & longitude > -5 & latitude < 50.5 & latitude > 49) %>% pull(haul_id) |> unique() # remove some outliers in 1995

haul_ids_del <- c(goa_del, gmex_del, neus_del, bits_del, nor_bts_del,
                  evhoe_del, ns_ibts_del, swc_ibts_del, nigfs_del) |> unique() # pull together bad hauls

data <- data %>% filter(!haul_id %in% haul_ids_del) # remove bad hauls

# fix missing or zero catch weights ---------------------------------------

data <- data %>%
  mutate(wgt_cpua = if_else(survey == "NEUS", wgt / area_swept, wgt_cpua)) %>% # recompute NEUS CPUE based on https://github.com/fishglob/FishGlob_data/issues/55
  filter(wgt_cpua != 0) # remove zero CPUE

# remove pelagic families -------------------------------------------------

pelagic_families <- c("Clupeidae", "Osmeridae", "Exocoetidae", "Atherinidae",
                      "Engraulidae", "Hemiramphidae", "Inermiidae", "Belonidae",
                      "Scomberesocidae", "Echeneidae", "Carangidae", "Bramidae",
                      "Scombridae", "Centrolophidae", "Istiophoridae", "Ammodytidae")

data <- data %>% filter(!family %in% pelagic_families)

# remove duplicates -------------------------------------------------------

group_by_cols <- c("survey", "source", "timestamp", "haul_id", "country", "sub_area",
                   "continent", "stat_rec", "station", "stratum", "year", "month", "day",
                   "quarter", "latitude", "longitude", "haul_dur", "area_swept", "gear",
                   "sbt", "sst", "depth", "accepted_name")

data <- data %>%
  group_by(across(all_of(group_by_cols))) %>%
  summarize(
    wgt_cpua = sum(wgt_cpua, na.rm = TRUE),
    wgt = sum(wgt, na.rm = TRUE),
    .groups = "drop"
  )

# add region info ---------------------------------------------------------

data$month <- as.numeric(data$month)

data <- data %>%
  mutate(region = case_when(
    survey %in% c('GOA') ~ "Gulf of Alaska",
    survey %in% c("DFO-HS", "DFO-QCS", "DFO-WCHG", "DFO-WCVI") ~ "British Columbia",
    survey %in% c('WCANN') ~ "U.S. West Coast",
    survey %in% c('BITS') ~ "Baltic Sea",
    survey %in% c('NS-IBTS') ~ "North Sea",
    survey %in% c('Nor-BTS') ~ "Barents Sea",
    survey %in% c("EBS") ~ "East Bering Sea",
    survey %in% c("GMEX") ~ "Gulf of Mexico",
    survey %in% c("EVHOE","IE-IGFS","SWC-IBTS", "NIGF") ~ "Celtic-Biscay Shelf",
    survey %in% c("SP-NORTH") ~ "North Iberian Coast",
    survey %in% c("SCS","NEUS") ~ "Northeast U.S. and Scotian Shelf",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(region)) %>%
  mutate(region_short = case_when(
    region == "Gulf of Alaska" ~ "GOA",
    region == "British Columbia" ~ "BC",
    region == "U.S. West Coast" ~ "USWC",
    region == "Baltic Sea" ~ "BAL",
    region == "North Sea" ~ "NS",
    region == "Barents Sea" ~ "BS",
    region == "East Bering Sea" ~ "EBS",
    region == "Gulf of Mexico" ~ "GOM",
    region == "Celtic-Biscay Shelf" ~ "CBS",
    region == "North Iberian Coast" ~ "NIC",
    region == "Northeast U.S. and Scotian Shelf" ~ "NEUS-SS",
    TRUE ~ NA_character_
  ))

# reorder columns ---------------------------------------------------------
data <- data %>%
  select(survey, region, region_short, source, timestamp, haul_id, country,
         sub_area, continent, stat_rec, year, month, day, quarter, latitude,
         longitude, haul_dur, area_swept, gear, sbt, sst, depth, accepted_name,
         wgt_cpua, wgt)

# complete dataset with all haul-species combos ---------------------------
#For each region, make sure every haul has a row for every species, even if it wasn’t caught, set its weight per unit area to zero if missing, fill in any other missing info, and keep only data from 1994 onward

metadata_cols <- c(
  "survey", "region", "region_short", "source", "timestamp", "country",
  "sub_area", "continent", "stat_rec", "year", "month", "day",
  "quarter", "latitude", "longitude", "haul_dur", "area_swept",
  "gear", "sbt", "sst", "depth"
) # get haul level data

haul_info <- data %>%
  select(all_of(metadata_cols), haul_id) %>%
  distinct(haul_id, .keep_all = TRUE) # ensures one row per haul

expand_data <- data %>%
  select(region, haul_id, accepted_name, wgt_cpua) %>%
  group_by(region) %>%
  complete(haul_id, accepted_name, fill = list(wgt_cpua = 0)) %>%
  ungroup() # get haul_id - species combinations and add 0 wgt_cpua for missing combos

data_complete <- expand_data %>%
  left_join(haul_info, by = c("region", "haul_id")) %>%  filter(year >= 1994) # takes a while!

# retain only common species -----------------------------------------------

# set thresholds
min_freq_threshold <- 15 # minimum % frequency of occurrence in a region
cumulative_biomass_threshold <- 0.99  # 99% of biomass per region

selected_species <- data_complete %>%
  # Add presence flag
  mutate(present = as.integer(wgt_cpua > 0)) %>%
  # Species metrics
  group_by(region, accepted_name) %>%
  summarise(
    years_present  = n_distinct(year[present == 1]),              # distinct years present
    positive_hauls = sum(present),                                # # of positive hauls
    mean_wgt_cpua  = mean(wgt_cpua[present == 1], na.rm = TRUE),  # mean biomass in positive hauls
    total_biomass  = sum(wgt_cpua, na.rm = TRUE),                 # total biomass density
    total_hauls_region = n_distinct(haul_id),                      # total hauls in this region
    .groups = "drop_last"
  ) %>%
  
  # Apply cumulative biomass and frequency filters
  arrange(region, desc(total_biomass)) %>%
  mutate(
    cumulative_biomass = cumsum(total_biomass) / sum(total_biomass),
    freq_occurrence = 100 * positive_hauls / total_hauls_region
  ) %>%
  filter(
    cumulative_biomass <= cumulative_biomass_threshold,
    freq_occurrence >= min_freq_threshold
  ) %>%
  ungroup() %>%
  # Remove genus-only taxa
  filter(str_detect(accepted_name, "\\s")) %>%
  select(region, accepted_name) %>%
  # Filter data_complete to only selected species
  semi_join(data_complete, by = c("region", "accepted_name"))

# filter dataset
filtered_data <- data_complete %>%
  semi_join(selected_species, by = c("region", "accepted_name"))

rm(list = setdiff(ls(), "filtered_data"))

# add depth if missing ----------------------------------------------------

# depth_raster <- terra::rast(here("data/environmental/depth/GEBCO_2023_sub_ice_topo.nc"))
# 
# filtered_data$depth_gebco <- terra::extract(depth_raster$elevation,
#                                             filtered_data %>% select(longitude, latitude))$elevation %>%
#   abs()
# 
# # number of trawls with missing depth info
# sum(is.na(filtered_data$depth))
# 
# filtered_data <- filtered_data %>%
#   mutate(depth = coalesce(depth, depth_gebco)) %>%
#   select(-depth_gebco)


# match temperature, this is needed only for the thermal envelope models -------------------------------------------------------

# need to add minimun month of the survey for matching temperature
filtered_data = filtered_data |> group_by(region) |> mutate(min_month = min(month,na.rm = TRUE)) |> ungroup() |> mutate(
  month_year_date = as.Date(paste0(year, "-", sprintf("%02d", min_month), "-01"))
)

# load temperature nc file ------------------------------------------------

hauls = filtered_data |> distinct(region,haul_id,longitude,latitude,year,month_year_date) # simplifies extraction

nc_path <- here('data/environmental/temperature/cmems_mod_glo_phy_my_0.083deg_P1M-m_bottomT_180.00W-179.92E_20.00N-87.00N_1993-01-01-2021-06-01.nc')

# source time to date function for .nc files
source(here('utils/nc_time_to_date.R'))

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

unique_dates <- sort(unique(hauls$month_year_date))

for (i in seq_along(unique_dates)) {
  
    this_date <- unique_dates[i]
    
    data_date <- filter(hauls, month_year_date == this_date)
    
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

hauls_temp <- bind_rows(temp_list)

# remove hauls with missing temperature
#hauls_temp <- hauls_temp %>% filter(!is.na(mean_temp)) # it's not needed 

# rejoin tempearture to the main dataset
data <- filtered_data %>% left_join(hauls_temp %>% select(haul_id, region, mean_temp, coldest_temp, warmest_temp), by = c('haul_id','region')) #%>% filter(!is.na(mean_temp))

filtered_data %>%
  anti_join(hauls_temp, by = c("haul_id", "region")) # get unmatched data

# final selection and save -----------------------------------------------
data <- data %>%
  select(survey, region, region_short, haul_id, country, continent,
         year, quarter, month, day, min_month, month_year_date, latitude, longitude, haul_dur,
         area_swept, gear, mean_temp, coldest_temp, warmest_temp, sbt, sst, depth, accepted_name, wgt_cpua)

write_rds(data, here("data/processed/fishglob_clean.rds"))

