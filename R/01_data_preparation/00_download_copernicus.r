library(here)
# code from https://help.marine.copernicus.eu/en/articles/8638253-how-to-download-data-via-the-copernicus-marine-toolbox-in-r
#first you need to get copernicus marine toolbox executable file https://help.marine.copernicus.eu/en/articles/10750437-copernicus-marine-toolbox-executable-no-installation

path_copernicusmarine <- "/Users/fedma/copernicus/copernicusmarine.cli" # change with you own path!!
# e.g. path_copernicusmarine <- "C:/Users/username/Downloads/copernicusmarine.exe"

file.exists(path_copernicusmarine) # check if it exists

output_directory <- here('R/data/environmental/temperature/')

command <- paste(
  shQuote(path_copernicusmarine),
  "subset",
  "--dataset-id cmems_mod_glo_phy_my_0.083deg_P1M-m",
  "--variable bottomT",
  "--start-datetime 1993-01-01T00:00:00",
  "--end-datetime 2024-01-01T00:00:00",
  "--minimum-longitude -180",
  "--maximum-longitude 179.9166717529297",
  "--minimum-latitude 20",
  "--maximum-latitude 87",
  "--minimum-depth 0.49402499198913574",
  "--maximum-depth 0.49402499198913574",
  "-o", shQuote(output_directory),
  sep = " "
)

system(command)

# updated version, no need for installation
# follow this guide https://help.marine.copernicus.eu/en/articles/10750437-copernicus-marine-toolbox-executable-no-installation
# https://help.marine.copernicus.eu/en/articles/8638253-how-to-download-data-via-the-copernicus-marine-toolbox-in-r

# dowload the binaries



