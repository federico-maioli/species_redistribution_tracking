# Data preparation

In this directory `data_preparation/` you can find the code to download and prepare the data used in this project. The processed data are saved in `data/processed/`. `data_preparation/00_download_copernicus.R` downloads the temperature data from Copernicus. `data_preparation/01_fishglob_data_prep.R` downloads, prepares the trawl survey data and matches model based temperature data to the trawl survey data. `data_preparation/02_create_prediction_grid.R` prepares the spatiotemporal grids used for predictions.

## Trawl survey data

The trawl survey data used in this project are sourced from [FishGlob](https://github.com/fishglob/FishGlob_data). The primary dataset and methodology are described in the following publication:

Maureaud, A.A., Palacios-Abrantes, J., Kitchel, Z. et al. FISHGLOB_data: an integrated dataset of fish biodiversity sampled with scientific bottom-trawl surveys. Sci Data 11, 24 (2024). <https://doi.org/10.1038/s41597-023-02866-w>

## Temperature data

This is needed for the thermal envelopes models only.
Temperature data are sourced from [Copernicus](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description) and can be downloaded using the [Copernicus Marine Toolbox](https://help.marine.copernicus.eu/en/articles/8638253-how-to-download-data-via-the-copernicus-marine-toolbox-in-r). These data are large and we therefore ...
