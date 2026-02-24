# Data preparation

The directory `01_data_preparation/` contains the scripts used to download and prepare all datasets required for the project. Because downloading the temperature data requires installation of the Copernicus Marine Toolbox and the preprocessing steps can be time-consuming, you may skip this directory if the processed files are already available in `data/processed/`.

- `01_data_preparation/00_download_copernicus.R` downloads bottom temperature data from Copernicus.
- `01_data_preparation/01_fishglob_data_prep.R` downloads and prepares the trawl survey data.
- `01_data_preparation/02_create_prediction_grid.R` creates the spatiotemporal prediction grids and matches temperature data to the survey data.

For a faster workflow, you may start directly from `03_bayesian_trends/`, which runs the hierarchical Bayesian trend analyses using the model-derived distribution metrics saved in `data/processed/`.

# Data

## Trawl survey data

The trawl survey data used in this project are sourced from [FishGlob](https://github.com/fishglob/FishGlob_data). The primary dataset and methodology are described in the following publication:

Maureaud, A.A., Palacios-Abrantes, J., Kitchel, Z. et al. FISHGLOB_data: an integrated dataset of fish biodiversity sampled with scientific bottom-trawl surveys. Sci Data 11, 24 (2024). <https://doi.org/10.1038/s41597-023-02866-w>

## Temperature data

Temperature data are sourced from [Copernicus](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description) and can be downloaded using the [Copernicus Marine Toolbox](https://help.marine.copernicus.eu/en/articles/8638253-how-to-download-data-via-the-copernicus-marine-toolbox-in-r).
