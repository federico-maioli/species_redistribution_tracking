# Project Overview

This repository contains all analyses and manuscript files for the study “Thermal niche warming is more con-
sistent than range shifts in marine species under climate change”.

*Authors: Federico Maioli (fedma@aqua.dtu.dk), Daniël van Denderen, Max Lindmark, Marcel Montanyès, Eric J. Ward, Sean C. Anderson, Martin Lindegren*

# Repository structure

 - `R/`
Contains all R scripts used for the data analysis. This folder includes code for data cleaning, statistical analysis, and generating outputs.

 - `tex/`
Contains the LaTeX files used for writing the manuscript. This includes the main text, sections, references, and any formatting needed for the final document.

 - `renv/`
Used to manage the project’s R environment. This ensures reproducibility by specifying the exact package versions used in the analysis. You shouldn't care about this.

 - `output/`
Contains all figures, tables, and numerical results produced by the analysis and reported in the text. This folder serves as the central location for outputs referenced in the manuscript.

# Data availability 

All datasets used in this project are publicly available from their respective sources: 

 - Fish biomass data from standardized, fishery-independent bottom-trawl surveys compiled by [FISHGLOB](https://fishglob.sites.ucsc.edu/)
 
 - Bottom temperature data from the [Copernicus Marine Environment Monitoring Service (CMEMS) Global Ocean Physics Reanalysis](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description)
 
 - Bathymetric data from the [GEBCO](https://www.gebco.net)
 
# Reproducibility

This project uses the `renv` package to ensure that all R package versions match those used to produce the results in the manuscript.

To reproduce the analysis:

1. Clone or download this repository.
2. Open the `.Rproj` file in RStudio.
3. Install `renv` (if not already installed):

```r
   install.packages("renv")
   renv::restore()
```

### **Note**

Downloading environmental data from Copernicus Marine Service may take some time, depending on connection speed. Fitting the species distribution models (`R/02_sdm_modeling/01_fit_sdm.R`) is computationally intensive and may require substantial runtime (approximately overnight on a MacBook M3, 16 GB RAM, 2023 model).

If your primary interest is the Bayesian trend analysis, you may start directly from:

`R/03_bayesian_trends/`

The required processed input data are available in the repository at:

`R/data/processed/derived_quantities_sdm.rds`

This allows reproduction of the Bayesian analyses without re-running the full SDM workflow.


A complete archived version of the repository is available on Zenodo [![DOI](https://zenodo.org/badge/1067761557.svg)](https://doi.org/10.5281/zenodo.18761947)
