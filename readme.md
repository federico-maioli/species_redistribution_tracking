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
