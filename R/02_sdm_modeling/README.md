# Fitting species distribution models to survey data

This markdown accompanies the scripts in the `sdm_modeling/` directory. It explains how we fit spatial-temporal species distribution models (SDMs) to bottom trawl survey data using the `sdmTMB` package in R.

We model species-specific distributions and relative abundance over space and time, including:

-   Depth (log-transformed and standardized)
-   Spatial and spatiotemporal random fields
-   Survey fixed effects and month (random intercept), when it's the case

These models produce predictions used to derive quantities for trend estimation and multivariate analyses.

## Workflow

1.  Load and filter survey data by species-region (our "populations")
2.  Standardize covariates (e.g. log-depth)
3.  Convert coordinates to UTM and construct spatial mesh
4.  Fit a GLMM with spatial and temporal components
5.  Check diagnostics — rerun if max gradient \> 0.001
6.  Save model only if it passes `sdmTMB::sanity()`

Models are fit in parallel across species and regions.

## Model formulation in `sdmTMB` terms

The general model structure:

-   Tweedie family (for zero-inflated biomass), with log link
-   Spatial random fields: `spatial = "on"`
-   Spatiotemporal random fields:
    -   `"iid"` for most regions
    -   `"rw"` (random walk) for BC

### Formula (simplified):

``` r
fit <- sdmTMB(
  formula = wgt_cpua ~ 0 + as.factor(year) + logdepth + logdepth2, # + as.factor(survey) / + (1 | month_f)
  mesh = spde,
  data = sub,
  time = "year",
  family = tweedie(link = "log"),
  spatial = "on",
  spatiotemporal = "iid",  # or "rw"
  share_range = TRUE
)
```
Fitted models are in `sdm_modeling/fitted/` directory, named by region and species (not on Github because of size).

In `sdm_modeling/02_derive_quantities.R`, we derive annual estimates of range centroid, depth niche and thermal niche per species. We save these in `data/processed/derived_quantities/` for further analyses.