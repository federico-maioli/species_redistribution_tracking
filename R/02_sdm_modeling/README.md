# Fitting species distribution models to survey data

This markdown accompanies the scripts in the `02_sdm_modeling/` directory. It explains how we fit spatial-temporal species distribution models (SDMs) to bottom trawl survey data using the `sdmTMB` package in R.

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
  formula = wgt_cpua ~ 0 + as.factor(year) + logdepth + logdepth2, # + as.factor(survey) / + as.factor(quarter)
  mesh = spde,
  data = sub,
  time = "year",
  family = tweedie(link = "log"),
  spatial = "on",
  spatiotemporal = "iid",  # or "rw"
  share_range = TRUE
)
```
Fitted models are in `02_sdm_modeling/fitted/` directory, named by region and species (not on Github because of size).

In `02_sdm_modeling/02_derive_quantities.R`, we derive annual estimates (with standard errors) of range centroid, depth centroid and mean occupied bottom temperature per population, as biomass-weighted averages over the prediction grid. We also compute `thermal_niche_constant_density`: the same weighted average but with the density field fixed at its across-year mean, so temperature changes alone drive it. Comparing the two separates redistribution from in-place warming (see `03_bayesian_trends/08_static_vs_dynamic_occupied_temp.R`).

Output: `data/processed/derived_quantities_sdm.rds`.
