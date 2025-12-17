# Bayesian trend analysis

---

This README accompanies the scripts in the `bayesian_trends/` directory. We use a hierarchical Bayesian framework to estimate temporal trends in latitudinal and longitudinal centroids, depth niche, and thermal niche. The models are fitted using the `brms` package, which interfaces with Stan for Bayesian inference.
We propagate standard errors from the spatiotemporal modeling into the trend analysis to account for observation-level uncertainty. The hierarchical structure allows us to estimate global, regional, and species-level trends within a single unified model.

**Data locations:**

- Input data: `R/data/processed/derived_quantities_sdm.rds`

- Fitted models: `R/bayesian_trends/fitted/`

- Estimated trends: `R/data/processed/estimated_trends/`

the main script are:

### `01_fit_bayesian_trends.R`

<details>
  <summary>Click here to expand</summary>

- Fits hierarchical Bayesian models for all response variables using a **Student-t likelihood**, robust to extreme values.
- Includes:
  - **Global slope** for time (decade effect)
  - **Region-level deviations**
  - **Species-level deviations nested within region**
- The term `|ID|gr(region:species, by = region)` forces slopes to be correlated within regions.
- Each response variable is modeled with its standard error (`se(...)`) from the spatiotemporal model.
- Example syntax:

```r
brm(
  formula = bf(
    cog_y_c | se(cog_y_se, sigma = TRUE) ~ 0 + year_c + 
      (0 + year_c | region) + 
      (0 + year_c | ID | gr(region:species, by = region))
  ) +
  bf(...) + ... # other responses
  set_rescor(FALSE),
  data = data,
  family = student()
)
```
 - Centered response variables and time predictor per species–region combination so that slopes represent deviations from species–region mean.
 - Time rescaled to decades for interpretable slope units (change per decade).


### `02_model_diagnostics.R`

<details>
  <summary>Click here to expand</summary>


-	Performs posterior predictive checks (PPCs) to evaluate fit.
- Checks convergence using:
	-	$\hat{R}$
	-	Bulk and tail effective sample size (ESS)
	-	Divergent transitions
	-	Visualizes diagnostics and confirms good enough sampling.


### `03_prior_diagnostics.R`

<details>
  <summary>Click here to expand</summary>
  
	- Examines priors for global, regional, and species-level slopes.
	- Compares priors against literature and empirical distributions from the data.



### `04_plot_model_slopes.R`

<details>
  <summary>Click here to expand</summary>

-	Generates figures of estimated slopes for species and regions.


### `05_create_macro_latex.R`

<details>
  <summary>Click here to expand</summary>

- Automates LaTeX Output
	-	Extracts posterior summaries.
	-	Produces LaTeX macros automatically for manuscript reporting.
Code is from Sean Anderson 




