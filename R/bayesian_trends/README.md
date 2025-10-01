# Fitting bayesian trend analysis

This markdown accompanies the scripts in the `bayesian_trends/` directory. It explains how we fit Bayesian mixed effects models to estimate trends in species distributional responses. The models are fitted using the `brms` package, which uses Stan for Bayesian inference.
Data are in `data/processed/derived_quantities.rds` and the models fitted are in `bayesian_trends/fitted`. The results are saved in `data/processed/trend.rds`.

We fitted different model but with the same basic structure: 
outcome are modelled jointly. They are mean-centered so they represent anomalies. Also year is centered so that we can interpret the slope as the change in the outcome per year and we can omit the intercept.



