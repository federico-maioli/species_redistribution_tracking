library(DiagrammeR)
grViz("
digraph pipeline {

  graph [rankdir = LR, fontsize = 12]

  # Data sources
  node [shape = parallelogram]
  FISHGLOB
  Copernicus
  GEBCO

  # Model
  node [shape = box]
  'Spatiotemporal GLMM'
  
  # Extra
  node [shape = parallelogram]
  'Prediction grid'

  # Output
  node [shape = parallelogram]
  Predictions

  # Arrows
  FISHGLOB -> 'Spatiotemporal GLMM' -> 'Prediction grid'
  Copernicus -> 'Prediction grid'
  GEBCO -> 'Prediction grid'
  'Prediction grid' -> 'Predictions'
  'Predictions' -> 'Latitudinal range centroid'
  'Predictions' -> 'Longitudinal range centroid'
  'Predictions' -> 'Depth niche'
  'Predictions' -> 'Thermal niche'
  'Latitudinal range centroid' -> 'Bayesian trend analysis'
  'Longitudinal range centroid' -> 'Bayesian trend analysis'
  'Depth niche' -> 'Bayesian trend analysis'
  'Thermal niche'-> 'Bayesian trend analysis'
  'Bayesian trend analysis' -> 'Latitudinal range centroid shifts'
  'Bayesian trend analysis' -> 'Longitudinal range centroid shifts'
  'Bayesian trend analysis' -> 'Depth niche shifts'
  'Bayesian trend analysis' -> 'Thermal niche shifts'
}
")
