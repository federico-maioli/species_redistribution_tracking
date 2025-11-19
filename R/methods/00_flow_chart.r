library(DiagrammeR)


grViz("
digraph pipeline {

  graph [fontsize = 16, rankdir = TB, splines=ortho]

  # ---------------------------
  # BLOCK 0: Raw data (outside clusters)
  node [shape = parallelogram, style = filled, fillcolor = lightblue, fontname = Helvetica]
  fish [label='FishGlob\\n(Maureaud et al., 2024)']
  temp [label='Temperature\\nCopernicus Global Ocean Physics Reanalysis\\n(EU-Copernicus, 2018)']
  depth [label='Bathymetry\\n(GEBCO, 2023)']
  grid [label='Regular grid\\n5 x 5 km grid cells']

  # ---------------------------
  # BLOCK 1: Species spatiotemporal modelling
  subgraph cluster_sdm {
    label = 'Species spatiotemporal modelling'
    style = filled
    color = lightyellow

    node [shape = box, style = filled, fillcolor = antiquewhite, fontname = Helvetica]
    sdm [label='Fit species spatiotemporal\\nmodels (GLMMs)']
    predict [label='Spatiotemporal predictions']
  }

  # ---------------------------
  # BLOCK 2: Spatial & thermal metrics
  subgraph cluster_metrics {
    label = 'Spatial and thermal metrics'
    style = filled
    color = khaki1

    node [shape = box, style = filled, fillcolor = khaki1, fontname = Helvetica]
    derived_quantities [label='Spatial & thermal metrics:\\n- Range centroids\\n- Depth niche\\n- Thermal niche']
    pop_abundance [label='Population abundance']
  }

  # ---------------------------
  # BLOCK 3: Bayesian trend analysis
  subgraph cluster_bta {
    label = 'Bayesian trend analysis'
    style = filled
    color = orange

    node [shape = box, style = filled, fillcolor = orange, fontname = Helvetica]
    bta [label='Bayesian trend analysis']
    rates_of_change [label='Rates of change:\\n- Range shifts\\n- Depth shifts\\n- Thermal shifts']
  }

  # ---------------------------
  # BLOCK 4: Redistribution expectations
  subgraph cluster_thermal {
    label = 'Redistribution expectations based on thermal envelopes'
    style = filled
    color = lightcyan

    node [shape = box, style = filled, fillcolor = lightblue, fontname = Helvetica]
    thermal_envelopes [label='Fit thermal envelopes models']
  }

  # ---------------------------
  # BLOCK 5: Emergent patterns
  subgraph cluster_patterns {
    label = 'Emergent patterns across and within regions'
    style = filled
    color = lightpink

    node [shape = box, style = filled, fillcolor = pink, fontname = Helvetica]
    pca [label='Principal Component Analysis']
  }

  # ---------------------------
  # BLOCK 6: Linking strategies to abundance
  subgraph cluster_abundance {
    label = 'Linking strategies to population abundance'
    style = filled
    color = antiquewhite

    node [shape = box, style = filled, fillcolor = lightgreen, fontname = Helvetica]
    correlation [label='Correlation analysis']
    pop_trends [label='Population trends']
  }

  # ---------------------------
  # Arrows between nodes
  fish -> sdm [label = 'Model species distributions over time']
  temp -> grid
  depth -> grid
  thermal_envelopes -> grid [label = 'Project thermal envelopes over regular grid']

  sdm -> predict
  predict -> derived_quantities [label = 'Calculate spatial and thermal metrics']
  predict -> pop_abundance [label = 'Calculate population abundance']
  derived_quantities -> bta -> rates_of_change [label = 'Calculate rates of change']
  rates_of_change -> pca
  rates_of_change -> correlation
  pop_abundance -> pop_trends [label = 'Calculate population trends']
  pop_trends -> correlation [label = 'Strategies associated with abundance changes?']
}
")
