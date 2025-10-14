source("renv/activate.R")

# Load necessary packages
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
library(ggplot2)

if (!requireNamespace("ggtext", quietly = TRUE)) stop("Please install ggtext")
library(ggtext)

# Define the custom ggplot theme
theme_manuscript <- function(base_size = 12, title = 10) {
  ggplot2::theme_light(base_size = base_size) +
    ggplot2::theme(
      #legend.position = "none",
      plot.title = ggtext::element_markdown(hjust = 0.5, size = title),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) 
}

# Set it globally for all ggplots in this session
ggplot2::theme_set(theme_manuscript())

message("theme_manuscript() set as default.")
