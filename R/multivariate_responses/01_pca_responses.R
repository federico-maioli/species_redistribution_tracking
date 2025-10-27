# libraries
library(broom)
library(brms)
library(FactoMineR)
library(ggcorrplot)
library(ggrepel)
library(here)
library(paletteer)
library(patchwork)
library(tidybayes)
library(tidyverse)
library(purrr)

set.seed(123)

# load slopes -----------------------------------------------------------

slopes <- read_rds(here('R/data/processed/bayesian_species_trends.rds'))

# prepare data for PCA -----------------------------------------------------------------

pca_input <- slopes %>%
  mutate(
    pretty_outcome = case_when(
      outcome == "lat_shift"         ~ "Poleward shift",
      outcome == "lon_shift"         ~ "Eastward shift",
      outcome == "depth_shift"   ~ "Deepening",
      outcome == "thermal_shift" ~ "Warming niche"
    )
  ) %>%
  select(region, species, outcome, median_slope) %>%
  pivot_wider(names_from = outcome, values_from = median_slope) %>%
  drop_na()


# PCA by region -----------------------------------------------------------

pca_by_region <- pca_input %>%
  group_by(region) %>%
  group_nest() %>%
  mutate(
    pca = map(data, ~ PCA(select(.x, -species), graph = FALSE, ncp = 3)),
    scores = map2(pca, data, ~ .x$ind$coord %>%
                    as_tibble() %>%
                    mutate(
                      species = .y$species,
                      dist_from_origin = sqrt(Dim.1^2 + Dim.2^2)
                    )),
    loadings = map(pca, ~ .x$var$coord %>%
                     as_tibble(rownames = "outcome"))
  )



# Variance explained --------------------------------------

variance_info <- pca_by_region %>%
  mutate(var_exp = map(pca, ~ as_tibble(.x$eig) %>%
                         slice(1:2) %>%
                         mutate(PC = paste0("PC", row_number())))) %>%
  select(region, var_exp) %>%
  unnest(var_exp)

total_variance = variance_info |> filter(PC == 'PC2') |>  summarise(
  min_var = min(`cumulative percentage of variance`),
  max_var = max(`cumulative percentage of variance`)
)

variance_summary <- variance_info %>%
  group_by(PC) %>%
  summarise(
    min_var = min(`percentage of variance`),
    max_var = max(`percentage of variance`)
  )

x_lab <- paste0("PC1 (", round(variance_summary$min_var[variance_summary$PC=="PC1"], 1),
                " – ", round(variance_summary$max_var[variance_summary$PC=="PC1"], 1), "%)")
y_lab <- paste0("PC2 (", round(variance_summary$min_var[variance_summary$PC=="PC2"], 1),
                " – ", round(variance_summary$max_var[variance_summary$PC=="PC2"], 1), "%)")

# Get loadings ----------------------------------------------------------

scores_all <- bind_rows(pca_by_region$scores, .id = "region_id") %>%
  mutate(region = pca_by_region$region[as.integer(region_id)]) %>%
  select(-region_id)

# save it 
write_rds(scores_all, here('R/data/processed/pca_scores.rds'))

loadings_all <- bind_rows(pca_by_region$loadings, .id = "region_id") %>%
  mutate(
    region = pca_by_region$region[as.integer(region_id)],
    pretty_outcome = case_when(
      outcome == "lat_shift"         ~ "Poleward shift",
      outcome == "lon_shift"         ~ "Eastward shift",
      outcome == "depth_shift"   ~ "Deepening",
      outcome == "thermal_shift" ~ "Warming niche"
    )
  ) %>%
  select(-region_id)

loadings_all <-  loadings_all %>% select(outcome, Dim.1,Dim.2, region, pretty_outcome)

loadings_all <- loadings_all %>%
  mutate(color_fill = case_when(
    pretty_outcome == "Poleward shift" & Dim.1 < 0        ~ "#762A83",
    pretty_outcome == "Poleward shift" & Dim.1 >= 0       ~ "#1B7837",
    pretty_outcome == "Eastward shift" & Dim.1 < 0        ~ "#543005",
    pretty_outcome == "Eastward shift" & Dim.1 >= 0       ~ "#003C30",
    pretty_outcome == "Deepening" & Dim.1 < 0  ~ "#8C510A",
    pretty_outcome == "Deepening" & Dim.1 >= 0 ~ "#01665E",
    pretty_outcome == "Warming niche" & Dim.1 < 0 ~ "#2166AC",
    pretty_outcome == "Warming niche" & Dim.1 >= 0 ~ "#B2182B",
    TRUE ~ "#DDDDDD"  
  )) # colors

# use pretty region names
loadings_all <- loadings_all %>% mutate(

  region_full = case_when(
    region == "GOA"      ~ "Gulf of Alaska (GOA)",
    region == "BC"       ~ "British Columbia (BC)",
    region == "COW"     ~ "U.S. West Coast (USWC)",
    region == "BAL"      ~ "Baltic Sea (BAL)",
    region == "NS"       ~ "North Sea (NS)",
    region == "BS"       ~ "Barents Sea (BS)",
    region == "EBS"      ~ "East Bering Sea (EBS)",
    region == "GMX"      ~ "Gulf of Mexico (GOM)",
    region == "CBS"      ~ "Celtic–Biscay Shelf (CBS)",
    region == "NIC"      ~ "North Iberian Coast (NIC)",
    region == "NEUS"  ~ "NE U.S. & Scotian Shelf (NEUS–SS)",
    TRUE ~ region
  ),
  region_full = factor(
    region_full,
    levels = c(
      "East Bering Sea (EBS)",
      "Gulf of Alaska (GOA)",
      "British Columbia (BC)",
      "U.S. West Coast (USWC)",
      "NE U.S. & Scotian Shelf (NEUS–SS)",
      "Gulf of Mexico (GOM)",
      "Barents Sea (BS)",
      "North Sea (NS)",
      "Baltic Sea (BAL)",
      "Celtic–Biscay Shelf (CBS)",
      "North Iberian Coast (NIC)"
    )))


ggplot(loadings_all, aes(x = Dim.1, y = pretty_outcome, fill = color_fill)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = round(Dim.1, 2)),
            hjust = ifelse(loadings_all$Dim.1 >= 0, -0.1, 1.1),
            size = 2) +
  facet_wrap(~region_full, ncol = 1, scales = "free_y") +
  labs(x = "PC1 Loading", y = NULL) +
  scale_fill_identity() +
  theme_minimal(base_size = 8) +
  theme(
    #axis.text.y = element_text(size = 10),
    #axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 8)
  )

ggsave(
  here('output/figures/supp/pca_loadings_supp.png'),
  width = 180,
  height = 190,
  dpi = 600,
  units = "mm"
)


#write_csv(loadings_all, here('data/processed/pca_loadings.csv'))


# plot pca ----------------------------------------------------------------

# get posterior correlaitons

m <- read_rds(here('R/bayesian_trends/fitted/m_stud.rds'))

rho <-
  posterior_summary(m) %>% 
  data.frame() %>% 
  rownames_to_column("param") %>% 
  filter(str_detect(param, "^cor_")) %>% 
  mutate(param = str_remove(param, "^cor_region:species__")) %>% 
  separate(param, into = c("left", "right"), sep = "__") %>% 
  mutate(
    region = str_extract(left, "region[A-Z]+"),
    left   = str_remove(left, "_?region[A-Z]+"),
    right  = str_remove(right, "_?region[A-Z]+"), 
    region_short = str_remove(region, "region"),  # e.g. regionBAL → BAL
    region_full = case_when(
      region_short == "GOA"      ~ "Gulf of Alaska (GOA)",
      region_short == "BC"       ~ "British Columbia (BC)",
      region_short == "COW"     ~ "U.S. West Coast (USWC)",
      region_short == "BAL"      ~ "Baltic Sea (BAL)",
      region_short == "NS"       ~ "North Sea (NS)",
      region_short == "BS"       ~ "Barents Sea (BS)",
      region_short == "EBS"      ~ "East Bering Sea (EBS)",
      region_short == "GMX"      ~ "Gulf of Mexico (GOM)",
      region_short == "CBS"      ~ "Celtic–Biscay Shelf (CBS)",
      region_short == "NIC"      ~ "North Iberian Coast (NIC)",
      region_short == "NEUS"  ~ "NE U.S. & Scotian Shelf (NEUS–SS)",
      TRUE ~ region_short
    ),
    region_full = factor(
      region_full,
      levels = c(
        "East Bering Sea (EBS)",
        "Gulf of Alaska (GOA)",
        "British Columbia (BC)",
        "U.S. West Coast (USWC)",
        "NE U.S. & Scotian Shelf (NEUS–SS)",
        "Gulf of Mexico (GOM)",
        "Barents Sea (BS)",
        "North Sea (NS)",
        "Baltic Sea (BAL)",
        "Celtic–Biscay Shelf (CBS)",
        "North Iberian Coast (NIC)"
      ))
  )

rho <- rho %>% 
  mutate(
    left = str_remove(left, "_year_c:"),
    right = str_remove(right, "_year_c:")
  )

# ok now from here we can get the correlations we need for thermal avoiding strategies
lowest_corr <- rho %>%
  filter(left == "thermalnichec" | right == "thermalnichec") %>%
  group_by(region) %>%
  slice_min(order_by = Estimate, n = 1) %>%
  ungroup() %>% mutate(region=region_short)

# merge with pca
loadings_all <- loadings_all %>%
  left_join(lowest_corr, by = c("region","region_full"))

loadings_all <- loadings_all %>% mutate(var = case_when(
  left == 'cogyc'~"lat_shift",
  left == 'cogxc'~"lon_shift",
  left == "depthnichec"~"depth_shift",
  left == 'thermalnichec'~"thermal_shift"
))

#  now plot ---------------------------------------------------------------

base_colors <- paletteer_d("tidyquant::tq_dark", n = 12) |> as.character()

# assign them explicitly to regions
region_colors <- c(
  "GOA"     = base_colors[1],  # Gulf of Alaska
  "BC"      = base_colors[2],  # British Columbia
  "COW"    = base_colors[3],  # U.S. West Coast
  "BAL"     = base_colors[4],  # Baltic Sea
  "NS"      = base_colors[5],  # North Sea
  "BS"      = base_colors[6],  # Barents Sea
  "EBS"     = base_colors[7],  # East Bering Sea
  "GMX"     = base_colors[8],  # Gulf of Mexico
  "CBS"     = base_colors[9],  # Celtic-Biscay Shelf
  "NIC"     = base_colors[10], # North Iberian Coast
  "NEUS" = base_colors[11]  # Northeast U.S. & Scotian Shelf
)

ggplot(scores_all, aes(x = Dim.1, y = Dim.2)) +
  geom_segment(
    data = loadings_all |> filter(var != outcome),
    aes(x = 0, y = 0, xend = Dim.1 * 4, yend = Dim.2 * 4),
    arrow = arrow(length = unit(0.15, "cm")),
    color = "grey80",
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = loadings_all |> filter(var == outcome),
    aes(x = 0, y = 0, xend = Dim.1 * 4, yend = Dim.2 * 4, color = Estimate),
    arrow = arrow(length = unit(0.15, "cm")),
    inherit.aes = FALSE,
    size = 1
  ) +
  geom_point(aes(fill = region), size = 1.8, alpha = .4, shape = 21) +
  geom_text_repel(
    data = loadings_all,
    aes(x = Dim.1 * 5.2, y = Dim.2 * 5.2, label = pretty_outcome),
    size = 2.5,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.4,
    segment.color = "grey60"
  ) +
  scale_color_gradient(
    low = "#08306b",
    high = "#d0e6fa",
    name = expression(rho ~ "with warming niche")
  ) +
  guides(fill = "none") +
  scale_x_continuous(expand = expansion(mult = 0.4)) +
  facet_wrap(~region, scales = "free", ncol = 3) +
  scale_fill_manual(values = region_colors) +
  labs(x = x_lab, y = y_lab) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(0.3, "lines"),
    strip.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    legend.position = c(0.95, -0.01),  # x, y in normalized coordinates (0 = left/bottom, 1 = top/right)
    legend.justification = c("right", "bottom"),
    #legend.background = element_rect(fill = "white", color = "grey80"),
    #legend.box.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  filename = "output/figures/main/pca.png",
  width = 180,
  height = 190,
  dpi = 600,
  units = "mm"
)




# #  correlation for text ---------------------------------------------------
# 
# # Plot heatmap
# ggplot(cor_long, aes(x = variable, y = region, fill = correlation)) +
#   geom_tile(color = "white") +
#   geom_text(aes(label = round(correlation, 2)), size = 4) +
#   scale_fill_gradient2(expression(rho),
#                        low = "#59708b", mid = "#FCF9F0", high = "#A65141", midpoint = 0,
#                        labels = c(-1, "", 0, "", 1), limits = c(-1, 1)) +
#   labs(x = "", y = "Region") +
#   theme(axis.text = element_text(size = 12),
#         axis.ticks = element_blank(),
#         legend.text = element_text(hjust = 1)) + scale_x_discrete(NULL, expand = c(0, 0), labels = ggplot2:::parse_safe, position = "top") +
#   scale_y_discrete(NULL, expand = c(0, 0), labels = ggplot2:::parse_safe)
# 
# ggsave(here('plots/main/corr_outcomes.png'), width = 6, height = 4, dpi = 600)



