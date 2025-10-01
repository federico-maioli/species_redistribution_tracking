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

slopes <- read_rds(here('data/processed/bayesian_species_trends.rds'))

# prepare data for PCA -----------------------------------------------------------------

pca_input <- slopes %>%
  mutate(
    pretty_outcome = case_when(
      outcome == "cogyc"         ~ "Poleward shift",
      outcome == "cogxc"         ~ "Eastward shift",
      outcome == "depthnichec"   ~ "Deepening",
      outcome == "thermalnichec" ~ "Warming niche"
    )
  ) %>%
  select(region, species, outcome, mean_slope) %>%
  pivot_wider(names_from = outcome, values_from = mean_slope) %>%
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
write_rds(scores_all, here('data/processed/pca_scores.rds'))

loadings_all <- bind_rows(pca_by_region$loadings, .id = "region_id") %>%
  mutate(
    region = pca_by_region$region[as.integer(region_id)],
    pretty_outcome = case_when(
      outcome == "cogyc"         ~ "Poleward shift",
      outcome == "cogxc"         ~ "Eastward shift",
      outcome == "depthnichec"   ~ "Deepening",
      outcome == "thermalnichec" ~ "Warming niche"
    )
  ) %>%
  select(-region_id)

loadings_all <-  loadings_all %>% select(outcome, Dim.1, region, pretty_outcome)

ggplot(loadings_all, aes(x = Dim.1, y = pretty_outcome, fill = pretty_outcome)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = round(Dim.1, 2)),
            hjust = ifelse(loadings_all$Dim.1 >= 0, -0.1, 1.1),
            size = 3) +
  facet_wrap(~region, ncol = 1, scales = "free_y") +
  labs(x = "PC1 Loading", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

ggsave(here('multivariate_responses/figures/loadings.png'),bg = 'white', width = 8, height = 12, dpi = 600)

write_csv(loadings_all, here('data/processed/pca_loadings.csv'))

# Correlation and dominant strategies -------------------------------------

wide_df <- slopes %>%
  select(species, region, outcome, mean_slope) %>%
  pivot_wider(names_from = outcome, values_from = mean_slope)

cor_by_region <- wide_df %>%
  group_by(region) %>%
  summarise(
    cor_lat = cor(cogyc, thermalnichec, use = "complete.obs"),
    cor_lon = cor(cogxc, thermalnichec, use = "complete.obs"),
    cor_depth = cor(depthnichec, thermalnichec, use = "complete.obs")
  ) # get correlation with thermal niche by region and outcome

# NEUS correlation
neus_corr <- wide_df %>% filter(region == 'NEUS') %>% group_by(region) |> 
  summarise(
    cor_lat_lon = cor(cogyc, cogxc, use = "complete.obs")
  )

# NS correlation
ns_corr <- wide_df %>% filter(region == 'NS') %>% group_by(region) |> 
  summarise(
    cor_lat_depth = cor(cogyc, depthnichec, use = "complete.obs")
  )

# BAL correlation
bal_corr <- wide_df %>% filter(region == 'BAL') %>% group_by(region) |>
  summarise(cor_lat_depth = cor(cogxc, depthnichec, use = "complete.obs"))

# get the dominant strategy for warming less
dominant_strategy <- cor_by_region %>%
  pivot_longer(cols = starts_with("cor_"),
               names_to = "strategy",
               values_to = "correlation") %>%
  mutate(strategy = recode(strategy,
                           cor_lat = "Poleward shift",
                           cor_lon = "Eastward shift",
                           cor_depth = "Deepening")) %>%
  group_by(region) %>%
  slice_min(order_by = correlation, n = 1, with_ties = FALSE) %>%
  ungroup()

loadings_all <- loadings_all %>%
  left_join(dominant_strategy, by = "region")

# Plot ---------------------------------------------------------------------

ggplot(scores_all, aes(x = Dim.1, y = Dim.2)) +
  geom_segment(
    data = loadings_all |> filter(pretty_outcome != strategy),
    aes(x = 0, y = 0, xend = Dim.1 * 4, yend = Dim.2 * 4),
    arrow = arrow(length = unit(0.15, "cm")),
    color = "grey80",
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = loadings_all |> filter(pretty_outcome == strategy),
    aes(x = 0, y = 0, xend = Dim.1 * 4, yend = Dim.2 * 4, color = correlation),
    arrow = arrow(length = unit(0.15, "cm")),
    inherit.aes = FALSE,
    size = 1
  ) +
  geom_point(aes(fill = region), size = 1.8, alpha = .8, color='black', shape = 21) +
  geom_text(
    data = loadings_all,
    aes(x = Dim.1 * 5, y = Dim.2 * 5, label = pretty_outcome),
    size = 3, fontface = "bold"
  ) +
  scale_color_gradient(
    low = "#08306b",    
    high = "#d0e6fa", 
    name =  expression(rho ~ "with warming niche")
  ) +
  guides(fill = "none") +
  scale_x_continuous(expand = expansion(mult = 0.4)) +
  facet_wrap(~region, scales = "free", ncol = 3) +
  scale_fill_paletteer_d("yarrr::basel") +
  labs(x = x_lab, y = y_lab) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(0.3, "lines"),
    strip.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(here('multivariate_responses/figures/pca.png'), width = 12, height = 9, dpi = 600)


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



