library(here)
library(tidyverse)
library(brms)

# read in model
m <- read_rds(here('R/bayesian_trends/fitted/m_stud.rds'))
summary(m)

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
    right = str_remove(right, "_year_c:"),
    label = formatC(Estimate, digits = 2, format = "f") %>% str_replace("0\\.", ".")
  )

# define nicer variable labels in the desired order
nice_labels <- c(
  "cogyc"          = "Latitudinal shift",
  "cogxc"          = "Longitudinal shift",
  "depthnichec"    = "Depth shift",
  "thermalnichec"  = "Thermal niche shift"
)

rho <- rho %>%
  mutate(
    left_label  = recode(as.character(left), !!!nice_labels),
    right_label = recode(as.character(right), !!!nice_labels),
    left_label  = factor(left_label, levels = nice_labels),
    right_label = factor(right_label, levels = rev(nice_labels))
  )

rho %>%
  full_join(rho,
            by = c("left_label", "right_label", "Estimate", "Est.Error", "Q2.5", "Q97.5", "label", "region", 'region_full')) %>%
  ggplot(aes(x = left_label, y = right_label)) +
  geom_tile(aes(fill = Estimate)) +
  geom_text(aes(label = sprintf("%.2f", Estimate)), size = 2) +
  scale_fill_gradient2(expression(rho),
                       low = "#59708b", mid = "#FCF9F0", high = "#A65141",
                       midpoint = 0, limits = c(-1, 1)) +
  scale_x_discrete(NULL, expand = c(0, 0)) +
  scale_y_discrete(NULL, expand = c(0, 0)) +
  facet_wrap(~ region_full, ncol = 3) +  # you can adjust ncol for layout
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(face = "bold", size = 7),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
  )

ggsave(
  filename = "output/figures/supp/corr_supp.png",
  width = 180,
  height = 160,
  dpi = 600,
  units = "mm"
)


# save rho values ---------------------------------------------------------

# rounding helper
mround <- function(x, digits) sprintf(paste0("%.", digits, "f"), round(x, digits))

# LaTeX macro writer
write_tex <- function(x, macro, append = TRUE) {
  paste0("\\newcommand{\\", macro, "}{", x, "}") |>
    write_lines("output/values/rho_values.tex", append = append)
}

# clean file before writing
unlink("output/values/rho_values.tex")

# loop through each row of rho table
for (i in seq_len(nrow(rho))) {
  
  left_var   <- rho$left[i]
  right_var  <- rho$right[i]
  region     <- rho$region_short[i]
  
  est        <- mround(rho$Estimate[i], 2)
  lower      <- mround(rho$Q2.5[i], 2)
  upper      <- mround(rho$Q97.5[i], 2)
  
  # value string
  value_str <- paste0(est, " (95\\% CI: ", lower, "--", upper, ")")
  
  # macro names (both directions)
  macro1 <- paste0("rho_", left_var, "_", right_var, "_", region)
  macro2 <- paste0("rho_", right_var, "_", left_var, "_", region)
  
  # write both macros with identical content
  write_tex(value_str, macro1)
  write_tex(value_str, macro2)
}

