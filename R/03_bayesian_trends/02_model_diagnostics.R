# libraries ---------------------------------------------------------------
library(brms)
library(here)
library(tidyverse)
library(broom)
library(knitr)
library(patchwork)
library(bayesplot)
library(posterior)

# add loo criterion -------------------------------------------------------

m_stud <- readRDS(here("R/03_bayesian_trends/fitted/m_stud.rds"))
m_mvn <- readRDS(here("R/03_bayesian_trends/fitted/m_mvn.rds"))

m_stud <- add_criterion(m_stud, c("loo"))
m_mvn <- add_criterion(m_mvn, c("loo"))

loo <- loo_compare(m_stud, m_mvn)

# print table
# Convert to data frame
loo_df <- as.data.frame(loo)
loo_df$Model <- rownames(loo_df)
rownames(loo_df) <- NULL

# make a lookup table of nicer labels
model_labels <- c(
  "m_stud" = "Student-t $+$ MVN",
  "m_mvn"  = "Gaussian $+$ MVN"
)

# map safely
loo_df <- loo_df %>%
  mutate(Model = model_labels[Model])


# Reorder and rename columns
loo_df <- loo_df[, c("Model", "elpd_diff", "se_diff")]
colnames(loo_df) <- c("Model", "ELPD difference", "SE difference")

# Replace with nicer labels
loo_df$Model <- model_labels
#loo_df$`Random effects` <- re_structures

loo_df$`ELPD difference` <- round(loo_df$`ELPD difference`, 2)
loo_df$`SE difference` <- round(loo_df$`SE difference`, 2)

# Export LaTeX table
table_elpd <- kable(
  loo_df, 
  format = "latex", 
  label = "elpd",
  booktabs = TRUE,
  escape = FALSE
  #caption = "Model comparison based on leave-one-out cross-validation (LOO).  MVN denotes a multivariate normal distribution, while Student-t $+$ MVN specifies a univariate Student-t distribution for each response combined with multivariate normal random effects. The expected log predictive density (ELPD) reflects out-of-sample predictive accuracy, with higher values indicating better performance. ELPD differences are shown relative to the best-fitting model, so negative values indicate worse fit. The standard error (SE) of the difference reflects uncertainty; differences large relative to SE suggest meaningful performance gaps."
)

writeLines(table_elpd, here('output/tables/supp/table_elpd.tex'))

# pp_checks ----------------------------------------------------------------

# extract observed data range for each variable
xlim_cogyc <- range(m_stud$data$cog_y_c, na.rm = TRUE)
xlim_cogxc <- range(m_stud$data$cog_x_c, na.rm = TRUE)
xlim_depthnichec <- range(m_stud$data$depth_niche_c, na.rm = TRUE)
xlim_thermalnichec <- range(m_stud$data$thermal_niche_c, na.rm = TRUE)

# now make posterior predictive checks using those limits
pp_cogyc <- pp_check(m_stud, resp = "cogyc", ndraws = 100,  size = 1.3) + xlim(xlim_cogyc) + ggtitle('**Latitudinal centroid (km)**')
pp_cogxc <- pp_check(m_stud, resp = "cogxc", ndraws = 100, size = 1.3) + xlim(xlim_cogxc) + ggtitle('**Longitudinal centroid (km)**') 
pp_depthnichec <- pp_check(m_stud, resp = "depthnichec", ndraws = 100,size = 1.3) + xlim(xlim_depthnichec) + ggtitle('**Depth niche (m)**')
pp_thermalnichec <- pp_check(m_stud, resp = "thermalnichec", ndraws = 100, size = 1.3) + xlim(xlim_thermalnichec) + ggtitle('**Thermal niche (\u00B0C)**')

# combine plots
pp_plot <- (pp_cogyc | pp_cogxc) / (pp_depthnichec | pp_thermalnichec) +
plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

pp_plot

ggsave(
  here('output/figures/supp/pp_supp.png'),
  #plot = pp_plot,  
  width = 180,  
  height = 140,
  dpi = 600,
  units = "mm"
)

# r hat and ESS check ------------------------------------------------------------
draws <- as_draws_df(m_stud)

summary_df <- summarise_draws(draws, mean, sd, rhat, ess_bulk, ess_tail)

# r hat
p_rhat <- ggplot(summary_df, aes(y = reorder(variable, rhat), x = rhat)) +
  # yellow shading for problematic area
  annotate("rect", xmin = 1.01, xmax = Inf, 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "yellow") +
  geom_point(color = "grey20", size = 1) +
  geom_vline(xintercept = 1.01, linetype = "dashed", color = "black", linewidth = 0.5) +
  coord_flip() +
  labs(
    x = NULL,
    y = 'Parameter',
    title = "**R-hat**"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# bulk ESS
p_bulk <- ggplot(summary_df, aes(y = reorder(variable, -ess_bulk), x = ess_bulk)) +
  annotate("rect", xmin = -Inf, xmax = 400, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "yellow") +
  geom_point(color = "grey20", size = 1) +
  geom_vline(xintercept = 400, linetype = "dashed", color = "black", linewidth = 0.5) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Parameter",
    title = "**Bulk ESS**"
  ) +
  #theme_manuscript(title=12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# tail ESS
p_tail <- ggplot(summary_df, aes(y = reorder(variable, -ess_tail), x = ess_tail)) +
  annotate("rect", xmin = -Inf, xmax = 400, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "yellow") +
  geom_point(color = "grey20", size = 1) +
  geom_vline(xintercept = 400, linetype = "dashed", color = "black", linewidth = 0.5) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Parameter",
    title = "**Tail ESS**"
  ) + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_diag <- p_rhat / p_bulk / p_tail + plot_layout(axes = "collect")

p_diag 

ggsave(
  here('output/figures/supp/rhat_supp.png'),
  #plot = pp_plot,  
  width = 180,  
  height = 190,
  dpi = 600,
  units = "mm"
)

# problem_params <- summary_df %>%
#   filter(rhat > 1.01 | ess_bulk < 400 | ess_tail < 400) %>%
#   pull(variable)
# 
# pretty_names <- c(
#   "sd_region__cogyc_year_c" = "sigma[region]^'lat centroid'",
#   "r_region__cogyc[CBS,year_c]" = "beta['region[CBS]']^'lat centroid'"
# )
# 
# # slighlty bad param
# mcmc_plot(m_stud, variable = problem_params, type = "trace") + xlab()
#   facet_wrap(.~parameter, labeller = as_labeller(pretty_names, label_parsed), scales = 'free') # oh they look good, come on!
# 
# ggsave(
#   here('output/figures/supp/bad_param_supp.png'),
#   #plot = pp_plot,  
#   width = 180,  
#   height = 80,
#   dpi = 600,
#   units = "mm"
# )



