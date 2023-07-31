
# Setup ------------------------------------------------------------------------

library(tidyverse)
library(ggpubr)
theme_set(theme_bw())

simres <- readRDS("simulation/psyc520/simres.rds")

simres_long <- simres %>%
  tibble::rowid_to_column("cond") %>%
  pivot_longer(
    bias.est_PI_meani:coverage.ci_FSR_vars, 
    names_to = c("criteria", "method", "stat", "par"), 
    values_to = "value", 
    names_pattern = "(.*)_(.*)_(.*)(.)"
  ) |>
  mutate(method = recode(method, "PI" = "JSEM", "TSPAR" = "'2S-PA'", 
                         "FSR" = "'FS-PA'"), 
         method = factor(method,
                         levels = c("JSEM", "'2S-PA'", "'FS-PA'", "TSPAB", "FSB")),
         kappa2 = paste0("kappa[2]==", kappa2)) |>
  filter(!(method %in% c("TSPAB", "FSB")))

# Bias -------------------------------------------------------------------------

bias_meani <- simres_long %>%
  filter(criteria == "bias.est", stat == "mean", 
         par == "i") %>%
  mutate(n = as.factor(n)) %>%
  ggplot(aes(x = as.factor(r_ni), y = value, 
             color = n, shape = n)) +
  geom_point(aes(col = n), position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .15) +
  facet_grid(kappa2~ method, labeller = label_parsed) +
  ylim(-.01, .01) +
  labs(x = expression(r[ni]/p[ni]), 
       y = "Bias", title = "Mean of Intercept", 
       color = "Sample Size", shape = "Sample Size")

bias_means <- simres_long %>%
  filter(criteria == "bias.est", stat == "mean", 
         par == "s") %>%
  mutate(n = as.factor(n)) %>%
  ggplot(aes(x = as.factor(r_ni), y = value, 
             color = n, shape = n)) +
  geom_point(aes(col = n), position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .15) +
  facet_grid(kappa2~ method, labeller = label_parsed) +
  ylim(-.01, .01) +
  labs(x = expression(r[ni]/p[ni]), 
       y = "Bias", title = "Mean of Slope", 
       color = "Sample Size", shape = "Sample Size")

bias_vari <- simres_long %>%
  filter(criteria == "bias.est", stat == "var", 
         par == "i") %>%
  mutate(n = as.factor(n)) %>%
  ggplot(aes(x = as.factor(r_ni), y = value, 
             color = n, shape = n)) +
  geom_point(aes(col = n), position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .15) +
  facet_grid(kappa2~ method, labeller = label_parsed) +
  ylim(-.07, .07) +
  labs(x = expression(r[ni]/p[ni]), 
       y = "Bias", title = "Variance of Intercept", 
       color = "Sample Size", shape = "Sample Size")

bias_vars <- simres_long %>%
  filter(criteria == "bias.est", stat == "var", 
         par == "s") %>%
  mutate(n = as.factor(n)) %>%
  ggplot(aes(x = as.factor(r_ni), y = value, 
             color = n, shape = n)) +
  geom_point(aes(col = n), position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .15) +
  facet_grid(kappa2~ method, labeller = label_parsed) +
  ylim(-.07, .07) +
  labs(x = expression(r[ni]/p[ni]), 
       y = "Bias", title = "Variance of Slope", 
       color = "Sample Size", shape = "Sample Size")

bias_plots <- ggarrange(bias_meani, bias_means, bias_vari, bias_vars, 
                        nrow = 2, ncol = 2, common.legend = TRUE, 
                        legend = "bottom")

ggsave(here::here("manuscript/psychometrika/proposal/figures/bias_plots.png"), 
       bias_plots, width = 3000, height = 2160, units = "px")
