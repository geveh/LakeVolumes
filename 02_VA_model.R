### Fit a lake Volume-Area (V-A) hierarchical regression model.

## Load packages that we will need in this script.

# If packages haven't been installed beforehand, 
# first use install.packages("packagename")

require(tidyverse)
require(scales)
require(brms)
require(modelr)
require(tidybayes)
require(see)
require(doParallel)
require(sf)
require(pbapply)
require(ggrepel)
require(cowplot)
require(dplyr)

# Set working directory

setwd("D:/nrc_user/veh/Zhang_glacial_lakes_global/")

# Add the function HDIofMCMC from John Kruschke to this directory. We will use
# this function to derive the highest density intervals of posterior 
# (predictive) distributions of lake volumes.

source("HDIofMCMC.R")

# Useful functions

scale_this <- function(x){
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}


################################################################################

# To fit a multi-level model volume-area regression model, we load the
# text file that includes empirically measured lake areas and volumes.

vat <- read.csv2("D:/nrc_user/veh/Literature/va.txt", 
                  header = T, 
                  sep = "\t", 
                  dec = ".") %>%
  as_tibble() 

# Determine unique references.

unique(vat$Original_Reference)

# Lakes with repeated surveys.

unique(vat$Lake) %>% sort()

length(which((table(vat$Lake) %>% unname()) >1))

# We create a new column that converts the lake area from sq-meters to 
# square kilometers.

vat$Area_km2 <- vat$Area_m2 / 10^6

# If lakes were sampled repeatedly, we only take the one with the largest
# lake area, as our focus in on improving the fit for large lakes.

vat.large <- vat %>%
  group_by(Lake) %>%
  slice_max(Area_km2) %>%
  ungroup() %>%
  group_by(Lake_type) %>%
  filter(n() > 1) %>%
  ungroup()

# Plot each lake type to get an overview on the distribution of V-A data using 
# a faceted plot.

ggplot(data = vat.large,
       mapping = aes(Area_km2,
                     Bathymetric_volume_mil_m3)) +
  geom_point() +
  facet_wrap(~Lake_type, scales = "free")  +
  geom_smooth(method = "lm") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 3),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 3),
                labels = trans_format("log10", math_format(10 ^ .x)))

# Convert A and V to log10-scale and standardise them (mean of zero, unit
# standard deviation)

vat.large <- vat.large %>%
  mutate(alog = log10(Area_km2),
         vlog = log10(Bathymetric_volume_mil_m3),
         alog_scale = scale_this(alog),
         vlog_scale = scale_this(vlog)) 

# Plot again.

ggplot(data = vat.large,
       mapping = aes(alog_scale,
                     vlog_scale)) +
  geom_point() +
  facet_wrap(~Lake_type, scales = "free")  +
  geom_smooth(method = "lm") 

# Define weakly informed prior for all model parameters. They refer to data
# on the log10-scale.

prior0 <- prior(normal(0, 1.5), class = "Intercept") +
  prior(normal(1, 1.5), class = "b") +
  prior(normal(0, 0.2), class = "sd", lb = 0) +  # Tighter prior on sd
  prior(normal(0, 0.2), class = "sigma") +
  prior(normal(2, 5), class = "nu")

# Fit the Bayesian hierachical linear model that predicts log10-V 
# from log10-A, distinguished by lake type.

fit0 <- brm(bf(vlog_scale ~ alog_scale + (alog_scale | Lake_type)), 
            data = vat.large, 
            family = "student",
            prior = prior0,
            iter = 4000, 
            warmup = 1000,
            chains = 4, 
            cores = 4,
            control = list(adapt_delta = 0.99, max_treedepth = 15))

saveRDS(fit0, "VA_model.RDS")
# fit0 <- readRDS("VA_model.RDS")

# Show the summary of the model. 

fit0
fit0_plot <- plot(fit0, ask = F, nvariables = 7, theme = theme_bw(base_size =  7)) 

ggsave(filename = "model_parameters.pdf",
       fit0_plot[[1]],
       height = 200,
       width = 160,
       units = "mm")

pp_check(fit0)

################################################################################

# After fitting model, we summarise the model on two hierarchical levels:
# - the grand mean (population level);
# - the group level (individual dam types);

# Grand mean: all lake types combined.
# Define a range of lake areas (on log10-scale) for which we would like to obtain 
# draws from the posterior predictive distribution. 

alog_scale <- seq(min(vat.large$alog_scale, na.rm = T), 
                  max(vat.large$alog_scale, na.rm = T), length.out = 200)

pred.grid.grand <- add_predicted_draws(
  object = fit0, 
  newdata = data.frame(alog_scale = alog_scale),
  value = "vlog_scale", 
  re_formula = NA) %>%
  ungroup() %>%
  mutate(Area_km2 = 10^((alog_scale * sd(vat.large$alog)) + mean(vat.large$alog)),
         Bathymetric_volume_mil_m3 = 10^((vlog_scale * sd(vat.large$vlog)) + mean(vat.large$vlog)))

# Plot the trend in the grand mean of all glacier lake types.

plot.trend.grand.mean <- pred.grid.grand %>%
  ggplot(aes(x = Area_km2, y = Bathymetric_volume_mil_m3)) +
  scale_fill_manual(name = "Posterior rate", values = "grey50") +
  stat_lineribbon(aes(y = Bathymetric_volume_mil_m3), .width = 0.68,
                  point_interval = median_qi) +
  geom_jitter2(data = vat.large,
               aes(x = Area_km2, y = Bathymetric_volume_mil_m3),
               alpha = 0.8) +
  theme_bw() +
  labs(x = "Lake area [km²]",
       y = expression("Lake volume [" ~ 10^6 ~ " m"^3 * "]")) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         aspect.ratio = 1,
         legend.position = "none")  +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))

################################################################################

# Group-level effects: assess the effect (i.e. the slope) of different dam types
# in the V-A regression model.

# Transformation of parameters needs to recover the mean and standard deviation
# of log10-transformed V and A.

mu_alog <- mean(vat.large$alog)
mu_vlog <- mean(vat.large$vlog)
sd_alog <- sd(vat.large$alog)
sd_vlog <- sd(vat.large$vlog)

# Model slope (i.e. trend of V with a change in A)

# We first extract posterior draws on the population level. 

fixef_samples_b1 <- fixef(fit0, pars = "alog_scale", summary = F) %>% 
  as_tibble() 

fixef_b1_orig <- (fixef_samples_b1 * (sd_vlog / sd_alog))  

fixef_summary_b1 <- fixef_b1_orig %>%
  summarise(lo = HDIofMCMC(alog_scale, 0.68)[1],
            md = median(alog_scale),
            hi = HDIofMCMC(alog_scale, 0.68)[2]) 

# Then, we extract random effects for the posterior samples for 
# the slope parameter ("alog") and the lake type. The random effects measure
# the deviation from the population level.

ranef_samples_types_b1 <- ranef(fit0, 
                               groups = "Lake_type", 
                               pars = "alog_scale",  
                               summary = F) 

ranef_samples_types_b1 <- ranef_samples_types_b1$Lake_type[, , "alog_scale"]

# Convert random effects to a tibble.

ranef_samples_types_b1 <- ranef_samples_types_b1 %>% 
  as_tibble() 

# As the random effects measure the deviations from the grand mean, we need to 
# add them to the fixed (population-level) effects.

fixef_and_ranef_b1 <- (ranef_samples_types_b1 + fixef_samples_b1$alog_scale) %>% 
  as_tibble() * (sd_vlog / sd_alog)

# Same now for the intercept. Note that that the recovery of the
# original model parameters works a bit different to that for the slope
# https://mc-stan.org/docs/stan-users-guide/efficiency-tuning.html#standardizing-predictors

fixef_samples_b0 <- fixef(fit0, pars = "Intercept", summary = F) %>% 
  as_tibble() 

fixef_summary_b0 <- ((fixef_samples_b0$Intercept * sd_vlog) + mu_vlog - (fixef_b1_orig * mu_alog)) %>%
  as_tibble() %>%
  rename("Intercept" = "alog_scale") %>%
  summarise(lo = HDIofMCMC(Intercept, 0.68)[1],
            md = median(Intercept),
            hi = HDIofMCMC(Intercept, 0.68)[2]) 

# Then, we extract random effects for the posterior samples for 
# the slope parameter ("alog") and the lake type. The random effects measure
# the deviation from the population level.

ranef_samples_types_b0 <- ranef(fit0, 
                                groups = "Lake_type", 
                                pars = "Intercept",  
                                summary = F) 

ranef_samples_types_b0 <- ranef_samples_types_b0$Lake_type[, , "Intercept"]

# As the random effects measure the deviations from the grand mean, we need to 
# add them to the fixed (population-level) effects.

fixef_and_ranef_b0 <- (ranef_samples_types_b0 + fixef_samples_b0$Intercept) 

fixef_and_ranef_b0 <- (((fixef_and_ranef_b0 * sd_vlog) + mu_vlog) - (as.matrix(fixef_and_ranef_b1) * mu_alog)) %>%
  as_tibble()

# We summarise the trend of in the V-A relationship.

type.summary <- rbind(
  fixef_and_ranef_b0 %>% mutate(param = "b0"),
  fixef_and_ranef_b1 %>% mutate(param = "b1")) %>%
  pivot_longer(cols = unique(vat.large$Lake_type),
               names_to = "Lake_type") %>%
  group_by(param, Lake_type) %>% 
  summarise(quantslow = HDIofMCMC(value, 0.68)[1], 
            median    = median(value, 0.68),
            quantsup  = HDIofMCMC(value, 0.68)[2])

# We would like to generate draws from the expected posterior distribution.
# Therefore we generate a sequence of area values for each lake type.

pred.grid.types <- lapply(seq_range(vat.large$alog_scale, n = 101), 
                           function (x) { (unique(vat.large[,  "Lake_type"])) %>%
                               mutate(alog_scale = x) })

pred.grid.types <- do.call(rbind, pred.grid.types) %>% 
  ungroup()

# We then draw from the expected posterior predictive distribution for
# each log10-transformed lake area and lake type. 

pred.grid.types.o <- add_predicted_draws(
  object = fit0, 
  newdata = pred.grid.types,
  value = "vlog_scale", 
  ndraws = 1000,
  re_formula = ~ alog_scale + (alog_scale | Lake_type)) %>%
  mutate(Lake_type = factor(Lake_type, 
                            levels = c("supraglacial", 
                                       "ice", 
                                       "moraine", 
                                       "moraine/bedrock", 
                                       "bedrock"))) %>%
  mutate(Area_km2 = 10^((alog_scale * sd(vat.large$alog)) + mean(vat.large$alog)),
         Bathymetric_volume_mil_m3 = 10^((vlog_scale * sd(vat.large$vlog)) + mean(vat.large$vlog)))

# We plot the posterior trend in of the V-A model for each dam type 
# with the posterior trend of the population mean below and the original data as
# scattered points on top.

plot.trend.types.one.figure <- pred.grid.types.o %>%
  ggplot(aes(x = Area_km2, y = Bathymetric_volume_mil_m3/10^3, color = Lake_type)) +
  stat_lineribbon(data = pred.grid.grand,
                  aes(x = Area_km2, 
                      y = Bathymetric_volume_mil_m3/10^3), 
                  .width = 0.68,
                  fill = "grey85",
                  color = "grey50",
                  point_interval = median_hdi) +
  geom_jitter2(data = vat.large,
               aes(x = Area_km2, 
                   y = Bathymetric_volume_mil_m3/10^3, 
                   colour = Lake_type),
               alpha = 0.8) +
  stat_lineribbon(aes(x = Area_km2, 
                      y = Bathymetric_volume_mil_m3/10^3), 
                  .width = 0,
                  linewidth = 1,
                  alpha = 0.8,
                  point_interval = median_qi) +
  scale_fill_viridis_d(option = "plasma") +
  scale_color_viridis_d(name = "Dam type", option = "plasma") +
  theme_bw() +
  labs(x = "Lake area [km²]",
       y = "Lake volume [km³]") +
  theme(text = element_text(size = 7),
        legend.position.inside = c(0.2, 0.7),
       # aspect.ratio = 1,
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank())  + 
  guides(fill   = "none",
         colour = guide_legend(position = "inside")) +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))

ggsave(filename = "va_model_one_panel.pdf",
       plot = plot.trend.types.one.figure, 
       width = 140,
       height = 70,
       units = "mm",
       device = cairo_pdf, 
       family = "Arial")

# Same plot as above, but now we plot all dam types as separate panels, including
# the posterior regression estimate for that dam type.

lab <- vat.large %>% 
  group_by(Lake_type) %>% 
  summarise(n = n(),
            lab = paste0("n = ", n)) %>%
  mutate(Lake_type = factor(Lake_type, 
                            levels = c("supraglacial", 
                                       "ice", 
                                       "moraine", 
                                       "moraine/bedrock", 
                                       "bedrock"))) 

# This plot adds facets.

plot.trend.types <- pred.grid.types.o %>%
  mutate(Lake_type = fct_relevel(Lake_type, 
                                 "supraglacial", 
                                 "ice", 
                                 "moraine", 
                                 "moraine/bedrock", 
                                 "bedrock")) %>%
  ggplot(aes(x = Area_km2, y = Bathymetric_volume_mil_m3, fill = Lake_type)) +
  stat_lineribbon(data = pred.grid.grand,
                  aes(x = Area_km2,
                      y = Bathymetric_volume_mil_m3),
                  .width = 0.68,
                  fill = "grey90",
                  color = "grey30",
                  point_interval = median_hdi,
                  linewidth = 1) +
  stat_lineribbon(aes(x = Area_km2, 
                      y = Bathymetric_volume_mil_m3), 
                  .width = 0.68,
                  color = "black",
                  alpha = 0.75,
                  point_interval = median_hdi,
                  linewidth = 1) +
  scale_fill_viridis_d(option = "plasma") +
  geom_text(data = lab,
            aes(x = 10^-4, y = 10^5, label = lab),
            inherit.aes = FALSE,    # use annotation data only
            color = "black",
            fontface = "italic",
            size = 7 * 0.35,
            hjust = 0) +
  facet_wrap(~Lake_type, nrow = 2) +
  geom_jitter2(data = vat.large,
               aes(x = Area_km2, y = Bathymetric_volume_mil_m3),
               alpha = 0.5,
               color = "black") +
  theme_bw() +
  labs(x = "Lake area [km²]",
       y = expression("Lake volume [" ~ 10^6 ~ " m"^3 * "]")) +
  theme(text = element_text(size = 7),
        strip.background = element_blank(),
        legend.position = "none", 
        aspect.ratio = 1)  + 
  scale_y_log10(labels = function(x) {
    sapply(x, function(val) {
      if (is.na(val)) {
        NA 
      } else if (val > 1) {
        format(round(val, 0), scientific = FALSE)
      } else {
        format(val, digits = 2)
      }
    })
  }) +
  scale_x_log10(labels = function(x) {
    sapply(x, function(val) {
      if (is.na(val)) {
        NA  
      } else if (val > 1) {
        format(round(val, 0), scientific = FALSE)
      } else {
        format(val, digits = 2)
      }
    })
  })

# We draw the the posterior slopes for the population-level (all dam types 
# combined) as vertical lines, and the regional model intercepts an slopes as 
# distributions.

vlines <- bind_rows(fixef_summary_b0, fixef_summary_b1) %>%
  mutate(param = c("beta[0]", "beta[1]"))

plot.post.types <- rbind(
  fixef_and_ranef_b0 %>% mutate(param = "b0"),
  fixef_and_ranef_b1 %>% mutate(param = "b1")) %>%
  mutate(param = recode(param, 
                        "b0" = "beta[0]", 
                        "b1" = "beta[1]")) %>%
  pivot_longer(cols = unique(vat.large$Lake_type),
               names_to = "Lake_type") %>%
  mutate(Lake_type = factor(Lake_type, 
                            levels = c("supraglacial", 
                                       "ice", 
                                       "moraine", 
                                       "moraine/bedrock", 
                                       "bedrock"))) %>%
  group_by(param, Lake_type) %>%
  filter(value > HDIofMCMC(value, 0.99)[1],
         value < HDIofMCMC(value, 0.99)[2]) %>%
  mutate(md_grp_alog = median(value)) %>% 
  ggplot(aes(x = value, 
             y = Lake_type, 
             fill = Lake_type)) +
  geom_vline(data = vlines,
             aes(xintercept = lo), 
             color = "darkgrey", linetype = "dashed") +
  geom_vline(data = vlines,
             aes(xintercept = hi), 
             color = "darkgrey", linetype = "dashed") +
  geom_vline(data = vlines,
             aes(xintercept = md), 
             color = "darkgrey") +
  stat_halfeye(.width = 0.68,
               interval_size = 0.5, 
               shape = 21,
               point_color = "red",
               point_fill = "white",
               point_size = 1.5,
               show.legend = FALSE) +
  scale_fill_viridis_d(option = "plasma") +
  labs(x = expression("Intercept" ~ beta[0] ~ "and slope" ~ beta[1] ~ 
                        "for" ~ log[10] ~ 
                        "-transformed lake areas and volumes"), 
       y = "Dam type") +
  facet_wrap(~param, scales = "free_x", labeller = label_parsed) +
  theme_bw() +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        text = element_text(size = 7))

# We combine plots showing the data with trends, and the posterior slope 
# distribution.

plot.types.post <- plot_grid(plot.trend.types, 
                            plot.post.types, 
                            nrow = 2,
                            axis = "l",
                            rel_heights = c(2,1),
                            labels = c('a', 'b'),
                            align = "v",
                            label_size = 8)

ggsave(filename = "Trend_types_and_posterior.pdf",
       plot = plot.types.post, 
       height = 170,
       width = 150,
       units = "mm")


################################################################################

# Predict the posterior distribution of lake volume for all in-situ surveyed
# glacial lakes. We use the 68% HDI (approx. 1 standard deviation in frequentist
# statistics).

preds <- posterior_predict(fit0,
                           ndraws = 1000,
                           newdata = vat.large,
                           re_formula = NULL)

preds <- ((preds * sd(log10(vat.large$Bathymetric_volume_mil_m3))) + mean(log10(vat.large$Bathymetric_volume_mil_m3))) 

preds <- tibble(lo = apply(preds, 2, HDIofMCMC, 0.68)[1, ],
                md = apply(preds, 2, median),
                hi = apply(preds, 2, HDIofMCMC, 0.68)[2, ]) 

vat.pred <- bind_cols(vat.large, preds) 

# Volume stored in Lake Hazen, the largest lake in our sample.

hazen <- vat.pred %>%
  filter(Lake == "Lake Hazen") %>%
  summarise(md = md,
            lo = lo-md,
            hi = hi-md)

################################################################################

# Error analysis

# On log-log-transformed scale.
# We compare median errors above and below suspected breakpoints in the V-A distribution,
# 0.5 km² or 5 km², as reported by Shugar et al. (2020) and Zhang et al. (2024)

m.res.gt05 <- bind_cols(vat.large, preds) %>%
  mutate(resid = md-vlog) %>%
  group_by(Lake_type) %>%
  filter(Area_km2 > 0.5) %>%
  summarise(m_res = median(resid))

m.res.sm05 <- bind_cols(vat.large, preds) %>%
  mutate(resid = md-vlog) %>%
  group_by(Lake_type) %>%
  filter(Area_km2 < 0.5) %>%
  summarise(m_res = median(resid))

m.res.gt5 <- bind_cols(vat.large, preds) %>%
  mutate(resid = md-vlog) %>%
  group_by(Lake_type) %>%
  filter(Area_km2 > 5) %>%
  summarise(m_res = median(resid))

m.res.sm5 <- bind_cols(vat.large, preds) %>%
  mutate(resid = md-vlog) %>%
  group_by(Lake_type) %>%
  filter(Area_km2 < 5) %>%
  summarise(m_res = median(resid))

# Plot the absolute errors (median predicted vs. reported V)

plot.abs.errors.log <- bind_cols(vat.large, preds) %>%
  ggplot(mapping = aes(x = Area_km2,
                       y = md-vlog)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "tomato") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "deepskyblue") +
  geom_point(alpha = 0.66) +
  geom_linerange(aes(ymin = lo-vlog,
                     ymax = hi-vlog),
                 alpha = 0.33) +
  geom_segment(data = m.res.sm05,
               mapping = aes(x = min(vat.large$Area_km2, na.rm = T),
                             y = m_res,
                             xend = 0.5,
                             yend = m_res),
               color = "tomato1",
               linewidth = 1) +
  geom_segment(data = m.res.gt05,
               mapping = aes(x = 0.5,
                             y = m_res,
                             xend = max(vat.large$Area_km2, na.rm = T),
                             yend = m_res),
               color = "deepskyblue1",
               linewidth = 1)+
  geom_segment(data = m.res.sm5,
               mapping = aes(x = min(vat.large$Area_km2, na.rm = T),
                             y = m_res,
                             xend = 5,
                             yend = m_res),
               color = "tomato4",
               linewidth = 1) +
  geom_segment(data = m.res.gt5,
               mapping = aes(x = 5,
                             y = m_res,
                             xend = max(vat.large$Area_km2, na.rm = T),
                             yend = m_res),
               color = "deepskyblue4",
               linewidth = 1.2)+
  facet_wrap(~Lake_type, scales = "free_y") +
  theme_bw() +
  labs(x = "Lake area [km²]",
       y = "Absolute errors of the exponents\nbetween predicted and observed lake volumes") +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         strip.background = element_blank(),
         legend.position = "none", 
         aspect.ratio = 0.5)  +
  scale_x_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))

################################################################################

# Errors back-transformed on original scale

# Same as above: Absolute error for predicted volumes of lakes greater than 
# 0.5 or 5 km².

m.res.gt05 <- bind_cols(vat.large, preds) %>%
  mutate(resid = 10^md-10^vlog) %>%
  group_by(Lake_type) %>%
  filter(Area_km2 > 0.5) %>%
  summarise(m_res = median(resid))

m.res.sm05 <- bind_cols(vat.large, preds) %>%
  mutate(resid = 10^md-10^vlog) %>%
  group_by(Lake_type) %>%
  filter(Area_km2 < 0.5) %>%
  summarise(m_res = median(resid))

m.res.gt5 <- bind_cols(vat.large, preds) %>%
  mutate(resid = 10^md-10^vlog) %>%
  group_by(Lake_type) %>%
  filter(Area_km2 > 5) %>%
  summarise(m_res = median(resid))

m.res.sm5 <- bind_cols(vat.large, preds) %>%
  mutate(resid = 10^md-10^vlog) %>%
  group_by(Lake_type) %>%
  filter(Area_km2 < 5) %>%
  summarise(m_res = median(resid))

# Plot absolute errors on original scale.

plot.abs.errors.orig <- bind_cols(vat.large, preds) %>%
  ggplot(mapping = aes(x = 10^alog,
                       y = 10^md-10^vlog)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "tomato") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "deepskyblue") +
  geom_point(alpha = 0.66) +
  geom_segment(data = m.res.sm05,
               mapping = aes(x = min(vat.large$Area_km2, na.rm = T),
                             y = m_res,
                             xend = 0.5,
                             yend = m_res),
               color = "tomato1",
               linewidth = 1) +
  geom_segment(data = m.res.gt05,
               mapping = aes(x = 0.5,
                             y = m_res,
                             xend = max(vat.large$Area_km2, na.rm = T),
                             yend = m_res),
               color = "deepskyblue1",
               linewidth = 1)+
  geom_segment(data = m.res.sm5,
               mapping = aes(x = min(vat.large$Area_km2, na.rm = T),
                             y = m_res,
                             xend = 5,
                             yend = m_res),
               color = "tomato4",
               linewidth = 1) +
  geom_segment(data = m.res.gt5,
               mapping = aes(x = 5,
                             y = m_res,
                             xend = max(vat.large$Area_km2, na.rm = T),
                             yend = m_res),
               color = "deepskyblue4",
               linewidth = 1)+
  facet_wrap(~Lake_type, scales = "free_y") +
  theme_bw() +
  labs(x = "Lake area [km²]",
       y = expression(atop("Absolute errors [10"^6 * " m"^3 * "]",
                           "predicted vs. observed lake volumes"))) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         strip.background = element_blank(),
         legend.position = "none", 
         aspect.ratio = 0.5) +
  scale_x_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))


# Show relative errors, i.e. as percentage of under- or overestimation,
# again for data smaller or greater than 0.5 or 5 km².

m.res.rel.gt05 <- bind_cols(vat.large, preds) %>%
  mutate(resid = (10^md - 10^vlog)/(10^vlog)) %>%
  group_by(Lake_type) %>%
  filter(Area_km2 > 0.5) %>%
  summarise(m_res = median(resid))

m.res.rel.sm05 <- bind_cols(vat.large, preds) %>%
  mutate(resid = (10^md - 10^vlog)/(10^vlog)) %>%
  group_by(Lake_type) %>%
  filter(Area_km2 < 0.5) %>%
  summarise(m_res = median(resid))


m.res.rel.gt5 <- bind_cols(vat.large, preds) %>%
  mutate(resid = (10^md - 10^vlog)/(10^vlog)) %>%
  group_by(Lake_type) %>%
  filter(Area_km2 > 5) %>%
  summarise(m_res = median(resid))

m.res.rel.sm5 <- bind_cols(vat.large, preds) %>%
  mutate(resid = (10^md - 10^vlog)/(10^vlog)) %>%
  group_by(Lake_type) %>%
  filter(Area_km2 < 5) %>%
  summarise(m_res = median(resid))

plot.rel.errors <- bind_cols(vat.large, preds) %>%
  ggplot(mapping = aes(x = 10^alog,
                       y = (10^md - 10^vlog)/(10^vlog))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "tomato") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "deepskyblue") +
    geom_point(alpha = 0.66) +
  geom_segment(data = m.res.rel.sm05,
               mapping = aes(x = min(vat.large$Area_km2, na.rm = T),
                             y = m_res,
                             xend = 0.5,
                             yend = m_res),
               color = "tomato1",
               linewidth = 1) +
  geom_segment(data = m.res.rel.gt05,
               mapping = aes(x = 0.5,
                             y = m_res,
                             xend = max(vat.large$Area_km2, na.rm = T),
                             yend = m_res),
               color = "deepskyblue1",
               linewidth = 1)+
  geom_segment(data = m.res.rel.sm5,
               mapping = aes(x = min(vat.large$Area_km2, na.rm = T),
                             y = m_res,
                             xend = 5,
                             yend = m_res),
               color = "tomato4",
               linewidth = 1) +
  geom_segment(data = m.res.rel.gt5,
               mapping = aes(x = 5,
                             y = m_res,
                             xend = max(vat.large$Area_km2, na.rm = T),
                             yend = m_res),
               color = "deepskyblue4",
               linewidth = 1)+
  facet_wrap(~Lake_type, scales = "free_y") +
  theme_bw() +
  labs(x = "Lake area [km²]",
       y = "Relative errors (ratio between\npredicted and observed lake volumes)") +
  theme(text = element_text(size = 7),
        strip.background = element_blank(),
        legend.position = "none", 
        aspect.ratio = 0.5) +
  scale_x_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))

plot.errors <- plot_grid(plot.abs.errors.log,
                         plot.abs.errors.orig, 
                         plot.rel.errors, 
                         nrow = 3,
                         align = "v",
                         axis = "l",
                         labels = c('a', 'b', 'c'),
                         label_size = 8)

ggsave(filename = "VA_model_errors.pdf",
       plot = plot.errors, 
       height = 200,
       width = 180,
       units = "mm")

saveRDS(vat.large, "VA_data.RDS")
