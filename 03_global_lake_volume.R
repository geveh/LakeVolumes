### Predict glacier lake volumes globally.

## Load packages that we'll need in the course of this script

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
require(sf)
require(lwgeom)
require(dplyr)

# Set working directory

setwd("D:/nrc_user/veh/Zhang_glacial_lakes_global/")

# Add the function HDIofMCMC from John Kruschke to this directory. We will use
# this function to derive the highest density intervals of posterior 
# (predictive) distributions of lake volumes.

source("HDIofMCMC.R")

# We also write a function to remove extreme values from the posterior predictions
# to reduce the influence of a few outliers on the posterior predicitive 
# distribution.

replace_extremes <- function(row) {
  
  HDIrange <- HDIofMCMC(row, credMass = 0.98)
  
  lower <- HDIrange[1]  # 1st percentile
  upper <-  HDIrange[2]  # 99th percentile
  
  # Identify extreme values
  extreme_idx <- which(row < lower | row > upper)
  
  if (length(extreme_idx) > 0) {
    # Replace with random values from the same row (excluding extremes)
    row[extreme_idx] <- sample(row[row >= lower & row <= upper], 
                               length(extreme_idx), 
                               replace = TRUE)
  }
  
  return(row)
}

# We need the glacial lake areas mapped in 1990 and 2020 to predict the local, 
# regional, and global glacial lake volumes. We also load the hierachical V-A
# model and the object containing the empirically surveyed V-A data.

all.lakes <- readRDS("Lakes19902020_damtype.RDS")
fit0 <- readRDS("VA_model.RDS")
vat.large <- readRDS("VA_data.RDS")

# How many lakes do we have for each year?

all.lakes %>%
  st_drop_geometry() %>%
  group_by(Year_new) %>%
  summarise(n())

################################################################################

# Obtain the posterior predictive distribution for every lake globally.
# We obtain 1000 posterior predictions for every lake in our sample

all.lakes <- all.lakes %>%
  mutate(alog_scale = (alog - mean(vat.large$alog)) / sd(vat.large$alog) )

pred.vol.all <- posterior_predict(fit0, 
                                  all.lakes,
                                  ndraws = 1000,
                                  re_formula = ~ alog_scale + (alog_scale | Lake_type)) %>%
  t()

# Apply function to remove outliers to each row of the matrix

pred.vol.all <- t(apply(pred.vol.all, 1, replace_extremes))

# Convert data to original scale

pred.vol.all <- 10^((pred.vol.all  * sd(vat.large$vlog)) + mean(vat.large$vlog))

colnames(pred.vol.all) <- paste0("Pred_", seq_len(ncol(pred.vol.all)))

pred.vol.all <- pred.vol.all %>%
  as_tibble()

# For each lake, obtain a summary of the posterior predictive distribution,
# including the median and the 68% HDI.

vol.per.lake <- apply(pred.vol.all, MARGIN = 1, function (x) {
  
  x.orig <- as.numeric(x)
  
  HDI <- HDIofMCMC(x.orig, credMass = 0.68)
  md <- median(x.orig)
  
  return(tibble(lo = HDI[1],
                md = md,
                hi = HDI[2])) } )


vol.per.lake <- bind_rows(vol.per.lake)

# Add this posterior estimate to the lake.

all.lakes <- bind_cols(all.lakes, pred.vol.all)

saveRDS(all.lakes, "all_lakes_with_volumes.RDS")
# all.lakes <- readRDS("all_lakes_with_volumes.RDS")

# Summary per lake

all.lakes.summary <- readRDS("Lakes19902020_damtype.RDS") %>%
  bind_cols(., vol.per.lake)

# Find lakes with median volumes larger than 1 km³.

all.lakes.summary %>%
  filter(Year_new == 2020,
         md > 1000)

# Individual lakes and their size distribution in 2020.

o1regions <- st_read("D:/nrc_user/veh/RGI2000-v7.0-regions/RGI2000-v7.0-o1regions.shp")
o1regions <- o1regions[!duplicated(o1regions$o1region), ]

# Calculate empirical exceedance probabilities of lake volumes per region.
# To this end, we also need to add the regional identifier of the RGI to 
# each lake.

all.lakes.sort <- all.lakes %>%
  bind_cols(., vol.per.lake) %>%
  st_drop_geometry() %>%
  filter(Year_new == 2020) %>%
  dplyr::left_join(., o1regions %>% st_drop_geometry(), 
                   by = "o1region", 
                   relationship = "many-to-one") %>%
  group_by(full_name.y)

# Sort lakes in descending order (i.e. from high to low exceedance probability)
# for the median predicted lake volume for each lake.

med.vol.sort <- all.lakes.sort %>% 
  dplyr::mutate(md = md *10^6) %>%
  dplyr::arrange(md, .by_group = TRUE) %>% 
  dplyr::mutate(p = seq(n(), 1, -1) / n())

# Same for the lower 68% HDI estimate for each lake...

lo.vol.sort <- all.lakes.sort %>% 
  dplyr::mutate(lo = lo *10^6) %>%
  dplyr::arrange(lo*10^6, .by_group = TRUE) %>% 
  dplyr::mutate(p = seq(n(), 1, -1) / n())

# ...  and the upper 68% HDI estimate for each lake.

hi.vol.sort <- all.lakes.sort %>% 
  dplyr::mutate(hi = hi *10^6) %>%
  dplyr::arrange(hi*10^6, .by_group = TRUE) %>% 
  dplyr::mutate(p = seq(n(), 1, -1) / n())

# We show the regional exceedance probabilities in a plot, and add a label that 
# describes how many lakes have volumes smaller than 1M m³ or larger than 1 km³.

label.vol.sort <- med.vol.sort %>%
  dplyr::summarise(n_sm1mil = sum(md <10^6),
                   p_sm1mil = round((n_sm1mil/ n())*100, digits = 1),
                   n_gt1km = sum(md >10^9),
                   p_gt1km = round((n_gt1km/ n())*100, digits = 1)) 

plot.regional.size.dist <- ggplot(med.vol.sort,
       mapping = aes(x = md,
                     y = p)) +
  geom_line(data = lo.vol.sort,
            mapping = aes(x = lo,
                          y = p),
            inherit.aes = F,
            color = "grey50",
            linewidth = 0.5) +
  geom_line(data = hi.vol.sort,
            mapping = aes(x = hi,
                          y = p),
            inherit.aes = F,
            color = "grey50",
            linewidth = 0.5) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 10^6,
             linetype = "dashed",
             color = "darkorange",
             linewidth = 0.4) +
  geom_text(data = label.vol.sort,
            aes(x = 10^2, y = 10^-3, 
               label = n_sm1mil),
            color = "darkorange2",
            size = 7 * 0.35,
            hjust = 0) +
  geom_text(data = label.vol.sort,
            aes(x = 10^2, y = 10^-3.5, 
                label = paste0(p_sm1mil, "%")),
            color = "darkorange2",
            size = 7 * 0.35,
            hjust = 0) +
  geom_vline(xintercept = 10^9,
             linetype = "dashed",
             color =  "navy",
             linewidth = 0.4) +
  geom_text(data = label.vol.sort,
            aes(x = 10^9.5, y = 10^-0.5, 
                label = n_gt1km),
            color = "navy",
            size = 7 * 0.35,
            hjust = 0) +
  geom_text(data = label.vol.sort,
            aes(x = 10^9.5, y = 10^-1, 
                label = paste0(p_gt1km, "%")),
            color = "navy",
            size = 7 * 0.35,
            hjust = 0) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 6),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 6),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  facet_wrap(~full_name.y, ncol = 4) +
  theme_bw() +
  labs(x = "Estimated lake volume [m³]",
       y = "Empirical exceedance probability") + 
  theme(text = element_text(size = 7),
        aspect.ratio = 1,
        panel.background = element_blank(),
        strip.background = element_blank())
  
# Save this plot with the exceedance probabilities to disk (Figure S9).

ggsave("Regional_size_distribution.pdf",
       plot.regional.size.dist,
       width = 180,
       height = 190,
       units = "mm",
       device = cairo_pdf, 
       family = "Arial")

# Assessing individual changes in glacial lake volumes.

# We remove supraglacial lakes as those are too volatile through time (they can
# move as the glacier flows down-valley) to be robustly detected between 1990 
# and 2020. 

all.lakes.summary.no.sup <- all.lakes.summary %>%
  filter(Lake_type != "supraglacial")

# Estimate the global glacial lake volume in 1990 and 2020. All posterior 
# predictions from all lakes contribute equally to this estimate.

volume.all.lakes <- all.lakes %>% 
  st_drop_geometry() %>% 
  dplyr::select("Year_new", "o1region", starts_with("Pred_")) %>%
  group_by(Year_new) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE))) %>%
  group_by(Year_new) %>%
  summarise(lo = HDIofMCMC(c_across(starts_with("Pred_")), 0.68)[1],
            md = quantile(c_across(starts_with("Pred_")), 0.5),
            hi = HDIofMCMC(c_across(starts_with("Pred_")), 0.68)[2])

# Global lake volume plus/ minus median.

volume.all.lakes  %>% 
  group_by(Year_new) %>%
  summarise(median = md/1000,
            qminus = (md-lo) /1000,
            qplus = (hi-md)/1000)

# Median increase in global lake volume.

(volume.all.lakes %>% filter(Year_new == 2020) %>% pull(md) -
volume.all.lakes %>% filter(Year_new == 1990) %>% pull(md)) /
volume.all.lakes %>% filter(Year_new == 1990) %>% pull(md)

# Volume regional in 1990 and 2020. Summary statistics are aggregated on
# regional level.

volume.regional.lakes <- all.lakes %>% 
  st_drop_geometry() %>% 
  dplyr::select("Year_new", "o1region", starts_with("Pred_")) %>%
  group_by(Year_new, o1region) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = starts_with("Pred_")) %>% 
  group_by(Year_new, o1region) %>%
  summarise(lo = HDIofMCMC(value, 0.68)[1],
            md = quantile(value, 0.5),
            hi = HDIofMCMC(value, 0.68)[2])

saveRDS(volume.regional.lakes, "volume_regional_lakes.RDS")

# Regional volume, if Arctic Canada North and South is aggregated.

volume.can.agg <- all.lakes %>% 
  st_drop_geometry() %>% 
  dplyr::select("Year_new", "o1region", starts_with("Pred_")) %>%
  mutate(o1region = case_when(
    o1region == "03" ~ "03, 04",    
    o1region == "04" ~ "03, 04",  
    TRUE ~ o1region          # Keep all other values unchanged
  )) %>%
  group_by(Year_new, o1region) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = starts_with("Pred_")) %>% 
  group_by(Year_new, o1region) %>%
  summarise(lo = HDIofMCMC(value, 0.68)[1],
            md = quantile(value, 0.5),
            hi = HDIofMCMC(value, 0.68)[2],
            mdkm3 = md/1000,
            lominuskm3 = (lo-md)/1000,
            hipluskm3 = (hi-md)/1000)

# Median regional share in lake volumes.

volume.regional.lakes %>% 
  filter(Year_new == 2020) %>%
  mutate(share_global = (md / volume.all.lakes %>% 
                           filter(Year_new == 2020) %>% 
                           pull(md))*100 )

# Calculate the share of RGI regions in global glacial lake volume.

volume.global.by.pred <- all.lakes %>% 
  st_drop_geometry() %>% 
  dplyr::select("Year_new", "o1region", starts_with("Pred_")) %>%
  filter(Year_new == 2020) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE))) %>%
  as_tibble() 

volume.regional.by.pred <- all.lakes %>% 
  st_drop_geometry() %>% 
  filter(Year_new == 2020) %>%
  dplyr::select("o1region", starts_with("Pred_")) %>%
  group_by( o1region) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE))) %>%
  as_tibble()

# Identify common columns between the two tibbles (excluding o1region if 
# it's not numeric)

common_vars <- intersect(names(volume.regional.by.pred), 
                         names(volume.global.by.pred))

# Convert the one-row tibble (volume.global.by.pred) to a named list

global_values <- as.list(volume.global.by.pred)

# Divide each common numeric column in volume.regional.by.pred by its 
# corresponding global value

regional.share <- volume.regional.by.pred %>%
  mutate(across(all_of(common_vars), ~ .x / global_values[[cur_column()]])) %>%
  group_by(o1region) %>%
  summarise(md = median(c_across(starts_with("Pred_"))) * 100, 
            lo = HDIofMCMC(c_across(starts_with("Pred_")), credMass = 0.68)[1] * 100,
            hi = HDIofMCMC(c_across(starts_with("Pred_")), credMass = 0.68)[2] * 100,
            minuslo = lo-md, 
            plushi = hi-md) 

# Aggregate volume for High Mountain Asia to make it comparable to the glacier 
# volumes reported in Millan et al. (2022).

regional.share.hma.agg <- all.lakes %>% 
  st_drop_geometry() %>% 
  filter(Year_new == 2020) %>%
  dplyr::select("o1region", starts_with("Pred_")) %>%
  mutate(o1region = case_when(
    o1region == "13" ~ "13, 14, 15",    
    o1region == "14" ~ "13, 14, 15",  
    o1region == "15" ~ "13, 14, 15",  
    TRUE ~ o1region          # Keep all other values unchanged
  )) %>%
  group_by(o1region) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE))) %>%
  as_tibble() %>%
  mutate(across(all_of(common_vars), ~ .x / global_values[[cur_column()]])) %>%
  group_by(o1region) %>%
  summarise(md = median(c_across(starts_with("Pred_"))) * 100, 
            lo = HDIofMCMC(c_across(starts_with("Pred_")), credMass = 0.68)[1] * 100,
            hi = HDIofMCMC(c_across(starts_with("Pred_")), credMass = 0.68)[2] * 100,
            minuslo = lo-md, 
            plushi = hi-md) 

# Calculate global change in glacial lake volume.

change.all.lakes <- all.lakes %>% 
  st_drop_geometry() %>% 
  dplyr::select("Year_new", "o1region", starts_with("Pred_")) %>%
  group_by(Year_new) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("Pred_")) %>% 
  group_by(name) %>%
  summarise(change = (value[Year_new == 2020]-value[Year_new == 1990])/ value[Year_new == 1990] ) %>%
  summarise(lo = HDIofMCMC(change, 0.68)[1]*100,
            md = quantile(change, 0.5)*100,
            hi = HDIofMCMC(change, 0.68)[2]*100)

# Regional absolute change in volume

change.regional.absolute <- all.lakes %>% 
  st_drop_geometry() %>% 
  dplyr::select("Year_new", "o1region", starts_with("Pred_")) %>%
  group_by(Year_new, o1region) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE)))  %>%
  pivot_longer(cols = starts_with("Pred_")) %>%
  pivot_wider(names_from = Year_new, values_from = value) %>%  # Reshape to wide format
  mutate(change = `2020` - `1990`) %>%  # Calculate absolute change
  group_by(o1region) %>%  # Group by region
  summarise(
    lo = HDIofMCMC(change, 0.68)[1],
    md = quantile(change, 0.5),
    hi = HDIofMCMC(change, 0.68)[2]
  )


# Regional relative change in volume

gc()

change.regional.relative <- all.lakes %>% 
  st_drop_geometry() %>% 
  dplyr::select("Year_new", "o1region", starts_with("Pred_")) %>%
  group_by(Year_new, o1region) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE)))  %>%
  pivot_longer(cols = starts_with("Pred_")) %>%
  pivot_wider(names_from = Year_new, values_from = value) %>%  # Reshape to wide format
  mutate(change_relative = ((`2020` - `1990`) / (`1990`)) *100) %>%  # Calculate absolute change
  group_by(o1region) %>%  # Group by region
  summarise(
    lo = HDIofMCMC(change_relative, 0.68)[1],
    md = quantile(change_relative, 0.5),
    hi = HDIofMCMC(change_relative, 0.68)[2],
    gt0 = sum(change_relative > 0) / n()
  )

# Cumulative volume and sorted

vol.cumsum.2020 <- apply(pred.vol.all[all.lakes$Year_new == 2020, ], 2, function (x) {
  
  vol.sort <- sort(x)
  
  o <-  cumsum(vol.sort)
  
  df <- tibble(volsort = vol.sort,
               frac_cumsum = o / max(o))
  
})

for(i in 1:length(vol.cumsum.2020)) {
  
  vol.cumsum.2020[[i]]$pred <- names(vol.cumsum.2020)[i]
  
}

vol.cumsum.2020 <- do.call(rbind, vol.cumsum.2020)

gc()

# Define fixed percentages of cumulative volume, from 10^-10 (almost zero) 
# to 0. Log10 at 0 is not defined, so we do not go simply from zero to 1 on the
# y-axis.

fixed.perc <- 10^seq(from = -10, to = 0, length.out = 150)

# For every prediction: Linearly interpolate between the cumulative

interp.perc.2020 <- vol.cumsum.2020 %>%
  group_by(pred) %>%
  do({
    data.frame(frac_cumsum = fixed.perc,
               volsort = approx(.$frac_cumsum, 
                                .$volsort, 
                                xout = fixed.perc, rule = 2:1)$y)
  }) %>%
  ungroup()

# Aggregate over all lakes: 

summary.cumsum.2020 <- interp.perc.2020 %>%
  ungroup() %>%
  group_by(frac_cumsum) %>%
  summarize(md = median(volsort),
            lo = HDIofMCMC(volsort, 0.68)[1],
            hi = HDIofMCMC(volsort, 0.68)[2])

vol.per.lake.stats.2020 <- vol.per.lake[all.lakes$Year_new == 2020, ] %>%
  arrange(md) %>%
  mutate(csum50 = cumsum(md),
         frac50 = csum50/csum50[n()],
         frac_of_tot50 = md/csum50[n()])

min.pr <- min(log10(pred.vol.all))
max.pr <- max(log10(pred.vol.all))

dens.2020 <- apply(pred.vol.all[all.lakes$Year_new == 2020, ], 
                   MARGIN = 2, function (x) {
                     
                     h2 <- density(log10(x), bw = 0.125, from = min.pr, to = max.pr)
                     h2$y
                     
                   })

brea <- density(log10(pred.vol.all$Pred_1), bw = 0.125,
                from = min.pr, 
                to = max.pr)$x

dens.summary.2020 <- apply(dens.2020, MARGIN = 1, 
                           function (x) {
                             
                             tibble(lo = HDIofMCMC(x, 0.68)[1],
                                    md = quantile(x, 0.5),
                                    hi = HDIofMCMC(x, 0.68)[2])
                           }) %>% 
  bind_rows() %>%
  mutate(brea = 10^brea)

## Same for 1990

vol.cumsum.1990 <- apply(pred.vol.all[all.lakes$Year_new == 1990, ], 2, function (x) {
  
  vol.sort <- sort(x)
  
  o <-  cumsum(vol.sort)
  
  df <- tibble(volsort = vol.sort,
               frac_cumsum = o / max(o))
  
})

for(i in 1:length(vol.cumsum.1990)) {
  
  vol.cumsum.1990[[i]]$pred <- names(vol.cumsum.1990)[i]
  
}

vol.cumsum.1990 <- do.call(rbind, vol.cumsum.1990)

gc()

# For every prediction: Linearly interpolate between the cumulative

interp.perc.1990 <- vol.cumsum.1990 %>%
  group_by(pred) %>%
  do({
    data.frame(frac_cumsum = fixed.perc,
               volsort = approx(.$frac_cumsum, 
                                .$volsort, 
                                xout = fixed.perc, rule = 2:1)$y)
  }) %>%
  ungroup()

# 

summary.cumsum.1990 <- interp.perc.1990 %>%
  ungroup() %>%
  group_by(frac_cumsum) %>%
  summarize(md = median(volsort),
            lo = HDIofMCMC(volsort, 0.68)[1],
            hi = HDIofMCMC(volsort, 0.68)[2])

vol.per.lake.stats.1990 <- vol.per.lake[all.lakes$Year_new == 1990, ] %>%
  arrange(md) %>%
  mutate(csum50 = cumsum(md),
         frac50 = csum50/csum50[n()],
         frac_of_tot50 = md/csum50[n()])

dens.1990 <- apply(pred.vol.all[all.lakes$Year_new == 1990, ], 
                   MARGIN = 2, function (x) {
                     
                     h2 <- density(log10(x), bw = 0.125, from = min.pr, to = max.pr)
                     
                     h2$y
                     
                   })

dens.summary.1990 <- apply(dens.1990, MARGIN = 1, 
                           function (x) {
                             
                             tibble(lo = HDIofMCMC(x, 0.68)[1],
                                    md = quantile(x, 0.5),
                                    hi = HDIofMCMC(x, 0.68)[2])
                           }) %>% 
  bind_rows() %>%
  mutate(brea = 10^brea)

summary_cumsum <- bind_rows(summary.cumsum.1990 %>% mutate(Year_new = 1990),
                            summary.cumsum.2020 %>% mutate(Year_new = 2020)) %>%
  mutate(Year_new = as_factor(Year_new))

dens_summary <- bind_rows(dens.summary.1990 %>% mutate(Year_new = 1990),
                          dens.summary.2020 %>% mutate(Year_new = 2020)) %>%
  mutate(Year_new = as_factor(Year_new))


# Plot with densities of lake volumes in 1990 and 2020, and their cumulated
# volume

density.cumul.plot <- ggplot(summary_cumsum,
                             aes(x = md,
                                 y = frac_cumsum,
                                 group = Year_new, 
                                 fill = Year_new, 
                                 color = Year_new)) +
  geom_lineribbon(aes(xmin = lo,
                      xmax = hi),
                  alpha = 0.5,
                  linewidth = 0) +
  geom_line(linewidth = 0.8) +  # Line for each Year_new
  geom_lineribbon(data = dens_summary,
                  mapping = aes(x = brea,
                                y = md*2,
                                ymin = lo*2,
                                ymax = hi*2,
                                group = Year_new, 
                                fill = Year_new, 
                                color = Year_new),
                  inherit.aes = F,
                  alpha = 0.5,
                  linewidth = 0) +
  geom_line(data = dens_summary,
            mapping = aes(x = brea,
                          y = md*2,
                          group = Year_new, 
                          color = Year_new),
            inherit.aes = F,
            linewidth = 0.8) +  # Line for each Year_new
  scale_fill_manual(name = "Year", values = c("lightblue", "darkorange")) +  # Both ribbons grey
  scale_color_manual(name = "Year", values = c("navy", "darkorange4")) +
  geom_rug(vol.per.lake.stats.2020 %>% filter(frac50 > 0.5),
           sides = "r",
           linewidth = 0.25,
           mapping = aes(y = frac50),
           color = "darkorange4", inherit.aes = F) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 3),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_continuous(sec.axis = sec_axis(transform=~., 
                                         name = "Cumulated volume")) +
  labs(y = "Density",
       x = expression("Lake volume [" * 10^6 ~ "m"^3 * "]")) +
  theme_bw() +
  theme(text = element_text(size = 7),
        legend.position.inside = c(0.15, 0.85),
        aspect.ratio = 1,
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank())  + 
  guides(fill   = guide_legend(position = "inside"),
         colour = guide_legend(position = "inside")) 

ggsave(filename = "density_and_cumulative_dist.pdf",
       plot = density.cumul.plot, 
       width = 90,
       height = 90,
       units = "mm")

### Same plot without cumulative distribution

density.plot <- ggplot(data = dens_summary,
                  mapping = aes(x = brea/10^3,
                                y = md,
                                ymin = lo,
                                ymax = hi,
                                group = Year_new, 
                                fill = Year_new, 
                                color = Year_new)) +
  geom_lineribbon(alpha = 0.5,
                  linewidth = 0) +
  geom_line(mapping = aes(x = brea/10^3,
                          y = md,
                          group = Year_new, 
                          color = Year_new),
            linewidth = 0.8) +  # Line for each Year_new
  scale_fill_manual(name = "Year", values = c("lightblue", "darkorange")) +  # Both ribbons grey
  scale_color_manual(name = "Year", values = c("navy", "darkorange4")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 3),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  labs(y = "Density",
       x = "Lake volume [km³]") +
  theme_bw() +
  theme(text = element_text(size = 7),
        legend.position.inside = c(0.15, 0.85),
        aspect.ratio = 1,
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank())  + 
  guides(fill   = guide_legend(position = "inside"),
         colour = guide_legend(position = "inside")) 

ggsave(filename = "density_volume.pdf",
       plot = density.plot, 
       width = 70,
       height = 70,
       units = "mm",
       device = cairo_pdf, 
       family = "Arial")

################################################################################

# Regional relative change in volume

gc()

test.mat <- all.lakes %>% 
  st_drop_geometry() %>% 
  select(starts_with("Pred_")) %>%
  as.matrix() %>%
  apply(., 1, function (x) {
    
    HDI <- HDIofMCMC(x, 0.68)
    
    out <- x[(x >= HDI[1]) & (x <= HDI[2])]
    
    if (length(out) == 682) {out <- out[-682]}
    
    out
    
  })
  
lakes.trim <- as_tibble(t(test.mat)) 
colnames(lakes.trim) <- paste0("Pred_", 1:ncol(lakes.trim))

change.size.classes <- lakes.trim %>%
  mutate(Area = all.lakes$Area,
         Year_new = all.lakes$Year_new) %>%
  mutate(interval = cut(Area, breaks = c(0, 0.1, 1, 10, 1000))) %>% 
  group_by(Year_new, interval) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE)))  %>%
  summarise( lo = HDIofMCMC(value, 0.68)[1],
             md = quantile(value, 0.5),
             hi = HDIofMCMC(value, 0.68)[2] )

change.size.classes <- all.lakes %>%
  st_drop_geometry() %>%
  mutate(interval = cut(Area, breaks = c(0, 0.1, 1, 10, 1000))) %>% 
  group_by(Year_new, interval) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE))) %>%
  rowwise() %>%
  mutate(
    lo = HDIofMCMC(c_across(starts_with("Pred")), 0.68)[1],
    md = quantile(c_across(starts_with("Pred")), 0.5),
    hi = HDIofMCMC(c_across(starts_with("Pred")), 0.68)[2]
    
  ) 

labs <- all.lakes %>%
  st_drop_geometry() %>%
  mutate(interval = cut(Area, breaks = c(0, 0.1, 1, 10, 1000))) %>% 
  group_by(Year_new, interval) %>% 
  count(interval)

plot_data <- change.size.classes %>%
  left_join(labs, by = c("Year_new", "interval"))

lakes.in.size.classes <- ggplot(plot_data, aes(x = interval, 
                      y = md/10^3, 
                      fill = factor(Year_new))) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_errorbar(aes(ymin = lo/10^3, 
                    ymax = hi/10^3),
                position = position_dodge(width = 0.9), 
                width = 0.2) +
  geom_text(aes(label = scales::comma(n), 
                y = hi/10^3 + 50,
                color = factor(Year_new)), 
            position = position_dodge(width = 0.9),
            vjust = 0, size = 7/.pt) +
  labs(x = "Lake area (km²)", y = "Volume (km³)", fill = "Year") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = c("navy", "darkorange")) +
  scale_color_manual(values = c("navy", "darkorange")) +
  theme_bw() +
  theme(text = element_text(size = 7),
        legend.position = "none",
        aspect.ratio = 1,
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank())

ggsave(filename = "lakes_in_size_classes.pdf",
       plot = lakes.in.size.classes, 
       width = 70,
       height = 70,
       units = "mm",
       device = cairo_pdf, 
       family = "Arial")

################################################################################

change.regional.relative <- all.lakes %>% 
  st_drop_geometry() %>% 
  dplyr::select("Year_new", "o1region", starts_with("Pred_")) %>%
  group_by(Year_new, o1region) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE)))  %>%
  pivot_longer(cols = starts_with("Pred_")) %>%
  pivot_wider(names_from = Year_new, values_from = value) %>%  # Reshape to wide format
  mutate(change_relative = ((`2020` - `1990`) / (`1990`)) *100) %>%  # Calculate absolute change
  group_by(o1region) %>%  # Group by region
  summarise(
    lo = HDIofMCMC(change_relative, 0.68)[1],
    md = quantile(change_relative, 0.5),
    hi = HDIofMCMC(change_relative, 0.68)[2],
    gt0 = sum(change_relative > 0) / n()
  )

# Cumulative volume and sorted

vol.cumsum.2020 <- apply(pred.vol.all[all.lakes$Year_new == 2020, ], 2, function (x) {
  
  vol.sort <- sort(x)
  
  o <-  cumsum(vol.sort)
  
  df <- tibble(volsort = vol.sort,
               frac_cumsum = o / max(o))
  
})


 
################################################################################

# Regional relative change in volume

# Cumulative volume and sorted

vol.cumsum.2020 <- apply(pred.vol.all[all.lakes$Year_new == 2020, ], 2, function (x) {
  
  vol.sort <- sort(x)
  
  o <-  cumsum(vol.sort)
  
  df <- tibble(volsort = vol.sort,
               frac_cumsum = o / max(o))
  
})


################################################################################

# Where is the peak of the distribution in lake volumes? 

dens_summary %>% 
  filter(Year_new == 2020) %>% 
  slice_max(md) %>% pull(brea)

# How many lakes impound 50% of the total volume?

vol.per.lake.stats.2020 %>%
  filter(frac50 > 0.50)

tail(vol.per.lake.stats.2020)

# Filter lakes only for 2020

l.2020 <- all.lakes.summary %>% filter(Year_new == 2020)

top.largest.2020 <- l.2020 %>%
  slice_max(order_by = md, n = nrow(vol.per.lake.stats.2020 %>%
                                      filter(frac50 > 0.50)))

all.lakes.1990 <- all.lakes %>%
  filter(Year_new == 1990)

int.largest <- st_intersects(top.largest.2020, all.lakes.1990, sparse = F)

change.of.the.largest <- apply(int.largest, 1, function(x)  {
  
  tf <- which(x)
  
  if (length(tf) >= 1) {
  
  all.lakes.1990[tf, ] %>%
    summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE))) %>%
    as_tibble() %>%
    summarise(md1990 = median(c_across(starts_with("Pred_"))) , 
              lo1990 = HDIofMCMC(c_across(starts_with("Pred_")), credMass = 0.68)[1] ,
              hi1990 = HDIofMCMC(c_across(starts_with("Pred_")), credMass = 0.68)[2] ) 
  
  } else { tibble(md1990 = 0,
                  lo1990 = 0,
                  hi1990 = 0)}
  
})

bind_cols(top.largest.2020, bind_rows(change.of.the.largest)) %>%
  filter(hi1990 < lo)

# Then, we split the data into 2020 and 1990.

l.2020 <- all.lakes.summary.no.sup %>% filter(Year_new == 2020)

# We find lakes that intersect in both years.

int.2020.1990 <- st_intersects(l.2020, all.lakes.1990, sparse = F)

# We deem a change credible, if the upper bound of the 68% HDI in 1990 is smaller
# than the lower bound of the 68% HDI in 2020. If more than one lake intersects, 
# we stack the volumes of the individual lakes and calculate the 68% HDI.

individual.changes <- apply(int.2020.1990, 1, function (x) {
  
  tf <- which(x)
  
  if (length(tf) >= 1) {
    
    all.lakes.1990[tf, ] %>%
      summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE))) %>%
      as_tibble() %>%
      summarise(md1990 = median(c_across(starts_with("Pred_"))) , 
                lo1990 = HDIofMCMC(c_across(starts_with("Pred_")), credMass = 0.68)[1] ,
                hi1990 = HDIofMCMC(c_across(starts_with("Pred_")), credMass = 0.68)[2] ) 
    
  } else { tibble(md1990 = 0,
                  lo1990 = 0,
                  hi1990 = 0)}
  }
)

individual.changes <- bind_rows(individual.changes)

# Number of lakes that existed in 1990 and 2020.

indiv.volumes.1990.2020 <- bind_cols(l.2020, individual.changes)

saveRDS(indiv.volumes.1990.2020, "individual_volumes_in_1990_2020_no_sup.RDS")

# Lakes that credibly gained volume. Lake must be present in 1990

cred.growth <- bind_cols(l.2020, individual.changes) %>%
  filter(md1990 != 0,
         lo > hi1990) 

o <- st_intersects(cred.growth, all.lakes.1990, sparse = F)


area.growth.1990.2020 <- sapply(1: nrow(o), function(i) {

  c(area1990 = sum(all.lakes.1990[o[i, ], ]$Area),
    area2020 = cred.growth$Area[i])
  
})

ag9020 <- t(area.growth.1990.2020) %>% 
  as_tibble() %>% 
  mutate(perc_ch = ((area2020-area1990)/ area1990) * 100) 

quantile(ag9020$perc_ch, c(0.16, 0.5, 0.84))


# Lakes that credibly gained volume until 2020. Can also be new lakes (NA in 1990).

l.2020.w.change <- bind_cols(l.2020, individual.changes) %>%
  filter(!is.na(md1990),
         hi1990 < lo)

# Only change in the median, without credibility interval

bind_cols(l.2020, individual.changes) %>%
  filter(!is.na(md1990),
         md1990 < md)

# Determine the volume of lakes that are new in 2020 (and not supraglacial)

l.2020.without.change <- bind_cols(l.2020, individual.changes) %>%
  filter(md1990 == 0) %>%
  st_drop_geometry() %>%
  summarise(medsum = sum(md))
  

################################################################################

# Obtain the regional share of lake volume from the global lake volume 

# First, calculate the global volume for each year (1990 and 2020) for each draw
# from the posterior predictive distribution.

df_global <- all.lakes %>% 
  st_drop_geometry() %>% 
  select("Year_new", "o1region", starts_with("Pred_")) %>%
  group_by(Year_new) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE))) %>%
  mutate(o1region = "Global")

# Do the same, calculate the regional glacial lake volume for both 1990 and 2020.

df_regional <- all.lakes %>% 
  st_drop_geometry() %>% 
  select("Year_new", "o1region", starts_with("Pred_")) %>%
  group_by(Year_new, o1region) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE)))

# Divide the regional volume by the global volume to obtain the
# regional share of global lake volume.

vol.stats <- lapply(c(1990, 2020), function(i) {
  
  global.year <- df_global %>% filter(Year_new == i)
  regional.year <- df_regional %>% filter(Year_new == i)
  
  # Iterate over all regions
  
  for (j in unique(df_regional$o1region)) {
    
    # Filter the predictions of the regional glacial lake volume
    
    reg <- regional.year %>% 
      ungroup() %>%
      select(starts_with("Pred_"))
    
    # Total regional volume: median and 68% HDI.
    
    vol <- reg %>% apply(., 1, function(x) {
      hdi <- HDIofMCMC(x, 0.68)
      md  <- quantile(x, 0.5)
      q <- c(hdi[1], md, hdi[2])
      
      names(q) <- c("lo", "md", "hi")
      q}) %>% t() %>% as_tibble()
    
    # Calculate the regional fraction compared to the total volume
    
    frac <- apply(reg, 1, function(x) {as.numeric(x) / global.year %>% 
        select(starts_with("Pred_")) %>% 
        as.numeric()}) %>% 
      apply(., 2, function(x) {
        hdi <- HDIofMCMC(x, 0.68)
        md  <- quantile(x, 0.5)
        q <- c(hdi[1], md, hdi[2])
        names(q) <- c("qfraclo", "qfracmd", "qfrachi")
        q}) %>%
      t() %>%
      as_tibble() %>% 
      mutate(Year_new = i,
             o1region = regional.year$o1region)
    
    vol.frac <- bind_cols(frac, vol)
    
    return(vol.frac)
  }
})

vol.stats <- bind_rows(vol.stats)


################################################################################

# Maps showing

sf_use_s2(FALSE)

date_line <- st_linestring(matrix(c(180, -90, 180, 90), ncol = 2, byrow = TRUE))
date_line_sf <- st_sfc(date_line, crs = 4326)

world_sf <- st_read("D:/nrc_user/veh/continent shapefile/continent_dissolve.shp") %>%
  st_simplify(dTolerance = 0.01) %>%
  st_split(., date_line_sf)

world_robin <- st_transform(world_sf, crs = "+proj=robin") %>% 
  st_make_valid()

grat_robin <- 
  st_graticule(lon = c(-179.9, seq(-150, 150, 25), 179.9),
               lat = seq(-75, 75, 25))  %>%
  st_transform_proj(crs =  "+proj=robin")

# vectors of latitudes and longitudes that go once around the 
# globe in 1-degree steps

lats <- c(90:-90, -90:90, 90)
longs <- c(rep(c(180, -180), each = 181), 180)

# turn into correctly projected sf collection

robin_outline <- 
  list(cbind(longs, lats)) %>%
  st_polygon() %>%
  st_sfc(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% 
  st_sf() %>%
  st_transform_proj(crs =  "+proj=robin") # transform to Robinson

# Example colors

land_color  <- "#2c5559" # "#779396"  # Dark green-gray
ocean_color <- "#EDF6F9" # Very light cyan

o1 <- st_read("D:/nrc_user/veh/RGI2000-v7.0-regions/RGI2000-v7.0-o1regions.shp") %>%
  st_make_valid()

o1_robin <-  st_transform(o1, crs = "+proj=robin")

region_pts_robin <- st_point_on_surface(o1_robin) %>%
  left_join(., vol.stats %>% filter(Year_new == 2020),
            by = 'o1region') 


################################################################################

df_points <- region_pts_robin %>% 
  cbind(st_coordinates(region_pts_robin)) %>% 
  filter(o1region != "19",
         o1region != "20") %>%
  mutate(hi = hi/1000,
         md = md/1000,
         lo = lo/1000,
         qfracmd = qfracmd*100,
         qfraclo = qfraclo*100,
         qfrachi = qfrachi*100,
         lab = paste0(full_name, ":\n", 
                      round(qfracmd, 1), "(-",
                      round((qfracmd-qfraclo), 1), "/+",
                      round((qfrachi-qfracmd), 1), ")")) 


df_points <- region_pts_robin %>% 
  cbind(st_coordinates(region_pts_robin)) %>% 
  filter(o1region != "19",
         o1region != "20") %>%
  mutate(hi = hi/1000,
         md = md/1000,
         lo = lo/1000,
         qfracmd = qfracmd*100,
         qfraclo = qfraclo*100,
         qfrachi = qfrachi*100,
         lab = paste0(full_name, ":\n", 
                      round(qfracmd, 1), "(-",
                      round((qfracmd-qfraclo), 1), "/+",
                      round((qfrachi-qfracmd), 1), ")")) 


map.vol <- df_points %>%
  ggplot(aes(x = X, y = Y, size = hi)) + 
  geom_point(stroke = 0.1) +
  geom_sf(data = robin_outline, fill = ocean_color, color = NA, inherit.aes = F) +
  geom_sf(data = world_robin, color = NA, fill = land_color, inherit.aes = F) +
  geom_sf(data = o1_robin  %>% filter(o1region != "19", 
                                      o1region != "20"), 
          color = "white",
          fill = fill_alpha("lightblue", 0.5), inherit.aes = F) +
  scale_size_continuous(name = "Regional glacier\nlake volume [km³]",
                        range = c(3, 12)) +
  geom_sf(shape = 21, fill = "orangered2", color = "orangered3", stroke = 0.2) +
  geom_point(aes(size = lo, fill = qfracmd), shape = 21, color = "orangered3", stroke = 0.2) +
  geom_point(aes(size = md), fill = NA, color = "navy", shape = 21) +
  scale_fill_gradient(name = "Regional share [%] of\nglobal glacier lake volume",
                      low = "grey95", high = "steelblue4") +
  geom_text_repel(
    aes(label = lab, 
        geometry = geometry),
    stat = "sf_coordinates",
    force = 10,
    color = "black",     # text color
    bg.r = 0.2,          # shadow radius
    max.overlaps = Inf,
    segment.color = "grey45",
    min.segment.length = 3,
    size = 20/.pt) +
  theme_map() +
  theme( legend.position = "bottom",  # Center at the bottom
         legend.justification = c(0.5, 0),  # Align the legend's center
         legend.title = element_text(size = 20), 
         legend.text = element_text(size = 20),
         legend.box = "vertical") +
  guides(color = guide_colorbar(order = 1),  # Colorbar appears in the first row
         size  = guide_legend(order = 2))      # Size legend appears in the second row

################################################################################

# Change in regional glacial lake volume between 1990 and 2020.

df_change <- st_point_on_surface(o1_robin) %>%
  left_join(., change.regional.relative,
            by = 'o1region') %>% 
  cbind(st_coordinates(.)) %>% 
  filter(o1region != "19",
         o1region != "20") %>%
  mutate(gt0 = gt0 * 100, 
         q50_simple = if_else(md > 100, 100, md),
         lab = paste0(full_name, ":\n", 
                      round(md, 1), "(-",
                      round((md-lo), 1), "/+",
                      round((hi-md), 1), ")")) 

map.change <- df_change %>%
  ggplot(aes(x = X, y = Y, size = q50_simple, fill = gt0)) + 
  geom_sf(data = robin_outline, fill = ocean_color, color = NA, inherit.aes = F) +
  geom_sf(data = world_robin, color = NA, fill = land_color, inherit.aes = F) +
  geom_sf(data = o1_robin  %>% filter(o1region != "19", 
                                      o1region != "20"), 
          color = "white",
          fill = fill_alpha("lightblue", 0.5), inherit.aes = F) +
  geom_point(shape = 21, stroke = 0.1) +
  scale_size_continuous(name = "Median change in\nglacier lake volume [%]\n1990-2020",
                        range = c(1, 12)) +
  geom_sf(aes(size = q50_simple, fill = gt0), shape = 21) +
  scale_fill_gradientn(name = "Proportion [%] of\nposterior mass >0",
    colors = c("lightblue", "white", "darkred"),
    values = scales::rescale(c(0.25, 0.5, 1), from = c(0.25, 1)),
    limits = c(25, 100),
    breaks = c(25, 50, 100)) +
  geom_text_repel(
    aes(label = lab, 
        geometry = geometry),
    stat = "sf_coordinates",
    force = 10,
    color = "black",     # text color
    bg.r = 0.2,          # shadow radius
    max.overlaps = Inf,
    segment.color = "grey45",
    min.segment.length = 3,
    size = 7/.pt) +
  theme_map() +
  theme( legend.position = "bottom",  # Center at the bottom
         legend.justification = c(0.5, 0),  # Align the legend's center
         legend.title=element_text(size = 7), 
         legend.text=element_text(size = 7),
         legend.box = "vertical") +
  guides(color = guide_colorbar(order = 1),  # Colorbar appears in the first row
         size  = guide_legend(order = 2))      # Size legend appears in the second row

ggsave(plot     = map.change, 
       filename = "Global_lake_change.pdf",
       width    = 180,
       height   = 120,
       units    = "mm")

volume.maps <- plot_grid(map.vol,
                         map.change, 
                         nrow = 2,
                         labels = c('a', 'b'),
                         label_size = 8)

ggsave(filename = "volume_maps.pdf",
       plot = volume.maps, 
       height = 200,
       width = 180,
       units = "mm")
