### Compare the regional glacial lake volume with the regional glacier volume;
### Compare the regional change in glacial lake volume with rates of glacial
### mass loss.
### Assess the flowpath distance and elevation distribution of glacial lakes.

## Load packages that we'll need in the course of this script

# If packages haven't been installed beforehand, 
# first use install.packages("packagename")

require(tidyverse)
require(scales)
require(brms)
require(tidybayes)
require(ggrepel)
require(cowplot)
require(sf)

# Set working directory

setwd("D:/nrc_user/veh/Zhang_glacial_lakes_global/")

# O1 regions from the RGI.

o1regions <- st_read("D:/nrc_user/veh/RGI2000-v7.0-regions/RGI2000-v7.0-o1regions.shp")

# Load lake data

all.lakes <- readRDS("all_lakes_with_volumes.RDS")

# Catchment stats for lakes

lake.stats.sf <- readRDS("lakes_with_catchment_and_landover.rds")

# Define the join key to join predicted lake volumes with catchments stats 
# via the common key UniqueID.

join_key <- "UniqueID"

# Determine the columns in df_right that are not in df_left (excluding the join key)
new_vars <- setdiff(names(lake.stats.sf), c(names(all.lakes), join_key))

# Perform the left join by including only the join key and the new variables from df_right
lakes.w.stats <- all.lakes %>%
  left_join(., dplyr::select(lake.stats.sf, all_of(c(join_key, new_vars))), by = join_key)

# Add the function HDIofMCMC from John Kruschke to this directory. We will use
# this function to derive the highest density intervals of posterior 
# (predictive) distributions of lake volumes.

source("HDIofMCMC.R")

# Regional comparison of lake volumes and lake volume change with 
# regional ice volume and ice volume change

# Comparison with the ice volumes estimated by Millan et al. (2022),
# https://www.nature.com/articles/s41561-021-00885-z
# Millan et al. aggregated Alaska and BC (RGI regions 01 and 02), and the three
# regions in High Mountain Asia (13, 14, 15)

lakes.millan <- lakes.w.stats %>% filter(Year_new == 2020)
lakes.millan$o1region[(lakes.millan$o1region == "01") | 
                        (lakes.millan$o1region == "02")] <- "01, 02"
lakes.millan$o1region[(lakes.millan$o1region == "13") | 
                        (lakes.millan$o1region == "14") |
                        (lakes.millan$o1region == "15")] <- "13, 14, 15"

# Calculate regional lake volumes based on these regions.

vol.stats.millan <- lakes.millan %>% 
  st_drop_geometry() %>% 
  dplyr::select( "o1region", starts_with("Pred_")) %>%
  group_by(o1region) %>%
  summarise(across(starts_with("Pred_"), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = starts_with("Pred_")) %>% 
  group_by(o1region) %>%
  summarise(qvollo = HDIofMCMC(value, 0.68)[1],
            qvolmd = quantile(value, 0.5),
            qvolhi = HDIofMCMC(value, 0.68)[2])

# Load Millan's data, stored in an Excel sheet.

vol.gl <- readxl::read_xlsx("Millan_glacier_volumes.xlsx", sheet = 1)

# Join glacier and ice volume data. Convert both volumes to cubic kilometers
# and log10-transform the data.

lake.gl.vol <- left_join(vol.stats.millan, vol.gl, by = "o1region") %>%
  mutate(Ice_volume_km3 = Ice_volume_1000km3 * 10^3,
         Ice_volume_error_km3 = Ice_volume_error *10^3,
         qvolmd_km3 = qvolmd/10^3,
         qvollo_km3 = qvollo/10^3,
         qvolhi_km3 = qvolhi/10^3,
         uncertainty = qvolhi_km3 - qvolmd_km3,
         ice_log = log10(Ice_volume_km3),
         vol_log = log10(qvolmd_km3))

# Calculate a trend between regional ice volume and lake volume.
# Regional ice volumes are assumed to be two orders of magnitude larger 
# than lake volumes (see prior for the intercept).

prior1 <- prior(normal(2, 2.5), class = "Intercept") +
  prior(normal(0, 2.5), class = "b") +
  prior(normal(2, 5), class = "nu") +
  prior(normal(0, 2.5), class = "sigma") 

fit1 <- brm(bf(ice_log ~ vol_log),
            data = lake.gl.vol, 
            family = "student",
            prior = prior1,
            iter = 4000,
            warmup  = 1000,
            chains  = 4,
            cores   = 4, 
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15))

### Plot model performance statistics.

fit1_plot <- plot(fit1, ask = F, nvariables = 4, theme = theme_bw(base_size =  7)) 

fit1.slope.post <- posterior_samples(fit1, pars = "b_vol_log")

median(fit1.slope.post$b_vol_log)
HDIofMCMC(fit1.slope.post$b_vol_log, 0.68)

### Predictions for new data.

# Generate a sequence of new log10-transformed lake volumes

vol_log <- seq(min(lake.gl.vol$vol_log, na.rm = T), 
               max(lake.gl.vol$vol_log, na.rm = T), length.out = 100)

# Draw from the posterior predictive distribution

pred.grid.ice_log <- add_predicted_draws(
  object = fit1, 
  newdata = data.frame(vol_log = vol_log),
  value = "ice_log", 
  ndraws = 4000,
  re_formula = NA) 

# Remove outliers from the prediction.

for (i in unique(pred.grid.ice_log$.row)) {
  
  row <- pred.grid.ice_log$ice_log[pred.grid.ice_log$.row == i]
  
  HDIrange <- HDIofMCMC(row, credMass = 0.99)
  
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
  
  pred.grid.ice_log$ice_log[pred.grid.ice_log$.row == i] <- row
  
}

# Plot regional lake volume versus regional glacier volume.

lake.vol.plot <- ggplot(data = lake.gl.vol,
                        aes(qvolmd_km3,
                            Ice_volume_km3,
                            size = uncertainty)) +
  stat_lineribbon(data = pred.grid.ice_log,
                  mapping = aes(x = 10^vol_log,
                                y = 10^ice_log), 
                  .width = 0.68,
                  point_interval = median_hdi,
                  fill = "grey66", 
                  inherit.aes = F) +
  geom_linerange(aes(ymin = Ice_volume_km3-Ice_volume_error_km3, 
                     ymax = Ice_volume_km3+Ice_volume_error_km3),
                 linewidth = 1,
                 color = "chocolate4") +
  geom_point(fill = "chocolate3", color = "white", shape = 21, stroke = 0.3) +
  geom_text_repel(aes(label = Region),
                  size = 7/.pt,
                  max.overlaps = 25) +
  scale_size_continuous(
    name = "Width of 68% HDI of\nlake volume [km設",
    range = c(1,10)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 3),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 3),
                labels = trans_format("log10", math_format(10 ^ .x)))  +
  labs(x = "Regional lake volume [km設 in 2020",
       y = "Regional glacier volume [km設\nin ~2017/2018") +
  theme_bw() +
  theme(text = element_text(size = 7),
        aspect.ratio = 0.8,
        legend.position = "bottom")

################################################################################

# Same analysis now for glacier volume change and lake volume change.

# Aggregate lake volume change on regional level between 1990 and 2020.

change.regional.absolute <- lakes.w.stats %>% 
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
    hi = HDIofMCMC(change, 0.68)[2])

# Load glacier mass loss data from Hugonnet et al. (2021),
# available at https://doi.org/10.6096/13

gl.loss <- read.csv("Hugonnet_glacier_loss.txt", 
                    sep = ",", header = T,
                    colClasses = c("o1region" = "character")) %>%
  as_tibble() %>%
  mutate(absolute_vol_km3 = abs(dvol/10^9),
         error_vol_km3 =  err_dvol/10^9)

# Join lake and glacier volume change.

lakes.gl.loss <- left_join(change.regional.absolute, gl.loss, by = "o1region") %>%
  mutate(qvolmd_km3 = md/10^3,
         qvollo_km3 = lo/10^3,
         qvolhi_km3 = hi/10^3, 
         uncertainty = qvolhi_km3-qvolmd_km3,
         vol_change_log = log10(qvolmd_km3),
         ice_loss_log = log10(absolute_vol_km3)) %>%
  filter(!is.na(vol_change_log))

# Add RGI O1 code

o1regions <- st_read("D:/nrc_user/veh/RGI2000-v7.0-regions/RGI2000-v7.0-o1regions.shp") %>%
  st_drop_geometry()

lakes.gl.loss <- left_join(lakes.gl.loss, o1regions, by = "o1region")

# Assess trend of glacier mass loss with lake volume change

prior2 <- prior(normal(2, 2.5), class = "Intercept") +
  prior(normal(0, 2.5), class = "b") +
  prior(normal(2, 5), class = "nu") +
  prior(normal(0, 2.5), class = "sigma") 

fit2 <- brm(bf(ice_loss_log ~ vol_change_log),
            data = lakes.gl.loss, 
            family = "student",
            prior = prior2,
            iter = 4000,
            warmup  = 1000,
            chains  = 4,
            cores   = 4, 
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15))

fit2.slope.post <- posterior_samples(fit2, pars = "b_vol_change_log")

median(fit2.slope.post$b_vol_change_log)
HDIofMCMC(fit2.slope.post$b_vol_change_log, 0.68)

### Plot model performance statistics

fit2_plot <- plot(fit2, ask = F, nvariables = 4, theme = theme_bw(base_size =  7)) 

# Combine both plots showing the performance of models fit1 and fit2.

plot.fit1.fit2 <- plot_grid(fit1_plot[[1]], 
                            fit2_plot[[1]], 
                            nrow = 2,
                            axis = "l",
                            labels = c('a', 'b'),
                            align = "v",
                            label_size = 8)

# Save the plot to disk.

ggsave(filename = "parameter_ice_vol_loss_fits.pdf",
       plot.fit1.fit2,
       width = 160,
       height = 200,
       units = "mm")

# Make predictions for glacier change with lake volume change.

vol_change_log <- seq(min(lakes.gl.loss$vol_change_log, na.rm = T), 
                      max(lakes.gl.loss$vol_change_log, na.rm = T), 
                      length.out = 100)

pred.grid.ice_loss_log <- add_predicted_draws(
  object = fit2, 
  newdata = data.frame(vol_change_log = vol_change_log),
  value = "ice_loss_log", 
  ndraws = 4000,
  re_formula = NA) 

for (i in unique(pred.grid.ice_loss_log$.row)) {
  
  row <- pred.grid.ice_loss_log$ice_loss_log[pred.grid.ice_loss_log$.row == i]
  
  HDIrange <- HDIofMCMC(row, 
                        credMass = 0.99)
  
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
  
  pred.grid.ice_loss_log$ice_loss_log[pred.grid.ice_loss_log$.row == i] <- row
  
}

# Plot the regional change in lake volume with glacier volume.
# Posterior trend is shown as a shade, with median trend on top of it, and
# point size scaled to the uncertainty in lake volume change.

lake.change.plot <- ggplot(lakes.gl.loss,
                           aes(x = qvolmd_km3,
                               y = abs(dvol)/10^9,
                               size = uncertainty)) +
  stat_lineribbon(data = pred.grid.ice_loss_log,
                  mapping = aes(x = 10^vol_change_log,
                                y = (10^ice_loss_log)/10^9), 
                  .width = 0.68,
                  point_interval = median_hdi,
                  fill = "grey66", 
                  inherit.aes = F) +
  geom_linerange(aes(ymin = (abs(dvol)/10^9) - (err_dvol/10^9),
                     ymax = (abs(dvol)/10^9) + (err_dvol/10^9)),
                 color = "mediumorchid4",
                 linewidth = 1) +
  geom_point(fill = "mediumorchid3", color = "white", shape = 21, stroke = 0.3) +
  scale_size_continuous(
    name = "Width of 68% HDI of lake\nvolume change [km設",
    range = c(1, 10)) +
  geom_text_repel(aes(label = full_name),
                  size = 7/.pt,
                  max.overlaps = 25) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 3),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 3),
                labels = trans_format("log10", math_format(10 ^ .x)))  +
  labs(x = "Regional gain in lake volume [km設 (1990 - 2020)",
       y = "Regional loss in glacier\nvolume [km設 (2000 - 2020)") +
  theme_bw() +
  theme(text = element_text(size = 7),
        aspect.ratio = 0.8,
        legend.position = "bottom")

# Combine both plots (lake volume and glacier volume; lake volume change and
# glacier volume change).

lake.gl.vol.change <- plot_grid(lake.vol.plot,
                                lake.change.plot , 
                                ncol = 2,
                                labels = c('a', 'b'),
                                label_size = 8)

ggsave(filename = "lake_and_glacier_volume_change.pdf",
       plot = lake.gl.vol.change, 
       height = 90,
       width = 180,
       units = "mm")

################################################################################

# Obtain the cumulated volume with elevation and distance to the ocean/ sink.

# We use the median elevation for every lake extracted from the COP30 elevation 
# dataset. Lakes are sorted by elevation in ascending order.

global.elevation <- lakes.w.stats %>% 
  filter(Year_new == 2020) %>%
  st_drop_geometry() %>%
  mutate(ElevationCOP = case_when(is.na(ElevationCOP) ~ Elevation,
                                  TRUE ~ ElevationCOP),
         ElevationCOP = case_when(ElevationCOP < 1 ~ 1,
                                  TRUE ~ ElevationCOP)) %>%
  arrange(ElevationCOP) 

# We now calculate the cumulative glacial lake volume for all predicted lake 
# volumes.

elevation.cumsum <- global.elevation %>% 
  select(starts_with(("Pred_"))) %>% 
  apply(., 2, function(x) {
    csum <- cumsum(x)
    
    return(csum/max(csum))
    
  })

# and extract the median and 68% HDI along the cumulative distribution.

vol.with.elev <- apply(elevation.cumsum, 1, 
                       function (x) {
                         tibble(lo = HDIofMCMC(x, 0.68)[1],
                                md = median(x), 
                                hi = HDIofMCMC(x, 0.68)[2])})

ve <- bind_rows(vol.with.elev) %>%
  mutate(x = global.elevation$ElevationCOP,
         var = "Elevation")

# Read flowpath data extracted from 90-m Merit DEM from Wolfgang Schwanghart.
# Data are a irregular txt file organised by lines. Every lake is identified 
# via its ID, and the X, Y, distance, and Z coordinate from the source towards
# the coast. So, after 5 lines, a new flow path starts.

fread <- read_lines("D:/nrc_user/veh/Zhang_glacial_lakes_global/flowpaths_detailed.txt")

idx <- seq(1, (nrow(lakes.w.stats %>% filter(Year_new == 2020))-1) *5, by = 5)

fp.list <- list()

for(i in idx) {
  
  # We only extract the last entry, which is the entails the cumulative flow path
  # length from the source.
  
  last.row <-  tibble(UniqueID = as.numeric(fread[i]),
                      X = fread[i+1] %>% str_split(" ") %>% unlist() %>% as.numeric() %>% na.omit(),
                      Y = fread[i+2] %>% str_split(" ") %>% unlist() %>% as.numeric() %>% na.omit(),
                      dist = fread[i+3] %>% str_split(" ") %>% unlist() %>% as.numeric() %>% na.omit(),
                      Z = fread[i+4] %>% str_split(" ") %>% unlist() %>% as.numeric() %>% na.omit()) %>%
    slice(n())
  
  fp.list[[i]] <- last.row
  
}

# Bind all flow paths to one object.

flowpaths <- fp.list[!sapply(fp.list, is.null)] %>% 
  do.call(rbind, .)

# How many lakes drain into endorheic basins?

flowpaths %>% 
  filter(Z > 500)

saveRDS(flowpaths, "flowpath_length.RDS")
# flowpaths <- readRDS("flowpath_length.RDS")

# We join the length of the flow paths with the lake volumes

lakes.w.flowpath <- left_join(lakes.w.stats, flowpaths, by = "UniqueID") 

# We sort the lake volumess by flow path distance.

global.flowpath <- lakes.w.flowpath %>% 
  filter(Year_new == 2020) %>%
  st_drop_geometry() %>%
  arrange(dist) 

# ... and calculate the cumulated lake volume with flow path distance.

flow.cumsum <- global.flowpath %>% 
  select(starts_with(("Pred_"))) %>%
  apply(., 2, function(x) {
    csum <- cumsum(x)
    
    return(csum/max(csum))
    
  })

# Take the median and 68% HDI of the cumulated volume with distance to the coast.

vol.with.dist <- apply(flow.cumsum, 1, 
                       function (x) {
                         tibble(lo = HDIofMCMC(x, 0.68)[1],
                                md = median(x), 
                                hi = HDIofMCMC(x, 0.68)[2])})

vd <- bind_rows(vol.with.dist) %>%
  mutate(x = global.flowpath$dist,
         var = "Distance upstream")

# We combine both the cumulated data (elevation and distance to coast) to
# show them in one plot

ve.vd <- bind_rows(ve, vd)

# What is the median elevation/ distance to coast?

median_data <- ve.vd %>%
  group_by(var) %>%
  filter(md >0.5) %>% 
  slice_min(md) %>%
  mutate(Z_label = paste0("italic(q)[50] == ", 
                          round(x, digits = 0), 
                          if_else(var == "Distance upstream", "~km", "~m"))) 

# Plot the cumulative distributions.

plot.cumul.dist.vol <- ggplot(data = ve.vd,
                              mapping = aes(x = x, 
                                            y = md,
                                            fill = var,
                                            color = var)) +
  geom_lineribbon(mapping = aes(ymin = lo,
                                ymax = hi)) +
  scale_color_manual(values = c("chocolate4", "mediumorchid4")) +
  scale_fill_manual(values = c("chocolate1", "mediumorchid1")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 5),
                limits = c(1, 10000),
                sec.axis = dup_axis(
                  name = "Flowpath distance [km] from\nthe lake to the ocean/sink", 
                  labels = function(x) round(x, digits = 0))) +
  geom_text(data = median_data,
            aes(x = x, y = 0.5, label = Z_label),
            parse = TRUE,
            size = 7/.pt,
            hjust = 0,
            nudge_x = 0.2,
            color = c("chocolate4", "mediumorchid4"),
            fontface = "italic") +
  theme_bw() +
  labs(x = "Elevation [m a.s.l.]",
       y = "Cumulated glacier lake volume")+ 
  theme(text = element_text(size = 7), 
        axis.text.x.bottom = element_text(color = "mediumorchid4"), 
        axis.title.x.bottom = element_text(color = "mediumorchid4"),
        axis.text.x.top = element_text(color = "chocolate4"),
        axis.title.x.top = element_text(color = "chocolate4"),
        legend.position = "none",
        aspect.ratio = 0.8)

# Save this plot to disk (Fig. 3)

ggsave(filename = "cumulative_distances_and_elevations.pdf", 
       plot = plot.cumul.dist.vol,
       width = 90,
       height = 90,
       units = "mm")

###############################################################################

# Regional cumulative lake volume with elevation. Generate a list of tibbles
# for each region, and sort the lake cumulatively by elevation.

regional.elevation <- lakes.w.stats %>% 
  filter(Year_new == 2020) %>%
  st_drop_geometry() %>%
  left_join(., o1regions %>% st_drop_geometry(), 
            by = "o1region", 
            relationship = "many-to-many") %>%
  mutate(ElevationCOP = case_when(is.na(ElevationCOP) ~ Elevation,
                                  TRUE ~ ElevationCOP),
         ElevationCOP = case_when(ElevationCOP < 1 ~ 1,
                                  TRUE ~ ElevationCOP)) %>%
  group_by(full_name.y) %>%
  group_split()

# Same approach as above. We only iterate over the list of regions, instead of making
# a global estimate.

regional.cumsum.z <- lapply(regional.elevation, function (x) { 
  
  x1 <- x %>%  
    arrange(ElevationCOP) 
  
  x2 <- x1 %>%  
    select(starts_with(("Pred_"))) %>% 
    apply(., 2, function(x) {
      csum <- cumsum(x)
      
      return(csum/max(csum))
      
    })
  
  
  x3 <- apply(x2, 1, 
              function (x) {
                
                hdi <- HDIofMCMC(x, 0.68)
                
                return(tibble(lo = hdi[1],
                              md = median(x),
                              hi = hdi[2]))})
  
  regional.cumsum <- bind_rows(x3) %>%
    mutate(Z = x1$ElevationCOP,
           var = "Elevation",
           region = x1$full_name.y)
  
  return(regional.cumsum)
})


ve.regional <- bind_rows(regional.cumsum.z) 

# Regional flow path distance

regional.flowpath <- lakes.w.flowpath %>% 
  filter(Year_new == 2020) %>%
  st_drop_geometry() %>%
  left_join(., o1regions %>% st_drop_geometry(), 
            by = "o1region", 
            relationship = "many-to-many") %>%
  group_by(full_name.y) %>%
  group_split()

# We also sort the regional lake volumes by flow path distance.

regional.cumsum.fp <- lapply(regional.flowpath, function (x) { 
  
  x1 <- x %>%  
    arrange(dist) 
  
  x2 <- x1 %>%  
    select(starts_with(("Pred_"))) %>% 
    apply(., 2, function(x) {
      csum <- cumsum(x)
      
      return(csum/max(csum))
      
    })
  
  
  x3 <- apply(x2, 1, 
              function (x) {
                
                hdi <- HDIofMCMC(x, 0.68)
                
                return(tibble(lo = hdi[1],
                              md = median(x),
                              hi = hdi[2]))})
  
  regional.cumsum <- bind_rows(x3) %>%
    mutate(dist = x1$dist,
           var = "Distance upstream",
           region = x1$full_name.y)
  
  return(regional.cumsum)
  
})

vd.regional <- bind_rows(regional.cumsum.fp) 

# Combine both tibbles.

ve.vd.regional <- bind_rows(ve.regional %>% rename("x" = "Z"), 
                            vd.regional %>% rename("x" = "dist"))

median_data.regional <- ve.vd.regional %>%
  group_by(var, region) %>%
  filter(md >0.5) %>% 
  slice_min(md) %>%
  mutate(Z_label = paste0("italic(q)[50] == ", 
                          round(x, digits = 0), 
                          if_else(var == "Distance upstream", "~km", "~m"))) 

# Generate a facet plot by region showing the cumulated volume with elevation and 
# flow path distance to the coast. 

plot.cumul.dist.vol.regional <- ggplot(data = ve.vd.regional,
                                       mapping = aes(x = x, 
                                                     y = md,
                                                     fill = var,
                                                     color = var)) +
  geom_lineribbon(mapping = aes(ymin = lo,
                                ymax = hi)) +
  scale_color_manual(values = c("chocolate4", "mediumorchid4")) +
  scale_fill_manual(values = c("chocolate1", "mediumorchid1")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 5),
                limits = c(1, 10000),
                sec.axis = dup_axis(
                  name = "Flowpath distance [km] from\nthe lake to the ocean/sink", 
                  labels = function(x) round(x, digits = 0))) +
  geom_text(data = median_data.regional %>% filter(var == "Distance upstream"),
            aes(x = 1, y = 0.9, label = Z_label),
            parse = TRUE,
            size = 7/.pt,
            hjust = 0,
            color = c("mediumorchid4"),
            fontface = "italic") +
  geom_text(data = median_data.regional %>% filter(var == "Elevation"),
            aes(x = 1, y = 0.75, label = Z_label),
            parse = TRUE,
            size = 7/.pt,
            hjust = 0,
            color = c("chocolate4"),
            fontface = "italic") +
  theme_bw() +
  labs(x = "Elevation [m a.s.l.]",
       y = "Cumulated glacier lake volume") +
  facet_wrap(~region, ncol = 4) +
  theme(text = element_text(size = 7), 
        axis.text.x.bottom = element_text(color = "mediumorchid4"), 
        axis.title.x.bottom = element_text(color = "mediumorchid4"),
        axis.text.x.top = element_text(color = "chocolate4"),
        axis.title.x.top = element_text(color = "chocolate4"),
        strip.background = element_rect(fill="white"),
        legend.position = "none")

# Save this plot to disk (Figure S6)

ggsave(filename = "cumulative_distances_and_elevations_regional.pdf", 
       plot = plot.cumul.dist.vol.regional,
       width = 160,
       height = 200,
       units = "mm")