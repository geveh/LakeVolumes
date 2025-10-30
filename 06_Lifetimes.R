### Estimate sedimentation-driven lifespan of glacial lakes on local,
### regional, and global level.
### Assess land cover in contributing catchments of glacial lakes.

## Load packages that we'll need in the course of this script

# If packages haven't been installed beforehand, 
# first use install.packages("packagename")

require(sf)
require(tidyverse)
require(scales)
require(brms)
require(tidybayes)
require(ggrepel)
require(cowplot)
require(readxl)
require(rnaturalearthdata)
require(rnaturalearth)
require(fitdistrplus)
require(doParallel)

# Set working directory

setwd("D:/nrc_user/veh/Zhang_glacial_lakes_global/")

# Add the function HDIofMCMC from John Kruschke to this directory. We will use
# this function to derive the highest density intervals of posterior 
# (predictive) distributions of lake volumes.

source("HDIofMCMC.R")

# Load the O1 regions from the Randolph Glacier Inventory V7.0.

o1regions <- st_read("D:/nrc_user/veh/RGI2000-v7.0-regions/RGI2000-v7.0-o1regions.shp")
o1regions <- o1regions[!duplicated(o1regions$o1region), ]

ymin <- rep(NA, nrow(o1regions))

for(i in 1:nrow(o1regions)) {ymin[i] <- st_bbox(o1regions[i, ])["ymin"] }

o1regions$Y <- ymin

# Load lake data with predicted volumes

all.lakes <- readRDS("all_lakes_with_volumes.RDS")

# Catchment stats for lakes

lake.stats.sf <- readRDS("lakes_with_catchment_and_landcover.rds")

# Define the join key to join predicted lake volumes with catchments stats 
# via the common key UniqueID.

join_key <- "UniqueID"

# Determine the columns in df_right that are not in df_left (excluding the join key)
new_vars <- setdiff(names(lake.stats.sf), c(names(all.lakes), join_key))

# Perform the left join by including only the join key and the new variables all.lakes

lakes.w.stats <- all.lakes %>%
  left_join(., dplyr::select(lake.stats.sf, all_of(c(join_key, new_vars))), by = join_key)

# Remove supraglacial lakes, as they are too volatile to be permanent features
# for infill

lakes.2020.no.sup <- lakes.w.stats %>% 
  filter(Year_new == 2020, 
         Lake_type != "supraglacial", 
         Lake_type != "ice") %>%
  filter(!is.na(glaciercov)) %>%
  dplyr::left_join(., o1regions %>% st_drop_geometry(), 
            by = "o1region", 
            relationship = "many-to-one")

# Load glacial erosion rates compiled in the paper by Wilner et al. (2024)
# Data S1 and S2 can be downloaded from the supplement of the original paper: 
# https://www.science.org/doi/10.1126/sciadv.adr2009

d.glac <- read_excel("D:/nrc_user/veh/Erosion_rates/sciadv.adr2009_data_s1_and_s2/adr2009_data_s1.xlsx",
                     sheet = 1)

# We only use 'contemporary' erosion rates (time interval < 500 years)

d.glac.filt <- d.glac %>% filter(`Time interval (yr)` < 500) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326)

d.fluv <- read_excel("D:/nrc_user/veh/Erosion_rates/sciadv.adr2009_data_s1_and_s2/adr2009_data_s2.xlsx",
                     sheet = 1)

d.fluv.filt <- d.fluv %>% filter(`Time interval (yr)` < 500) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326)

# Fit a normal distributin to log10-transformed fluvial and glacial erosion rates

glac.norm <- fitdist(log10(d.glac.filt$`Erosion rate (mm/yr)`), "norm")
fluv.norm <- fitdist(log10(d.fluv.filt$`Erosion rate (mm/yr)`), "norm") 

# Generate the curve for the normal distribution

x_vals <- seq(-5, 5, length.out = 1000)
y.glac <- dnorm(x_vals, mean = glac.norm$estimate["mean"], sd = glac.norm$estimate["sd"])
y.fluv <- dnorm(x_vals, mean = fluv.norm$estimate["mean"], sd = fluv.norm$estimate["sd"])


plot(density(log10(d.glac.filt$`Erosion rate (mm/yr)`)), xlim = c(-5,5), col = "blue")
lines(density(log10(d.fluv.filt$`Erosion rate (mm/yr)`)), col = "darkorange")

# Add the normal distribution curve to the plot
lines(x_vals, y.glac, col = "blue", lwd = 2, lty = "dashed")
lines(x_vals, y.fluv, col = "darkorange", lwd = 2, lty = "dashed")

# Derive the HDI of predicted glacial lake volumes

lakes.hdi <- lakes.2020.no.sup %>% 
  st_drop_geometry() %>% 
  filter(!is.na(glaciercov)) %>%
  dplyr::select(starts_with("Pred_")) %>%
  data.matrix() %>%
  apply(., 1, function(x) { 
    hdi <- HDIofMCMC(x, 0.68)
    return(unname(x[(x > hdi[1]) & (x < hdi[2])]))
    })

# We model the expected lifetime of each lake until complete infill by simulating
# glacial and fluvial erosion rates, rock and sediment bulk densities, and trapping
# efficiencies, including their uncertainties. The outcome are 5000 simulated
# lifetimes, i.e. the period until complete infill.

sed.prod.sys <- list()
  
for(i in 1: 5000) {
  
  # From each lake, sample one predicted volume
  
  vols <- lakes.hdi %>% 
    sapply(., function(y) sample(y, size = 1, replace = F))
  
  # Multiply catchment area (in meters, so area x 10^6) with erosion rate (in
  # meters, so x 0.001), and rock density (2.6)
  
  # Simulate a glacial erosion rate for each lake and truncate extreme ones.
  
  erosion.rate.glac <- 10^rnorm(n = length(lakes.2020.no.sup$catchmnolake),
                           mean = glac.norm$estimate["mean"], 
                           sd = glac.norm$estimate["sd"])
  
  hdi.er.glac <- HDIofMCMC(erosion.rate.glac, 0.98)

  # Identify extreme values
  extreme.er.glac <- which(erosion.rate.glac < hdi.er.glac[1] | erosion.rate.glac > hdi.er.glac[2])
  
  if (length(extreme.er.glac) > 0) {

    erosion.rate.glac[extreme.er.glac] <- sample(erosion.rate.glac[erosion.rate.glac >= hdi.er.glac[1] & erosion.rate.glac <= hdi.er.glac[2]], 
                                        length(extreme.er.glac), 
                                        replace = TRUE)
  }
  
  # Obtain a fluvial erosion rate for each lake and truncate extreme ones.
  
  erosion.rate.fluv <- 10^rnorm(n = length(lakes.2020.no.sup$catchmnolake),
                                mean = fluv.norm$estimate["mean"], 
                                sd = fluv.norm$estimate["sd"])
  
  hdi.er.fluv <- HDIofMCMC(erosion.rate.fluv, 0.98)
  
  # Identify extreme values
  
  extreme.er.fluv <- which(erosion.rate.fluv < hdi.er.fluv[1] | erosion.rate.fluv > hdi.er.fluv[2])
  
  if (length(extreme.er.fluv) > 0) {

    erosion.rate.fluv[extreme.er.fluv] <- sample(erosion.rate.fluv[erosion.rate.fluv >= hdi.er.fluv[1] & erosion.rate.fluv <= hdi.er.fluv[2]], 
                                                 length(extreme.er.fluv), 
                                                 replace = TRUE)
  }
  
  # Simulate a rock density in the catchment feeding the lake
  
  rock.density <- rnorm(n = length(lakes.2020.no.sup$catchmnolake),
                        mean = 2.6,
                        sd = 0.2)
  
  hdi.rd <- HDIofMCMC(rock.density, 0.98)
  
  # Identify extreme values
  
  extreme_rd <- which(rock.density < hdi.rd[1] | rock.density > hdi.rd[2])
  
  if (length(extreme_rd) > 0) {
  
    rock.density[extreme_rd] <- sample(rock.density[rock.density >= hdi.rd[1] & rock.density <= hdi.rd[2]], 
                                       length(extreme_rd), 
                                       replace = TRUE)
  }
  
  # Correct for changing density (i.e. compaction) from catchment-wide erosion 
  # to deposition in a given lake
  
  deposition.density <- rnorm(n = length(lakes.2020.no.sup$catchmnolake),
                              mean = 1.6,
                              sd = 0.2) 
  
  hdi.dd <- HDIofMCMC(deposition.density, 0.98)
  
  # Identify extreme values
  
  extreme_dd <- which(deposition.density < hdi.dd[1] | deposition.density > hdi.dd[2])
  
  if (length(extreme_dd) > 0) {
  
    deposition.density[extreme_dd] <- sample(deposition.density[deposition.density >= hdi.dd[1] & deposition.density <= hdi.dd[2]], 
                                       length(extreme_dd), 
                                       replace = TRUE)
  }
  
  
  # Simulate a trapping efficiency for each lake. We assume 
  # a beta distribution with mode 0.8 to reflect that trapping efficiency is 
  # usually very high (~80%)
  
  trap.eff <- rbeta(length(lakes.2020.no.sup$catchmnolake), shape1 = 9, shape2 = 3)
  
  # Multiply the catchment area (in sq meters) with the erosion rate
  # and the rock density to determine the annual sediment production.
  
  ann.sediment.prod.glac <- lakes.2020.no.sup$catchmnolake * 10^6 * 
    (0.001 * erosion.rate.glac) * rock.density
  
  # Divide production rate by deposition density
  
  ann.sed.depo.glac <- ann.sediment.prod.glac / deposition.density

  # The time to infill is the initially available lake volume (in cubic meters) 
  # as of 2020, divided by the amount that is deposited every year, and 
  # corrected for a trapping efficiency that is lower than 100%. 
  
  time.to.fill.glac <- (vols*10^6) / (ann.sed.depo.glac * trap.eff )
  
  # Same for fluvial
  
  ann.sediment.prod.fluv <- lakes.2020.no.sup$catchmnolake * 10^6 * 
    (0.001 * erosion.rate.fluv) * rock.density
  
  ann.sed.depo.fluv <- ann.sediment.prod.fluv / deposition.density
  
  time.to.fill.fluv <- (vols*10^6) / (ann.sed.depo.fluv * trap.eff )
  
  # Weighted by glacier cover: 
  
  time.to.fill.weighted <- ifelse(runif(length(lakes.2020.no.sup$glaciercov)) < lakes.2020.no.sup$glaciercov,
                   time.to.fill.fluv, time.to.fill.glac)
  
  sed.prod.sys[[i]] <- tibble(draw = i,
                              glac_er = time.to.fill.glac,
                              fluv_er = time.to.fill.fluv,
                              weighted_er = time.to.fill.weighted)
  
} 

saveRDS(sed.prod.sys, "simulated_infill_times.RDS")
# sed.prod.sys <- readRDS("simulated_infill_times.RDS")

lifetime.all <- bind_rows(sed.prod.sys)

################################################################################

# Modelling the year until a given percentage of the initial lake volume (as of
# 2020) will be lost. Here we do not provide an estimate weighted by glacial
# cover. The code itself is similar to that above and left largely without
# comments.

fill.lakes <- lakes.2020.no.sup %>% 
  st_drop_geometry() %>% 
  filter(!is.na(glaciercov))

registerDoParallel(20)

# For each lake, we make 1500 simulations to approximate the year until complete
# infill under either fluvial or glacial erosion rates.

sims <- 1500

year.to.fill <- foreach(i = 1:nrow(fill.lakes),
        .packages = c("dplyr")) %dopar% {
  
  # First, take a specific lakes and sample predicted volumes from the 68% HDI.        
          
  vol.fill.lakes <- fill.lakes[i, ] %>%
    dplyr::select(starts_with("Pred_")) %>%
      data.matrix() %>%
      apply(., 1, function(x) {
        hdi <- HDIofMCMC(x, 0.98)
        return(unname(x[(x > hdi[1]) & (x < hdi[2])]))
      })

  # Convert the initial lake volume to cubic meters.
  
  vols <- sample(vol.fill.lakes, size = sims, replace = T) * 10^6
    
  # Multiply catchment area (in meters, so area x 10^6) with erosion rate (in
  # meters, so x 0.001), and rock density (2.6)
  
  # Obtain glacier erosion rates and truncate extreme ones.
  
  erosion.rate.glac <- 10^rnorm(n = sims,
                                mean = glac.norm$estimate["mean"], 
                                sd =   glac.norm$estimate["sd"]) * 0.001
    
  hdi.er.glac <- HDIofMCMC(erosion.rate.glac, 0.98)
    
  # Identify extreme values
 
  extreme.er.glac <- which(erosion.rate.glac < hdi.er.glac[1] | erosion.rate.glac > hdi.er.glac[2])
  
  if (length(extreme.er.glac) > 0) {
    
    erosion.rate.glac[extreme.er.glac] <- sample(erosion.rate.glac[erosion.rate.glac >= hdi.er.glac[1] & erosion.rate.glac <= hdi.er.glac[2]], 
                                                 length(extreme.er.glac), 
                                                 replace = TRUE)
  }
  
  # Obtain fluvial erosion rates and truncate extreme ones.
  
  erosion.rate.fluv <- 10^rnorm(n = sims,
                                mean = fluv.norm$estimate["mean"], 
                                sd =   fluv.norm$estimate["sd"]) * 0.001
  
  hdi.er.fluv <- HDIofMCMC(erosion.rate.fluv, 0.98)
  
  # Identify extreme values
  
  extreme.er.fluv <- which(erosion.rate.fluv < hdi.er.fluv[1] | erosion.rate.fluv > hdi.er.fluv[2])
  
  if (length(extreme.er.fluv) > 0) {
    
    erosion.rate.fluv[extreme.er.fluv] <- sample(erosion.rate.fluv[erosion.rate.fluv >= hdi.er.fluv[1] & erosion.rate.fluv <= hdi.er.fluv[2]], 
                                                 length(extreme.er.fluv), 
                                                 replace = TRUE)
  }
  
  # Rock densities
  
  rock.density <- rnorm(n = sims,
                        mean = 2.6,
                        sd = 0.2)
  
  hdi.rd <- HDIofMCMC(rock.density, 0.98)
  
  # Identify extreme values
  
  extreme_rd <- which(rock.density < hdi.rd[1] | rock.density > hdi.rd[2])
  
  if (length(extreme_rd) > 0) {
    
    rock.density[extreme_rd] <- sample(rock.density[rock.density >= hdi.rd[1] & rock.density <= hdi.rd[2]], 
                                       length(extreme_rd), 
                                       replace = TRUE)
  }
  
  # Correct for changing density (i.e. compaction) from erosion to deposition
  
  deposition.density <- rnorm(n = sims,
                              mean = 1.6,
                              sd = 0.2) 
  
  hdi.dd <- HDIofMCMC(deposition.density, 0.98)
  
  # Identify extreme values
  
  extreme_dd <- which(deposition.density < hdi.dd[1] | deposition.density > hdi.dd[2])
  
  if (length(extreme_dd) > 0) {

    deposition.density[extreme_dd] <- sample(deposition.density[deposition.density >= hdi.dd[1] & deposition.density <= hdi.dd[2]], 
                                             length(extreme_dd), 
                                             replace = TRUE)
  }
  
  # Trapping efficiency according to beta distribution with mode 0.8
  
  trap.eff <- rbeta(sims, shape1 = 9, shape2 = 3)
  
  # Multiply the catchment area (in sq meters) with the erosion rate
  # and the rock density to determine the annual sediment production.
  
  ann.sediment.prod.glac <- fill.lakes$catchmnolake[i] * 10^6 * 
    erosion.rate.glac * rock.density
  
  ann.sed.depo.glac <- ann.sediment.prod.glac / deposition.density
  
  time.to.fill.glac <- vols / (ann.sed.depo.glac * trap.eff)

  # Same for fluvial
  
  ann.sediment.prod.fluv <- fill.lakes$catchmnolake[i] * 10^6 * 
     erosion.rate.fluv * rock.density
  
  ann.sed.depo.fluv <- ann.sediment.prod.fluv / deposition.density
  
  time.to.fill.fluv <- vols / (ann.sed.depo.fluv * trap.eff)

  filled.lake <- tibble(UniqueID = fill.lakes$UniqueID[i],
                        full_name = fill.lakes$full_name.y[i],
                        Area = fill.lakes$Area[i],
                        glac_lt =  time.to.fill.glac,
                        fluv_lt =  time.to.fill.fluv)
  
  return(filled.lake)
  
  } 

# Combine these estimates until complete infill to one large data frame.
# We assign each lake to a given size class.

ytf <- do.call(rbind, year.to.fill) %>%
  mutate(interval = cut(Area, breaks = c(0, 0.1, 1, 10, 1000)))

# Estimate the year when 10, 25, 50, 75 and 100% of the initial lake volume
# will be lost.

# For glacial erosion rates...

remaining.storage.glac <- lapply(c(0.10, 0.25, 0.50, 0.75, 1), function(x) {
  
  o <- ytf %>%
   mutate(perc_lost = (ytf$glac_lt* x) + 2020) %>%
   group_by(interval) %>%
   dplyr::summarise(lo = HDIofMCMC(perc_lost, 0.68)[1],
                    md = quantile(perc_lost, 0.5),
                    hi = HDIofMCMC(perc_lost, 0.68)[2]) %>%
   mutate(removed = x)
  }) %>% do.call(rbind, .) %>%
  dplyr::mutate(
    interval = factor(interval, levels = unique(interval)),  # preserve order
    removed = factor(paste0(removed*100, "%")),
    removed = factor(removed, levels = c("10%", "25%", "50%", "75%", "100%")),
    er = "Glacial erosion rate"
  )

# ... and for fluvial erosion rates

remaining.storage.fluv <- lapply(c(0.10, 0.25, 0.50, 0.75, 1), function(x) {
  
  o <- ytf %>%
    mutate(perc_lost = (ytf$fluv_lt* x) + 2020) %>%
    group_by(interval) %>%
    dplyr::summarise(lo = HDIofMCMC(perc_lost, 0.68)[1],
                     md = quantile(perc_lost, 0.5),
                     hi = HDIofMCMC(perc_lost, 0.68)[2]) %>%
    dplyr:: mutate(removed = x)
  
  } ) %>% do.call(rbind, .)  %>%
  dplyr:: mutate(
    interval = factor(interval, levels = unique(interval)),  
    removed = factor(paste0(removed*100, "%")),
    removed = factor(removed, levels = c("10%", "25%", "50%", "75%", "100%")),
    er = "Fluvial erosion rate"
  )

# Combine infill periods from both glacial and fluvial erosion estimates

remaining.storage <- rbind(remaining.storage.glac,
                           remaining.storage.fluv) %>%
  dplyr::mutate(er = factor(er, levels = c("Glacial erosion rate",
                                           "Fluvial erosion rate")))

# Plot the year until a given lake storage capacity (10, 25, 50, 75, 100%) has 
# been infilled under these two erosion scenarios as a faceted plot.

storage.plot <- ggplot(remaining.storage, aes(x = interval, 
                                              y = md, 
                                              color = removed, 
                                              group = removed)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  geom_linerange(aes(ymin = lo, ymax = hi),
                 position = position_dodge(width = 0.7), size = 0.8) +
  labs(x = "Lake area [km²]", 
       y = "Year of capacity loss (median and 68% HDI)",
       color = "Capacity loss") +
  geom_text(aes(label = round(md), 
                y = hi + (hi*0.2),
                color = removed), 
            position = position_dodge(width = 0.7),
            vjust = 0.5, 
            hjust = 0, 
            size = 7/.pt,
            angle = 90) +
  labs(x = "Lake area (km²)", y = "Volume (km³)", fill = "Year") +
  facet_wrap(~er) +
  scale_colour_viridis_d() +
  scale_y_log10(
    breaks = c(2000,    2100,   2300,   2500,   5000,   10000,   20000,   50000,   100000),
    labels = c("2000", "2100", "2300", "2500", "5000", "10000", "20000", "50000", "100000")
  ) +
  theme_bw() +
  theme(text = element_text(size = 7), 
        strip.background = element_blank())

ggsave(filename = "remaining_storage.pdf",
       plot = storage.plot, 
       width = 140,
       height = 90,
       units = "mm",
       device = cairo_pdf, 
       family = "Arial")

################################################################################

# What is the median volume of each lake?

med.vols <- lakes.2020.no.sup %>% 
  st_drop_geometry() %>% 
  dplyr::filter(!is.na(glaciercov)) %>%
  dplyr::select(starts_with("Pred_")) %>%
  data.matrix() %>% 
  apply(., 1, median)

large.idx <- sort.int(med.vols, decreasing = T, index.return = T)$ix[1:40]

# What is the lifetime of the 40 largest lake?

lifetime.largest <- lapply(sed.prod.sys, function (x) {
  
  x[large.idx, ] %>% 
    mutate(ID = 1:40)
  
})

lifetime.largest <- bind_rows(lifetime.largest) %>%
  pivot_longer(cols = c("glac_er", "fluv_er", "weighted_er")) %>%
  group_by(ID, name) %>%
  dplyr::summarise(md = median(value)) 

lifetime.largest.glac <- lifetime.largest %>% 
  filter(name == "glac_er") 
lifetime.largest.fluv <-lifetime.largest %>% 
  filter(name == "fluv_er") 

# We also generate a plot showing the probability densities of stacked individual 
# lifetimes, distinguished by glacial and fluvial erosion rates, and an estimate 
# weighted by glacier cover.

min.pr <- min(log10(lifetime.all))
max.pr <- max(log10(lifetime.all))

# The density is derived using a fixed bandwidth for each of the 1500 
# simulated lifetimes per lake and each erosion scenario.

dens.2020 <- lapply(sed.prod.sys, function (x) {
                     
                    h1 <- density(log10(x$glac_er), bw = 0.125, from = min.pr, to = max.pr)
                    h2 <- density(log10(x$fluv_er), bw = 0.125, from = min.pr, to = max.pr)
                    h3 <- density(log10(x$weighted_er), bw = 0.125, from = min.pr, to = max.pr)
                     
                    densities <-  tibble(dens_glac = h1$y,
                                         dens_fluv = h2$y,
                                         dens_weighted = h3$y)
                     
                    return(densities)
                    
                   })

brea <- density(log10(sed.prod.sys[[1]]$glac_er), bw = 0.125,
                from = min.pr, 
                to = max.pr)$x

# The uncertainty in the density is also summarised by the 68% HDI.

# For the glacial erosion rate scenario

dens.summary.glac <- sapply(dens.2020, function(x) x$dens_glac) %>% 
  apply(., MARGIN = 1, function(y) {
    
    hdi <- HDIofMCMC(y, 0.68)
    md = median(y)
    
    tibble(lo = hdi[1],
           md = md,
           hi = hdi[2],
           er = "glacial")
  }) %>% bind_rows()

# the fluvial erosion rate scenario

dens.summary.fluv <- sapply(dens.2020, function(x) x$dens_fluv) %>% 
  apply(., MARGIN = 1, function(y) {
    
    hdi <- HDIofMCMC(y, 0.68)
    md = median(y)
    
    tibble(lo = hdi[1],
           md = md,
           hi = hdi[2],
           er = "fluvial")
  }) %>% bind_rows()

# and the erosion rate scenario weighted by glacial cover in the feeding 
# catchment.

dens.summary.weighted <- sapply(dens.2020, function(x) x$dens_weighted) %>% 
  apply(., MARGIN = 1, function(y) {
    
    hdi <- HDIofMCMC(y, 0.68)
    md = median(y)
    
    tibble(lo = hdi[1],
           md = md,
           hi = hdi[2],
           er = "weighted")
  }) %>% bind_rows()

# Combine all summarised densities of stacked individual lifetimes.

dens_summary <- bind_rows(dens.summary.glac,
                          dens.summary.fluv,
                          dens.summary.weighted) %>%
  mutate(brea = rep(brea, 3))

# We also want to know when a given fraction of the initial lake volume will 
# be lost. To this end, we obtain the empirical quantiles of modelled infill
# times in 1%-steps.

er.quant <- lapply(sed.prod.sys, function (x) {
  x %>%
   reframe(glacial = quantile(glac_er, seq(0, 1, by = 0.01)),
           fluvial = quantile(fluv_er, seq(0, 1, by = 0.01)),
           weighted = quantile(weighted_er, seq(0, 1, by = 0.01))) %>%
    mutate(q = seq(0, 1, by = 0.01))})

er.cumulative.summary <- bind_rows(er.quant) %>%
  pivot_longer(cols = c("glacial", "fluvial", "weighted"),
               names_to = "er") %>%
  dplyr::rename("brea" = "q") %>%
  group_by(brea, er) %>%
  dplyr::summarise(md = median(value),
            lo = HDIofMCMC(value, 0.68)[1],
            hi = HDIofMCMC(value, 0.68)[2])


med.lifetime.glacial <- lifetime.all %>% 
  summarise(md_glac = round(median(glac_er)),
            lo = round(md_glac-HDIofMCMC(glac_er, 0.68)[1]),
            hi = round(HDIofMCMC(glac_er, 0.68)[2]-md_glac),
            label_glac = paste0("glacial :",scales::comma(md_glac),
                                "(+", 
                                scales::comma(hi),
                                "/-", 
                                scales::comma(lo),
                                ") yr"))

med.lifetime.fluvial <- lifetime.all %>% 
   summarise(md_fluv = round(median(fluv_er)),
             lo = round(md_fluv- HDIofMCMC(fluv_er, 0.68)[1]),
             hi = round(HDIofMCMC(fluv_er, 0.68)[2]-md_fluv),
             label_fluv = paste0("fluvial :", scales::comma(md_fluv),
                                 "(+", 
                                 scales::comma(hi),
                                 "/-", 
                                 scales::comma(lo),
                                 ") yr"))

med.lifetime.weighted <- lifetime.all %>% 
  summarise(md_weighted = round(median(weighted_er)),
            lo = round(md_weighted - HDIofMCMC(weighted_er, 0.68)[1]),
            hi = round(HDIofMCMC(weighted_er, 0.68)[2]-md_weighted),
            label_weighted = paste0("weighted :", scales::comma(md_weighted),
                               "(+", 
                               scales::comma(hi),
                               "/-", 
                               scales::comma(lo),
                               ") yr"))

# Plot the stacked probability densities of projected sedimentation-driven 
# lifetimes of all glacial lakes as of 2020.

density.lifetime.plot <- ggplot(dens_summary, 
         mapping = aes(x = 10^brea,
                       y = md,
                       ymin = lo,
                       ymax = hi,
                       group = er, 
                       fill = er, 
                       color = er)) +
  geom_lineribbon(alpha = 0.5,
                  linewidth = 0) +
  geom_line(mapping = aes(x = 10^brea,
                          y = md,
                          group = er, 
                          color = er),
            linewidth = 0.8) +  # Line for each Year_new
  geom_text(data = med.lifetime.fluvial,
            mapping = aes(x = 10^-2.5,
                          y = 0.25,
                          label = label_fluv),
            color = "navy",
            size = 7/.pt, hjust = 0,  inherit.aes = F) +
  geom_text(data = med.lifetime.glacial,
            mapping = aes(x = 10^-2.5,
                          y = 0.35,
                          label = label_glac),
            color = "darkorange2",
            size = 7/.pt, hjust = 0,  inherit.aes = F) +
  geom_text(data = med.lifetime.weighted,
            mapping = aes(x = 10^-2.5,
                          y = 0.3,
                          label = label_weighted),
            color = "mediumvioletred",
            size = 7/.pt, hjust = 0, inherit.aes = F) +
  scale_fill_manual(name = "Erosion rate", 
                    values = c("lightblue", "darkorange", "maroon1")) + 
  scale_color_manual(name = "Erosion rate", 
                     values = c("navy", "darkorange2", "mediumvioletred")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 6),
                labels = trans_format("log10", math_format(10 ^ .x)),
                limits = c(10^-3, 10^8)) +
  labs(y = "Density",
       x = expression("Estimated life time [years]")) +
  theme_bw() +
  theme(text = element_text(size = 7),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank())

# Save this plot (Figure 4a).

ggsave(filename = "density_lifetime.pdf",
       plot = density.lifetime.plot, 
       width = 70,
       height = 90,
       units = "mm",
       device = cairo_pdf, 
       family = "Arial")

# Global lifetime of 10%, 50% and 90% of all glaciers

er.cumulative.summary %>% 
  filter(brea == 0.1)

er.cumulative.summary %>% 
  filter(brea == 0.5)

er.cumulative.summary %>% 
  filter(brea == 0.9)

HDIofMCMC(lifetime.all$weighted_er, 0.68)
median(lifetime.all$weighted_er)

################################################################################

# We also model the year when a given fraction, here in 1%-steps, of the 
# regionally available glacial lake volume could be entirely filled with sediments.

sed.prod.sys.reg <- lapply(sed.prod.sys, function (x) {
  y <- x %>% mutate(o1region = lakes.2020.no.sup$full_name.y)
  return(y)})

reg.lifetime <- sed.prod.sys.reg %>%
  bind_rows() %>%
  group_by(draw, o1region) %>%
  dplyr::reframe(glac = quantile(glac_er, seq(0, 1, 0.01)),
                 fluv = quantile(fluv_er, seq(0, 1, 0.01)),
                 weighted = quantile(weighted_er, seq(0, 1, 0.01)),
                 q = seq(0, 1, 0.01))

# We add labels that show the year when a a given regional fraction (10, 25 and   
# 50%) of glacial lakes will be theoretically filled with sediments.

reg.lifetime.summary <- reg.lifetime %>%
  pivot_longer(cols = c("glac", "fluv", "weighted")) %>%
  group_by(o1region, name, q) %>%
  dplyr::reframe(md = median(value),
          lo = HDIofMCMC(value, 0.68)[1],
          hi =  HDIofMCMC(value, 0.68)[2])

label.lifetime <- reg.lifetime.summary %>%
  filter(q == 0.1 | q == 0.25 | q == 0.5) %>%
  mutate(md = md + 2020) %>%
  dplyr::select(o1region, name, q, md) %>%
  pivot_wider(names_from = q, values_from = md, names_prefix = "q") %>%
  mutate(label = paste0(round(q0.1), " | ",  round(q0.25), " | ", round(q0.5)))

# We only show the trajectory of sedimentation-based infill until the year 
# 10000 AD.

plot.regional.lifetimes <- reg.lifetime.summary %>%
  filter(lo < 9780) %>%
  arrange(name, q) %>%  
  ggplot(aes(x = md + 2020,
             y = q, 
             color = name,
             fill = name)) +
  geom_lineribbon(mapping = aes(xmin = lo+2020,
                                xmax = hi+2020),
                  linewidth = 0.6) +
  scale_fill_manual(name = "Erosion rate", 
                    values = c("lightblue", "darkorange", "maroon1")) +  # Both ribbons grey
  scale_color_manual(name = "Erosion rate", 
                     values = c("navy", "darkorange2", "mediumvioletred")) +
  geom_text(data = label.lifetime %>% filter(name == "glac"),
            aes(x = 4000, y = 0.3, label = label),
            inherit.aes = FALSE,    # use annotation data only
            color = "darkorange2",
            size = 7 * 0.35,
            hjust = 0) +
  geom_text(data = label.lifetime %>% filter(name == "weighted"),
            aes(x = 4000, y = 0.17, label = label),
            inherit.aes = FALSE,    # use annotation data only
            color = "mediumvioletred",
            size = 7 * 0.35,
            hjust = 0) +
  geom_text(data = label.lifetime %>% filter(name == "fluv"),
            aes(x = 4000, y = 0.02, label = label),
            inherit.aes = FALSE,    # use annotation data only
            color = "navy",
            size = 7 * 0.35,
            hjust = 0) +
  geom_hline(yintercept = c(0.1, 0.25, 0.5),
             linetype = "dashed",
             color = "grey33",
             linewidth = 0.4) +
  labs(y = "Empirical cumulative distribution",
       x = "Estimated year [AD] of complete infill") +
  theme_bw() +
  xlim(c(2020, 10000)) +
  facet_wrap(~o1region, ncol = 4) +
  theme(text = element_text(size = 7),
        legend.position = "none",
        aspect.ratio = 1,
        strip.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank())  

# Save the regional lake lifetimes to disk (Figure S8).

ggsave("Regional_lake_lifetimes.pdf",
       plot.regional.lifetimes,
       width = 180,
       height = 190,
       units = "mm",
       device = cairo_pdf, 
       family = "Arial")

# Small plots for supporting figure

# Size distribution of contributing catchments

plot.catchm <- lakes.2020.no.sup %>% 
  ggplot(mapping = aes(x = catchmnolake)) +
  geom_density(fill = "lightblue") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 3),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  theme_bw() +
  labs(x = NULL,
       y = "Density") +
  theme(text = element_text(size = 7),
        aspect.ratio = 0.8)

# Figure S7a

ggsave("small_panel_catchment_size.pdf",
       plot.catchm,
       width = 45,
       height = 40,
       units = "mm")

# Erosion rates

erosion.rate.fluv <- 10^rnorm(n = 10^6,
                              mean = fluv.norm$estimate["mean"], 
                              sd = fluv.norm$estimate["sd"])

erosion.rate.glac <- 10^rnorm(n = 10^6,
                              mean = glac.norm$estimate["mean"], 
                              sd = glac.norm$estimate["sd"])

er.sim <- bind_rows(
  tibble(rate = erosion.rate.glac,
         Type = "glacial"),
  tibble(rate = erosion.rate.fluv,
         Type = "fluvial"))

er.emp <- bind_rows(
  tibble(rate = d.glac.filt$`Erosion rate (mm/yr)`,
         Type = "glacial"),
  tibble(rate = d.fluv.filt$`Erosion rate (mm/yr)`,
         Type = "fluvial"))
  
plot.er <- er.emp %>%  
  ggplot(mapping = aes(x = rate, fill = Type)) +
  geom_density( alpha=.5, color = "grey75") +
  scale_fill_manual(values = c("lightblue", "darkorange")) +
  geom_density(data = er.sim,
               mapping = aes(x = rate,
                             group = Type), 
               inherit.aes = F) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  theme_bw() +
  labs(x = NULL,
       y = "Density") +
  theme(text = element_text(size = 7),
        aspect.ratio = 0.8,
        legend.position = "none")

ggsave("small_panel_erosion_rate.pdf",
       plot.er,
       width = 45,
       height = 40,
       units = "mm")

# Plot empirical probability density of glacier cover in the contribution
# catchments of glacial lakes. 

plot.glaccov <- lakes.2020.no.sup %>% 
  ggplot(mapping = aes(x = glaciercov)) +
  geom_density(fill = "cyan") + 
  theme_bw() +
  labs(x = NULL,
       y = "Density")+
  theme(text = element_text(size = 7),
        aspect.ratio = 0.8,
        legend.position = "none")

# Fig. 7c

ggsave("small_panel_glac_cov.pdf",
       plot.glaccov,
       width = 45,
       height = 40,
       units = "mm")

# Rock density

x_vals <- seq(1.8, 3.4, length.out = 500)

df <- data.frame(x = x_vals,
                 y = dnorm(x_vals, mean = 2.6, sd = 0.2))

# Plot probability distribution of rock density in the contributing catchment

plot.rd <- ggplot(df, aes(x = x, y = y)) +
  geom_area(fill = "burlywood1", alpha = 0.5) +
  geom_line()  +
  theme_bw() +
  labs(x = NULL,
       y = "Density") +
  theme(text = element_text(size = 7),
        aspect.ratio = 0.8,
        legend.position = "none")

# Fig. 7a

ggsave("small_panel_rock_density.pdf",
       plot.rd,
       width = 45,
       height = 40,
       units = "mm")

# Annual sediment production

sed.prod <- bind_rows(
                tibble(sp = ann.sediment.prod.glac,
                       Type = "glacial"),
                tibble(sp = ann.sediment.prod.fluv,
                       Type = "fluvial"))

plot.sed.prod <- sed.prod %>%  
  ggplot(mapping = aes(x = sp, fill = Type, color = Type)) +
  geom_density( alpha = .5) +
  scale_fill_manual(values = c("lightblue", "darkorange")) +
  scale_color_manual(values = c("navy", "darkorange4")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  theme_bw() +
  labs(x = NULL,
       y = "Density") +
  theme(text = element_text(size = 7),
        aspect.ratio = 0.8,
        legend.position = "none")

# Fig. S7b

ggsave("small_panel_sediment_prod.pdf",
       plot.sed.prod,
       width = 45,
       height = 40,
       units = "mm")

# Annual sediment deposition

sed.depo <- bind_rows(
  tibble(sd = ann.sed.depo.glac,
         Type = "glacial"),
  tibble(sd = ann.sed.depo.fluv,
         Type = "fluvial"))

plot.sed.depo <- sed.depo %>%  
  ggplot(mapping = aes(x = sd, fill = Type, color = Type)) +
  geom_density( alpha = .5) +
  scale_fill_manual(values = c("lightblue", "darkorange")) +
  scale_color_manual(values = c("navy", "darkorange4")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  theme_bw() +
  labs(x = NULL,
       y = "Density") +
  theme(text = element_text(size = 7),
        aspect.ratio = 0.8,
        legend.position = "none")

# Fig. S7b

ggsave("small_panel_sediment_depo.pdf",
       plot.sed.depo,
       width = 45,
       height = 40,
       units = "mm")

# Bulk sediment density

x_vals <- seq(0.8, 2.4, length.out = 500)
df <- data.frame(
  x = x_vals,
  y = dnorm(x_vals, mean = 1.6, sd = 0.2)
)

plot.sd <- ggplot(df, aes(x = x, y = y)) +
  geom_area(fill = "grey75", alpha = 0.5) +
  geom_line()  +
  theme_bw() +
  labs(x = NULL,
       y = "Density")  +
  theme(text = element_text(size = 7),
        aspect.ratio = 0.8,
        legend.position = "none")

ggsave("small_panel_bulk_density.pdf",
       plot.sd,
       width = 45,
       height = 40,
       units = "mm")

# Trapping efficiency

x_vals <- seq(0, 1, length.out = 500)
df <- data.frame(
  x = x_vals,
  y = dbeta(x_vals, 9, 3)
)

plot.trapeff <- ggplot(df, aes(x = x, y = y)) +
  geom_area(fill = "steelblue", alpha = 0.5) +
  geom_line()  +
  theme_bw() +
  labs(x = NULL,
       y = "Density") +
  theme(text = element_text(size = 7),
        aspect.ratio = 0.8,
        legend.position = "none")

# Fig. S7c

ggsave("small_panel_trapping_efficiency.pdf",
       plot.trapeff,
       width = 45,
       height = 40,
       units = "mm")

# Lake volumes

lake.vol <- lakes.2020.no.sup %>%
  st_drop_geometry() %>%
  dplyr::select(starts_with("Pred_"))

med.vol <- apply(lake.vol, 1, median) 
med.vol <- tibble(medvol = med.vol)

plot.lakevol <- med.vol %>% 
  ggplot(mapping = aes(x = medvol*10^6)) +
  geom_density(fill = "blue") + 
  theme_bw()  +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  theme_bw() +
  labs(x = NULL,
       y = "Density") +
  theme(text = element_text(size = 7),
        aspect.ratio = 0.8,
        legend.position = "none")

# Fig. S7c

ggsave("small_panel_lake_volume.pdf",
       plot.lakevol,
       width = 45,
       height = 40,
       units = "mm")

###############################################################################

# Show how lake volume influences lake lifetime (i.e. theoretical year of 
# complete infill)

# From the estimated lifetimes per lake, just select 1000 samples
# of glacial and fluvial erosion rates to make it faster.

lifetime.by.col <- lapply(sed.prod.sys[seq(1, 5000, by = 5)], function (x) {
  
  y <-  x %>% dplyr::select("glac_er", "fluv_er")
  return(y)

  }) %>% bind_cols()

# combine both fluvial and glacial erosion rates for a given lake, and 
# estimate the median and 68% HDI of estimated lifetimes.

lifetime.by.lake <- apply(lifetime.by.col, 1, function (x) {
  hdi <- HDIofMCMC(x, 0.68)
  
  tibble(md = median(x),
         lo = hdi[1],
         hi = hdi[2])
})

lifetime.by.lake <- bind_rows(lifetime.by.lake)

# Analogous approach for lake volumes.

volume.by.lake <- lakes.2020.no.sup %>% 
  st_drop_geometry() %>% 
  dplyr::select(starts_with("Pred_")) %>%
  apply(., 1, function (x) {
    hdi <- HDIofMCMC(x, 0.68)
    
    tibble(md = median(x),
           lo = hdi[1],
           hi = hdi[2])
  })

volume.by.lake <- bind_rows(volume.by.lake)

# Given that resulting dataframes of estimated lifetimes and volumes have the 
# same structure, we can just join them side-by-side.

summary.sed.vol <- bind_cols(
  lifetime.by.lake %>% 
    dplyr::rename(lo_life = "lo",
                  md_life = "md",
                  hi_life = "hi"),
  volume.by.lake %>% 
    dplyr::rename(lo_vol = "lo",
                  md_vol = "md",
                  hi_vol = "hi"))

# Plot lake volume versus estimated lifetime until complete infill of all 
# non-supraglacial lakes as of 2020. We show the median in lake volume
# and lifetime as dots, and uncertainties as lines.

p <- ggplot(summary.sed.vol[sample(1:nrow(summary.sed.vol), 5000), ],
            mapping = aes(x = md_vol,
                          y = md_life)) +
  geom_linerange(aes(ymin = lo_life, 
                     ymax = hi_life),
                 alpha = 0.05, linewidth = 0.5) +
  geom_linerange(aes(xmin = lo_vol, 
                     xmax = hi_vol),
                 alpha = 0.05, linewidth = 0.5) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 3),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 3),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  theme_bw() +
  labs(x = expression("Lake volume [" * 10^6 ~ "m"^3 * "]"),
       y = "Estimated life time of glacier lakes") +
  theme(text = element_text(size = 7))

# Save this figure (Figure S10)

ggsave("volume_vs_lifetime.pdf",
       p,
       width = 170,
       height = 100,
       units = "mm")

################################################################################

# Show the land cover in the catchments feeding glacial lakes.

lakes.w.regions <- lakes.w.stats %>% 
  st_drop_geometry() %>%
  filter(Year_new == 2020,
         !is.na(glaciercov),
         !is.na(domin_catchm_lc)) %>%
  dplyr::left_join(., o1regions %>% st_drop_geometry(), 
                   by = "o1region", 
                   relationship = "many-to-one") 

regional.landcover <- lakes.w.regions %>% 
  group_by(full_name.y, domin_catchm_lc, .drop = F) %>%
  dplyr::summarise(n_type = n()) %>% 
  ungroup() %>%
  group_by(full_name.y) %>%
  dplyr::mutate(group_total = sum(n_type, na.rm = TRUE),
         prop = (n_type / group_total) *100) %>%
  ungroup() %>%
  dplyr::rename("Landcover" = "domin_catchm_lc")

glacier.coupled.region <- lakes.w.regions %>% 
  group_by(full_name.y) %>% 
  summarise(perc_glac = round((sum(glaciercoup == "Yes")/ n())*100) )

# Every landcover class has a different color.
# Those are the original colors in the WorldCover dataset, but the reviewer
# wanted to have higher contrast between "herbaceous wetland" and 
# "permanent water bodies" 


landcover.colors <- tibble(
  Landcover = c("Tree cover",
                "Shrubland",
                "Grassland",
                "Cropland",
                "Built-up",
                "Bare/sparse vegetation",
                "Snow and ice",
                "Permanent water bodies",
                "Herbaceous wetland",
                "Mangroves",
                "Moss and lichen"),
  hex = c("#006400",
          "#FFBB22",
          "#FFFF4C",
          "#F096FF",
          "#FA0000",
          "#B4B4B4",
          "#F0F0F0",
          "#0064C8",
          "#013c40",
          "#00CF75",
          "#FAE6A0"))

order.by.lat <- c("Arctic Canada North", 
                  "Arctic Canada South",
                  "Svalbard and Jan Mayen", 
                  "Russian Arctic",
                  "Greenland Periphery",
                  "Scandinavia",
                  "Iceland",
                  "Alaska",
                  "Western Canada and USA",
                  "Central Europe",
                  "Caucasus and Middle East",
                  "North Asia",
                  "Central Asia",
                  "South Asia West",   
                  "South Asia East",              
                  "Low Latitudes",           
                  "New Zealand",                                      
                  "Southern Andes")          
 
regional.landcover <- left_join(regional.landcover, 
                                landcover.colors, 
                                by = "Landcover")

regional.landcover$full_name.y <- factor(regional.landcover$full_name.y, 
                                         levels = rev(order.by.lat))

labels_df <- regional.landcover %>%
  filter(Landcover %in% c("Bare/sparse vegetation", "Snow and ice")) %>%
  mutate(label = paste0(round(prop), "%"))

# The plot shows the number of lakes in a given region and the regional share
# of lakes being in glacier contact.

sample_sizes <- regional.landcover %>%
  group_by(full_name.y) %>%
  summarise(group_total = unique(group_total)) %>%
  arrange(full_name.y) %>%
  mutate(group_total_comma = scales::comma(group_total)) %>%
  left_join(., glacier.coupled.region, by = "full_name.y") %>%
  mutate(samp_and_glac = paste0(group_total_comma, " (", perc_glac, "%)"))

# The plot is a bar plot with horizontal bars showing the study regions and
# the fraction of landcover in the catchments feeding lakes in a given region.

landcover.plot <- ggplot(regional.landcover, 
                         aes(x = full_name.y, 
                             y = prop , 
                             fill = Landcover)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = setNames(regional.landcover$hex, 
                                      regional.landcover$Landcover),
                    guide = guide_legend(nrow = 2)) +
  coord_flip(clip = "off") +
  geom_text(data = labels_df,
            aes(label = label),
            position = position_fill(vjust = 0.5),
            size = 6/(14/5)  , 
            color = "black") +
  geom_text(
    data = sample_sizes %>% distinct(full_name.y, samp_and_glac),
    aes(x = full_name.y, y = 1, label = samp_and_glac),  # place just outside bar
    inherit.aes = FALSE,
    hjust = 0,
    size = 7 / (14/5) 
  ) +
  scale_y_continuous(expand = c(0, 0),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Regions sorted by latitude [North to South]", 
       y = "Dominant land cover class in contributing catchments", 
       fill = "Land cover") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 7, color = "black"),
    plot.margin = unit(c(5, 55, 5, 5), "pt"),
    axis.text.y = element_text(colour = "black", size = 7),
    axis.text.x = element_text(colour = "black", size = 7)
  )

# Save plot showing the regional landcover in the contributing catchments of 
# glacial lakes (Figure 5).

ggsave("landcover.pdf",
       landcover.plot,
       width = 150,
       height = 120,
       units = "mm",
       device = cairo_pdf, 
       family = "Arial")
