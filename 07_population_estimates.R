require(pbapply)
require(parallel)
require(sf)
require(tidyverse)
require(terra)
require(exactextractr)
require(cowplot)

setwd("D:/nrc_user/veh/Zhang_glacial_lakes_global")

catch.dirs <- list.dirs(path = "D:/nrc_user/veh/All_lakes", full.names = T, recursive = F)

cl <- makeCluster(28) # Adapt to the number of available cores on YOUR machine.

# Export R-packages to the cluster, which will be used in the apply loop. Make sure
# you installed these packages before.

clusterEvalQ(cl = cl, c(library("sf"),
                        library("tidyverse"),
                        library("stringr")))

# Export the file list of the glacier shapefiles and the extent of the study
# regions to the clusters.

clusterExport(cl = cl, list("catch.dirs"))

catch.stats <- pblapply(catch.dirs, cl = cl, function (x) {
  
  catch.files <- list.files(pattern = "_catchment.gpkg",
                            path = x,
                            full.names = T, 
                            include.dirs = T,
                            recursive = T)
  
  stats.in.dir <- lapply(catch.files, function (y) {
    
    dir.str <- strsplit(dirname(y), "/")[[1]]
    uni.id <- dir.str[length(dir.str)]
    o1region <- str_sub(dir.str[length(dir.str)-1], 1,2)
    
    result <- tryCatch({
      st_read(y, quiet = T) %>%
       mutate(UniqueID = uni.id,
              o1region = o1region) %>%
        st_transform(4326)
    }, error = function(e) {
      NULL
    })
    
  } )
  
  stats.in.dir <- stats.in.dir[!sapply(stats.in.dir, is.null)]
  
  catch.bind <- do.call(rbind, stats.in.dir)
  
  return(catch.bind)
  
})

stopCluster(cl)

# Bind all individual catchments

catchm.all <- do.call(rbind, catch.stats)

st_write(catchm.all, "all_catchments.gpkg")
saveRDS(catchm.all, "all_catchments.rds")

# Dissolving and deleting holes done in QGIS

catchm.reg <- st_read(dsn = "all_catchments_dissolved_no_holes.gpkg") %>%
  st_make_valid()

# Flow paths from MERIT DEM

all.lakes <- readRDS("all_lakes_with_volumes.RDS")

idx <- seq(1, (nrow(all.lakes %>% filter(Year_new == 2020))-1) *5, by = 5)

# Read flowpath data extracted from the 90 m Merit DEM.
# Data are a irregular txt file organised by lines. Every lake is identified 
# via its UniqueID, and the X, Y, distance, and Z coordinate from the source towards
# the coast. So, after 5 lines, a new flow path starts.

fread <- read_lines("D:/nrc_user/veh/Zhang_glacial_lakes_global/flowpaths_detailed.txt")

registerDoParallel(12)

bufs <- foreach(i = idx,
       .packages = c("dplyr",
                     "sf",
                     "stringr")
       ) %dopar% {

  
  # We only extract the last entry, which is the entails the cumulative flow path
  # length from the source.
  
  fp <- tibble(UniqueID = as.numeric(fread[i]),
               dist = fread[i+1] %>% str_split(" ") %>% unlist() %>% as.numeric() %>% na.omit(),
               X = fread[i+2] %>% str_split(" ") %>% unlist() %>% as.numeric() %>% na.omit(),
               Y = fread[i+3] %>% str_split(" ") %>% unlist() %>% as.numeric() %>% na.omit(),
               Z = fread[i+4] %>% str_split(" ") %>% unlist() %>% as.numeric() %>% na.omit()) %>%
    filter(dist <= 50)
  
  if (nrow(fp) <= 2 ) {return(NULL)}
  
  line <- fp %>% 
    group_by(UniqueID) %>% 
    dplyr::summarise(
      geometry = st_sfc(st_linestring(as.matrix(dplyr::select(cur_data(), Y, X)))), 
      .groups = "drop"
    ) %>%
    st_as_sf(crs = 4326) %>%
    st_make_valid()
   
  buf <- st_buffer(line, dist = 1000, endCapStyle = "ROUND") %>%
    st_make_valid()
  
  return(buf)
    
}

# Bind all flow paths to one object.

flowpath.bufs <- bufs[!sapply(bufs, is.null)] 

all.lakes.region <- all.lakes %>% 
  filter(Year_new == 2020) %>% 
  st_drop_geometry() %>% 
  as_tibble() %>%
  dplyr::select(c("UniqueID", "o1region", "o2region", "full_name", "long_code"))

flowpath.w.region <- lapply(flowpath.bufs, function (x) {left_join(x, all.lakes.region, by = "UniqueID")})

uni.reg <- unique(all.lakes.region$o1region)

registerDoParallel(18)

bufs <- foreach(i = uni.reg,
                .packages = c("dplyr",
                              "sf") ) %dopar% {
  
  
  # We only extract the last entry, which is the entails the cumulative flow path
  # length from the source.
  
  idx <- sapply(flowpath.w.region, function(x) {x$o1region == i} )                              
  
  reg.bufs <- lapply(flowpath.w.region[idx], st_make_valid)
  
  regional.buffers <- do.call(rbind, reg.bufs) 
  return(regional.buffers)
  
}

flowpath.bufs <- do.call(rbind, bufs)

st_write(flowpath.bufs, "flowpath_buffers.gpkg")

# Dissolving in QGIS, calculating difference between catchments and dissolved
# flowpath buffers.

flowpath.dissolved <- st_read("diff.gpkg") %>%
  st_make_valid() 

# Read four different gridded population data sets, Landscan, WorldPop, GPW, and
# GHS. For each dataset, first extract the population upstream, then downstream
# of glacial lakes.

# Landscan2020 data (year when lakes were mapped): https://landscan.ornl.gov/

ls <- rast("D:/nrc_user/veh/Population_data/landscan-global-2020-assets/")

ls.catch.extr <- exact_extract(ls, catchm.reg, include_cols = "o1region")
ls.catch.all <- do.call(rbind, ls.catch.extr) %>%
  dplyr::mutate(pop_weight = value * coverage_fraction) %>%
  group_by(o1region) %>%
  dplyr::summarise(pop_sum_ls = sum(pop_weight, na.rm = T))

ls.ds.extr <- exact_extract(ls, flowpath.dissolved, include_cols = "o1region")
ls.ds.all <- do.call(rbind, ls.ds.extr) %>%
  dplyr::mutate(pop_weight = value * coverage_fraction) %>%
  group_by(o1region) %>%
  dplyr::summarise(pop_sum_ls = sum(pop_weight, na.rm = T))

# Worldpop 1km global: https://hub.worldpop.org/geodata/listing?id=64

wpop <- rast("D:/nrc_user/veh/Population_data/ppp_2020_1km_Aggregated.tif")

wpop.catch.extr <- exact_extract(wpop, catchm.reg, include_cols = "o1region")
wpop.catch.all <- do.call(rbind, wpop.catch.extr) %>%
  dplyr::mutate(pop_weight = value * coverage_fraction) %>%
  group_by(o1region) %>%
  dplyr::summarise(pop_sum_wpop = sum(pop_weight, na.rm = T))

wpop.ds.extr <- exact_extract(wpop, flowpath.dissolved, include_cols = "o1region")
wpop.ds.all <- do.call(rbind, wpop.ds.extr) %>%
  dplyr::mutate(pop_weight = value * coverage_fraction) %>%
  group_by(o1region) %>%
  dplyr::summarise(pop_sum_wpop = sum(pop_weight, na.rm = T))

# NASA's Gridded Population of the World (GPW), version 4
# https://www.earthdata.nasa.gov/data/projects/gpw/data-access-tools

gpw <- rast("D:/nrc_user/veh/Population_data/gpw-v4-population-count-rev11_2020_30_sec_tif/gpw_v4_population_count_rev11_2020_30_sec.tif")

gpw.catch.extr <- exact_extract(gpw, catchm.reg, include_cols = "o1region")
gpw.catch.all <- do.call(rbind, gpw.catch.extr) %>%
  dplyr::mutate(pop_weight = value * coverage_fraction) %>%
  group_by(o1region) %>%
  dplyr::summarise(pop_sum_gpw = sum(pop_weight, na.rm = T))

gpw.ds.extr <- exact_extract(gpw, flowpath.dissolved, include_cols = "o1region")
gpw.ds.all <- do.call(rbind, gpw.ds.extr) %>%
  dplyr::mutate(pop_weight = value * coverage_fraction) %>%
  group_by(o1region) %>%
  dplyr::summarise(pop_sum_gpw = sum(pop_weight, na.rm = T))

# Global Human Settlement Layer (GHSL)
# https://human-settlement.emergency.copernicus.eu/download.php

ghs <- rast("D:/nrc_user/veh/Population_data/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0.tif")

ghs.catch.extr <- exact_extract(ghs, catchm.reg, include_cols = "o1region")
ghs.catch.all <- do.call(rbind, ghs.catch.extr) %>%
  dplyr::mutate(pop_weight = value * coverage_fraction) %>%
  group_by(o1region) %>%
  dplyr::summarise(pop_sum_ghs = sum(pop_weight, na.rm = T))

ghs.ds.extr <- exact_extract(ghs, flowpath.dissolved, include_cols = "o1region")
ghs.ds.all <- do.call(rbind, ghs.ds.extr) %>%
  dplyr::mutate(pop_weight = value * coverage_fraction) %>%
  group_by(o1region) %>%
  dplyr::summarise(pop_sum_ghs = sum(pop_weight, na.rm = T))

gc()
closeAllConnections()

# Bind all regional estimates of population counts

reg.pop.counts.catch <- left_join(ls.catch.all, wpop.catch.all, by = "o1region") %>%
  left_join(., gpw.catch.all, by = "o1region") %>%
  left_join(., ghs.catch.all, by = "o1region") %>%
  mutate(location = "Upstream")

all.catch <- reg.pop.counts.catch %>%
  dplyr::summarise(across(starts_with("pop_"), sum, na.rm = TRUE)) %>%
  mutate(location = "Upstream")

# HMA alone

hma.catch <- reg.pop.counts.catch %>%
  filter(o1region == "13" | o1region == "14" | o1region == "15") %>%
  dplyr::summarise(across(starts_with("pop_"), sum, na.rm = TRUE))

# Arctic alone

arctic.catch <- reg.pop.counts.catch %>%
  filter(o1region == "03" | o1region == "04" | o1region == "05" | o1region == "09" ) %>%
  dplyr::summarise(across(starts_with("pop_"), sum, na.rm = TRUE))

arctic.catch/ all.catch %>% dplyr::select(!location) * 100

reg.pop.counts.ds <- left_join(ls.ds.all, wpop.ds.all, by = "o1region") %>%
  left_join(., gpw.ds.all, by = "o1region") %>%
  left_join(., ghs.ds.all, by = "o1region") %>%
  mutate(location = "Downstream")

all.ds <- reg.pop.counts.ds %>%
  dplyr::summarise(across(starts_with("pop_"), sum, na.rm = TRUE)) %>%
  mutate(location = "Downstream")

# HMA alone

hma.ds <- reg.pop.counts.ds %>%
  filter(o1region == "13" | o1region == "14" | o1region == "15") %>%
  dplyr::summarise(across(starts_with("pop_"), sum, na.rm = TRUE))

# Arctic alone

arctic.ds <- reg.pop.counts.ds %>%
  filter(o1region == "03" | o1region == "04" | o1region == "05" | o1region == "09" ) %>%
  dplyr::summarise(across(starts_with("pop_"), sum, na.rm = TRUE))

hma.ds

hma.ds/all.ds
arctic.ds/all.ds*100

require(scales)

# O1 regions from the RGI.

o1regions <- st_read("D:/nrc_user/veh/RGI2000-v7.0-regions/RGI2000-v7.0-o1regions.shp")

# Combine both population counts upstream and downstream of lakes. We take
# the median and the min and max of these estimates per region.

regional.pop <- rbind(reg.pop.counts.catch, reg.pop.counts.ds) %>% 
  rowwise() %>%
  mutate(
    pop_min = min(c_across(starts_with("pop_")), na.rm = TRUE),
    pop_min = if_else(pop_min < 1, 1, pop_min),
    pop_median = median(c_across(starts_with("pop_")), na.rm = TRUE),
    pop_median = if_else(pop_median < 1, 1, pop_median),
    pop_max = max(c_across(starts_with("pop_")), na.rm = TRUE),
    pop_max = if_else(pop_max < 1, 1, pop_max),
    location = factor(location, levels = c("Upstream", "Downstream"))
  ) %>% 
  left_join(., o1regions, 
            by = "o1region",
            relationship = "many-to-many") %>%
  mutate(o1region2 = paste0(o1region, ": ", full_name))

# Do the same also for the global population

global.pop <- rbind(all.ds, all.catch) %>% 
  rowwise() %>%
  mutate(
    pop_min = min(c_across(starts_with("pop_")), na.rm = TRUE),
    pop_min = if_else(pop_min < 1, 1, pop_min),
    pop_median = median(c_across(starts_with("pop_")), na.rm = TRUE),
    pop_median = if_else(pop_median < 1, 1, pop_median),
    pop_max = max(c_across(starts_with("pop_")), na.rm = TRUE),
    pop_max = if_else(pop_max < 1, 1, pop_max),
    location = factor(location, levels = c("Upstream", "Downstream")),
    o1region2 = "Global"
  ) 

# Produce bar plots of regional population counts.

pop.catch.ds.plot <- rbind(regional.pop %>% 
                             dplyr::select(location, pop_min, pop_median, pop_max, o1region2), 
      global.pop %>% 
        dplyr::select(location, pop_min, pop_median, pop_max, o1region2)) %>%
  ggplot(mapping = aes(x = location,
                       y = pop_median,
                       fill = location)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.9)) +
  geom_linerange(mapping = aes(ymin = pop_min,
                               ymax = pop_max)) +
  facet_wrap(~o1region2, scales = "free_y") + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 3),
                labels = trans_format("log10", math_format(10 ^ .x)),
                limits = c(1, 10^7))  +
  scale_fill_manual(values = c("lightgreen", "darkgreen")) +
  labs(x = "",
       y = "") +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    legend.position = "none",
    legend.key = element_blank(),
    aspect.ratio = 1.5,
    plot.background = element_blank(),
    strip.background = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(2, "lines")
  ) + coord_cartesian(clip = "off")



ggsave(filename = "population_upstream_downstream.pdf",
       plot = pop.catch.ds.plot, 
       width = 180,
       height = 120,
       units = "mm",
       device = cairo_pdf, 
       family = "Arial")

################################################################################

sf_use_s2(FALSE)

require(sf)
require(lwgeom)

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


world_map <- 
  ggplot() + 
  geom_sf(data = robin_outline, fill = ocean_color, color = NA, inherit.aes = F) +
  geom_sf(data = world_robin, color = NA, fill = land_color, inherit.aes = F) +
  geom_sf(data = o1_robin  %>% filter(o1region != "19", 
                                      o1region != "20"), 
          color = "white",
          fill = fill_alpha("lightblue", 0.5), inherit.aes = F) +
  
  theme_map() +
  theme( legend.position = "bottom",  # Center at the bottom
         legend.justification = c(0.5, 0),  # Align the legend's center
         legend.title=element_text(size = 20), 
         legend.text=element_text(size = 20),
         legend.box = "vertical") 

ggsave(plot     = world_map, 
       filename = "World_map.pdf",
       width    = 120,
       height   = 90,
       units    = "mm")


################################################################################

volume.regional.lakes <- readRDS("volume_regional_lakes.RDS") %>%
  filter(Year_new == 2020)

pop.vol <- left_join(regional.pop, volume.regional.lakes, by = "o1region")

plot.vol.vs.pop <- ggplot(pop.vol,
       aes(x = md/1000,
       y = pop_median,
       fill = location,
       color = location)) +
  geom_point(color = "white", shape = 21, stroke = 0.3) +
  geom_linerange(mapping = aes(xmin = lo/1000,
                               xmax = hi/1000)) +
  geom_linerange(mapping = aes(ymin = pop_min,
                               ymax = pop_max)) +
  geom_text_repel(aes(label = o1region),
                  color = "black",
                  size = 7/.pt,
                  max.overlaps = 25) +
  scale_fill_manual(values = c("lightgreen", "darkgreen")) +
  scale_color_manual(values = c("lightgreen", "darkgreen")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 5),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ x, n = 3),
                labels = trans_format("log10", math_format(10 ^ .x)))  +
  labs(x = "Regional lake volume [km³] in 2020",
       y = "Regional population count in ~2020") +
  theme_bw() +
  facet_wrap(~location) +
  theme(text = element_text(size = 7),
        aspect.ratio = 0.5,
        legend.position = "none",
        panel.background = element_blank(),
        strip.background = element_blank())

ggsave("Regional_volume_vs_population.pdf",
       plot.vol.vs.pop,
       width = 140,
       height = 80,
       units = "mm",
       device = cairo_pdf, 
       family = "Arial")
