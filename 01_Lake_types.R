### Preprocess lakes and assign dam types.

## Load packages that we'll need in the course of this script

# If packages haven't been installed beforehand, 
# first use install.packages("packagename")

require(tidyverse)
require(scales)
require(modelr)
require(doParallel)
require(sf)

# Set YOUR working directory. 

setwd("D:/nrc_user/veh/Zhang_glacial_lakes_global/")

# Download and unzip all glacial lake shapefiles from Zhang et al. (2024)
# (https://doi.org/10.11888/Cryos.tpdc.300938)into  your working directory.

# Download and unzip all glaciers mapped in the Randolph Glacier Inventory (RGI) 
# V7.0 (https://doi.org/10.5067/f6jmovy5navz)

################################################################################
## Read all lake shapefiles to disk ############################################

# Load outlines of RGI regions. We slightly modified these regions, e.g.
# making several subregions for Greenland, and adding more vertices along
# straight lines.

o2_new <- st_read("D:/nrc_user/veh/RGI2000-v7.0-regions/RGI2000-v7.0-o2regions_modified_densified.gpkg") %>%
  st_make_valid()

# List the locations of all lake shapefiles.

lake.list <- list.files(pattern = "Glacial_lake_2020.shp$|Glacial_lake_1990.shp$|Glacial_lake_2000.shp$",
                        recursive = T,
                        include.dirs = T)

# Remove Antarctica from the list.

lake.list <- grep("Antarctic", lake.list, invert = T, value = T)

# Register cluster here, export sf and tidyverse

cl <- makeCluster(30) # Adapt to the number of available cores on YOUR machine.

# Export R-packages to the cluster, which will be used in the apply loop. Make sure 
# you installed these packages before.

clusterEvalQ(cl = cl, c(library("sf"),
                        library("tidyverse")))

# Export the file list of the glacier shapefiles and the extent of the study 
# regions to the clusters.

clusterExport(cl = cl, list("lake.list", "o2_new"))

all.lakes <- pblapply(lake.list, cl = cl, function(x) {
  
  shp <- st_read(x) %>%
    st_make_valid() %>%
    st_transform(4326) %>%
    st_join(., o2_new, suffix = c("", ""))
  
  # Not all lakes in the first time step (~1990) were mapped in that year.
  # For simplicity, we assign all lakes that were mapped prior to 2010 the
  # 'New_year' the new value 1990.
  
  if(mean(shp$Year, na.rm = T) < 2010) {
    shp %>% mutate(Year_new = 1990)} else { shp %>% mutate(Year_new = 2020) }
  
})

stopCluster(cl)

# Make one large sf object from the regional time slices.

all.lakes <- do.call(rbind, all.lakes)

# Some geometries are empty. We remove those.

all.lakes <- all.lakes[!st_is_empty(all.lakes), , drop=FALSE]

# Correct RGI7.0 glacier geometries. Some are invalid. 

all.glaciers <- list.files(path = "D:/nrc_user/veh/RGI2000-v7.0-G-global",
                           pattern = ".shp$", 
                           recursive = T,
                           full.names = T, 
                           include.dirs = T)

for(i in 1: length(all.glaciers)) {
  
  st_read(all.glaciers[i]) %>%
    st_make_valid() %>%
    st_make_valid() %>%
    st_write(dsn = paste0(dirname(all.glaciers[i]), "/",
                          str_sub(basename(all.glaciers[i]), 1, -5),
                          "_valid.gpkg"), delete_dsn = T)
} 

# List corrected glacier inventories

all.glaciers <- list.files(path = "D:/nrc_user/veh/RGI2000-v7.0-G-global",
                           pattern =  "_valid.gpkg$", 
                           recursive = T,
                           full.names = T, 
                           include.dirs = T)

# Iterate over each region and check if polygons are entirely within a glacier. 
# These lakes are deemed to be supraglacial lakes.

uni.reg <- unique(all.lakes$o1region)

all.lakes.list <- list()

for(i in uni.reg) {
  
  sf_use_s2(FALSE)
  
  # Filter only lakes from a given O1 region.
  
  lakes.reg <- all.lakes %>% 
    filter(o1region == i) %>%
    st_make_valid()
  
  # The same filter operation for all glaciers.
  
  glaciers.reg <- all.glaciers[str_sub(basename(all.glaciers), 16, 17) %in% i] %>%
    st_read() %>%
    st_make_valid()
  
  # Lakes that are fully within the RGI outlines are deemed to be supraglacial.
  
  lakes.in.glaciers <- st_within(lakes.reg, glaciers.reg, sparse = F)
  
  supras.idx <- apply(lakes.in.glaciers, 1, any)
  
  # We add a new empty column called 'Lake_type'.
  # Lakes that are Type1 == 1 become moraine or bedrock, the Type1 == 3 are
  # ice-dammed lakes according to Zhang et al., and those that are fully 
  # within the glacier are supraglacial lakes, according to our definition.
  
  lakes.reg$Lake_type <- NA
  lakes.reg$Lake_type[lakes.reg$Type1 == 1] <- "moraine/bedrock"
  lakes.reg$Lake_type[lakes.reg$Type1 == 3] <- "ice"
  lakes.reg$Lake_type[supras.idx] <- "supraglacial"
  
  all.lakes.list[[i]] <- lakes.reg
  
}

# Now, all lakes have dam types

all.lakes.types <- do.call(rbind, all.lakes.list)

# A few were actually non-glacier fed lakes (Type1 = 2), which we remove from this
# shapefile. Then, we add a unique ID such that we can identify each lake 
# indivudally.

all.lakes.types <- all.lakes.types %>%
  filter(Type1 != 2) %>%
  mutate(alog = log10(Area), 
         UniqueID = 1: n())


saveRDS(all.lakes.types, "all_lakes_including_dam_type.rds")
# all.lakes.types <- readRDS("all_lakes_including_dam_type.rds")

################################################################################

# Select only lakes from 2020

lakes.2020 <- all.lakes.types %>%
  filter(Year_new == 2020)

# Export centroid coordinates to calculate flow paths to the ocean.

# lakes.2020 %>%
#   st_point_on_surface() %>%
#   st_write(., dsn = "V:/xchange/wolfgang/glacier_lakes_2020_centroids.shp")

lakes.1990 <- all.lakes.types %>%
  filter(Year_new == 1990) 

# Given that we do not know the lake type for those lakes being not supraglacial
# or ice-dammed (it could be either moraine, bedrock, or moraine/bedrock), 
# we randomly assign any of these types. We could also play more with that later.

lakes.2020$Lake_type[!((lakes.2020$Lake_type == "ice") | (lakes.2020$Lake_type == "supraglacial"))] <- 
  sample(c("moraine", "bedrock", "moraine/bedrock"),
         size = length(which(!((lakes.2020$Lake_type == "ice") | (lakes.2020$Lake_type == "supraglacial")))), 
         replace = T)

int.2020.1990 <- st_intersects(lakes.2020, lakes.1990, sparse = F)

# Assign the same dam type in 1990 as of 2020.

for(i in 1:ncol(int.2020.1990)) {
  
  # If the lake was ice or supraglacial in 1990 according to Zhang et al. (2024), 
  # then we keep that dam type.
  
  if((lakes.1990$Lake_type[i] == "ice") | (lakes.1990$Lake_type[i] == "supraglacial")) {next}
  
  int.idx <- which(int.2020.1990[ ,i])
  
  # If the lake did exist in 1990, but not in 2020, we randomly sample a new dam 
  # type.
  
  if(length(int.idx) == 0) {
    
    lakes.1990$Lake_type[i] <- sample(c("moraine", 
                                        "bedrock", 
                                        "moraine/bedrock"), 
                                      1) } else if (length(int.idx) == 1) {
                                        
                                        lakes.1990$Lake_type[i] <- lakes.2020$Lake_type[int.idx]
                                        
                                      } else if (length(int.idx) > 1) {
                                        
                                        lakes.1990$Lake_type[i] <-  sample(lakes.2020$Lake_type[int.idx], 1)
                                        
                                      }
  
}

# The intersection object is quite large, but is not necessary anymore.
# rm(int.2020.1990)

# Combine lakes from both years including their dam types.

all.lakes <- rbind(lakes.1990, lakes.2020)

# Save this object to disk. This will be our 'master object' for further analysis.

saveRDS(all.lakes, "Lakes19902020_damtype.RDS")
