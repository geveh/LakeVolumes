### Download the Copernicus GLO30 DEM and land cover maps. 
### Derive the catchment of each lake, and summarise catchment-wide statistics 
### including glacial cover, land cover, and relief.

### The script makes use of SAGA GIS, called through the package RGSAGA.
### Thus, you need to pre-install SAGA GIS (here V9.6.1) before executing this code.

## Load packages that we'll need in the course of this script

# If packages haven't been installed beforehand, 
# first use install.packages("packagename")

require(tidyverse)
require(sf)
require(rstac)
require(stars)
require(pbapply)
require(doParallel)
require(dplyr)
require(RSAGA)
require(units)
require(terra)
require(gdalUtilities)

# Define useful functions

# Get the corresponding UTM zone from a given longitude

get_utm_zone <- function(lon) {
  return((floor((lon + 180) / 6) %% 60) + 1)
}

# Get the EPSG code of a UTM zone from a given latitude and longitude

get_epsg_code <- function(lat, lon) {
  utm_zone <- get_utm_zone(lon)
  if (lat >= 0) {
    epsg_code <- 32600 + utm_zone
  } else {
    epsg_code <- 32700 + utm_zone
  }
  return(epsg_code)
}

# We will download all COP-30 DEM tiles for the regions in the Randolph Glacier
# Inventory, and deposit the data in the same folder structure

o2 <- st_read("D:/nrc_user/veh/RGI2000-v7.0-regions/RGI2000-v7.0-o2regions.shp") %>%
  st_make_valid()

o2_new <- st_read("D:/nrc_user/veh/RGI2000-v7.0-regions/RGI2000-v7.0-o2regions_modified_densified.gpkg") %>%
  st_make_valid()

# Establish a stac endpoint to Microsoft Planetary Computer.

s <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")

# We iterate over the RGI regions and make a stac query for 
# Copernicus-30m DEM data covering the extent of the RGI region.
# We then sign in to planetary computer, and download the data to a new directory.

for (i in 1:nrow(o2)) {

  bbox.grid <- st_bbox(o2[i, ])

  stac_query <- s %>%
    stac_search(collections = "cop-dem-glo-30",
                bbox = c(bbox.grid$xmin,
                         bbox.grid$ymin,
                         bbox.grid$xmax,
                         bbox.grid$ymax),
                limit = 999)


  executed_stac_query <- rstac::get_request(stac_query)

  signed_stac_query <- rstac::items_sign(
    executed_stac_query,
    rstac::sign_planetary_computer()
  )

  target.dir <- paste0("D:/nrc_user/veh/COPs/", o2$o2region[i])

  dir.create(target.dir)

  rstac::assets_download(signed_stac_query, "data",
                         output_dir = target.dir,
                         overwrite = T)

}

# Some DEMs overlap with several RGI regions, and we delete the duplicate ones.

all.dem.paths <- list.files(path = "D:/nrc_user/veh/COPs/",
                            pattern = "DEM.tif$",
                            recursive = T,
                            full.names = T,
                            include.dirs = T)

dems.to.delete <- all.dem.paths %>%
  basename() %>%
  duplicated()

file.remove(all.dem.paths[dems.to.delete])

# Some DEMs seem to be empty (only a few kilobytes large). We delete
# them, and instead download the 90 m version of the DEM.

file.remove(all.dem.paths[file.size(all.dem.paths) < 1000])

# Generate a lookup shapefile that contains the extent of all DEM tiles

remaining.dems <- list.files(path = "D:/nrc_user/veh/COPs/",
                             pattern = "DEM.tif$",
                             recursive = T,
                             full.names = T,
                             include.dirs = T)

# We manually identify regions where we need to download additional data (this
# largely occurs in Greenland)

s <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")

for (i in c(39,41,47,48)) {

  bbox.grid <- st_bbox(o2_new[i, ])

  stac_query <- s %>%
    stac_search(collections = "cop-dem-glo-90",
                bbox = c(bbox.grid$xmin,
                         bbox.grid$ymin,
                         bbox.grid$xmax,
                         bbox.grid$ymax),
                limit = 999)

  executed_stac_query <- rstac::get_request(stac_query)

  available.dems <- sapply(executed_stac_query$features, function (x) basename(x$assets$data$href))

  tf <- !(str_sub(available.dems, -22, -1) %in% str_sub(basename(remaining.dems), -22, -1))

  if (all(tf == F)) next

  filtered_stac_query <- executed_stac_query$features[tf]

  # Rebuild the stac_response object
  filtered_stac_response <- executed_stac_query
  filtered_stac_response$features <- filtered_stac_query

  # Confirm the class is preserved
  class(filtered_stac_response)

  signed_stac_query <- rstac::items_sign(
    filtered_stac_response,
    rstac::sign_planetary_computer()
  )

  target.dir <- paste0("D:/nrc_user/veh/COPs/", o2_new$o2region[i])

  if (!dir.exists(target.dir)) {

    dir.create(target.dir)

  }

  rstac::assets_download(signed_stac_query, "data",
                         output_dir = target.dir,
                         overwrite = T)

}

# Delete duplicate COP DEMs after download of the additional tiles.

all.dem.paths <- list.files(path = "D:/nrc_user/veh/COPs/",
                            pattern = "DEM.tif$",
                            recursive = T,
                            full.names = T,
                            include.dirs = T)

dems.to.delete <- all.dem.paths %>%
  basename() %>%
  str_sub(., -22, -1) %>%
  duplicated()

file.remove(all.dem.paths[dems.to.delete])

# We obtain the extent of each DEM and generate a simple feature of it.

cl <- makeCluster(35) 

clusterEvalQ(cl = cl, c(library("sf"),
                        library("tidyverse"),
                        library("stars")))

all.dem.paths <- list.files(path = "D:/nrc_user/veh/COPs",
                            pattern = "DEM.tif$",
                            recursive = T,
                            full.names = T,
                            include.dirs = T)

dem.extents <- pblapply(all.dem.paths, cl = cl, function (x) {
  
  f <- read_stars(x) %>% 
    st_bbox() %>% 
    st_as_sfc() %>% 
    st_sf()
  
  f$File <- x
  
  f
  
})

stopCluster(cl)

dem.ext <- do.call(rbind, dem.extents)


# Download landcover maps. We again use Microsoft's Planetary Computer API to
# download the ESA Worldcover Version 2 data.
# The data are again deposited in the structure of the O2 region of the RGI.

s <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/",
          force_version = 1)

for (i in 1:nrow(o2)) {

  bbox.grid <- st_bbox(o2[i, ])

  stac_query <- s %>%
    stac_search(collections = "esa-worldcover",
                bbox = c(bbox.grid$xmin,
                         bbox.grid$ymin,
                         bbox.grid$xmax,
                         bbox.grid$ymax),
                limit = 999)
  
  executed_stac_query <- rstac::get_request(stac_query) %>% 
    items_filter(properties$`esa_worldcover:product_version` == "2.0.0")
  
  signed_stac_query <- rstac::items_sign(
    executed_stac_query,
    rstac::sign_planetary_computer()
  )
  
  target.dir <- paste0("D:/nrc_user/veh/Landcover/", o2$o2region[i])

  dir.create(target.dir)

  rstac::assets_download(signed_stac_query, "map",
                         output_dir = target.dir,
                         overwrite = T)
  
}

# Delete duplicate landcover maps.

all.lc.paths <- list.files(path = "D:/nrc_user/veh/Landcover/",
                            pattern = "_Map.tif$",
                            recursive = T,
                            full.names = T,
                            include.dirs = T)

lc.to.delete <- all.lc.paths %>%
  basename() %>%
  duplicated()

file.remove(all.lc.paths[lc.to.delete])

# Generate a lookup shapefile that contains the extent of all land cover tiles.

remaining.lcs <- list.files(path = "D:/nrc_user/veh/Landcover/",
                            pattern = "Map.tif$",
                            recursive = T,
                            full.names = T,
                            include.dirs = T)


cl <- makeCluster(35) 

clusterEvalQ(cl = cl, c(library("sf"),
                        library("tidyverse"),
                        library("stars")))


lc.extents <- pblapply(remaining.lcs, cl = cl, function (x) {

  f <- read_stars(x) %>%
    st_bbox() %>%
    st_as_sfc() %>%
    st_sf()

  f$File <- x

  f

})

stopCluster(cl)

# Bind all extent polygons to one large polygon.

lc.ext <- do.call(rbind, lc.extents)

################################################################################

# Catchment-wide analysis.

# Load the combined shapefile of all lakes in 1990 and 2020, including their 
# predicted volumes.

all.lakes <- readRDS("D:/nrc_user/veh/Zhang_glacial_lakes_global/Lakes19902020_damtype.RDS")

# Filter only lakes as of 2020.

lakes.2020 <- all.lakes %>%
  filter(Year_new == 2020)

# We will derive catchments, assess relief, glacier and cover per glacial lake.
# We will deposit these data for each lake individually, structured by RGI region.

main.dir <- "D:/nrc_user/veh/All_lakes/" 

dir.create(main.dir)

# Correct invalid glacier geometries in the Randolph Glacier Inventory

all.glaciers <- list.files(path = "D:/nrc_user/veh/RGI2000-v7.0-G-global",
                           pattern = ".shp$", 
                           recursive = T,
                           full.names = T, 
                           include.dirs = T)

for(i in 1: length(all.glaciers)) {
  
  st_read(all.glaciers[i]) %>%
    st_make_valid() %>%
    st_write(dsn = paste0(dirname(all.glaciers[i]), "/",
                          str_sub(basename(all.glaciers[i]), 1, -5),
                          "_valid.gpkg"), delete_dsn = T)
} 

# Catchments are derived using SAGA GIS. Define SAGA environment, a pointer
# for RSAGA towards the installed SAGA version.

env <- rsaga.env()

registerDoParallel(35)

foreach(i = 1:nrow(lakes.2020),
       .packages = c("tidyverse",
                     "sf",
                     "RSAGA",
                     "gdalUtilities",
                     "stars",
                     "terra",
                     "units")) %dopar% {                           
  
  spec.lake <- lakes.2020[i, ] %>%
    st_make_valid()
  
  # Get the EPSG code of lake
  
  epsg_code <- get_epsg_code(spec.lake$Lat, spec.lake$Lon)
  
  # Transform lake to UTM
  
  spec.lake.utm <- spec.lake %>%
    st_transform(epsg_code)
  
  # Write shapefile to disk in the corresponding RGI region.
  
  region.dir <- paste0(main.dir, spec.lake.utm$o2region)
  
  if(!dir.exists(region.dir)) {dir.create(region.dir)}
  
  lake.dir <- paste0(region.dir, "/", spec.lake.utm$UniqueID)
  
  # Create a separate folder for each lake.
  
  dir.create(lake.dir, recursive = T)
  
  shape.loc <- paste0(lake.dir, "/", spec.lake.utm$UniqueID, ".shp")
  
  st_write(spec.lake.utm, dsn = shape.loc, delete_dsn = T)
  
  # Find overlapping COP30-DEMs with the lake.
  
  lake.dem.int <- st_intersects(spec.lake, dem.ext, sparse = F)
  
  suit.dem <- dem.ext[lake.dem.int[1, ], ]
  
  catchment.int <- nrow(suit.dem) + 1
  int.catch.dem <- lake.dem.int[1, ]
  
  # Within the lake directory, we generate a new DEM directory, in which
  # we will process the catchments
  
  dem.dir <- paste0(lake.dir, "/DEM/")
  dir.create(dem.dir)
  
  # If several DEMs cover a lake, they need to be stitched to one larger DEM,
  # in SAGA GIS format (.sgrd/ .sdat).
  # We then fill each DEM, derive the upstream contributing area, both in UTM,
  # and in WGS84.
  # The while-loop makes sure that the process of finding new adjacent DEM tiles
  # is repeated, if the upstream contributing catchment hits the boundary of 
  # the DEM tile.
  
  while(nrow(suit.dem) < catchment.int) {
    
    suit.dem <- dem.ext[int.catch.dem, ]
    
    raw.dem.sdat <- paste0(dem.dir, "mosaic.sdat")
    raw.dem.sgrd <- paste0(dem.dir, "mosaic.sgrd")
    
    if (nrow(suit.dem) > 1) {
      
      lake.dem <- stars::st_mosaic(suit.dem$File,
                            dst = paste0(dem.dir, "temp"),
                            options = c("-resolution", "highest",
                                        "-vrtnodata", "-99999", 
                                        "-srcnodata", "-32767",
                                        "-r", "cubic")) 
      
      gdalUtilities::gdalwarp(srcfile = paste0(dem.dir, "temp"),
               dstfile =  raw.dem.sdat,
               t_srs = paste0("epsg:",epsg_code), 
               tr = c(25,25), 
               r = "cubic",
               dstnodata = -99999,
               wm = 500,
               overwrite = T)
      
    } else { lake.dem <- suit.dem$File
    
    gdalUtilities::gdalwarp(srcfile = lake.dem,
                             dstfile = raw.dem.sdat,
                             t_srs = paste0("epsg:",epsg_code), 
                             tr = c(25,25), 
                             r = "cubic",
                             dstnodata = -99999,
                             wm = 500,
                             overwrite = T)
    
    }
    
    # Fill DEM
    
    filled.dem.sgrd <- paste0(dem.dir, "filled.sgrd")
    filled.dem.sdat <- paste0(dem.dir, "filled.sdat")
    
    rsaga.geoprocessor(lib = "ta_preprocessor", 
                       module = 5, 
                       env = env,
                       param = list(ELEV =   raw.dem.sgrd,
                                    FILLED = filled.dem.sgrd,
                                    MINSLOPE = 0.05))
    
    # Generate Target Grid
    
    lake.rasterized <- paste0(dem.dir, "lake_rasterized.sgrd")
    
    rsaga.geoprocessor(lib = "grid_gridding", 
                       env = env,
                       module = 0, 
                       param = list(INPUT = shape.loc,
                                    OUTPUT = 0,
                                    TARGET_DEFINITION = 1,
                                    TARGET_TEMPLATE = filled.dem.sgrd,
                                    GRID = lake.rasterized))
    
    # Calculate upslope contributing area
    
    upslope.area.sgrd <- paste0(dem.dir, "upslope_area.sgrd")
    upslope.area.sdat <- paste0(dem.dir, "upslope_area.sdat")
    
    rsaga.geoprocessor(lib = "ta_hydrology", 
                       module = 4, 
                       env = env,
                       param = list(TARGET = lake.rasterized,
                                    ELEVATION = filled.dem.sgrd,
                                    AREA = upslope.area.sgrd,
                                    METHOD = 0))
    
    # Convert the SAGA grid showing the upslope area into a polygon
    
    catchment.utm <- read_stars(upslope.area.sdat) %>%
      st_as_sf(
        as_points = FALSE,
        merge = T,
        na.rm = TRUE,
        use_integer = T) %>%
      st_union()
    
    catchment.wgs <- catchment.utm %>% 
      st_union() %>%
      st_buffer(dist = 80) %>%
      st_transform(4326)
    
    int.catch.dem <- st_intersects(catchment.wgs, dem.ext, sparse = F)[1, ]
    
    catchment.int <- length(which(int.catch.dem))
    
  }
  
  # Write catchment to disk
  
  catchment.loc <- paste0(lake.dir, "/", spec.lake.utm$UniqueID, "_catchment.gpkg")
  
  st_write(catchment.utm, catchment.loc, delete_dsn = T)
  
  # Obtain the catchment area with and without the lake in the catchment
  
  catchment.area <- st_area(catchment.utm) %>% 
    set_units("km^2")  
  
  catchment.area.nolake <- st_difference(catchment.utm, spec.lake.utm) %>%
    st_area(catchment.utm) %>% 
    set_units("km^2")
  
  # Convert the DEM to terra object.
  
  dem.utm.rast <- rast(raw.dem.sdat)
  
  # Extract median lake elevation
  
  lake.elev <- terra::extract(dem.utm.rast, 
                              vect(spec.lake.utm), 
                              fun = median, 
                              na.rm = T)
  
  # We also do some side-analysis that is not shown in the paper, including
  # relief analysis. Here, we create a 1-km buffer and mask the lake to extract  
  # the relief (largest elevation difference) around the lake.
  
  buffer.1km <- st_buffer(spec.lake.utm, 1000) %>%
    st_difference(., spec.lake.utm)
  
  surrounding.dem <- terra::mask(dem.utm.rast, vect(buffer.1km))
  
  max.elev.1km <- terra::extract(dem.utm.rast, 
                                 vect(buffer.1km), 
                                 fun = max, 
                                 na.rm = T)
  
  relief.1km <- max.elev.1km[ ,2] - lake.elev[ , 2]
  
  # Maximum elevation of the feeding catchment
  
  max.catchment <- terra::extract(dem.utm.rast, 
                                  vect(catchment.utm), 
                                  fun = max, 
                                  na.rm = T)
  
  relief.catchment <- max.catchment[ ,2] - lake.elev[ , 2]
  
  # Slope within 1 km and in the entire catchment
  
  slope <- terrain(dem.utm.rast, v = "slope")
  
  slope.1km <- mask(slope, vect(buffer.1km))
  
  slope.catchm <- mask(slope, vect(catchment.utm))
  
  slope.1km.val <- global(slope.1km, function(x) quantile(x, 0.5, na.rm = T)) %>% 
                    unname() %>%
                    as_vector()
  
  slope.catchm.val <- global(slope.catchm, function(x) quantile(x, 0.5, na.rm = T)) %>% 
                   unname() %>%
                   as_vector()
  
  # Angle of reach for all cells >30° slope.
  
  slope.gt30 <- slope < 30

  # Rasterize the lake polygon
  
  lake_raster <- rasterize(spec.lake.utm, rast(filled.dem.sdat), field = 1)
  
  # Create a distance raster (distance from lake boundary)
  
  dist_raster <- distance(lake_raster)
  
  # Extract lake boundary elevations
  
  lake_boundary_elevation <- mask(dem.utm.rast, lake_raster)
  
  # Elevation difference (DEM - lake boundary elevation)
  
  elevation_difference <- dem.utm.rast - lake.elev$mosaic

  # Calculate the angle of reach (in degrees)
  
  angle_of_reach <- atan(elevation_difference / dist_raster) * (180 / pi)
  
  angle_crop <- mask(angle_of_reach, slope.gt30, maskvalues = T) %>%
                 mask(spec.lake.utm, inverse = T) %>%
                 mask(vect(catchment.utm))
  
  # Sum of the angle of reach values for cells with slopes >30° according
  # to Allen et al. (2019)
  
  angle_sum <- global(angle_crop, 
                      fun = function(x) sum(atan(x), na.rm = T)) %>% 
                unname() %>%
                as_vector()
  
  # Glacier cover in the catchment
  # Intersect the catchment with RGI regions.
  
  rgi.reg.int <- st_intersects(catchment.wgs, o2_new, sparse = F)[1, ]
  
  o2.reg <- o2_new[rgi.reg.int, ]
  
  all.glaciers <- list.files(path = "D:/nrc_user/veh/RGI2000-v7.0-G-global",
                             pattern = "_valid.gpkg$", 
                             recursive = T,
                             full.names = T, 
                             include.dirs = T)
  
  # Find the glacier shapefiles in the RGI region.
  
  reg.glaciers <- all.glaciers[grep(o2.reg$o1region, basename(all.glaciers))]
  
  uni.regions <- unique(o2.reg$o1region)
  
  # Read these glacier shapefiles to memory.
  
  if (length(uni.regions) > 1) {
    
    several.regions <- lapply(uni.regions, function (x) {
      
      all.glaciers[grep(x, basename(all.glaciers))] %>%
        st_read() 
      
    })
    
    reg.glaciers <- do.call(rbind, several.regions)
    
    
  } else {
    
    matching.reg <- grep(uni.regions, basename(all.glaciers))
    
    reg.glaciers <- all.glaciers[matching.reg] %>%
      st_read() }
  
  reg.glaciers <- reg.glaciers %>% 
    st_make_valid() 
  
  sf_use_s2(FALSE)
  
  # Clip the glaciers with the catchment to obtain the glacier cover in the
  # catchment.
  
  glacier.int <- st_intersects(catchment.wgs, reg.glaciers, sparse = F)[1, ]
  
  if(length(which(glacier.int)) > 0) {
    
    # Find glaciers that intersect with the catchment
    
    reg.glaciers.utm <- reg.glaciers[glacier.int, ] %>%
      st_transform(st_crs(catchment.utm)) %>%
      st_make_valid()
    
    glaciers.clip <- st_intersection(catchment.utm, reg.glaciers.utm) %>% 
      st_make_valid()
    
    # If there is more than one glacier intersecting the catchment area...
    
    if(!(is.null(length(glaciers.clip)) | (length(glaciers.clip) == 0))) {
      
      # Write glaciers to disk
      
      glaciers.loc <- paste0(lake.dir, "/", spec.lake.utm$UniqueID, "_glaciers.gpkg")
      
      st_write(glaciers.clip, glaciers.loc, delete_dsn = T)
      
      # Calculate total glacier area.
      
      tot.glacier.area <- sum(st_area(glaciers.clip)) %>% 
        set_units("km^2")  
      
      glacier.cover <- (tot.glacier.area/ catchment.area) %>% 
        drop_units() 
      
      # Check if lake is in contact with glacier
      
      lake.glacier.int <- st_intersects(spec.lake.utm, glaciers.clip, sparse = F)
      
      if(any(lake.glacier.int[1, ])) {
        
        glacier.coupled <- "Yes"
        
        # Assess Euclidean distance between lake and glacier. We don't account
        # for the flow distance.
        
        lake.glacier.dist <- st_distance(spec.lake.utm, 
                                         glaciers.clip %>% 
                                           st_union(),
                                         which = "Euclidean") %>%
          drop_units() %>% as.numeric()
        
        if(length(as.numeric(lake.glacier.dist)) == 0) { lake.glacier.dist <- NA }
        
      }  else { glacier.coupled <- "No"
                 tot.glacier.area <- 0
                 glacier.cover <- 0
                 lake.glacier.dist <- NA }  } else {  glacier.coupled <- "No"
                                                      tot.glacier.area <- 0
                                                      glacier.cover <- 0
                                                      lake.glacier.dist <- NA } }
                      
  sf_use_s2(TRUE)
  
  # Assess landcover in catchment
  
  lc.int <- st_intersects(catchment.wgs %>% st_buffer(dist = 500), 
                          lc.ext, sparse = F)[1, ]
  
  map <- lc.ext$File[lc.int]
  
  # If the catchment is covered by several land cover maps, we need to generate
  # a land cover mosaic.
  
  if (length(map) > 1) {
  
    lc.mosaic <- stars::st_mosaic(map,
                                  dst = paste0(lake.dir, "/temp"),
                                  options = c("-resolution", "highest",
                                              "-r", "near"))
    lc.map <- rast(lc.mosaic)
  
  } else { lc.map <- rast(map) }
  
  # Save the original color table
  
  color_table <- coltab(lc.map)
  
  # Define the bbox of the catchment and the lake surrounding as an extent object.
  
  bbox.lake.buf <- st_bbox(buffer.1km %>% st_transform(4326))
  bbox.catchment <- st_bbox(catchment.utm %>% st_transform(4326))
  
  crop.extent.catchment <- ext(bbox.catchment)
  crop.extent.lake.buf <- ext(bbox.lake.buf)
  
  # Crop the raster and assign color table
  
  cropped.catchment <- crop(lc.map, crop.extent.catchment)
  levels(cropped.catchment) <- color_table

  masked.catchment <- mask(cropped.catchment,
                           spec.lake,
                           inverse = T) %>%
    mask(., catchment.utm %>% st_transform(4326) %>% vect())

  cropped.lake <- crop(lc.map, crop.extent.lake.buf)
  levels(cropped.lake) <- color_table

  masked.lake <- mask(cropped.lake, spec.lake, inverse = T) %>%
    mask(., buffer.1km %>% st_transform(4326) %>% vect())

  # Land cover legend according to the WorldCover metadata.

  legend <- data.frame(
    class_id = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100),
    land_cover = c(
      "Tree cover", "Shrubland", "Grassland", "Cropland",
      "Built-up", "Bare/sparse vegetation", "Snow and ice",
      "Permanent water bodies", "Herbaceous wetland",
      "Mangroves", "Moss and lichen"
    )
  )
  
  allNAlake  <- global(masked.lake, fun = function(x) all(is.na(x)))
  
  gc()
  
  allNAcatch <- global(masked.catchment, fun = function(x) all(is.na(x)))

  gc()
  
  # a few dozen lakes). If the surrounding area of a lake is NA, we write NA 
  # into the catchment stats.
  
  if(! (allNAlake$global | allNAcatch)) {
  
  max.lc.class.catchment <- freq(masked.catchment) %>%
    as_tibble() %>%
    mutate(frac = count/ sum(count)) %>%
    slice_max(frac,
              with_ties = F)

  lc.catchment <- legend$land_cover[legend$class_id %in% max.lc.class.catchment$value]

  max.lc.class.lake <- freq(masked.lake) %>%
    as_tibble() %>%
    mutate(frac = count/ sum(count)) %>%
    slice_max(frac,
              with_ties = F)

  lc.lake <- legend$land_cover[legend$class_id %in% max.lc.class.lake$value]

  # Save the cropped raster

  outfile.catchment <- paste0(lake.dir, "/",
                              basename(lake.dir),
                              "_landcover_masked_catchment.tif")

  outfile.lake <- paste0(lake.dir, "/",
                         basename(lake.dir),
                         "_landcover_masked_lake.tif")

  writeRaster(masked.catchment,
              outfile.catchment, overwrite = TRUE)

  writeRaster(masked.lake,
              outfile.lake, overwrite = TRUE) 
  
  } else { 
    
    max.lc.class.catchment <- data.frame(frac = NA)
    max.lc.class.lake <- data.frame(frac = NA)
    lc.lake <- NA
    lc.catchment <- NA
    
    }

  # Populate the lake shapefile with all statistics regarding catchment area,
  # relief, elevation, glacier area and cover, distance to and coupling with
  # glaciers, and land cover.
  
  spec.lake.with.stats <- spec.lake.utm %>%
    mutate(catchmarea   = catchment.area,
           catchmnolake = catchment.area.nolake,
           catchmrelief = relief.catchment,
           ElevationCOP = lake.elev[1, 2],
           relief1km    = relief.1km,
           glacierarea  = tot.glacier.area,
           glaciercov   = glacier.cover,
           glaciercoup  = glacier.coupled,
           dist_to_gl   = lake.glacier.dist,
           slope1km     = slope.1km.val,
           slopecatchm  = slope.catchm.val,
           aorsum       = angle_sum,
           domin_lake_lc = lc.lake,
           domin_lake_frac = max.lc.class.lake$frac,
           domin_catchm_lc = lc.catchment,
           domin_catchm_frac = max.lc.class.catchment$frac)

  # Write this shapefile to disk.
  
  stats.loc <- paste0(lake.dir, "/", spec.lake.utm$UniqueID, "_with_stats.gpkg")
  
  st_write(spec.lake.with.stats, dsn = stats.loc, delete_dsn = T)
  
  # Clip all SAGA grid files to the outline of the catchment to reduce the
  # amount of data.
  
  sdats <- list.files(dem.dir, pattern = ".sdat$",
                      full.names = T)
  
  for (j in 1: length(sdats)) {
    
    outfile <- paste0(str_sub(sdats[j], 1, -6), "_crop.sdat")
    
    f <- rast(sdats[j]) %>%
      crop(., catchment.utm)
    
    terra::writeRaster(f,  outfile, overwrite = T)
    
  }
  
  all.dem.files <- list.files(dem.dir, full.names = T)
  
  gc()
  
  files.to.delete <- grep(pattern = "crop", all.dem.files, invert = T, value = T)
  
  file.remove(files.to.delete)
  
  closeAllConnections()
  
}

doParallel::stopImplicitCluster()
gc()

# Catchments in Greenland can be covered both by peripherial glaciers and those 
# descending from the ice sheets. However, the ice sheets in Greenland are not 
# part of the RGI and therefore not part of the previous analysis. We beed to 
# update the glacier cover for Greenland accordingly.

# Load GrIS shapefile

gris <- st_read("D:/nrc_user/veh/glaciers_cci_gi_rgi05_TM-ETM_1994-2009_v170525/glaciers_cci_gi_greenland_gis-cl2_2000.shp") %>%
  st_transform(4326) %>%
  st_make_valid()

# List all lake shapefiles in the Greenland directory starting with RGI code "05"

greenland.dir <- list.dirs(path = "D:/nrc_user/veh/All_lakes", 
                           full.names = T, recursive = F) %>% 
                 grep(pattern = "/05-", value = T)

greenland.files <- list.files(pattern = "_with_stats.gpkg",
                               path = greenland.dir,
                               full.names = T, 
                               include.dirs = T,
                               recursive = T) 

# We iterate over all lakes in Greenland. 
# We load the catchment shapefile and clip the GrIS shapefile with
# the catchment shapefile, if they overlap.

registerDoParallel(45)

foreach(i = 1: length(greenland.files),
       .packages = c("tidyverse",
                     "sf",
                     "units")
       ) %dopar% {
         
  spec.lake.utm <- greenland.files[i] %>%
    st_read() %>%
    st_make_valid()
  
  # Write catchment to disk
  
  catchment.utm <- list.files(pattern = "catchment.gpkg", 
                              path = dirname(greenland.files[i]), 
                              full.names = T) %>%
    st_read() %>%
    st_make_valid()
  
  # change to WGS84
  
  catchment.wgs <- st_transform(catchment.utm, 4326) %>%
    st_make_valid()
  
  # Obtain the catchment area with and without the lake in the catchment
  
  catchment.area <- st_area(catchment.utm) %>% 
    set_units("km^2")  

  # Glacier cover in the catchment.
  # Intersect the catchment with RGI regions and the GrIS shapefile.
  
  gris.int <- st_intersects(catchment.wgs, gris, sparse = F)[1, ]
  
  gris.reg <- gris[gris.int, ]
  
  sf_use_s2(FALSE)
  
  # Clip the GrIS with the catchment to obtain the GrIS cover in the
  # catchment.
  
  if(nrow(gris.reg) > 0) {
    
    # Find glaciers that intersect with the catchment
    
    reg.gris.utm <- gris.reg %>%
      st_transform(st_crs(catchment.utm)) %>%
      st_make_valid()
    
    gris.clip <- st_intersection(catchment.utm, reg.gris.utm) %>% 
      st_make_valid()
    
    # If there is more than one glacier intersecting the catchment area...
    
    if(!(is.null(nrow(gris.clip)) | (nrow(gris.clip) == 0))) {
      
      # Write glaciers to disk
      
      gris.loc <- paste0(dirname(greenland.files[i]), 
                             "/", spec.lake.utm$UniqueID, "_gris.gpkg")
      
      st_write(gris.clip, gris.loc, delete_dsn = T)
      
    } 
  }
       }
    
doParallel::stopImplicitCluster()
gc()

# We then iterate over all lakes, and merge glaciers and GrIS shapefiles,
# if needed.

registerDoParallel(45)

foreach(i = 1 : length(greenland.files),
        .packages = c("tidyverse",
                      "sf",
                      "units")) %dopar% {
        
          # Read the lake shapefile              
                        
          stats.shp <- greenland.files[i] %>% st_read()
            
          gris.file <- list.files(pattern = "_gris.gpkg", 
                             path = dirname(greenland.files[i]), 
                             full.names = T)
          
          # Jump to next shapefile, if we hadn't found any overlap with the GrIS 
          # and its catchment in the preceding loop.
          
          if(length(gris.file) == 0) {
            
            return(NULL)} else {
              
              gris.shp <- st_read(gris.file) %>%
                st_make_valid()}
          
          # Also load the glacier shapefile, if available.
          
          glaciers.file <- list.files(pattern = "_glaciers.gpkg", 
                                  path = dirname(greenland.files[i]), 
                                  full.names = T)
          
          # Merge glacier and GrIS shapefile.
          
          if(length(glaciers.file) == 1) {
            
            glaciers.shp <- glaciers.file %>% st_read()
            
            if(nrow(glaciers.shp) >0) {gris.shp <- st_union(gris.shp, glaciers.shp)}}
          
          gris.shp <- gris.shp %>% summarise() 
          
          # Calculate the area of the GrIS and the glaciers.
          
          tot.gris.area <- sum(st_area(gris.shp)) %>% 
              set_units("km^2")  
            
          # Assess glacier cover in the catchment
          
          gris.cover <- (tot.gris.area/ stats.shp$catchmarea) %>% 
              drop_units() 
          
          # and the distance between glacier and the lake.
          
          lake.gris.dist <- st_distance(stats.shp, 
                                        gris.shp,
                                        which = "Euclidean") %>%
            drop_units() %>% as.numeric()
            
          if(length(as.numeric(lake.gris.dist)) == 0) { lake.gris.dist <- NA }
          
          # Check if lake is in contact with glacier
          
          lake.gris.int <- st_intersects(stats.shp, gris.shp, sparse = F)
          
          if(any(lake.gris.int[1, ])) { gris.coupled <- "Yes"
            
           }  else { gris.coupled <- "No"} 
            
          spec.lake.with.stats <- stats.shp %>%
            mutate(glacierarea  = tot.gris.area %>% drop_units(),
                   glaciercov   = gris.cover,
                   glaciercoup  = gris.coupled,
                   dist_to_gl   = lake.gris.dist)
          
          st_write(spec.lake.with.stats, greenland.files[i], delete_dsn = T)
          
          }

doParallel::stopImplicitCluster()
gc()

# Bind all shapefiles with stats to one large shapefile. To this end, first
# list all files, then read them to memory in parallel.

lake.stats.files <- list.files(pattern = "with_stats.gpkg",
                               path = "D:/nrc_user/veh/All_lakes",
                               full.names = T, 
                               include.dirs = T,
                               recursive = T)

cl <- makeCluster(40) # Adapt to the number of available cores on YOUR machine.

# Export R-packages to the cluster, which will be used in the apply loop. Make sure
# you installed these packages before.

clusterEvalQ(cl = cl, c(library("sf"),
                        library("tidyverse")))

# Export the file list of the glacier shapefiles and the extent of the study
# regions to the clusters.

clusterExport(cl = cl, list("lake.stats.files"))

lake.stats <- pblapply(lake.stats.files, cl = cl, function (x) {
  
  f <- st_read(x) %>% st_drop_geometry()
  
  return(f)
  
})

stopCluster(cl)

lake.stats.sf <- do.call(rbind, lake.stats)

saveRDS(lake.stats.sf, "D:/nrc_user/veh/Zhang_glacial_lakes_global/lakes_with_catchment_and_landcover.rds")
