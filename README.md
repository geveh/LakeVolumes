# **LakeVolumes**: Code for *"Evolving resource potential of glacial lakes with ongoing deglaciation"*

## Overview

**This repository contains eight scripts to estimate the local, regional, and global volume of glacial lakes from a given lake area. In addition, we provide protocols to assess the resource potential of glacial lakes, including their distance to, and elevation above, the coast; their potential lifetimes associated with sedimentation; the number of people living above and below catchments hosting glacial lakes.**

- [01_Lake_types.R](#01_lake_typesr)
- [02_VA_model.R](#02_va_modelr)
- [03_global_lake_volume.R](#03_global_lake_volumer)
- [04_Catchment_statistics.R](#04_catchment_statisticsr)
- [05_Lake_volumes_glaciers_elev_dist.R](#05_lake_volumes_glaciers_elev_distr)
- [06_Lifetimes.R](#06_lifetimer)
- [07_population_estimates.R](#07_population_estimatesr)
- [process_flowpaths_detailed.m](#process_flowpaths_detailedm)


The codes are written in the statistical programming language **R** (https://www.r-project.org/), Version 4.2.2, and called within
the Graphical User Interface **RStudio** (https://posit.co/download/rstudio-desktop/) under a Microsoft Windows Server 2019 operating system. 
Please install both R and RStudio on your machine to successfully run the codes and produce figures and R data objects.

The R codes depend on a number of packages, listed at the beginning of all scripts. Please install those packages before running the scripts. 
The comments within the scripts provide further details on model dependencies and usage of functions. 

Each script will call one or more input data object(s), which are available via ***Zenodo***.  
We also use freely available digital elevation models (DEMs), land cover maps, glaciological data (e.g., outlines of glaciers and ice sheets), and glacial lake outlines. Sources for these datasets are given in the script where they are needed. Please download the data before executing the scripts.  
Please put all input files into the same folder, and change the working directory (set at the beginning of each script) to your folder structure. The scripts can be executed one after another, with the each script generating output that is used as input for the next script.
The scripts (and parts thereof) can also be run independent of each other using the input files (in most cases *.RDS* files) from Zenodo.
Each script will produce output in form of a figure (displayed in the associate manuscript and supplementary figures) or R-data objects.

## Scripts

### 01_Lake_types.R

**Script to assign a dam type to each lake in the glacial lake inventory of Zhang et al. (2024) as of 1990 and 2020.**

*Mandatory input data*: 
- All regional glacial lake shapefiles from Zhang et al. (2024), available at https://doi.org/10.11888/Cryos.tpdc.300938
- The Randolph Glacier Inventory (RGI) V7.0, available at https://doi.org/10.5067/f6jmovy5navz
- "RGI2000-v7.0-o2regions_modified_densified.gpkg" (Modified outlines of the O2 regions in the RGI V7.0)

*Main outputs*: 
- "all_lakes_including_dam_type.rds" (R-object of all glacial lakes with an assigned dam type)
- "glacier_lakes_2020_centroids.shp" (R-object of all glacial lake centroids for the year 2020, used to derive downstream flow paths)
- "Lakes19902020_damtype.RDS" (R-object of all glacial lakes in 1990 and 2020 with a dam type assigned)

---

### 02_VA_model.R

**Script to fit a Bayesian hierarchical linear regression model of lake area *A* versus volume *V* distinguished by dam type.**

*Mandatory input data*: 
- "va.txt" (Text file of bathymetrically surveyed lakes including their volume *V* and area *A* (converted from the Excel Sheet)) 
- "HDIofMCMC.R" (R-function to estimate the highest density interval for a given distribution, written by John K. Kruschke, available in this zip folder: https://drive.google.com/file/d/1rkpJC148LB4Hi7K-1IRKkqyED3PGWKeo/view)

*Main outputs*: 
- "VA_model.RDS" (R-object with linear regression model of *V* versus *A* distinguished by dam type)
- "model_parameters.pdf" (Summary of model output as a figure; Extended Data Figure 2)
- "va_model_one_panel.pdf" (PDF figure showing all V-A data with posterior regression estimates; Figure 1a)
- "Trend_types_and_posterior.pdf" (PDF figure showing the trend in *V*-*A* per lake type, and the posterior regression)
- "VA_model_errors.pdf" (PDF figure showing prediction errors of the *V*-*A* model compared to the original data that entered the model; Extended Data Figure 3)
- "VA_data.RDS" (R-object of all non-repetitively surveyed glacial lakes)

---

### 03_global_lake_volume.R

**Script to estimate the local, regional, and global volume of glacial lakes.**

*Mandatory input data*: 
- "RGI2000-v7.0-o1regions.shp" (O1 regions of the RGI to aggregate lake volumes on regional level from https://doi.org/10.5067/f6jmovy5navz)
- "VA_model.RDS" (R-object with linear regression model of *V* versus *A* distinguished by dam type)
- "VA_data.RDS" (R-object of all non-repetitively surveyed glacial lakes)
- "continent_dissolve.shp" (Shapefile of dissolved continent outlines from ArcGIS Hub: https://hub.arcgis.com/datasets/esri::world-continents/about)
- "HDIofMCMC.R" (R-function to estimate the highest density interval for a given distribution, written by John K. Kruschke, available in this zip folder: https://drive.google.com/file/d/1rkpJC148LB4Hi7K-1IRKkqyED3PGWKeo/view)

*Output*: 
- "Lakes19902020_damtype.RDS" (R-object of glacial lakes mapped by Zhang et al. (2024) in 1990 and 2020 with a posterior median and 68% HDI estimate of their volume)
- "all_lakes_with_volumes.RDS" (R-object with posterior median and 68% volume estimate for each lake in 1990 and 2020)
- "Regional_size_distribution.pdf" (PDF figure showing the empirical exceedance probabilities of lake volumes; Extended Data Figure 9)
- "volume_regional_lakes.RDS" (R-object with a posterior estimate of regional glacial lake volumes, including the median and 68% HDI of lake volumes)

---

### 04_Catchment_statistics.R

**Script to download the Copernicus GLO30 DEM and land cover maps. The script then derives the catchment of each lake, and summarises catchment-wide statistics including glacial cover, land cover, and relief.**

*Mandatory input data*: 
- "RGI2000-v7.0-o2regions.shp" (Original outlines of the O2 regions in the RGI V7.0, https://doi.org/10.5067/f6jmovy5navz)
- "RGI2000-v7.0-o2regions_modified_densified.gpkg" (Modified outlines of the O2 regions in the RGI V7.0)
- "Lakes19902020_damtype.RDS" (R-object of glacial lakes mapped by Zhang et al. (2024) in 1990 and 2020 with a posterior median and 68% HDI estimate of their volume)
- RGI2000-v7.0-G-global (All glaciers in the RGI V7.0, https://doi.org/10.5067/f6jmovy5navz)
- Copernicus GLO30 and GLO90 DEMs (30- and 90-m resolution DEMs derived from the TanDEM-X mission, available through Microsoft Planetary Computer: https://planetarycomputer.microsoft.com/api/stac/v1/)
- "glaciers_cci_gi_greenland_gis-cl2_2000.shp" (outlines of the Greenland ice sheet in 2000, available at https://glaciers-cci.enveo.at/crdp2/index.html)

*Output*:
- "lakes_with_catchment_and_landcover.rds" (R-object of all glacial lakes as of 2020 including statistics on relief, as well as glacier and land cover)

---

### 05_Lake_volumes_glaciers_elev_dist.R

**Script to estimate trends of glacier volume (change) with glacial lake volume (change), and to assess flow path distances from glacial lakes to the coast.**

*Mandatory input data*: 
- "Millan_glacier_volumes.xlsx" (Excel-file of regional glacier volumes, extracted from Millan et al. (2022): https://www.nature.com/articles/s41561-021-00885-z
- "Hugonnet_glacier_loss.txt" (Text-file glacier with mass loss data from Hugonnet et al. (2021), available at https://doi.org/10.6096/13)
- "RGI2000-v7.0-o1regions.shp" (O1 regions of the RGI to aggregate lake volumes on regional level)
- "all_lakes_with_volumes.RDS" (R-object with posterior median and 68% volume estimate for each lake in 1990 and 2020)
- "lakes_with_catchment_and_landcover.rds" (R-object with catchment-wide statistics on relief, glacier and land cover for all lakes in 2020)
- "flowpaths_detailed.txt" (Irregular txt file provided by Wolfgang Schwanghart. The file is organised by lines, with every lake is identified by 5 lines: its ID, and the X, Y, distance, and Z coordinate from the source towards the coast. Every fifth line, a new flow path starts from a new lake starts.)

*Output*:
- "lake_and_glacier_volume_change.pdf" (PDF-figure showing trends of lake volume (change) with glacier mass (change); Extended Data Figure 5)
- "parameter_ice_vol_loss_fits.pdf" (PDF-figure showing the performance of the two models estimating glacier volume (change) from glacial lake (change)
- "cumulative_distances_and_elevations.pdf" (PDF-figure showing the global cumulated volume of glacier lakes with distance to the coast/ elevation; Figure 3)
- "cumulative_distances_and_elevations_regional.pdf" (PDF-figure showing the regional cumulated volume of glacier lakes with distance to the coast/ elevation; Extended Data Figure 6)

---

### 06_Lifetime.R

**Script to approximate the local, regional, and global lifetime (longevity) of glacial lakes.**

*Mandatory input data*: 
- "RGI2000-v7.0-o1regions.shp" (O1 regions of the RGI to aggregate regional lake lifetimes)
- "lakes_with_catchment_and_landcover.rds" (R-object with catchment-wide statistics on relief, glacier and land cover for all lakes in 2020)
- "VA_data.RDS" (R-object of all non-repetitively surveyed glacial lakes)
- "all_lakes_with_volumes.RDS" (R-object with posterior median and 68% volume estimate for each lake in 1990 and 2020)
- "continent_dissolve.shp" (Shapefile of dissolved continent outlines)
- "HDIofMCMC.R" (R-function to estimate the highest density interval for a given distribution, written by John K. Kruschke, available in this zip folder: https://drive.google.com/file/d/1rkpJC148LB4Hi7K-1IRKkqyED3PGWKeo/view)
- "sciadv.adr2009_data_s1_and_s2/adr2009_data_s1.xlsx": glacial erosion rates from Wilner et al. (2024), available at https://www.science.org/doi/10.1126/sciadv.adr2009
- "sciadv.adr2009_data_s1_and_s2/adr2009_data_s2.xlsx": fluvial erosion rates from Wilner et al. (2024), available at https://www.science.org/doi/10.1126/sciadv.adr2009

*Output*: 
- "simulated_infill_times.RDS" (R-object with simulated lifetimes of each glacial lake under fluvial and glacial erosion rates, and a scenario weighted by glacial cover in a given catchment)
- "density_lifetime.pdf" (PDF figure showing the stacked density of estimated individual lake lifetimes; Figure 4a)
- "remaining_storage.pdf" (PDF figure showing the sedimentation-driven storage loss of glacial lakes; Figure 4b)
- "Regional_lake_lifetimes.pdf" (PDF figure showing the cumulative estimated regional lifetime, i.e. year when x% of the regional glacial lake volume will be lost; Extended Data Figure 8)
- 7 small panels to show the workflow to estimate the sedimentation-driven lifetime of glacial lakes for Extended Data Figure 7 ("small_panel_catchment_size.pdf"; "small_panel_erosion_rate.pdf"; "small_panel_glac_cov.pdf"; "small_panel_rock_density.pdf"; "small_panel_sediment_prod.pdf"; "small_panel_sediment_depo.pdf"; "small_panel_lake_volume.pdf")
- "volume_vs_lifetime.pdf" (PDF figure showing the median volume + 68% HDI vs. simulated median sedimentation-driven lifetime + 68% HDI of all glacial lakes)  
- "landcover.pdf" (PDF figure showing the regional landcover in catchments feeding glacial lakes; Figure 5)

---

### 07_population_estimates.R

**Script to extract population in the catchment feeding glacial lakes, and within the first 50 km downstream of glacial lakes**

*Mandatory input data*: 
- "all_catchments_dissolved_no_holes.gpkg" (Geopackage with all catchments feeding glacial lakes, dissolved, derived from COP30 DEM data)
- "all_lakes_with_volumes.RDS" (R-object with posterior median and 68% volume estimate for each lake in 1990 and 2020)
- "flowpaths_detailed.txt" (Flowpaths from all lakes as of 2020 from the lake to the coast/ endorheic basin)
- "volume_regional_lakes.RDS" (R-object with a posterior estimate of regional glacial lake volumes, including the median and 68% HDI of lake volumes)
- "diff.gpkg" (Geopackage that contains the difference between catchments and a dissolved 1-km buffer around flowpaths originating from glacial lakes)
- "RGI2000-v7.0-o1regions.shp" (O1 regions of the RGI to aggregate lake volumes on regional level from https://doi.org/10.5067/f6jmovy5navz)
- "continent_dissolve.shp" (Shapefile of dissolved continent outlines from ArcGIS Hub: https://hub.arcgis.com/datasets/esri::world-continents/about)
- "landscan-global-2020.tif" (Landscan2020 data corresponding to the year when lakes were mapped: https://landscan.ornl.gov/)
- "ppp_2020_1km_Aggregated.tif" (Worldpop 1km global: https://hub.worldpop.org/geodata/listing?id=64)
- "gpw_v4_population_count_rev11_2020_30_sec.tif" (NASA's Gridded Population of the World (GPW), version 4; https://www.earthdata.nasa.gov/data/projects/gpw/data-access-tools)
- "GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0.tif" (Global Human Settlement Layer (GHSL); https://human-settlement.emergency.copernicus.eu/download.php)

*Output*:
- "flowpath_buffers.gpkg" (1km buffer along the first 50-km flow paths downstream of glacial lakes)
- "all_catchments.gpkg" / "all_catchments.gpkg" (Geopackage / R-object containing all catchments feeding glacial lakes worldwide in 2020)
- "population_upstream_downstream.pdf" (PDF figure showing little barplots of the regional population living in the upstream catchment and along a 1-km wide, 50-km long buffer downstream of lakes)  
- "World_map.pdf" (PDF figure showing the major RGI regions, used in Figure 6a)  
- "Regional_volume_vs_population.pdf" (PDF figure showing a comparison of regional population versus regional lake volume upstream and downstream of lakes; Figure 6b,c)

---

### process_flowpaths_detailed.m

**Matlab Script by Wolfgang Schwanghart to calculate flow path from glacial lake centroids to the coast or endorheic basins. Please contact Wolfgang Schwanghart (wolfgang.schwanghart@uni-potsdam.de), if you have further questions.**

*Mandatory input data*: 
- "merit500/DEM.mat": a composite of the Merit DEM (https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/) at 500-m pixel resolution.
- "merit500/FD.mat": Flow directions calculated from this Merit DEM
- "lakes_2020.gpkg": a point vector layer with the centroids of all lakes mapped by Zhang et al. (2024) in 2020.

*Output*: 
- "flowpaths_detailed.txt" (Unstructured Textfile with coordinates, length, and elevation of edges along the flow path from a lake to the sea/ endorheic basin)

## Input data

Please visit the repository on Zenodo to obtain the input files.


## References

Georg Veh, Wolfgang Schwanghart, Oliver Korup, and Jonathan L. Carrivick: *Evolving resource potential of glacial lakes with ongoing deglaciation*. Nature Water, In Revision.

## Contact

**Georg Veh**  
Postdoctoral researcher in the research group on natural hazards  
Institute of Environmental Sciences and Geography  
University of Potsdam  
georg.veh@uni-potsdam.de  
https://www.uni-potsdam.de/de/umwelt/forschung/ag-naturgefahren.html
