# Project Description and Code Overview ====
#
# Project:      Eelgrass species distribution modelling
# Description:  Creating a high spatial resolution ensemble species distribution
#               model for eelgrass (Zostera marina) to support conservation
#               and marine spatial planning in Canadian Maritimes 

# Contact:      e-mail: John.Obrien@dfo-mpo.gc.ca | tel: +1.782.640.1522
# Publication:  
#
# Code Overview:
# 1. Create fishnet grid for study domain from Cape Sable Island to Cape North
# 2. Aggregate direction observations of eelgrass occurrence and environmental
#    conditions to spatial resolution of predictor layers
# 3. Update, stack, and crop environmental predictor layers
# 4. Format data for biomod2
# 5. Split data into cross validation fold using spatial blocking
#
# Requirements:
# R version 3.6.3

# House Keeping ====

# Load required packages

library(dplyr)
library(readr)
library(sf)
library(raster)
library(fasterize)
library(blockCV)
library(biomod2)
library(usdm)
library(ggplot2)
library(terra)
library(spThin)

# set raster processing options
rasterOptions(chunksize = 1e+05, maxmemory = 1e+09)





# 1. Load fishnet grid ====
fishnet <- st_read("data/shapefile/sg_sdm_fishnet.shp")





# 2. Aggregate eelgrass occurrence observations within cells and then thin spatially ====

# Read in eelgrass occurrence records
sg_records <- read.csv("data/csv/eelgrass_occurrence_records_20240117.csv") %>%
  dplyr::select(-starts_with("sampling")) %>% 
  mutate(SG = as.numeric(seagrass)) %>% # convert species column to numeric
  filter(!is.na(latitude) & !is.na(longitude)) %>% # filter out any points with no lat/lon
  filter(id != 4692 & id != 4710) %>% # filter out two problematic observations where the point is over a shallow bank in the corner of the cell but the centroid is in a deep channel
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(26920)

# Join eelgrass records with fishnet grid
joined_grid <- st_join(fishnet, sg_records) %>%
  filter(!is.na(seagrass)) 

# Aggregate multiple observations in the same grid cell
join_gridAgg <- joined_grid %>% 
  group_by(GridID) %>% # group by grid cell ID
  summarise(subObs_int = raster::modal(subObs_JT, na.rm = TRUE), # modal substrate value
            seagrass = sum(seagrass), # any SG presence = present
            cover = mean(percent.cover, na.rm = TRUE), 
            depth = mean(bathy_m*-1, na.rm = TRUE), # invert depths to be consistent with bathymetry raster
            freq = n()) %>% # frequency of observations per grid cell
  mutate(seagrass = if_else(seagrass > 0, 1, 0)) # rescale P-A to 0-1
str(join_gridAgg)

# The next few steps are setting up a dataframe for spatial thinning

# Covert aggregated grid back to lat/lon, create columns for centroid lat and lon, and add a species column as required by spThin
join_gridAgg2 <- st_transform(join_gridAgg, crs = 4326) %>%
  st_centroid(.) %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  mutate(species = rep("eelgrass", length(join_gridAgg$GridID)))
str(join_gridAgg2)

# Add column with row names (case ID) to facilitate left join after thinning
join_gridAgg2$id <- as.numeric(rownames(join_gridAgg2))
str(join_gridAgg2)

# Spatially thin grid cells using NND of 200 m
thinned <-
  thin(loc.data = join_gridAgg2,
       lat.col = "lat", long.col = "lon",
       spec.col = "species", 
       thin.par = 0.2, reps = 100, # set 200 m nearest neighbour distance, repeat process 100 times because there is a random element; proceed with data set that retains the most observations within the NND constraint
       locs.thinned.list.return = TRUE,
       write.files = FALSE,
       out.dir = "output/",
       write.log.file = TRUE,
       log.file = "spatial_thin_log.txt",
       verbose = TRUE)

# Save list of retained cells as data frame
thinned2 <- as.data.frame(thinned[1]) %>%
  mutate(id = as.integer(rownames(.))) %>% # add row names as id column
  rename(., lon = Longitude, lat = Latitude) %>%
  dplyr::select(-lat, -lon)

str(thinned2) # 1064 cases retained
thinned2

write.csv(thinned2, "output/thinned_locs_200m_model12_20240117.csv")

# Left join list of retained cells with variables in join_gridAgg2
join_gridAgg2 <- left_join(thinned2, join_gridAgg2, by = "id")
join_gridAgg2  
str(join_gridAgg2)

# Covert back to sf object 
join_gridAgg2 <- st_as_sf(join_gridAgg2)
str(join_gridAgg2)

# Transform to UTM
join_gridAgg2 %>% 
 st_transform(26920)
str(join_gridAgg2)

# Write joined grid to shapefile
st_write(obj = join_gridAgg2,
         dsn = "data/shapefile",
         layer = "12_SG_pop_grid",
         driver = "ESRI Shapefile",
         append = FALSE)

# Read joined grid
join_gridAgg <- st_read("data/shapefile/12_SG_pop_grid.shp") %>%
  st_transform(26920)





# 3. Update, stack, and crop environmental predictor layers==== 

# Bathymetry and derivitives

# Bathymetry
# dem <- raster("data/raster/bathy_yw_orig_v2_adj0.tif")

# Replace cell values with observed depth values if available
# cells <- join_gridAgg %>%
#  filter(!is.na(depth)) %>% 
#  st_centroid() %>% 
#  st_coordinates() %>% 
#  cellFromXY(dem, .)
#  dem_vals <- values(dem)
#  new_vals <- filter(join_gridAgg, !is.na(depth)) %>% 
#  pull(depth)
#  dem_vals[cells] <- new_vals
#  dem <- setValues(dem, dem_vals)

# Write new raster to file
# writeRaster(dem, filename = "data/raster/bathy_yw_orig_v2_adj0_fieldvals.tif")

dem <- raster("data/raster/bathy_yw_orig_v2_adj0_fieldvals.tif")

# Seabed slope
slope <- raster("data/raster/slope_20230920_JGS.tif") 

# Broad scale bathymetric position index
bpi_broad <- raster("data/raster/bathy_Broad_BPI_v2_20230920_JGS.tif") 

# Fine scale bathymetric position index
bpi_fine <- raster("data/raster/bathy_Fine_BPI_20230920_JGS.tif")

# Relative wave exposure index
rei <- raster("data/raster/rei_scaled.tif")

# Stack predictors, crop to study domain, and rename
predictors <- stack(dem, slope, bpi_broad, bpi_fine, rei) 
names(predictors) <- c("bathy_m", "slope", "bpi_broad", "bpi_fine", "rei")
str(predictors)





# 4. Format data for biomod2 ====

# Name of the study species
myRespName <- "ZosteraMarina"

# Combine values of predictors for each sampling point into a dataframe

coords <- st_centroid(join_gridAgg) # coordinates of grid cells with eelgrass data

myExpl <- coords %>% 
  raster::extract(predictors, ., df = TRUE) %>% 
  dplyr::select(-ID)
dim(myExpl) # 1064

# Fill missing predictor values
slope_NA <- is.na(myExpl$slope)
rei_NA <- is.na(myExpl$rei)

myExpl[slope_NA, "slope"] <- raster::extract(predictors$slope, coords[slope_NA,], buffer = 100, fun = mean)
myExpl[rei_NA, "rei"] <- raster::extract(predictors$rei, coords[rei_NA,], buffer = 100, fun = mean)
dim(myExpl) # 1064

# Presence/absences data for our species 
myResp <- join_gridAgg$seagrass[!is.na(myExpl$slope) & !is.na(myExpl$rei) & !is.na(myExpl$bathy_m)]
length(myResp) # 1046
str(myResp)

# XY coordinates of species data
myRespXY <- st_centroid(join_gridAgg) %>% 
  st_coordinates()
myRespXY <- myRespXY[!is.na(myExpl$slope) & !is.na(myExpl$rei) & !is.na(myExpl$bathy_m),]
length(myRespXY) # 2092
str(myRespXY)

# biomod2 format
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = filter(myExpl, !is.na(slope) & !is.na(rei) & !is.na(bathy_m)),
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)
str(myBiomodData) # 1046

# Save biomod2 formatted data
saveRDS(myBiomodData, file = "data/rds/12_biomod2_data.rds")

# Check for collinearity among predictors

# Pairwise correlations - |rho| < 0.70
cor_spear <- cor(x = myBiomodData@data.env.var, 
                 use = "pairwise.complete",
                 method = "spearman") # Spearman's Correlation

diag(cor_spear) <- NA
apply(abs(cor_spear), MARGIN = 2, FUN = max, na.rm = TRUE)

# Variance Inflation Factors - VIF < 3
vif_pred <- vif(myBiomodData@data.env.var)
range(vif_pred$VIF) # 1.28 - 1.65





# 5. Split data into CV fold using spatial blocking ====

# Examine effective range of spatial autocorrelation in predictors

# original function, now deprecated
# sac <- spatialAutoRange(rasterLayer = dropLayer(predictors, 5), 
#                        sampleNumber = 5000,
#                        doParallel = TRUE,
#                        showPlots = TRUE)
# summary(sac) # mean = 29240

# attempt at using cv_spatial_autocor, as the deprecation message directed
sac <- cv_spatial_autocor(r = predictors, 
                          num_sample = 5000,
                          progress = TRUE, 
                          plot = FALSE) 
sac$range # 20399

# spatial blocking by range of spatial autocorrelation with random assignment
sb <- spatialBlock(speciesData = join_gridAgg[!is.na(myExpl$slope) & !is.na(myExpl$rei) & !is.na(myExpl$bathy_m),],
                   species = "seagrass",
                   rasterLayer = predictors,
                   theRange = sac$range, # size of the blocks
                   k = 5,
                   selection = "random",
                   iteration = 100, # find evenly dispersed folds
                   biomod2Format = TRUE,
                   xOffset = 0, # shift the blocks horizontally
                   yOffset = 0) 
DataSplitTable <- sb$biomodTable
head(DataSplitTable)
dim(DataSplitTable) # 1483

# Save CV partition
saveRDS(DataSplitTable, file = "data/rds/12_CV_folds_biomod.rds")

# Plot spatial blocks
sb_plots <- ggplot() +
  sb$plots$layers[[1]] + 
  geom_sf(data = st_centroid(join_gridAgg),
          alpha = 0.5,
          size = 3,
          shape = 21,
          col = "black",
          fill = "darkgoldenrod1") +
  sb$plots$layers[[2]] +
  geom_sf_label(data = sb$plots$layers[[3]]$data,
               aes(label = folds),
               # label.r = unit(0.75,"lines"),
               alpha = 0.5,
               colour = NA) +
  geom_sf_text(data = sb$plots$layers[[3]]$data,
                aes(label = folds),
                fontface = "bold") +
  labs(x = "", y = "") +
  theme_bw() +
  scale_fill_viridis_c(guide = "none") +
  theme(text = element_text(colour = "red"),
        axis.text = element_text(colour = "black"))

sb_plots

ggsave(sb_plots, filename = "Output/Spatial_blocks.tiff",
       compression = "lzw", dpi = 600)

# Save spatial blocks
saveRDS(sb, file = "Data/SpatialBlocks.rds")

# Move on to model fitting and evaluation
cat("Proceed to 2ModelBuilding.R")