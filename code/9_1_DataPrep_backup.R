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
library(tidyverse)

# set raster processing options
rasterOptions(chunksize = 1e+05, maxmemory = 1e+09)





# 1. Load fishnet grid ====

# fishnet <- st_read("data/shapefile/sg_sdm_fishnet.shp")





# 2. Aggregate eelgrass occurrence observations within cells and then thin spatially ====

# Read in eelgrass occurrence records
# sg_records <- read.csv("data/csv/eelgrass_occurrence_records_20230918.csv") %>%
#  dplyr::select(-starts_with("sampling")) %>% 
#  mutate(SG = as.numeric(seagrass)) %>% # convert species column to numeric
#  filter(!is.na(latitude) & !is.na(longitude)) %>%
#  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
#  st_transform(26920)

# Join eelgrass records with fishnet grid
# joined_grid <- st_join(fishnet, sg_records) %>%
#  filter(!is.na(seagrass)) 

# Aggregate multiple observations in the same grid cell
# join_gridAgg <- joined_grid %>% 
#  group_by(GridID) %>% # group by grid cell ID
#  summarise(subObs_int = raster::modal(subObs_JT, na.rm = TRUE), # modal substrate value
#            seagrass = sum(seagrass), # any SG presence = present
#            cover = mean(percent.cover, na.rm = TRUE), 
#            freq = n()) %>% # frequency of observations per grid cell
#  mutate(seagrass = if_else(seagrass > 0, 1, 0)) # rescale P-A to 0-1
# str(join_gridAgg)

# Covert back to lat/lon
# join_gridAgg2 <- st_transform(join_gridAgg, crs = 4326) %>%
#  st_centroid(.) %>%
#  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
#                lat = sf::st_coordinates(.)[,2]) %>%
#  mutate(species = rep("eelgrass", 3396))
# str(join_gridAgg2)

# Add column with row names (case ID) to facilitate left join after spatial thinning
# join_gridAgg2$id <- as.numeric(rownames(join_gridAgg2))
# str(join_gridAgg2)





# Spatially thin grid cells using NND of 200 m
# thinned <-
#  thin(loc.data = join_gridAgg2,
#       lat.col = "lat", long.col = "lon",
#       spec.col = "species", 
#       thin.par = 0.2, reps = 100, # set 200 m nearest neighbour distance, repeat process 100 times because there is a random element; proceed with data set that retains the most observations within the NND constraint
#       locs.thinned.list.return = TRUE,
#       write.files = FALSE,
#       out.dir = "output/",
#       write.log.file = TRUE,
#       log.file = "spatial_thin_log.txt",
#       verbose = TRUE)

# Save list of retained cells as data frame
# thinned2 <- as.data.frame(thinned[1]) %>%
#  mutate(id = as.integer(rownames(.))) %>% # add row names as id column
#  rename(., lon = Longitude, lat = Latitude) %>%
#  dplyr::select(-lat, -lon)

# str(thinned2) # 1062 cases retained
# thinned2

# write.csv(thinned2, "output/thinned_locations_200m_20231206.csv")

# Left join list of retained cells with variables in join_gridAgg2
# join_gridAgg2 <- left_join(thinned2, join_gridAgg2, by = "id")
# join_gridAgg2  
# str(join_gridAgg2)

# Covert back to sf object 
# join_gridAgg2 <- st_as_sf(join_gridAgg2)
# str(join_gridAgg2)

# Transform to UTM
# join_gridAgg2 %>% 
#  st_transform(26920)
# str(join_gridAgg2)
# View(join_gridAgg2)

# Write joined grid to shapefile
# st_write(obj = join_gridAgg2,
#         dsn = "data/shapefile",
#         layer = "9_SG_pop_grid",
#         driver = "ESRI Shapefile",
#         append = FALSE)





# 3. Stack eelgrass data, environmental predictor layers, habitat suitability and coordinates ==== 

# Bathymetry
dem <- raster("data/raster/bathy_yw_orig_v2_aligned_clipped_recalc_ext_20230810.tif")

# Seabed slope
slope <- raster("data/raster/slope_20230920_JGS.tif") 

# Broad scale bathymetric position index
bpi_broad <- raster("data/raster/bathy_Broad_BPI_v2_20230920_JGS.tif") 

# Fine scale bathymetric position index
bpi_fine <- raster("data/raster/bathy_Fine_BPI_20230920_JGS.tif")

# REI
rei <- raster("data/raster/rei_scaled.tif")

# substrate <- raster("data/raster/substrate_reclassified_20230911_mergeCHS_manualupdate.tif")
substrate <- raster("data/raster/substrate.tif") # this is a raster created from the in situ observation point data in the eelgrass occurrence records file (soft, mixed, rock)

# Ratify substrate layer
# substrate <- as.factor(substrate)
# substrate <- ratify(substrate)
# levels(substrate)
# rat <- levels(substrate)[[1]]
# rat$substrate <- c("NA", "soft", "mixed", "rocky")
# levels(substrate) <- rat

# Median of daily median temperature
temp1m_median <- raster("data/raster/temp1m_median.tif") # median of daily median temperature, June 1st - Sep 15th, 2021

# Median of daily difference between 90th and 10th quantiles
temp1m_qrange <- raster("data/raster/temp1m_qrange.tif") # mean of daily difference between 90th and 10th quantiles of temperature, June 1st - Sep 15th, 2021

# Median of daily 90th quantile of water speed
ws1m_q90 <- raster("data/raster/ws1m_q90.tif") # mean of daily 90th quantile of water velocity, June 1st - Sep 15th, 2021

# Stack predictors, crop to study domain, and rename
# predictors <- stack(dem, slope, bpi_broad, bpi_fine, substrate, rei) 
# names(predictors) <- c("bathy_m", "slope", "bpi_broad", "bpi_fine", "sub_type", "rei")
# str(predictors)

# Load eelgrass PA and % cover rasters (rasters created via spatial join of attributes from join_gridAgg2 with fishnet grid in QGIS, then conversion to raster)
# join_gridAgg2 was created by 1) aggregating all data in fishnet grid cells (PA = any presence, % cover = mean) and then spatially thinning using a NND threshold of 200 m and 100 repititions
seagrass_PA <- raster("data/raster/eelgrass_PA_raster_model9.tif") 
seagrass_cover <- raster("data/raster/eelgrass_cover_raster_model9.tif")

# Load ensemble model predicted habitat suitability layer
suitability <- raster("data/raster/proj_8_EM_CurrentZostera_projection_ZosteraMarina_ensemble.tif")
suitability10 <- raster("data/raster/proj_10_EM_CurrentZostera_projection_ZosteraMarina_ensemble.tif")

# Create raster stack
all_data <- stack(seagrass_PA, seagrass_cover, suitability, suitability10, dem, slope, bpi_broad, bpi_fine, substrate, rei, temp1m_median, temp1m_qrange, ws1m_q90)
names(all_data) <- c("seagrass_PA", "seagrass_cover", "suitability", "suitability10", "bathy_m", "slope", "bpi_broad", "bpi_fine", "sub_type", "rei", "temp1m_median", "temp1m_qrange", "ws1m_q90")
str(all_data)

# Convert all_data stack to dataframe
# read join_gridAgg2 from shapefile because spatial thinning has already been done and has a random element to it (starting point), which means it would differ slightly if run every time
join_gridAgg2 <- st_read("data/shapefile/9_SG_pop_grid.shp")

# Extract grid cell centroid coordinates
coords <- st_centroid(join_gridAgg2)
str(coords)

# Extract data from all_data_df associated with each set of coordinates and put into dataframe
all_data_df <- coords %>% 
  raster::extract(all_data, ., df = TRUE) %>% 
  dplyr::select(-ID) %>%
  cbind(., coords$lat, coords$lon)
str(all_data_df)

# write.csv(all_data_df, "data/csv/all_data_df.csv")

# Create data frame for cases with % cover and complete data
cover_df <- drop_na(all_data_df)
str(cover_df)

cover_df$biomass <- 23.549243 + 0.118732*cover_df$seagrass_cover + 0.028087*(cover_df$seagrass_cover^2)
str(cover_df)

write.csv(cover_df, "data/csv/cover_dataframe.csv") # export and exclude data points with presence but 0 % cover; this was a mistake in creation of cover data frame (NAs ended up being included as 0s)
### cover_dataframe2.csv, loaded below, excludes NAs for % cover

cover_df2 <- read.csv("data/csv/cover_dataframe2.csv") # load dataframe with only positive % cover data
str(cover_df2)





##### Data exploration
hist(cover_df2$seagrass_cover)
hist(cover_df2$biomass)
hist(cover_df2$suitability)
hist(cover_df2$suitability10)

# Plot biomass-cover relationship just to verify that function looks correctly used
plot(biomass ~ seagrass_cover, data = cover_df2)



# 1. Two-stage SAM: OLS using habitat suitability as a predictor variable


par(mfrow = c(2,1))
plot(seagrass_cover ~ suitability10, data = cover_df2)
plot(seagrass_cover ~ suitability, data = cover_df2)
plot(biomass ~ suitability10, data = cover_df2)

hist(cover_df2$seagrass_cover)

m1 <- lm(seagrass_cover ~ suitability10, data = cover_df2)
summary(m1)
par(mfrow = c(2,2))
plot(m1)

m1 <- lm(biomass ~ suitability10, data = cover_df2)
summary(m1)
par(mfrow = c(2,2))
plot(m1)

library(betareg)
m2 <- betareg(seagrass_cover_prop ~ suitability10, data = cover_df2)
summary(m2)


# Histograms of % cover when present and habitat suitability
hist(pos$seagrass_cover)
hist(pos$suitability)

par(mfrow = c(2,2))
plot(seagrass_cover ~ suitability, data = pos, xlab = "Habitat suitability", ylab = "% cover")
plot(seagrass_cover ~ temp1m_median, data = pos, xlab = "Median temperature", ylab = "% cover")
plot(seagrass_cover ~ temp1m_qrange, data = pos, xlab = "90th-10th temperature quantile", ylab = "% cover")
plot(seagrass_cover ~ ws1m_q90, data = pos, xlab = "90th quantile water speed", ylab = "% cover")

par(mfrow = c(2,2))
plot(biomass ~ suitability, data = pos, xlab = "Habitat suitability", ylab = "% cover")
plot(biomass ~ temp1m_median, data = pos, xlab = "Median temperature", ylab = "% cover")
plot(biomass ~ temp1m_qrange, data = pos, xlab = "90th-10th temperature quantile", ylab = "% cover")
plot(biomass ~ ws1m_q90, data = pos, xlab = "90th quantile water speed", ylab = "% cover")

m1 <- lm(seagrass_cover ~ suitability, data = pos)

summary(m1)

par(mfrow = c(2,2))
plot(m1)

m2 <- lm(biomass ~ suitability, data = pos)
summary(m2)



# still need to apply the grouped function



# 2. Abundance when present SAM: beta regression



# 3. Abundance-absence model
str(cover_df)

m1 <- glm(seagrass_cover ~ suitability10, data = )
summary(m1)

par(mfrow = c(2,2))
plot(m1)




# 4. Format data for biomod2 ====

# Name of the study species
myRespName <- "ZosteraMarina"

# Combine values of predictors for each sampling point into a dataframe

coords <- st_centroid(join_gridAgg) # coordinates of grid cells with eelgrass data

myExpl <- coords %>% 
  raster::extract(predictors, ., df = TRUE) %>% 
  dplyr::select(-ID)
dim(myExpl) # 1049

# Fill missing predictor values
slope_NA <- is.na(myExpl$slope)
rei_NA <- is.na(myExpl$rei)
# sub_type_NA <- is.na(myExpl$sub_type)

myExpl[slope_NA, "slope"] <- raster::extract(predictors$slope, coords[slope_NA,], buffer = 100, fun = mean)
myExpl[rei_NA, "rei"] <- raster::extract(predictors$rei, coords[rei_NA,], buffer = 100, fun = mean)
# myExpl[sub_type_NA, "sub_type"] <- raster::extract(predictors$sub_type, coords[sub_type_NA,], buffer = 80, fun = raster::modal)
dim(myExpl) # 1049

# Presence/absences data for our species 
myResp <- join_gridAgg$seagrass[!is.na(myExpl$slope) & !is.na(myExpl$rei) & !is.na(myExpl$bathy_m) & !is.na(myExpl$sub_type)]
length(myResp) # 1030
str(myResp)

# XY coordinates of species data
myRespXY <- st_centroid(join_gridAgg) %>% 
  st_coordinates()
myRespXY <- myRespXY[!is.na(myExpl$slope) & !is.na(myExpl$rei) & !is.na(myExpl$bathy_m) & !is.na(myExpl$sub_type),]
length(myRespXY) # 2060
str(myRespXY)

# biomod2 format
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = filter(myExpl, !is.na(slope) & !is.na(rei) & !is.na(bathy_m) & !is.na(myExpl$sub_type)),
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)
str(myBiomodData) # 1030

# Save biomod2 formatted data
saveRDS(myBiomodData, file = "data/rds/8_biomod2_data.rds")

# Check for collinearity among predictors

# Pairwise correlations - |rho| < 0.70
cor_spear <- cor(x = myBiomodData@data.env.var, 
                 use = "pairwise.complete",
                 method = "spearman") # Spearman's Correlation

diag(cor_spear) <- NA
apply(abs(cor_spear), MARGIN = 2, FUN = max, na.rm = TRUE)

# Variance Inflation Factors - VIF < 3
vif_pred <- vif(myBiomodData@data.env.var)
range(vif_pred$VIF) # 1.28 - 1.92





# 5. Split data into CV fold using spatial blocking ====

# Examine effective range of spatial autocorrelation in predictors

# original function, now deprecated
# sac <- spatialAutoRange(rasterLayer = dropLayer(predictors, 5), 
#                        sampleNumber = 5000,
#                        doParallel = TRUE,
#                        showPlots = TRUE)
# summary(sac) # mean = 29240

# attempt at using cv_spatial_autocor, as the deprecation message directed
predictors.sub <- stack(dem, slope, bpi_broad, bpi_fine, rei) # using only continuous predictors with complete raster layers

sac <- cv_spatial_autocor(r = predictors.sub, 
                          num_sample = 5000,
                          progress = TRUE, 
                          plot = FALSE) 
sac$range # 16871

# spatial blocking by range of spatial autocorrelation with random assignment
sb <- spatialBlock(speciesData = join_gridAgg[!is.na(myExpl$slope) & !is.na(myExpl$rei) & !is.na(myExpl$bathy_m) & !is.na(myExpl$sub_type),],
                   species = "seagrass",
                   rasterLayer = predictors,
                   theRange = sac$range, # size of the blocks; 41772 m
                   k = 5,
                   selection = "random",
                   iteration = 100, # find evenly dispersed folds
                   biomod2Format = TRUE,
                   xOffset = 0, # shift the blocks horizontally
                   yOffset = 0) 
DataSplitTable <- sb$biomodTable
head(DataSplitTable)
dim(DataSplitTable) # 1030

# Save CV partition
saveRDS(DataSplitTable, file = "data/rds/8_CV_folds_biomod.rds")

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