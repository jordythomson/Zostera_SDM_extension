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
fishnet <- st_read("data/shapefile/sg_sdm_fishnet.shp")





# 2. Aggregate eelgrass occurrence observations within cells and then thin spatially ====

# Read in eelgrass occurrence records
sg_records <- read.csv("data/csv/eelgrass_occurrence_records_20240118.csv") %>%
  dplyr::select(-starts_with("sampling")) %>% 
  mutate(SG = as.numeric(seagrass)) %>% # convert species column to numeric
  filter(!is.na(latitude) & !is.na(longitude)) %>% # filter out any points with no lat/lon
  filter(id != 4692 & id != 4710 & id != 4705 & id != 4700 & id != 4699 & id != 4691 & id != 4692 & id != 4669 & id != 4668 & id != 4667 & id != 4666 & id != 4665) %>% # filter out two problematic observations where the point is over a shallow bank in the corner of the cell but the centroid is in a deep channel
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(26920)

# Join eelgrass records with fishnet grid
joined_grid <- st_join(fishnet, sg_records) %>%
  filter(!is.na(percent.cover)) # filter out those without % cover (abundance) data

# Aggregate multiple observations in the same grid cell
join_gridAgg <- joined_grid %>% 
  group_by(GridID) %>% # group by grid cell ID
  summarise(subObs_int = raster::modal(subObs_JT, na.rm = TRUE), # modal substrate value
            seagrass = sum(seagrass), # any SG presence = present
            cover = mean(percent.cover, na.rm = TRUE), 
            freq = n()) %>% # frequency of observations per grid cell
  mutate(seagrass = if_else(seagrass > 0, 1, 0)) # rescale P-A to 0-1
str(join_gridAgg) # 2258 % cover observations

# Write joined grid to shapefile
st_write(obj = join_gridAgg,
         dsn = "data/shapefile",
         layer = "SG_cover_pop_grid_all",
         driver = "ESRI Shapefile",
         append = FALSE)

# The next few steps are setting up a dataframe for spatial thinning

# Covert aggregated grid back to lat/lon, create columns for centroid lat and lon, and add a species column as required by spThin
join_gridAgg2 <- st_transform(join_gridAgg, crs = 4326) %>%
  st_centroid(.) %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  mutate(species = rep("eelgrass", length(join_gridAgg$GridID)))

# Add column with row names (case ID) to facilitate left join after thinning
join_gridAgg2$id <- as.numeric(rownames(join_gridAgg2))
str(join_gridAgg2)

# Spatially thin grid cells using NND of 200 m
thinned <-
  thin(loc.data = join_gridAgg2,
       lat.col = "lat", long.col = "lon",
       spec.col = "species", 
       thin.par = 0.2, reps = 100, # set 200 m nearest neighbour distance, repeat process 100 times because there is a random element; proceed with data set that retains the most observations
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

str(thinned2) # 1067 cases retained
thinned2

write.csv(thinned2, "output/thinned_locs_cover_200m_20240118.csv")

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
         layer = "SG_cover_pop_grid_200m",
         driver = "ESRI Shapefile",
         append = FALSE)

# At this stage, I went to QGIS to join SG_cover_pop_grid_200m to fishnet grid by location
# From there, I rasterized the layer using % cover data and exported
# That produces separate raster layers for % cover that can be loaded and stacked with predictors
# This was repeated with 100 m NND distance threshold





# 3. Update, stack, and crop environmental predictor layers==== 

# Bathymetry
dem <- raster("data/raster/bathy_yw_orig_v2_adj0_fieldvals.tif")

# Seabed slope
slope <- raster("data/raster/slope_20230920_JGS.tif") 

# Broad scale bathymetric position index
bpi_broad <- raster("data/raster/bathy_Broad_BPI_v2_20230920_JGS.tif") 

# Fine scale bathymetric position index
bpi_fine <- raster("data/raster/bathy_Fine_BPI_20230920_JGS.tif")

# Relative wave exposure index
rei <- raster("data/raster/rei_scaled.tif")

# Substrate type
substrate <- raster("data/raster/substrate.tif") # this is a raster created from the in situ observation point data in the eelgrass occurrence records file (soft, mixed, rock)

# Ratify substrate layer
# substrate <- as.factor(substrate)
# substrate <- ratify(substrate)
# levels(substrate)
# rat <- levels(substrate)[[1]]
# rat$substrate <- c("NA", "soft", "mixed", "rocky")
# levels(substrate) <- rat

temp1m_median <- raster("data/raster/temp1m_median.tif") # median of daily median temperature, June 1st - Sep 15th, 2021

temp1m_qrange <- raster("data/raster/temp1m_qrange.tif") # mean of daily difference between 90th and 10th quantiles of temperature, June 1st - Sep 15th, 2021

ws1m_q90 <- raster("data/raster/ws1m_q90.tif") # mean of daily 90th quantile of water velocity, June 1st - Sep 15th, 2021

# Stack predictors, crop to study domain, and rename
predictors <- stack(dem, slope, bpi_broad, bpi_fine, rei, substrate) 
names(predictors) <- c("bathy_m", "slope", "bpi_broad", "bpi_fine", "rei", "substrate")
str(predictors)

# Load eelgrass PA and % cover rasters (rasters created via spatial join of attributes from join_gridAgg2 with fishnet grid in QGIS, then conversion to raster)
# join_gridAgg2 was created by 1) aggregating all data in fishnet grid cells (PA = any presence, % cover = mean) and then spatially thinning using a NND threshold of 200 m and 100 repititions
PA <- raster("data/raster/eelgrass_PA_raster_200m_20240118.tif") 
cover <- raster("data/raster/eelgrass_cover_raster_200m_20240118.tif")

# Load ensemble model predicted habitat suitability layer
suitability10 <- raster("data/raster/proj_10_EM_CurrentZostera_projection_ZosteraMarina_ensemble.tif")
suitability12 <- raster("data/raster/proj_12_EM_CurrentZostera_projection_ZosteraMarina_ensemble.tif")

# Create raster stack
all_data <- stack(PA, cover, suitability10, suitability12, dem, slope, bpi_broad, bpi_fine, substrate, rei, temp1m_median, temp1m_qrange, ws1m_q90)
names(all_data) <- c("PA", "cover", "suitability9", "suitability12", "bathy_m", "slope", "bpi_broad", "bpi_fine", "sub_type", "rei", "temp1m_median", "temp1m_qrange", "ws1m_q90")
str(all_data)

# Convert all_data stack to dataframe
# read join_gridAgg2 from shapefile because spatial thinning has already been done and has a random element to it (starting point), which means it would differ slightly if run every time
join_gridAgg2 <- st_read("data/shapefile/SG_cover_pop_grid_200m.shp")

# Extract grid cell centroid coordinates
coords <- st_centroid(join_gridAgg2)
str(coords)

# Extract data from all_data_df associated with each set of coordinates and put into dataframe
all_data_df <- coords %>% 
  raster::extract(all_data, ., df = TRUE) %>% 
  dplyr::select(-ID) %>%
  cbind(., coords$lat, coords$lon)
str(all_data_df)

write.csv(all_data_df, "data/csv/PA_cover_suitability_predictors_200m.csv")

# At this point, I went to the csv for some quick quality control to check for any PA 1s that have 0 % cover or 0s that have some cover
# Deleted two points
# Need to automate this step!

# Re-load data
all_data_df <- read.csv("data/csv/PA_cover_suitability_predictors_200m.csv")

# Create data frame for cases with % cover and complete data
cover_df <- drop_na(all_data_df)
str(cover_df)
View(cover_df)

cover_df$biomass <- 23.549243 + 0.118732*cover_df$cover + 0.028087*(cover_df$cover^2)
str(cover_df)

write.csv(cover_df, "data/csv/cover_dataframe_200m_20240118.csv")

# Create dataframe for presence-only cover/abundance
pres_cover_df <- filter(cover_df, cover > 0)
pres_cover_df

write.csv(pres_cover_df, "data/csv/pres_cover_dataframe_200m_20240118.csv") # export to modify for beta regression by adding cover.prop (cover/100) and cover.prop.br (1s replaced with 0.99s)
# Need to automate this step!

# Re-load data
pres_cover_df <- read.csv("data/csv/pres_cover_dataframe_200m_20240118.csv")





# Data exploration --------------------------------------------------------

# histograms of abundance indices and habitat suitability for all data including zeros
par(mfrow = c(2,2))
hist(cover_df$cover, xlab = "% cover", main = "")
hist(cover_df$biomass, xlab = "biomass (g DW m-2)", main = "")
hist(cover_df$suitability9, xlab = "habitat suitability [model 9]", main = "")
hist(cover_df$suitability12, xlab = "habitat suitability [model 9]", main = "")


# histograms of abundance indices and habitat suitability for presence data only
par(mfrow = c(2,2))
hist(pres_cover_df$cover, xlab = "% cover", main = "")
hist(pres_cover_df$biomass, xlab = "biomass (g DW m-2)", main = "")
hist(pres_cover_df$suitability9, xlab = "habitat suitability [model 8]", main = "")
hist(pres_cover_df$suitability12, xlab = "habitat suitability [model 8]", main = "")


# Plot biomass-cover relationship just to verify that function looks correctly used
dev.off()
plot(biomass ~ cover, data = cover_df)


# PA and cover
str(pres_cover_df)
corr.data1 <- dplyr::select(pres_cover_df, cover, biomass, suitability9) 

library(psych)
pairs.panels(corr.data1,
             smooth = TRUE,       # If TRUE, draws loess smooths
             scale = FALSE,       # If TRUE, scales the correlation text font
             density = FALSE,     # If TRUE, adds density plots and histograms
             ellipses = FALSE,    # If TRUE, draws ellipses
             method = "pearson",  # Correlation method (also "spearman" or "kendall")
             pch = 16,            # pch symbol
             lm = FALSE,          # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,          # If TRUE, reports correlations
             jiggle = FALSE,      # If TRUE, data points are jittered
             factor = 2,          # Jittering factor
             hist.col = 4,        # Histograms color
             stars = TRUE,        # If TRUE, adds significance level with stars
             ci = TRUE)  

str(pres_cover_df)
corr.data2 <- dplyr::select(pres_cover_df, cover, biomass, sub_type, bathy_m, rei, slope, bpi_broad, bpi_fine) 

pairs.panels(corr.data2,
             smooth = TRUE,       # If TRUE, draws loess smooths
             scale = FALSE,       # If TRUE, scales the correlation text font
             density = FALSE,     # If TRUE, adds density plots and histograms
             ellipses = FALSE,    # If TRUE, draws ellipses
             method = "pearson",  # Correlation method (also "spearman" or "kendall")
             pch = 16,            # pch symbol
             lm = FALSE,          # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,          # If TRUE, reports correlations
             jiggle = FALSE,      # If TRUE, data points are jittered
             factor = 2,          # Jittering factor
             hist.col = 4,        # Histograms color
             stars = TRUE,        # If TRUE, adds significance level with stars
             ci = TRUE)  

str(pres_cover_df)
corr.data3 <- dplyr::select(pres_cover_df, cover, biomass, temp1m_median, temp1m_qrange, ws1m_q90) 

pairs.panels(corr.data3,
             smooth = TRUE,       # If TRUE, draws loess smooths
             scale = FALSE,       # If TRUE, scales the correlation text font
             density = FALSE,     # If TRUE, adds density plots and histograms
             ellipses = FALSE,    # If TRUE, draws ellipses
             method = "pearson",  # Correlation method (also "spearman" or "kendall")
             pch = 16,            # pch symbol
             lm = FALSE,          # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,          # If TRUE, reports correlations
             jiggle = FALSE,      # If TRUE, data points are jittered
             factor = 2,          # Jittering factor
             hist.col = 4,        # Histograms color
             stars = TRUE,        # If TRUE, adds significance level with stars
             ci = TRUE)  



# 1. Two-stage abundance-based SDM ----------------------------------------
str(pres_cover_df)

par(mfrow = c(2,2))
plot(cover ~ suitability9, data = cover_df)
plot(biomass ~ suitability9, data = cover_df)
plot(cover ~ suitability12, data = cover_df)
plot(biomass ~ suitability12, data = cover_df)

par(mfrow = c(2,2))
plot(cover ~ suitability9, data = pres_cover_df, pch = 16, xlab = "habitat suitability [model 9]", ylab = "% cover")
plot(biomass ~ suitability9, data = pres_cover_df)
plot(cover ~ suitability12, data = pres_cover_df, pch = 16, xlab = "habitat suitability [model 9]", ylab = "% cover")
plot(biomass ~ suitability12, data = pres_cover_df)


ggplot(pres_cover_df, aes(x = suitability9, y = cover, col = bathy_m, label = X)) +
  geom_point() +
  geom_text(nudge_y = 1) +
  labs(x = "Habitat suitability [model 9]", y = "% cover")

ggplot(pres_cover_df, aes(x = suitability12, y = cover, col = bathy_m, label = X)) +
  geom_point() +
  geom_text(nudge_y = 1) +
  labs(x = "Habitat suitability [model 12]", y = "% cover")



library(betareg)

m1 <- betareg(cover.prop.br ~ suitability9, link = "log", data = pres_cover_df)
summary(m1)
plot(m1)

m2 <- betareg(cover.prop.br ~ suitability12, link = "log", data = pres_cover_df)
summary(m2)
plot(m2)






# 2. Abundance when present -----------------------------------------------
# Beta regression
plot(cover ~ temp1m_qrange, data = pres_cover_df)
plot(cover ~ temp1m_median, data = pres_cover_df)
plot(cover ~ ws1m_q90, data = pres_cover_df)
plot(cover ~ ws1m_q90, data = subset(pres_cover_df, ws1m_q90 < 0.8))

m1 <- betareg(cover.prop.br ~ bathy_m*rei + slope + sub_type + temp1m_qrange + bathy_m + bpi_broad + bpi_fine + temp1m_median + temp1m_qrange + ws1m_q90, link = "log", data = pres_cover_df)
summary(m1)

ggplot(pres_cover_df, aes(x = bathy_m, y = cover.prop.br, col = rei)) +
  geom_point() +
  labs(x = "Depth (m)", y = "% cover")

m1b <- update(m1, cover.prop ~ . - bpi_broad)
summary(m1b)

m1c <- update(m1b, cover.prop ~ . - slope)
summary(m1c)

m1d <- update(m1c, cover.prop ~ . - sub_type)
summary(m1d)

m1e <- update(m1d, cover.prop ~ . - temp1m_qrange)
summary(m1e)

m1f <- update(m1e, cover.prop ~ . - bpi_fine)
summary(m1f)

m1g <- update(m1f, cover.prop ~ . - temp1m_median)
summary(m1g)

m1h <- update(m1g, cover.prop ~ . - ws1m_q90)
summary(m1h)

m1i <- update(m1h, cover.prop ~ . - rei)
summary(m1i)


plot(cover.prop ~ bathy_m, data = pres_cover_df, pch = 16, xlab = "depth (m)", ylab = "% cover")


# Linear model on log-transformed abundance data
m1 <- lm(log(biomass) ~ bathy_m + slope + bpi_broad + bpi_fine + rei + temp1m_median + temp1m_qrange + ws1m_q90, data = pres_cover_df)
summary(m1)

par(mfrow = c(2,2))
plot(m1)

# Random forest
# https://favtutor.com/blogs/random-forest-in-r#:~:text=Random%20Forest%20for%20Regression%20in%20R%201%20Step,...%207%20Step%207%3A%20Evaluate%20the%20Model%20





# 3. Abundance-absence ----------------------------------------------------
# glm

glm1 <- glm(biomass ~ bathy_m + slope + bpi_broad + bpi_fine + rei + temp1m_median + temp1m_qrange + ws1m_q90, family = Gamma(link = "inverse"), data = pres_cover_df)
summary(glm1)

glm2 <- glm(cover ~ bathy_m + slope + bpi_broad + bpi_fine + rei + temp1m_median + temp1m_qrange + ws1m_q90, family = Gamma(link = "inverse"), data = pres_cover_df)
summary(glm2)

# hurdle model
library(pscl)
hurdle1 <- hurdle(cover ~ bathy_m + slope + bpi_broad + bpi_fine + rei + temp1m_median + temp1m_qrange + ws1m_q90, data = cover_df, dist = c("poisson"), zero.dis = "negbin")

# random forest
library(randomForest)








