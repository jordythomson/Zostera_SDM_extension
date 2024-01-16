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
# 1. Fit 7 individual sdm model types to data
# 2. Combine individual sdm's to generate ensemble model
# 3. Calculate variable importance
# 4. Project predictions of individual models over study domain
# 5. Generate ensemble prediction over study domain

# Requirements:
# R version 3.6.3
# biomod2 formatted eelgrass occurrence and predictor data
# CV folds of data
# Environmental predictor layers

# House Keeping ====

# Load required packages
library(biomod2)
library(dplyr)
library(tidyr)
library(data.table)
library(forcats)
library(purrr)
library(raster)
library(sf)
library(terra)
library(doBy)

# list of categorical predictors
facVars <- c("sub_type")

# set raster processing options
rasterOptions(chunksize = 1e+05, maxmemory = 1e+09)





# 1. Fit individual sdm's ====

# Import formatted biomod2 data and CV folds

myBiomodData <- readRDS("data/rds/11_biomod2_data.rds") # biomod2 data
str(myBiomodData) # 823

DataSplitTable <- readRDS("data/rds/11_CV_folds_biomod.rds")
dim(DataSplitTable) # 823

# Define Models Options using default options

myBiomodOption <- BIOMOD_ModelingOptions(GLM = list(type = "quadratic",
                                                    interaction.level = 0),
                                         ANN = list(maxit = 500),
                                         MARS = list(type = "simple"))

myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData, 
                                     models = c("RF", "ANN", "GBM", "FDA", "GLM", "CTA", "MARS"), 
                                     bm.options = myBiomodOption, 
                                     data.split.table = DataSplitTable,
                                     var.import = 10, 
                                     metric.eval = c("ROC","TSS"),
                                     do.full.models = FALSE,
                                     modeling.id = "AtlCoast_AUC_0.7")

# Playing around with output summary functions in biomod2
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c("expl.var", "algo", "run"))
bm_PlotResponseCurves(bm.out = myBiomodModelOut)
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo'), dataset = "validation")
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo'), dataset = "calibration")
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'), dataset = "validation")
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'), dataset = "calibration")
bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = "validation")
bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = "calibration")

# get cv evaluations

# John's original code
# eval <- get_evaluations(myBiomodModelOut,as.data.frame=T) %>% 
#  mutate(Model.type = strsplit(x = as.character(Model.name), split = "_"),
#         Model.type = map_chr(Model.type, first)) %>%  
#  group_by(Model.type, Eval.metric) %>%   
#  summarise(across(c("Testing.data", "Sensitivity", "Specificity"), 
#                   list(mean = mean, sd = sd), 
#                   .names = "{col}_{fn}",
#                   na.rm = TRUE))

# Updated code
eval <- get_evaluations(myBiomodModelOut) %>% 
  group_by(algo, metric.eval) %>%   
  summarise(across(c("validation", "sensitivity", "specificity"), # Verify with John that "Testing.data" from orig code should be "validation", not "calibration"
               list(mean = mean, sd = sd), 
               .names = "{col}_{fn}",
               na.rm = TRUE))
eval

saveRDS(eval, file = "Output/11_AllMods_evalStats.rds")

# Get variable importance
var_impALL <- get_variables_importance(myBiomodModelOut)
var_impALL

# Summarize
imp_df <- summaryBy(var.imp ~ expl.var + algo, FUN = c(min, mean, max), data = var_impALL)
imp_df
# Re-order for easy comparison with John's table
imp_df[order(imp_df$var.imp.mean),]

# Save importance values
saveRDS(imp_df, file = "Output/11_VarImp_ALL.rds")





# 2. Doing Ensemble Modelling ====
myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                       models.chosen = "all",
                                       em.by = "PA+run", # changed from em.by = "PA_dataset+repet"
                                       em.algo = "EMwmean",
                                       metric.select = c("ROC"),
                                       metric.select.thresh = c(0.7),
                                       metric.eval = c("TSS","ROC"),
                                       # prob.mean = FALSE,
                                       # prob.cv = FALSE,
                                       # prob.ci = FALSE,
                                       # prob.median = FALSE,
                                       # committee.averaging = FALSE,
                                       # prob.mean.weight = TRUE,
                                       # prob.mean.weight.decay = 'proportional',
                                       EMwmean.decay = "proportional",
                                       var.import = 10)   

# print summary
myBiomodEM

# Updated code
eval_em <- get_evaluations(myBiomodEM) %>% 
  mutate(name = "Ensemble") %>% 
  group_by(name, metric.eval) %>%   
  summarise(across(c("validation", "sensitivity", "specificity"), # changed Testing.data to validation again
                   list(mean = mean, sd = sd), 
                   .names = "{col}_{fn}",
                   na.rm = TRUE))
eval_em

# Save evaluations
saveRDS(eval_em, file = "Output/11_Ensemble_evalStats.rds")

# Updated code
var_impEM <- get_variables_importance(myBiomodEM)
var_impEM

table_var_impEM <- summaryBy(var.imp ~ expl.var + merged.by.run, FUN = c(mean), data = var_impEM) %>%
  pivot_wider(., names_from = merged.by.run, values_from = var.imp.mean) %>%
  rowwise() %>% 
  mutate(.,
         max = max(c_across(starts_with("RUN"))),
         mean = mean(c_across(starts_with("RUN"))),
         min = min(c_across(starts_with("RUN"))))

table_var_impEM

# Save importance values
saveRDS(table_var_impEM, file = "Output/11_VarImp_Ensemble.rds")





# 4. Project current distribution ====

# a. Load predictors
dem <- raster("data/raster/bathy_yw_orig_v2_aligned_clipped_recalc_ext_20230810.tif")
# substrate <- raster("data/raster/substrate_reclassified_20230911_mergeCHS_manualupdate.tif")
# substrate <- as.factor(substrate)
# substrate <- ratify(substrate)
# rat <- levels(substrate)[[1]]
# levels(substrate)
# rat$substrate <- c("NA", "soft", "mixed", "rock")
# levels(substrate) <- rat
slope <- raster("data/raster/slope_20230920_JGS.tif")
bpi_broad <- raster("data/raster/bathy_Broad_BPI_v2_20230920_JGS.tif")
bpi_fine <- raster("data/raster/bathy_Fine_BPI_20230920_JGS.tif")
rei <- raster("data/raster/rei_scaled.tif")


# land polygon with 5-km buffer
land_buffer <- st_read("data/shapefile/coast50k_buff_5km.shp")

# Mask for Bay of Fundy
BoF_mask <- st_read("data/shapefile/BoF_mask.shp")

# Create mask for study area
atl_mask <- raster("data/raster/NS_12m_contour.tif") %>% 
  mask(., land_buffer) %>% # exclude cells > 5 km from shore
  mask(., BoF_mask, inverse = TRUE)

# Create fishnet grid from reclassified raster template
predictors <- stack(dem, slope, bpi_broad, bpi_fine, rei) %>% 
  mask(., mask = atl_mask)

names(predictors) <- c("bathy_m","slope","bpi_broad", # rename layers
                       "bpi_fine", "rei")

# writeRaster(predictors, filename = "data/raster/8_PredictorStack_mask.tif", overwrite = FALSE)





# b. Project current distribution for all individual models and CV runs
myBiomodProjection <- BIOMOD_Projection(myBiomodModelOut,
                                        new.env = stack(predictors),
                                        proj.name = "CurrentZostera_projection",
                                        models.chosen = "all",
                                        metric.binary = NULL,
                                        build.clamping.mask = FALSE, 
                                        compress = TRUE)

get_built_models(myBiomodModelOut) # use this command to get list of fitted model names

# Outputs predictions into CurrentZostera sub-folder, file = proj_CurrentZostera_ZosteraMarina.tif
# This .tif file has 35 bands (one for each model); these are added to QGIS project and through symbology can display one at a time
# However, the individual_projections sub-folder is empty... this is where I expected individual model .tif files to show up

str(myBiomodProjection)
myBiomodProjection





# 5. Project current distribution predicted from ensemble ====
myEnsembleProjection <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                                   bm.proj = myBiomodProjection,
                                                   proj.name = "10_EM_CurrentZostera_projection",
                                                   models.chosen = "all",
                                                   metric.binary = c("ROC","TSS"),
                                                   compress = TRUE)

get_built_models(myBiomodEM)

?BIOMOD_EnsembleForecasting

cat("Proceed to 3FigurePrep.R")