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

myBiomodData <- readRDS("data/rds/3_biomod2_data.rds") # biomod2 data
str(myBiomodData) # 2602

DataSplitTable <- readRDS("data/rds/3_CV_folds_biomod.rds")
dim(DataSplitTable) # 2602

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

saveRDS(eval, file = "Output/AllMods_evalStats.rds")

# Load John's original output for comparison
John_eval <- readRDS("Output/AllMods_evalStats_JohnOrig.rds")
John_eval
eval

# Get variable importance
var_impALL <- get_variables_importance(myBiomodModelOut)
var_impALL

# imp_list <- list()
# for (i in 1:5) {
#   imp_list[[i]] <- var_impALL[,,i,1]
# }
# imp_df <- map(imp_list, as.data.frame) %>% 
#  map(., ~cbind(predictor = row.names(.), .)) %>% 
#  rbindlist(., idcol = "CV_run") %>% 
#  pivot_longer(., cols = RF:MARS, values_to = "var_imp", names_to = "Model.type") %>%  
#  group_by(Model.type, predictor) %>% 
#  summarise(across(var_imp, 
#                   list(max = max, 
#                        mean = mean, 
#                        min = min),
#                   .names = "{fn}")) %>% 
#  ungroup() %>% 
#  mutate(predictor = forcats::as_factor(predictor)) %>%
#  arrange(mean) %>% 
#  mutate(predictor = fct_reorder(predictor, mean))

# Save importance values
# saveRDS(imp_df, file = "Output/VarImp_ALL.rds")

# John's original variable importance summary table
VarImp_ALL_JohnOrig <- readRDS("Output/VarImp_ALL_JohnOrig.rds")
print(VarImp_ALL_JohnOrig, n = 42)

# My attempt to recreate it
imp_df <- summaryBy(var.imp ~ expl.var + algo, FUN = c(min, mean, max), data = var_impALL)
imp_df
# Re-order for easy comparison with John's table
imp_df[order(imp_df$var.imp.mean),]

# Save importance values
saveRDS(imp_df, file = "Output/3_VarImp_ALL.rds")







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

# get evaluation scores

# John's original code
# eval_em <- get_evaluations(myBiomodEM, as.data.frame = TRUE) %>% 
#  mutate(name = "Ensemble") %>% 
#  group_by(name, Eval.metric) %>%   
#  summarise(across(c("Testing.data", "Sensitivity", "Specificity"), 
#                   list(mean = mean, sd = sd), 
#                   .names = "{col}_{fn}",
#                   na.rm = TRUE))

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
saveRDS(eval_em, file = "Output/3_Ensemble_evalStats.rds")






# 3. Get variable importance ====

# John's original code
# var_impEM <- get_variables_importance(myBiomodEM)
# imp_list <- list()
# for (i in 1:5) {
# imp_list[[i]] <- var_impEM[,,i]
# }
# imp_df <- map(imp_list, as.data.frame) %>% 
#   map(., ~cbind(predictor = row.names(.), .)) %>% 
#   map(., function(.data){
#     .data %>% 
#       rowwise() %>% 
#      mutate(var_imp = mean(c_across(rand1:rand10)), .keep = "unused") %>% 
#       arrange(var_imp)
#  }) %>% 
#  rbindlist(., idcol = "CV_run") %>% 
#  pivot_wider(., values_from = var_imp, names_from = CV_run, names_prefix = "RUN") %>% 
#  rowwise() %>% 
#  mutate(.,
#         max = max(c_across(starts_with("RUN"))),
#         mean = mean(c_across(starts_with("RUN"))),
#         min = min(c_across(starts_with("RUN")))) %>%
#  ungroup() %>% 
#  mutate(predictor = forcats::as_factor(predictor)) %>% 
#  mutate(predictor = fct_reorder(predictor, mean))

# Read in John's variable importance ensemble summary for comparison
readRDS("Output/VarImp_Ensemble_JohnOrig.rds")

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
saveRDS(table_var_impEM, file = "Output/3_VarImp_Ensemble.rds")





# 4. Project current distribution ====

# a. Load predictors
dem <- raster("Data/Rasters/bathy_m.tif")
substrate <- raster("Data/Rasters/substrate_reclassified.tif")
substrate <- as.factor(substrate)
substrate <- ratify(substrate)
rat <- levels(substrate)[[1]]
rat$substrate <- c("Muddy", "Sand & Mud", "Sandy", "Rocky")
levels(substrate) <- rat
slope <- raster("Data/Rasters/MAR_coastal_slope.tif")
bpi_broad <- raster("Data/Rasters/bpi_broad.tif")
bpi_fine <- raster("Data/Rasters/bpi_fine.tif")
rei <- raster("Data/Rasters/rei_scaled.tif")

# land polygon with 5-km buffer
land_buffer <- st_read("Data/Shapefiles/coast50k_buff_5km.shp")

# Mask for Bay of Fundy
BoF_mask <- st_read("Data/Shapefiles/BoF_mask.shp")

# Create mask for study area
atl_mask <- raster("Data/Rasters/NS_12m_contour.tif") %>% 
  mask(., land_buffer) %>% # exclude cells > 5 km from shore
  mask(., BoF_mask, inverse = TRUE)

# Create fishnet grid from reclassified raster template
predictors <- stack(dem, slope, bpi_broad, bpi_fine, substrate, rei) %>% 
  mask(., mask = atl_mask)

names(predictors) <- c("bathy_m","slope","bpi_broad", # rename layers
                       "bpi_fine","sub_type","rei")

writeRaster(predictors, filename = "Data/Rasters/PredictorStack_mask.tif", overwrite = TRUE)




  
# b. Project current distribution for all individual models and CV runs
myBiomodProjection <- BIOMOD_Projection(myBiomodModelOut,
                                        new.env = stack(predictors),
                                        proj.name = "CurrentZostera2_test",
                                        models.chosen = c("ZosteraMarina_allData_RUN1_RF", "ZosteraMarina_allData_RUN1_ANN"),
                                        metric.binary = NULL,
                                        build.clamping.mask = FALSE, 
                                        compress = TRUE)

get_built_models(myBiomodModelOut) # use this command to get list of fitted model names

?BIOMOD_Projection

# Outputs predictions into CurrentZostera sub-folder, file = proj_CurrentZostera_ZosteraMarina.tif
# This .tif file has 35 bands (one for each model); these are added to QGIS project and through symbology can display one at a time
# However, the individual_projections sub-folder is empty... this is where I expected individual model .tif files to show up

str(myBiomodProjection)
myBiomodProjection





# 5. Project current distribution predicted from ensemble ====
myEnsembleProjection <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                                   bm.proj = myBiomodProjection,
                                                   proj.name = "EM_Curr_Projection",
                                                   models.chosen = c("ZosteraMarina_EMwmeanByROC_allData_RUN1_mergedAlgo", "ZosteraMarina_EMwmeanByROC_allData_RUN2_mergedAlgo"),
                                                   metric.binary = c("ROC","TSS"),
                                                   compress = TRUE)

get_built_models(myBiomodEM)

?BIOMOD_EnsembleForecasting

cat("Proceed to 3FigurePrep.R")