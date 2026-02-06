# ==============================================================================
# SCRIPT: TESTING THE ASSOCIATION BETWEEN ENVIRONMENTAL FACTORS AND THE DISPERSAL
# LOCATIONS OF INFERRED VIRAL LINEAGES 
# ==============================================================================
# Description: Tests the association between environmental variables (landscape 
# factors) and the inferred dispersal locations of viral lineages using the
# seraphim framework. This is a generalized framework applicable to all analysed 
# epidemics. Users can adapt the script for different periods by updating the
# directory paths and environmental rasters.
# ==============================================================================

# 1. SETUP WORKING ENVIRONMENT AND LIBRARIES
# ------------------------------------------------------------------------------
setwd("./")

library(lubridate)
library(seraphim)
library(sf)
library(terra)
library(raster)

localTreesDirectory = "Tree_extractions"
nberOfExtractionFiles = 1000


# 3. ENVIRONMENTAL RASTER PREPARATION
# ------------------------------------------------------------------------------
# Loads .tif files, ensures no negative values exist, and crops/masks them 
# to the study area (Italy) defined by a shapefile.
rast.names <- list.files("./data/env_raster")
rast.names <- rast.names[endsWith(rast.names, ".tif")]
envt.list <- lapply(rast.names, function(x) raster(paste("./raster_epidemic_2016-2017/",x, sep = "")))
names(envt.list) <- gsub(".tif", "", rast.names)

for(i in 1:length(envt.list)) {
  names(envt.list[[i]]) <- names(envt.list)[i]
}

# Cleanup: Avoid negative values and mask to Italy
envt.list <- lapply(envt.list, function(x) {x[x[]<0] = 0; return(x)})
ita.sh <- read_sf("./italy.shp")
ita.sh <- st_transform(ita.sh, crs("+proj=longlat +datum=WGS84 +no_defs"))
envt.list <- lapply(envt.list, function(x) crop(x, ita.sh))
envt.list <- lapply(envt.list, function(x) raster::mask(x, ita.sh))

# Specific fix for farm density raster NAs
envt.list$numero_farms_2017[is.na(envt.list$numero_farms_2017)] <- 0
envt.list$numero_farms_2017 <- raster::mask(envt.list$numero_farms_2017, ita.sh)


# 4. NULL MODEL GENERATION (TREE RANDOMISATION)
# ------------------------------------------------------------------------------
# Generates a null distribution of dispersal by randomising branch locations.
# This serves as the "background" to compare against observed dispersal.
nberOfRandomisations = nberOfExtractionFiles
randomProcedure = 3
nberOfCores = 5
set.seed(123)
treesRandomisation(localTreesDirectory, nberOfRandomisations, envt.list, randomProcedure, nberOfCores)


# 5. SPATIAL DATA EXTRACTION: OBSERVED VS. RANDOM
# ------------------------------------------------------------------------------
# Iterates through every tree (observed and randomised) and extracts the 
# environmental values (e.g., farm density) at the end-point of every branch.
envVariableNames = rep(NA, length(envt.list))
for(i in 1:length(envt.list)) envVariableNames[i] = names(envt.list[[i]])

for(i in 1:nberOfRandomisations) {
  print(paste("Processing extraction:", i))
  obs = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), header=T)
  ran = read.csv(paste0(localTreesDirectory,"/TreeRandomisation_",i,".csv"), header=T)
  envValues_obs = matrix(nrow=dim(obs)[1], ncol=length(envt.list))
  envValues_ran = matrix(nrow=dim(ran)[1], ncol=length(envt.list))
  colnames(envValues_obs) = envVariableNames
  colnames(envValues_ran) = envVariableNames
  
  for(j in 1:length(envt.list)) {
    envValues_obs[,j] = raster::extract(envt.list[[j]], SpatialPoints(obs[,c("endLon","endLat")]))
    envValues_ran[,j] = raster::extract(envt.list[[j]], SpatialPoints(ran[,c("endLon","endLat")]))
  }
  write.csv(envValues_obs, paste0(localTreesDirectory,"/EnvValues_obs_",i,".csv"), row.names=F, quote=F)
  write.csv(envValues_ran, paste0(localTreesDirectory,"/EnvValues_ran_",i,".csv"), row.names=F, quote=F)
}


# 6. STATISTICAL TESTING: BAYES FACTORS (BF)
# ------------------------------------------------------------------------------
# Calculates Bayes Factors to determine if the observed mean environmental values 
# are significantly higher or lower than the null expectation. 
BFs = matrix(nrow=length(envVariableNames), ncol=2)
colnames(BFs) = c("lower","higher")
row.names(BFs) = envVariableNames
meanEnvValues_obs_list = list()
meanEnvValues_ran_list = list()

for(i in 1:length(envVariableNames)) {
  print(paste("Testing variable:", envVariableNames[i]))
  meanEnvValues_obs = rep(NA, nberOfRandomisations)
  meanEnvValues_ran = rep(NA, nberOfRandomisations)
  lowerEnvValues_randomisations = 0
  envVariableName = envVariableNames[i]
  
  for(j in 1:nberOfRandomisations) {
    envValues_obs1 = read.csv(paste0(localTreesDirectory,"/EnvValues_obs_",j,".csv"))[,envVariableName]
    envValues_ran1 = read.csv(paste0(localTreesDirectory,"/EnvValues_ran_",j,".csv"))[,envVariableName]
    
    indices = which(!is.na(envValues_obs1))
    envValues_obs2 = envValues_obs1[indices]
    meanEnvValues_obs[j] = mean(envValues_obs2, na.rm=T)
    
    indices = which(!is.na(envValues_ran1))
    envValues_ran2 = envValues_ran1[indices]
    meanEnvValues_ran[j] = mean(envValues_ran2, na.rm=T)
    
    if(meanEnvValues_obs[j] < meanEnvValues_ran[j]) lowerEnvValues_randomisations = lowerEnvValues_randomisations + 1
  }
  
  # Calculate BF for lower than expected
  p = lowerEnvValues_randomisations/nberOfRandomisations
  BFs[i,"lower"] = round((p/(1-p))/(0.5/(1-0.5)),1)
  
  # Calculate BF for higher than expected
  p = (1-lowerEnvValues_randomisations/nberOfRandomisations)
  BFs[i,"higher"] = round((p/(1-p))/(0.5/(1-0.5)),1)
  
  meanEnvValues_obs_list[[i]] = meanEnvValues_obs
  meanEnvValues_ran_list[[i]] = meanEnvValues_ran
}

write.csv(BFs, paste0("./dispersal_location_H5N8_2016-2017.csv"), quote=F)

