# ==============================================================================
# SCRIPT: VIRAL LINEAGE DISPERSAL VELOCITY
# ==============================================================================
# Description: Tests the impact of environmental factors on diffusion velocity.
# The script can be used for (1) Initial Correlation screening or (2) BF estimation.
# Reference: https://github.com/sdellicour/seraphim/blob/master/tutorials/impact_on_diffusion_velocity.pdf
# ==============================================================================

# 1. SETUP WORKING ENVIRONMENT AND LIBRARIES
# ------------------------------------------------------------------------------
# setwd("./")

library(seraphim)
library(raster)
library(sf)

localTreesDirectory = "Tree_extractions"


# 2. ENVIRONMENTAL RASTER PREPARATION & MASKING
# ------------------------------------------------------------------------------
# Load environmental .tif files, ensure non-negative values, and mask to Italy.
rast.names <- list.files("./data/env_raster")
rast.names <- rast.names[endsWith(rast.names, ".tif")]

envt.list <- lapply(rast.names, function(x) raster(paste("./raster_epidemic_2016-2017/", x, sep = "")))
names(envt.list) <- gsub(".tif", "", rast.names)

for(i in 1:length(envt.list)) {
  names(envt.list[[i]]) <- names(envt.list)[i]
}

# Clean data and apply spatial constraints
envt.list <- lapply(envt.list, function(x) {x[x[] < 0] = 0; return(x)}) 
ita.sh <- read_sf("./data/shp/ita_country/gadm41_ITA_0.shp")
ita.sh <- st_transform(ita.sh, crs("+proj=longlat +datum=WGS84 +no_defs"))
envt.list <- lapply(envt.list, function(x) raster::crop(x, ita.sh))
envt.list <- lapply(envt.list, function(x) raster::mask(x, ita.sh))

# Handle specific farm density NAs
envt.list$numero_farms_2017[is.na(envt.list$numero_farms_2017)] <- 0
envt.list$numero_farms_2017 <- raster::mask(envt.list$numero_farms_2017, ita.sh)


# 3. RASTER TRANSFORMATION (SCALING FACTORS k)
# ------------------------------------------------------------------------------
# Create three sets of rasters transformed by k=10, 100, and 1000.
# Transformation formula: 1 + k * (v / vmax)
envt.list.k10 <- envt.list.k100 <- envt.list.k1000 <- list()

# Helper loop for transformations
for(i in 1:length(envt.list)) {
  temp_raster <- envt.list[[i]]
  vmax <- max(values(temp_raster), na.rm = TRUE)
  
  envt.list.k10[[i]]   <- overlay(temp_raster, fun=function(vo) {return(1 + 10 * (vo / vmax))})
  envt.list.k100[[i]]  <- overlay(temp_raster, fun=function(vo) {return(1 + 100 * (vo / vmax))})
  envt.list.k1000[[i]] <- overlay(temp_raster, fun=function(vo) {return(1 + 1000 * (vo / vmax))})
}

# Apply names to the new lists
names(envt.list.k10) <- names(envt.list.k100) <- names(envt.list.k1000) <- names(envt.list)


# 4. VARIABLE SELECTION BASED ON PREVIOUS CORRELATION
# ------------------------------------------------------------------------------
# This block requires 'var.sel.RDS', which contains variables with p(Q) > 0.90.
# Requirement: Run STEP 1 of 'impact_dispersal_velocity_processing.R' first.
# Note: If you are currently running the initial correlation analysis, skip this part.

var.sel <- readRDS("./var.sel.RDS")

envt.list.k10.cond <- envt.list.k10[var.sel[var.sel$dataset == "k10_cond", "var_names"]]
envt.list.k10.res  <- envt.list.k10[var.sel[var.sel$dataset == "k10_res", "var_names"]]

envt.list.k100.cond <- envt.list.k100[var.sel[var.sel$dataset == "k100_cond", "var_names"]]
envt.list.k100.res  <- envt.list.k100[var.sel[var.sel$dataset == "k100_res", "var_names"]]

envt.list.k1000.cond <- envt.list.k1000[var.sel[var.sel$dataset == "k1000_cond", "var_names"]]
envt.list.k1000.res  <- envt.list.k1000[var.sel[var.sel$dataset == "k1000_res", "var_names"]]


# 5. SERAPHIM ANALYSIS
# ------------------------------------------------------------------------------
# Runs the analyses for each scaling factor (k) and each mode (Resistance vs Conductance).
# We iterate through different movement models and randomization settings:
# Path Models: 2 = Least-cost path, 3 = Circuitscape
# Randomizations: 0 = Correlation only, 1 = Bayes Factor estimation

# Define the parameter combinations
pathModels <- c(2, 3) # 2 = Least-cost, 3 = Circuitscape
randomSettings <- c(0, 1) # 0 = Correlation, 1 = Bayes Factor

# Fixed parameters
nberOfExtractionFiles <- 1000
randomProcedure <- 3
fourCells <- FALSE
showingPlots <- FALSE

# Define datasets to iterate over
datasets <- list(
  list(env = envt.list.k10.res, res = TRUE,  suffix = "k10_res"),
  list(env = envt.list.k10.cond, res = FALSE, suffix = "k10_cond"),
  list(env = envt.list.k100.res, res = TRUE,  suffix = "k100_res"),
  list(env = envt.list.k100.cond, res = FALSE, suffix = "k100_cond"),
  list(env = envt.list.k1000.res, res = TRUE,  suffix = "k1000_res"),
  list(env = envt.list.k1000.cond, res = FALSE, suffix = "k1000_cond"))

for (pm in pathModels) {
  for (rnd in randomSettings) {
    
    model_lab <- ifelse(pm == 2, "least-cost", "circuitscape")
    type_lab  <- ifelse(rnd == 0, "corr", "BF")
    
    for (ds in datasets) {
      current_outputName <- paste(model_lab, ds$suffix, type_lab, sep = "_")
      
      message(paste(">>> Processing:", current_outputName))
      
      tryCatch({
        spreadFactors(
          localTreesDirectory = localTreesDirectory, 
          nberOfExtractionFiles = nberOfExtractionFiles, 
          envVariables = ds$env, 
          pathModel = pm, 
          resistances = rep(ds$res, length(ds$env)), 
          avgResistances = rep(ds$res, length(ds$env)),
          fourCells = fourCells, 
          nberOfRandomisations = rnd, 
          randomProcedure = randomProcedure, 
          outputName = current_outputName, 
          showingPlots = showingPlots
        )
      }, error = function(e) {
        message(paste("ERROR in", current_outputName, ":", e$message))
      })
      
      gc()
    }
  }
}

#-------------------------------------------------------------------------------------------------------------------------------------------
### NOTES
## Scale of interpretation defined by Jeffreys
# BF values     log10(BF)     Strength of evidence 
# 3.16 - 10     0.5 - 1       substantial
# 10 - 31.62    1 - 1.5       strong 
# 31.62 - 100   1.5 - 2       very strong 
# >100          >2            decisive

## Scale of Kass & Raftery
# BF values       Strength of evidence
# 3 - 20          positive
# 20 - 150        strong
# >150            very strong

