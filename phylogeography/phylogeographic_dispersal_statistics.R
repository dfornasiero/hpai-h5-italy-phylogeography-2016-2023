# ==============================================================================
# SCRIPT: PHYLOGEOGRAPHIC DISPERSAL STATISTICS
# ==============================================================================
# Description: This script estimates several dispersal statistics (e.g., velocity, 
# weighted diffusion coefficient, spatial wavefront distance over time, etc.)
# The general dispersal statistics can be inferred for any of the three epidemics.
# However, the spatial wavefront analysis in this version is tailored for a single
# monophyletic group analysis (e.g. 2016-2017 epidemic).
# ==============================================================================

# 1. SETUP WORKING ENVIRONMENT AND LIBRARIES
# ------------------------------------------------------------------------------
setwd("./")

library(seraphim)
library(doMC)
library(dplyr)
library(ggplot2)
library(lubridate)

source("./spreadStatistics.R") # Available from: https://github.com/sdellicour/seraphim/blob/master/windows/R/spreadStatistics.r

# 2. ESTIMATION OF DISPERSAL STATISTICS
# ------------------------------------------------------------------------------
# Configuration for calculating weighted diffusion coefficients and velocity.
localTreesDirectory = "Tree_extractions"
nberOfExtractionFiles = 1000
timeSlices = 100
onlyTipBranches = FALSE
outputName = "summary"
slidingWindow = 1/24 # Two-week sliding window
nberOfCores = 10
showingPlots = TRUE

set.seed(123)
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices,
                 onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
