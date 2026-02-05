# ==============================================================================
# SCRIPT: BEAST POSTERIOR TREE EXTRACTION
# ==============================================================================
# Description: Extracts 1,000 posterior trees from BEAST .trees files for 
# downstream spatial analysis in seraphim. This script is designed for the
# 2016-2017 epidemic period, characterised by a single monophyletic group analysis. 
# ==============================================================================

setwd("./")

library(seraphim)
library(lubridate)

# 1. POSTERIOR TREES INFORMATION
# ------------------------------------------------------------------------------
# Extracting spatio-temporal information embedded in posterior trees
localTreesDirectory = "Tree_extractions"
allTrees = scan(file="./main_clade_2016-2017.trees", what="", sep="\n", quiet=T)
burnIn = ceiling(length(allTrees) * 0.1)
randomSampling = TRUE
nberOfTreesToSample = 1000
mostRecentSamplingDatum = lubridate::decimal_date(ymd("2017-12-11"))
coordinateAttributeName = "location"

set.seed(123)
treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample,
                mostRecentSamplingDatum, coordinateAttributeName)