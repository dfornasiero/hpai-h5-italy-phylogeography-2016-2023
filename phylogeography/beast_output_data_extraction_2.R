# ==============================================================================
# SCRIPT: BEAST POSTERIOR TREE EXTRACTION AND SPATIAL FILTERING
# ==============================================================================
# Description: Extracts 1,000 posterior trees from BEAST .trees files for 
# downstream spatial analysis in seraphim. This script is designed for the
# 2021-2022 and 2022-2023 epidemics, characterised by multiple monophyletic groups.
# ==============================================================================

setwd("./")

library(lubridate)
library(seraphim)

source("./post_tree_extraction.r")

# 1. POSTERIOR TREES INFORMATION
# ------------------------------------------------------------------------------
# Extracting spatio-temporal information embedded in posterior trees
allTrees = readAnnotatedNexus("./RRW_2021-2022.trees")
mostRecentSamplingDatum = decimal_date(ymd("2022-03-15"))
nberOfTreesToSample = 1000

burnIn = ceiling(length(allTrees) * 0.1)
indices = (burnIn+1):length(allTrees)
set.seed(123)
smpl <- sample(indices, nberOfTreesToSample, replace=F)  # sample 1000 random trees from the beast output
smpl <- sort(smpl)
sampledTrees = allTrees[smpl]

# --- EXTRACTION LOOP ---
for(i in 1:length(sampledTrees)) {
    print(i)
    
    # Extracting tree information
    tab = post_tree_extraction(post_tre=sampledTrees[[i]], mostRecentSamplingDatum)
    write.csv(tab, paste0("./RRW_2021-2022_ext1/TreeExtractions_",i,".csv"), row.names=F, quote=F)

    # Filtering for study area (Italy)
    tab2 = tab
    tab2 = tab2[which((tab2[,"startStudyArea"]=="Italy")&(tab2[,"endStudyArea"]=="Italy")),]
    write.csv(tab2, paste0("./RRW_2021-2022_ext2/TreeExtractions_",i,".csv"), row.names=F, quote=F)
}
