# ==============================================================================
# SCRIPT: SPATIAL WAVEFRONT STATISTICS PER MONOPHYLETIC GROUP
# ==============================================================================
# Description: This script calculates spatial wavefront distances for multiple 
# independently analyzed monophyletic groups. Designed for complex epidemic
# periods (e.g., 2021-2022 adn 2022-2023 epidemics) where multiple groups of
# sequences are analysed separately. The script iterates through group-specific
# directories to generate individual wavefront summary statistics.
# ==============================================================================

setwd("./")

library(seraphim)
library(doMC)
library(dplyr)
library(ggplot2)
library(lubridate)
library(reshape)

# Source the custom statistics script
source("./spreadStatistics.R") # Available from: https://github.com/sdellicour/seraphim/blob/master/windows/R/spreadStatistics.r

# ------------------------------------------------------------------------------
# PART 1: SPATIAL WAVEFRONT EXTRACTION (2021-2022)
# ------------------------------------------------------------------------------
# Pre-processing step:
# 1. Use FigTree to identify monophyletic groups.
# 2. Create a CSV/table with columns: "sequences" and "group".
# 3. Map each sequence name to its respective group number (1 to n).
ds.groups <- readxl::read_excel("./groups_2021-2022.xlsx")

get_path_to_root <- function(node, parent_lookup) {
  path <- c(node)
  while(as.character(node) %in% names(parent_lookup)) {
    node <- parent_lookup[[as.character(node)]]
    path <- c(path, node)
  }
  return(path)
}

# Loop for Monophyletic Group Extraction
for(i in 1:length(unique(ds.groups$group))) { 
  print(paste("Extracting group", i))
  temp.cl <- ds.groups[ds.groups$group == i, ]
  
  for(j in 1:1000) { 
    file_path <- paste0("./TreeExtractions_", j, ".csv")
    if(!file.exists(file_path)) next
    
    temp.tre <- read.csv(file_path)
    parent_lookup <- setNames(temp.tre$node1, temp.tre$node2)
    
    my_tips <- temp.cl$sequences
    tip_rows <- temp.tre[temp.tre$tipLabel %in% my_tips, ]
    if(nrow(tip_rows) < 2) next
    
    tip_nodes <- tip_rows$node2
    tip_labels <- tip_rows$tipLabel
    
    # Trace paths and find MRCA
    paths <- lapply(tip_nodes, get_path_to_root, parent_lookup = parent_lookup)
    names(paths) <- tip_labels
    
    find_mrca <- function(p1, p2) {
      intersect_nodes <- intersect(p1, p2)
      if(length(intersect_nodes) == 0) return(NA)
      return(intersect_nodes[1])
    }
    
    tip_combinations <- combn(names(paths), 2, simplify = FALSE)
    mrca_list <- sapply(tip_combinations, function(x) find_mrca(paths[[x[1]]], paths[[x[2]]]))
    mrca_table <- sort(table(mrca_list), decreasing = TRUE)
    if(length(mrca_table) == 0) next 
    
    best_mrca <- names(mrca_table)[1]
    kept_tips <- names(paths)[sapply(paths, function(p) best_mrca %in% p)]
    if(length(kept_tips) < 2) next
    
    # Rebuild and Save Subtree
    tip_nodes_sub <- temp.tre$node2[temp.tre$tipLabel %in% kept_tips]
    nodes_to_check <- tip_nodes_sub
    visited_rows <- data.frame()
    while(length(nodes_to_check) > 0) {
      matching_rows <- temp.tre[temp.tre$node2 %in% nodes_to_check, ]
      visited_rows <- rbind(visited_rows, matching_rows)
      nodes_to_check <- setdiff(matching_rows$node1, visited_rows$node2)
    }
    selected_nodes <- unique(c(visited_rows$node1, visited_rows$node2))
    sub.temp.tre <- temp.tre[temp.tre$node1 %in% selected_nodes & temp.tre$node2 %in% selected_nodes, ]
    
    out_dir <- paste0("./group", i)
    if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    write.csv(sub.temp.tre, paste0(out_dir, "/TreeExtractions_", j, ".csv"), row.names = FALSE)
  }
}

# ------------------------------------------------------------------------------
# PART 2: ESTIMATION OF SPREAD STATISTICS FOR EACH MONOPHYLETIC GROUP
# ------------------------------------------------------------------------------
for(i in 1:length(unique(ds.groups$group))) {
  message(paste("Calculating Spread Stats for group", i))
  localTreesDirectory = paste0("./group", i, "/")
  outputDir = paste0("./s", i, "/")
  if(dir.exists(localTreesDirectory) && length(list.files(localTreesDirectory)) > 0) {
    if(!dir.exists(outputDir)) {
      message(paste("Creating directory:", outputDir))
      dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
    }
    outputName = paste0(outputDir, "s", i, "_summary")
    
    spreadStatistics(localTreesDirectory = localTreesDirectory, 
                     nberOfExtractionFiles = 1000, 
                     timeSlices = 100,
                     onlyTipBranches = FALSE, 
                     showingPlots = FALSE, 
                     outputName = outputName, 
                     nberOfCores = 10, 
                     slidingWindow = 1/24)
    gc()
  } else {
    warning(paste("Directory", localTreesDirectory, "is missing or empty. Skipping group", i))
  }
}

# ------------------------------------------------------------------------------
# PART 3: WAWEFRONT DISTANCE DATA PREPARATION
# ------------------------------------------------------------------------------
hpd.swd <- data.frame(); median.swd <- data.frame(); mean.swd <- data.frame()

for(i in 1:length(unique(ds.groups$group))) {
  path_prefix <- paste0("./s", i, "/s", i, "_summary_")
  
  if(file.exists(paste0(path_prefix, "95%HPD_spatial_wavefront_distance.txt"))) {
    # HPD
    hpd.temp <- read.table(paste0(path_prefix, "95%HPD_spatial_wavefront_distance.txt"), header = T)
    hpd.temp$time <- as.Date(lubridate::date_decimal(hpd.temp$time))
    hpd.temp$X95.HPD_lower_value[is.na(hpd.temp$X95.HPD_lower_value)] <- 0
    hpd.temp$X95.HPD_higher_value[is.na(hpd.temp$X95.HPD_higher_value)] <- 0
    hpd.temp$group <- paste0("s", i)
    hpd.swd <- rbind(hpd.swd, hpd.temp)
    
    # Median
    median.temp <- read.table(paste0(path_prefix, "median_spatial_wavefront_distance.txt"), header = T)
    median.temp$time <- as.Date(lubridate::date_decimal(median.temp$time))
    median.temp$distance[is.na(median.temp$distance)] <- 0
    median.temp$group <- paste0("s", i)
    median.swd <- rbind(median.swd, median.temp)
    
    # Mean
    mean.temp <- read.table(paste0(path_prefix, "mean_spatial_wavefront_distance.txt"), header = T)
    mean.temp$time <- as.Date(lubridate::date_decimal(mean.temp$time))
    mean.temp$distance[is.na(mean.temp$distance)] <- 0
    mean.temp$group <- paste0("s", i)
    mean.swd <- rbind(mean.swd, mean.temp)
  }
}

# Filter for specific groups of interest (i.e., the ones containing >= 10 sequences)
target_groups <- c("s2", "s4", "s6")
hpd.swd <- hpd.swd[hpd.swd$group %in% target_groups, ]
median.swd <- median.swd[median.swd$group %in% target_groups, ]
mean.swd <- mean.swd[mean.swd$group %in% target_groups, ]

# Final save
saveRDS(mean.swd, "./mean.swd.RDS")
saveRDS(median.swd, "./median.swd.RDS")
saveRDS(hpd.swd, "./hpd.swd.RDS")

