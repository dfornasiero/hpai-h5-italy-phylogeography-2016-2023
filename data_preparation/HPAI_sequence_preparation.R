# ==============================================================================
# SCRIPT: HA SEQUENCE DATA PREPARATION AND CLEANING
# ==============================================================================
# Description: Initial data cleaning and preparation pipeline for HA sequences.
# This script works for all three epidemic periods analysed.
# ==============================================================================

# 1. SETUP WORKING ENVIRONMENT AND LIBRARIES
# ------------------------------------------------------------------------------
# Set working directory to the project root
# setwd("./")

library(readxl)
library(xlsx)
library(bioseq)
library(ape)
library(dplyr)
library(stringr)
library(tidyverse)
library(phylotools)

# 2. METADATA LOADING AND PRE-PROCESSING
# ------------------------------------------------------------------------------
# Loading European metadata from GISAID and standardizing subtype naming
meta <- read_excel("sequences_europe_metadata.xls", sheet = 1)
meta$Subtype <- gsub("A / ","", meta$Subtype)

# Loading Italian epidemiological and viral data
ita.seq <- readRDS("epi_vir_ds.RDS")
ita.seq <- merge(ita.seq, meta[c("Isolate_Name","Isolate_Id")], by = "Isolate_Name", all.x = T)
ita.seq <- ita.seq[order(ita.seq$n),]
ita.seq <- ita.seq[order(ita.seq$epidemic),]

# 3. LOADING EUROPEAN BACKGROUND FASTA SEQUENCES
# ------------------------------------------------------------------------------
# Loading European background sequences for each of the three epidemic periods
gisaid.list <- list()
seg <- c("2016-2017","2021-2022","2022-2023")

for(i in 1:length(seg)) {
  gisaid.list[[i]] <- read_fasta(paste("sequences_", seg[[i]], "_europe.fasta", sep =""), type = "DNA")
}
names(gisaid.list) <- paste("outgroup_", seg, sep = "")

# 4. SEQUENCE EXPLORATION AND TIBBLE TRANSFORMATION
# ------------------------------------------------------------------------------
# Inspecting sequence lengths and splitting sequences into separate nucleotide positions
lapply(gisaid.list, function(x) range(seq_nchar(x, gaps = T)))
seqSplit <- lapply(gisaid.list, function(x) stringr::str_split(x, pattern = "", simplify = FALSE))
seqTibble <- lapply(gisaid.list, as_tibble.bioseq_dna)

nucl <- 1:1791 # Number of positions according to Aliview
nucl <- as.character(nucl)
seqTibble <- lapply(seqTibble, function(x) separate(x, sequence, nucl, sep = "")) 

# Basic structural checks
str(seqTibble[[1]][1:50])
rowSums(seqTibble[[1]][c(1:1791)] != c('A', 'T', 'C', 'G', '-'), na.rm=T)

# 5. OUTGROUP DATA CLEANING AND FILTERING
# ------------------------------------------------------------------------------
# Filtering outgroups: removing Italian sequences, duplicates, and short sequences (<1700bp)
list.otg.c <- list()
for(i in 1:length(gisaid.list)) {
  set.seed(890)
  otg <- as.data.frame(gisaid.list[[i]])
  otg$Isolate_Id <- gsub("\r","", names(gisaid.list[[i]]))
  names(otg) <- c("HA","Isolate_Id")
  otg <- otg[c(2,1)]
  rownames(otg) <- NULL
  
  # Merge with metadata for filtering
  otg <- merge(otg, meta[c("Isolate_Id","Isolate_Name","Subtype","Host","Collection_Date")], by = "Isolate_Id", all.x = T)
  
  # Filter: Remove Italian sequences found in the outgroup set
  otg <- otg[-grep("Italy|italy|Italia|ITALIA", otg$Isolate_Name), ]
  
  # Filter: Remove exact duplicates
  otg <- otg[!duplicated(otg$Isolate_Name, otg$HA), ]
  
  # Filter: Keep high-quality sequences only
  otg$HA_length <- str_count(otg$HA, "A|T|C|G")
  otg <- otg[otg$HA_length >= 1700, ]
  
  # Standardize isolate names for file compatibility
  otg$Isolate_Name <- gsub(" ", "_", otg$Isolate_Name)
  otg$Isolate_Name <- gsub("[()]", "", otg$Isolate_Name)
  
  list.otg.c[[i]] <- otg
}
names(list.otg.c) <- names(gisaid.list)

# Verify cleaning results
lapply(list.otg.c, function(x) x[duplicated(x$Isolate_Id),])
lapply(list.otg.c, function(x) length(unique(x$Isolate_Id)))

# 6. EXPORTING FINAL FASTA FILES (ITALY + EUROPEAN OUTGROUPS)
# ------------------------------------------------------------------------------
# Creating combined FASTA files for each epidemic period
for(i in 1:length(list.otg.c)) {
  
  # Process and write Italian sequences
  temp.ita <- ita.seq[ita.seq$epidemic == unique(ita.seq$epidemic)[[i]], ]
  seqID.ita <- with(temp.ita, paste(Isolate_Id,
                                    "A",
                                    species_original,
                                    "Italy",
                                    RDP, # IZSVe Virology lab code  
                                    format(collection_date, "%Y"),
                                    collection_date,
                                    sep="/")) %>% gsub(" ","_",.)
  
  # Output file will be saved in the same directory (e.g., ita_outgroup_2016-2017.fasta)
  out_filename <- paste("ita_", names(list.otg.c[i]), ".fasta", sep = "")
  
  write.dna(setNames(strsplit(temp.ita$HA, ""), seqID.ita),
            file = out_filename,
            format = "fasta", nbcol = -1, colsep = "")
  
  # Append European outgroup sequences to the same file
  temp.eu <- list.otg.c[[i]]
  seqID.eu <- with(temp.eu, paste(Isolate_Id,
                                  Isolate_Name,
                                  Collection_Date,
                                  sep="/"))
  seqID.eu <- gsub(" ","_", seqID.eu)
  
  write.dna(setNames(strsplit(temp.eu$HA, ""), seqID.eu),
            file = out_filename,
            format = "fasta", nbcol = -1, colsep = "", append = T)
}
