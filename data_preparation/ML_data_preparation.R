# ==============================================================================
# SCRIPT: HPAI SEQUENCE CLEANING AND BEAST DTA FILE PREPARATION
# ==============================================================================
# Description: Standardized pipeline for cleaning sequence data and formatting 
# files for Discrete Trait Analysis (DTA) in BEAST.This script works for all 
# three epidemic periods analysed.
# ==============================================================================


# 1. SETUP WORKING ENVIRONMENT AND LIBRARIES
# ------------------------------------------------------------------------------
# setwd("./")

library(readxl)
library(xlsx)
library(bioseq)
library(ape)
library(dplyr)
library(stringr)
library(tidyverse)
library(phylotools)

# 2. LOAD METADATA
# ------------------------------------------------------------------------------
meta <- read_excel("sequences_europe_metadata.xls", sheet = 1)
meta$Subtype <- gsub("A / ","", meta$Subtype)
meta$Country <-  sub('.*Europe / ', '', meta$Location)
meta$Country <- sub(' /.*', '', meta$Country)
meta$Country <- ifelse(meta$Country == "Isle of Man", "United Kingdom", meta$Country)

# 3. LOAD ALIGNED AND TRIMMED SEQUENCES
# ------------------------------------------------------------------------------
# Loading the aligned fasta files
geno.list <- list()
seg <- c("2016-2017","2021-2022","2022-2023")

for(i in 1:length(seg)) {
  geno.list[[i]] <- read_fasta(paste("ita_outgroup_", seg[[i]], "_al.fasta", sep =""), type = "DNA")
  print(paste("Epidemic", seg[i], "- sequences n°", length(geno.list[[i]])))
}
names(geno.list) <- paste("outgroup_", seg, sep = "")

# 4. LOAD EPIDEMIOLOGICAL DATA
# ------------------------------------------------------------------------------
epi.ds <- readRDS("epi_vir_ds.RDS")
epi.ds$Isolate_Name <- paste(epi.ds$Isolate_Name, epi.ds$collection_date, sep = "/")

# 5. DEFINE OUTLIER SEQUENCES
# This is the list of outliers previously identified in TempEst, per epidemic wave.
# ------------------------------------------------------------------------------
list.outliers <- list( 
  # 2016-2017:
  c("EPI_ISL_239435","EPI_ISL_239436","EPI_ISL_320574","EPI_ISL_240678","EPI_ISL_297234","EPI_ISL_4070128","EPI_ISL_297235"),
  # 2021-2022:
  c("EPI_ISL_7050532","EPI_ISL_10942814","EPI_ISL_13955171","EPI_ISL_11561592","EPI_ISL_5530613","EPI_ISL_5804789",
    "EPI_ISL_8769034","EPI_ISL_10913846","EPI_ISL_6761005","EPI_ISL_18245790","EPI_ISL_11007541","EPI_ISL_11007722",
    "EPI_ISL_11007528","EPI_ISL_11007721","EPI_ISL_17584739","EPI_ISL_17584742","EPI_ISL_17584745","EPI_ISL_14391865",
    "EPI_ISL_18445703","EPI_ISL_9616212","EPI_ISL_10497305","EPI_ISL_17033245","EPI_ISL_11406402","EPI_ISL_13370696",
    "EPI_ISL_13370697","EPI_ISL_13990743","EPI_ISL_14864652","EPI_ISL_8768966","EPI_ISL_11560325","EPI_ISL_13370893",
    "EPI_ISL_18245597","EPI_ISL_8338456","EPI_ISL_17786288","EPI_ISL_10911953","EPI_ISL_11007527"),
  # 2022-2023:
  c("EPI_ISL_18750016","EPI_ISL_18946720","EPI_ISL_16618976","EPI_ISL_16618973","EPI_ISL_16384126","EPI_ISL_18946718",
    "EPI_ISL_16618978","EPI_ISL_15234362","EPI_ISL_16384170","EPI_ISL_15579543","EPI_ISL_17716077","EPI_ISL_15234362",
    "EPI_ISL_16384190","EPI_ISL_18946719","EPI_ISL_16507103","EPI_ISL_18001861","EPI_ISL_15878548","EPI_ISL_16507101",
    "EPI_ISL_15585888","EPI_ISL_15579544","EPI_ISL_16618979","EPI_ISL_15234361","EPI_ISL_15878549","EPI_ISL_15234358",
    "EPI_ISL_16618980","EPI_ISL_16618981","EPI_ISL_16618982","EPI_ISL_16618972","EPI_ISL_16618975","EPI_ISL_16618974")
  )
names(list.outliers) <- names(geno.list)

# 6. CLEANING AND GENETIC CLUSTERING
# ------------------------------------------------------------------------------
list.seq <- list()
for(i in 1:length(geno.list)) {
  set.seed(890)
  otg <- as.data.frame(geno.list[[i]])
  otg$Isolate_Id <- sub("/.*", "", names(geno.list[[i]]))
  otg$Isolate_Name <- gsub("\r","", names(geno.list[[i]]))
  otg$Isolate_Name <- sub("^.*?/(A/.*)", "\\1", otg$Isolate_Name)
  names(otg) <- c("HA","Isolate_Id","Isolate_Name")
  
  # Split into Italian and European subsets
  otg.ita <- otg[grep("Italy|italy|Italia|ITALIA", otg$Isolate_Name), ]
  otg.eu <- otg[!otg$Isolate_Name %in% otg.ita$Isolate_Name, ]
  
  # Process Italian sequence headers for merging with epi data
  otg.ita$RDP <- gsub(".*?/([^/]*VIR[^/]*)/.*", "\\1", otg.ita$Isolate_Name)
  otg.ita$collection_date <- sub(".*([0-9]{4}-[0-9]{2}-[0-9]{2})$", "\\1", otg.ita$Isolate_Name)
  otg.ita$collection_date <- as.Date(otg.ita$collection_date, origin = "1970-01-01")
  otg.ita <- merge(otg.ita, epi.ds[c("RDP","subtype","adis","species_original","collection_date","type")],
                   by = c("RDP","collection_date"), all.x = T)
  
  print(paste("Epidemic", seg[i]))
  
  #--- EUROPEAN SEQUENCES CLEANING (Clustering at 0.5%)
  otg.eu <- merge(otg.eu, meta[c("Isolate_Id","Subtype","Country","Host","Collection_Date")],
                  by = c("Isolate_Id"), all.x = T)
  otg.eu <- otg.eu[!duplicated(otg.eu$Isolate_Name, otg.eu$HA), ]
  otg.eu <- otg.eu[!otg.eu$Isolate_Id %in% list.outliers[[i]], ]
  
  otg.eu <- otg.eu %>%
    mutate(Cluster = bioseq::seq_cluster(bioseq::dna(HA), threshold = 0.005, method = "complete")) %>%
    group_by(Country, Cluster) %>%
    slice_sample(n = 1) %>% ungroup()
  
  #--- ITALIAN SEQUENCES CLEANING (Clustering at 0.5%)
  otg.ita <- otg.ita %>%
    mutate(Cluster = bioseq::seq_cluster(bioseq::dna(HA), threshold = 0.005, method = "complete")) %>%
    group_by(adis, species_original, collection_date, Cluster) %>%
    slice_sample(n = 1) %>% ungroup()
  
  # Reporting statistics
  print(paste("outgroup sequences after =", length(otg.eu$Isolate_Id)))
  print(paste("total final sequences =", length(otg.eu$Isolate_Id) + length(otg.ita$Isolate_Id)))
  print("-----------------------------------")
  
  # Combine subsets
  final.seq <- rbind(otg.ita[c("HA","Isolate_Id","Isolate_Name")],
                     otg.eu[c("HA","Isolate_Id","Isolate_Name")])
  list.seq[[i]] <- final.seq
}
names(list.seq) <- names(geno.list)

# 7. EXPORT CLEANED FASTA FILES
# ------------------------------------------------------------------------------
for(i in 1:length(list.seq)) {
  temp.eu <- list.seq[[i]]
  seqID.eu <- with(temp.eu, paste(Isolate_Id, Isolate_Name, sep="/"))
  
  # Export
  out_fasta <- paste("seq_", names(list.seq[i]), ".fasta", sep = "")
  write.dna(setNames(strsplit(temp.eu$HA, ""), seqID.eu),
            file = out_fasta, format = "fasta", nbcol = -1, colsep = "", append = FALSE)
}

# 8. PREPARE TRAIT FILES FOR BEAST ANALYSIS (DTA)
# ------------------------------------------------------------------------------
# Creating a discrete trait file mapping taxa to Country and State (Italy vs Non-Italy)
trait.list <- list()
for(i in 1:length(list.seq)) {
  temp <- list.seq[[i]]
  temp <- merge(temp, meta[c("Isolate_Id","Country")], by = "Isolate_Id", all.x = T)
  temp$Country <- ifelse(is.na(temp$Country), "Italy", temp$Country)
  trait.list[[i]] <- temp
}
names(trait.list) <- names(list.seq)

# Final formatting and file export
trait.list <- lapply(trait.list, function(x) {x$State <- ifelse(x$Country == "Italy","Italy","Non-Italy"); return(x)})

for(i in 1:length(trait.list)) {
  temp <- trait.list[[i]]
  temp$name <- paste(temp$Isolate_Id, temp$Isolate_Name, sep="/")
  temp <- temp[c("name","Country","State")]
  
  # Exporting .txt trait files
  out_trait <- paste(names(trait.list[i]), ".txt", sep = "")
  caroline::write.delim(temp, out_trait, quote = FALSE, row.names = FALSE, sep = "\t")
}
