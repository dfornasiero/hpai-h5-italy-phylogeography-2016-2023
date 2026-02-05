# ==============================================================================
# SCRIPT: DISPERSAL VELOCITY POST-PROCESSING
# ==============================================================================
# Description: Aggregates results from SERAPHIM 'spreadFactors' runs. This script
# parses correlation coefficients, p(Q) values, and Bayes Factors (BF) to generate
# a consolidated results table.The script can be adapted and used for all epidemic
# periods. This is the final step in the environmental analysis pipeline.
# ==============================================================================

setwd("./")

library(dplyr)
library(reshape)

# 1. HELPER FUNCTION: MEAN & 95% CONFIDENCE INTERVAL
# ------------------------------------------------------------------------------
mean_ci <- function(x) {
  n <- length(x)
  mean_x <- round(mean(x, na.rm = TRUE), 4)
  stderr <- sd(x, na.rm = TRUE) / sqrt(n)
  error_margin <- qnorm(0.975) * stderr  # Z-score for 95% confidence level
  
  lower_ci <- round(mean_x - error_margin, 4)
  upper_ci <- round(mean_x + error_margin, 4)
  
  return(c(mean = mean_x, lower_ci = lower_ci, upper_ci = upper_ci))
}

# 2. DEFINE MODELS TO PROCESS
# ------------------------------------------------------------------------------
methods <- c("least-cost", "circuitscape")

for (method in methods) {
  
  # Set dynamic paths based on the method
  method_path <- paste0("./dispersal_velocity/", method, "_path_method/")
  corr_subdir <- paste0(method, "_corr/")
  rand_subdir <- paste0(method, "_rand/")
  
  message(paste("Processing results for:", method))
  
  #-----------------------------------------------------------------------------
  # STEP 1: CORRELATION RESULTS
  #-----------------------------------------------------------------------------
  
  tab.names <- list.files(paste0(method_path, corr_subdir))
  tab.names <- tab.names[grep("linear_regression_results", tab.names)]
  tab.list <- lapply(tab.names, function(x) read.table(paste0(method_path, corr_subdir, x), header = T))
  names(tab.list) <- gsub(".txt", "", tab.names)
  tab.coef <- lapply(tab.list, function(x) lapply(x, function(y) sum(y>0)))
  
  res.list <- list()
  for(i in 1:length(tab.coef)) {
    tryCatch({
      tab <- tab.coef[[i]]
      tab <- as.data.frame((tab))
      tab1 <- reshape::melt(tab)
      if(grepl("_cond_", names(tab.coef[i]))) {
        tab1$variable <- sub("_C$", "", tab1$variable)
        tab1$var_names <- paste(sub(".*(coefficients_|R2_)(.*)", "\\2", tab1$variable), "(C)")
      }
      if(grepl("_res_", names(tab.coef[i]))) {
        tab1$variable <- sub("_R$", "", tab1$variable)
        tab1$var_names <- paste(sub(".*(coefficients_|R2_)(.*)", "\\2", tab1$variable), "(R)")
      }
      
      tab1$sts <- ifelse(grepl("Univariate_LR_coefficients",tab1$variable), "Univariate_LR_coefficients",
                         ifelse(grepl("Univariate_LR_R2",tab1$variable), "Univariate_LR_R2",
                                ifelse(grepl("Univariate_LR_delta_R2",tab1$variable), "Univariate_LR_delta_R2", NA)))
      tab1$variable <- NULL
      tab2 <- reshape(tab1, direction = "wide", 
                      idvar = "var_names",
                      timevar = "sts")
      names(tab2) <- sub("value.", "", names(tab2))
      
      # Statistics for Coefficients
      tab.rc <- tab.list[[i]]
      tab.rc <- tab.rc[ , grep("_LR_coefficients_", colnames(tab.rc))]
      tab.mci.rc <- tab.rc %>% summarise_all(~ mean_ci(.)) %>% t() %>% as.data.frame()
      colnames(tab.mci.rc) <- c("Coeff_mean", "Coeff_lower_95ci", "Coeff_upper_95ci")
      
      if(grepl("_cond_", names(tab.list[i]))) {
        tab.mci.rc$variable <- sub("_C$", "", rownames(tab.mci.rc))
        tab.mci.rc$var_names <- paste(sub(".*(coefficients_|R2_)(.*)", "\\2", tab.mci.rc$variable), "(C)")
      }
      if(grepl("_res_", names(tab.list[i]))) {
        tab.mci.rc$variable <- sub("_R$", "", rownames(tab.mci.rc))
        tab.mci.rc$var_names <- paste(sub(".*(coefficients_|R2_)(.*)", "\\2", tab.mci.rc$variable), "(R)")
      }
      tab.mci.rc$variable <- NULL
      rownames(tab.mci.rc) <- NULL
      tab3 <- merge(tab2, tab.mci.rc, by = "var_names", all = T)
      
      # Statistics for Q (delta R2)
      tab.q <- tab.list[[i]]
      tab.q <- tab.q[ , grep("_delta_R2_", colnames(tab.q))]
      tab.mci.q <- tab.q %>% summarise_all(~ mean_ci(.)) %>% t() %>% as.data.frame()
      colnames(tab.mci.q) <- c("Q_mean", "Q_lower_95ci", "Q_upper_95ci")
      
      if(grepl("_cond_", names(tab.list[i]))) {
        tab.mci.q$variable <- sub("_C$", "", rownames(tab.mci.q))
        tab.mci.q$var_names <- paste(sub(".*(coefficients_|R2_)(.*)", "\\2", tab.mci.q$variable), "(C)")
      }
      if(grepl("_res_", names(tab.list[i]))) {
        tab.mci.q$variable <- sub("_R$", "", rownames(tab.mci.q))
        tab.mci.q$var_names <- paste(sub(".*(coefficients_|R2_)(.*)", "\\2", tab.mci.q$variable), "(R)")
      }
      tab.mci.q$variable <- NULL
      rownames(tab.mci.q) <- NULL
      
      tab4 <- merge(tab3, tab.mci.q, by = "var_names", all = T)
      res.list[[i]] <- tab4
    }, error = function(e) { message(paste("Error processing correlation index", i, ":", e$message)) })
  }
  names(res.list) <- names(tab.coef)
  
  
  # --- VARIABLE SELECTION FOR RANDOMIZATION ---
  # This section identifies variables with p(Q) > 0.90 for follow-up Bayes Factor (BF) estimation.
  # Skip this step if you have already performed the randomization/BF analysis
  # in the 'impact_dispersal_velocity.R' script.

  final.temp <- data.frame()
  for(i in 1:length(res.list)) {
    temp <- res.list[[i]]
    temp$dataset <- gsub(paste0(method, "_|_linear_regression_results"), "", names(res.list)[i])
    final.temp <- rbind(final.temp, temp)
  }
  
  var.sel <- final.temp[final.temp$Univariate_LR_delta_R2 >= 900,]
  var.sel <- var.sel[!is.na(var.sel$var_names),]
  var.sel <- var.sel[c("var_names","dataset")]
  var.sel$var_names <- gsub(" \\(C\\)| \\(R\\)", "", var.sel$var_names)
  saveRDS(var.sel, paste0(method_path, "var.sel.RDS"))
  
  
  #-----------------------------------------------------------------------------
  # STEP 2: RANDOMISATION RESULTS
  #-----------------------------------------------------------------------------
  
  bf.names <- list.files(paste0(method_path, rand_subdir))
  bf.names <- bf.names[grep("_randomisation_", bf.names)]
  
  if(length(bf.names) > 0) {
    bf.list <- lapply(bf.names, function(x) read.table(paste0(method_path, rand_subdir, x), header = T))
    names(bf.list) <- gsub(".txt", "", bf.names)
    
    final.df <- data.frame()
    for(i in 1:length(res.list)) {
      tryCatch({
        tab.res <- res.list[[i]]
        tab.bf <- bf.list[[i]] 
        
        if(grepl("_cond_", names(bf.list[i]))) {
          tab.bf$variable <- sub("_C", "", rownames(tab.bf))
          tab.bf$var_names <- paste(tab.bf$variable, "(C)")
        } 
        if(grepl("_res_", names(bf.list[i]))) {
          tab.bf$variable <- sub("_R$", "", rownames(tab.bf))
          tab.bf$var_names <- paste(tab.bf$variable, "(R)")
        }
        rownames(tab.bf) <- NULL
        tab.final <- merge(tab.res, tab.bf, by = "var_names", all.x = T)
        tab.final$dataset <- gsub(paste0(method, "_|_linear_regression_results"), "", names(res.list)[i])
        final.df <- rbind(final.df, tab.final)
      }, error = function(e) { message(paste("Error merging randomization index", i, ":", e$message)) })
    }
    
    # Save the result for each method
    write.csv(final.df, paste0(method_path, method, "_final_results.csv"), row.names = F)
  }
  
  message(paste("Finished processing", method))
  message("-------------------------------------------")
}