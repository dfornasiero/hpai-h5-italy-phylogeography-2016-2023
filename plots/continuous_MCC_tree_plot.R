# ==============================================================================
# SCRIPT: CONTINUOUS PHYLOGEOGRAPHY MCC TREE PLOT
# DESCRIPTION: Visualizing MCC Trees with Geographic and Posterior Annotations
# ==============================================================================

library(treeio)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtree)
library(ggpubr)
library(lubridate)
library(stringr)

# 1. CONFIGURATION & UTILITIES -------------------------------------------------
base_path <- "./"
epi_ds_path <- file.path(base_path, "./epi_vir_ds.RDS")

# Function to convert BEAST numeric heights to Date objects
calc_hpd_dates <- function(range, mrsd) {
  if (all(is.na(range))) return(NA)
  date_max <- mrsd - (range[1] * 365)
  date_min <- mrsd - (range[2] * 365)
  return(c(date_min, date_max))
}

# 2. CORE PLOTTING FUNCTION ----------------------------------------------------
generate_epidemic_tree <- function(tree_path, mrsd, title, main_color, is_complex = FALSE) {
  
  # Load Tree
  tr <- read.beast(tree_path)
  mcc_data <- fortify(tr)
  
  # Extract RDP and Merge Metadata
  mcc_data$RDP <- str_match(mcc_data$label, "/([^/]*VIR[^/]*)/")[ , 2]
  mcc_data <- merge(mcc_data, epi_ds[c("RDP","species_original","type","lon","lat")], all.x = TRUE, by = "RDP")
  
  # Logic for internal node filtering
  mcc_data <- mcc_data %>% 
    mutate(
      node_prob_above_threshold = if(is_complex) {
        (posterior >= 0.9 & studyArea.prob >= 0.9)
      } else {
        posterior >= 0.9
      },
      node_type = ifelse(is.na(label), "Internal node", "Tip node")
    )
  
  # Calculate HPD Dates
  mcc_data$height_0.95_HPD_dates <- lapply(mcc_data$height_0.95_HPD, calc_hpd_dates, mrsd = ymd(mrsd))
  
  # Base Tree Plot
  p <- ggtree(tr, mrsd = mrsd, as.Date = TRUE, size = 0.3, 
              aes(color = if(is_complex) studyArea else NULL)) +
    geom_rootedge(rootedge = 15, color = if(is_complex) "grey80" else main_color) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    ggtitle(title) +
    theme(
      legend.position = "top",
      legend.background = element_rect(fill = "transparent"),
      legend.title = element_blank(),
      axis.ticks.x = element_line(colour = "black", size = 0.5),
      axis.line.x = element_line(colour = "grey70"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = "grey70")
    )
  
  # Add Data Layers
  p <- p %<+% mcc_data +
    geom_tippoint(aes(shape = node_type), size = 1.2, color = main_color) +
    geom_nodepoint(aes(subset = (node_prob_above_threshold == TRUE), shape = node_type), 
                   size = 1.5, colour = "grey30") +
    geom_range(range = "height_0.95_HPD_dates", color = "grey80", alpha = 0.5, size = 2) +
    scale_shape_manual(values = c("Internal node" = 18, "Tip node" = 19),
                       labels = c(paste0("Internal node (", if(is_complex) "PP/ASP" else "PP", " ≥ 90%)"), "Tip node")) +
    guides(shape = guide_legend(nrow = 2), color = guide_legend(nrow = 2))
  
  # Apply studyArea colors if applicable
  if(is_complex) {
    p <- p + scale_color_manual(values = c("Italy" = main_color, "other" = "grey80"),
                                labels = c("Italy" = "Italy", "other" = "Other EU"))
  } else {
    p <- p + scale_color_manual(values = main_color)
  }
  
  return(p)
}

# 3. DATA PROCESSING -----------------------------------------------------------
epi_ds <- readRDS(epi_ds_path)

# Tree 1: 2016-2017
t1_path <- file.path(base_path, "./MCC_main_clade_2016-2017.trees")
fig1 <- generate_epidemic_tree(t1_path, "2017-12-11", "Epidemic 2016–2017", "darkcyan") +
  geom_vline(xintercept = as.Date("2017-01-01"), linetype = "dashed", color = "grey50")

# Tree 2: 2021-2022
t2_path <- file.path(base_path, "./MCC_2021-2022.trees")
fig2 <- generate_epidemic_tree(t2_path, "2022-03-15", "Epidemic 2021–2022", "orange2", is_complex = TRUE) +
  geom_vline(xintercept = as.Date("2022-01-01"), linetype = "dashed", color = "grey50")

# Tree 3: 2022-2023
t3_path <- file.path(base_path, "./MCC_2022-2023.trees")
fig3 <- generate_epidemic_tree(t3_path, "2023-07-29", "Epidemic 2022–2023", "navyblue", is_complex = TRUE) +
  geom_vline(xintercept = as.Date("2023-01-01"), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = as.Date("2022-01-01"), linetype = "dashed", color = "grey50")

# 4. FINAL ASSEMBLY & SAVE -----------------------------------------------------
final_plot <- ggarrange(fig1, fig2, fig3, nrow = 1, common.legend = FALSE)

ggsave("./MCC_continuous_trees.png", plot = final_plot, device = "png", dpi = 300, units = "cm", width = 40, height = 30)

