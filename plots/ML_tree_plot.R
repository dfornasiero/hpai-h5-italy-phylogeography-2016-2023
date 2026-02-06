# ==============================================================================
# SCRIPT: ML TREE PLOT
# Description: This script loads Maximum Likelihood (ML) trees for multiple 
# epidemic waves, identifies Italian clades, calculates node support
# (SH-aLRT/UFBoot), and generates a multi-panel figure.
# ==============================================================================

library(ggtree)
library(ggplot2)
library(ggpubr)
library(treeio)
library(dplyr)
library(lubridate)
library(sf)
library(phangorn)

# 1. Configuration & Data Loading
# ------------------------------------------------------------------------------
setwd("./")

epdm_labels <- c("Epidemic 2016–2017", "Epidemic 2021–2022", "Epidemic 2022–2023")
seg <- c("2016-2017", "2021-2022", "2022-2023")

# 2. Function to Process and Plot Trees
# ------------------------------------------------------------------------------
plot_ml_tree <- function(folder, title) {
  
  # Construct Path
  tree_path <- paste0("./epi_", folder, "/sequences_", folder, ".fasta.treefile") # ML tree
  
  # Load and Prepare Tree
  raw_tree <- read.tree(tree_path)
  nodes <- raw_tree$tip.label
  
  # Group Italy vs others and midpoint root
  tree_grouped <- groupOTU(raw_tree, nodes[grep("Italy", nodes)], group_name = "origin")
  tree_mid <- midpoint(tree_grouped)
  
  # Fortify data to extract SH-aLRT and UFBoot from labels (format: "SH/UF")
  tree_data <- fortify(tree_mid) %>%
    mutate(
      SHaLRT = as.numeric(sub("/.*", "", label)),
      UFBoot = as.numeric(sub(".*/", "", label)),
      both_supp = ifelse(UFBoot >= 90 & SHaLRT >= 80, TRUE, FALSE))
  
  # Generate Plot
  p <- ggtree(tree_mid, aes(color = origin), size = 0.5, ladderize = FALSE) +
    geom_tippoint(aes(color = origin, size = origin)) +
    geom_rootedge(rootedge = 0.0005, color = 'grey80') +
    scale_size_manual(values = c("0" = 1, "1" = 1.5), guide = "none") +
    scale_color_manual(values = c("0" = "grey80", "1" = "#ae017e"),
                       labels = c("Other EU countries", "Italy")) +
    scale_shape_manual(values = c("TRUE" = 18), 
                       labels = c("Internal node (SH-aLRT ≥ 80%, UFBoot ≥ 90%)")) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    geom_treescale(fontsize = 2.5, offset = 4, y = -25) +
    theme_tree_white() + 
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "grey70"),
          plot.title = element_text(size = 12, hjust = 0.5),
          legend.text = element_text(size = 8)) +
    ggtitle(title)
  
  # Add high-support node markers
  p <- p %<+% tree_data + 
    geom_nodepoint(aes(subset = (both_supp == TRUE), shape = both_supp), 
                   size = 1.5, colour = "grey50")
  
  return(p)
}

# 3. Generate All Plots
# ------------------------------------------------------------------------------
# Iterate through folders and labels to create a list of plots
all_plots <- mapply(plot_ml_tree, seg, epdm_labels, SIMPLIFY = FALSE)

# 4. Arrange and Save
# ------------------------------------------------------------------------------
final_plot <- ggarrange(plotlist = all_plots, 
                        nrow = 1, ncol = 3, 
                        common.legend = TRUE, 
                        legend = "bottom") +
  bgcolor("white") + 
  border("white")

print(final_plot)
ggsave("sm_figure1.png", plot = final_plot, device = "png", dpi = 300, width = 4768, height = 2890, units = "px")
