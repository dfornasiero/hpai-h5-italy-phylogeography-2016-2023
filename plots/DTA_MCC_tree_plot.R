# ==============================================================================
# SCRIPT: DTA TREE PLOT (MCC TREES)
# Description: Processes BEAST MCC trees to visualize the discrete phylogenetic
# reconstruction, incorporating Posterior Probabilities and Ancestral State Probs.
# ==============================================================================

library(ggtree)
library(ggplot2)
library(ggpubr)
library(treeio)
library(dplyr)
library(lubridate)
library(ggrepel)

# 1. Configuration & Global Data
# ------------------------------------------------------------------------------
setwd("./")
epi_ds <- readRDS("./epi_vir_ds.RDS")

# Define specific parameters for each wave
wave_params <- list(
  "2016-2017" = list(mrsd = "2017-12-11", label = "Epidemic 2016–2017", nodes = c(238), x_lim = NULL, year_sep = "2017-01-01", years = c("2016", "2017")),
  "2021-2022" = list(mrsd = "2022-03-15", label = "Epidemic 2021–2022", nodes = c(531, 769, 951, 523), x_lim = "2021-02-01", year_sep = "2022-01-01", years = c("2021", "2022")),
  "2022-2023" = list(mrsd = "2023-07-29", label = "Epidemic 2022–2023", nodes = c(615, 524, 844), x_lim = "2021-11-01", year_sep = "2023-01-01", years = c("2022", "2023"))
)

# 2. Functional Plotting Logic
# ------------------------------------------------------------------------------
plot_dta_tree <- function(wave_id, p) {
  
  # Paths
  base_path <- glue::glue("/dta_{wave_id}")
  tree_path <- file.path(base_path, glue::glue("MCC_{wave_id}.trees")) # MCC tree
  
  # Load Tree
  mcc_tree <- read.beast(tree_path)
  
  # Data processing
  mcc_data <- fortify(mcc_tree) %>%
    mutate(
      RDP = stringr::str_match(label, "/([^/]*VIR[^/]*)/")[, 2],
      node_type = ifelse(is.na(label), "Internal node", "Tip node"),
      # Threshold logic: PP >= 0.9 and ASP >= 0.9
      node_prob_above_threshold = ifelse(posterior >= 0.9 & location.prob >= 0.9, TRUE, FALSE),
      node_pp_asp = paste0("PP: ", round(as.numeric(posterior), 2), "\n", 
                           "ASP: ", round(as.numeric(location.prob), 2))
    ) %>%
    left_join(epi_ds[c("RDP", "type")], by = "RDP")
  
  # Base Tree Plot
  plt <- ggtree(mcc_tree, mrsd = p$mrsd, as.Date = TRUE, aes(color = location), size = 0.3) +
    geom_rootedge(rootedge = 15, color = 'grey80') +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = if(!is.null(p$x_lim)) as.Date(c(p$x_lim, NA)) else NULL) +
    scale_color_manual(values = c("Italy" = "#ce1256", "other" = "grey80"),
                       breaks = c("Italy", "other"),
                       labels = c("Italy", "Other EU countries")) +
    scale_shape_manual(values = c("Internal node" = 18, "Tip node" = 19),
                       labels = c("Internal node (PP ≥ 90% and ASP ≥ 90%)", "Tip node")) +
    scale_size_manual(values = c("D" = 1.2, "W" = 1.2, "Italy" = 1.2), guide = "none") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", colour = "grey70")) +
    ggtitle(p$label)
  
  # Add layers with data
  plt <- plt %<+% mcc_data +
    geom_tippoint(aes(shape = node_type, size = type), color = "#ce1256") +
    geom_nodepoint(aes(subset = (node_prob_above_threshold == TRUE), shape = node_type), size = 1.5, colour = "grey50") +
    geom_label_repel(aes(label = ifelse(node %in% p$nodes, node_pp_asp, NA)), 
                     size = 2.5, alpha = 0.7, nudge_x = -100, direction = "y") +
    # Annual separators
    geom_vline(xintercept = as.Date(p$year_sep), linetype = "dashed", color = "grey30") +
    annotate("label", x = as.Date(p$year_sep), y = Inf, label = p$years[1], hjust = 1.2, vjust = 1.5, size = 3) +
    annotate("label", x = as.Date(p$year_sep), y = Inf, label = p$years[2], hjust = -0.2, vjust = 1.5, size = 3)
  
  return(plt)
}

# 3. Final plot
# ------------------------------------------------------------------------------
plot_list <- mapply(plot_dta_tree, names(wave_params), wave_params, SIMPLIFY = FALSE)

combined_fig <- ggarrange(plotlist = plot_list, nrow = 1, common.legend = TRUE, legend = "bottom") +
  bgcolor("white") + border("white")

ggsave("./mcc_tree.png", plot = combined_fig, device = "png", dpi = 300, width = 4768, height = 2890, units = "px")
