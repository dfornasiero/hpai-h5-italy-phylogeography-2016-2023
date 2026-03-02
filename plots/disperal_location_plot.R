# ==============================================================================
# SCRIPT: DISPERSAL LOCATIONS ANALYSIS (BF FACTORS)
# Description: Processes Bayes Factor (BF) results across epidemics to visualize 
# associations between environmental variables and HPAI lineage locations.
# ==============================================================================

require(xlsx)
require(ggplot2)
library(purrr)
library(dplyr)
library(reshape2)

# 1. Configuration
# ------------------------------------------------------------------------------
setwd("./")

epdm <- c("2016-2017", "2021-2022", "2022-2023")
subtype <- c("H5N8", "H5N1", "H5N1")

custom_labels <- c("Anas platyrhynchos" = "italic(Anas)*plain(' ')*italic(platyrhynchos)",
                   "Bubulcus ibis" = "italic(Bubulcus)*plain(' ')*italic(ibis)",
                   "Egretta garzetta" = "italic(Egretta)*plain(' ')*italic(garzetta)",
                   "Phasianus colchicus" = "italic(Phasianus)*plain(' ')*italic(colchicus)",
                   "Larus ridibundus" = "italic(Larus)*plain(' ')*italic(ridibundus)",
                   "Larus melanocephalus" = "italic(Larus)*plain(' ')*italic(melanocephalus)")

# 2. Data Loading and pre-processing
# ------------------------------------------------------------------------------
# Load dispersal location results
ds.list <- list()
for(i in 1:length(epdm)) {
  ds.list[[i]] <- read.csv(paste0("./epidemic_", epdm[i], "/dispersal_location/Results_dispersal_location_", subtype[i], "_", epdm[i], ".csv"))
}
names(ds.list) <- paste("Epidemic", epdm)
ds <- map_df(ds.list, ~as.data.frame(.x), .id = "epidemic")

# Load and format variable labels
var.lab <- read.xlsx("./var_labs.xlsx", sheetIndex = 1)
var.lab$var_type <- factor(var.lab$var_type, levels = c("Landscape", "Poultry-related", "Species abundance"))
var.lab$label1 <- factor(var.lab$label1, levels = c("Anas platyrhynchos", "Bubulcus ibis", "Egretta garzetta", "Phasianus colchicus", "Larus ridibundus", "Larus melanocephalus",
  "Elevation", "CLC 1 — Artificial surfaces", "CLC 2 — Agricultural areas", "CLC 3 — Forest and semi-natural areas", "CLC 4 — Wetlands", "CLC 5 — Water bodies",
  "Wind speed", "Distance from wetland", "Distance from Important Bird Areas", "Start-of-season date", "End-of-season date", "Season length",
  "Human footprint", "Fine particulate matter (PM2.5)", "Landscape fragmentation","no. poultry farms", "Simpson pop. farms per species"))

# Merge and reshape
var.lab$label <- NULL
ds <- merge(ds, var.lab, all.x = T, by = "variable")
ds$variable <- as.factor(ds$variable)
ds$epidemic <- as.factor(ds$epidemic)
ds$variable <- NULL

ds <- reshape2::melt(ds, id.vars = c("epidemic", "label1", "var_type", "year"), variable.name = "BF")

# Filter variables with BF > 20 in at least one epidemic
ds <- ds[ds$value > 3, ]
ds <- ds %>%
  group_by(label1) %>%
  filter(any(value > 20)) %>%
  ungroup()

# 3. Visualization logic
# ------------------------------------------------------------------------------
plt <- ggplot(ds) +
  facet_grid(var_type ~ epidemic, scales = "free_y", space = "free", switch = "y", 
             labeller = as_labeller(c("Epidemic 2016-2017" = "Epidemic 2016–2017",
                                      "Epidemic 2021-2022" = "Epidemic 2021–2022",
                                      "Epidemic 2022-2023" = "Epidemic 2022–2023",
                                      "Poultry-related" = "Poultry-\nrelated",
                                      "Landscape" = "Landscape",
                                      "Species abundance" = "Species\nabundance"))) +
  geom_col(aes(y = label1, x = value, fill = as.factor(BF)), width = 0.7) + 
  scale_y_discrete(limits = rev,
                   labels = function(x) ifelse(x %in% names(custom_labels), parse(text = custom_labels[x]), x)) +
  scale_x_continuous(breaks = c(3, 20, 150), labels = c(3, 20, 150)) +
  coord_cartesian(xlim = c(0, 160)) +
  scale_fill_manual(breaks = c("lower", "higher"),
                    labels = c("Negative association\nwith the dispersal locations\nof inferred HPAI lineages",
                               "Positive association\nwith the dispersal locations\nof inferred HPAI lineages"),
                    values = c("#D6604D", "#4393C3")) +
  labs(x = "Bayes factor", y = NULL) +
  theme_minimal() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_text(angle = 0, size = 13),
        strip.text.x = element_text(size = 15, margin = margin(0, 0, 8, 0)),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "gray80", linetype = "dotted", linewidth = 0.1),
        panel.grid.major.x = element_line(colour = "gray20", linetype = "dashed", linewidth = 0.1),
        panel.spacing = unit(0.2, "lines"),
        axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill = "white", colour = "grey70"))

# 4. Final plot
# ------------------------------------------------------------------------------
windows()
print(plt)
ggsave("./figure4.png", device   = "png", dpi = 600, width = 7679, height = 5211, units = "px")
