# ==============================================================================
# SCRIPT: MULTI-PANEL SPATIOTEMPORAL SYNTHESIS (FIGURE 1)
# ==============================================================================
# Epidemic Period: 2016-2017
# Description: Generates the primary composite figure for the 2016-2017 epidemic.
#              This script integrates four output obtained with the previous
#              R/BEAST analyses:
#              A. Epidemic Curves
#              B. MCC Continuous Tree
#              C. Spatial Wavefront Distance
#              D. Continuous Spatial Reconstruction
# ==============================================================================

# 1. SETUP WORKING ENVIRONMENT AND LIBRARIES
# ------------------------------------------------------------------------------
setwd("./")

library(devtools)
library(seraphim)
library(diagram)
library(rgdal)
library(cowplot)
library(treeio)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(sf)
library(ape)
library(glue)
library(readr)
library(ggtree)
library(ggpubr)
library(ggspatial)
library(RColorBrewer)
library(readxl)
library(FluMoDL)
library(aweek)
library(ggspatial)

# 2. DATA EXTRACTION: POSTERIOR TREES AND MCC TREE INFORMATION
# ------------------------------------------------------------------------------
### Continuous reconstruction on map
# 1. Extracting spatio-temporal information embedded in the MCC tree
# For the 2016-2017 epidemic, use mccTreeExtractions.r
source("./mccTreeExtractions.r") # Available from: https://github.com/sdellicour/seraphim/blob/master/windows/R/mccTreeExtractions.R
# For the 2021-2022 and 2022-2023 epidemics, use mcc_tree_extraction.r
source("./mcc_tree_extraction.r") # Available from: https://github.com/dfornasiero/hpai-h5-italy-phylogeography-2016-2023/blob/main/plots/mcc_tree_extraction.r

localTreesDirectory = "Tree_extractions"
nberOfTreesToSample = 1000
mostRecentSamplingDatum = lubridate::decimal_date(ymd("2017-12-11"))

mcc_tre = readAnnotatedNexus("MCC_main_clade_2016-2017.trees")
mcc_data = as.treedata(mcc_tre)
mcc_tab = mccTreeExtractions(mcc_tre, mostRecentSamplingDatum) # 2016-2017 epidemic
mcc_tab = mcc_tree_extraction(mcc_tre, mostRecentSamplingDatum) # 2021-2022 and 2022-2023 epidemics
write.table(mcc_tab, "./epidemic_2016-2017_mcc.csv", head=T,row.names=F, quote=F, sep=",")  


# 2. ESTIMATING HPD REGIONS FOR TIME SLICES
# ------------------------------------------------------------------------------
prob = 0.95
precision = 1/12 #1 month
nberOfExtractionFiles = nberOfTreesToSample
minYears = rep(NA, nberOfExtractionFiles)
for(i in 1:nberOfExtractionFiles) {
  tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
  minYears[i] = min(tab[ ,"startYear"])
}
startDatum = HDInterval::hdi(minYears)[1] # to get the lower bound of the 95% HPD interval for the starting year to consider
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
lapply(polygons, function(x) x@bbox) # check for negative/Inf polygon estimates -> to remove before the next step


# 3. DATA TRANSFORMATION: CUMULATIVE TIME PERIODS FOR MAPPING
# ------------------------------------------------------------------------------
mcc = read.csv("./epidemic_2016-2017_mcc.csv", head=T)  
mcc = mcc[order(mcc[,"endYear"],mcc[,"startYear"],decreasing=F),]
mcc$type_nodes <- ifelse(mcc$node2 %in% mcc$node1, "Internal node", "Tip node") # for internal nodes (dots)
polygons.bind <- lapply(polygons, function(x) {proj4string(x) <- CRS('+proj=longlat +datum=WGS84 +no_defs'); return(x)})
polygons.bind <- do.call(bind, polygons.bind) 
polygons.sf <- st_as_sf(polygons.bind)
polygons.sf$dec_dates <- grep("[0-9]", colnames(polygons.sf), value = TRUE)
polygons.sf$dec_dates <- as.numeric(polygons.sf$dec_dates)
dates.labels <- lubridate::date_decimal(as.numeric(grep("[0-9]", colnames(polygons.sf), value = TRUE)))
dates.labels <- format(dates.labels, '%Y-%b')
dates.labels[duplicated(dates.labels)] <- ""

# Splitting the time periods
# Choose the best time periods to display
dmin <- min(min(mcc$endYear), min(as.numeric(polygons.sf$dec_dates)))
d1 <- decimal_date(as.POSIXct("2016-12-28 23:59:59", tz = "UTC"))
d2 <- decimal_date(as.POSIXct("2017-03-25 23:59:59", tz = "UTC"))
d3 <- decimal_date(as.POSIXct("2017-08-31 23:59:59", tz = "UTC"))
dmax <- max(max(mcc$endYear), max(as.numeric(polygons.sf$dec_dates)))

# Create the labels and cumulative datasets
mcc <- mcc %>%
  mutate(time_range = case_when(
    endYear >= dmin & endYear <= d1 ~ paste(format(date_decimal(dmin), "%b %Y"), "—", format(date_decimal(d1), "%b %Y")),
    endYear > d1 & endYear <= d2 ~ paste(format(date_decimal(d1) %m+% months(1), "%b %Y"), "—", format(date_decimal(d2), "%b %Y")),
    endYear > d2 & endYear <= d3 ~ paste(format(date_decimal(d2) %m+% months(1), "%b %Y"), "—", format(date_decimal(d3), "%b %Y")),
    endYear > d3 & endYear <= dmax ~ paste(format(date_decimal(d3) %m+% months(1), "%b %Y"), "—", format(date_decimal(dmax), "%b %Y")),
    TRUE ~ "Other")) %>%
  mutate(time_range = factor(time_range, levels = c(
    paste(format(date_decimal(dmin), "%b %Y"), "—", format(date_decimal(d1), "%b %Y")),
    paste(format(date_decimal(d1) %m+% months(1), "%b %Y"), "—", format(date_decimal(d2), "%b %Y")),
    paste(format(date_decimal(d2) %m+% months(1), "%b %Y"), "—", format(date_decimal(d3), "%b %Y")),
    paste(format(date_decimal(d3) %m+% months(1), "%b %Y"), "—", format(date_decimal(dmax), "%b %Y")),
    "Other"), 
    ordered = TRUE))

polygons.sf <- polygons.sf %>%
  mutate(time_range = case_when(
    dec_dates <= dmin ~ paste(format(date_decimal(dmin), "%b %Y"), "—", format(date_decimal(d1), "%b %Y")),
    dec_dates > dmin & dec_dates <= d1 ~ paste(format(date_decimal(dmin), "%b %Y"), "—", format(date_decimal(d1), "%b %Y")),
    dec_dates > d1 & dec_dates <= d2 ~ paste(format(date_decimal(d1) %m+% months(1), "%b %Y"), "—", format(date_decimal(d2), "%b %Y")),
    dec_dates > d2 & dec_dates <= d3 ~ paste(format(date_decimal(d2) %m+% months(1), "%b %Y"), "—", format(date_decimal(d3), "%b %Y")),
    dec_dates > d3 & dec_dates <= dmax ~ paste(format(date_decimal(d3) %m+% months(1), "%b %Y"), "—", format(date_decimal(dmax), "%b %Y")),
    dec_dates > dmax ~ paste(format(date_decimal(d3) %m+% months(1), "%b %Y"), "—", format(date_decimal(dmax), "%b %Y")))) %>%
  mutate(time_range = factor(time_range, levels = c(
    paste(format(date_decimal(dmin), "%b %Y"), "—", format(date_decimal(d1), "%b %Y")),
    paste(format(date_decimal(d1) %m+% months(1), "%b %Y"), "—", format(date_decimal(d2), "%b %Y")),
    paste(format(date_decimal(d2) %m+% months(1), "%b %Y"), "—", format(date_decimal(d3), "%b %Y")),
    paste(format(date_decimal(d3) %m+% months(1), "%b %Y"), "—", format(date_decimal(dmax), "%b %Y"))), 
    ordered = TRUE))

mcc.p1 <- mcc[mcc$time_range == levels(mcc$time_range)[1],]; mcc.p1$time_cum <- "t1"
mcc.p2 <- mcc[mcc$time_range %in% levels(mcc$time_range)[1:2], ]; mcc.p2$time_cum <- "t2"
mcc.p3 <- mcc[mcc$time_range %in% levels(mcc$time_range)[1:3], ]; mcc.p3$time_cum <- "t3"
mcc.p4 <- mcc[mcc$time_range %in% levels(mcc$time_range)[1:4], ]; mcc.p4$time_cum <- "t4"
mcc.cum <- rbind(mcc.p1,mcc.p2,mcc.p3,mcc.p4)

poly.p1 <- polygons.sf[polygons.sf$time_range == levels(polygons.sf$time_range)[1],]; poly.p1$time_cum <- "t1"
poly.p2 <- polygons.sf[polygons.sf$time_range %in% levels(polygons.sf$time_range)[1:2], ]; poly.p2$time_cum <- "t2"
poly.p3 <- polygons.sf[polygons.sf$time_range %in% levels(polygons.sf$time_range)[1:3], ]; poly.p3$time_cum <- "t3"
poly.p4 <- polygons.sf[polygons.sf$time_range %in% levels(polygons.sf$time_range)[1:4], ]; poly.p4$time_cum <- "t4"
poly.cum <- rbind(poly.p1,poly.p2,poly.p3,poly.p4)

t_labs <- c("t1" = paste(format(date_decimal(dmin), "%b %Y"), "—", strsplit(as.character(levels(mcc$time_range)[1]), " — ")[[1]][2]),
            "t2" = paste(format(date_decimal(dmin), "%b %Y"), "—", strsplit(as.character(levels(mcc$time_range)[2]), " — ")[[1]][2]),
            "t3" = paste(format(date_decimal(dmin), "%b %Y"), "—", strsplit(as.character(levels(mcc$time_range)[3]), " — ")[[1]][2]),
            "t4" = paste(format(date_decimal(dmin), "%b %Y"), "—", strsplit(as.character(levels(mcc$time_range)[4]), " — ")[[1]][2]))


# 4. PLOTTING: GEOGRAPHIC MAP (MCC TREE AND HPD)
# ------------------------------------------------------------------------------7
# Co-plotting the HPD regions and MCC tree
clean_time_range_levels <- function(time_range_levels) {
  sapply(time_range_levels, function(x) {
    parts <- strsplit(x, " — ")[[1]]
    if (length(parts) == 2 && parts[1] == parts[2]) {
      return(parts[1])
    } else {
      return(x)
    }
  })
}

shp.ita <- read_sf("./italy.shp")
shp.ita <- st_transform(shp.ita, crs = CRS('+proj=longlat +datum=WGS84 +no_defs'))
shp.eu <- read_sf("./europe.shp")
shp.eu <- st_transform(shp.eu, crs = CRS('+proj=longlat +datum=WGS84 +no_defs'))

map.mcc <- ggplot(mcc.cum) +
  facet_wrap(~ time_cum, nrow = 1, labeller = as_labeller(t_labs)) +
  geom_sf(data = shp.ita, fill = "snow", color = "gray90") +
  geom_sf(data = shp.eu[shp.eu$CNTRYNAME != "Italy",], fill = "gray95", color = "gray70") +
  geom_sf(data = shp.eu[shp.eu$CNTRYNAME == "Italy",], fill = "transparent") +
  geom_sf(aes(fill = time_range), alpha = 0.3, data = poly.cum, color = "transparent") +
  geom_point(aes(x = startLon, y = startLat, shape = type_nodes), color = "darkcyan", size = 1.2) + 
  geom_point(aes(x = endLon, y = endLat, shape = type_nodes), color = "darkcyan", size = 1.2) + 
  geom_curve(aes(x = startLon, y = startLat, xend = endLon, yend = endLat), linewidth = 0.2, curvature = -0.2, col = "gray40") +
  theme_map() +
  scale_shape_manual(breaks = c("Internal node", "Tip node"), values = c(21,19))  +
  scale_fill_manual("95% Highest Posterior Density region",
                    breaks = levels(poly.cum$time_range),
                    labels = clean_time_range_levels(levels(poly.cum$time_range)),
                    values = c("#FDE725FF","#67CC5CFF","#2B748EFF","#440154FF")) +
  theme(strip.text.x = element_text(size = 8),
        legend.title = element_text(size=8),
        legend.position = "bottom",
        legend.justification = "center",
        legend.box="horizontal",
        legend.box.just = "center",
        legend.spacing.x = unit(1.5, "cm"),
        legend.key.spacing.y = unit(0, "cm"),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size = 8)) +
  guides(shape = guide_legend(title = NULL, position = "bottom", order = 1, nrow = 2),
         fill = guide_legend(title.position = "top", order = 2, nrow = 1)) +
  coord_sf(crs = CRS('+proj=longlat +datum=WGS84 +no_defs'), xlim = c(6.8, 18), ylim = c(37, 47), expand = T) +
  annotation_scale(location = "bl", width_hint = 0.2, height = unit(0.1, "cm"), bar_cols = c("gray40", "white"),
                   data = subset(mcc.cum, time_cum == "t1")) +
  ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                    style = north_arrow_fancy_orienteering(text_col = "gray40", text_size = 6, fill = c("white", "gray40"),),
                                    height = unit(0.7, "cm"), width = unit(0.7, "cm"),
                                    data = subset(mcc.cum, time_cum == "t4"))


# 5. PLOTTING: SPATIAL WAVEFRONT STATISTICS
# ------------------------------------------------------------------------------
hpd.sw <- read.table("./summary_95%HPD_spatial_wavefront_distance.txt", header = T)
median.sw <- read.table("./summary_median_spatial_wavefront_distance.txt", header = T)
hpd.sw$time <- as.Date(lubridate::date_decimal(hpd.sw$time))
median.sw$time <- as.Date(lubridate::date_decimal(median.sw$time))
median.sw$group <- "H5N8-A/wild duck/Poland/82 A/2016-like"

sw <- ggplot() +
  geom_ribbon(aes(x=time, ymin=X95.HPD_lower_value, ymax = X95.HPD_higher_value), 
              fill = "slategray2", alpha = 0.5, data = hpd.sw) +
  geom_line(aes(x=time, y=distance, linetype = group), 
            data=median.sw, col = "steelblue4", linewidth = 0.5) +
  annotate("segment", x = as.Date("2017-01-01"), xend = as.Date("2017-01-01"), 
           y = 0, yend = Inf, linetype = "dashed", color = "grey30") +
  annotate("label", x = as.Date("2017-01-01"), y = Inf, label = "2016", 
           color = "grey30", size = 3, vjust = 1, hjust = 1, nudge_x = -1) +
  annotate("label", x = as.Date("2017-01-01"), y = Inf, label = "2017", 
           color = "grey30", size = 3, vjust = 1, hjust = 0, nudge_x = 1) +
  scale_x_date(limits = c(ymd("2016-12-20"), ymd("2017-12-15")),
               date_breaks = "1 month", date_labels = "%b") +
  scale_linetype_manual("Monophyletic group", values = "solid") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        legend.position = c(0.80, 0.15),
        legend.background = element_blank()) +
  ylab("Spatial wavefront distance (km)")

# ----------------------------------------------------------------------------
# MULTI-CLADE VISUALIZATION (OPTIONAL)
# In case of multiple monophyletic groups (see: spatial_wavefront_distance_per_monophyletic_group.R),
# replace in the plot the following block to map group-specific linetypes and labels.
# ----------------------------------------------------------------------------
# geom_line(aes(x=time, y=distance, group = group, linetype = group), 
#           data=median.sw, col = "steelblue4", linewidth = 0.5) +
# scale_linetype_manual("Monophyletic group", 
#                       values = c("dotdash","solid","dotted"),
#                       breaks = c("s4","s6","s2"),
#                       labels = c("s2"="(c) EA-2021-AD", "s4"="(a) EA-2020-C", "s6"="(b) EA-2020-C")) +
# guides(linetype=guide_legend(ncol=1, byrow=TRUE)) +


# 6. PLOTTING: CONTINUOUS MCC TREE VISUALIZATION
# ------------------------------------------------------------------------------
epi.ds <- readRDS("./epi_vir_ds.RDS")
mcc.list.1 <- read.beast("./MCC_main_clade_2016-2017.trees")
mcc.data.1 <- fortify(mcc.list.1)
mcc.data.1$RDP <- stringr::str_match(mcc.data.1$label, "/([^/]*VIR[^/]*)/")[ , 2]
mcc.data.1 <- merge(mcc.data.1, epi.ds[c("RDP","species_original","type","lon","lat")], all.x = T, by = "RDP")
mcc.data.1 <- mcc.data.1 %>% mutate(node_prob_above_threshold = case_when(posterior >= 0.9 ~ T, T ~ F))
mcc.data.1$node_type <- ifelse(is.na(mcc.data.1$label), "Internal node", "Tip node")

numeric_date_range_to_date <- function(range, mrsd) {
  if (all(is.na(range))) return(NA)
  date_max <- mrsd - (range[1] * 365)
  date_min <- mrsd - (range[2] * 365)
  return(c(date_min, date_max))
}

mcc.data.1$height_0.95_HPD_dates <- lapply(X = mcc.data.1$height_0.95_HPD, FUN = numeric_date_range_to_date, mrsd = ymd("2017-12-11"))

tree_figure.1 <- ggtree(temp.tre.1, mrsd = "2017-12-11", color = "darkcyan", as.Date = T, size = 0.3) +
  geom_rootedge(rootedge = 20, color="darkcyan") +
  # Static markers for year transition
  annotate("segment", x = as.Date("2017-01-01"), xend = as.Date("2017-01-01"), 
           y = -Inf, yend = Inf, linetype = "dashed", color = "grey30") +
  annotate("label", x = as.Date("2017-01-01"), y = Inf, label = "2016", 
           color = "grey30", size = 3, vjust = 1, hjust = 1) +
  annotate("label", x = as.Date("2017-01-01"), y = Inf, label = "2017", 
           color = "grey30", size = 3, vjust = 1, hjust = 0) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  theme_tree2() %<+% mcc.data.1 + 
  geom_tippoint(aes(shape = node_type), size = 1.2, color = "darkcyan") +
  geom_nodepoint(aes(subset = node_prob_above_threshold == T, shape = node_type), 
                 size = 1.5, colour = "grey30") + 
  scale_shape_manual(values = c("Internal node" = 18, "Tip node" = 19),
                     labels = c("Internal node (PP ≥ 90%)", "Tip node")) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", colour = "grey70"))

# 7. PLOTTING: EPIDEMIC CURVE DATA PREPARATION
# ------------------------------------------------------------------------------
epi.dom <- read_excel("./epi_vir_dom.xlsx", sheet = 1) # Domestic outbreak dataset
epi.wb <- read_excel("./epi_vir_wb.xlsx", sheet = 1) # Wild bird outbreak dataset
epi.all <- rbind(epi.dom, epi.wb)
epi.all$week <- isoweek(epi.all$dt_conf, type = "week")
epi.all$year <- format(as.Date(epi.all$dt_conf, format="%Y-%m-%d"),"%Y")
epi.all <- epi.all[c("epidemic","adis","type","lon","lat","dt_conf","dt_ext","week","year")]
epi.all <- epi.all[!duplicated(epi.all),]
epi.all$num <- 1

geno.all <- readRDS("./epi_vir_ds.RDS")
geno.all$n_seq <- 1
geno.agg <- aggregate(n_seq ~ adis + type, geno.all, sum)

epi.all <- merge(epi.all, geno.agg, all.x = T, by = c("adis", "type")) 
epi.all$n_seq <- ifelse(is.na(epi.all$n_seq),0, epi.all$n_seq)
epi1 <- epi.all[epi.all$epidemic == "Epidemic 2016-2017", ]
epi1$dt_conf <- as.Date(epi1$dt_conf)
epi.agg1 <- epi1 %>%
  mutate(week = floor_date(dt_conf, unit="week", week_start = 1)) %>%
  group_by(week, type) %>%
  summarise(num = sum(num),
            n_seq = sum(n_seq))
all.weeks <- tibble(week = seq(min(epi.agg1$week), max(epi.agg1$week), by = "1 week"))
expanded.weeks <- all.weeks %>%
  crossing(type = c("D", "W"))
expanded.epi.agg1 <- expanded.weeks %>%
  left_join(epi.agg1, by = c("week", "type")) %>%
  mutate(num = replace_na(num, 0),
         n_seq = replace_na(n_seq, 0)) %>%
  arrange(week, type)
epi.dates <- range(epi.agg1$week)
all.weeks <- tibble(week = seq(min(epi.agg1$week), max(epi.agg1$week), by = "1 week"))


# 8. PLOTTING: EPIDEMIC CURVE BARPLOT
# ------------------------------------------------------------------------------
epi.curve <- ggplot(expanded.epi.agg1) +
  geom_col(aes(x = week, y = n_seq, fill = type, group = week), alpha = 0.6, position = position_dodge2(preserve = "single")) +
  geom_col(aes(x = week, y = num, color = type, group = week), fill = "transparent", linewidth = 0.5, position = position_dodge2(preserve = "single")) +
  annotate("segment", x = as.Date("2017-01-01"), xend = as.Date("2017-01-01"), 
           y = 0, yend = Inf, linetype = "dashed", color = "grey30") +
  annotate("label", x = as.Date("2017-01-01"), y = Inf, label = "2016", 
           color = "grey30", size = 3, vjust = 1, hjust = 1) +
  annotate("label", x = as.Date("2017-01-01"), y = Inf, label = "2017", 
           color = "grey30", size = 3, vjust = 1, hjust = 0) +
  scale_x_date(limits = c(ymd("2016-12-20"), ymd("2017-12-14")), date_breaks = "1 month", date_labels = "%b") +
  scale_fill_manual("No. Sequences", breaks = c("D", "W"), values = c("D" = "navyblue", "W" = "violetred1"), labels = c("Poultry farms", "Non-poultry birds")) +
  scale_color_manual("No. Outbreaks", breaks = c("D", "W"), values = c("D" = "navyblue", "W" = "violetred1"), labels = c("Poultry farms", "Non-poultry birds")) +
  theme(axis.title.x = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 11),
        panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(2, "mm"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10), axis.title.y = element_text(margin = margin(0,5,0,0), size = 10),
        legend.box = "horizontal", legend.background = element_rect(fill = "transparent"),
        legend.position = c(0.35,0.85), legend.title = element_text(size=8),
        legend.key.spacing.y = unit(0, "cm"), legend.key = element_blank(),
        legend.key.height = unit(0.4, 'cm'), legend.key.width = unit(0.4, 'cm'),
        legend.text = element_text(size = 8), panel.background = element_rect(fill = "white", colour = "grey70")) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE), color=guide_legend(nrow=2,byrow=TRUE)) +
  ylab("No. outbreaks")


# 9. FINAL ASSEMBLY AND OUTPUT SAVING
# ------------------------------------------------------------------------------
windows(12.62, 10.1)

# Combined plot using plot_grid
plot_grid(plot_grid(plot_grid(epi.curve, sw,
                              ncol=1, align="hv", axis = "lr",label_fontface = "plain",
                              labels = c("A.","C."),label_x = -0.015),
                    tree_figure.1, labels = c("","B."), align="v", axis = "bt", rel_widths = c(1,0.7),label_x = -0.035, label_fontface = "plain"),
          map.mcc, ncol=1, nrow=2, rel_heights = c(1,0.8), align="h", axis = "rl", labels = c("","D."), label_x = 0.005, label_fontface = "plain") +
  bgcolor("white") +
  border("white")

ggsave("./figure1.png", device = "png", dpi = 300, units = "cm", width = 32, height = 25.6)

