# ==============================================================================
# SCRIPT: WEIGTHED DIFFUSION COEFFIENT PLOT
# ==============================================================================
# Description: This script processes post-BEAST diffusion statistics for three
# HPAI epidemic waves in Italy (2016-2023). It extracts median coefficients and
# 95% HPD intervals, converts annual rates to weekly values, and generates a
# multi-panel plot.
# ==============================================================================

library(lubridate)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

# 1. Configuration: Define base paths and wave-specific parameters
# ------------------------------------------------------------------------------
base_path <- "./"
waves <- list(
  "Epidemic 2016—2017" = list(folder = "epidemic_2016-2017", start = "2016-12-20", end = "2017-12-15"),
  "Epidemic 2021—2022" = list(folder = "epidemic_2021-2022", start = "2021-10-05", end = "2022-04-01"),
  "Epidemic 2022—2023" = list(folder = "epidemic_2022-2023", start = "2022-09-01", end = "2023-08-31")
)

# 2. Functional Processing Suite
# ------------------------------------------------------------------------------
# Function to load, merge, and clean diffusion data for a specific wave
process_diffusion <- function(label, params) {
  path_root <- file.path(base_path, params$folder, "./diffusion_statistics/")
  
  # Load raw files
  df_med <- read.table(file.path(path_root, "summary_median_weighted_diffusion_coefficient.txt"), header = TRUE)
  df_hpd <- read.table(file.path(path_root, "summary_95%HPD_weighted_diffusion_coefficient.txt"), header = TRUE)
  
  # Processing pipeline: Merge -> Date Convert -> Normalize to Week -> Filter
  df <- inner_join(df_med, df_hpd, by = "time") %>%
    mutate(
      time = as.Date(date_decimal(time)),
      difc_week = diffusion_coefficient / 52,
      HPD_low_week = X95.HPD_lower_value / 52,
      HPD_hig_week = X95.HPD_higher_value / 52,
      epidemic = label
    ) %>%
    filter(time > ymd(params$start) & time < ymd(params$end)) %>%
    replace(is.na(.), 0) # Replace NAs with 0
  
  return(df)
}

# 3. Data Integration
# ------------------------------------------------------------------------------
difc_all <- map2_dfr(names(waves), waves, process_diffusion)

# 4. Comparative Visualization
# ------------------------------------------------------------------------------
# Define specific year markers for the dashed vertical lines in each facet
year_markers <- data.frame(
  epidemic = names(waves),
  year_date = as.Date(c("2017-01-01", "2022-01-01", "2023-01-01")),
  label = c("2017", "2022", "2023")
)

plot_diffusion <- ggplot(difc_all, aes(x = time)) +
  facet_wrap(~ epidemic, nrow = 1, scales = "free_x") +
  # Draw 95% HPD Area
  geom_ribbon(aes(ymin = HPD_low_week, ymax = HPD_hig_week), alpha = 0.3, fill = "slategray2") +
  # Draw Median Line
  geom_line(aes(y = difc_week), linewidth = 0.6, color = "steelblue4") +
  # Add Year Indicators
  geom_vline(data = year_markers, aes(xintercept = year_date), linetype = "dashed", color = "grey50") +
  geom_label(data = year_markers, aes(x = year_date, y = Inf, label = label), 
             color = "grey30", size = 3, vjust = 1, nudge_x = 1) +
  # Axis and Labels
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(breaks = seq(0, 6000, 500)) +
  labs(y = expression("Weighted Diffusion Coefficient (km"^2*"/week)"), x = NULL) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "plain"),
    strip.background = element_rect(fill = "white", colour = "grey70"),
    panel.background = element_rect(fill = "white", colour = "grey70"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "none"
  )

# ggsave("all_diffusion_coef.png", plot = plot_diffusion, device = "png", dpi = 300, width = 12, height = 7)

