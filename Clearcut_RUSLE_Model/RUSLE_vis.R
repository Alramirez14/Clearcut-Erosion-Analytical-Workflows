###############################################################
######### This script visualizes RUSLE results             ####
#########                                                  ####
######### Created 05/25/2025 by Clearcut Erosion Modeling  ####
###############################################################

library(ggplot2)
library(ggpmisc)
library(mgcv)

rusle_results <- read.csv("rusle_summary_by_year.csv")

rusle_results_pixel <- read.csv("rusle_summary_by_year_pixelwise.csv")

w# Define pretty variable names for facet labels
pretty_names <- c(
  Mean_C = "Cover Factor (C) ",
  Delta_NDVI = "ΔNDVI (year-over-year)",
  Mean_R = "Rainfall Factor (R)"
)

# Prepare long-form data
long_data <- rusle_results %>%
  select(Year, Mean_A, Mean_C, Delta_NDVI, Mean_R) %>%
  pivot_longer(cols = -c(Year, Mean_A), names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = recode(Variable, !!!pretty_names))



# Plot
ggplot(long_data, aes(x = Value, y = Mean_A)) +
  geom_point(size = 3, aes(color = Variable), alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, aes(color = Variable), linewidth = 1.2) +
  facet_wrap(~ Variable, scales = "free", nrow = 1) +
  stat_poly_eq(
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    parse = TRUE,
    label.x.npc = "right", label.y.npc = "top",
    size = 5,
    color = "black",
    inherit.aes = TRUE
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "#E0E0E0", color = NA),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(),
    axis.ticks.x = element_line()
  ) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Drivers of Soil Loss (RUSLE A)",
    x = NULL,
    y = expression(paste("Mean Soil Loss (A) [Tons·acre"^{-1}~"·yr"^{-1}*"]"))
  )

#Plot pixelwise

# Pivot longer for these variables

plot_vars <- c("Delta_NDVI", "LS", "C", "R")

long_pixel_data <- rusle_results_pixel %>%
  select(all_of(c("A", plot_vars))) %>%
  pivot_longer(
    cols = all_of(plot_vars),
    names_to = "Variable",
    values_to = "Value"
  )


# Plot with facets, points, and linear trendline + R^2 equation
ggplot(long_pixel_data, aes(x = Value, y = A)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#1f77b4") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  stat_poly_eq(
    aes(label = ..rr.label..),
    formula = y ~ x, parse = TRUE, size = 3,
    label.x.npc = "right", label.y.npc = 0.1
  ) +
  facet_wrap(~ Variable, scales = "free_x") +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Pixel-level Drivers of Soil Loss (RUSLE A)",
    x = NULL,
    y = expression(paste("Soil Loss (A) [Mg·ha"^{-1}~"·yr"^{-1}*"]"))
  )

###Time series

ggplot(long_data, aes(x = Year, y = Value, group = 1)) +
  geom_line(color = "#2c7fb8", size = 1) +
  geom_point(color = "#2c7fb8", size = 2) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1,
             labeller = labeller(Variable = c(
               Mean_A = "Mean Soil Loss (A)",
               Mean_R = "Mean Rainfall Factor (R)",
               Mean_C = "Mean Cover Factor (C)",
               Delta_NDVI = "Change in NDVI (ΔNDVI)"
             ))) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    panel.grid.minor = element_blank()
  ) +
  labs(
    y = NULL,
    title = "Time Series of RUSLE Components and Vegetation Change"
  )
