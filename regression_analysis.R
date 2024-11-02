# Author: Nathan Malamud
# Date: 2024.10.25
#
# Goal: use regression analysis to identify signficant correlations
# across spectrophotometric, structural, and physiological leaf measurements.
#

# Analysis Packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(GGally)
library(dplyr)
library(readxl)
library(reshape2)
library(scales)
library(LeafArea)
library(stringr)
library(RColorBrewer)

# Load data for traits
# REMINDER: Set Working Directory -> Source File Location
data <- read_csv("./data/traits.csv")

# TODO: move prospect_rtm_output to data
# Load prospect measurements from spec curves
prospect <- read_csv("./data/molecular_content.csv") %>%
  select(-c(sampleID, species, treatment_mmol))
data <- merge(data, prospect, by = "barcodeID")

# Define Factor Levels (Treatment and Species)
data$treatment_mmol <- data$treatment_mmol
data$species <- as.factor(data$species)
levels(data$species) <- c("R. sativus", "B. officinalis", "H. vulgare")
data$LDMC <- as.numeric(data$LDMC)

# Create a new column for aggregated treatment levels
data <- data %>%
  mutate(treatment_level = case_when(
    treatment_mmol <= 5 ~ "0 - 5 mmol",
    treatment_mmol <= 15 ~ "10 - 15 mmol",
    treatment_mmol <= 25 ~ "20 - 25 mmol",
    treatment_mmol <= 35 ~ "30 - 35 mmol",
    TRUE ~ "Other"  # Optional: catch any values above 35
  ))

# Define a custom theme for all plots
custom_theme <- theme_minimal(base_family = "sans") + 
  theme(
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.25, color = "grey80"),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5, color = "black"),
    strip.text = element_text(hjust = 0, size = 10, face = "italic"),
    legend.position = "none"
  )

# Define colors for each species
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c", "grey" = "#d3d3d3")

# Function to calculate significant regressions with debug print statements
calculate_significant_regressions <- function(data, x, y) {
  data %>%
    group_by(species, treatment_level) %>%
    summarize(
      fit = list(lm(reformulate(x, y), data = cur_data())),
      r_squared = summary(fit[[1]])$r.squared,
      p_value = summary(fit[[1]])$coefficients[2, 4]
    ) %>%
    filter(p_value < 0.05) %>%  # Only retain significant results
    ungroup() %>%  # Ungroup to avoid nested list issues
    mutate(variable_pair = paste(x, y, sep = "_"))  # Track variable pairs for clarity
}

# Centralized plot function with added print statements to inspect data structure
plot_correlation <- function(data, x, y, x_label, y_label) {
  print(paste("Plotting", x, "vs", y))
  print("Data preview:")
  print(head(data))  # Inspect data structure
  
  ggplot(data, aes_string(x = x, y = y)) +
    geom_point(aes(color = color), size = 1.5) +
    scale_color_identity() +
    geom_smooth(data = subset(data, !is.na(p_value) & species == "R. sativus"),
                aes(color = josef_colors["R. sativus"]), method = "lm", se = TRUE, size = 0.8) +
    geom_smooth(data = subset(data, !is.na(p_value) & species == "B. officinalis"),
                aes(color = josef_colors["B. officinalis"]), method = "lm", se = TRUE, size = 0.8) +
    geom_smooth(data = subset(data, !is.na(p_value) & species == "H. vulgare"),
                aes(color = josef_colors["H. vulgare"]), method = "lm", se = TRUE, size = 0.8) +
    labs(x = x_label, y = y_label) +
    stat_cor(
      data = subset(data, !is.na(p_value) & species == "R. sativus"),
      aes(label = paste(..rr.label.., ..p.label.., sep = " ~ "), color = josef_colors["R. sativus"]),
      label.x.npc = 0.45, label.y.npc = 0.85, size = 2.5, method = "pearson"
    ) +
    stat_cor(
      data = subset(data, !is.na(p_value) & species == "B. officinalis"),
      aes(label = paste(..rr.label.., ..p.label.., sep = " ~ "), color = josef_colors["B. officinalis"]),
      label.x.npc = 0.45, label.y.npc = 0.75, size = 2.5, method = "pearson"
    ) +
    stat_cor(
      data = subset(data, !is.na(p_value) & species == "H. vulgare"),
      aes(label = paste(..rr.label.., ..p.label.., sep = " ~ "), color = josef_colors["H. vulgare"]),
      label.x.npc = 0.45, label.y.npc = 0.65, size = 2.5, method = "pearson"
    ) +
    custom_theme
}

# Generate plots with corrected color assignment logic for significant regressions
generate_plots <- function(relationships, data) {
  plots <- list()
  for (rel in relationships) {
    # Calculate significant regressions for the current relationship
    p_values <- calculate_significant_regressions(data, rel$x, rel$y)
    
    # Merge significant regressions back into the main data
    data_rel <- data %>%
      left_join(p_values, by = c("species", "treatment_level")) %>%
      mutate(
        color = case_when(
          !is.na(p_value) & species == "R. sativus" ~ josef_colors["R. sativus"],
          !is.na(p_value) & species == "B. officinalis" ~ josef_colors["B. officinalis"],
          !is.na(p_value) & species == "H. vulgare" ~ josef_colors["H. vulgare"],
          TRUE ~ josef_colors["grey"]  # Non-significant points in grey
        )
      ) %>% 
      mutate(p_value = round(p_value, 2))  # Round p-values for display
    
    # Check if color assignment worked
    print("Color assignment preview:")
    print(data_rel %>% select(species, treatment_level, r_squared, p_value, color) %>% head(10))
    
    # Generate plot for the current relationship
    plot <- plot_correlation(data_rel, rel$x, rel$y, rel$x_label, rel$y_label) +
      facet_grid(~ treatment_level)
    
    plots[[length(plots) + 1]] <- plot
  }
  return(plots)
}

# Ensure the figures directory exists
if (!dir.exists("./figures")) {
  dir.create("./figures")
}

# TODO: clean up annotations
# - - - 
# Correlation of traits with growth
# Define relationships to plot
growth_relationships <- list(
  list(x = "LMA", y = "dry_whole_g", x_label = expression(LMA~(g~m^{-2})), y_label = expression(Shoot~Biomass~(g))),
  list(x = "LMA", y = "LDMC", x_label = expression(LMA~(g~m^{-2})), y_label = expression(LDMC~(mg~g^{-1}))),
  list(x = "CHL", y = "LMA", x_label = expression(CHL~(mu*g~cm^{-2})), y_label = expression(LMA~(g~m^{-2})))
)

# Generate and save plots for growth relationships
growth_plots <- generate_plots(growth_relationships, data)
pdf("./figures/growth_relationships.pdf", width = 10, height = 10)
print(ggarrange(plotlist = growth_plots, ncol = 1, nrow = length(growth_relationships), common.legend = TRUE, legend = "bottom"))
dev.off()

# - - - 
# Correlation of molecular contents
# Define relationships to plot
molecular_relationships <- list(
  list(x = "CHL", y = "CAR", x_label = expression(CHL~(mu*g~cm^{-2})), y_label = expression(CAR~(mu*g~cm^{-2}))),
  list(x = "CHL", y = "ANT", x_label = expression(CHL~(mu*g~cm^{-2})), y_label = expression(ANT~(mu*g~cm^{-2}))),
  list(x = "CAR", y = "ANT", x_label = expression(CAR~(mu*g~cm^{-2})), y_label = expression(ANT~(mu*g~cm^{-2})))
)

# Generate and save plots for molecular content relationships
molecular_plots <- generate_plots(molecular_relationships, data)
pdf("./figures/molecular_relationships.pdf", width = 10, height = 10)
print(ggarrange(plotlist = molecular_plots, ncol = 1, nrow = length(molecular_relationships), common.legend = TRUE, legend = "bottom"))
dev.off()

# - - - 
# Correlation of leaf structural traits
# To mesophyll structure
mesophyll_relationships <- list(
  list(x = "N", y = "CHL", x_label = "N layers", y_label = expression(CHL~(mu*g~cm^{-2}))),
  list(x = "CHL", y = "LDMC", x_label = expression(CHL~(mu*g~cm^{-2})), y_label = expression(LDMC~(mg~g^{-1}))),
  list(x = "LDMC", y = "LMA", x_label = expression(LDMC~(mg~g^{-1})), y_label = expression(LMA~(g~m^{-2})))
)

# Generate and save plots for molecular content relationships
mesophyll_plots <- generate_plots(mesophyll_relationships, data)
pdf("./figures/mesophyll_relationships.pdf", width = 10, height = 10)
print(ggarrange(plotlist = mesophyll_plots, ncol = 1, nrow = length(mesophyll_relationships), common.legend = TRUE, legend = "bottom"))
dev.off()
