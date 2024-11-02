# Examine random effects of treatment levels
# across the nitrogen gradient
# Examine relationships between CHL and growth rate
# across treatments as well as ETR and growth rate
#
# author: Nathan Malamud
# date: 2024.10.29

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

# Fit random effects model to examine relationships
library(lme4)
library(visreg)
z <- lmer(dry_whole_g ~ LMA*CHL + species + (1|treatment_level), data = data)
summary(z)

# Plot N vs LMA across treatments
p1 <- ggplot(data, aes(x = N, y = LMA, color = treatment_level)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Number of layers", y = "Leaf Mass Area (g / m^2)") +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.25, color = "grey80"),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5, color = "black"),
    strip.text = element_text(hjust = 0, size = 10, face = "italic"),
    legend.position = "bottom",
    legend.title = element_blank()
  ) + facet_wrap(~species) + scale_color_brewer(palette = "Greens")

# Plot CHL vs AB (g) across treatments
p2 <- ggplot(data, aes(x = CHL, y = LDMC, color = treatment_level)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "CHL (ug)", y = "Aboveground Biomass (g)") +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.25, color = "grey80"),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5, color = "black"),
    strip.text = element_text(hjust = 0, size = 10, face = "italic"),
    legend.position = "bottom",
    legend.title = element_blank()
  ) + facet_wrap(~species) + scale_color_brewer(palette = "Greens")

# Plot LMA vs AB (g) across treatments
p3 <- ggplot(data, aes(x = N, y = LDMC)) +
  geom_point(color = data$treatment_level) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Leaf Dry Matter Content (mg / g)", y = "Aboveground Biomass (g)") +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.25, color = "grey80"),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5, color = "black"),
    strip.text = element_text(hjust = 0, size = 10, face = "italic"),
    legend.position = "bottom",
    legend.title = element_blank()
  ) + facet_wrap(~species) + scale_color_brewer(palette = "Greens")

# Save plots to PDF
# TODO: investigate warning messages
pdf("./figures/random_effects_plots.pdf", width = 8, height = 8)
print(ggarrange(p1, p2, p3, ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom"))
dev.off()
