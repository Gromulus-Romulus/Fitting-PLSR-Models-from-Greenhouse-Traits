# Author: Nathan Malamud
# Date: 2024.05.24

# Analysis Packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readxl)
library(reshape2)
library(scales)
library(LeafArea)
library(dplyr)
library(stringr)
library(RColorBrewer)

# REMINDER: Set Working Directory -> Source File Location
data <- read_excel("./leaf_mass_data.xlsx")

# Define Factor Levels (Treatment and Species)
data$treatment_mmol <- as.factor(data$treatment_mmol)
data$species <- as.factor(data$species)
levels(data$species) <- c("radish", "borage", "barley")
data$LDMC <- as.numeric(data$LDMC)

# Run ImageJ Software on Leaf Scans
# This line of code initiates a pop-up window.
# It will take a while to run, that's expected.
# Afterwards, merge values for image scans with the dataframe.
scans <- run.ij(set.directory = './scans')
names(scans) <- c("barcodeID", "area_cm2")
data <- merge(data, scans, by = "barcodeID")

# Calculate LMA and SLA from Weight and Area Measurements
data$SLA <- (data$area_cm2 / 100) / data$dry_leaf_g

# Load Fluorometry Data from Porometer and Rename Columns
# Use Fs and Fm' to Calculate Quantum Yield of Fluorescence
# Merge Porometer Measurements with Data Afterwards
# Source: https://www.licor.com/env/products/LI-600/
fluor <- read_csv("./fluor/fluor_data.csv") %>%
  select(c(unique_id, gsw, Fs, `Fm'`, PhiPS2))
names(fluor) <- c("barcodeID", "gsw", "Fs", "Fm_prime", "Phi_PS2")
data <- merge(data, fluor, by = "barcodeID")

# Load Spectroscopy Data from Porometer and Aggregate Data
# N Bison Mentioned There Are 2-3 Measurements per Individual
# TODO: Investigate Warnings Here
spec <- read_csv("./spec/spec_data.csv") %>%
  aggregate(by = list(.$barcodeID), FUN = mean)
names(spec) <- c("barcodeID", "delete_me", names(spec[3:1026]))
spec <- spec %>% select(-delete_me)

# Calculate Chlorophyll Index (CI)
# Source: https://eos.com/make-an-analysis/chlorophyll-index/
# CI Green = ρNIR / ρGreen – 1 = ρ730 / ρ530 – 1
# CI Red-Edge = ρNIR / ρRed-Edge – 1 = ρ850 / ρ730 – 1
# Merge CI Green and CI Red with Data
CIs <- spec %>% select(`529.8`, `729.9`, `849.9`, barcodeID)
CIs$CI_green <- (CIs$`729.9` / CIs$`529.8`) - 1
CIs$CI_red <- (CIs$`849.9` / CIs$`729.9`) - 1
data <- merge(data, CIs, by = "barcodeID")

# Box Plot of Leaf CI Green in Response to Treatment
# TODO: Consider Using Loess Smooth Curves Instead of Boxplots?
# TODO: Ensure Colors Are Consistent with Josef's Manuscript
josef_colors <- c("#299680", "#7570b2", "#ca621c")

data.long <- data %>%
  melt(id = c("species", "treatment_mmol", "barcodeID", "sampleID")) %>%
  select(species, treatment_mmol, variable, value)

# TODO: Investigate Warning Message
data.long$value <- as.numeric(data.long$value)

# Boxplot Visualization
data.long %>%
  subset(variable %in% c("CI_green", "SLA", "dry_leaf_g")) %>%
  ggplot(aes(y = value, x = treatment_mmol, color = species)) +
  geom_boxplot(aes(alpha = 0.25)) +
  facet_wrap(variable ~ species, scales = "free", nrow = 3) +
  scale_color_manual(values = josef_colors) +
  ylab("") +
  theme_classic() +
  theme(legend.position = "none")

# Load Color Percentages from Images
image_colors <- read_csv("./leaf_color_percentages.csv")
data <- merge(data, image_colors, by = "barcodeID")

data.colors.long <- data %>%
  melt(id = c("species", "treatment_mmol", "barcodeID", "sampleID")) %>%
  select(species, treatment_mmol, variable, value)
data.colors.long$value <- as.numeric(data.colors.long$value)

TOT_PIX_RES <- 2480 * 3507

# Boxplot of CI Green Across Treatment Levels
data.colors.long %>%
  subset(variable %in% c("CI_green")) %>%
  ggplot(aes(y = value, x = as.factor(treatment_mmol), color = species)) +
  geom_boxplot(aes(alpha = 0.15)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "loess", alpha = 0.35, linetype = "dashed") +
  facet_wrap(variable ~ species, ncol = 3) +
  scale_color_manual(values = josef_colors) +
  ylab("") +
  theme_classic() +
  theme(legend.position = "none")

# Boxplot of Dry Whole Mass by Treatment
data %>%
  ggplot(aes(y = dry_whole_g, x = treatment_mmol, color = species)) +
  geom_boxplot(aes(alpha = 0.25)) +
  facet_wrap(~species, ncol = 3) +
  scale_color_manual(values = josef_colors) +
  theme_classic() +
  theme(legend.position = "none")

# Histogram of Size Distributions for Species
den1 <- data %>%
  ggplot(aes(dry_whole_g, group = treatment_mmol, fill = as.numeric(treatment_mmol))) +
  geom_density(alpha = 0.25) +
  theme_classic() +
  facet_wrap(.~species) +
  theme(legend.position = "none")

den2 <- data %>%
  ggplot(aes(dry_leaf_g, group = treatment_mmol, fill = as.numeric(treatment_mmol))) +
  geom_density(alpha = 0.25) +
  theme_classic() +
  facet_wrap(.~species) +
  theme(legend.position = "none")

den3 <- data %>%
  ggplot(aes(SLA, group = treatment_mmol, fill = as.numeric(treatment_mmol))) +
  geom_density(alpha = 0.25) +
  theme_classic() +
  facet_wrap(.~species) +
  theme(legend.position = "none")

ggarrange(den1, den2, den3, nrow = 3)

# Plot Spectra Data
spec <- merge(spec, data[c("barcodeID", "sampleID", "species", "treatment_mmol")], by = "barcodeID")

# Average Spectral Signature per Species, per Treatment
spec <- spec %>% select(-c(barcodeID, sampleID))

spec <- spec %>%
  melt(id = c("species", "treatment_mmol"))
names(spec) <- c("species", "treatment_mmol", "wavelength", "reflectance")
spec$wavelength <- as.numeric(spec$wavelength)
spec$reflectance <- as.numeric(spec$reflectance)

# Reflectance Curves per Species and Treatment
# How Reflectance Curves Change Across Treatments for Each Species
spec %>%
  subset(treatment_mmol %in% c(5, 25, 35)) %>%
  ggplot(aes(x = wavelength, y = reflectance, color = as.numeric(treatment_mmol), alpha = as.numeric(treatment_mmol))) +
  geom_point() +
  facet_wrap(species ~ ., ncol = 1) +
  theme_classic() +
  xlim(530, 800)
