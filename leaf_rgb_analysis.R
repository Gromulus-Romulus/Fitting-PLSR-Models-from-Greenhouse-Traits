# See if leaf greenness index (calculated w/ MatLab)
# correlates with growth rate or leaf index
##' @author Nathan Malamud
##' @date 2024-09-24
##' 
##' Indices calculated here:
##' - NGI: Normalized Greenness Index
##' TODO: - green_difference: TODO define
##' TODO:- green_difference2: TODO defne
##' - veg_index: Vegetation index
##' - greenness: Greenness index
##' 
##' Source:
##'   https://www.biorxiv.org/content/10.1101/2023.08.23.554481v1.full

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(GGally)
library(readxl)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load data from RDS files
traits <- readRDS("./proc/traits.RDS")
spec <- readRDS("./proc/spec.RDS")

# Leaf color percentages calculated from MatLab
leaf_color <- read_excel("./raw/image_rgb.xlsx")

# merge with traits data
data <- merge(traits, leaf_color, by="barcodeID")

# CHL, CAR, EWT content from PROSPECT-PRO
RT_output <- read_csv("./prospect_rtm_output.csv")
data <- merge(data, RT_output %>% select(-c(sampleID, species, treatment_mmol)), by="barcodeID")

# TODO: Organize dataframe columns, investigate warning
metadata <- c("barcodeID", "species", "treatment_mmol", "sampleID")
data <- data %>% select(metadata, everything())

# Calculate Chlorophyll Index (CI)
# Source: https://eos.com/make-an-analysis/chlorophyll-index/
# CI Green = ρNIR / ρGreen – 1 = ρ730 / ρ530 – 1
# CI Red-Edge = ρNIR / ρRed-Edge – 1 = ρ850 / ρ730 – 1
# Merge CI Green and CI Red with Data
CIs <- spec %>% select(`529.8`, `729.9`, `849.9`, barcodeID)
CIs$CI_green <- (CIs$`729.9` / CIs$`529.8`) - 1
CIs$CI_red <- (CIs$`849.9` / CIs$`729.9`) - 1
data <- merge(data, CIs, by = "barcodeID")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Visualization
# Create two plots, one relating percentGreen to total dry biomass
# another relating leaf area / total dry biomass (La)
# External validation fit stats
# Box Plot of Leaf CI Green in Response to Treatment
# TODO: Consider Using Loess Smooth Curves Instead of Boxplots?
# TODO: Ensure Colors Are Consistent with Josef's Manuscript
josef_colors <- c("#299680", "#7570b2", "#ca621c")

# Pariwise feature plot of Chlorophyll indices green with CHL (RT estimated)
data %>%
  select(CHL, CI_green, CI_red, NGI, greenness, veg_index, green_difference) %>%
  ggpairs(
    lower = list(continuous = wrap("points", alpha = 0.5, color = josef_colors[1])),
    upper = list(continuous = wrap("cor", size = 5, color = josef_colors[2])),
    diag = list(continuous = wrap("barDiag", fill = josef_colors[3]))
  ) + 
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend if it's not needed

# Pariwise feature plot of CI green with traits
data %>%
  select(CI_green, dry_whole_g, area_cm2, LDMC, LMA, CHL) %>%
  ggpairs(
    lower = list(continuous = wrap("points", alpha = 0.5, color = josef_colors[1])),
    upper = list(continuous = wrap("cor", size = 5, color = josef_colors[2])),
    diag = list(continuous = wrap("barDiag", fill = josef_colors[3]))
  ) + 
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend if it's not needed

# Box plots across species and treatments
data %>%
  ggplot(aes(y = CHL, x = treatment_mmol, color = species)) +
  geom_boxplot(aes(alpha = 0.25)) +
  facet_wrap(~species, ncol = 3) +
  scale_color_manual(values = josef_colors) +
  theme_classic() +
  theme(legend.position = "none")

# Scatterplot of CI Green vs. Dry Biomass with R² and equation
data %>%
  ggplot(aes(x = NGI, y = dry_whole_g, color = species)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = josef_colors) +
  theme_classic() +
  facet_wrap(~species, nrow = 3) +
  theme(legend.position = "bottom")  # Remove legend if it's not needed
  
# summary of mixed effects model
library(lmerTest)
z <- lmer(dry_whole_g ~ CI_green + (1|species), data = data)
summary(z)
  
# Histogram of Size Distributions for Species
den1 <- data %>%
  ggplot(aes(dry_whole_g, group = treatment_mmol, fill = as.numeric(treatment_mmol))) +
  geom_density(alpha = 0.25) +
  theme_classic() +
  facet_wrap(.~species) +
  theme(legend.position = "e")

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
  
