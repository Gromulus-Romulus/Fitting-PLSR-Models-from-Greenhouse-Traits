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

# Create a lookup table for species to Latin names
latin_names <- c("borage" = "B. officinalis",
                 "barley" = "H. vulgare",
                 "radish" = "R. sativus")

# Rename species column using the lookup table
data <- data %>%
  mutate(species = recode(species,
                          "borage" = latin_names["borage"],
                          "barley" = latin_names["barley"],
                          "radish" = latin_names["radish"]))

# Pariwise feature plot of Chlorophyll indices green with CHL (RT estimated)
data %>%
  select(CHL, CI_green, CI_red, NGI, veg_index) %>%
  ggpairs(
    lower = list(continuous = wrap("points", alpha = 0.5, color = josef_colors[1])),
    upper = list(continuous = wrap("cor", size = 5, color = josef_colors[2])),
    diag = list(continuous = wrap("barDiag", fill = josef_colors[3]))
  ) + 
  theme(text = element_text(family = "Helvetica", face = "bold"), legend.position = "none") +
  theme_minimal()  # Remove legend if it's not needed

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
  ggplot(aes(y = CI_green, x = treatment_mmol, color = species)) +
  geom_boxplot(aes(alpha = 0.25)) +
  facet_wrap(~species, ncol = 3) +
  scale_color_manual(values = josef_colors) +
  theme_minimal() +
  labs(x = "Ammoniacal N (mmol)", y = "CI (ρ730 / ρ530 – 1)") +
  theme(
    text = element_text(family = "Helvetica", face = "plain"),  # Regular font for all text
    axis.title = element_text(size = 12, face = "bold"),  # Axis titles in regular
    axis.text = element_text(size = 10, face = "plain"),   # Axis text in regular
    legend.title = element_text(size = 12, face = "bold"),# Legend title in regular
    legend.text = element_text(size = 10, face = "bold.italic"), # Legend text in regular
    legend.position = "None",
    strip.text = element_text(size = 12, face = "bold.italic")  # Species names in italic in facet labels
  )

# summary of mixed effects model
library(lmerTest)
z <- lmer(CHL ~ treatment_mmol + (1|species), data = data)
summary(z)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Scatterplot of CI Green vs. Dry Biomass with R² and equation

data %>%
  ggplot(aes(x = CI_green, y = dry_whole_g, color = species, shape=species)) +
  geom_point(alpha=0.75) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = josef_colors) +
  labs(x = "Chlorophyll Index (ρ730 / ρ530 – 1)", y = "Dry Biomass (g)") +
  theme_minimal() +
  #facet_wrap(~species, nrow = 3) +
  theme(
    text = element_text(family = "Helvetica", face = "plain"),  # Regular font for all text
    axis.title = element_text(size = 12, face = "bold"),  # Axis titles in regular
    axis.text = element_text(size = 10, face = "plain"),   # Axis text in regular
    legend.title = element_text(size = 12, face = "bold"),# Legend title in regular
    legend.text = element_text(size = 10, face = "bold.italic"), # Legend text in regular
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold.italic")  # Species names in italic in facet labels
  )
  
# summary of mixed effects model
library(lmerTest)
z <- lmer(dry_whole_g ~ CI_green + (1|species), data = data)

z <- lm(dry_whole_g ~ CI_green + species, data=data)
summary(z)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Density plots of structural traits (CHL, LDMC, N)

# Convert treatment_mmol to a factor with levels ordered by nitrate concentration
data <- data %>%
  mutate(treatment_mmol = factor(treatment_mmol, levels = c(0, 5, 10, 15, 20, 25, 30, 35)))

# Create a custom color palette for the discrete categories
color_palette <- c("lightblue", "skyblue", "deepskyblue", "dodgerblue", 
                   "blue", "darkblue", "navy", "midnightblue")

# Create density plots using log density for LMA
den1 <- data %>%
  ggplot(aes(CHL, group = treatment_mmol, fill = treatment_mmol)) +
  geom_density(aes(y = ..density..), alpha = 0.25) +
  theme_classic() +
  facet_wrap(.~species) +
  scale_fill_manual(values = color_palette) +
  theme(legend.position = "e")

den2 <- data %>%
  ggplot(aes(LDMC, group = treatment_mmol, fill = treatment_mmol)) +
  geom_density(aes(y = ..density..), alpha = 0.25) +
  theme_classic() +
  facet_wrap(.~species) +
  scale_fill_manual(values = color_palette) +
  theme(legend.position = "none")

# Use log density for LMA plot
den3 <- data %>%
  ggplot(aes(N, group = treatment_mmol, fill = treatment_mmol)) +
  geom_density(aes(y = ..density..), alpha = 0.25) +
  theme_classic() +
  facet_wrap(.~species) +
  scale_fill_manual(values = color_palette) =
  theme(legend.position = "none")

# Arrange the plots in a grid
ggarrange(den1, den2, den3, nrow = 3)

