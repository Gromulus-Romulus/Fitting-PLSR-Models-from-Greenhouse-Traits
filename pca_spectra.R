##' Objective: Create a PCA plot of spectral signatures
##' across low and high treatments
##'
##' Adapted from code written for 7th
##' annual plant functional traits course.
##' 
##' @author [Nicole Bison, Nathan Malamud]
##' @date 2024.10.03
##' 
##' Source:
##'   https://github.com/MichaletzLab/pftc7_spectroscopy
##'   

library(FactoMineR)
library(ggfortify)
library(factoextra)
library(scatterplot3d)
library(tidyverse)
library(MetBrewer)
library(dplyr)

# - - - - - 
# Read in input files and remove "date" column
# as it is irrelevant for now.
spectra_wide <- read.csv("./data/spec_data_wide.csv") |>
  select(-date)

# Remove "X" characters from beginning of wavelengths
spec_names <- colnames(spectra_wide)[5:ncol(spectra_wide)]
colnames(spectra_wide)[5:ncol(spectra_wide)] <- gsub("X", "", spec_names)

# Assign wavelength ranges
# Don't include below < 400 and above 2400
# due to increased noise at range edges
# Ranges set using SVC i-series Field Spectroscopy guide
FULL = seq(400, 2400, by = .1)
VIS = seq(400, 700, by = .1)
NIR = seq(700, 1000, by = .1)
SWIR = seq(1000, 2400, by = .1)

# - - - - - 
# Categorize treatments into "LO" and "HI" based on set threshold for mmol
THRESHOLD <- 15
spectra_wide <- spectra_wide %>%
  mutate(treatment_category = ifelse(treatment_mmol <= THRESHOLD, "LO", "HI"))

# Sort values
spectra_wide <- spectra_wide %>%
  select(barcodeID, sampleID, treatment_mmol, treatment_category, species, everything())

spec_names <- colnames(spectra_wide)[6:ncol(spectra_wide)]

# Filter out spectral matrix
X <- spectra_wide %>% select(spec_names) %>%
  as.matrix()

# Filter spec_names for the SWIR range (700 to 1000 nm)
FULL_names <- spec_names[spec_names %in% FULL]

# Filter the spectral matrix for FULL range and convert to matrix
X_FULL <- spectra_wide %>%
  select(all_of(FULL_names)) %>%
  as.matrix()

# PCA
# specs <- c('LO', "HI")
# topca <- spectra_wide %>% filter(species %in% specs)

# - - - - - 
# PCA by species
# TODO: look into PC loadings for prcomp
# Loadings tell you what the PCA stands for
pca <- prcomp(X_FULL, scale=TRUE)

pc <- as.data.frame(pca$x)
pc$species <- spectra_wide$species
pc$barcode  <- spectra_wide$barcodeID
pc$treatment_category <- spectra_wide$treatment_category

pcaplot <- ggplot(data = pc,
                  aes(x = PC1, y=PC2, color = factor(treatment_category))) +
  geom_point() + 
  stat_ellipse() +
  #scale_color_met_d("Lakota") +
  theme_minimal() + facet_wrap(~species)
pcaplot

# How many compoenents needed
eig <- fviz_eig(pca, 
                addlabels = TRUE, 
                ylim = c(0, 70),
                main="Scree Plot")
eig

# - - - - -
library(ggplot2)
library(ggrepel)

# TODO: look into PC loadings for prcomp
# Extract PC loadings from the PCA object
loadings <- as.data.frame(pca$rotation)

# View the first few rows of loadings
head(loadings)

# Add variable names to the loadings data frame
loadings$wavelength <- rownames(loadings)

# Plot the loadings for PC1 and PC2
loading_plot <- ggplot(loadings, aes(x = PC1, y = PC2, label = wavelength)) +
  geom_point() +
  geom_text_repel(size = 3) +
  theme_minimal() +
  labs(title = "PC Loadings Plot",
       x = "PC1 Loadings",
       y = "PC2 Loadings")

loading_plot