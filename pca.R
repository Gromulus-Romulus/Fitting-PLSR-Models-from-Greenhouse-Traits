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
# Load spectral matrices for PLSR model (made by signal_cleaning.R)
signal.matrices <- readRDS("./data/signal.matrices.rds")

X_raw <- signal.matrices$X_raw
X_avg <- signal.matrices$X_avg
X_d1 <- signal.matrices$X_d1
X_d2 <- signal.matrices$X_d2

# Assign wavelength ranges
# Don't include below < 400 and above 2400
# due to increased noise at range edges
# Ranges set using SVC i-series Field Spectroscopy guide
FULL = seq(400, 2400, by = .1)
VIS = seq(400, 700, by = .1)
NIR = seq(700, 1000, by = .1)
SWIR = seq(1000, 2400, by = .1)

pca <- prcomp(X_raw, center = TRUE, scale = TRUE)

# How many compoenents needed
eig <- fviz_eig(pca, 
                addlabels = TRUE, 
                ylim = c(0, 70),
                main="Scree Plot")
eig

# - - - - -
library(ggplot2)
library(ggrepel)

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

