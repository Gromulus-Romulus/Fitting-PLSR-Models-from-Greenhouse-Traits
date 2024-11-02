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

# - - - - - 
# Structural trait PCA
data <- read.csv("./data/traits.csv")

# Merge with prospect output
prospect <- read.csv("./data/molecular_content.csv") %>%
  select(barcodeID, CHL, N)

data <- merge(data, prospect, by = "barcodeID")

# Create a new column for aggregated treatment levels
data <- data %>%
  mutate(treatment_level = case_when(
    treatment_mmol <= 5 ~ "0 - 5 mmol",
    treatment_mmol <= 15 ~ "10 - 15 mmol",
    treatment_mmol <= 25 ~ "20 - 25 mmol",
    treatment_mmol <= 35 ~ "30 - 35 mmol",
    TRUE ~ "Other"  # Optional: catch any values above 35
  ))

# Calculate "theoretical" Electron Transport Rate (ETR)
data$ETR <- (data$Phi_PS2 / data$Qamb) * 500

pca_data <- data %>%
  select(barcodeID, species, treatment_level, ETR, CHL, N, LMA, LDMC, dry_whole_g)

pca <- prcomp(pca_data %>% select(ETR, CHL, N, LMA, LDMC, dry_whole_g), center = TRUE, scale = TRUE)

# Biplot to visualize PCA results
pca_biplot <- fviz_pca_biplot(pca, 
                              geom.ind = "point", # Use points for individuals
                              col.ind = pca_data$species, # Color by species
                              addEllipses = TRUE, # Add concentration ellipses
                              legend.title = "Species",
                              repel = TRUE) # Avoid label overlap

pca_biplot

# Define colors for each species
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c", "grey" = "#d3d3d3")

# Extract PCA scores and loadings
pca_scores <- as.data.frame(pca$x)
pca_scores$barcodeID <- pca_data$barcodeID
pca_scores$species <- pca_data$species
pca_scores$treatment_level <- pca_data$treatment_level

pca_loadings <- as.data.frame(pca$rotation)
pca_loadings$Variable <- rownames(pca_loadings)

# Define colors for each species
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c", "grey" = "#d3d3d3")

# Load required libraries
library(ggplot2)
library(dplyr)
library(factoextra)

# Define colors for each species
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c", "grey" = "#d3d3d3")

# Split data by species
species_list <- unique(pca_data$species)
pca_results <- list()

# Perform separate PCA for each species
for (species in species_list) {
  # Filter data for each species
  species_data <- pca_data %>% filter(species == !!species)
  
  # Perform PCA on the subset of data
  pca <- prcomp(species_data %>% select(ETR, CHL, N, LMA, LDMC, dry_whole_g), center = TRUE, scale = TRUE)
  
  # Store PCA results in a list
  pca_results[[species]] <- pca
}

# Create PCA biplots for each species
pca_plots <- list()
# Adjusted code to avoid row mismatch

# Perform separate PCA for each species
for (species in species_list) {
  # Filter data for each species
  species_data <- pca_data %>% filter(species == !!species)
  
  # Perform PCA on the subset of data
  pca <- prcomp(species_data %>% select(ETR, CHL, N, LMA, LDMC, dry_whole_g), center = TRUE, scale = TRUE)
  
  # Extract PCA scores and add treatment_level information
  pca_scores <- as.data.frame(pca$x)
  pca_scores$treatment_level <- species_data$treatment_level
  
  # Extract PCA loadings
  pca_loadings <- as.data.frame(pca$rotation)
  pca_loadings$Variable <- rownames(pca_loadings)
  
  # Plot PCA
  pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = treatment_level)) +
    geom_point(size = 2) +
    geom_segment(data = pca_loadings, 
                 aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5), # Adjust multiplier for arrow length
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "black") +
    geom_text(data = pca_loadings, 
              aes(x = PC1 * 5.5, y = PC2 * 5.5, label = Variable), # Adjust multiplier for label position
              color = "black", 
              size = 3, 
              hjust = 0.5, 
              vjust = 0.5) +
    labs(title = paste("PCA for", species), x = "PC1", y = "PC2") +
    theme_minimal()
  
  # Add plot to list
  pca_plots[[species]] <- pca_plot
}

# Display plots for each species
pca_plots


# Biplot to visualize PCA results
pca_biplot <- fviz_pca_biplot(pca, 
                              geom.ind = "point", # Use points for individuals
                              col.ind = pca_data$treatment_level, # Color by species
                              addEllipses = TRUE, # Add concentration ellipses
                              legend.title = "Species",
                              repel = TRUE) # Avoid label overlap

pca_biplot

# How many compoenents needed
eig <- fviz_eig(pca, 
                addlabels = TRUE, 
                ylim = c(0, 70),
                main="Scree Plot")
eig





