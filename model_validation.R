##' Model cross validation - calculate PRESS statistics.
##'
##' @author Nathan Malamud
##' @date 2024.10.04
##' 

library(dplyr)
library(tidyverse)
library(prospectr)
library(caret)
library(spectratrait)
library(pls)
library(ggplot2)
library(ggpubr)  # Load ggpubr for ggarrange

# Load traits data
traits <- read.csv("./data/traits.csv")

# Load spectral matrices for PLSR model (made by signal_cleaning.R)
signal.matrices <- readRDS("./data/signal.matrices.rds")

# - - - - -
# 1) TODO: See which matrices (X_raw, X_avg, X_d1, X_d2)
# Reduce RMSEP error multiple traits (e.g. LMA, LDMC, EWT)
# ... by using PLSR model with different signal matrices
X_raw <- signal.matrices$X_raw
X_avg <- signal.matrices$X_avg
X_d1 <- signal.matrices$X_d1
X_d2 <- signal.matrices$X_d2

# Important - the uniq_ids array
# helps us translate between matrix row numbers
# and barcodeIDs for retrieving sample metadata
barcodeID <- signal.matrices$uniq_ids
spectral_columns <- colnames(X_raw) %>% as.character()

# Extract prediction matrix Y
Y <- traits %>%
  filter(barcodeID %in% signal.matrices$uniq_ids) %>%
  select("LMA", "LDMC", "EWT") %>% as.matrix()

# Fit PLSR model using Shawn Serbin's demonstration code:
# Documentation guide (Serbin et al. 2022): 
#   https://github.com/plantphys/spectratrait/blob/main/spectratrait_1.2.5.pdf 
plsr_data <- cbind(barcodeID, Y, X_d2) %>%
  as.data.frame()

# Traits of interest
trait_ids <- c("LMA", "LDMC", "EWT")
species_ids <- traits$species %>% unique()

# Create a list to store plots
p <- list()

# TODO: this code does not work due oto a misuse of the PLSR function
# y fitted values are not being extracted correctly

# Iterate: for every model, one predicted trait
for (t in trait_ids) {
  
  # TODO: Determine optimal number of components using PCA analysis
  n = 3
  
  # Extract y and x matrices for model fitting
  y_obs = plsr_data[[t]]
  x = plsr_data[, spectral_columns] %>% as.matrix()
  
  plsr_model <- plsr(y_obs ~ x, ncomp=3)
  
  # TODO: why are we outputing matrices, and not vectors???
  # Extract fitted values and observed values
  y_hat <- plsr_model$fitted.values |> as.vector()
  
  # Create a data frame for ggplot
  plot_data <- data.frame(observed = y_obs, predicted = y_hat)
  
  # TODO: calculate R2 and RMSEP values

  # Create observed vs predicted regression plot
  plot <- ggplot(plot_data, aes(x = predicted, y = observed)) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(title = paste(t),
         x = "Predicted Values",
         y = "Observed Values") +
    theme_minimal() +
    coord_fixed(ratio = 1) +  # Keep the aspect ratio 1:1
    xlim(range(plot_data$observed)) +  # Set x limits based on observed data
    ylim(range(plot_data$observed))    # Set y limits based on observed data

  # Add plot to list of plots
  p[[t]] <- plot
}

# Combine the plots into a 3 x 3 grid (or as many as needed)
grid_plot <- ggarrange(plotlist = p, ncol = 3, nrow = 1, common.legend = TRUE)

# Display the grid plot
print(grid_plot)
