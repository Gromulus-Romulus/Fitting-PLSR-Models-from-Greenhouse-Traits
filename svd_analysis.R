##' Objective: Rank spectral frequencies by scores and loadings
##' obtained from fitted PLSR models using pls R package.
##' 
##' @author Nathan Malamud
##' @date 2024.10.10
##' 
##'   

library(dplyr)
library(tidyverse)
library(prospectr)
library(caret)
library(pls)
library(ggplot2)
library(ggpubr)

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

# Load traits data, only include barcodes that we
# have spectral curves for
traits <- read.csv("./data/traits.csv") %>%
  filter(barcodeID %in% signal.matrices$uniq_ids)

# Save barcode metadata in separate dataframe
# Extract prediction matrix Y
Y <- traits %>%
  select("LMA", "LDMC", "EWT") %>% as.matrix()

# Save metadata in separate dataframe
metadata <- traits %>%
  select(barcodeID, sampleID, treatment_mmol, species)

# Fit PLSR model using Shawn Serbin's demonstration code:
# Documentation guide (Serbin et al. 2022): 
#   https://github.com/plantphys/spectratrait/blob/main/spectratrait_1.2.5.pdf 
plsr_data <- cbind(barcodeID, Y, X_avg) %>%
  as.data.frame()

# Traits and species of interest
trait_ids <- c("LMA", "LDMC", "EWT")
units_lookup <- c("LMA" = "g/m^2", "LDMC" = "mg/g", "EWT" = "g/m^2") 
species_ids <- traits$species %>% unique()

# Create a list to store plots
# And save the models as well!
p <- list()
mods <- list()

# Max iterations for SVD
P_DIM = 10 

# Create a PDF device to save the plots
pdf("Figures/PLSR_Loadings_Xavg.pdf", width = 10, height = 8.5)

# Iterate: for every model, one predicted trait
for (t in trait_ids) {
  
  # Extract y and x matrices for model fitting
  y_obs = plsr_data[[t]]
  x = plsr_data[, spectral_columns] %>% as.matrix()
  
  # TODO: experiment with different cross validation methods
  plsr_model <- plsr(y_obs ~ x, ncomp=P_DIM, validation="LOO")
  
  # Save model to mods
  mods[[t]] <- plsr_model
  
  # Which number of components worked best?
  nComps <- selectNcomp(plsr_model, plot=F)
  
  # Extract fitted values and observed values
  y_hat <- predict(plsr_model, x, ncomp=nComps) %>% as.vector()
  
  # Calculate RMSEP and R2 values
  # Rounded to 3 decimal places
  # TODO: Look at PRESS Statistics (Shawn Serbin paper)
  resids = y_obs - y_hat
  y_bar = mean(y_obs)
  RSS = sum(resids**2)
  TSS = sum((y_obs - y_bar)**2)
  r2 = round(1.0 - (RSS/TSS), 3)
  rmsep = round(sqrt(mean(resids**2)), 3)
  
  # TODO: this can only be calculated if validation method is LOO
  press <- round(sum(plsr_model$validation$PRESS), 3)
  
  # Create a data frame for ggplot
  plot_data <- data.frame(observed = y_obs, predicted = y_hat, species = metadata$species)
  
  # Plot with color by species
  josef_colors <- c("#7570b2", "#ca621c", "#299680")  # Custom color palette for species
  
  plot <- ggplot(plot_data, aes(x = observed, y = predicted, color = species)) +  # Color points by species
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(title = bquote(atop(.(paste(t, "(comps =", nComps, ")")),
                             atop("RMSEP:" ~ .(rmsep) ~ ", R2:" ~ .(r2) ~ ", PRESS:" ~ .(press)))),
         x = "Observed Values",
         y = "Predicted Values") +
    theme_minimal() +
    coord_fixed(ratio = 1) +  # Keep the aspect ratio 1:1
    xlim(range(plot_data$observed)) +  # Set x limits based on observed data
    ylim(range(plot_data$observed)) +  # Set y limits based on observed data
    scale_color_manual(values = josef_colors)  # Apply custom colors to species
  
  # Add plot to list of plots
  p[[t]] <- plot
  
  # - - - use separately fit models to extract scores and loadings
  scores <- plsr_model$scores %>% as.matrix()
  loadings <- plsr_model$loadings %>% as.matrix()
  
  plot(plsr_model, plottype = "loadings", comps = 1:3, legendpos ="topleft", labels="numbers", xlab="nm")
  title(paste(t, "Loadings"))
  
  abline(h=0)
}

# Combine the plots into a 3 x 3 grid (or as many as needed)
grid_plot <- ggarrange(plotlist = p, ncol = 3, nrow = 1, common.legend = TRUE)

# Display the grid plot
# TODO: look into plotting warnings
print(grid_plot)

dev.off()



