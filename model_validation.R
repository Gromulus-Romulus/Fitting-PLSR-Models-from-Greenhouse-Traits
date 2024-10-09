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
plsr_data <- cbind(barcodeID, Y, X_d2) %>%
  as.data.frame()

# Traits of interest
trait_ids <- c("LMA", "LDMC", "EWT")
species_ids <- traits$species %>% unique()

# Create a list to store plots
p <- list()

# TODO: this code does not work due to a misuse of the PLSR function
# y fitted values are not being extracted correctly

# Iterate: for every model, one predicted trait
for (t in trait_ids) {
  
  # TODO: Determine optimal number of components using PCA analysis
  # Select 5 for now.
  # (Q): number of principal components different for each trait?
  P_DIM = 5 
  
  # Extract y and x matrices for model fitting
  y_obs = plsr_data[[t]]
  x = plsr_data[, spectral_columns] %>% as.matrix()
  
  # TODO: experiment with different cross validation methods
  plsr_model <- plsr(y_obs ~ x, ncomp=P_DIM)
  
  # TODO: why are we outputing matrices, and not vectors???
  # Extract fitted values and observed values
  y_hat <- predict(plsr_model, x, ncomp=P_DIM) %>% as.vector()
  
  # TODO: use metrics package
  # Calculate RMSEP and R2 values
  # Rounded to 3 decimal places
  resids = y_obs - y_hat
  y_bar = mean(y_obs)
  RSS = sum(resids**2)
  TSS = sum((y_obs - y_bar)**2)
  r2 = round(1.0 - (RSS/TSS), 3)
  rmsep = round(sqrt(mean(resids**2)), 3)
  
  # Create a data frame for ggplot
  plot_data <- data.frame(observed = y_obs, predicted = y_hat)

  # Create observed vs predicted regression plot
  plot <- ggplot(plot_data, aes(x = observed, y = predicted)) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(title = paste(t, "|", "RMSEP:", rmsep, ",", "R2:", r2, ",", "comp:", P_DIM),
         x = "Observed Values",
         y = "Predicted Values") +
    theme_minimal() +
    coord_fixed(ratio = 1) +  # Keep the aspect ratio 1:1
    xlim(range(plot_data$observed)) +  # Set x limits based on observed data
    ylim(range(plot_data$observed))    # Set y limits based on observed data
  # Add plot to list of plots
  p[[t]] <- plot
}

# Combine the plots into a 3 x 3 grid (or as many as needed)
grid_plot <- ggarrange(plotlist = p, ncol = 3, nrow = 1, common.legend = TRUE)
    # TODO: ggtitle(paste("PLSR Model Performance by Trait", "n=", length(plsr_data), "comps=", P_DIM))

# Display the grid plot
print(grid_plot)

# - - - - -
# 2) TODO: Comparitive statistics using cross-validation
# Calculate PRESS statistics, RMSEP, R2, for each model
THRESH = 15 # mmol
metadata$treatment_category <- ifelse(metadata$treatment_mmol > THRESH, "high", "low")

plsr_data_hi <- plsr_data %>%
  filter(metadata$treatment_category == "high") # 52 obs > 15 mmol

plsr_data_lo <- plsr_data %>%
  filter(metadata$treatment_category == "low") # 45 obs <= 15 mmol

# Define the traits of interest
trait_ids <- c("LMA", "LDMC", "EWT")
P_DIM <- 5 # Number of components for PLSR

# Create a list to store plots for all traits
p_list <- list()

# Loop over each trait and perform PLSR for "FULL", "LO", and "HI" data partitions
for (t in trait_ids) {
  
  # 1) Define observed and spectral matrices for each partition
  y_obs <- plsr_data[[t]]
  y_obs_hi <- plsr_data_hi[[t]]
  y_obs_lo <- plsr_data_lo[[t]]
  
  x <- plsr_data[, spectral_columns] %>% as.matrix()
  x_hi <- plsr_data_hi[, spectral_columns] %>% as.matrix()
  x_lo <- plsr_data_lo[, spectral_columns] %>% as.matrix()
  
  # 2) Fit PLSR models for each partition (calibrate to separate data sets)
  plsr_model_full <- plsr(y_obs ~ x, ncomp = P_DIM)
  plsr_model_hi <- plsr(y_obs_hi ~ x_hi, ncomp = P_DIM)
  plsr_model_lo <- plsr(y_obs_lo ~ x_lo, ncomp = P_DIM)
  
  # 3) Predict and extract fitted values for each model (fit to WHOLE data)
  y_hat_full <- predict(plsr_model_full, x, ncomp = P_DIM) %>% as.vector()
  y_hat_hi <- predict(plsr_model_hi, x, ncomp = P_DIM) %>% as.vector()
  y_hat_lo <- predict(plsr_model_lo, x, ncomp = P_DIM) %>% as.vector()
  
  # 4) Calculate RMSEP, R², and PRESS statistics for each partition
  # Use a helper function to encapsulate the series of calculations
  calc_stats <- function(y_obs, y_hat) {
    resids <- y_obs - y_hat
    RSS <- sum(resids**2)
    y_bar <- mean(y_obs)
    TSS <- sum((y_obs - y_bar)**2)
    r2 <- round(1.0 - (RSS / TSS), 3)
    rmsep <- round(sqrt(mean(resids**2)), 3)
    return(list(rmsep = rmsep, r2 = r2))
  }
  
  stats_full <- calc_stats(y_obs, y_hat_full)
  stats_hi <- calc_stats(y_obs, y_hat_hi)
  stats_lo <- calc_stats(y_obs, y_hat_lo)
  
  # 5) Create data frames for plotting
  plot_data_full <- data.frame(observed = y_obs, predicted = y_hat_full)
  plot_data_hi <- data.frame(observed = y_obs, predicted = y_hat_hi)
  plot_data_lo <- data.frame(observed = y_obs, predicted = y_hat_lo)
  
  # 6) Create observed vs predicted plots for each partition with titles
  # Create observed vs predicted plots for each partition with improved titles
  plot_full <- ggplot(plot_data_full, aes(x = observed, y = predicted)) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(title = bquote(atop(.(t), atop("FULL", 
                                        "RMSEP:" ~ .(stats_full$rmsep) ~ ", R2:" ~ .(stats_full$r2) ~ ", comp:" ~ .(P_DIM)))),
         x = "Observed Values", y = "Predicted Values") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, hjust = 0.5)) +  # Adjust font size and center alignment
    coord_fixed(ratio = 1) +
    xlim(range(plot_data_full$observed)) +
    ylim(range(plot_data_full$observed))
  
  plot_hi <- ggplot(plot_data_hi, aes(x = observed, y = predicted)) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(title = bquote(atop(.(t), atop("HI", 
                                        "RMSEP:" ~ .(stats_hi$rmsep) ~ ", R2:" ~ .(stats_hi$r2) ~ ", comp:" ~ .(P_DIM)))),
         x = "Observed Values", y = "Predicted Values") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, hjust = 0.5)) +  # Adjust font size and center alignment
    coord_fixed(ratio = 1) +
    xlim(range(plot_data_hi$observed)) +
    ylim(range(plot_data_hi$observed))
  
  plot_lo <- ggplot(plot_data_lo, aes(x = observed, y = predicted)) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(title = bquote(atop(.(t), atop("LO", 
                                        "RMSEP:" ~ .(stats_lo$rmsep) ~ ", R2:" ~ .(stats_lo$r2) ~ ", comp:" ~ .(P_DIM)))),
         x = "Observed Values", y = "Predicted Values") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, hjust = 0.5)) +  # Adjust font size and center alignment
    coord_fixed(ratio = 1) +
    xlim(range(plot_data_lo$observed)) +
    ylim(range(plot_data_lo$observed))
  
  # 7) Store all plots in a list for the current trait
  p_list[[t]] <- list(plot_full, plot_hi, plot_lo)
}

# 8) Arrange all plots into a 3x3 grid and display
p_combined <- gridExtra::grid.arrange(
  p_list$LMA[[1]], p_list$LMA[[2]], p_list$LMA[[3]],
  p_list$LDMC[[1]], p_list$LDMC[[2]], p_list$LDMC[[3]],
  p_list$EWT[[1]], p_list$EWT[[2]], p_list$EWT[[3]],
  nrow = 3, ncol = 3
)

print(p_combined)