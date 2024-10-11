##' Model cross validation across traits, treatments, and species.
##' TODO: clean up documentation with GPT.
##'
##' @author Nathan Malamud
##' @date 2024.10.08
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
plsr_data <- cbind(barcodeID, Y, X_raw) %>%
  as.data.frame()

# Traits and species of interest
trait_ids <- c("LMA", "LDMC", "EWT")
units_lookup <- c("LMA" = "g/m^2", "LDMC" = "mg/g", "EWT" = "g/m^2") 
species_ids <- traits$species %>% unique()

# - - - - -
# Comparative statistics using cross-validation
# Calculate RMSEP, R2, PRESS for each model and cross-compare
THRESH = 15 # mmol
metadata$treatment_category <- ifelse(metadata$treatment_mmol > THRESH, "high", "low")

plsr_data_hi <- plsr_data %>%
  filter(metadata$treatment_category == "high") # 52 obs > 15 mmol

plsr_data_lo <- plsr_data %>%
  filter(metadata$treatment_category == "low") # 45 obs <= 15 mmol

# Define the traits of interest
trait_ids <- c("LMA", "LDMC", "EWT")
P_DIM <- 10

# Plot with color by species (Borage, Barley, Radish)
josef_colors <- c("#7570b2", "#ca621c", "#299680")

# Create a list to store plots for all traits
p_list <- list()

# TODO : export ggplot to PDF
# Loop over each trait and perform PLSR for "FULL", "LO", and "HI" data partitions
for (t in trait_ids) {
  
  # 1) Define observed and spectral matrices for each partition
  y_obs <- (plsr_data[[t]])
  y_obs_hi <- plsr_data_hi[[t]]
  y_obs_lo <- plsr_data_lo[[t]]
  
  x <- plsr_data[, spectral_columns] %>% as.matrix()
  x_hi <- plsr_data_hi[, spectral_columns] %>% as.matrix()
  x_lo <- plsr_data_lo[, spectral_columns] %>% as.matrix()
  
  # 2) Fit PLSR models for each partition
  plsr_model_full <- plsr(y_obs ~ x, ncomp = P_DIM, validation = "LOO")
  plsr_model_hi <- plsr(y_obs_hi ~ x_hi, ncomp = P_DIM, validation = "LOO")
  plsr_model_lo <- plsr(y_obs_lo ~ x_lo, ncomp = P_DIM, validation = "LOO")
  
  nComps_full <- selectNcomp(plsr_model_full, plot = F)
  
  # 3) Predict and extract fitted values
  y_hat_full <- predict(plsr_model_full, x, ncomp = nComps_full) %>% as.vector()
  y_hat_hi <- predict(plsr_model_hi, x, ncomp = nComps_full) %>% as.vector()
  y_hat_lo <- predict(plsr_model_lo, x, ncomp = nComps_full) %>% as.vector()
  
  # 4) Calculate RMSEP and R²
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
  plot_data_full <- data.frame(observed = y_obs, predicted = y_hat_full, Crop = metadata$species)
  plot_data_hi <- data.frame(observed = y_obs, predicted = y_hat_hi, Crop = metadata$species)
  plot_data_lo <- data.frame(observed = y_obs, predicted = y_hat_lo, Crop = metadata$species)
  
  # 6) Create observed vs predicted plots for each partition with species-colored points
  plot_full <- ggplot(plot_data_full, aes(x = observed, y = predicted, color = Crop)) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(title = bquote(atop(bold(.(paste(t, ": all data ( n =", nrow(plsr_data), ")"))),
                             atop("RMSEP:" ~ .(stats_full$rmsep) ~ ", R2:" ~ .(stats_full$r2) ~ ", comp:" ~ .(nComps_full)))),
         x = "Observed", y = "Predicted") +
    scale_color_manual(values = josef_colors) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, hjust = 0.5)) +
    coord_fixed(ratio = 1) +
    xlim(range(plot_data_full$observed)) +
    ylim(range(plot_data_full$observed))
  
  plot_hi <- ggplot(plot_data_hi, aes(x = observed, y = predicted, color = Crop)) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(title = bquote(atop(bold(.(paste(t, ": >", THRESH, "mmol ( n =", nrow(plsr_data_hi), ")"))),
                             atop("RMSEP:" ~ .(stats_hi$rmsep) ~ ", R2:" ~ .(stats_hi$r2) ~ ", comp:" ~ .(nComps_full)))),
         x = "Observed", y = "Predicted") +
    scale_color_manual(values = josef_colors) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, hjust = 0.5)) +
    coord_fixed(ratio = 1) +
    xlim(range(plot_data_hi$observed)) +
    ylim(range(plot_data_hi$observed))
  
  plot_lo <- ggplot(plot_data_lo, aes(x = observed, y = predicted, color = Crop)) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(title = bquote(atop(bold(.(paste(t, ": <=", THRESH, "mmol ( n =",  nrow(plsr_data_lo), ")"))),
                             atop("RMSEP:" ~ .(stats_lo$rmsep) ~ ", R2:" ~ .(stats_lo$r2) ~ ", comp:" ~ .(nComps_full)))),
         x = "Observed", y = "Predicted") +
    scale_color_manual(values = josef_colors) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, hjust = 0.5)) +
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
