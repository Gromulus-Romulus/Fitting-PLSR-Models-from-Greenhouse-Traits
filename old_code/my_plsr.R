# Proof of Concept for Fitting PLSR Model
# Predicting Traits from Spectra Data
##' @author Nathan Malamud
##' @date 2024.09.19
#
# TODO: why running script more than once causes fatal R error (abort?)

library(tidyverse)
library(readxl)
library(LeafArea)
library(prospectr)
library(pls)
library(gridExtra)

# Load traits and spectroscopy data
traits <- readRDS("./proc/traits.RDS")
spec <- readRDS("./proc/spec.RDS")

# Use prospectr for PLSR: Spectra -> Traits Factored by N Treatment
# prospectr is widely used in spectroscopic and remote-sensing applications
#   Source: https://cran.r-project.org/web/packages/prospectr/vignettes/prospectr.html

# Note: At this point, we have a matrix of trait vectors (Y) and a matrix of spectral vectors (X).
# The goal is to fit X into a parameterized PLSR model and compare predicted Y with measured Y.

# Important: Not all leaves measured have spec readings, and vice versa.
# We filter only the barcode IDs that exist in both the X and Y matrices.
barIDs <- intersect(spec$barcodeID, traits$barcodeID)

# Isolate Spectral Matrix (X)
X <- spec %>%
  filter(barcodeID %in% barIDs) %>%
  select(-c(barcodeID, sampleID, treatment_mmol, species)) %>%
  as.matrix()

# Burnett et al. 2021 reccomend filtering out wavelengths < 400 nm and > 2500 nm
# Due to increased sensor noise and lower signal-to-noise ratio

# Isolate Trait Matrix (Y)
Y <- traits %>%
  filter(barcodeID %in% barIDs) %>%
  select(LDMC, LMA) %>%
  as.matrix()

# Labels for N_mmol treatment
L <- traits %>%
  filter(barcodeID %in% barIDs) %>%
  select(treatment_mmol) %>%
  as.matrix()

# Step 1: Clean Signal Noise for Spectral Matrix
# Opted for computationally lighter method for calibration - SNV
# TODO: Using Multiplicative Scatter Correction (MSC)
snv_spc <- standardNormalVariate(X %>% scale())

# Step 2: Fit PLSR Model
# The number of components (ncomp) can be tuned based on cross-validation (CV)
plsr_model <- plsr(Y ~ snv_spc, ncomp = 5, validation = "CV")

# Assess model performance by plotting the validation results
validationplot(plsr_model, val.type = "MSEP")

# Step 4: Predictions
predicted_Y <- predict(plsr_model, newdata = snv_spc)

# Step 5: Compare Predicted Y with Measured Y (Optional)
# TODO: look into warnings here

# Combine measured and predicted values into a single dataframe
Measured = as.data.frame(Y)
Predicted = as.data.frame(predicted_Y)
colnames(Measured) <- c("Measured_LMA", "Measured_LDMC")
colnames(Predicted) <- c("Predicted_LMA", "Predicted_LDMC")

#  Create a new column indicating whether values are Measured or Predicted
# Combine measured and predicted into one data frame
comparison <- cbind(Measured, Predicted)
comparison$Nmmol <- as.vector(L) 

# Plot for LMA: Measured vs Predicted
plot_LMA <- ggplot(comparison, aes(x = Measured_LMA, y = Predicted_LMA, color = as.numeric(Nmmol))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Perfect prediction line
  scale_color_continuous(name = "N Treatment (mmol)") +  # Color by Nmmol
  labs(title = "Measured vs Predicted LMA",
       x = "Measured LMA",
       y = "Predicted LMA") +
  theme_classic() + xlim(75, 200) + ylim(75, 200)

# Plot for LDMC: Measured vs Predicted
plot_LDMC <- ggplot(comparison, aes(x = Measured_LDMC, y = Predicted_LDMC, color = as.numeric(Nmmol))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Perfect prediction line
  scale_color_continuous(name = "N Treatment (mmol)") +  # Color by Nmmol
  labs(title = "Measured vs Predicted LDMC",
       x = "Measured LDMC",
       y = "Predicted LDMC") +
  theme_classic() + xlim(0.03, 0.07) + ylim(0.03, 0.07)

# Display the two plots side by side
grid.arrange(plot_LMA, plot_LDMC, ncol = 2)

# Model Summary
summary(plsr_model)
