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
  select(-c(barcodeID, sampleID, treatment_mmol)) %>%
  as.matrix()

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
plsr_model <- plsr(Y ~ snv_spc, ncomps = 10, validation = "CV")

# Assess model performance by plotting the validation results
validationplot(plsr_model, val.type = "MSEP")

# Step 4: Predictions
predicted_Y <- predict(plsr_model, newdata = snv_spc)

# Step 5: Compare Predicted Y with Measured Y (Optional)

# TODO: look into warnings here
merged_data <- spec %>%
  filter(barcodeID %in% barIDs) %>%
  select(barcodeID, sampleID, treatment_mmol) %>%
  inner_join(traits %>% filter(barcodeID %in% barIDs), by = "barcodeID")

comparison <- data.frame(
  Measured = as.vector(Y),
  Predicted = as.vector(predicted_Y),
  Nmmol = as.vector(L)
)

# Plot the Comparison
ggplot(comparison, aes(x = Measured, y = Predicted, color=as.numeric(Nmmol))) +
  geom_point() +
  geom_abline(color="red", slope = 1, intercept = 0) +
  labs(title = "Measured vs Predicted Traits",
       x = "Measured Y",
       y = "Predicted Y", ) + 
  scale_color_continuous(name="N Treatment (mmol)") +
  theme_classic()

# Model Summary
summary(plsr_model)
