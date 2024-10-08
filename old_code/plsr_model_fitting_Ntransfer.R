#' @author Nathan Malamud
#' @date 2024.09.23

# Clear namespace
rm(list=ls())

# Load dependencies
library(pls)
library(dplyr)
library(tidyverse)

# Define useful dplyr "macros"
`%notin%` <- Negate(`%in%`)

# PLSR options
pls.options(plsralg = "oscorespls")

# Load data
traits <- readRDS("./proc/traits.RDS")
spec <- readRDS("./proc/spec.RDS")

# Filter matching barcode IDs
barIDs <- intersect(traits$barcodeID, spec$barcodeID)
traits <- traits %>% filter(barcodeID %in% barIDs)
spec <- spec %>% filter(barcodeID %in% barIDs)

# Select and filter wavelengths
wv <- seq(500, 2400, 0.1)
spec_long <- spec %>%
  pivot_longer(cols = c(5:1028), names_to = "wavelength", values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),  reflectance = as.numeric(reflectance))

spec <- spec_long |>
  filter(as.numeric(wavelength) %in% wv) %>%
  pivot_wider(names_from = "wavelength", values_from = "reflectance")

# Merge traits and spec data
plsr_data <- left_join(traits, spec, by = c("barcodeID", "sampleID", "species", "treatment_mmol")) %>%
  select(barcodeID, sampleID, LMA, treatment_mmol)

# Categorize treatments
plsr_data <- plsr_data %>% mutate(treatment_category = ifelse(treatment_mmol %notin% c("15", "20", "25", "30", "35"), "Low", "High"))

# Split data into calibration and validation sets
cal.plsr.data <- filter(plsr_data, treatment_category == "Low")
val.plsr.data <- filter(plsr_data, treatment_category == "High")

# Fit PLSR model
plsr_model <- plsr(LMA ~ ., data = cal.plsr.data, ncomp = 10, validation = "CV")

# Predict and calculate RMSE
predictions <- predict(plsr_model, newdata = val.plsr.data[, -c(1:4)], ncomp = 10)
rmse <- sqrt(mean((val.plsr.data$LMA - predictions)^2))

cat("RMSE:", rmse, "\n")
