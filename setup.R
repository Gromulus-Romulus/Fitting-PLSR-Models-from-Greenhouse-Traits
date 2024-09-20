# Objective: prepare raw data for PLSR analysis.
#
# Sections:
#   * Read in leaf mass Excel sheet
#   * Run imageJ software on leaf scans
#   * Calculate LMA and LDMC from weight and area measurements
#   * Load fluorometry data from porometer and rename columns
#   * Load spectroscopy data from SVC and aggregate data
#   * Write traits and spectral data to separate R data files
#
# Output from R script will be fed into plsr.R for model fitting and cross-validation.
# Dependencies: requires imageJ to be downloaded, may need devtools for R package installations.
#
# Author: Nathan Malamud
# Date: 2024.09.20

# Pre-processing packages
library(tidyverse)
library(readxl)
library(reshape2)
library(LeafArea)
library(ggplot2)
library(ggpubr)

# REMINDER: Set Working Directory -> Source File Location
mass_data <- read_excel("./raw/leaf_mass_data.xlsx")

# Define Factor Levels (treatment, species, barcode)
mass_data$treatment_mmol <- as.factor(mass_data$treatment_mmol)
mass_data$species <- as.factor(mass_data$species)
mass_data$barcodeID <- as.factor(mass_data$barcodeID)
levels(mass_data$species) <- c("radish", "borage", "barley")
mass_data$LDMC <- as.numeric(mass_data$LDMC)

# Run ImageJ Software on Leaf Scans
# This line of code initiates a pop-up window.
# It will take a while to run, that's expected.
# Afterwards, merge values for image scans with the dataframe.
scans <- run.ij(set.directory = './raw/scans')
names(scans) <- c("barcodeID", "area_cm2")
mass_data <- merge(mass_data, scans, by = "barcodeID")

# Calculate LMA from Weight and Area Measurements (units = g / m2)
mass_data$LMA <- ((mass_data$dry_leaf_g) / (mass_data$area_cm2)) * 100

# Load Fluorometry Data from Porometer and Rename Columns
# Use Fs and Fm' to Calculate Quantum Yield of Fluorescence
# Merge Porometer Measurements with Data Afterwards
#   Source: https://www.licor.com/env/products/LI-600/
fluor <- read_csv("./raw/fluor/fluor_data.csv") %>%
  select(c(unique_id, gsw, Fs, `Fm'`, PhiPS2))
names(fluor) <- c("barcodeID", "gsw", "Fs", "Fm_prime", "Phi_PS2")

# Write traits (only those of interest) data to R data file.
# Also merge measured traits with LI-COR fluorometry data
traits <- mass_data |>
  select("barcodeID",
         "sampleID",
         "treatment_mmol",
         "species",
         "dry_whole_g",
         "LDMC", "LMA") %>%
  merge(fluor, by="barcodeID")

# Remove duplicate values
traits <- traits %>%
  distinct(barcodeID, .keep_all = TRUE)
  
saveRDS(traits, file = ("./proc/traits.RDS"))

# Load Spectroscopy Data from Porometer and Aggregate Data
# N Bison Mentioned There Are 2-3 Measurements per Individual
# TODO: Investigate 50+ Warnings Here
spec <- read_csv("./spec/spec_data.csv") %>%
  aggregate(by = list(.$barcodeID), FUN = mean)
names(spec) <- c("barcodeID", "delete_me", names(spec[3:1026]))
spec$barcodeID <- as.factor(spec$barcodeID)
spec <- spec %>% select(-delete_me)

# Write spec data to R data file, join with "barcodeID", "sampleID", "treatment_mmol"
# Filter out junk columns and remove problematic Barley sample (NA values for ID 38450)
spec <- left_join(traits %>% select("barcodeID", "sampleID", "treatment_mmol", "species"), spec, by="barcodeID") %>%
  filter(barcodeID != "38450")
saveRDS(spec, file = ("./proc/spec.RDS"))
