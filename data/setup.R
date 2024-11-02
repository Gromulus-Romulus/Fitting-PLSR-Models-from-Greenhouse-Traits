##' Objective: prepare raw traits data for analysis.
##'
##' Sections:
##'   * Read in leaf mass Excel sheet
##'   * Run imageJ software on leaf scans
##'   * Calculate LMA and LDMC from weight and area measurements
##'   * Load fluorometry data from porometer and rename columns
##'   * Write traits and fluorometry to CSV file
##'
##' Output from R script will be fed into plsr.R for model fitting and cross-validation.
##' Dependencies: requires ImageJ to be downloaded, may need devtools for R package installations.
##'
##' @author: Nathan Malamud
##' @date: 2024.09.20

# Pre-processing packages
library(tidyverse)
library(readxl)
library(reshape2)
library(LeafArea)

# REMINDER: Set Working Directory -> Source File Location
mass_data <- read_excel("./leaf_mass_data.xlsx")

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
scans <- run.ij(set.directory = './scans/raw')
names(scans) <- c("barcodeID", "area_cm2")
mass_data <- merge(mass_data, scans, by = "barcodeID")

# Calculate LMA from Weight and Area Measurements (units = g / m2)
mass_data$LMA <- ((mass_data$dry_leaf_g) / (mass_data$area_cm2)) * (100^2)

# Calculate EWT from weight and Area Measurements (units = g / m2)
mass_data$EWT <- ((mass_data$wet_leaf_g - mass_data$dry_leaf_g) / (mass_data$area_cm2)) * (100^2)

# Load Fluorometry Data from Porometer and Rename Columns
# Use Fs and Fm' to Calculate Quantum Yield of Fluorescence
# Merge Porometer Measurements with Data Afterwards
#   Source: https://www.licor.com/env/products/LI-600/
fluor <- read_csv("./fluor/fluor_data.csv") %>%
  select(c(unique_id, gsw, Fs, `Fm'`, PhiPS2, ETR, Qamb, VPleaf))
names(fluor) <- c("barcodeID", "gsw", "Fs", "Fm_prime", "Phi_PS2", "ETR", "Qamb", "vpdl")

# Write traits (only those of interest) data to R data file.
# Also merge measured traits with LI-COR fluorometry data
traits <- mass_data |>
  select("barcodeID",
         "sampleID",
         "treatment_mmol",
         "species",
         "dry_whole_g", "dry_leaf_g", "area_cm2",
         "LDMC", "LMA", "EWT") %>%
  merge(fluor, by="barcodeID")

# Remove duplicate values
traits <- traits %>%
  distinct(barcodeID, .keep_all = TRUE)

# Create a lookup table for species to Latin names
latin.names <- c("borage" = "B. officinalis",
                 "barley" = "H. vulgare",
                 "radish" = "R. sativus")

# Rename species column using the lookup table
traits <- traits %>%
  mutate(species = recode(species,
                          "borage" = latin.names["borage"],
                          "barley" = latin.names["barley"],
                          "radish" = latin.names["radish"]))

write_csv(traits, file = ("./traits.csv"))
