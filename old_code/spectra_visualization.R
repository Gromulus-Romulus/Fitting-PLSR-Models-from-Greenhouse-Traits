# Show differences in spectra between species and treatments
#
# Documentation guide (Serbin et al. 2022): 
#   https://github.com/plantphys/spectratrait/blob/main/spectratrait_1.2.5.pdf
#
##' @author Nathan Malamud
##' @date 2024.09.23

##' @author Nathan malamud
##' @date 2024.09.23

# clear namespace
rm(list=ls())

# load dependencies
library(dplyr)
library(ggplot2)
library(spectratrait)
library(tidyverse)

# Define useful dplyr "macros"
`%notin%` <- Negate(`%in%`)

# Script options
output_dir <- "tempdir"

Start.wave <- 500
End.wave <- 2400
wv <- seq(Start.wave, End.wave, .1)

# Load traits and spectroscopy data
traits <- readRDS("./proc/traits.RDS")
spec <- readRDS("./proc/spec.RDS")

# Filter only samples with both spec and trait data
barIDs <- intersect(spec$barcodeID, traits$barcodeID)
traits <- traits %>% filter(barcodeID %in% barIDs)
spec <- spec %>% filter(barcodeID %in% barIDs)

# Remove NA rows from spec dataframe
spec <- spec %>% filter(complete.cases(.))

# Filter spectra data frame to get wavelengths in specified range
spec <- spec |>
  pivot_longer(cols = c(5:1028), names_to="wavelength", values_to="reflectance") %>%
  filter(as.numeric(wavelength) >= Start.wave & as.numeric(wavelength) <= End.wave) %>%
  pivot_wider(names_from="wavelength", values_from="reflectance")

# Remove duplicate items
spec <- spec %>% distinct(barcodeID, .keep_all = TRUE) %>% distinct(sampleID, .keep_all = TRUE)

# Update wavelength vector
wv <- as.numeric(colnames(spec[5:length(colnames(spec))]))

# Rename columns of spectra dataframe
names(spec)[5:length(names(spec))] <- paste0("Wave_", names(spec)[5:length(names(spec))])

# Focus on relevant traits
traits <- traits %>% select(barcodeID, sampleID, species, treatment_mmol, LDMC, LMA, Phi_PS2, Fm_prime, dry_whole_g)

# Merge traits and spectral data
plsr_data <- merge(traits, spec, by=c("barcodeID", "sampleID", "species", "treatment_mmol"))

# Create a lookup table for species to Latin names
latin_names <- c("borage" = "B. officinalis",
                 "barley" = "H. vulgare",
                 "radish" = "R. sativus")

# Rename species column using the lookup table
plsr_data <- plsr_data %>%
  mutate(species = recode(species,
                          "borage" = latin_names["borage"],
                          "barley" = latin_names["barley"],
                          "radish" = latin_names["radish"]))

# Filter dataset for LO and HI treatments
# Categorize treatments based on numeric values
plsr_data <- plsr_data %>%
  mutate(treatment_category = ifelse(treatment_mmol %in% c(0, 5, 10),
                                     "0-10 mmol (low)", ">15 mmol (high)"))

lo_hi_data <- plsr_data %>%
  filter(treatment_category %in% c("0-10 mmol (low)", ">15 mmol (high)"))

# Pivot data for plotting
long_spectra <- lo_hi_data %>%
  pivot_longer(cols = starts_with("Wave_"), names_to = "wavelength", values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(gsub("Wave_", "", wavelength)))

# Plot spectral differences
# Plot spectral differences with italic species names and regular face for everything else
ggplot(long_spectra, aes(x = wavelength, y = reflectance, color = treatment_category)) +
  geom_line(aes(group = interaction(barcodeID, species)), alpha = 0.4) +
  facet_wrap(~ species, scales = "free") +
  labs(x = "Wavelength (nm)", y = "Prop. Reflectance", color = "Treatment") +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", face = "plain"),  # Regular font for all text
    axis.title = element_text(size = 12, face = "bold"),  # Axis titles in regular
    axis.text = element_text(size = 10, face = "plain"),   # Axis text in regular
    legend.title = element_text(size = 12, face = "bold"),# Legend title in regular
    legend.text = element_text(size = 10, face = "bold"), # Legend text in regular
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold.italic")  # Species names in italic in facet labels
  )
