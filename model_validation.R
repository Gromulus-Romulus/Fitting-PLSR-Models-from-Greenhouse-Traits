##' Model cross validation - calculate PRESS statistics.
##'
##' @author Nathan Malamud
##' @date 2024.10.04
##' 

library(dplyr)
library(tidyverse)
library(prospectr)
library(caretr)
library(spectratrait)

spectra_wide <- read.csv("./data/spec_data_wide.csv") |>
  select(-date)

# Remove "X" characters from beginning of wavelengths
spec_names <- colnames(spectra_wide)[5:ncol(spectra_wide)]
colnames(spectra_wide)[5:ncol(spectra_wide)] <- gsub("X", "", spec_names)

# - - - - -
# TODO: use signal_cleaning.R to clean up the data
# Export derived data to a new file

# TODO: import signal_output.csv, merge with traits data
# FIT PLSR model to assigned calibration vs validation datasets
# Use Sean Serbin's package for fitting PLSR models

# TODO: calculate PRESS statistics ... evaluate model transferability
# For multiple traits (e.g. LMA, LDMC)