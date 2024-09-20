# Objective: merge all csv files together from raw data
# to have a single csv file of spectroscopy measurements
# Reference: https://github.com/MichaletzLab/Excision-Photosynthesis-Methods/blob/main/Read.and.Modify.Data.R
# Author: nathan malamud
# date: 2024.05.24
#

# load libraries
library(dplyr)
library(purrr)
library(readxl)
library(tidyverse)

# WARNING: this library is a deprecated build
# source documentation: https://github.com/meireles/spectrolab
library(spectrolab)

# Load spectroscopy files.
# Remove metadata file from paths variable
spec.files <- list.files(
  path = "./raw",
  full.names = TRUE)
spec.files <- spec.files[spec.files != "./raw/metadata.xlsx"]

processed_data <- lapply(spec.files, function(file) {
  # read data and import as matrix
  spec <- read_spectra(file, format="sig")
  x <- as.matrix(x=spec) %>% as.data.frame()
  
  # extract barcodeID with regex and append to matrix
  pattern <- "(?<prelim>[0-9]{8})_(?<barcodeID>.*?)_[0-9]{4}\\.sig"
  match <- str_match(rownames(x)[1], pattern) %>% as.data.frame()
  x$barcodeID <- match$barcodeID[1] %>% as.character()
  
  return(x)
} )

# Check the column names of all data frames in processed_data
# Find the common column names across all data frames
# Filter the data frames to keep only the common columns
all_col_names <- lapply(processed_data, colnames)
common_col_names <- Reduce(intersect, all_col_names)
processed_data_filtered <- lapply(processed_data, function(df) df[, common_col_names, drop = FALSE])

# Make into a single data frame called full.exin.data
full.exin.data <- do.call(rbind, processed_data_filtered) %>%
  select(barcodeID, everything())

# write to current directory
write_csv(full.exin.data, "./spec_data.csv")
