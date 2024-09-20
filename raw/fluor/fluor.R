# Objective: merge all csv files together from raw data
# to have a single csv file of fluorometry measurements
# Reference: https://github.com/MichaletzLab/Excision-Photosynthesis-Methods/blob/main/Read.and.Modify.Data.R
# Author: nathan malamud
# date: 2024.05.16
#

# load libraries
library(dplyr)
library(purrr)
library(readxl)
library(tidyverse)

# REMEMBER: set working directory -> source file location
# Get the files in the directory
files <- list.files(path = "./to_merge", full.names = TRUE)

processed_data <- lapply(files, function(file) {
  # read in the file and skip the first header
  data <- read_csv(file, skip = 1)
  
  # remove the first row containing units
  data <- data[-1, ]
  
  # Return the processed data frame
  return(data)
})

# Check the column names of all data frames in processed_data
all_col_names <- lapply(processed_data, colnames)

# Find the common column names across all data frames
common_col_names <- Reduce(intersect, all_col_names)

# Filter the data frames to keep only the common columns
processed_data_filtered <- lapply(processed_data, function(df) df[, common_col_names, drop = FALSE])

# Make into a single data frame called full.exin.data
full.exin.data <- do.call(rbind, processed_data_filtered)

# write to current directory
write_csv(full.exin.data, "./fluor_data.csv")
