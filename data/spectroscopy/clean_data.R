##' Objective: merge all csv files together from raw data
##' to have a single csv file of spectroscopy measurements
##'
##' Adapted from code written for 7th
##' annual plant functional traits course.
##' 
##' @author [Nicole Bison, Nathan Malamud]
##' @date 2024.10.02
##' 
##' Source:
##'   https://github.com/MichaletzLab/pftc7_spectroscopy
##'   

# - - - - -
# Libraries
library(readxl)
library(tidyverse) 
library(dplyr)
library(readr)
library(lubridate)

# (!!!) This script depends on the "spectrolab" library.
#       This library is a deprecated build. Can't install with CRAN.
#   Source: https://github.com/meireles/spectrolab
#
# Download using the following code:
#   devtools::install_github("meireles/spectrolab")
#
library(spectrolab)

# Define useful dplyr "macros" here
`%notin%` <- Negate(`%in%`)

# - - - - -
# Find the files to read
# IMPORTANT: set working directory > source file location
# Ensure working directory contains folder with the raw data files.
INPUT_DIR <- "./raw_spec" 
OUTPUT_DIR <- "../"

if (!file.exists(INPUT_DIR)) {
  print("Please set working directory to source file location in RStudio.")
  stop("Script terminated.")
}


# - - - - -
# Load spectroscopy files.
all.files = list.files(path=INPUT_DIR,
                      full.names = T,
                      recursive = T,
                      pattern = 'sig')

# Remove metadata file from paths variable
all.files <- all.files[
  all.files != paste(INPUT_DIR, "/", "/metadata.xlsx")
  ]

# Loop through, read file, correct sensor overlap, return as list of data.frames
spectra_list <- lapply(all.files, function(file) {
     # Read data and import as matrix
     x <- read_spectra(file, format="sig")
    
     # (!!!) Nicole B. says: important to correct for sensor overlap
     # in 990-1900 range due to feature of SVC equipment.
     x2 <- match_sensors(x, c(990,1900)) %>% as.data.frame()
     
     # Convert x to data frame so filename is accessible
     x <- as.matrix(x) %>% as.data.frame()
     
     # Match filename of spectra to extract info
     pattern <- "(?<date>[0-9]{8})_(?<barcodeID>.*?)_(?<scan>[0-9]{4})\\.sig"
     match <- str_match(rownames(x)[1], pattern) %>% as.data.frame()
     
     # Extract barcodeID with Reg-Ex and append to matrix
     x2$date <- match$date[1] %>% as.character()
     x2$barcodeID <- match$barcodeID[1] %>% as.character()
     x2$scan <- match$scan[1] %>% as.character()
     
     return(x2)
   }
 )
 
# - - - - -
# Check the column names of all data frames in processed_data
# Find the common column names across all data frames
# Filter the data frames to keep only the common columns
all_col_names <- lapply(spectra_list, colnames)
common_col_names <- Reduce(intersect, all_col_names)
processed_data_filtered <- lapply(spectra_list, function(df) df[, common_col_names, drop = FALSE])
 
# Make into a single data frame called full.exin.data
all.spectra <- do.call(rbind, processed_data_filtered) %>%
  select(barcodeID, scan, date, everything())

# Remove white references and practice test scans
all.spectra$barcodeID <- toupper(all.spectra$barcodeID)
all.spectra <- all.spectra %>% filter(barcodeID %notin% c("WR", "WR_", "TEST", "SIG", "BARCODE"))

# Convert dates to standard format using lubridate
all.spectra$date <- ymd(all.spectra$date)

# (!!!) Metadata file in ./raw/metadata.xlsx
# says to delete id number 27826
all.spectra <- all.spectra %>% filter(barcodeID != "27826")

# - - - - - 
# Average the spectra for all scans for each barcode ID
all.spectra <- all.spectra %>%
  select(-c(scan, sample_name)) %>%              # Remove unnecessary columns
  group_by(barcodeID, date) %>%                  # Group by barcodeID and date
  summarize(across(everything(), mean)) %>%      # Calculate mean for all other columns
  ungroup()                                      # Ungroup to finalize the result

# Merge spectra with metadata csv
metadata <- read.csv("./metadata.csv")
all.spectra <- merge(metadata, all.spectra, by="barcodeID")

# Write to current directory
all.spectra.wide <- all.spectra
all.spectra.long <- all.spectra %>% pivot_longer(
  cols = -c(colnames(metadata), date),
  names_to = "wavelength",
  values_to = "reflectance"
)

write_csv(all.spectra.wide, paste(OUTPUT_DIR, "spec_data_wide.csv", sep=""))
write_csv(all.spectra.long, paste(OUTPUT_DIR, "spec_data_long.csv", sep=""))
