##' Objective: apply numerical techniques to clean and preprocess spectral data.
##' Reduces SNR (signal-to-noise) ratio using rectangular moving-average and Savitzky-Golay filters.
##' 
##' Written with Chat-GPT4
##'   Source: https://chatgpt.com/share/6705888f-c220-8012-88c2-0f4aa09fc319
##'   
##' @author [Nathan Malamud]
##' @date [2024.10.08]
##' 

library(dplyr)
library(tidyverse)
library(prospectr)
library(zoo)
library(ggplot2)

spectra_wide <- read.csv("./data/spec_data_wide.csv") |>
  select(-date)

# Clean the column names to remove 'X' prefix after loading csv
spec_names <- colnames(spectra_wide)[5:ncol(spectra_wide)]
colnames(spectra_wide)[5:ncol(spectra_wide)] <- gsub("X", "", spec_names)

# Extract column names registered from SVC spectrophotometer
spectral_columns <- as.numeric(gsub("X", "", colnames(spectra_wide)[5:ncol(spectra_wide)]))

# - - - - -
# Define helper functions - these are moving window filters
# that smooth spectral data at an interval of w discrete bands.

# Function to apply a moving average filter using rollapply
movav <- function(X, w) {
  t(apply(X, 1, function(row) zoo::rollapply(row, width = w, FUN = mean, fill = NA, align = "center")))
}

# Function to apply Savitzky-Golay polynomial smooth using prospectr functions
polyav <- function(X, w) {
  # Clean the column names to remove 'X' prefix
  spectral_columns <- colnames(X)
  
  # Apply Savitzky-Golay smoothing and generate first and second derivatives
  sg_first_deriv <- savitzkyGolay(X = X, p = 3, w = w, m = 1) # first derivative
  sg_second_deriv <- savitzkyGolay(X = X, p = 3, w = w, m = 2) # second derivative
  
  # Initialize matrices to hold the derivatives with the same dimensions as X, filled with NA
  X_d1 <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
  X_d2 <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
  
  # Calculate the number of bands lost due to the window size
  n_bands_lost <- floor(w / 2)
  
  # Fill X_d1 and X_d2 with the corresponding values from sg_first_deriv and sg_second_deriv
  X_d1[, (n_bands_lost + 1):(ncol(X) - n_bands_lost)] <- sg_first_deriv
  X_d2[, (n_bands_lost + 1):(ncol(X) - n_bands_lost)] <- sg_second_deriv
  
  # Set column names to match the original X
  colnames(X_d1) <- spectral_columns
  colnames(X_d2) <- spectral_columns
  
  # Double-triple check they are matrix types
  X_d1 <- as.matrix(X_d1)
  X_d2 <- as.matrix(X_d2)
  
  # Return the first and second derivative matrices
  return(list(X_d1 = X_d1, X_d2 = X_d2))
}


# - - - - -
# Set filter size for moving average
FILTER <- 11 

# Assign wavelength ranges
# Don't include below < 400 and above 2400
# due to increased noise at range edges
# Ranges set using SVC i-series Field Spectroscopy guide
#   Note: SVC i-series manual says full range goes to 2400,
#   Nicole says common practice is to cut off at 2500
FULL_R = seq(400, 2500, by = .1)
VIS_R = seq(400, 700, by = .1)
NIR_R = seq(700, 1000, by = .1)
SWIR_R = seq(1000, 2400, by = .1)

# Clean the column names to remove 'X' prefix
spectral_columns <- as.numeric(gsub("X", "", colnames(spectra_wide)[5:ncol(spectra_wide)]))

# Extract spectral data and convert to matrix
X <- as.matrix(spectra_wide[, 5:ncol(spectra_wide)])

# Apply moving average filters (rectangular and polynomial convolutions)
X_avg <- movav(X = X, w = FILTER)
sg <- polyav(X = X, w = FILTER)

X_d1 <- sg$X_d1
X_d2 <- sg$X_d2

# Trim the spectral matrices (X, X_filtered, sg_first_deriv, sg_second_deriv) to the FULL range
# by identifying which spectral_columns are within the FULL range
trim_range <- spectral_columns %in% FULL_R

# Trim the spectral_columns and matrices
spectral_columns_trimmed <- spectral_columns[trim_range]
X <- X[, trim_range]
X_avg <- X_avg[, trim_range]
X_d1 <- X_d1[, trim_range]
X_d2 <- X_d2[, trim_range]

# - - - - -
# Unique barcode IDs
uniq_ids <- unique(spectra_wide$barcodeID)

# Extract species and treatments for all unique IDs
species <- spectra_wide$species[match(uniq_ids, spectra_wide$barcodeID)]
treatments <- spectra_wide$treatment[match(uniq_ids, spectra_wide$barcodeID)]

# Open PDF for output to ./Figures directory
pdf("./Figures/SNR_spectra_plots.pdf")  # Change file name as needed

# Loop through each unique barcode ID
for (i in 1:length(uniq_ids)) {  # Correct iteration
  # Extract metadata
  id <- uniq_ids[i]
  s <- species[i]  # Correct indexing
  t <- treatments[i]  # Correct indexing
  
  # Extract reflectance vectors
  x <- X[i, ]
  x_avg <- X_avg[i, ] 
  x_d1 <- X_d1[i, ] 
  x_d2 <- X_d2[i, ]
  
  # Combine raw, moving average, and derivative data into a data frame for ggplot
  df <- data.frame(
    wavelength = rep(spectral_columns_trimmed, 4),  # Replicate trimmed wavelengths
    reflectance = c(x, x_avg, x_d1, x_d2),  # Reflectance values
    signal_type = rep(c("raw", "avg", "1d", "2d"), 
                   each = length(spectral_columns_trimmed))  # Signal processing type
  )
  
  # Plot using ggplot2 with facet_wrap allowing free y-axis scales
  p <- ggplot(df, aes(x = wavelength, y = reflectance, color = signal_type)) +
    geom_line(size = 1.5) +
    labs(x = "lambda", y = "% R", color = "Signal Processing Type") +
    theme_minimal() +
    scale_color_manual(values = c("black", "grey", "grey", "black")) +
    ggtitle(paste(id, ": ", s, " | ", t, "mmol")) +  # Correct title formatting
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~signal_type, labeller = label_both, scales = "free")  # Allow free y-axis ranges
  
  # Print the plot to the PDF
  print(p)  # Important to print the plot
  
}

# Close the PDF
dev.off()

# Save spectral matrices to RDS so we can work with them for PLSR fitting
signal.matrices <- list( 
  uniq_ids = uniq_ids,
  spectral_columns = spectral_columns,
  X_raw = X,
  X_avg = X_avg,
  X_d1 = X_d1,
  X_d2 = X_d2
)

saveRDS(signal.matrices, file = "./data/signal.matrices.RDS")
