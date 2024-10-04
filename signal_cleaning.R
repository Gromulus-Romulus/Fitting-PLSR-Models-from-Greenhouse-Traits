##' Clean up spectroscopy data
##' prior to model fitting.
##'
##' @author Nathan Malamud
##' @date 2024.10.04
##' 

library(dplyr)
library(tidyverse)
library(prospectr)

spectra_wide <- read.csv("./data/spec_data_wide.csv") |>
  select(-date)

# Remove "X" characters from beginning of wavelengths
spec_names <- colnames(spectra_wide)[5:ncol(spectra_wide)]
colnames(spectra_wide)[5:ncol(spectra_wide)] <- gsub("X", "", spec_names)

# - - - - -
# Plot of 2nd derivative transformation
# calculate absorbance
# https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter?useskin=vector
# TODO: smooth derivatives, review this paper:
#   Source: https://academic.oup.com/jxb/article/63/1/489/560215
# m = order of the derivative
# w = gap size
# s = segment size
# first derivative with a gap of 5 bands
library(prospectr)

# Apply Savitzky-Golay smoothing to the spectral matrix
sg_smooth <- savitzkyGolay(X = X, p = 3, w = 11, m = 0)

# Apply the first derivative
sg_first_deriv <- savitzkyGolay(X = X, p = 3, w = 11, m = 1)

# Apply the second derivative
sg_second_deriv <- savitzkyGolay(X = X, p = 3, w = 11, m = 2)

# Plot the first and second derivatives for the first observation
plot(as.numeric(colnames(sg_first_deriv)), 
     sg_first_deriv[1, ], 
     type = "l", 
     lwd = 1.5, 
     xlab = "Wavelength", 
     ylab = "")

lines(as.numeric(colnames(sg_second_deriv)), sg_second_deriv[1, ], lwd = 1.5, col = "red")
grid()
legend("topleft", 
       legend = c("1st der", "2nd der"), 
       lty = c(1, 1),
       col = c("black", "red"))

# TODO: rewrite this plot using ggplot2
# TODO: output pdfs to figures/ directory

# TODO: see if Savitzky-Golay Smoothing + Normalization reduces RMSEP in model_validation.R


