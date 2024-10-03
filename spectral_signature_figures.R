##' Objective: Plot differences in spectral signatures
##' across low and high treatments
##'
##' Adapted from code written for 7th
##' annual plant functional traits course.
##' 
##' TODO: Do a wavelength-by-wavelength T test
##' and add p values to each chart
##'   Reference: https://journals.ashs.org/hortsci/view/journals/hortsci/53/5/article-p669.xml
##' 
##' @author [Nicole Bison, Nathan Malamud]
##' @date 2024.10.02
##' 
##' Source:
##'   https://github.com/MichaletzLab/pftc7_spectroscopy
##'   

library(dplyr)
library(tidyverse)
library(rcartocolor)
library(MetBrewer)

# - - - - - 
# Read in input files and remove "date" column
# as it is irrelevant for now.
spectra_long <- read.csv("./data/spec_data_long.csv") |>
  select(-date)

spectra_wide <- read.csv("./data/spec_data_wide.csv") |>
  select(-date)

# Remove "X" characters from beginning of wavelengths
spec_names <- colnames(spectra_wide)[5:ncol(spectra_wide)]
colnames(spectra_wide)[5:ncol(spectra_wide)] <- gsub("X", "", spec_names)

# - - - - - 
# Categorize treatments into "LO" and "HI" based on set threshold for mmol
THRESHOLD <- 15
spectra_long <- spectra_long %>%
  mutate(treatment_category = ifelse(treatment_mmol <= THRESHOLD, "LO", "HI"))

# - - - - -
# Average measurements per-species and treatment
per_species_treatment <- spectra_long %>% group_by(species, treatment_category, wavelength) %>% 
  summarize(mean_r = mean(reflectance),
            sd_r = sd(reflectance),
            n = n(),
            se_r = sd_r / sqrt(n)) %>% ungroup()

# Assign wavelength ranges
# Don't include below < 400 and above 2400
# due to increased noise at range edges
# Ranges set using SVC i-series Field Spectroscopy guide
FULL = seq(400, 2400, by = .1)
VIS = seq(400, 700, by = .1)
NIR = seq(700, 1000, by = .1)
SWIR = seq(1000, 2400, by = .1)

spectra_long <- spectra_long %>%
  mutate(wave_type = case_when(
    wavelength %in% FULL ~ "Full",
    wavelength %in% VIS ~ "Visible",
    wavelength %in% NIR ~ "NIR",
    wavelength %in% SWIR ~ "SWIR"
  ))

# - - - - -
# Helper function for making spectral plots
spectra_plot <- function(data) {
  ggplot(data, aes(x = wavelength, color = treatment_category, fill = treatment_category)) +
    geom_ribbon(aes(ymin = 100*(mean_r - 2*se_r), ymax = 100*(mean_r + 2*se_r), fill = treatment_category), 
                alpha = 0.3, color = NA) +
    geom_line(aes(y = 100*mean_r), linewidth = 0.75) +  # Thicker lines for better visibility
    ylab("% Reflectance") + 
    xlab("Wavelength (nm)") +
    scale_fill_manual(name = "Treatment Category", 
                      values = c("LO" = "darkgrey", "HI" = "red"),
                      labels = c("LO" = "low (0 - 10 mmol)", "HI" = "high (15 - 35 mmol)")) +
    scale_color_manual(name = "Treatment Category", 
                       values = c("LO" = "darkgrey", "HI" = "red"),
                       labels = c("LO" = "low (0 - 10 mmol)", "HI" = "high (15 - 35 mmol)")) +
    theme_minimal() + 
    facet_wrap(~species) +
    theme(
      text = element_text(family = "Helvetica", face = "plain"),  # Regular font for all text
      axis.title = element_text(size = 14, face = "bold"),  # Larger axis titles
      axis.text = element_text(size = 12, face = "plain"),   # Larger axis text
      legend.title = element_text(size = 14, face = "bold"), # Larger legend title
      legend.text = element_text(size = 12, face = "bold"),  # Larger legend text
      legend.position = "bottom",
      legend.spacing.x = unit(0.5, 'cm'),                   # Increase spacing in the legend
      strip.text = element_text(size = 14, face = "bold.italic"),  # Larger facet labels
      panel.grid.major = element_line(size = 0.5, color = "lightgrey"),  # Subtle major gridlines
      panel.grid.minor = element_line(size = 0.3, color = "lightgrey")   # Subtle minor gridlines
    )
}

# Function to filter data by wavelength range and generate the spectral plot
OUTPUT_DIR <- "figures"
plot_by_range <- function(data, range_name, wavelength_range) {
  data_filtered <- data %>% filter(wavelength %in% wavelength_range)
  
  plot <- spectra_plot(data_filtered) +
    ggtitle(paste("Spectral Plot -", range_name))  # Adding a title for each plot
  
  # Save plot to a PDF file
  ggsave(filename = paste0(OUTPUT_DIR, "/", "spectral_plot_", range_name, ".pdf"), 
         plot = plot, device = "pdf", width = 10, height = 6)
}

# Iterate through the ranges and create plots
plot_by_range(per_species_treatment, "FULL", FULL)
plot_by_range(per_species_treatment, "VIS", VIS)
plot_by_range(per_species_treatment, "NIR", NIR)
plot_by_range(per_species_treatment, "SWIR", SWIR)

# - - - - -
# Plot of 2nd derivative transformation
# calculate absorbance
# TODO: smooth derivatives, review this paper:
#   Source: https://academic.oup.com/jxb/article/63/1/489/560215
spec_names <- colnames(spectra_wide)[6:ncol(spectra_wide)]

# Filter out spectral matrix
X <- spectra_wide %>% select(spec_names) %>%
  as.matrix()

d1 <- t(diff(t(X), differences = 1)) # first derivative
d2 <- t(diff(t(X), differences = 2)) # second derivative
plot(as.numeric(colnames(d1)), 
     d1[1,], 
     type = "l", 
     lwd = 1.5, 
     xlab = "Wavelength", 
     ylab = "")

lines(as.numeric(colnames(d2)), d2[1,], lwd = 1.5, col = "red")
grid()
legend("topleft", 
       legend = c("1st der", "2nd der"), 
       lty = c(1, 1),
       col = c("black", "red"))
