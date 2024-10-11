##' Objective: Fit prospect model to "low" and "high" nitrate treatments.
##' 
##' @author Nathan Malamud
##' @date 2024.09.27
##' 
##' References:
##'   https://jbferet.gitlab.io/prospect/index.html
##'   https://www-sciencedirect-com.ezproxy.library.ubc.ca/science/article/pii/003442579090100Z

library(tidyverse)
library(prospect)
library(dplyr)

# - - - - -
# Filter out wavelengths outside of 500-2400 nm
FULL <- seq(500.0, 2400.0, .01)

# Assume relationship between leaf reflectance and absorptance
spectra_long <- read.csv("./data/spec_data_long.csv") |>
  select(-date) %>% filter(wavelength %in% FULL)

# Aggregate lambda values via rounding -> mean reflectance
spectra_long$wavelength <- round(spectra_long$wavelength, 0)
spectra_long_agg <- aggregate(
  reflectance ~ barcodeID + species + treatment_mmol + wavelength,
  data = spectra_long, FUN=mean, na.rm = T
)

# - - - 
# Define parameters to estimate for PROSPECT-PRO engine
#   Source: https://jbferet.gitlab.io/prospect/articles/prospect.html
Parms2Estimate <- c('CHL', 'CAR', 'ANT', 'EWT', 'PROT', 'CBC')

# Create an empty list to store results for each leaf
results_list <- list()

# Get all unique barcodeIDs
unique_barcodes <- unique(spectra_long_agg$barcodeID)

# Loop over each unique barcodeID
for (barcode in unique_barcodes) {
  
  # Subset data for the current leaf
  leaf_spec <- spectra_long_agg %>% filter(barcodeID == barcode)
  
  # Adjust spectral domain
  SubData <- FitSpectralData(lambda = leaf_spec$wavelength,
                             Refl = leaf_spec$reflectance, 
                             Tran = NULL)
  
  # Invert PROSPECT to estimate parameters
  OutPROSPECTPRO <- Invert_PROSPECT(SpecPROSPECT = SubData$SpecPROSPECT, 
                                    Refl = SubData$Refl, 
                                    Tran = NULL,
                                    Parms2Estimate = Parms2Estimate, 
                                    PROSPECT_version = 'PRO')
  
  # Extract estimated parameters and store them along with the barcodeID
  estimated_params <- data.frame(
    barcodeID = barcode,
    CHL = OutPROSPECTPRO$CHL,
    CAR = OutPROSPECTPRO$CAR,
    ANT = OutPROSPECTPRO$ANT,
    EWT = OutPROSPECTPRO$EWT,
    PROT = OutPROSPECTPRO$PROT,
    CBC = OutPROSPECTPRO$CBC
  )
  
  # Append results to the list
  results_list[[as.character(barcode)]] <- estimated_params
}

# Combine all results into a single dataframe
final_results <- do.call(rbind, results_list)

# Merge with original data
final_results <- merge(final_results, spectra_wide %>%
                         select(barcodeID, sampleID, species, treatment_mmol), by = "barcodeID")

# - - -
# Create graph of how molecular content varies with nitrate treatment
# and across species
# Define colors for each species
species_colors <- c("R. sativus" = "#299680",      # Green
                    "B. officinalis" = "#7570b2", # Blue
                    "H. vulgare" = "#ca621c")     # Orange

final_results %>%
  pivot_longer(cols = c(CHL, CAR, ANT, EWT, PROT, CBC), names_to = "molecule", values_to = "value") %>%
  ggplot(aes(x = as.factor(treatment_mmol), y = value, color = species)) +
  geom_boxplot() +
  facet_wrap(~molecule, scales = "free_y") +
  scale_color_manual(values = species_colors) +  # Use the defined species colors
  labs(title = "Molecular content vs. nitrate treatment",
       x = "Nitrate treatment (mmol)",
       y = "Molecular content") +
  theme_minimal()



# - - - 
# Write prospect model output to csv
write.csv(final_results, file = "./prospect_rtm_output.csv", row.names = FALSE)
