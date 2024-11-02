##' Objective: Fit prospect model to "low" and "high" nitrate treatments.
##' 
##' @author Nathan Malamud
##' @date 2024.09.27
##' 
##' Reference:
##'   https://jbferet.gitlab.io/prospect/index.html
##'   
##' Chat-GPT4 Log:
##'  https://chatgpt.com/share/66f7258c-28a4-8012-84b4-4f00c74ff8f6

library(tidyverse)
library(prospect)
library(dplyr)

# Assume relationship between leaf reflectance and absorptance
# R/T ratio is equivalent to R / (1 - R) and (1 - T) / T
spec_wide <- read_csv("./data/spec_data_wide.csv")

# Reorder columns
spec_wide <- spec_wide %>% 
  select(barcodeID, sampleID, species, treatment_mmol, everything())

# Filter out wavelengths outside of 500-2400 nm
FULL <- seq(500.0, 2500.0, .01)
NIR <- seq(700.0, 1100.0, .01)
SWIR <- seq(1100.0, 2500.0, .01)

# Remove duplicate items 
spec_wide <- spec_wide %>% distinct(barcodeID, .keep_all = T) %>% distinct(sampleID, .keep_all = T)

# Pivot dataframe to get "lambda", "reflectance", and "transmittance" values
spec_long <- spec_wide |>
  pivot_longer(cols = c(5:1028), names_to = "lambda", values_to = "reflectance") |>
  mutate(
    lambda = as.numeric(lambda),       # Convert 'lambda' to numeric
    reflectance = as.numeric(reflectance) # Convert 'reflectance' to numeric
  ) %>% filter(lambda %in% FULL)

# Aggregate lambda values via rounding -> mean reflectance
spec_long$lambda <- round(spec_long$lambda, 0)
spec_long_agg <- aggregate(reflectance ~ barcodeID + species + treatment_mmol + lambda, data = spec_long, FUN=mean, na.rm = T)
spec_long_agg$transmittance = 1 - spec_long_agg$reflectance

# Estimate all parameters for PROSPECT-D
# Source: https://jbferet.gitlab.io/prospect/articles/prospect5.html
# define set of parameters to be assessed
# Create an empty list to store results for each leaf
results_list <- list()

# Get all unique barcodeIDs
unique_barcodes <- unique(spec_long_agg$barcodeID)

# Define parameters to estimate for PROSPECT-D (or PRO)
Parms2Estimate <- c('N', 'CHL', 'CAR', 'EWT')

# Loop over each unique barcodeID
for (barcode in unique_barcodes) {
  
  # Subset data for the current leaf
  leaf_spec <- spec_long_agg %>% filter(barcodeID == barcode)
  
  # Adjust spectral domain
  SubData <- FitSpectralData(lambda = leaf_spec$lambda,
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
    N = OutPROSPECTPRO$N,
    CHL = OutPROSPECTPRO$CHL,
    CAR = OutPROSPECTPRO$CAR
  )
  
  # Append results to the list
  results_list[[as.character(barcode)]] <- estimated_params
}

# Combine all results into a single dataframe
final_results <- do.call(rbind, results_list)

# Merge with original data
final_results <- merge(final_results, spec_wide %>%
                         select(barcodeID, sampleID, species, treatment_mmol), by = "barcodeID")

# Write prospect model output to csv
# Units: 
write.csv(final_results, file = "./prospect_rtm_output.csv", row.names = FALSE)

# View the final results dataframe
# Box plots across species and reatments
josef_colors <- c("#299680", "#7570b2", "#ca621c")

final_results %>%
  ggplot(aes(y = N, x = treatment_mmol, color = species)) +
  geom_boxplot(aes(alpha = 0.25)) +
  facet_wrap(~species, ncol = 3) +
  scale_color_manual(values = josef_colors) +
  theme_classic() +
  theme(legend.position = "none")
