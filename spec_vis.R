##' Visualize data from spectral curves with multiple samples
##' and add averaged reflectance with confidence intervals
##' @author Nathan Malamud
##' @date 2024.09.20

library(tidyverse)

# Load spectral data
spec <- readRDS("./proc/spec.RDS")

# Aggregate data according to species, treatment, and individual samples
spec <- spec %>%
  select(-c(barcodeID, sampleID)) %>%
  melt(id = c("species", "treatment_mmol"))
names(spec) <- c("species", "treatment_mmol", "wavelength", "reflectance")
spec$wavelength <- as.numeric(spec$wavelength)
spec$reflectance <- as.numeric(spec$reflectance)

# Calculate mean and confidence intervals for each species, treatment, and wavelength
spec_summary <- spec %>%
  group_by(species, treatment_mmol, wavelength) %>%
  summarise(mean_reflectance = mean(reflectance),
            sd_reflectance = sd(reflectance),
            n = n()) %>%
  mutate(se_reflectance = sd_reflectance / sqrt(n),
         lower = mean_reflectance - 1.96 * se_reflectance,
         upper = mean_reflectance + 1.96 * se_reflectance)

# Create plot
spec %>%  subset(treatment_mmol != 0) |>  # Exclude control treatment
  ggplot(aes(x = wavelength, y = reflectance)) +
  geom_point(alpha=0.5) +
  facet_wrap(species ~ treatment_mmol, nrow = 3) +
  theme_classic() +
  labs(title = "Plot of spectral measurements for all species and treatments",
       x = "Wavelength (nm)",
       y = "Reflectance") +
  xlim(330, 1000)

# Export to PDF
ggsave("spectral_curves_overlay.pdf", plot = p, device = "pdf", width = 10, height = 8)
