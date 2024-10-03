#------ ------ Plotting average spectra by species ------------------ #
library(rcartocolor)
library(MetBrewer)

Start.wave <- 500
End.wave <- 2400
wv <- seq(Start.wave, End.wave, .1)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Nate's Code (spectra_visualization.R)
# Load traits and spectroscopy data
traits <- readRDS("../proc/traits.RDS")
spec <- readRDS("../proc/spec.RDS")

# Filter only samples with both spec and trait data
barIDs <- intersect(spec$barcodeID, traits$barcodeID)
traits <- traits %>% filter(barcodeID %in% barIDs)
spec <- spec %>% filter(barcodeID %in% barIDs)

# Remove NA rows from spec dataframe
spec <- spec %>% filter(complete.cases(.))

# Filter spectra data frame to get wavelengths in specified range
spec <- spec |>
  pivot_longer(cols = c(5:1028), names_to="wavelength", values_to="reflectance") %>%
  filter(as.numeric(wavelength) >= Start.wave & as.numeric(wavelength) <= End.wave) %>%
  pivot_wider(names_from="wavelength", values_from="reflectance")

# Remove duplicate items
spec <- spec %>% distinct(barcodeID, .keep_all = TRUE) %>% distinct(sampleID, .keep_all = TRUE)

# Update wavelength vector
wv <- as.numeric(colnames(spec[5:length(colnames(spec))]))

# Rename columns of spectra dataframe
names(spec)[5:length(names(spec))] <- paste0("Wave_", names(spec)[5:length(names(spec))])

# Focus on relevant traits
traits <- traits %>% select(barcodeID, sampleID, species, treatment_mmol, LDMC, LMA, Phi_PS2, Fm_prime, dry_whole_g)

# Merge traits and spectral data
plsr_data <- merge(traits, spec, by=c("barcodeID", "sampleID", "species", "treatment_mmol"))

# Create a lookup table for species to Latin names
latin_names <- c("borage" = "B. officinalis",
                 "barley" = "H. vulgare",
                 "radish" = "R. sativus")

# Rename species column using the lookup table
plsr_data <- plsr_data %>%
  mutate(species = recode(species,
                          "borage" = latin_names["borage"],
                          "barley" = latin_names["barley"],
                          "radish" = latin_names["radish"]))

# Filter dataset for LO and HI treatments
# Categorize treatments based on numeric values
plsr_data <- plsr_data %>%
  mutate(treatment_category = ifelse(treatment_mmol %in% c(0, 5, 10),
                                     "0-10 mmol (low)", ">15 mmol (high)"))

lo_hi_data <- plsr_data %>%
  filter(treatment_category %in% c("0-10 mmol (low)", ">15 mmol (high)"))

# Pivot data for plotting
long_spectra <- lo_hi_data %>%
  pivot_longer(cols = starts_with("Wave_"), names_to = "wavelength", values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(gsub("Wave_", "", wavelength)))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# # ---- per-site ---- #
# # rename site 6
# sspectra <- spectra %>% mutate(siteID = replace(siteID, siteID == 'HS', '6'))
# # remove missing sites
# sspectra <- sspectra %>% filter(!is.na(siteID))
# per_site <- sspectra %>% group_by(siteID,  lambda) %>% 
#   summarize(mean_r = mean(reflectance),
#             sd_r = sd(reflectance),
#             n = n(),
#             se_r = sd_r / sqrt(n)) %>% ungroup()

# ---- per-aspect ---- #
# rename site 6
# sspectra <- spectra %>% mutate(siteID = replace(siteID, siteID == 'HS', '6'))
# # remove missing sites
# sspectra <- sspectra %>% filter(!is.na(siteID))
# per_aspect <- sspectra %>% group_by(siteID, aspect, lambda) %>% 
#   summarize(mean_r = mean(reflectance),
#             sd_r = sd(reflectance),
#             n = n(),
#             se_r = sd_r / sqrt(n)) %>% ungroup()

# ------ all per-species ------ # 
# calculate mean, sd per-species, per-wavelength
# filter to focal species

# all_per_species <- long_spectra %>% group_by(species, wavelength) %>% 
#   summarize(mean_r = mean(reflectance),
#             sd_r = sd(reflectance),
#             n = n(),
#             se_r = sd_r / sqrt(n)) %>% ungroup()


# ------ per-species ------ # 

# calculate mean, sd per-species, per-wavelength
# # filter to focal species
# sspectra <- long_spectra %>% mutate(species = replace(species, species == 'Helichrysum nudifolium', 'Helichrysum nudifolium/pallidum')) %>% 
#   mutate(species = replace(species, species == 'Helichrysum pallidum', 'Helichrysum nudifolium/pallidum')) %>%
#   mutate(species = replace(species, species == 'Senecio glaberimous', 'Senecio glaberrimus'))
# 
# focals <- c( 'Helichrysum nudifolium/pallidum', "Senecio glaberrimus", "Helichrysum ecklonis")
# sspectra <- sspectra %>% filter(species %in% focals)
# 
# per_species <- sspectra %>% group_by(species, lambda) %>% 
#   summarize(mean_r = mean(reflectance),
#             sd_r = sd(reflectance),
#             n = n(),
#             se_r = sd_r / sqrt(n)) %>% ungroup()
  
# ------ per-species and treatment ------ # 

per_species_treatment <- long_spectra %>% group_by(species, treatment_mmol, wavelength) %>% 
  summarize(mean_r = mean(reflectance),
            sd_r = sd(reflectance),
            n = n(),
            se_r = sd_r / sqrt(n)) %>% ungroup()

per_species_treatment_category <- long_spectra %>% group_by(species, treatment_category, wavelength) %>% 
  summarize(mean_r = mean(reflectance),
            sd_r = sd(reflectance),
            n = n(),
            se_r = sd_r / sqrt(n)) %>% ungroup()

# ------ per-barcode ------ # 

# # include response variable columns
# #  merge in dry weight column, fresh weight column
# per_barcode <- long_spectra %>% group_by(species, treatment_mmol, barcode, lambda) %>% 
#   summarize(mean_l = mean(reflectance), n = n()) %>% ungroup()

# filter out scans with no dry weight
# per_barcode <- per_barcode %>% filter(!is.na(dry_mass))

# save out to give to nicole
# write.csv(per_barcode, "./plsr_per_leaf.csv")

# # keep only not messy wavelengths
# per_site <- per_site %>% filter(lambda > 400)
# per_species <- per_species %>% filter(lambda > 400)
# all_per_species <- all_per_species %>% filter(lambda > 400)
# per_species_site <- per_species_site %>% filter(lambda > 400)


#plot focal species
speciesplot <- ggplot(data = per_species_treatment_category,
                      aes(x = wavelength, color = treatment_category, fill = treatment_category)) +
  geom_ribbon(aes(ymin = 100*(mean_r-2*se_r), ymax = 100*(mean_r+2*se_r), fill = treatment_category), alpha = 0.3, color = NA) +
  ylab("% Reflectance") + 
  xlab("Wavelength (nm)") +
  geom_line(aes(y = 100*mean_r)) + 
  scale_fill_manual(name = "Treatment Category", values = c("0-10 mmol (low)" = "darkgrey", ">15 mmol (high)" = "red")) +
  scale_color_manual(name = "Treatment Category", values = c("0-10 mmol (low)" = "darkgrey", ">15 mmol (high)" = "red")) +
  theme_minimal() + 
  facet_wrap(~species) +
  theme(
    text = element_text(family = "Helvetica", face = "plain"),  # Regular font for all text
    axis.title = element_text(size = 12, face = "bold"),  # Axis titles in regular
    axis.text = element_text(size = 10, face = "plain"),   # Axis text in regular
    legend.title = element_text(size = 12, face = "bold"),# Legend title in regular
    legend.text = element_text(size = 10, face = "bold"), # Legend text in regular
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold.italic")  # Species names in italic in facet labels
  )
speciesplot
# ggsave("per_focal_species.pdf", speciesplot, units = "mm", width = 200, height = 100)

library(dplyr)

# Create a hash table as a named numeric vector (0-35 in intervals of 5)
treatment_map <- c("0" = 0, "5" = 5, "10" = 10, "15" = 15, "20" = 20, 
                   "25" = 25, "30" = 30, "35" = 35)

# Convert the factor to numeric using the hash table
per_species_treatment <- per_species_treatment %>%
  mutate(treatment_mmol_numeric = as.numeric(treatment_map[as.character(treatment_mmol)]))

# Now 'treatment_mmol_numeric' is a numeric variable with the correct numeric values


# Plot treatments across species
treatmentplot <- ggplot(data = per_species_treatment,
                        aes(x = wavelength)) +
  geom_ribbon(aes(ymin = 100*(mean_r - 2*se_r), ymax = 100*(mean_r + 2*se_r)), fill = "grey80", alpha = 0.3) +
  ylab("% Reflectance") + 
  xlab("Wavelength (nm)") +
  geom_line(aes(y = 100*mean_r, color = treatment_mmol_numeric)) + 
  scale_color_viridis_c(name = "Nitrate Treatment (mmol)", option = "D") +
  theme_minimal() + 
  facet_wrap(~species)
treatmentplot

# 
# all_per_species <- all_per_species %>% filter(n > 5)
# all_per_species <- data.table(all_per_species)
# all_per_species[, c("genus") := tstrsplit(species, " ", keep = c(1))]
# #plot all species
# allspeciesplot <- ggplot(data = all_per_species,
#                       aes(x = lambda, color = genus, fill=species)) +
# #  geom_ribbon(aes(ymin = 100*(mean_r-2*se_r), ymax = 100*(mean_r+2*se_r), fill = genus), alpha = 0.3, color = NA) +
#   ylab("% Reflectance") + 
#   xlab("Wavelength (nm)") +
#   geom_line(aes(y = 100*mean_r), linewidth=1.5, alpha=.8) + 
#   # scale_fill_carto_d(name = "genus", palette='Vivid') +
#   scale_color_carto_d(name = "genus", palette='Vivid') +
#   # scale_fill_met_d("Signac") +
#   scale_color_met_d("Signac") +
#   theme_minimal()
# allspeciesplot
# ggsave("per_all_species.pdf", allspeciesplot, units = "mm", width = 200, height = 100)
# 
# #plot per-site species
# siteplot <- ggplot(data = per_site,
#                       aes(x = lambda, color = siteID, fill = siteID)) +
#   geom_ribbon(aes(ymin = 100*(mean_r-2*se_r), ymax = 100*(mean_r+2*se_r), fill = siteID), alpha = 0.3, color = NA) +
#   ylab("% Reflectance") + 
#   xlab("Wavelength (nm)") +
#   geom_line(aes(y = 100*mean_r)) + 
#   scale_fill_carto_d(name = "Site") +
#   scale_color_carto_d(name = "Site") +
#   scale_fill_met_d("Lakota") +
#   scale_color_met_d("Lakota") +
#   theme_minimal()
# siteplot
# ggsave("per_site.pdf", siteplot, units = "mm", width = 200, height = 100)
# 
# #plot per-site and per-species
# sitespecplot$siteID <- as.character(sitespecplot$siteID)
# per_species_site <- per_species_site %>% filter(!is.na(siteID))
# sitespecplot <- ggplot(data = per_species_site,
#                    aes(x = lambda, color = factor(siteID), fill = factor(siteID))) +
#   facet_grid(.~species) +
#   geom_ribbon(aes(ymin = 100*(mean_r-2*se_r), ymax = 100*(mean_r+2*se_r), fill = factor(siteID)), alpha = 0.3, color = NA) +
#   ylab("% Reflectance") + 
#   xlab("Wavelength (nm)") +
#   geom_line(aes(y = 100*mean_r)) +
#   scale_fill_carto_d(name = "Site") +
#   scale_color_carto_d(name = "Site") +
#   scale_fill_met_d("Lakota") +
#   scale_color_met_d("Lakota") +
#   theme_minimal()
# 
# sitespecplot
# ggsave("per_site_per_spec.pdf", sitespecplot, units = "mm", width = 200, height = 100)
# 
# 
# 
# # plot per-aspect??
# 
# per_aspect <- per_aspect %>% filter(!is.na(aspect))
# aspectplot <- ggplot(data = per_aspect,
#                    aes(x = lambda, color = factor(siteID), fill = factor(siteID))) +
#   geom_ribbon(aes(ymin = 100*(mean_r-2*se_r), ymax = 100*(mean_r+2*se_r), fill = factor(siteID)), alpha = 0.3, color = NA) +
#   facet_grid(.~aspect) +
#   ylab("% Reflectance") + 
#   xlab("Wavelength (nm)") +
#   geom_line(aes(y = 100*mean_r)) + 
#   scale_fill_carto_d(name = "Site") +
#   scale_color_carto_d(name = "Site") +
#   scale_fill_met_d("Lakota") +
#   scale_color_met_d("Lakota") +
#   theme_minimal()
# aspectplot
# ggsave("per_site_aspect.pdf", aspectplot, units = "mm", width = 200, height = 100)
