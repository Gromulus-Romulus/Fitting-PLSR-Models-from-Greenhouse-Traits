# Use random forest for feature selection, QC,
# and validation of measurements
#
# author: Nathan Malamud
# date: 2024.10.29

# Analysis Packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(GGally)
library(dplyr)
library(readxl)
library(reshape2)
library(scales)
library(LeafArea)
library(stringr)
library(RColorBrewer)
library(caret)

# Load data
# REMINDER: Set Working Directory -> Source File Location
data <- read_csv("./data/traits.csv")

# Load prospect measurements from spec curves
prospect <- read_csv("./data/molecular_content.csv") %>%
  select(-c(sampleID, species, treatment_mmol))
names(prospect) <- c("barcodeID", "CHL", "CAR", "ANT", "N_layers", "PROT", "CBC")
data <- merge(data, prospect, by = "barcodeID")

# Define Factor Levels (Treatment and Species)
data$treatment_mmol <- data$treatment_mmol
data$species <- as.factor(data$species)
levels(data$species) <- c("R. sativus", "B. officinalis", "H. vulgare")
data$LDMC <- as.numeric(data$LDMC)

# Create a new column for aggregated treatment levels
data <- data %>%
  mutate(treatment_level = case_when(
    treatment_mmol <= 5 ~ "0 - 5 mmol", 
    treatment_mmol <= 15 ~ "10 - 15 mmol",
    treatment_mmol <= 25 ~ "20 - 25 mmol", 
    treatment_mmol <= 35 ~ "30 - 35 mmol",
    TRUE ~ "Other"  # Optional: catch any values above 35
  ))

# Define colors for each species
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c")

# TODO clean up plots
# Verify that ETR measurements aren't due to variation in vpdl and Qamb
# Make boxplot of Qamb across treatment levels and species, color by species
p1 <- ggplot(data, aes(x = treatment_level, y = Qamb, fill = species)) +
  geom_boxplot() +
  scale_fill_manual(values = josef_colors) + theme_minimal() + facet_wrap(~species) +
  labs(x = "Treatment Level (mmol)", y = "Qamb (μmol photons m-2 s-1)")

# Regression plot showing Qamb confounds ETR measurements
p2 <- ggplot(data, aes(x = Qamb, y = ETR)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color="red") +
  labs(x = "Qamb (μmol photons m-2 s-1)", y = "ETR (μmol electrons m-2 s-1)") +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.25, color = "grey80"),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5, color = "black"),
    strip.text = element_text(hjust = 0, size = 10, face = "italic"),
    legend.position = "none"
  )

# Verify that gsw measurements aren't due to variation in vpdl
# Make boxplot of gsw across treatment levels and species, color by species
p3 <- ggplot(data, aes(x = treatment_level, y = vpdl, fill = species)) +
  geom_boxplot() +
  scale_fill_manual(values = josef_colors) + theme_minimal() + facet_wrap(~species) +
  labs(x = "Treatment Level (mmol)", y = "vpdl")

# Regression plot showing vpdl confounds gsw measurements
p4 <- ggplot(data, aes(x = vpdl, y = gsw)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color="red") +
  labs(x = "vpdl", y = "gsw") +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.25, color = "grey80"),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5, color = "black"),
    strip.text = element_text(hjust = 0, size = 10, face = "italic"),
    legend.position = "none"
  )

# We have shown that vpdl and qamb are confounding factors for
# physiology measurements, therefore they can't be used.

# - - - - - - - - -
# GGgally plot showing feature correlations
# Select relevant columns and include species
correlation_data <- data %>%
  select(treatment_level, N_layers, CHL, LMA, LDMC, dry_whole_g)

# Plot pairwise correlations with species color
# R Color brewer themes: display.brewer.all()
library(GGally)
library(RColorBrewer)

# Define the color palette manually
# Define a custom color palette
#   See palettes: display.brewer.all()
# Define a custom theme for all plots
custom_theme <- theme_minimal(base_family = "sans") + 
  theme(
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.25, color = "grey80"),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5, color = "black"),
    strip.text = element_text(hjust = 0, size = 10, face = "italic"),
    legend.position = "none"
  )

correlation_data <- data %>%
  select(treatment_level, N_layers, CHL, LMA, LDMC, dry_whole_g)
  

custom_colors <- brewer.pal(6, "YlGn")[c(3:6)]
names(custom_colors) <- levels(correlation_data$treatment_level)

# Remove the treatment_level column from the data
correlation_data_subset <- correlation_data[, !colnames(correlation_data) %in% "treatment_level"]

# Generate the plot without the treatment_level column
ggpairs(correlation_data_subset,
        mapping = aes(color = correlation_data$treatment_level, fill = correlation_data$treatment_level),
        upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal()

# Now, physiology and molecular content
correlation_data <- data %>%
  select(treatment_level, CHL, N_layers, CAR, ANT, Phi_PS2)

# Remove Phi_PS2 outliers (< 0.5)
correlation_data <- correlation_data %>%
  filter(Phi_PS2 > 0.5)

custom_colors <- brewer.pal(6, "YlGn")[c(3:6)]
names(custom_colors) <- levels(correlation_data$treatment_level)

# Remove the treatment_level column from the data
correlation_data_subset <- correlation_data[, !colnames(correlation_data) %in% "treatment_level"]

# Generate the plot without the treatment_level column
ggpairs(correlation_data_subset,
        mapping = aes(color = correlation_data$treatment_level, fill = correlation_data$treatment_level),
        upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal()

# Look at changes in mesophyll structure and nitrate absorption
# Across species and treatment levels
correlation_data <- data %>%
  select(species, treatment_level, N_layers, CHL, LMA, LDMC, dry_whole_g)

p1 <- ggplot(correlation_data, aes(x = treatment_level, y = N_layers, color = species)) +
  geom_boxplot() +
  scale_color_manual(values = josef_colors) + theme_minimal() + facet_wrap(~species) +
  labs(x = "Treatment Level (mmol)", y = "N_layers")

p2 <- ggplot(correlation_data, aes(x = treatment_level, y = CHL, color = species)) +
  geom_boxplot() +
  scale_color_manual(values = josef_colors) + theme_minimal() + facet_wrap(~species) +
  labs(x = "Treatment Level (mmol)", y = "CHL")

p3 <- ggplot(correlation_data, aes(x = treatment_level, y = LMA, color = species)) +
  geom_boxplot() +
  scale_color_manual(values = josef_colors) + theme_minimal() + facet_wrap(~species) +
  labs(x = "Treatment Level (mmol)", y = "LMA")

ggarrange(p1, p2, p3, ncol = 1)

# - - - - - - - - 
# Ok, let's build a random forest model with candidate features
forest_data <- data %>%
  select(N_layers, CHL, LMA, LDMC, dry_whole_g)

# Treat species as numeric variable
#forest_data$species <- as.numeric(forest_data$species)

# Set up training control without ROC and for regression
control <- trainControl(method = "repeatedcv", number = 10, repeats = 3)

# Fit the random forest model for regression
set.seed(1234)
rf <- train(
  dry_whole_g ~ ., 
  data = forest_data,
  method = "rf",
  metric = "RMSE",  # Use RMSE for regression
  tuneGrid = expand.grid(mtry = 1:10),
  trControl = control
)

# Plot variable importance
plot(varImp(rf), main = "Variable Importance with Random Forest")

