##' Build functional group LDA classifier
##' using thesis greenhouse data.
##' 
##' Based on Anna Scheiger's "speciesID" code:
##'   https://github.com/annakat/speciesID/
##'   
##' Associated Paper: https://doi.org/10.1111/1365-2745.13972
##' 
##' @author [Nathan Malamud]
##' @date 2024.10.17
##' 

# Load libraries
library(readr)
library(pls)
library(MASS)
library(caret) # For confusion matrix and evaluation metrics
library(readxl)
library(reshape)
library(tidyverse) 
library(dplyr)
library(ggplot2)
library(gridExtra)

# Load spectral matrices for PLSR model (made by signal_cleaning.R)
signal.matrices <- readRDS("./data/signal.matrices.rds")

# - - - - -
# 1) TODO: See which matrices (X_raw, X_avg, X_d1, X_d2)
# Reduce RMSEP error multiple traits (e.g. LMA, LDMC, EWT)
# ... by using PLSR model with different signal matrices
X_raw <- signal.matrices$X_raw
X_avg <- signal.matrices$X_avg
X_d1 <- signal.matrices$X_d1
X_d2 <- signal.matrices$X_d2

# Important - the uniq_ids array
# helps us translate between matrix row numbers
# and barcodeIDs for retrieving sample metadata
barcodeID <- signal.matrices$uniq_ids
spectral_columns <- colnames(X_raw) %>% as.character()

# Load traits data, only include barcodes that we
# have spectral curves for
traits <- read.csv("./data/traits.csv") %>%
  filter(barcodeID %in% signal.matrices$uniq_ids)

# Save metadata in separate dataframe
# Note: naming conflict with MASS library - "select"
metadata <- traits %>%
  dplyr::select(barcodeID, sampleID, treatment_mmol, species)

# Save barcode metadata in separate dataframe
# Extract prediction matrix Y
Y <- traits %>%
  dplyr::select("LMA", "LDMC", "EWT") %>% as.matrix()

#species = as.factor(metadata$species)
plsr_data <- cbind(species, X_d2) %>%
  as.data.frame()

# ensure this is a factor
plsr_data$species <- as.factor(plsr_data$species)

# Partition data into LO and HI categories
# Comparative statistics using cross-validation
# Calculate RMSEP, R2, PRESS for each model and cross-compare
THRESH = 15 # mmol
metadata$treatment_category <- ifelse(metadata$treatment_mmol > THRESH, "high", "low")

plsr_data_hi <- plsr_data %>%
  filter(metadata$treatment_category == "high") # 52 obs > 15 mmol

plsr_data_lo <- plsr_data %>%
  filter(metadata$treatment_category == "low") # 45 obs <= 15 mmol

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 2) Build LDA classifier
# remove near zero variance predictors

# TODO: FIX SPECIES CLASS LABELS
# Step 1: Ensure the species column in both datasets have consistent factor levels
species_levels <- levels(plsr_data_lo$species)  # Assuming plsr_data_lo has correct species levels
plsr_data_hi$species <- factor(plsr_data_hi$species, levels = species_levels)

# Fit LDA models for both low and high datasets
lda_lo <- MASS::lda(species ~ ., data=plsr_data_lo)
lda_hi <- MASS::lda(species ~ ., data=plsr_data_hi)

# Classify the data using the LDA models
prediction_lohi <- predict(lda_lo, newdata=plsr_data_hi)
prediction_lolo <- predict(lda_lo, newdata=plsr_data_lo)
prediction_hihi <- predict(lda_hi, newdata=plsr_data_hi)
prediction_hilo <- predict(lda_hi, newdata=plsr_data_lo)

# Rename the predicted classes based on actual species names for all predictions
levels(prediction_lohi$class) <- species_levels
levels(prediction_lolo$class) <- species_levels
levels(prediction_hihi$class) <- species_levels
levels(prediction_hilo$class) <- species_levels

# Create data frames for all predictions with LD1, LD2, and species names
df_lohi <- data.frame(LD1 = prediction_lohi$x[,1], LD2 = prediction_lohi$x[,2], species = prediction_lohi$class)
df_lolo <- data.frame(LD1 = prediction_lolo$x[,1], LD2 = prediction_lolo$x[,2], species = prediction_lolo$class)
df_hihi <- data.frame(LD1 = prediction_hihi$x[,1], LD2 = prediction_hihi$x[,2], species = prediction_hihi$class)
df_hilo <- data.frame(LD1 = prediction_hilo$x[,1], LD2 = prediction_hilo$x[,2], species = prediction_hilo$class)

# Plot function with ellipsoids and fixed aspect ratio
plot_lda <- function(data, title) {
  ggplot(data, aes(x = LD1, y = LD2, color = species)) +
    geom_point(size = 2) +
    stat_ellipse(level = 0.95, aes(fill = species), alpha = 0.1, geom = "polygon") +
    scale_color_manual(values = josef_colors) +
    scale_fill_manual(values = josef_colors) +
    labs(title = title) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    coord_fixed(ratio=2, xlim=c(-15, 30), ylim=c(-6.5, 6.5))  # Ensure fixed aspect ratio
}

# Create individual plots for each LDA prediction
p1 <- plot_lda(df_lohi, "LDA: lo (trained) -> hi (predicted)")
p2 <- plot_lda(df_lolo, "LDA: lo (trained) -> lo (predicted)")
p3 <- plot_lda(df_hihi, "LDA: hi (trained) -> hi (predicted)")
p4 <- plot_lda(df_hilo, "LDA: hi (trained) -> lo (predicted)")

# Arrange plots in a 2x2 grid with fixed aspect ratio
grid.arrange(p1, p2, p3, p4, ncol = 2)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 1: Combine low and high datasets
plsr_data_combined <- rbind(plsr_data_lo, plsr_data_hi)

# Step 2: Fit LDA model to the full PLSR data
lda_combined <- MASS::lda(species ~ ., data = plsr_data_combined)

# Step 3: Classify the data using the LDA model fit on combined data
prediction_combined_lohi <- predict(lda_combined, newdata = plsr_data_hi)
prediction_combined_lolo <- predict(lda_combined, newdata = plsr_data_lo)

# Step 4: Rename the predicted classes based on actual species names
levels(prediction_combined_lohi$class) <- species_levels
levels(prediction_combined_lolo$class) <- species_levels

# Create data frames for predictions with LD1, LD2, and species names
df_combined_lohi <- data.frame(LD1 = prediction_combined_lohi$x[, 1], LD2 = prediction_combined_lohi$x[, 2], species = prediction_combined_lohi$class)
df_combined_lolo <- data.frame(LD1 = prediction_combined_lolo$x[, 1], LD2 = prediction_combined_lolo$x[, 2], species = prediction_combined_lolo$class)

# Create individual plots for each LDA prediction using combined model
p5 <- plot_lda(df_combined_lohi, "LDA : combined -> hi (predicted)")
p6 <- plot_lda(df_combined_lolo, "LDA : combined -> lo (predicted)")

# Arrange the new plots with the existing ones in a 2x3 grid
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TODO: Confusion matrices
# Generate the confusion matrix
conf_matrix <- confusionMatrix(plsr_data_lo$species, prediction_hilo$class)

# Extract the table from the confusion matrix
conf_table <- as.data.frame(conf_matrix$table)

# Rename columns for clarity
colnames(conf_table) <- c("Reference", "Prediction", "Count")

# Create a ggplot of the confusion matrix
ggplot(data = conf_table, aes(x = Reference, y = Prediction, fill = Count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#F2F2F2", high = "#4D4D4D") +  # Use Josef colors for the fill
  geom_text(aes(label = Count), color = "black") +
  theme_minimal() +
  labs(title = "Confusion Matrix", x = "True Class", y = "Predicted Class") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# - - - -

# Define colors for each species and grey for non-significant points
josef_colors <- c("R. sativus" = "#299680", "B. officinalis" = "#7570b2", "H. vulgare" = "#ca621c", "grey" = "#d3d3d3")

# Function to calculate significant regressions
calculate_significant_regressions <- function(data, x, y) {
  data %>%
    group_by(species, treatment_level) %>%
    summarize(
      fit = list(lm(reformulate(x, y), data = cur_data())),
      r_squared = summary(fit[[1]])$r.squared,
      p_value = summary(fit[[1]])$coefficients[2, 4]
    ) %>%
    filter(p_value < 0.05)  # Keep only significant results
}

# Custom theme with italicized treatment labels for consistency
custom_theme <- theme_minimal(base_family = "Arial") +
  theme(
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.25, color = "grey80"),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5, color = "black"),
    strip.text = element_text(hjust = 0, size = 10, face = "italic"),  # Italicize treatment levels
    legend.position = "none"
  )

# Standardized plot function
plot_correlation <- function(data, x, y, x_label, y_label) {
  ggplot(data, aes_string(x = x, y = y)) +
    geom_point(aes(color = color), size = 1.5) +
    scale_color_identity() +
    geom_smooth(data = subset(data, !is.na(p_value) & species == "R. sativus"),
                aes(color = "#299680"), method = "lm", se = TRUE, size = 0.8) +
    geom_smooth(data = subset(data, !is.na(p_value) & species == "B. officinalis"),
                aes(color = "#7570b2"), method = "lm", se = TRUE, size = 0.8) +
    geom_smooth(data = subset(data, !is.na(p_value) & species == "H. vulgare"),
                aes(color = "#ca621c"), method = "lm", se = TRUE, size = 0.8) +
    labs(x = x_label, y = y_label) +
    # Color-coded R^2 and p-value annotations for each species, consistent annotation position
    stat_cor(
      data = subset(data, !is.na(p_value) & species == "R. sativus"),
      aes(label = paste(..rr.label.., ..p.label.., sep = "~~~"), color = "#299680"),
      label.x.npc = 0.65, label.y.npc = 0.75,
      size = 2.5, method = "pearson"
    ) +
    stat_cor(
      data = subset(data, !is.na(p_value) & species == "B. officinalis"),
      aes(label = paste(..rr.label.., ..p.label.., sep = "~~~"), color = "#7570b2"),
      label.x.npc = 0.65, label.y.npc = 0.65,
      size = 2.5, method = "pearson"
    ) +
    stat_cor(
      data = subset(data, !is.na(p_value) & species == "H. vulgare"),
      aes(label = paste(..rr.label.., ..p.label.., sep = "~~~"), color = "#ca621c"),
      label.x.npc = 0.65, label.y.npc = 0.55,
      size = 2.5, method = "pearson"
    ) +
    custom_theme
}

# Prepare data for each new relationship in the specified order
new_relationships <- list(
  list(x = "CHL", y = "LMA", x_label = "CHL", y_label = "LMA"),
  list(x = "CAR", y = "LMA", x_label = "CAR", y_label = "LMA"),
  list(x = "ANT", y = "LMA", x_label = "ANT", y_label = "LMA")
)

new_plots <- list()

# Generate plots for each new relationship
for (rel in new_relationships) {
  p_values <- calculate_significant_regressions(data, rel$x, rel$y)
  data_rel <- data %>%
    left_join(p_values, by = c("species", "treatment_level")) %>%
    mutate(color = case_when(
      !is.na(p_value) & species == "R. sativus" ~ "#299680",
      !is.na(p_value) & species == "B. officinalis" ~ "#7570b2",
      !is.na(p_value) & species == "H. vulgare" ~ "#ca621c",
      TRUE ~ "#d3d3d3"
    )) %>% mutate(p_value = round(p_value, 2))
  
  plot <- plot_correlation(data_rel, rel$x, rel$y, rel$x_label, rel$y_label) +
    facet_grid(~ treatment_level)
  
  new_plots[[length(new_plots) + 1]] <- plot
}

# Arrange the new set of plots in a 3-row layout with a common legend at the bottom
ggarrange(plotlist = new_plots, ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

