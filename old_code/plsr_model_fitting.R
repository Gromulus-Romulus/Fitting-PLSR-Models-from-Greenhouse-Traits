# Playing around with Shawn Serbin's github
# repo for learning PLSR pipelines
#
# Documentation guide (Serbin et al. 2022): 
#   https://github.com/plantphys/spectratrait/blob/main/spectratrait_1.2.5.pdf
#
##' @author Nathan malamud
##' @date 2024.09.23

# clear namespace
rm(list=ls())

# load dependencies
library(pls)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(spectratrait)
library(tidyverse)

# - - - - - - - - - - - - - - - - - - - -
# Installation Guide
# devtools::install_github(
# repo = "plantphys/spectratrait",
# dependencies=TRUE)
# - - - - - - - - - - - - - - - - - - - -

# Define useful dplyr "macros"
`%notin%` <- Negate(`%in%`)

# Script options
pls::pls.options(plsralg = "oscorespls")
pls::pls.options("plsralg")
output_dir <- "tempdir"

Start.wave <- 500
End.wave <- 2400
wv <- seq(Start.wave, End.wave, .1)

# Load traits and spectroscopy data
# From ./proc (processed raw data) directory
traits <- readRDS("./proc/traits.RDS")
spec <- readRDS("./proc/spec.RDS")

# Important: Not all leaves measured have spec readings, and vice versa.
# We filter only the barcode IDs that exist in both the X and Y matrices.
# Let's do that quality check here.
barIDs <- intersect(spec$barcodeID, traits$barcodeID)

traits <- traits %>% filter(barcodeID %in% barIDs)
spec <- spec %>% filter(barcodeID %in% barIDs)

# Remove NA rows from spec dataframe
spec <- spec %>% filter(complete.cases(.))

# Filter spectra data frame so we only get wavelengths in specified range
spec <- spec |>
  pivot_longer(cols = c(5:1028), names_to="wavelength", values_to="reflectance") %>%
  filter(wavelength %in% wv) %>% pivot_wider(names_from="wavelength", values_from="reflectance")

# Remove duplicate items
spec <- spec %>% distinct(barcodeID, .keep_all = T) %>% distinct(sampleID, .keep_all = T)

# Update wv to only include intervals from svc data
wv <- as.numeric(colnames(spec[5:length(colnames(spec))]))

# Rename columns of spectra dataframe to Wave_<wavelength_value>
names(spec)[5:861] <- paste0("Wave_", names(spec[5:861]))

# For now... only focus on traits of interest
of_interest <- c("LDMC", "LMA", "Phi_PS2", "Fm_prime", "dry_whole_g")
traits <- traits %>% select(barcodeID, sampleID, species, treatment_mmol, of_interest)

# Remove duplicate items
traits <- traits %>% distinct(barcodeID, .keep_all = T) %>% distinct(sampleID, .keep_all = T)

# Merge into single dataframe
# TODO: deprecated merging syntax
plsr_data <- merge(traits, spec, by=c("barcodeID", "sampleID", "species", "treatment_mmol"))

# Filter out duplicate barcodeIDs and sampleIDs (just remove for now)
plsr_data <- plsr_data %>% distinct(barcodeID, .keep_all = T) %>% distinct(sampleID, .keep_all = T)

# View dataframe traits of Full PLSR dataset
head(plsr_data)[, 1:6]

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Create cal/val datasets
## Make a stratified random sampling in the strata species and treatment_mmol

# What is the target variable?
inVar <- "Phi_PS2"

method <- "base" #base/dplyr
# base R - a bit slow
# dplyr - much faster
split_data <- spectratrait::create_data_split(dataset=plsr_data, approach=method, 
                                              split_seed=23452135, prop=0.7, 
                                              group_variables="species")

cal.plsr.data <- split_data$cal_data
val.plsr.data <- split_data$val_data

rm(split_data)

##' @graph Visualize cal / val plots for traits
##' TODO: redo in GGplot and output to pdf file
cal_hist_plot <- qplot(cal.plsr.data[,paste0(inVar)],geom="histogram",
                       main = paste0("Cal. Histogram for ",inVar),
                       xlab = paste0(inVar),ylab = "Count",fill=I("grey50"),col=I("black"),
                       alpha=I(.7))
val_hist_plot <- qplot(val.plsr.data[,paste0(inVar)],geom="histogram",
                       main = paste0("Val. Histogram for ",inVar),
                       xlab = paste0(inVar),ylab = "Count",fill=I("grey50"),col=I("black"),
                       alpha=I(.7))

histograms <- grid.arrange(cal_hist_plot, val_hist_plot, ncol=2)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Format PLSR data for model fitting 
cal_spec <- as.matrix(cal.plsr.data[, which(names(cal.plsr.data) %in% paste0("Wave_",wv))])
cal.plsr.data <- data.frame(cal.plsr.data[, which(names(cal.plsr.data) %notin% paste0("Wave_",wv))],
                            Spectra=I(cal_spec))
head(cal.plsr.data)[1:5]

val_spec <- as.matrix(val.plsr.data[, which(names(val.plsr.data) %in% paste0("Wave_",wv))])
val.plsr.data <- data.frame(val.plsr.data[, which(names(val.plsr.data) %notin% paste0("Wave_",wv))],
                            Spectra=I(val_spec))
head(val.plsr.data)[1:5]

##' @graph Visualize cal / val spectral signatures
## TODO: look into `value` argument deprecation warning
par(mfrow=c(1,2)) # B, L, T, R
spectratrait::f.plot.spec(Z=cal.plsr.data$Spectra,wv=wv,plot_label="Calibration")
spectratrait::f.plot.spec(Z=val.plsr.data$Spectra,wv=wv,plot_label="Validation")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Use permutation to determine optimal number of components
if(grepl("Windows", sessionInfo()$running)){
  pls.options(parallel = NULL)
} else {
  pls.options(parallel = parallel::detectCores()-1)
}

method <- "pls" #pls, firstPlateau, firstMin
random_seed <- 1245565
seg <- 50
maxComps <- 16
iterations <- 80
prop <- 0.70

# TODO: look into number of items != multiple of replacement length
# Running Serbin's demo code on my greenhouse data shows optimal number of components is 8
if (method=="pls") {
  nComps <- spectratrait::find_optimal_components(dataset=cal.plsr.data, targetVariable=inVar,
                                                  method=method, 
                                                  maxComps=maxComps, seg=seg, 
                                                  random_seed=random_seed)
  print(paste0("*** Optimal number of components: ", nComps))
} else {
  nComps <- spectratrait::find_optimal_components(dataset=cal.plsr.data, targetVariable=inVar,
                                                  method=method,
                                                  maxComps=maxComps, iterations=iterations, 
                                                  seg=seg, prop=prop,
                                                  random_seed=random_seed)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Fit final model

## TODO: Two types of validation: leave one out (LOO), CV (cross-validation)
plsr.out <- plsr(as.formula(paste(inVar,"~","Spectra")), scale=FALSE, ncomp=nComps,
                 validation="CV", trace=FALSE, data=cal.plsr.data)
fit <- plsr.out$fitted.values[,1,nComps]
pls.options(parallel = NULL)

# External validation fit stats
par(mfrow=c(1,2)) # B, L, T, R
pls::RMSEP(plsr.out, newdata = val.plsr.data)

##' @graph RMSEP plot vs number of model components
plot(pls::RMSEP(plsr.out,estimate=c("test"),newdata = val.plsr.data), main="MODEL RMSEP",
     xlab="Number of Components",ylab="Model Validation RMSEP",lty=1,col="black",cex=1.5,lwd=2)
box(lwd=2.2)

pls::R2(plsr.out, newdata = val.plsr.data)

##' @graph R2 plot vs number of model components
plot(pls::R2(plsr.out,estimate=c("test"),newdata = val.plsr.data), main="MODEL R2",
     xlab="Number of Components",ylab="Model Validation R2",lty=1,col="black",cex=1.5,lwd=2)
box(lwd=2.2)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## PLSR fit observed vs predicted plot data
#calibration

cal.plsr.output <- data.frame(cal.plsr.data[, which(names(cal.plsr.data) %notin% "Spectra")],
                              PLSR_Predicted=fit,
                              PLSR_CV_Predicted=as.vector(plsr.out$validation$pred[,,nComps]))
cal.plsr.output <- cal.plsr.output %>%
  mutate(PLSR_CV_Residuals = PLSR_CV_Predicted-get(inVar))
head(cal.plsr.output)

# R2 and RMSEP values
cal.R2 <- round(pls::R2(plsr.out,intercept=F)[[1]][nComps],2)
cal.RMSEP <- round(sqrt(mean(cal.plsr.output$PLSR_CV_Residuals^2)),2)

val.plsr.output <- data.frame(val.plsr.data[, which(names(val.plsr.data) %notin% "Spectra")],
                              PLSR_Predicted=as.vector(predict(plsr.out, 
                                                               newdata = val.plsr.data, 
                                                               ncomp=nComps, type="response")[,,1]))
val.plsr.output <- val.plsr.output %>%
  mutate(PLSR_Residuals = PLSR_Predicted-get(inVar))
head(val.plsr.output)


# Plots and Residuals
val.R2 <- round(pls::R2(plsr.out,newdata=val.plsr.data,intercept=F)[[1]][nComps],2)
val.RMSEP <- round(sqrt(mean(val.plsr.output$PLSR_Residuals^2)),2)

rng_quant <- quantile(cal.plsr.output[,inVar], probs = c(0.001, 0.999))
cal_scatter_plot <- ggplot(cal.plsr.output, aes(x=PLSR_CV_Predicted, y=get(inVar))) + 
  theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1, color="dark grey", 
                                          linetype="dashed", size=1.5) + xlim(rng_quant[1], 
                                                                              rng_quant[2]) + 
  ylim(rng_quant[1], rng_quant[2]) +
  labs(x=paste0("Predicted ", paste(inVar), ""),
       y=paste0("Observed ", paste(inVar), ""),
       title=paste0("Calibration: ", paste0("Rsq = ", cal.R2), "; ", paste0("RMSEP = ", 
                                                                            cal.RMSEP))) +
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

cal_resid_histogram <- ggplot(cal.plsr.output, aes(x=PLSR_CV_Residuals)) +
  geom_histogram(alpha=.5, position="identity") + 
  geom_vline(xintercept = 0, color="black", 
             linetype="dashed", size=1) + theme_bw() + 
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

rng_quant <- quantile(val.plsr.output[,inVar], probs = c(0.001, 0.999))
val_scatter_plot <- ggplot(val.plsr.output, aes(x=PLSR_Predicted, y=get(inVar))) + 
  theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1, color="dark grey", 
                                          linetype="dashed", size=1.5) + xlim(rng_quant[1], 
                                                                              rng_quant[2]) + 
  ylim(rng_quant[1], rng_quant[2]) +
  labs(x=paste0("Predicted ", paste(inVar), ""),
       y=paste0("Observed ", paste(inVar), ""),
       title=paste0("Validation: ", paste0("Rsq = ", val.R2), "; ", paste0("RMSEP = ", 
                                                                           val.RMSEP))) +
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

val_resid_histogram <- ggplot(val.plsr.output, aes(x=PLSR_Residuals)) +
  geom_histogram(alpha=.5, position="identity") + 
  geom_vline(xintercept = 0, color="black", 
             linetype="dashed", size=1) + theme_bw() + 
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

# plot cal/val side-by-side
scatterplots <- grid.arrange(cal_scatter_plot, val_scatter_plot, cal_resid_histogram, 
                             val_resid_histogram, nrow=2,ncol=2)