##' Recreate estimation of chlorophyll content using leafDB
##' @author Nathan Malamud
##' @date 2024.09.27
##' 

library(tidyverse)
library(prospect)

# download ANGERS dataset
LeafDB <- download_LeafDB(dbName = 'ANGERS')
# Prior estimation of N using R only
Nprior_R <- Get_Nprior(lambda = LeafDB$lambda, 
                       Refl = LeafDB$Refl)
# Prior estimation of N using T only
Nprior_T <- Get_Nprior(lambda = LeafDB$lambda, 
                       Tran = LeafDB$Tran)

# Estimate all parameters for PROSPECT-D
Parms2Estimate  <- 'ALL'
InitValues <- data.frame(CHL = 40, CAR = 10, ANT = 0.1, BROWN = 0, 
                         EWT = 0.01, LMA = 0.01, N = 1.5)

# Adjust spectral domain for SpecPROSPECT to fit leaf optical properties 
SubData <- FitSpectralData(lambda = LeafDB$lambda,
                           Refl = LeafDB$Refl, 
                           Tran = LeafDB$Tran)
print('PROSPECT inversion using full spectral range')
res_Ronly <- Invert_PROSPECT(SpecPROSPECT = SubData$SpecPROSPECT, 
                             Refl = SubData$Refl, Tran = NULL, 
                             PROSPECT_version = 'D', 
                             Parms2Estimate = Parms2Estimate, 
                             InitValues = InitValues)

res_Tonly <- Invert_PROSPECT(SpecPROSPECT = SubData$SpecPROSPECT, 
                             Refl = NULL, Tran = SubData$Refl, 
                             PROSPECT_version = 'D',
                             Parms2Estimate = Parms2Estimate, 
                             InitValues = InitValues)

