# This script addresses sensitivity and exploratory analyses for
# structure-function coupling, from Reviewer 2 for "Canonical neurodevelopmental
# trajectories of structural and functional manifolds", in eLife (2025). Written
# by Alicja Monaghan, MRC Cognition and Brain Sciences Unit, July 2025.

rm(list = ls())
library(mgcv)

setwd('Desktop/neurodevelopmental_gradients/')
source('code/GAMM.functions.v3.R')

# Load the structure-function coupling estimates for neuro-typical CALM. 
calm_neurotypical_coupling = read.csv('data/calm/structure.function/calm_neurotypical_coupling.csv')
divisions = c("global", "Cont", "Default", "DorsAttn", "Limbic", "SalVentAttn", "SomMot", "Vis")
# Initialize an array to hold the p-value, F statistic, and EDF for the smooth
# scan age term! The final slice is for the FDR-corrected p-value.
scan_age_effects = array(NA, dim=c(length(divisions), 5))
# Use age, sex, and head motion (mean across structural and functional scans) to
# predict structure-function coupling. We're not including ID here, as all
# participants are from the baseline study.
for (division_idx in 1:length(divisions)){
  gam_formula = as.formula(sprintf('%s ~ s(scan_age, k=4) + meanFWD + Sex', divisions[division_idx]))
  gamm_output = gam(formula=gam_formula, data=calm_neurotypical_coupling)
  gamm_summary = summary.gam(gamm_output)
  # Collect the p-value, F statistic and EDF for the smooth scan age term
  scan_age_effects[division_idx, 1] = gamm_summary$s.table[4]
  scan_age_effects[division_idx, 2] = gamm_summary$s.table[3]
  scan_age_effects[division_idx, 3] = gamm_summary$s.table[2]
  # Get the adjusted R-squared for the whole model
  scan_age_effects[division_idx, 4] = gamm_summary$r.sq * 100
}
# Correct for multiple corrections by controlling the BH rate
scan_age_effects[,5] = p.adjust(scan_age_effects[,1], method = "BH")
# Convert into a data frame for reporting
scan_age_effects_df = data.frame(scan_age_effects, row.names = divisions)
colnames(scan_age_effects_df) = c("raw_p", "F", "EDF", "r_sq", "corrected_p")
write.csv(scan_age_effects_df, 'data/structure.function/calm_neurotypical_age_coupling.csv')
