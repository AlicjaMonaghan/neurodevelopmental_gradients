# This script conducts repeated 10-fold cross-validation on generalised additive
# mixed models (GAMMs) to assess the relationship between structure-function
# coupling in manifold space and the working memory dimension of cognition. This
# addresses the "prediction vs association analysis" eLife R1 comment. Written 
# by Alicja Monaghan at the MRC CBU, June 2025.

rm(list = ls())
library(mgcv)
library(splitstackshape)
library(dplyr)
source('code/GAMM.functions.v3.R')

setwd('Desktop/neurodevelopmental_gradients/')
# Read in the network-level structure-function coupling estimates 
coupling_df = read.csv('data/structure.function/coupling.data.df.csv')
# Load in the cognitive dimension loadings
cognitive_coupling_df = read.csv(
  'data/structure.function/cognitive.participant.loadings.csv') %>%
  # Select all columns apart from age_in_months: these appeared to be rounded 
  # down to the closest age. To match the precise values given by the NKI meta-
  # data sheet, use scan_age from the coupling data frame.
  select(!age_in_months) %>%
  # Merge with the structure-function coupling data frame
  left_join(coupling_df, by = c("dataset", "id", "timepoint"))
# Specify in which networks we're going to investigate structure-function 
# coupling: Figure 5 in the eLife paper shows that coupling globally, and 
# within somato-motor, dorsal attention, and default-mode networks.
divisions = c("global", "SomMot", "DorsAttn", "Default")
# Conduct 5-fold cross-validation 50 times to derive a distribution of 
# prediction accuracies. 
nrep = 50
# Create an array to hold the prediction accuracy estimates 
pred_corr_array = array(dim=c(nrep))
for (i in 1:nrep){
  # Split data into train and test using stratified sampling based on group (NKI
  # or CALM) and scan age. Shuffle the data to remove associations with age.
  # cognitive_coupling_df$scan_age = sample(cognitive_coupling_df$Default)
  # split_data = stratified(indt = cognitive_coupling_df, group = c("dataset", "scan_age"), size=0.50, bothSets = TRUE)
  split_data = stratified(indt = cognitive_coupling_df, size=0.90, bothSets = TRUE, group=c("id"))
  training_data = split_data[["SAMP1"]]
  test_data = split_data[["SAMP2"]]
  # Within each significant division, use structure-function coupling to predict
  # the second cognitive dimension (working memory): Correct in manuscript that
  # we're using 4 degrees of freedom, not 3!
  coupling.gamm = fit.gamm.tensor.interaction.with.random.effects(
    df=training_data, outcome = "Default", knots = 4, fx = FALSE, parametric = c("sex", "meanfwd"),
    tensor.interaction = list(c("scan_age", "Factor1"), c("scan_age", "Factor2")))
    # Use the fitted GAM to predict coupling based on the test data
    coupling.prediction = predict.gam(coupling.gamm[[2]]$gam, newdata = test_data)
    # Correlate the predicted values with the empirical values!
    pred_corr = cor.test(coupling.prediction, test_data$Default, method = "spearman")
    pred_corr_array[i] = pred_corr$estimate
    print(sprintf('Finished rep %d for default coupling', i))
}



