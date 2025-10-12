# This script conducts repeated 5-fold cross-validation on generalized additive
# mixed models (GAMMs) to assess the relationship between structure-function
# coupling in manifold space and the working memory dimension of cognition. This
# addresses the "prediction vs association analysis" eLife R1 comment. Written 
# by Alicja Monaghan at the MRC CBU, July-August 2025.

rm(list = ls())
library(splitstackshape)
library(dplyr)
library(rsample)
library(doParallel)
library(Metrics)
library(effsize)
library(mgcv)
library(stats)
library(corrr)
setwd('/Users/alicjamonaghan/Desktop/neurodevelopmental_gradients/')

# Find the number of cores we have and set them up for parallelization
cluster = makeCluster(detectCores() - 1)
registerDoParallel(cluster)

# This function will run the K-fold cross-validation, using split_df which is 
# the stratified K-fold data frame produced by the group_vfold_cv function in
# the rsample package. It also incorporates a null-testing procedure, where we 
# shuffle the second cognitive factor according to the number of observations
# each participant has. 
run_kfold_cv = function(split_df, outcome_var){
  # Initialize an array for the correlations between the empirical coupling in
  # the test data and predicted coupling from the model fit on the training data
  nfold = dim(split_df)[1]
  corrs = numeric(nfold)
  for (fold_idx in 1:nfold){
    train_data = analysis(split_df$splits[[fold_idx]])
    test_data = assessment(split_df$splits[[fold_idx]])
    # Set 3 knots to reduce model complexity
    predictors = "ti(scan_age, k=3) + ti(Factor2, k=3) + ti(scan_age, Factor2, k=3) + 
    ti(scan_age, Factor1, k=3) + ti(Factor1, k=3) + meanfwd + sex + s(id, bs='re')"
    gam_formula = as.formula(paste(outcome_var, "~", predictors))
    bam_output = bam(gam_formula, data=train_data, select=TRUE, discrete=TRUE)
    # Note, sometimes a warning will appear which flags that IDs in the training
    # data set are not present in the prediction dataset: this is normal as we
    # specifically grouped by ID so that each participant's data is restricted 
    # to either the train or testing data. 
    preds = predict.bam(bam_output, newdata=test_data)
    corrs[fold_idx] = cor(preds, test_data[[outcome_var]], method="spearman")
  }
  return(list(corrs))
}

set.seed(123)
# Read in the network-level structure-function coupling estimates 
coupling_df = read.csv('data/structure.function/coupling.data.df.csv')
# Load in the cognitive dimension loadings
cognitive_coupling_df = read.csv(
  'data/phenotypic/cognitive.participant.loadings.august25.csv') %>%
  # Merge with the structure-function coupling data frame
  left_join(coupling_df, by = c("dataset", "id", "timepoint")) %>%
  # Specify id as a factor
  mutate(id = as.factor(id))
# To examine the degree of collinearity between networks, conduct pairwise 
# Spearman correlations between all network-level coupling!
coupling_pairwise_corr = correlate(cognitive_coupling_df[,11:17], method="spearman")
# Specify which networks we're going to investigate: these are where 
# structure-function coupling is significantly predicted by cognition.
divisions = c("global", "Vis", "SomMot", "DorsAttn", "Default")
# Initialize an array which will hold the empirical and null prediction 
# accuracies for each division: we'll visualize these later on. Else, these can
# be loaded (with the same names in the data/coupling.phenotype sub-directory). 
n_rep = 100
kfold = 5
empirical_prediction_accuracies_all = array(NA, dim=c(n_rep*kfold, length(divisions)))
null_prediction_accuracies_all = array(NA, dim=c(n_rep*kfold, length(divisions)))
cohen_effsize_array = array(NA, dim=length(divisions))

for (division_idx in 1:length(divisions)){
  # Create cross-validation 5-fold splits for the data, where we group based on ID
  # and data-set, so that all observations related to a specific participant are
  # constrained to a single fold, to prevent data leakage, and the relative
  # proportions of CALM and NKI are similar. This code takes around 5 minutes to
  # run for each network-specific GAMM. 
  empirical_results_list = foreach(i=1:n_rep, .combine='c', .packages=c('rsample', 'mgcv', 'Metrics')) %dopar% {
    split_df = group_vfold_cv(data=cognitive_coupling_df, group="id",v=5, strata = "dataset")
    run_kfold_cv(split_df, outcome_var=divisions[division_idx])
  }
  # Every third entry in empirical_results_list is the prediction accuracy, 
  # starting with 1, whilst every even entry is the MAE. Collate these two!
  empirical_prediction_accuracy = unlist(empirical_results_list)
  empirical_prediction_accuracies_all[,division_idx] = empirical_prediction_accuracy
  print(sprintf('Empirical prediction accuracies for %s coupling computed', divisions[division_idx]))
  # To compare empirical prediction accuracies with what we'd expect by chance,
  # we shall permute the coupling data frame. As we're assessing the 
  # relationship between structure-function coupling and cognition, we shall
  # permute Factor2. We will group by participant ID and then shuffle!
  null_results_list = foreach(i=1:n_rep, .combine='c', .packages=c('rsample', 'mgcv', 'Metrics', 'dplyr')) %dopar% {
    # Find how many observations each participant has, and keep track of the 
    # order of each observation within each participant block using obs_index.
    cognitive_coupling_df_with_obs = cognitive_coupling_df %>% group_by(id) %>%
      mutate(obs_count = n(), obs_index = row_number()) %>% ungroup()
    # Shuffle participant IDs within each group defined by the number of 
    # observations (e.g. 320 with one observation, 105 with 2 observations, and
    # 33 with 3 observations). The shuffled_ids data frame will have 3 columns:
    # id, obs_count, and shuffled_id
    shuffled_ids = cognitive_coupling_df_with_obs %>% distinct(id, obs_count) %>%
      group_by(obs_count) %>% mutate(shuffled_id = sample(id)) %>% ungroup()
    # Join the shuffled IDs back to original data and replace Factor 2
    cognitive_coupling_df_shuffled = cognitive_coupling_df_with_obs %>%
      left_join(shuffled_ids, by = "id") %>%
      left_join(cognitive_coupling_df_with_obs %>% select(id, obs_index, Factor2) %>%
          rename(shuffled_id = id, Factor2_shuffled = Factor2),
        by = c("shuffled_id", "obs_index")) %>%
      mutate(Factor2 = Factor2_shuffled)
    # Now conduct the stratified 5-fold cross-validation
    split_df = group_vfold_cv(data=cognitive_coupling_df_shuffled, group="id",v=5, strata = "dataset")
    run_kfold_cv(split_df, outcome_var=divisions[division_idx])
  }
  null_prediction_accuracy = unlist(null_results_list)
  null_prediction_accuracies_all[,division_idx] = null_prediction_accuracy
  print(sprintf('Null prediction accuracies for %s coupling computed', divisions[division_idx]))
  # Report the mean and standard deviation of empirical prediction accuracy
  mean_empirical = mean(empirical_prediction_accuracy)
  std_empirical = sd(empirical_prediction_accuracy)
  print(sprintf('Empirical prediction accuracy for %s is %.3f +/- %.3f', divisions[division_idx], mean_empirical, std_empirical))
  # Test the null hypothesis that there is no significant difference in
  # prediction accuracies between the empirical and permuted distributions.
  prediction_accuracy_ttest = t.test(empirical_prediction_accuracy, null_prediction_accuracy)
  print(sprintf('Two-sample t-test of empirical (M = %.3f, SD = %.3f) versus null (M = %.3f, SD = %.3f) prediction accuracies',
                mean_empirical, std_empirical, mean(null_prediction_accuracy), sd(null_prediction_accuracy)))
  print(sprintf('T-statistic of %.3f, with %.3f degrees of freedom, and p of %.3f',
                prediction_accuracy_ttest$statistic, prediction_accuracy_ttest$parameter,
                prediction_accuracy_ttest$p.value))
  # Compute the effect size for the difference between the distributions!
  cohen_effsize = cohen.d(empirical_prediction_accuracy, null_prediction_accuracy)
  cohen_effsize_array[division_idx] = cohen_effsize$estimate
  print(sprintf('Cohen D of %.3f' , cohen_effsize$estimate))
}

write.table(empirical_prediction_accuracies_all, file='data/coupling.phenotypic/empirical_prediction_accuracies_all.csv',
          col.names = divisions, row.names = FALSE, sep = ",", quote=FALSE)
write.table(null_prediction_accuracies_all, file='data/coupling.phenotypic/null_prediction_accuracies_all.csv',
            col.names = divisions, row.names = FALSE, sep=",", quote=FALSE)
write.table(cohen_effsize_array, file = "data/coupling.phenotypic/cohen_effsize_array.csv",
            row.names = divisions, sep = ",", quote=FALSE)

# For the empirical prediction accuracies, order the networks by decreasing mean
# accuracy and conduct a series of independent sample t-tests to test for
# significant differences in mean accuracy. By ordering the networks by mean, we
# reduce the number of multiple comparisons (FDR).
empirical_prediction_accuracy_df = data.frame(empirical_prediction_accuracies_all)
colnames(empirical_prediction_accuracy_df) = divisions
empirical_prediction_accuracy_df = empirical_prediction_accuracy_df[
  , names(sort(colMeans(empirical_prediction_accuracy_df), decreasing=TRUE))]
raw_pvalues = array(NA, dim=length(divisions)-1)
for (i in 1:(length(divisions)-1)){
  network_comparison = t.test(empirical_prediction_accuracy_df[,i], empirical_prediction_accuracy_df[,i+1])
  network_effsize = cohen.d(empirical_prediction_accuracy_df[,i], empirical_prediction_accuracy_df[,i+1])
  print(sprintf('Comparison between %s and %s: t-test statistic of %.3f, degrees of freedom as %.3f and d of %.3f',
                colnames(empirical_prediction_accuracy_df)[i], colnames(empirical_prediction_accuracy_df)[i+1],
                network_comparison$statistic, network_comparison$parameter, network_effsize$estimate))
  raw_pvalues[i] = network_comparison$p.value
}
# Correct for multiple comparisons!
fdr_pvalues = p.adjust(raw_pvalues, method="BH")
# Test that the distribution of empirical prediction accuracies is significantly
# higher than 0 using a series of one-tailed t-tests.
for (i in 1:length(divisions)){
  ttest_result = t.test(empirical_prediction_accuracy_df[,i], mu=0, alternative = "greater")
  print(sprintf("One-sample t-test of %s empirical prediction accuracies has p-value of %.3f",
                divisions[i], ttest_result$p.value))
}


