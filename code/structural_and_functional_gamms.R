# This script runs generalized additive mixed models (GAMMs) for structural and 
# functional connectivity manifold eccentricity, using age as a smooth predictor
#, an interaction between age and data set as a factor-smooth interaction, and 
# mean frame-wise displacement, sex, and data set as categorical co-variates. We
# fit one model per modality, across both data sets. Written by Alicja 
# Monaghan, Alicja.Monaghan@mrc-cbu.cam.ac.uk

### PART 1 - Setting up the Work Space ####
rm(list = ls())
library(raveio)
library(gamm4)
library(ggplot2)
library(cowplot)
reticulate::py_run_string("import sys")
setwd('/Users/alicjamonaghan/Desktop/neurodevelopmental_gradients')
# Add the custom GAMM function code.
source('code/GAMM.functions.v3.R')
# Specify key DME variables
modalities = c("structural","functional")
datasets = c("calm","nki")
calm_timepoints = c("baseline","followup")
nnode = 200
ncomp = 3
# And load the region labels for Schaefer 200-node 7-network atlas.
schaefer200x7_metadata = read_mat("data/schaefer200x7_1mm_info.mat")
schaefer200x7_labels = unlist(schaefer200x7_metadata[["schaefer200x7.1mm.info"]][[1]])
nroi = length(schaefer200x7_labels)
# And set the labels for Yeo's (2011) resting-state functional connectivity networks
yeo_7_networks = c("Vis","SomMot","DorsAttn","SalVentAttn","Limbic","Cont","Default")
# Load the group-level and individual-level DME outputs for each data set and
# modality...
calm_dme_and_metadata_structural = 
  read.csv("data/calm/dme/dme.and.metadata.structural.connectivity.csv")[,-1]
calm_dme_and_metadata_functional = 
  read.csv("data/calm/dme/dme.and.metadata.functional.connectivity.csv")[,-1]
nki_dme_and_metadata_structural = 
  read.csv("data/nki/dme/dme.and.metadata.structural.connectivity.csv")[,-1]
nki_dme_and_metadata_functional = 
  read.csv("data/nki/dme/dme.and.metadata.functional.connectivity.csv")[,-1]

### PART 2 - Assessing Age Effects on Global and Network-Level Structural and Functional Manifold Eccentricity ####
ICN_modality_gamm_list = vector("list", 2)
for (modality_idx in 1:length(modalities)){
  # Initialize a list to hold the global and network-level plots
  modality_plot_list = vector("list", length = length(yeo_7_networks) + 1)
  modality = modalities[modality_idx]
  # Pull the meta-data and DME values for each data set - Initialize a vector
  # to hold the number of participants.
  nsub = array(NA,length(datasets))
  metadata_df_list = vector("list",length(datasets))
  for (dataset_idx in 1:length(datasets)){
    dataset = datasets[dataset_idx]
    metadata_df_list[[dataset_idx]] = get(paste0(dataset,'_dme_and_metadata_',modality))
    nsub[dataset_idx] = nrow(metadata_df_list[[dataset_idx]])
  }
  # Bind the rows of the two meta-data data frames
  metadata_df_list[[1]]$id = as.character(metadata_df_list[[1]]$id)
  dme_and_metadata = bind_rows(metadata_df_list)
  # Calculate network-level manifold eccentricity, for each of Yeo's 7 resting-
  # state functional connectivity networks
  dme_and_metadata_df = calculate.yeo.7.network.manifold.eccentricity(dme_and_metadata)
  # Add factor column for data sets
  dme_and_metadata_df$dataset = factor(c(rep("calm",times = nsub[1]),rep("nki",times = nsub[2])))
  # Save meta-data...
  write.csv(dme_and_metadata_df, sprintf('data/%s/calm_nki_dme_and_metadata_%s.csv', modality, modality), row.names = FALSE) 
  # Set the parametric, smooth, and interaction predictors to be used for all
  # GAMMs. 
  parametric = c("sex","meanfwd","dataset")
  smooth ="scan_age"
  interaction = c("scan_age","dataset")
  all_predictors = c(parametric,smooth,interaction)
  # Initialise lists to hold the p-values (FDR-un-corrected and corrected) for
  # all predictors across global and ICN-level analyses, as well as GAMM 
  # summary tables.
  manifold_eccentricity_pval = array(NA, dim=c(length(all_predictors)-1,length(yeo_7_networks)+1, 2))
  manifold_eccentricity_summary_tables = vector("list", length = length(yeo_7_networks)+1)
  manifold_eccentricity_gamm_output_list = vector("list",length(yeo_7_networks)+1)
  # Now run a GAMM to predict global manifold eccentricity!
  global_manifold_eccentricity_modality_gamm = 
    fit.gamm.with.random.effects(df = dme_and_metadata_df, parametric=parametric,
                                 smooth=smooth, interaction=interaction,
                                 outcome = "global_manifold_eccentricity", fx = FALSE, knots = 4,
                                 report_stats = TRUE)
  # Print the output summary table
  print(summary(global_manifold_eccentricity_modality_gamm[[2]]$gam))
  # Assign uncorrected p-values to the output array
  manifold_eccentricity_pval[,1,1] = global_manifold_eccentricity_modality_gamm[[1]][,3]
  # Assign the results table to the output 
  manifold_eccentricity_summary_tables[[1]] = global_manifold_eccentricity_modality_gamm[[1]]
  # And assign the models themselves to the output list
  manifold_eccentricity_gamm_output_list[[1]] = global_manifold_eccentricity_modality_gamm[[2]]
  # Visualize the normalized age effect
  manifold_partial_age_effect_plot = plot.manifold.eccentricity.partial.age.effect(
    model = global_manifold_eccentricity_modality_gamm[[2]], modality = modality,
    alpha = .95)
  # Assign to list
  modality_plot_list[[1]] = manifold_partial_age_effect_plot
  # Now run a GAMM to predict network-level manifold eccentricity, using the 
  # same co-variates as above. 
  for (network_idx in 1:length(yeo_7_networks)){
    network_level_manifold_eccentricity_gamm = fit.gamm.with.random.effects(
      df=dme_and_metadata_df, parametric=parametric, smooth=smooth,
      interaction=interaction, outcome=yeo_7_networks[network_idx], fx=FALSE,
      knots=4, report_stats = FALSE)
    # Assign un-corrected p-values for each predictor to output array. Note that
    # we add one to the network_idx to account for how we assigned the global
    # results to the first index.
    manifold_eccentricity_pval[,network_idx+1,1] = network_level_manifold_eccentricity_gamm[[1]][,3]
    # Assign the results table summary to the output
    manifold_eccentricity_summary_tables[[network_idx + 1]] = network_level_manifold_eccentricity_gamm[[1]]
    # And assign the models themselves to the output list
    manifold_eccentricity_gamm_output_list[[network_idx + 1]] = network_level_manifold_eccentricity_gamm[[2]]
  }
  # Correct for multiple comparisons for network-level manifold eccentricity 
  # GAMMs using false-discovery rate. We subtract 1 from the length of 
  # all_predictors because the interaction is provided as two separate terms.
  manifold_eccentricity_pval[,,2] = p.adjust(manifold_eccentricity_pval[,,1], method="fdr")
  # Visualize the normalized age x data set interaction for each ICN. Again, we
  # add one to the network_idx to account for how we assigned the global GAMM
  # outputs to the first index.
  for (network_idx in 1:length(yeo_7_networks)){
    network_partial_age_effect_plot = 
      plot.manifold.eccentricity.partial.age.effect(
        model=manifold_eccentricity_gamm_output_list[[network_idx + 1]],
        modality = modality, alpha=.05)
    # Assign to output list
    modality_plot_list[[1 + network_idx]] = network_partial_age_effect_plot
  }
  # For each network, create a summary table of results. Ensure to include the 
  # FDR-corrected p-values for the main effects and FDR-corrected bootstrapped 
  # p-values for the interaction term. We're also including global-level 
  # analyses, hence adding 1 to the number of Yeo networks.
  for (idx in 1:(length(yeo_7_networks)+1)){
    # Select the correct network summary table
    summary_table = manifold_eccentricity_summary_tables[[idx]]
    # The parametric predictors are reported first. Subset by these, and format.
    parametric_coefs = summary_table[1:length(parametric),]
    # Find the indices of parametric coefficients in all_predictors.
    parametric_idx = which(all_predictors == parametric)[c(1:length(parametric))]
    # Replace the final column with the FDR-corrected p-values
    parametric_coefs$X3 = manifold_eccentricity_pval[parametric_idx,idx,2]
    colnames(parametric_coefs) = c("estimate", "t", "FDR-corrected p")
    if (idx == 1){
      print(sprintf('GAMM with %s global eccentricity as an outcome...', modality))
    } else{
      print(sprintf('GAMM with %s eccentricity in the %s network as an outcome...',
                    modality, yeo_7_networks[idx-1]))
    }
    print(parametric_coefs)
    # Now find the smooths!
    smooth_coefs = summary_table[(length(parametric) + 1): nrow(summary_table), ]
    # Replace the smooth p-value from the GAMM with the FDR-corrected p-value
    smooth_coefs[['X3']] = manifold_eccentricity_pval[4:5, idx, 2]
    colnames(smooth_coefs) = c("edf", "F", "FDR-corrected p")
    print(smooth_coefs)
    # Bind rows of the parametric and smooth coefficients, and then replace the 
    # previous summary table. To ensure smooth binding, temporarily assign the 
    # column names of parametric_coefs to smooth_coefs, and then rename 
    # appropriately.
    colnames(parametric_coefs) = colnames(smooth_coefs)
    df = bind_rows(parametric_coefs, smooth_coefs)
    colnames(df) = c("estimate/edf", "t/F", "FDR-corrected p")
    manifold_eccentricity_summary_tables[[idx]] = df
    # Save the summary data frames!
    if (idx == 1){
      save_filename = sprintf('data/%s/global.eccentricity.gamm.csv', modality)
    } else{
      save_filename = sprintf('data/%s/%s.eccentricity.gamm.csv', modality, yeo_7_networks[idx-1])
    }
    write.csv(df, file = save_filename)
  }
  # Loop through each network again and global analysis again...
  all_analysis_levels = c("global", yeo_7_networks)
  for (idx in 1:length(all_analysis_levels)){
    print(sprintf('Extracting %s %s statistics', modality, all_analysis_levels[idx]))
    # Extract the GAMM summary table for this network
    summary_table = manifold_eccentricity_summary_tables[[idx]]
    # Extract the p value for the main age effect 
    main_age_effect = summary_table['s(scan_age)', 'FDR-corrected p']
    # If the age x data set interaction or age main effects are significant, 
    # plot significant age derivatives
    if (manifold_eccentricity_pval[5,idx,2] < .05 | main_age_effect < .05){
      # Visualize the significant first age derivatives for each data set
      calm_significant_age_deriv = plot.sensitive.age.developmental.period(
        model = manifold_eccentricity_gamm_output_list[[idx]], tensor=FALSE, dataset = "calm")
      nki_significant_age_deriv = plot.sensitive.age.developmental.period(
        model = manifold_eccentricity_gamm_output_list[[idx]], tensor=FALSE, dataset = "nki")
      # Create list
      significant_age_deriv = list(calm_significant_age_deriv, nki_significant_age_deriv)
      # We'll now combine the age trajectory and derivative plots together. 
      # Check if there are any periods of significant change.
      if (anyNA(significant_age_deriv) == TRUE){
        # Find the index of the missing value
        missing_idx = which(is.na(significant_age_deriv))
        if (missing_idx == 1){
          # This is when CALM age derivatives aren't statistically significant
          allplots = list(modality_plot_list[[idx]], NULL, significant_age_deriv[[2]])
        } else if (missing_idx == 2){
          # When NKI age derivatives aren't statistically significant
          allplots = list(modality_plot_list[[idx]], significant_age_deriv[[1]], NULL)
        } else {
          # When both CALM and NKI are not statistically significant
          allplots = list(modality_plot_list[[idx]], NULL, NULL)
        }} else{
    allplots = list(modality_plot_list[[idx]], significant_age_deriv[[1]], significant_age_deriv[[2]])}} else{
      # If the age-dataset interaction and main age effects are both statistically
      # non-significant, then assign NULLs in their places.
      allplots = list(modality_plot_list[[idx]], NULL, NULL)
    } 
    combined_plot = plot_grid(plotlist = allplots, align = "v", axis = "lr", ncol = 1, rel_heights = c(16, 2, 2))
    ggsave(filename = sprintf('visualisation/developmental.trends/%s.%s.manifold.eccentricity.png',modality, all_analysis_levels[idx]),
           combined_plot, height = 8, width = 10, dpi = 700)
}}

