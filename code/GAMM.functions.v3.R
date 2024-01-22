# This function fits generalized additive mixed models (GAMMs), using the gamm4
# package, to structural and functional eccentricities from two data sets: the 
# Center for Attention, Learning, and Memory (CALM), and the Nathan Kline Rock-
# land Institute (NKI). 
rm(list = ls())
library(dplyr)
library(gratia)
library(itsadug)
library(scales)
library(pbkrtest)
library(rlang)
library(parallel)
library(gamm4)
library(mgcv)
reticulate::use_miniconda('r-reticulate')
reticulate::py_run_string("import sys")
setwd('U:/gradients_open_access')
Sys.setenv(`_R_USE_PIPEBIND_` = TRUE)

# This function fits GAMMs with subject-specific intercepts as random effects. 
# The output are the statistics and p-values for parametric coefficients and 
# smooths, as well as the GAMM summary. We also provide a written summary of 
# each GAMM's results table. For parametric coefficients, we report the 
# estimate, t statistic, and p-value. For smooth terms, we report the effective
# degrees of freedom, F statistic, and p value.
fit.gamm.with.random.effects = function(df, parametric, smooth, outcome, interaction, knots, fx, report_stats){
  # Order each discrete parametric predictor
  for (predictor_idx in 1:length(parametric)){
    if (is.null(levels(df[[parametric[predictor_idx]]])) == FALSE){
      df[[parametric[predictor_idx]]] = ordered(df[[parametric[predictor_idx]]])
    }
  }
  # Also check if there are any categorical variables as interactions
  for (interaction_idx in 1:length(interaction)){
    if (is.null(levels(df[[interaction[interaction_idx]]])) == FALSE){
      df[[interaction[interaction_idx]]] = ordered(df[[interaction[interaction_idx]]])
    }
  }
  # Create the gamm formula. Note that we specify random effects and data set at
  # a later stage.
  gamm4_formula = as.formula(sprintf("%s ~ s(%s, k=%f, fx=%s) + %s + s(%s, by=%s, k=%f, fx=%s)",
                                     outcome,smooth,knots,fx,paste(parametric,collapse=" + "),interaction[1],interaction[2],knots,fx))
  # Fit the GAMM and extract summary
  gamm_output = gamm4(formula=gamm4_formula, random = ~(1|id), data = df)
  gamm_summary = summary(gamm_output$gam)
  # Initialize an output array for test statistics. We collect estimates and t
  # statistics for parametric predictors, and estimated degrees of freedom, 
  # F-statistic and p-values for smooths. We subtract one from the length of the
  # interaction input to account for how both predictors in the interaction term
  # are incorporated into a single line in the output table which we'll be 
  # indexing.
  output_array = array(NA,dim=c(length(parametric)+length(smooth)+length(interaction)-1,3))
  # Note that the ordered parametric predictors have different names now, to 
  # account for the levels of the ordering e.g. sex is now sexM. 
  parametric_names = rownames(gamm_summary$p.table)
  # Remove intercept from parametric_names.
  parametric_names = parametric_names[!parametric_names %in% "(Intercept)"]
  # Extract summary statistics for parametric predictors!
  for (predictor_idx in 1:length(parametric_names)){
    output_array[predictor_idx,1] = gamm_summary$p.table[c(parametric_names[predictor_idx]), "Estimate"]
    output_array[predictor_idx,2] = gamm_summary$p.table[c(parametric_names[predictor_idx]),"t value"]
    output_array[predictor_idx,3] = gamm_summary$p.table[c(parametric_names[predictor_idx]),"Pr(>|t|)"]
    if (report_stats == TRUE){
      # Report back to user!
      print(sprintf('Effect of %s: estimate of %.3f, t-statistic of %.3f, FDR uncorrected p-value of %.3f',
                    parametric_names[predictor_idx], round(output_array[predictor_idx, 1],3),
                    round(output_array[predictor_idx, 2],3), round(output_array[predictor_idx, 3],3)))
    }
  }
  # To extract summary statistics for the smooths, note that these names will 
  # also have changed if we included an interaction term with a factor, as that
  # factor will've been ordered, and thus its name changed. 
  smooth_names = rownames(gamm_summary$s.table)
  for (predictor_idx in 1:length(smooth_names)){
    output_array[predictor_idx+length(parametric_names),1] = gamm_summary$s.table[c(smooth_names[predictor_idx]), "edf"]
    output_array[predictor_idx+length(parametric_names),2] = gamm_summary$s.table[c(smooth_names[predictor_idx]),"F"]
    output_array[predictor_idx+length(parametric_names),3] = gamm_summary$s.table[c(smooth_names[predictor_idx]),"p-value"]
    if (report_stats == TRUE){
      print(sprintf('Effect of %s: estimate of %.3f, t-statistic of %.3f, FDR-uncorrected p-value of %.3f',
                    smooth_names[predictor_idx], round(output_array[predictor_idx + length(parametric_names), 1], 3),
                    round(output_array[predictor_idx + length(parametric_names), 2], 3),
                    round(output_array[predictor_idx + length(parametric_names), 3], 3)))
    }
  }
  # Format output array as a data frame
  output_df = data.frame(output_array)
  rownames(output_df) = c(parametric_names, smooth_names)
  return(list(output_df, gamm_output))
  # Check concurvity (~collinearity) of smooth effects in the model. This looks
  # at whether the smooth effects can be approximated by one or more other 
  # smooth term in the model. --> CHECK WITH DUNCAN, AS AGE AND SEX HAVE AN 
  # OBSERVED CONCURVITY OF ~.67
}
# This function fits GAMMs using tensors for main effects and interactions, when
# the interactions are between two continuous variables. Predictors not included
# in the continuous*continuous interaction are specified as parametric or smooth
# terms, like in the function above. The output are the statistics and p-values
# for parametric terms and smooth/tensor terms. The function accepts one tensor
# interaction. 
fit.gamm.tensor.interaction.with.random.effects = function(df, tensor.interaction, outcome, knots, fx, smooth, parametric){
  # Order each discrete parametric predictor
  for (predictor_idx in 1:length(parametric)){
    if (is.null(levels(df[[parametric[predictor_idx]]])) == FALSE){
      df[[parametric[predictor_idx]]] = ordered(df[[parametric[predictor_idx]]])
    }
  }
  # Specify the part of the GAMM formula which deals with parametric terms.
  nparametric = length(parametric)
  if (nparametric > 1){
    parametric.specification = paste0(parametric, collapse = " + ")
  } else{
    parametric.specification = parametric
  }
  # Specify the part of the GAMM formula which deals with tensor interactions. 
  # These interactions are inputted as a list.
  ninteraction = length(tensor.interaction)
  if (ninteraction > 1){
    # Initialize list to hold the formatted interactions
    interaction_list = vector("list", length = ninteraction)
    for (interaction_idx in 1:ninteraction){
      interaction = tensor.interaction[[interaction_idx]]
      # Note that we include a main effect tensor for the variable which is not
      # unique to each interaction e.g. for two age x Factor interactions, 
      # specify the Fcator main effect and the whole interaction.
      interaction_list[[interaction_idx]] = 
        sprintf('ti(%s, k=%f, fx=%s) + ti(%s, %s, k=%f, fx=%s)',
                interaction[[2]], knots, fx, interaction[[1]], interaction[[2]], knots, fx)
    }
    # Now specify the main effect of the common predictor across interactions
    common.main.effect = sprintf('ti(%s, k=%f, fx=%s)', interaction[[1]], knots, fx)
    interaction.specification = paste0(c(paste0(interaction_list, collapse = " + "), common.main.effect), collapse = " + ")
  } else{
    interaction.specification = 
      sprintf('ti(%s, k=%f, fx=%s) + ti(%s, k=%f, fx=%s) + ti(%s, %s, k=%f, fx=%s)',
              interaction[[1]], knots, fx, interaction[[2]], knots, fx, 
              interaction[[1]], interaction[[2]], knots, fx)
  }
  # Specify the GAMM formula. Note that we have to use the gamm function from
  # the MGCV package, not the gamm4 package.
  gamm.formula = as.formula(sprintf('%s ~ %s + %s', outcome, interaction.specification, parametric.specification)) 
  # Fit the GAMM and extract summary. Note that adding a control parameter can
  # overcome convergence issues if they exist --> control=lmeControl(
  # msMaxIter = 1000, msMaxEval = 1000)
  gamm_output = gamm(formula=gamm.formula, random = list(id=~1), data = df, method = "REML")
  gamm_summary = summary(gamm_output$gam)
  # Create an output array of test statistics for the parametric and smooth/
  # tensor terms. 
  output_array = array(NA,dim=c(length(parametric)+length(unlist(tensor.interaction))+1,3))
  # Note that the ordered parametric predictors have different names now, to 
  # account for the levels of the ordering e.g. sex is now sexM. 
  parametric_names = rownames(gamm_summary$p.table)
  # Remove intercept from parametric_names.
  parametric_names = parametric_names[!parametric_names %in% "(Intercept)"]
  # Extract summary statistics for parametric predictors!
  for (predictor_idx in 1:length(parametric_names)){
    output_array[predictor_idx,1] = gamm_summary$p.table[c(parametric_names[predictor_idx]), "Estimate"]
    output_array[predictor_idx,2] = gamm_summary$p.table[c(parametric_names[predictor_idx]),"t value"]
    output_array[predictor_idx,3] = gamm_summary$p.table[c(parametric_names[predictor_idx]),"Pr(>|t|)"]
  }
  # Extract summary statistics for smooths/tensors.
  smooth_names = rownames(gamm_summary$s.table)
  for (predictor_idx in 1:length(smooth_names)){
    output_array[predictor_idx+length(parametric_names),1] = gamm_summary$s.table[c(smooth_names[predictor_idx]), "edf"]
    output_array[predictor_idx+length(parametric_names),2] = gamm_summary$s.table[c(smooth_names[predictor_idx]),"F"]
    output_array[predictor_idx+length(parametric_names),3] = gamm_summary$s.table[c(smooth_names[predictor_idx]),"p-value"]
  }
  # Format output array as a data frame
  output_df = data.frame(output_array)
  rownames(output_df) = c(parametric_names, smooth_names)
  return(list(output_df, gamm_output))
}
# This function generates predicted data for a GAMM for variables X, whilst 
# holding all other variables constant, which can then be used to visualize 
# partial effects for X and to extract confidence intervals.
gamm.predictions = function(model, predictors, probability, rm.ranef, control.continuous){
  # First argument is the GAMM model, including both GAM and LME parts. Extract
  # the data used to fit the GAMM, and extract the predictors and the continuous
  # control variables. Note that we extract the continuous control variables 
  # separately from the discrete controls, as the prediction algorithm will 
  # automatically choose a level. This will be the same across both data sets, 
  # and therefore won't make a difference to predictions, whilst continuous
  # control values may differ between data sets/ categorical predictors.
  df = model$gam$model %>% dplyr::select(all_of(c(predictors, control.continuous)))
  # Using the itsadug package, create a prediction data frame comprised of 
  # the continuous predictor, such as age, from smallest to largest. We want to 
  # predict our outcome, such as global manifold eccentricity using predictors 
  # of interest only, whilst keeping all others constant, in both the LME and 
  # GAM parts of the GAMM. Thus, loop through each predictor. 
  cond = vector("list", length = length(c(predictors, control.continuous)))
  # Check if all predictors are continuous, factors, or a combination. If all
  # predictors are continuous, then predict using the minumum and maximum. 
  if (!("factor" %in% unlist(lapply(df,class))) == TRUE){
    for (predictor_idx in 1:length(predictors)){
      # If there's a continuous predictor and no factor, then predict data using 
      # the minimum and maximum, with number of repetitions being equal to the 
      # length of the original data frame.
      predictor_name = predictors[predictor_idx]
      cond[[predictor_idx]] = seq(min(df[[predictor_name]]),max(df[[predictor_name]]),length.out=nrow(df))
    }
    # As we won't be predicting based off of factors, the continuous control 
    # predictor will just be the median value.
    median.control.continuous = df %>% select(all_of(control.continuous)) %>% 
      apply(., 2, median) %>% as.numeric()
    for (control.idx in 1:length(control.continuous)){
      cond[[control.idx + 1]] = median.control.continuous[control.idx]
    }
    # Name the lists 
    names(cond) = c(predictors, control.continuous)
  } 
  # If all predictors are factors...
  if (!("numeric") %in% unlist(lapply(df,class)) == TRUE){
    # For each categorical variable, find the frequency of each level
    for (predictor_idx in 1:length(predictors)){
      predictor_name = predictors[predictor_idx]
      levels_frequency = summary(df[[predictor_name]])
      # And attach the levels to the condition list
      cond[[predictor_idx]] = levels(df[[predictor_name]])
    }
    # Set the continuous control predictor to the median of the entire data set
    median.control.continuous = df %>% select(all_of(control.continuous)) %>% 
      apply(., 2, median) %>% as.numeric()
    # Add this to the cond list, and name accordingly
    for (control.idx in 1:length(control.continuous)){
      cond[[control.idx + 1]] = median.control.continuous[control.idx]
    }
    # Name each of the factors in cond
    names(cond) = c(predictors, median.control.continuous)
  }
  # If there's a combination of continuous and discrete variables, then create a
  # prediction array of the continuous variable for each level of each discrete 
  # variable.
  if ("numeric" %in% unlist(lapply(df,class)) == TRUE & "factor" %in% unlist(lapply(df,class)) == TRUE){
    # First, find the categorical predictor. 
    categorical_predictor_idx = which(sapply(df,is.factor))[[1]]
    # Find the frequency of each level and add the level names to the cond list
    levels_frequency = summary(df[,categorical_predictor_idx])
    cond[[categorical_predictor_idx]] = levels(df[,categorical_predictor_idx])
    # Now find the indices of the continuous predictor
    continuous_predictor_idx = which(sapply(df[predictors],is.numeric))
    # For each level of the factor, extract the associated values for the 
    # continuous predictor. 
    continuous_Values_for_factor_list = vector("list",length(levels_frequency))
    # Also initialize a vector to hold the median values of the continuous 
    # control variable for each level of the factor.
    continuous_control_array = array(NA, dim = length(levels_frequency))
    # Loop across each factor level...
    for (factor_level_idx in 1:length(levels_frequency)){
      factor_level = levels(df[,categorical_predictor_idx])[factor_level_idx]
      continuous_values_for_factor = df[which(df[,categorical_predictor_idx] == factor_level),continuous_predictor_idx]
      # As we're splitting the data frame by a categorical predictor, we'll need
      # different values of the continuous control predictor. First, subset df
      # by the factor level we're interested in.
      df_subset = df[which(df[,predictors[categorical_predictor_idx]] == factor_level),]
      # Find the median value of the continuous predictor, and assign to output
      continuous_control_array[factor_level_idx] = median(df_subset[,control.continuous])
      # Create a prediction matrix ranging from minimum to maximum of this 
      # array, with length equal to the number of instances of that factor level.
      # Combine this with as a list with each factor level.
      continuous_Values_for_factor_list[[factor_level_idx]] =
        # Adding the factor level...
        list(factor_level,
             # Adding the continuous values for this factor level
             seq(min(continuous_values_for_factor),max(continuous_values_for_factor),
                 length.out=levels_frequency[[factor_level_idx]]),
             # Adding the median continuous control value for this factor level
             continuous_control_array[factor_level_idx])
      # Assign labels for each variable in this list
      names(continuous_Values_for_factor_list[[factor_level_idx]]) = 
        c(predictors[categorical_predictor_idx], predictors[continuous_predictor_idx],
          control.continuous)
    }
    # Combine all new data into a single matrix!
    continuous_values = vector("list", 2)
    for (idx in 1:length(levels_frequency)){
      continuous_values[[idx]] = continuous_Values_for_factor_list[[idx]]$scan_age
    }
    newdata = unlist(continuous_values)
  }
  # Now create the prediction matrices! 
  # probabilty dictates the % confidence intervals to
  if (probability == .95){
    f = 1.96
  } else if (probability == .99){
    f = 2.58
  }
  # remove_random_effects is a boolean argument
  # If the predictors are a combination of categorical and continuous variables,
  # then produce one prediction matrix per factor level.
  if ("numeric" %in% unlist(lapply(df,class)) == TRUE & "factor" %in% unlist(lapply(df,class)) == TRUE){
    categorical_predictor_idx = which(sapply(df,is.factor))[[1]]
    categorical_levels = levels(df[,categorical_predictor_idx])
    predictions_list = vector("list",length(categorical_levels))
    for (categorical_level_idx in 1:length(categorical_levels)){
      predictions_list[[categorical_level_idx]] = 
        get_predictions(model=model$gam, cond=continuous_Values_for_factor_list[[categorical_level_idx]],rm.ranef=rm.ranef,f=f, se=TRUE)
    }
    # Combine the two data frames of predicted data, one for each data set, into
    # a single data frame.
    predictions = bind_rows(predictions_list)
  } else {
    # If all predictors are either categorical or continuous, we only create one 
    # prediction matrix. The new categorical data is the same as the input. The 
    # new continuous data spans the minimum and maximum with N number of equally-
    # spaced breaks, where N is the number of rows of the input data frame. 
    predictions = get_predictions(model=model$gam, cond=cond, rm.ranef=rm.ranef, f=f, se=TRUE)
  }
  # Find credible intervals associated with X% probability, as inputted by the 
  # user. 
  predicted.smooth = fitted_values(object = model$gam, data = predictions, ci_level = probability, scale = "response")
  return(list(predicted.smooth, predictions))
}
# This function calculates manifold eccentricity across each of Yeo's (2011) 7
# resting-state functional connectivity networks. 
calculate.yeo.7.network.manifold.eccentricity = function(df){
  # The first input is contains subject meta-data and columns representing node-
  # wise manifold eccentricity for one modality. Specify the networks!
  yeo_7_networks = c("Vis","SomMot","DorsAttn","SalVentAttn","Limbic","Cont","Default")
  # Extract all meta-data from the input data frame
  metadata = df[,!grepl("7Networks_",names(df))]
  # Initialize an array to hold network-level manifold eccentricity
  network_level_manifold_eccentricity = array(NA,dim=c(nrow(metadata),length(yeo_7_networks)))
  # Subset the data frame to include only regional values, not meta-data
  network_df = df[,which(grepl('LH|RH', colnames(df)))]
  # For each participant, calculate their network-level manifold eccentricity
  for (network_idx in 1:length(yeo_7_networks)){
    network_level_manifold_eccentricity[,network_idx] = 
      rowMeans(network_df[,grepl(yeo_7_networks[network_idx],colnames(network_df))])
  }
  # Add network labels as column names, merge with meta-data, and return
  metadata_with_eccentricity = network_level_manifold_eccentricity %>% 
    as.data.frame %>% `colnames<-`(yeo_7_networks) %>% bind_cols(metadata)
  return(metadata_with_eccentricity)
}
# This function plots the partial age effects on manifold eccentricity, at a 
# global or intrinsic connectivity network (ICN) level. 
plot.manifold.eccentricity.partial.age.effect = function(model,modality,alpha){
  # Extract the outcome variable
  outcome_var = names(model$gam$model[1])
  # Load the data used to fit the model.
  df = model$gam$model
  # Predict manifold eccentricity as an interaction between data set and age,
  # and extract gaussian simultaneous confidence intervals for the predicted 
  # values. The function below does this both...
  pred.with.gaussian.ci = 
    gaussian.simultaneous.confidence.intervals(
      model = model, seed = 100, num_sim = 1000, probability = .95, 
      num_age_breaks = 200, categorical_group = "dataset", control_continuous = "meanfwd")
  # Add a dummy outcome variable
  pred.with.gaussian.ci[[outcome_var]] = 1
  # To compare the normalized age effect on manifold eccentricity between data
  # sets, we need to create a common scale, and thus visualize the interaction
  # between age and cohort minus the main cohort effect. For each data set, 
  # subtract the mean fit from the individual fit values, and assign to the
  # first slice of the second dimension of normalized_effects_array. Normalize 
  # the confidence intervals, and assign to the last two slices of the array.
  normalized_effects_array = array(NA,dim=c(nrow(pred.with.gaussian.ci),3))
  for (dataset_idx in 1:length(datasets)){
    dataset = datasets[dataset_idx]
    # Find the indices of rows in the predicted_manifold_eccentricity data frame
    # corresponding to this data set
    dataset_df_idx = which(pred.with.gaussian.ci$dataset == dataset)
    # For this data set, subtract mean fit from individual fit values
    normalized_effects_array[dataset_df_idx,1] = 
      pred.with.gaussian.ci$fit[dataset_df_idx] - mean(pred.with.gaussian.ci$fit[dataset_df_idx])
    # In the non-normalized space, calculate the difference between the lower 
    # CI and the original fit
    lower.ci.fit.difference = pred.with.gaussian.ci$fit[dataset_df_idx] - 
      pred.with.gaussian.ci$lwrS[dataset_df_idx]
    # Subtract this value from the normalized fit vector to give the normalized 
    # lower CI vector
    normalized_effects_array[dataset_df_idx, 2] = 
      normalized_effects_array[dataset_df_idx,1] - lower.ci.fit.difference
    # Add the value to the normalized fit vector to get the normalized upper CI
    # vector
    normalized_effects_array[dataset_df_idx, 3] =
      normalized_effects_array[dataset_df_idx, 1] + lower.ci.fit.difference
  }
  # Assign to predicted_manifold_eccentricity. Multiply by 1000 for easier 
  # visualisation.
  pred.with.gaussian.ci$normalized_age_effect = normalized_effects_array[,1]
  pred.with.gaussian.ci$normalized_lower_ci = normalized_effects_array[,2]
  pred.with.gaussian.ci$normalized_upper_ci = normalized_effects_array[,3]
  # Specify the breaks to be used for the plot's X axis (age)
  age_breaks = seq(min(pred.with.gaussian.ci$scan_age), max(pred.with.gaussian.ci$scan_age), by=3)
  # Also specify the Y axis breaks and limits (normalised age effect) - use the
  # same scale for all plots.
  if (outcome_var == "global_manifold_eccentricity" & modality == "structural" | modality == "functional"){
    y_axis_breaks = seq(round(min(pred.with.gaussian.ci$normalized_lower_ci),4),
                        round(max(pred.with.gaussian.ci$normalized_upper_ci), 4)+.0001, length.out = 4)
  } else if (outcome_var != "global_manifold_eccentricity" & modality == "structural"){
    y_axis_breaks = seq(-.0044, .0028, length.out = 4)
  } else{
    y_axis_breaks = seq(-.019, .034, length.out = 4)
  }
  # We only include axis and labels for the Y axis for the first ICN, as we're
  # plotting all trends on a common scale. For the global manifold eccentricity
  # plots, include all labels, as these will be on a different scale to the ICNs.
  if (outcome_var == "global_manifold_eccentricity"){
    specified.x.axis.text = element_text(size=40, colour = "black")
    specified.y.axis.text = element_text(size=40, colour = "black")
    specified.x.axis.tick.length = unit(.2,"cm")
    specified.y.axis.tick.length = unit(.2,"cm")
    specified.x.axis.title = element_text(size = 30, face = "bold")
    specified.y.axis.title = element_text(size = 30, face = "bold")
    # Use the same sized text, but colour white. This allows all plots to have
    # the same dimensions.
  } else if (outcome_var == yeo_7_networks[1]){
    specified.x.axis.text = element_text(size = 40, colour = "white")
    specified.y.axis.text = element_text(size = 40, colour = "black")
    specified.x.axis.tick.length = unit(0,"cm")
    specified.y.axis.tick.length = unit(0.2, "cm")
    specified.x.axis.title = element_text(size = 30, colour = "white")
    specified.y.axis.title = element_text(size = 30, face = "bold")
  } else{
    specified.x.axis.text = element_text(size = 40, colour = "white")
    specified.y.axis.text = element_text(size = 40, colour = "white")
    specified.x.axis.tick.length = unit(0, "cm")
    specified.y.axis.tick.length = unit(0, "cm")
    specified.x.axis.title = element_text(size = 30, colour = "white")
    specified.y.axis.title = element_text(size = 30, colour = "white")
  }
  # Plot the age trajectories!
  partial_normalized_age_effect_manifold_eccentricity_plot = 
    ggplot(pred.with.gaussian.ci, mapping = aes(
    x = scan_age, y = normalized_age_effect, group = dataset)) + 
    # Plot the smooth...
    geom_smooth(method="gam", linewidth = 3, aes(colour = dataset)) +
    # And the Gaussian simultaneous confidence intervals...
    geom_ribbon(mapping = aes(ymin = normalized_lower_ci, ymax = normalized_upper_ci,
                              fill = dataset), alpha = .2) +
    labs(x = "Age (Years)", y = "Normalized Age Effect \non Manifold Eccentricity") +
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1),
          legend.position = "none", axis.text.x = specified.x.axis.text,
          axis.text.y = specified.y.axis.text,
          axis.ticks.length.x = specified.x.axis.tick.length,
          axis.ticks.length.y = specified.y.axis.tick.length,
          axis.title.x = specified.x.axis.title,
          axis.title.y = specified.y.axis.title) +
    scale_x_continuous(breaks = age_breaks, labels = age_breaks) +
    # Visualize the Y axis labels to 2 decimal places
    scale_y_continuous(breaks = y_axis_breaks, labels = formatC(round(y_axis_breaks, 4), digits = 4, format = "f")) +
    # Use the cartesian coordinate system to ensure that the lines of best fit
    # do not expand beyond the confidence intervals
    coord_cartesian(ylim = c(min(y_axis_breaks), max(y_axis_breaks))) +
    scale_fill_manual(values = c("#440154", "#1F968B")) +
    scale_colour_manual(values = c("#440154", "#1F968B"))
  partial_normalized_age_effect_manifold_eccentricity_plot
  return(partial_normalized_age_effect_manifold_eccentricity_plot)

}
# This function plots the partial age effect on structure-function coupling, at
# a global or intrinsic connectivity network (ICN) level.
plot.structure.function.coupling.partial.age.effect = function(model, alpha, axis.labels){
  # Extract the outcome variable
  outcome_var = names(model$gam$model[1])
  # Load the data used to fit the model.
  df = model$gam$model
  # Predict coupling as an interaction between data set and age, and extract 
  # gaussian simultaneous confidence intervals for the predicted values. The 
  # function below does this both...
  pred.with.gaussian.ci = 
    gaussian.simultaneous.confidence.intervals(
      model = model, seed = 100, num_sim = 1000, probability = .95, 
      num_age_breaks = 200, control_continuous = "meanfwd", interaction.with.dataset = TRUE)
  # Add a dummy outcome variable
  pred.with.gaussian.ci[[outcome_var]] = 1
  # To compare the normalized age effect on coupling between data sets, we need
  # to create a common scale, and thus visualize the interaction between age and
  # cohort minus the main cohort effect. For each data set, subtract the mean 
  # fit from the individual fit values, and assign to the first slice of the
  # second dimension of normalized_effects_array. Normalize the confidence
  # intervals, and assign to the last two slices of the array.
  normalized_effects_array = array(NA,dim=c(nrow(pred.with.gaussian.ci),3))
  for (dataset_idx in 1:length(datasets)){
    dataset = datasets[dataset_idx]
    # Find the indices of rows in the predicted_manifold_eccentricity data frame
    # corresponding to this data set
    dataset_df_idx = which(pred.with.gaussian.ci$dataset == dataset)
    # For this data set, subtract mean fit from individual fit values
    normalized_effects_array[dataset_df_idx,1] = 
      pred.with.gaussian.ci$fit[dataset_df_idx] - mean(pred.with.gaussian.ci$fit[dataset_df_idx])
    # In the non-normalized space, calculate the difference between the lower 
    # CI and the original fit
    lower.ci.fit.difference = pred.with.gaussian.ci$fit[dataset_df_idx] - 
      pred.with.gaussian.ci$lwrS[dataset_df_idx]
    # Subtract this value from the normalized fit vector to give the normalized 
    # lower CI vector
    normalized_effects_array[dataset_df_idx, 2] = 
      normalized_effects_array[dataset_df_idx,1] - lower.ci.fit.difference
    # Add the value to the normalized fit vector to get the normalized upper CI
    # vector
    normalized_effects_array[dataset_df_idx, 3] =
      normalized_effects_array[dataset_df_idx, 1] + lower.ci.fit.difference
  }
  # Assign to output variable
  pred.with.gaussian.ci$normalized_age_effect = normalized_effects_array[,1]
  pred.with.gaussian.ci$normalized_lower_ci = normalized_effects_array[,2]
  pred.with.gaussian.ci$normalized_upper_ci = normalized_effects_array[,3]
  # Specify the breaks to be used for the plot's X axis (age)
  age_breaks = seq(min(pred.with.gaussian.ci$scan_age), max(pred.with.gaussian.ci$scan_age), by=3)
  # Also specify the Y axis breaks and limits (normalised age effect) - use the
  # same scale for all plots.
  y_axis_breaks = seq(min(pred.with.gaussian.ci$normalized_lower_ci), max(pred.with.gaussian.ci$normalized_upper_ci), length.out = 4)
  # We only include axis and labels for the Y axis for the first ICN, as we're
  # plotting all trends on a common scale. For the global manifold eccentricity
  # plots, include all labels, as these will be on a different scale to the ICNs.
  if (axis.labels == TRUE){
    specified.x.axis.text = element_text(size=40, colour = "black")
    specified.y.axis.text = element_text(size=40, colour = "black")
    specified.x.axis.tick.length = unit(.2,"cm")
    specified.y.axis.tick.length = unit(.2,"cm")
    specified.x.axis.title = element_text(size = 30, face = "bold")
    specified.y.axis.title = element_text(size = 30, face = "bold")
  } else{
    specified.x.axis.text = element_text(size = 40, colour = "white")
    specified.y.axis.text = element_text(size = 40, colour = "white")
    specified.x.axis.tick.length = unit(0, "cm")
    specified.y.axis.tick.length = unit(0, "cm")
    specified.x.axis.title = element_text(size = 30, colour = "white")
    specified.y.axis.title = element_text(size = 30, colour = "white")
  }
  # Plot the age trajectories!
  partial_normalized_age_effect_plot = 
    ggplot(pred.with.gaussian.ci, mapping = aes(
      x = scan_age, y = normalized_age_effect, group = dataset)) + 
    # Plot the smooth...
    geom_smooth(method="gam", linewidth = 3, aes(colour = dataset)) +
    # And the Gaussian simultaneous confidence intervals...
    geom_ribbon(mapping = aes(ymin = normalized_lower_ci, ymax = normalized_upper_ci,
                              fill = dataset), alpha = .2) +
    labs(x = "Age (Years)", y = "Normalized Age Effect on \nStructure-Function Coupling") +
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1),
          legend.position = "none", axis.text.x = specified.x.axis.text,
          axis.text.y = specified.y.axis.text,
          axis.ticks.length.x = specified.x.axis.tick.length,
          axis.ticks.length.y = specified.y.axis.tick.length,
          axis.title.x = specified.x.axis.title,
          axis.title.y = specified.y.axis.title) +
    scale_x_continuous(breaks = age_breaks, labels = age_breaks) +
    # Visualize the Y axis labels to 2 decimal places
    scale_y_continuous(breaks = y_axis_breaks, labels = formatC(round(y_axis_breaks, 4), digits = 4, format = "f")) +
    # Use the cartesian coordinate system to ensure that the lines of best fit
    # do not expand beyond the confidence intervals
    coord_cartesian(ylim = c(min(y_axis_breaks), max(y_axis_breaks))) +
    scale_fill_manual(values = c("#440154", "#1F968B")) +
    scale_colour_manual(values = c("#440154", "#1F968B"))
  partial_normalized_age_effect_plot
  return(partial_normalized_age_effect_plot)
}
# This function adds simultaneous confidence intervals to GAMM-predicted data,
# by sampling from a Gaussian posterior. Adapted from 
# https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/ 
gaussian.simultaneous.confidence.intervals = function(model, seed, num_sim, probability, num_age_breaks, interaction.with.dataset, control_continuous){
  # In brief, for each age sample, we find the maximum absolute standardized 
  # deviation of the fitted function from the true function, and then find the 
  # probability quantiles of these standardized deviations! We sample from the 
  # Gaussian posterior to get an estimate of the likelihood that the true 
  # population mean at any point is held a specific limit (confidence intervals).
  # Extract the Bayesian covariance matrix of the model coefficients. 
  Vb = vcov(model$gam)
  # This function will generate random values from a multivariate normal.
  rmvn <- function(n, mu, sig) { ## MVN random deviates
    L <- mroot(sig)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m*n), m, n))
  }
  # Predict manifold eccentricity using num_age_breaks equally-spaced ages, 
  # spanning from the minimum to the maximum ages in the inputted data frame, 
  # whilst keeping all other variables the same! Do this separately for the 
  # frequency of each category. 
  if (interaction.with.dataset == TRUE){
    gamm_predictions_output = 
      gamm.predictions(model = model, predictors = c("scan_age", "dataset"), 
                       probability = probability, rm.ranef = F, control.continuous = control_continuous)
  } else{
    gamm_predictions_output = 
      gamm.predictions(model = model, predictors = "scan_age", probability = probability, control.continuous = control_continuous, rm.ranef = FALSE)
  }
  # Since we're predicting data based on a continuous and categorical variable,
  # gamm_predictions_output is a list of length 2 - first being the predicted 
  # data frame, and second being the new data on which predictions were made.
  # Use the GAM to predict new values using this new data set. Note that we 
  # already have the predicted values from the first entry of the 
  # gamm_predictions_output list, but we need the standard error for the fit!
  newd = gamm_predictions_output[[2]]
  pred = predict(model$gam, newd, se.fit = TRUE)
  se.fit = pred$se.fit
  # Set the seed!
  set.seed(seed)
  # Take N draws from an approximately distributed multivariate normal with mean
  # of 0 and covariance matrix from the model.
  BUdiff = rmvn(num_sim, mu = rep(0, nrow(Vb)), sig = Vb)
  # Calculate the difference between the estimates and true function
  Cg = predict(model$gam, newd, type="lpmatrix")
  simDev = Cg %*% t(BUdiff)
  # Find the absolute values of the standardised deviations from the true model
  absDev = abs(sweep(simDev, 1, se.fit, FUN = "/"))
  # Find the maximum of the absolute standardised deviations from the true model
  masd = apply(absDev, 2L, max)
  # Find the critical value for the specified probability quantile
  crit = quantile(masd, prob = probability, type = 8)
  predicted = transform(cbind.data.frame(pred, newd),
                        uprP = fit + (2 * se.fit),
                        lwrP = fit - (2 * se.fit),
                        uprS = fit + (crit * se.fit),
                        lwrS = fit - (crit * se.fit))
  return(predicted)
}
# This function tests the significance of a GAMM interaction by comparing it to
# a main effect only model, adapted from code by Bart Larsen:
# (https://bart-larsen.github.io/GAMM-Tutorial/). When setting tensor == TRUE,
# specify the interaction we'd like to test, using the same syntax as the model.
gamm.interaction.significance = function(model, num_sim, tensor, interaction){
  # The first argument is the full GAMM model, with the interaction. The second
  # argument is the number of bootstrap simulations. Extract data used to fit
  # the GAMM.
  df = model$gam$model
  # Set the grouping variable to ID (random effects)
  group_var_name = "id"
  # And the outcome...
  resp = as.character(model$gam$terms[[2]])
  # Get the formula
  full_model = model$gam$formula
  # Get the terms
  theseVars = attr(terms(full_model), "term.labels")
  if (tensor == FALSE){
    # To form the main effects model, remove the interaction term from the formula
    main_effects_model = reformulate(theseVars[0:(length(theseVars)-1)],response = resp)
    # When the model is not a tensor, we have a single interaction term, which 
    # is at the end of the formula.
  } else{
    # When we have a tensor GAMM, the ordering of the interaction terms within 
    # the model may differ, as we tend to have several interactions. Therefore, 
    # the user can specify the interaction they'd like to test. 
    interaction_idx = which(theseVars %in% interaction)
    main_effects_model = reformulate(setdiff(theseVars, theseVars[interaction_idx]), response = resp)
  }
  # Fit a GAM to each model
  g1 = gam(full_model, data = df)
  g2 = gam(main_effects_model, data = df)
  # Get the design matrices
  mat1 = model.matrix(g1)
  mat2 = model.matrix(g2)
  # Add on the response variable and grouping variable
  group_var = df[,group_var_name]
  y = df[,resp]
  # For the original interaction model and main effects model, conduct a GAM, 
  # and extract the model matrix for this. The model design matrices contain all 
  # possible combinations of predictors, including interactions where 
  # applicable, with categorical predictors dummy-coded as constants. Use this
  # alongside the random effects variable to predict the outcome as a linear 
  # mixed-effects model.
  m1 = lmer(formula = y ~ -1 + mat1 + (1|group_var), data = df)
  m2 = lmer(formula = y ~ -1 + mat2 + (1|group_var), data = df)
  # To speed up the simulation, do computations with multiple processors.
  (nc <- detectCores())
  cl = makeCluster(rep("localhost", nc))
  # Create a bootstrap distribution with user-inputted number of simulations.
  ref_dist = PBrefdist(m1, m2, nsim = num_sim, seed = 100, cl=cl)
  # Compare the models, and extract the p-value for the bootstrap
  pb = PBmodcomp(m1,m2,nsim=numsim,ref = ref_dist, seed = 100, cl=cl)
  stopCluster(cl)
  pb_stats = pb$test["PBtest","p.value"]
  print(sprintf("Parametric bootstrapped likelihood ratios has a p-value of %.3f.",
                 round(pb_stats,3)))
  return(round(pb_stats, digits = 3))
}
# This function identifies significant periods of change in GAMMs, using 
# confidence intervals for the first derivative of the GAM part of GAMM. The 
# output is a vector of ages across which the change in the independent variable
# was significant. We do this separately for each data set. Adapted from
# https://fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-
# with-gams/, and https://bart-larsen.github.io/ GAMM-Tutorial/.
significant.developmental.change = function(model, probability, parametric, tensor){
  # Get the names of the data sets
  datasets = levels(model$model$dataset)
  # Get the name of the smooth, and find whether there's a tensor or not!
  if (tensor == FALSE){
    term = "s(scan_age)"
  } else{
    term = "ti(scan_age)"
  }
  # If we're examining confidence intervals across data sets...
  if (is.null(datasets) == FALSE){
    # Create an output list for the significant derivatives for each data set
    significant_deriv_list = vector("list", length(datasets))
    # Derive derivatives for each data set separately
    for (dataset in datasets){
      # Extract the data for this data set on which the model was fitted
      subset_df = model$model[which(model$model$dataset == dataset),]
      # Keep the columns which are the model variables (parametric and smooth terms)
      subset_df = subset_df[,c(parametric, "scan_age")]
      # Calculate the derivatives for each age. Confidence intervals are provided by
      # default. 
      model.d = gratia::derivatives(model, level = probability, term = term, data = subset_df)
      # Order by ascending age
      model.d = model.d[order(model.d$data), ]
      # Search where the confidence intervals do not contain 0
      significant_change_idx = c(which(model.d$lower < 0 & model.d$upper > 0))
      if (length(significant_change_idx) == 0){
        print(sprintf('No significant periods of change in %s.',dataset))
        # Append NA to the output list
        significant_deriv_list[[match(dataset, datasets)]] = NA
      } else{
        print(sprintf('Significant changes between %.2f and %.2f years old in %s, at %.1f intervals.', 
                      model.d$data[significant_change_idx[[1]]], 
                      model.d$data[[significant_change_idx[[length(significant_change_idx)]]]], dataset,
                      length(significant_change_idx)))
        # Set the non-significant derivatives (with confidence intervals crossing 
        # zero) to NA
        model.d$derivative[setdiff(1:nrow(model.d), significant_change_idx)] = NA
        # Append to output list
        significant_deriv_list[[match(dataset, datasets)]] = model.d
      }
    }
  } else{
    if (length(significant_change_idx) == 0){
      print('No significant periods of change in across datasets.')
    } else{
      print(sprintf(
        'Significant changes between %.2f and %.2f years old, at %.1f intervals.', 
        model.d$data[significant_change_idx[[1]]], model.d$data[[significant_change_idx[[length(significant_change_idx)]]]],
        length(significant_change_idx)))
      # Calculate the derivatives for each age, across data sets.
      model.d = gratia::derivatives(model, level = probability, term = "ti(scan_age)", n = nrow(gamm.model$gam$model))
      # Search where the confidence intervals do not contain 0
      significant_change_idx = c(which(model.d$lower < 0 & model.d$upper > 0))
      if (length(significant_change_idx) == 0){
        print('No significant change in rate of development across data sets.')
      } else{
        print(sprintf('Significant changes between %.2f and %.2f years old, at %.1f intervals.', 
                      model.d$data[significant_change_idx[[1]]], 
                      model.d$data[[significant_change_idx[[length(significant_change_idx)]]]],
                      length(significant_change_idx)))
      }
      # Set the non-significant derivatives (with confidence intervals crossing 
      # zero) to NA
      model.d$derivative[setdiff(1:nrow(model.d), significant_change_idx)] = NA
      return(model.d)
    }
  }
  return(significant_deriv_list)
}
# This function plots the age range across which age effects on manifold 
# eccentricity are significant in the GAMM i.e. having non-zero confidence 
# intervals for the age derivatives. We do this separately for each data set. 
# This code is adapted from Valerie Sydnor's 2023 Nature Neuroscience paper
# 'Intrinsic activity development unfolds along a sensorimotor-association 
# cortical axis in youth'. 
plot.sensitive.age.developmental.period = function(model, smooth, parametric, tensor, by_dataset){
  # Find where the age derivatives have non-zero confidence intervals. Note that
  # 'model' must be the GAM model only, not a list of GAMM outputs. The output 
  # of this will be two data frames - one for each data frame.
  significant_age_derivatives = significant.developmental.change(
    model = model, probability = .95, parametric = parametric, tensor = tensor)
  if (by_dataset == TRUE){
    # Plot each set of significant derivatives for each data set, using different
    # colors.
    dataset_cols = c("#440154", "#1F968B")
    # Initialize list to hold the plots for each data set
    dataset_plot_list = vector("list", length = 2)
    for (dataset_idx in 1:length(dataset_cols)){
      # First, check that this data set has periods of significant change in the 
      # first derivative of age
      if (length(significant_age_derivatives[[dataset_idx]]) == 1){
        # If there's no significant change, assign NA to dataset_plot_list
        dataset_plot_list[[dataset_idx]] = NA
      } else{
        # Select the data frame for this data set
        age_deriv_df = significant_age_derivatives[[dataset_idx]]
        # Find which derivatives are significant (non-zero confidence intervals)
        significant_age_deriv_df = age_deriv_df[which(!is.na(age_deriv_df$derivative)),] 
        # Sort by increasing age (labelled 'data' in data frame)
        significant_age_deriv_df = significant_age_deriv_df[order(significant_age_deriv_df$data), ]
        # Find the smallest and largest significant derivative, and the ages at 
        # which they occur. 
        earliest_age_deriv = significant_age_deriv_df$data[1]
        latest_age_deriv = significant_age_deriv_df$data[nrow(significant_age_deriv_df)]
        if (nrow(significant_age_deriv_df) != 0){
          # Now plot!
          derivative_plot = 
            ggplot(data = age_deriv_df, mapping = aes(x = data, y = .1)) +
            geom_rect(mapping = aes(xmin = earliest_age_deriv, xmax = latest_age_deriv,
                                    ymin = -Inf, ymax = Inf), 
                      colour = alpha(dataset_cols[dataset_idx], .5),
                      fill = alpha(dataset_cols[dataset_idx], .5)) +
            theme_classic() +
            theme(legend.position = "none", axis.title = element_blank(),
                  axis.text = element_blank(), axis.line = element_blank(),
                  axis.ticks = element_blank())
          dataset_plot_list[[dataset_idx]] = derivative_plot
        }
      }
    }
    return(dataset_plot_list)
  } else{
    # Check that there are significant age derivatives!
    if (sum(is.na(significant_age_derivatives$derivative)) == nrow(significant_age_derivatives)){
      return(NA)
    } else{
      # Find which derivatives are significant (non-zero confidence intervals)
      significant_age_deriv_df = significant_age_derivatives[which(!is.na(significant_age_derivatives$derivative)),] 
      # Find the smallest and largest significant derivative, and the ages at 
      # which they occur. 
      earliest_age_deriv = significant_age_deriv_df$data[1]
      latest_age_deriv = significant_age_deriv_df$data[nrow(significant_age_deriv_df)]
      # Now plot!
      derivative_plot = 
        ggplot(data = significant_age_deriv_df, mapping = aes(x = data, y = .1)) +
        geom_rect(mapping = aes(xmin = earliest_age_deriv, xmax = latest_age_deriv,
                                ymin = -Inf, ymax = Inf), colour = "blue", fill = "blue") +
        theme_classic() +
        theme(legend.position = "none", axis.title = element_blank(),
              axis.text = element_blank(), axis.line = element_blank(),
              axis.ticks = element_blank())
      return(derivative_plot)
     }
  }
}

