### PART 1 - Set Up the Work Space ####
rm(list = ls())
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(stringr)
library(colorRamp2)
library(circlize)
library(cowplot)
# Add the GAMM functions and set the working directory
setwd('U:/gradients_open_access/')
source('code/GAMM.functions.v3.R')
# Import the loadings of each cognitive/psychopathological measure onto each 
# component derived from principal components analysis.
psychopathology.participant.loadings = read.csv('data/phenotypic/psychopathology.participant.loadings.csv')
cognition.participant.loadings = read.csv('data/phenotypic/cognitive.participant.loadings.csv')
# And import the equivalent for individual measures
psychopathology.measure.loadings = read.csv('data/phenotypic/psychopathology.measure.loadings.csv')
cognition.measure.loadings = read.csv('data/phenotypic/cognitive.measure.loadings.csv')
# Load the structure-function coupling measure
structure.function.coupling.df = read.csv('data/structure.function/coupling.data.df.csv')[,-1]
# Rename the global.coupling column
names(structure.function.coupling.df)[names(structure.function.coupling.df) == 'global.coupling'] <- 'Global'

### PART 2 - Visualize Distribution of Participant Loadings for Psychopathology and Cognition ####
# Concatenate the two participant loading data frames and create a new variable
# which indexes where each data frame came from
cognition.psychopathology.participant.loadings = 
  bind_rows(psychopathology.participant.loadings, cognition.participant.loadings, .id = "domain") %>%
  mutate(domain = factor(domain, labels = c("psychopathology", "cognition"))) %>%
  mutate(dataset = factor(dataset)) %>%
  select(!age_in_months) 

# Loop across each factor, and plot the distribution of participant loadings as
# a function of domain and data set.
domains = c("cognition", "psychopathology")
for (domain_name in domains){
  if (domain_name == "cognition"){ 
    factors = c("Factor1", "Factor2")
  } else{
    factors = c("Factor1", "Factor2", "Factor3")
  }
  for (factor_name in factors){
    # Set the x-axis limits...
    domain.participant.loadings = filter(cognition.psychopathology.participant.loadings, domain == domain_name)
    xaxis_breaks = seq(from = min(domain.participant.loadings[[factor_name]]),
                       to = max(domain.participant.loadings[[factor_name]]),
                       length.out = 5)
    # And plot!
    participant.loadings.plot = 
      cognition.psychopathology.participant.loadings %>% filter(domain == domain_name) %>%
      ggplot(., mapping = aes(x= .data[[factor_name]], fill=dataset)) +
      geom_boxplot(mapping = aes(colour = dataset), lwd = 1) +
      labs(x = paste0("PC", str_sub(factor_name, start = -1), " Loading")) +
      theme(panel.background = element_blank(), legend.position = "none",
            axis.text.x = element_text(size=20), axis.title.x = element_text(size = 20, face="bold"),
            axis.ticks.length.y = unit(0, "cm"), axis.line.x = element_line(colour = "black"), 
            axis.text.y = element_blank(), axis.title.y = element_blank(),
            axis.ticks.length.x = unit(0.2, "cm")) +
      scale_fill_manual(values = c("calm"="white", "nki" = "white")) +
      scale_colour_manual(values = c("calm" = "#440154FF", "nki" = "#1F968BFF")) +
      scale_x_continuous(breaks = xaxis_breaks, labels = formatC(xaxis_breaks, digits = 2, format = "f")) +
      coord_cartesian(xlim = c(min(xaxis_breaks), max(xaxis_breaks)))
    ggsave(filename = sprintf("phenotypic/%s.component.%s.individual.loadings.png",domain_name,str_sub(factor_name, start = -1)),
           participant.loadings.plot, width = 5, height = 5, dpi = 700)
  }
}

### PART 3 - Visualize Measure Loadings onto Each Component ####
# Set the column names of each data frame as describing the measures used
names(psychopathology.measure.loadings)[names(psychopathology.measure.loadings) == "X"] <- "Measure"
psychopathology.measure.loadings$Measure = factor(psychopathology.measure.loadings$Measure)
names(cognition.measure.loadings)[names(cognition.measure.loadings) == "X"] <- "Measure"
cognition.measure.loadings$Measure = factor(cognition.measure.loadings$Measure)
# Loop across domains and components... Plot the distribution of the loadings of
# each measure onto each component...
for (domain_name in domains){
  domain.df = get(paste0(domain_name,".measure.loadings"))
  factors = colnames(domain.df)[-1]
  for (factor_name in factors){
    # Specify the y-axis limits
    if (all(domain.df[[factor_name]] > 0)){
      yaxis_breaks = seq(from = 0, to = max(domain.df[[factor_name]]), length.out=5)
    } else{
      yaxis_breaks = seq(from = min(domain.df[[factor_name]]), to = max(domain.df[[factor_name]]), length.out=5)
    }
    # If 0 is included in the y axis-breaks, specify this separately i.e. not to
    # 2 decimal places, but as an integer
    if (0 %in% yaxis_breaks){
      formatted_y_axis_labels = c(0, formatC(yaxis_breaks, digits = 2, format = "f")[-1])
    } else{
      formatted_y_axis_labels = formatC(yaxis_breaks, digits = 2, format = "f")
    }
    measure.loading.plot = ggplot(domain.df, mapping = aes(x = Measure, y = .data[[factor_name]])) +
      geom_bar(stat = "identity", fill = "white", lwd = 1, colour = "#2D708EFF") +
      ylab(paste0("PC ", str_sub(factor_name, start=-1))) +
      theme(panel.background = element_blank(), axis.ticks.length.y = unit(0, "cm"), 
            axis.text.y = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_text(size=20), axis.line.x = element_line(colour="black"), 
            axis.title.x = element_text(size = 20), axis.ticks.length.x = unit(0.2, "cm"), legend.position = "none") +
      scale_y_continuous(breaks = yaxis_breaks, labels = formatted_y_axis_labels) +
      coord_flip()
    measure.loading.plot
    ggsave(sprintf('phenotypic/%s.component.%s.measure.loadings.plot.png',domain_name,str_sub(factor_name, start=-1)),
           measure.loading.plot, height = 5, width = 5, dpi = 700)
  }
}

### PART 4 - Linking Structure-Function Coupling with Dimensions of Psychopathology and Cognition ####
# Specify the domains we'll examine and the networks
domains = c("psychopathology", "cognition")
outcome_array = c("Global", "Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default")
cognition_coupling_gamm_output_list = vector("list",length(outcome_array))
psychopathology_coupling_gamm_output_list = vector("list",length(outcome_array))
for (domain in domains){
  domain.df = get(paste0(domain,".participant.loadings"))
  # Merge with the structure-function coupling data frame!
  domain.structure.function.df = merge(structure.function.coupling.df, domain.df, by = c("dataset", "id", "timepoint"), all.x = TRUE)
  # Find the number of factors
  factors = colnames(domain.df)[which(grepl("Factor", colnames(domain.df)))]
  nfactors = length(factors)
  # Fit a GAMM with subject-specific random intercepts and the following 
  # predictors for global and network-level structure-function coupling: 
  # parametric terms of sex and mean frame-wise displacement. Additionally, we
  # include interactions between scan age and the domain measure (principal 
  # component loading for each participant). These are specified as tensors. 
  # Any main effects within the interaction are also specified as tensors. 
  parametric = c("sex", "meanfwd")
  # Specify the tensor interaction
  tensor.interaction = vector("list", length = nfactors)
  for (factor_idx in 1:nfactors){
    factor = factors[factor_idx]
    tensor.interaction[[factor_idx]] = c("scan_age", factor)
  }
  # Construct this GAMM with the outcome being coupling in one of the 8 levels 
  # of analysis in outcome_array. We use all factors as main effects and 
  # interactions with age as predictors, alongside co-variates of mean framewise
  # displacement, sex, and age. Initialize a list to hold the GAMM outputs and
  # an array to hold the p-values.
  coupling_summary_tables = vector("list", length = length(outcome_array))
  npredictors = length(parametric) + length(unlist(tensor.interaction)) + 1
  pval_array = array(numeric(), dim = c(npredictors, length(outcome_array), 2))
  # Loop across each outcome measure, with the outcome being coupling in one of
  # the 8 levels of analysis in outcome_array.
  for (outcome_idx in 1:length(outcome_array)){
    outcome_measure = outcome_array[outcome_idx]
    print(sprintf('%s structure-function: Examining link with dimensions of %s.',
            outcome_measure, domain))
    # Construct the GAMM
    coupling.gamm = fit.gamm.tensor.interaction.with.random.effects(
      df = domain.structure.function.df, tensor.interaction = tensor.interaction,
      outcome = outcome_measure, knots = 4, fx = FALSE, parametric = parametric)
    # Assign un-corrected p-values for each predictor to output array. 
    pval_array[,outcome_idx,1] = coupling.gamm[[1]][,3]
    # Assign the results table summary to the output
    coupling_summary_tables[[outcome_idx]] = coupling.gamm[[1]]
    # And assign the models themselves to the output list
    if (domain == "psychopathology"){
      psychopathology_coupling_gamm_output_list[[outcome_idx]] = coupling.gamm[[2]]
    } else{
      cognition_coupling_gamm_output_list[[outcome_idx]] = coupling.gamm[[2]]
    }
  }
  # Correct for multiple comparisons across networks and globally by controlling
  # the false discovery rate.
  pval_array[,,2] = p.adjust(pval_array[,,1], method="fdr")
  # For each outcome, replace the un-corrected p-values in the summary table 
  # with the FDR-corrected version.
  for (outcome_idx in 1:length(outcome_array)){
    summary_table = coupling_summary_tables[[outcome_idx]]
    summary_table$X3 = pval_array[,outcome_idx,2]
    # Format the table nicely
    colnames(summary_table) = c("Estimate/EDF", "t/F", "FDR-corrected p-value")
    write.csv(summary_table, file = sprintf("data/coupling.phenotypic/%s.%s.coupling.gamm.csv", outcome_array[outcome_idx], domain))
  }
}
saveRDS(cognition_coupling_gamm_output_list, file = "data/coupling.phenotypic/cognition.coupling.gamm.output.list.RData")
saveRDS(psychopathology_coupling_gamm_output_list, file = "data/coupling.phenotypic/psychopathology.coupling.gamm.output.list.RData")

### PART 5 - Visualizing Age * Factor Interaction when Predicting Cognition ####
# Our results show a significant interaction between age and cognitive factor 2
# when predicting coupling within the default-mode, dorsal attention, global, 
# and somato-motor networks. Therefore, for these networks, we shall split the
# participants according to cognitive factor, conduct a GAMM within each sub-
# group, and visualise age effects. 
cognition_coupling_gamm_output_list = readRDS("data/coupling.phenotypic/cognition.coupling.gamm.output.list.RData")
idx_to_plot = c(1, 3, 4, 8)
decomposed_interaction_effect_list = vector("list", length(idx_to_plot))
for (idx in idx_to_plot){
  interaction.plot = plot.continuous.tensor.interaction.with.age(
    gamm_model = cognition_coupling_gamm_output_list[[idx]],
    grouping_var = "Factor2", decompose_developmental_effects = TRUE)
  # Assign to output list and save individual plot...
  decomposed_interaction_effect_list[[idx]] = interaction.plot
  ggsave(filename = sprintf('data/coupling.phenotypic/decomposed.age.interaction.cognition.%s.png', outcome_array[idx]),
         plot = interaction.plot, height = 10, width = 10, dpi = 700)
  rm(interaction.plot)
}
