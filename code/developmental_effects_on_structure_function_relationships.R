### PART 1 - Set up the Work Space ####
rm(list = ls())
library(raveio)
library(dplyr)
library(ggplot2)
library(stringr)
library(rhdf5)
library(reshape2)
library(cowplot)
setwd('//cbsu/data/Imaging/astle/am10/diffusion.map.embedding/')
# Add the custom GAMM function code.
reticulate::py_run_string("import sys")
source('C:/Users/am10/diffusion_map_embedding_code/GAMM.functions.v3.R')
# Load the structure-function coupling metrics for CALM, across both time points.
calm.baseline.coupling = read.csv('data/calm/structure.function/euclidean.distance.spearman.baseline.schaefer200x7.csv')[,-1]
calm.followup.coupling = read.csv('data/calm/structure.function/euclidean.distance.spearman.followup.schaefer200x7.csv')[,-1]
# Load the coupling metrics for NKI, across the three time points
nki.bas1.coupling = read.csv('data/nki/structure.function/euclidean.distance.spearman.bas1.schaefer200x7.csv')[,-1]
nki.flu1.coupling = read.csv('data/nki/structure.function/euclidean.distance.spearman.flu1.schaefer200x7.csv')[,-1]
nki.flu2.coupling = read.csv('data/nki/structure.function/euclidean.distance.spearman.flu2.schaefer200x7.csv')[,-1]
# Load the labels for Yeo's (2011) 7 intrinsic functional connectivity networks.
yeo_7_networks = c("Vis","SomMot","DorsAttn","SalVentAttn","Limbic","Cont","Default")
# And set the data sets available
datasets = c("calm", "nki")
# Load the meta-data for each data set. 
calm.metadata.sc = read.csv('data/calm/dme/dme.and.metadata.structural.connectivity.csv')[,-1]
calm.metadata.fc = read.csv('data/calm/dme/dme.and.metadata.functional.connectivity.csv')[,-1]
nki.metadata.sc = read.csv('data/nki/dme/dme.and.metadata.structural.connectivity.csv')[,-1]
nki.metadata.fc = read.csv('data/nki/dme/dme.and.metadata.functional.connectivity.csv')[,-1]
### PART 2 - Calculate mean frame-wise displacement across structural and functional scans ####
# We'll use this as the estimate of head motion when assessing age effects on
# structure-function coupling. Start with CALM...
for (dataset in datasets){
  if (dataset == "calm"){
    timepoints = c("baseline", "followup")
    timepoint_var = "timepoint"
  } else{
    timepoints = c("bas1", "flu1", "flu2")
    timepoint_var = "session"
  }
  for (timepoint_name in timepoints){
    # Get the coupling data frame
    coupling.df = get(paste0(dataset, ".", timepoint_name, ".coupling"))
    # Get the SC meta-data sheet and subset by the time-point 
    sc = get(paste0(dataset, ".metadata.sc"))
    subset_sc = sc[which(sc[[timepoint_var]] == timepoint_name),]
    # Get the mean FWD for SC
    meanfwd_sc = subset_sc[which(subset_sc$id %in% coupling.df$id),"meanfwd"]
    # Repeat for FC...
    fc= get(paste0(dataset, ".metadata.fc"))
    subset_fc = fc[which(fc[[timepoint_var]] == timepoint_name),]
    meanfwd_fc = subset_fc[which(subset_fc$id %in% coupling.df$id),"meanfwd"]
    # Calculate the mean FWD and assign to appropriate data frame
    meanfwd_average = (meanfwd_sc + meanfwd_fc)/2
    if (dataset == "calm" & timepoint_name == "baseline"){
      calm.baseline.coupling$meanfwd_average = meanfwd_average
    } else if (dataset == "calm" & timepoint_name == "followup"){
      calm.followup.coupling$meanfwd_average = meanfwd_average
    } else if (dataset == "nki" & timepoint_name == "bas1"){
      nki.bas1.coupling$meanfwd_average = meanfwd_average
    } else if (dataset == "nki" & timepoint_name == "flu1"){
      nki.flu1.coupling$meanfwd_average = meanfwd_average
    } else{
      nki.flu2.coupling$meanfwd_average = meanfwd_average
    }
  }
}
# Create a single CALM data frame for structure-function coupling
calm.coupling = bind_rows(calm.baseline.coupling, calm.followup.coupling, .id = "Timepoint") %>%
  mutate(Timepoint = factor(Timepoint, labels = c("baseline", "followup"))) %>%
  # Convert ID to a character column
  mutate(id = as.character(id))
# Repeat for NKI, and merge with CALM
coupling.df = bind_rows(nki.bas1.coupling, nki.flu1.coupling, nki.flu2.coupling, .id = "Timepoint") %>%
  mutate(Timepoint = factor(Timepoint, labels = c("bas1", "flu1", "flu2"))) %>%
  # Convert ID to a character column
  mutate(id = as.character(id)) %>%
  bind_rows(., calm.coupling, .id = "Dataset") %>%
  mutate(Dataset = factor(Dataset, labels = c("nki", "calm")))
# Assign each column to an ICN, and average the coupling within each ICN. Start
# by initializing an array to hold the ICN averages for each participant.
ICN_coupling_array = array(NA, dim = c(nrow(coupling.df),length(yeo_7_networks)))
for (yeo_idx in 1:length(yeo_7_networks)){
  ICN_coupling_array[,yeo_idx] = rowMeans(coupling.df[,grepl(yeo_7_networks[yeo_idx], colnames(coupling.df))])
}
ICN.coupling.df = data.frame(ICN_coupling_array)
colnames(ICN.coupling.df)[1:length(yeo_7_networks)] = yeo_7_networks
ICN.coupling.df$Timepoint = coupling.df$Timepoint
ICN.coupling.df$id = coupling.df$id
ICN.coupling.df$Dataset = coupling.df$Dataset
# Melt into a single long data frame
ICN.coupling.df.melted = melt(ICN.coupling.df) %>%
  rename(., all_of(c(ICN = "variable", coupling = "value")))

### PART 3 - Structure-Function Coupling Generalized Additive Mixed Models (GAMMs) ####
# Subset the coupling data frame by CALM, and merge with the meta-data
calm.coupling.df = subset(ICN.coupling.df, Dataset == "calm") %>%
  # Rename columns to allow for easier merging with meta-data
  rename(., all_of(c(timepoint = "Timepoint"))) %>%
  merge(., calm.metadata.sc, by = c("id", "timepoint"), all.x = TRUE) %>%
  # Subset by the columns we want
  select(all_of(c("id", "sex", "meanfwd", "scan_age", "timepoint", yeo_7_networks)))
# Repeat with NKI!
nki.coupling.df = subset(ICN.coupling.df, Dataset == "nki") %>%
  rename(., all_of(c(session = "Timepoint"))) %>%
  merge(., nki.metadata.sc, by = c("id", "session")) %>%
  select(all_of(c("id", "session", "sex", "meanfwd", "scan_age", yeo_7_networks))) %>%
  # Rename session to timepoint for easier merging later on
  rename(., all_of(c(timepoint = "session")))
# Combine both data frames into one
coupling.df.with.metadata = bind_rows(calm.coupling.df, nki.coupling.df, .id = "dataset") %>%
  mutate(dataset = factor(dataset, labels = c("calm", "nki"))) 
# Calculate global coupling i.e. average coupling across all ICNs
coupling.df.with.metadata$global.coupling = rowMeans(coupling.df.with.metadata[,yeo_7_networks])
coupling.df.with.metadata$global = rowMeans(coupling.df.with.metadata[,yeo_7_networks])
# Ensure that all categorical variables are factors
coupling.df.with.metadata$sex = factor(coupling.df.with.metadata$sex)
coupling.df.with.metadata$id = factor(coupling.df.with.metadata$id)
# And save!
write.csv(x = coupling.df.with.metadata, file = "data/structure.function/coupling.data.df.csv")
# Report the median and interquartile range structure-function coupling for each
# level of analysis, and additionally as a function of dataset.
level_of_analysis = c("global", yeo_7_networks)
for (analysis_level in level_of_analysis){
  print(sprintf('Across datasets, %s coupling had a median of %.2f, and IQR of %.2f.',
          analysis_level, median(coupling.df.with.metadata[[analysis_level]]),
          IQR(coupling.df.with.metadata[[analysis_level]])))
  # Now repeat, but stratified by dataset!
  for (dataset in datasets){
    subset_df = coupling.df.with.metadata[which(coupling.df.with.metadata$dataset == dataset),]
    print(sprintf('In %s, %s coupling had a median of %.2f and IQR of %.2f.',
            dataset, analysis_level, median(subset_df[[analysis_level]]),
            IQR(subset_df[[analysis_level]])))
  }
  
}

### PART 3A - Summary statistics! ####
for (dataset in datasets){
  subset_df = coupling.df.with.metadata[which(coupling.df.with.metadata$dataset == dataset), ]
  if (dataset == "calm"){
    timepoints = c("baseline", "followup")
  } else{
    timepoints = c("bas1", "flu1", "flu2")
  }
  for (timepoint_idx in 1:length(timepoints)){
    timepoint.subset.df = subset_df %>% group_by(id) %>% filter(n() == timepoint_idx)
    df.summary = timepoint.subset.df %>% ungroup %>% summarise(mean = mean(scan_age), sd = sd(scan_age))
    # Get the percentage of male participants
    percent.male = (length(which(timepoint.subset.df$sex == "M")) / nrow(timepoint.subset.df)) *100
    print(sprintf('For %.1f %s participants with %.1f structure-function coupling value(s): 
                  mean age of %.2f, sd of %.2f, and %.2f percent male.', nrow(timepoint.subset.df), dataset,
                  timepoint_idx, round(df.summary$mean, 2), round(df.summary$sd, 2), round(percent.male, 2)))
  }
}

### PART 3B - Running the GAMMs ####
# Now run a GAMM to predict network-level manifold eccentricity, using the 
# same co-variates as above. Initialize an array to hold the corrected (for 
# multiple comparisons) and un-corrected p-values. 
coupling_summary_tables = vector("list", length = length(yeo_7_networks) + 1)
# Also initialize a list to hold the model outputs!
coupling_gamm_output_list = vector("list",length(yeo_7_networks) + 1)
# Initialize the GAMM predictors 
smooth = "scan_age"
parametric = c("dataset", "meanfwd", "sex")
interaction = c("scan_age", "dataset")
all_predictors = c(smooth, parametric, interaction)
# We correct for multiple comparisons across the 7 ICNs and global SF-coupling.
# Therefore, initialize an array to hold the un-corrected and corrected 
# p-values. We subtract 1 from the length of all_predictors as we'll be testing
# the significance of the interaction term separately. We add one onto the 
# length of yeo_7_networks as we'll place global SF-coupling statistics in this
# array too.
pval_array = array(numeric(), dim = c(length(all_predictors)-2, length(yeo_7_networks)+1, 2))
# Run the GAMMs at a global level and for each ICN...
for (idx in 1:(length(yeo_7_networks) + 1)){
  if (idx == 1){
    coupling.gamm = fit.gamm.with.random.effects(
      df = coupling.df.with.metadata, parametric = parametric, smooth = smooth,
      interaction = interaction, knots = 4, fx = FALSE, outcome = "global.coupling",
      report_stats = FALSE)
  } else{
    coupling.gamm = fit.gamm.with.random.effects(
      df=coupling.df.with.metadata, parametric=parametric, smooth=smooth, 
      interaction=interaction, outcome=yeo_7_networks[idx-1], 
      fx=FALSE, knots=4, report_stats = FALSE)
  }
  # Assign un-corrected p-values for each predictor to output array. This holds
  # the p-values for single predictors, not interactions.
  pval_array[,idx,1] = coupling.gamm[[1]][1:length(c(parametric, smooth)),3]
  # Assign the results table summary to the output
  coupling_summary_tables[[idx]] = coupling.gamm[[1]]
  # And assign the models themselves to the output list
  coupling_gamm_output_list[[idx]] = coupling.gamm[[2]]
}
### PART 3C - Correcting for multiple comparisons ####
# Correct for multiple comparisons for global and network-level SF-coupling 
# GAMMs using false-discovery rate. We subtract 2 from the length of 
# all_predictors because we will be testing the significance of the 
# interaction separately. 
for (predictor_idx in 1:length(all_predictors)-2){
  pval_array[predictor_idx,,2] = p.adjust(pval_array[predictor_idx,,1], method="fdr")
}
# Initialize an array to hold the bootstrapped FDR p-values for testing the 
# interaction between age and data set.
interaction_pval = array(NA, dim = length(yeo_7_networks) + 1)
# Calculate the p-value for the interaction between data set and age.
for (idx in 1:(length(yeo_7_networks) + 1)){
  interaction_output = gamm.interaction.significance(
    model = coupling_gamm_output_list[[idx]], num_sim = 10000, tensor = FALSE)
  # Keep the p-value...
  interaction_pval[[idx]] = interaction_output
}
# Correct for multiple comparisons on the interaction p-value by controlling 
# the false discovery rate.
adjusted_interaction_pval = p.adjust(interaction_pval, method = "fdr")
### PART 3D - Summary table of GAMM results ####
# For each network and globally, create a summary table of results. Ensure to 
# include the FDR-corrected p-values for the main effects and FDR-corrected 
# bootstrapped p-values for the interaction term. 
for (idx in 1:(length(yeo_7_networks) + 1)){
  # Select the correct network summary table
  summary_table = coupling_summary_tables[[idx]]
  # The parametric predictors are reported first. Subset by these, and format.
  parametric_coefs = summary_table[1:length(parametric),]
  # Find the indices of parametric coefficients in all_predictors.
  parametric_idx = which(all_predictors %in% parametric)[c(1:length(parametric))]
  # Replace the final column with the FDR-corrected p-values
  parametric_coefs$X3 = pval_array[parametric_idx,idx,2]
  colnames(parametric_coefs) = c("estimate", "t", "FDR-corrected p")
  if (idx == 1){
    print('GAMM with global SF-coupling as an outcome...')
  } else{
    print(sprintf('GAMM with %s network SF-coupling as an outcome...', yeo_7_networks[idx -1]))
  }
  print(parametric_coefs)
  # Now find the smooths!
  smooth_coefs = summary_table[(length(parametric) + 1): nrow(summary_table), ]
  # Replace the interaction p-value from the GAMM with the boostrapped version
  smooth_coefs[nrow(smooth_coefs), ncol(smooth_coefs)] = adjusted_interaction_pval[idx]
  # Replace the smooth p-value from the GAMM with the FDR-corrected p-value
  smooth_coefs[nrow(smooth_coefs)-1, ncol(smooth_coefs)] = pval_array[length(all_predictors)-2,idx,2]
  colnames(smooth_coefs) = c("edf", "F", "FDR-corrected p")
  print(smooth_coefs)
  # Bind rows of the parametric and smooth coefficients, and then replace the 
  # previous summary table. To do so, we'll temporarily use the same column
  # names across both parametric and smooth predictors, and then rename the 
  # columns appropriately.
  colnames(smooth_coefs) = colnames(parametric_coefs)
  updated.summary = bind_rows(parametric_coefs, smooth_coefs)
  colnames(updated.summary) = c("estimate/EDF", "t/F", "FDR-corrected p-value")
  coupling_summary_tables[[idx]] = updated.summary
  # Save this updated summary!
  if (idx == 1){
    save_filename = 'data/structure.function/summary/global.coupling.gamm.csv'
  } else{
    save_filename = sprintf('data/structure.function/summary/%s.coupling.gamm.csv', yeo_7_networks[idx-1])
  }
  write.csv(updated.summary, save_filename)
}
# Format the corrected p-values as a data frame, and save to CSV
pval_array[,,2] %>% as.data.frame() %>% `colnames<-`(c("global",yeo_7_networks)) %>%
  # Add the FDR-corrected bootstrapped p-value for the interaction
  rbind(., adjusted_interaction_pval) %>%
  `rownames<-`(c(parametric,smooth,paste0(interaction[1],":",interaction[2]))) %>%
  write.csv(., 'data/structure.function/coupling.age.gamm.csv')
# Loop through each network again...
### PART 4 - Visualizing trajectories of coupling in the default mode network ####
# Since coupling in the default mode network had the strongest effect of an 
# interaction between data set and age, we'll plot this interaction...
default.mode.interaction.plot = 
  plot.structure.function.coupling.partial.age.effect(
    model = coupling_gamm_output_list[[8]], alpha = .05, axis.labels = TRUE)
# Extract any sensitive periods of development.
significant.age.deriv = plot.sensitive.age.developmental.period(
  model = coupling_gamm_output_list[[8]]$gam, smooth = smooth, parametric = parametric)
# Check if there are any NA entries in the significant_age_deriv plot list
# , which indicates that there were no significant changes in the first
# derivatives of age. 
allplots = list(default.mode.interaction.plot, significant.age.deriv[[1]], significant.age.deriv[[2]])
combined_plot = plot_grid(plotlist = allplots, align = "v", axis = "lr", ncol = 1, rel_heights = c(16, 2, 2))
ggsave(filename = 'data/structure.function/developmental.trends/dmn.coupling.age.effect.png',
       combined_plot, height = 10, width = 10, dpi = 700)

### PART 5 - Explanatory Schematic of Structure-Function Coupling ####
# To explain how structure-function coupling works, we'll create a schematic for
# the first baseline NKI participant, in which we show the within-manifold 
# euclidean distance matrices for communicability and functional connectivity, 
# respectively. First, load these matrices.
nki.baseline.within.manifold.euclid = 
  h5read('data/nki/structure.function/within.manifold.euclidean.distance.bas1.h5', 
         '/within.manifold.euclidean.distances')
# Extract the euclidean distances for the first baseline participant.
nki.bas1.pp1.euclid = nki.baseline.within.manifold.euclid[,,,1]
# Extract the Spearman-rank correlation coefficient for corresponding rows of 
# the structural and functional euclidean distance matrices - this is our 
# coupling measure. 
nroi = 200
nki.bas1.coupling = array(data=NA, dim = c(nroi, 1))
for (roi in 1:nroi){
  nki.bas1.coupling[roi, ] = cor(nki.bas1.pp1.euclid[1,roi,], nki.bas1.pp1.euclid[2,roi,], method = "spearman")
}
# Plot the euclidean distances for each modality and the coupling!
measures = c("structural", "functional", "coupling")
for (measure in measures){
  if (measure == "structural"){
    measure_idx = 1
    visualisation_colour = "#440154FF"
  } else if (measure == "functional"){
    measure_idx = 2
    visualisation_colour = "#35608DFF"
  } 
  df = nki.bas1.pp1.euclid[measure_idx,,] %>% melt()
  colnames(df) = c("X", "Y", "value")
  euclidean.matrix.plot = 
    ggplot(df, aes(x=X, y=Y, fill=value)) + geom_tile() + coord_equal(expand = FALSE) +
    scale_fill_gradient(low = "white", high = visualisation_colour) +
    labs(x = "Regions", y = "Regions") +
    theme(panel.background = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), legend.text = element_blank(),
          axis.title = element_text(size = 15), legend.title = element_blank()) +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    guides(fill = guide_colorbar(ticks.colour = NA))
  ggsave(filename = sprintf('data/structure.function/schematic/%s.euclidean.distance.png', measure),
         plot = euclidean.matrix.plot, dpi = 700, height = 4, width = 4)
}

# Now create a 3D embedding visualization of one candidate region in the 
# structural and functional spaces of the first baseline NKI participant. Load
# the gradients (embeddings) for baseline NKI.
nki.baseline.gradients = read_mat('data/nki/dme/schaefer200x7_bas1.mat')
nki.bas1.pp1.structural.gradient = nki.baseline.gradients[['structural.eigenvectors']][,,1]
nki.bas1.pp1.functional.gradient = nki.baseline.gradients[['functional.eigenvectors']][,,1]
modalities = c("structural", "functional")
for (modality in modalities){
  if (modality == "structural"){
    modality_idx = 1
  } else{
    modality_idx = 2
  }
  df = data.frame(get(sprintf('nki.bas1.pp1.%s.gradient', modality)))
  colnames(df) = c("G1", "G2", "G3")
  # And get the euclidean distance matrix for the ROI of interest (75)
  df$euclid = nki.bas1.pp1.euclid[modality_idx,75,]
  noax = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, 
              showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)
  scope <- kaleido()
  embedding.plot = plot_ly(df, x=~G1, y=~G2, z=~G3, color = I("grey"),
                           type="scatter3d", mode = "markers") %>%
    # Highlight one ROI (50)
    add_markers(data = df, x=~df$G1[75], y=~G2[75], z=~G3[75],
                symbol=I('circle'), inherit=FALSE, color = I("#3CBB75FF")) %>%
    layout(scene = list(xaxis = noax, yaxis=noax, zaxis=noax, aspectmode = 'cube')) %>%
    hide_legend() %>% hide_colorbar()
  scope$transform(embedding.plot,sprintf('data/structure.function/schematic/%s.nki.bas1.pp1.embedding.png',modality)
                  ,width=4,height=4)
}

# Visualise correlation between embedding in 3D structural and functional space 
# for a sample region (75).
df = data.frame(structural.embedding = nki.bas1.pp1.euclid[1,75,], 
                functional.embedding = nki.bas1.pp1.euclid[2,75,])
# Get the X and Y axis breaks!
xaxis = seq(from = min(df$structural.embedding), to = round(max(df$structural.embedding),2), length.out = 6)
yaxis = seq(from = min(df$functional.embedding), to = round(max(df$functional.embedding),2), length.out = 6)
# Format them nicely
xaxis_formatted = formatC(xaxis, digits = 2, format = "f")
yaxis_formatted = formatC(yaxis, digits = 2, format = "f")
# Ensure that 0 is not to 3 digits, but 0!
xaxis_formatted[1] = "0"
yaxis_formatted[1] = "0"
sample.region.coupling.plot = ggplot(data=df, mapping = aes(x = structural.embedding, y = functional.embedding)) +
  geom_point(size = 3, colour = "grey") + geom_smooth(method = "lm", colour = "#3CBB75FF", se = FALSE, linewidth = 2) +
  labs(x = "Structural Embedding", y = "Functional Embedding") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 1),
        axis.text = element_text(size = 20), axis.title = element_text(size = 20, face = "bold"),
        axis.ticks.length = unit(0.2, "cm")) +
  scale_x_continuous(breaks = xaxis, labels = xaxis_formatted) +
  scale_y_continuous(breaks = yaxis, labels = yaxis_formatted)
ggsave(filename = 'data/structure.function/schematic/sample.region.coupling.correlation.plot.png',
       sample.region.coupling.plot, dpi = 700, height = 5, width = 5)

### PART 6 - Structure-Function Coupling Box Plot! ####
# Since structure-function coupling was significantly associated with data set
# across all ICNs and at a global level, we shall plot the distribution of 
# coupling values across data sets and ICNs. First, add global coupling to the 
# summary data frame.
ICN.coupling.df$Global = rowMeans(ICN.coupling.df[,1:7])
ICN.coupling.df.melted = melt(ICN.coupling.df) %>%
  rename(., all_of(c(ICN = "variable", coupling = "value")))
# Now, find the x axis breaks!
xaxis_breaks = seq(from = round(min(ICN.coupling.df.melted$coupling), 2),
                   to = round(max(ICN.coupling.df.melted$coupling), 2), length.out = 5)
coupling.grouped.box.plot = 
  ggplot(data = ICN.coupling.df.melted, mapping = aes(x = ICN, y = coupling, fill = Dataset)) +
  geom_boxplot(mapping = aes(colour = Dataset), fill = "white") +
  labs(y = "Structure-Function Coupling") +
  coord_flip() +
  theme(panel.background = element_blank(), axis.line.x = element_line(colour = "black"),
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.title.x = element_text(size = 20), axis.text.x = element_text(size=20),
        axis.ticks.length.y = unit(0, "cm"), axis.ticks.length.x = unit(0.2, "cm"),
        legend.position = "none") +
  scale_colour_manual(values = c("calm" = "#440154FF", "nki" = "#1F968BFF")) +
  scale_y_continuous(breaks = xaxis_breaks, labels = formatC(xaxis_breaks, digits = 2, format = "f"))
coupling.grouped.box.plot
ggsave(filename = 'data/structure.function/coupling.grouped.box.plot.png',
       plot = coupling.grouped.box.plot, dpi = 700, height = 5, width = 5)

# Report the mean and standard deviation in coupling for each level of analysis...
levels_of_analysis = c("Global", yeo_7_networks)
for (analysis_level in levels_of_analysis){
  subset.df = ICN.coupling.df[[analysis_level]]
  print(sprintf('%s structure-function coupling across both datasets: Mean of %.2f, standard deviation of %.2f.',
          analysis_level, round(mean(subset.df),2), round(sd(subset.df),2)))
}
