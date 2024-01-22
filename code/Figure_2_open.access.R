### PART 1 - Setting up the Work Space ####
rm(list = ls())
library(raveio)
library(plotly)
library(viridis)
library(scales)
library(dplyr)
library(reshape2)
setwd('U:/gradients_open_access/')
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

### PART 2 - Manifold Eccentricity 3D Schematic ####
# Load the group gradients for NKI, and then the structural eccentricity for the
# first participant in NKI (sub-A00018030), and structural gradients
nki_group_gradients = read_mat("data/calm.nki.group.gradients.mat")[["nki"]]
# Find the group-level manifold origin, by calculating the intersection of the 
# first 3 gradients for each hemisphere, and then the mean of that. Format into
# a data frame.
nki_structural_manifold_origin = 
  c((colMeans(nki_group_gradients[1,,]) + colMeans(nki_group_gradients[2,,])) /2) %>%
  t() %>% as.data.frame() %>% `colnames<-`(paste0("G",1:3))
# Load the first 3 structural eigen-vectors for the first participant in NKI 
# (sub-A00018030), merge with their eccentricity values, and format into a data
# frame.
nki_baseline_structural_gradients = 
  read_mat("data/nki/dme/schaefer200x7_bas1.mat")[["structural.eigenvectors"]][,,1] %>%
  as.data.frame() %>% 
  bind_cols(read_mat("data/nki/dme/schaefer200x7_bas1.mat")[["structural.manifold.eccentricity"]][1,]) %>%
  `colnames<-`(c(paste0("G",1:3),"eccentricity")) %>%
  # And add regional brain labels
  mutate(region = schaefer200x7_labels) %>%
  # Create a continuous colour scale to visualise the eccentricities.
  mutate(col = viridis(n=200, begin = min(eccentricity), end = max(eccentricity)))
# Plot on a 3D axis!
noax = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, 
            showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)
plot_ly(nki_baseline_structural_gradients, x=~G1, y=~G2, z=~G3, 
        color=~eccentricity, type="scatter3d", mode = "markers") %>%
  # Add the group manifold origin
  add_markers(data = nki_structural_manifold_origin, x=~G1, y=~G2, z=~G3,
            symbol=I('square'), inherit=FALSE, color = I("red")) %>%
  layout(scene = list(xaxis = noax, yaxis=noax, zaxis=noax, aspectmode = 'cube')) %>%
  hide_legend() %>% hide_colorbar()

### PART 3 - Nodal Coefficient of Variation in Eccentricity ####
# Start by initializing an array to hold the coefficient of variation for nodal
# eccentricity across data sets and modalities. 
cv_array = array(NA,dim=c(length(datasets), length(modalities), nnode))
for (dataset_idx in 1:length(datasets)){
  dataset = datasets[dataset_idx]
  for (modality_idx in 1:length(modalities)){
    modality = modalities[modality_idx]
    # Get the nodal eccentricities for this data set and modality
    dme_and_metadata_df = get(paste0(dataset,'_dme_and_metadata_',modality))
    for (node in 1:nnode){
      nodal_eccentricity = dme_and_metadata_df[,sub('7Networks_*','',schaefer200x7_labels[node])]
      cv_array[dataset_idx,modality_idx,node] = sd(nodal_eccentricity)/mean(nodal_eccentricity)
    }
  }
}

# Now plot the distributions of nodal coefficients of variation (CVs). Start by 
# creating one large data frame with CVs across data sets and modalities. 
cv_df = cv_array %>% melt(value.name = "eccentricity.cv") %>%
  rename(., all_of(c(dataset = "Var1", modality = "Var2", node = "Var3"))) %>%
  mutate_at(c("dataset", "modality"), factor) %>% 
  # Create new group variable to represent modality and data set combinations. 
  mutate(Group = case_when(dataset == 1 & modality == 1 ~ "CALM Structural",
                           dataset == 1 & modality == 2 ~ "CALM Functional",
                           dataset == 2 & modality == 1 ~ "NKI Structural",
                           dataset == 2 & modality == 2 ~ "NKI Functional")) %>%
  # And convert this to a factor
  mutate_at("Group", factor)

# We will plot the coefficient of variation for CALM as purple, and NKI as
# green. A darker colour represents structure, and lighter represents function. 
group_colours = c(scales::alpha("#440154FF", c(0.5, 0.25)), scales::alpha("#1F968BFF", c(0.5, 0.25)))
# Set the x axis breaks
breaks = seq(from = 0, to = 0.55, by = .11)
# To create a vertical line to indicate the mean for each group, terminating at
# where the density plot peaks, we need to format the data. 
group_levels = unique(cv_df$Group)
group_eccentricity_list = vector("list", length=length(group_levels))
for (level_idx in 1:length(group_levels)){
  level = group_levels[level_idx]
  # Subset cv_df by this grouping level
  cv_df_subset = filter(cv_df, Group == level)
  # Calculate the mean eccentricity CV for this level 
  mean_eccentricity_cv = mean(cv_df_subset$eccentricity.cv)
  # Find the X index in the density smooth corresponding most closely to the 
  # mean non-smoothed eccentricity
  x_idx = which.min(abs(density(cv_df_subset$eccentricity.cv)$x - mean_eccentricity_cv))
  # And find the corresponding Y value in the density smooth
  y_val = density(cv_df_subset$eccentricity.cv)$y[x_idx]
  # Create a data frame holding these values, and append to output list.
  group_eccentricity_cv_df = data.frame(
    'Group' = level, 'mean.cv' = mean_eccentricity_cv, 'corresponding_density_y_val' = y_val)
  group_eccentricity_list[[level_idx]] = group_eccentricity_cv_df
}
group_eccentricity = bind_rows(group_eccentricity_list)

# Plot the coefficient of variation for communicability as purple, and 
# functional connectivity as green. Use a lighter shade for CALM, and darker for
# NKI.
nodal_manifold_cv_distributions = 
  ggplot(data=cv_df, mapping=aes(x=eccentricity.cv, fill=Group, color=Group)) + 
  geom_density(position = "identity") +
  labs(x = "Coefficient of variation for nodal manifold eccentricity") +
  theme(panel.background = element_blank(), axis.line.x = element_line(colour="black"),
        axis.text.x = element_text(size = 25), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 25, face="bold", hjust=0.5),
        legend.position = "none", axis.ticks.length = unit(.2, "cm")) +
  # Remove space between histogram bars and x axis
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(limits = c(0, 0.55), breaks = breaks,
                     labels = formatC(breaks, digits = 2, format = "f")) +
  # Add a vertical line for each group to indicate the mean CV
  geom_segment(data = group_eccentricity, mapping = aes(
    x = mean.cv, xend = mean.cv, y = 0, yend = corresponding_density_y_val, group=Group, colour=Group),
    linetype="dashed", lwd=1.5) +
  scale_fill_manual(values = group_colours) + scale_color_manual(values = group_colours)
nodal_manifold_cv_distributions
ggsave(filename = 'visualisation/nodal.manifold.eccentricity.cv.distribution.png',
       plot = nodal_manifold_cv_distributions, width = 15, height = 8, dpi = 700)
# Also save the distribution of coefficients of variation so that we can perform
# statistical tests on them in Python!
write.csv(cv_df, 'data/coefficient.of.variation.csv', row.names = FALSE)
### PART 4 - Relationship Between Manifold Eccentricity and Graph Theory Metrics ####
# To assess what aspect of brain organisation manifold eccentricity measures, we
# calculated several graph theory metrics for each connectome. We plot the 
# nodal relationships, with each scatter plot having one line per data set, and 
# individual points coloured by data set.  

# Load the nodal graph theory metrics for both modalities and data sets.
nodal_graph_theory_list = vector("list", 2)
nroi = 200
for (modality_idx in 1:length(modalities)){
  dataset_df_list = vector("list", 2)
  for (dataset_idx in 1:length(datasets)){
    dataset_df_list[[dataset_idx]] = read.csv(
      sprintf('data/%s/connectomes/nodal.%s.graph.theory.metrics.csv', 
              datasets[dataset_idx],modalities[modality_idx]))[,-1]
  }
  datasets_df = bind_rows(dataset_df_list)
  datasets_df$dataset = factor(rep(datasets, each=nroi))
  nodal_graph_theory_list[[modality_idx]] = datasets_df
}
nodal_graph_theory_df = bind_rows(nodal_graph_theory_list)
nodal_graph_theory_df$modality = factor(rep(modalities, each=nroi*2))
# Find the metrics!
metrics = setdiff(colnames(nodal_graph_theory_df), c("modality", "dataset", "Mean.Manifold.Eccentricity"))
# Loop across each modality, and plot the relationship between each metric and
# mean (nodal) manifold eccentricity as a function of data set.
for (modality in modalities){
  # Subset by modality
  modality_nodal_graph_theory_df = nodal_graph_theory_df[which(
    nodal_graph_theory_df$modality == modality),]
  # Set Y axis breaks to be the same across graphs for the same modality. 
  modality_breaks = 
    seq(min(modality_nodal_graph_theory_df$Mean.Manifold.Eccentricity),
        max(modality_nodal_graph_theory_df$Mean.Manifold.Eccentricity),
        length.out = 5)
  # Loop across metrics...
  for (metric_idx in 1:length(metrics)){
    metric = metrics[metric_idx]
    # Set X axis breaks to be the same across graphs for the same metric.
    metric_mean_manifold_eccentricity = nodal_graph_theory_df[[metric]]
    metric_breaks = 
      seq(min(metric_mean_manifold_eccentricity),
          max(metric_mean_manifold_eccentricity), length.out = 5)
    # Set the correct X axis label
    if (metric == "Local.Modularity"){
      xaxis_label = "Local Modularity"
    } else if (metric == "Clustering.Coef"){
      xaxis_label = "Clustering Coefficient"
    } else if (metric == "Participation.Coef"){
      xaxis_label = "Participation Coefficient"
    } else {
      xaxis_label = expression(paste(bold("Within-Module Degree "), bolditalic("Z"), bold("-Score")))
    }
    # For the first metric for each modality, we'll include axis labels and 
    # ticks. For all others, we'll omit. 
    if (metric_idx == 1 & modality == "structural" | metric_idx == 1 & modality == "functional"){
      xaxis.text = element_text(size = 30, colour = "black")
      yaxis.text = element_text(size = 30, colour = "black")
      xaxis.title = element_text(size = 25, face = "bold")
      yaxis.title = element_text(size = 25, face = "bold")
      xaxis.ticks.length = unit(.2, "cm")
      yaxis.ticks.length = unit(.2, "cm")
    } else if (metric_idx != 1 & modality == "functional"){
      xaxis.text = element_text(size = 30, colour = "black")
      yaxis.text = element_text(size = 30, colour = "white")
      xaxis.title = element_text(size = 25, face = "bold")
      yaxis.title = element_text(size = 25, colour = "white")
      xaxis.ticks.length = unit(.2, "cm")
      yaxis.ticks.length = unit(0, "cm")
    } else {
      xaxis.text = element_text(size = 30, colour = "white")
      yaxis.text = element_text(size = 30, colour = "white")
      xaxis.title = element_text(size = 25, colour = "white")
      yaxis.title = element_text(size = 25, colour = "white")
      xaxis.ticks.length = unit(0, "cm")
      yaxis.ticks.length = unit(0, "cm")
    }
    # Subset the modality data frame by this metric
    nodal_df = modality_nodal_graph_theory_df %>% 
      dplyr::select(all_of(c(metric, "Mean.Manifold.Eccentricity", "dataset")))
    # Now plot!
    nodal_metric_manifold_eccentricity_relationship_plot = 
      ggplot(nodal_df, aes(x=.data[[metric]], y=Mean.Manifold.Eccentricity, group=dataset)) + 
      geom_point(aes(color=dataset), alpha = .5, size = 5) + 
      geom_smooth(aes(group=dataset, colour=dataset), method=lm, se=FALSE, fullrange = FALSE, linewidth=2) +
      labs(x = xaxis_label, y = "Mean Manifold Eccentricity") +
      theme(legend.position = "none", panel.background = element_blank(), 
            axis.line = element_line(colour="black", linewidth = 1), 
            axis.text.x = xaxis.text, axis.text.y = yaxis.text, 
            axis.title.x = xaxis.title, axis.title.y = yaxis.title,
            axis.ticks.length.x = xaxis.ticks.length, axis.ticks.length.y = yaxis.ticks.length) +
      scale_color_manual(values = c("#440154FF", "#1F968BFF")) +
      scale_y_continuous(breaks = modality_breaks, labels = formatC(modality_breaks, 2, format="f")) +
      scale_x_continuous(breaks = metric_breaks, labels = formatC(metric_breaks, 2, format="f")) +
      coord_cartesian(xlim = c(min(metric_breaks), max(metric_breaks)))
    ggsave(file = sprintf('visualisation/%s.nodal.%s.manifold.eccentricity.relationship.png', metric, modality),
           plot = nodal_metric_manifold_eccentricity_relationship_plot, height = 9, width = 9, dpi = 700)
  }
}
