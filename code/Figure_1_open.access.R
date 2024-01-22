### PART 1 - Set up the work space ####
rm(list = ls())
library(plotly)
library(raveio)
library(cowplot)
library(reshape2)
library(ggseg3d)
library(ggsegSchaefer)
library(dplyr)
library(ggplot2)
library(reticulate)
setwd('U:/gradients_open_access')
# Add the custom GAMM function code.
source('code/GAMM.functions.R')
modalities = c("structural","functional")
datasets = c("calm","nki")
calm_timepoints = c("baseline","followup")
nnode = 200
ncomp = 3
# And load the region labels for Schaefer 200-node 7-network atlas.
schaefer200x7_metadata = read_mat("data/schaefer200x7_1mm_info.mat")
schaefer200x7_labels = unlist(schaefer200x7_metadata[["schaefer200x7.1mm.info"]][[1]])
nroi = length(schaefer200x7_labels)
# Load the group-level DME outputs and meta-data!
calm_dme_and_metadata_structural = 
  read.csv("data/calm/dme/dme.and.metadata.structural.connectivity.csv")
calm_dme_and_metadata_functional = 
  read.csv("data/calm/dme/dme.and.metadata.functional.connectivity.csv")
nki_dme_and_metadata_structural = 
  read.csv("data/nki/dme/dme.and.metadata.structural.connectivity.csv")
nki_dme_and_metadata_functional = 
  read.csv("data/nki/dme/dme.and.metadata.functional.connectivity.csv")
# And the individual-level outputs...
nki_bas1_dme = read_mat('data/nki/dme/schaefer200x7_bas1.mat')
nki_flu1_dme = read_mat('data/nki/dme/schaefer200x7_flu1.mat')
nki_flu2_dme = read_mat('data/nki/dme/schaefer200x7_flu2.mat')
nki_ses_names = c("bas1","flu1","flu2")

### PART 2 - ICN and parcellation schematics ####
# Create a data frame from the X, Y, and Z coordinates of the Schaefer 200-node
# 7-network parcellation.
schaefer200x7_coord = data.frame(
  X = schaefer200x7_metadata[["schaefer200x7.1mm.info"]][[2]],
  Y = schaefer200x7_metadata[["schaefer200x7.1mm.info"]][[3]],
  Z = schaefer200x7_metadata[["schaefer200x7.1mm.info"]][[4]])
# Add a factor column highlighting which ICN each region belongs to.
icn_assignment = array(NA, nroi)
schaefer200x7_icn = unlist(schaefer200x7_metadata[["schaefer200x7.1mm.info"]][[1]])
for (idx in 1:nroi){
  icn_assignment[idx] = strsplit(schaefer200x7_icn[idx],'_')[[1]][3]
}
schaefer200x7_coord$ICN = factor(icn_assignment)
# Create a vector of colours for plots of each ICN, where each ICN is coloured 
# and the remaining colours are grey. 
icn_colour_vec = c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF", 
                   "#453781FF", "#287D8EFF", "#3CBB75FF")
yeo_7_networks = c("Vis","SomMot","DorsAttn","SalVentAttn","Limbic","Cont","Default")
yeo_7_networks_colour_list = vector("list", 7)
for (icn_idx in 1:7){
  ICN = yeo_7_networks[icn_idx]
  # Create a vector of grey colours
  colour_vec = rep("grey", 7)
  # Find the index of the ICN in the factor levels
  idx = which(levels(schaefer200x7_coord$ICN) %in% ICN)
  # Assign colour to that index
  colour_vec[idx] = icn_colour_vec[idx]
  # Append to output list
  yeo_7_networks_colour_list[[icn_idx]] = colour_vec
}
# Create a 3D plot of the coordinates! Note that we looped through each of the 
# 7 ICNs manually. 
plot_ly(schaefer200x7_coord, x =~X, y=~Y, z=~Z, color = ~ICN, 
        colors = icn_colour_vec) %>% 
  add_markers(showlegend=F, sizes = 10) %>% remove_axes()


### PART 3 - Variability in Structural and Functional Gradients ####
# For each modality, across data sets, create a scree plot with component number
# on the X axis, and proportion of variance explained on the Y axis. Initialize
# a list to hold the data frames with individual-level DME outputs across data
# sets and modalities.
individual_dme_modality_df_list = vector("list", length(modalities))
variability_boxplot_list = vector("list", 2)
for (modality_idx in 1:length(modalities)){
  # Format modality-specific individual-level DME outputs across both data sets,
  # starting with NKI. Average the variance explained over both hemispheres.
  modality = modalities[modality_idx]
  nki_individual_dme_list = vector("list", length(nki_ses_names))
  nki_ses_frequency = vector("numeric", length(nki_ses_names))
  for (nki_ses_idx in 1:length(nki_ses_names)){
    variance_explained = apply(get(paste0("nki_",nki_ses_names[nki_ses_idx],"_dme"))[[paste0(modality,".variance.explained")]],c(1,3), mean)
    nki_individual_dme_list[[nki_ses_idx]] = data.frame(variance_explained)
    nki_ses_frequency[nki_ses_idx] = dim(variance_explained)[1]
  }
  nki_individual_dme = bind_rows(nki_individual_dme_list) %>% 
    rename(., all_of(c("G1_var"="X1", "G2_var"="X2", "G3_var"="X3"))) %>%
    mutate(timepoint = factor(rep(nki_ses_names, nki_ses_frequency)))
  # Combine with CALM into a single data set...
  individual_dme_modality_df = 
    get(paste0("calm_dme_and_metadata_",modality))[c("G1_var", "G2_var", "G3_var", "timepoint")] %>%
    bind_rows(., nki_individual_dme, .id = "dataset") %>%
    mutate(dataset = factor(dataset, labels = c("calm", "nki"))) %>% melt() %>%
    rename(., c("variance" = "value")) %>%
    # Multiply by 100 to get the percentage
    mutate(variance = variance*100)
  # Append to output list!
  individual_dme_modality_df_list[[modality_idx]] = individual_dme_modality_df
  modality = modalities[modality_idx]
  if (modality == "structural"){
    formatted_modality_name = "Structural"
  } else{
    formatted_modality_name = "Functional"
  }
  # Find the minimum and maximum values to be plotted
  yaxis_breaks = seq(from = min(individual_dme_modality_df_list[[modality_idx]]$variance),
                     to = max(individual_dme_modality_df_list[[modality_idx]]$variance),
                     length.out = 5)
  # Now plot as a grouped box plot...
  variance_explained_boxplot = 
    ggplot(data = individual_dme_modality_df_list[[modality_idx]], 
           mapping = aes(x = variable, y = variance, fill = dataset)) +
    geom_boxplot(alpha = 0.5, mapping = aes(color = dataset)) +
    labs(x = "Component", y = "Percentage variance explained") +
    scale_y_continuous(breaks = yaxis_breaks, labels = formatC(yaxis_breaks, digits = 2, format="f")) +
    theme(panel.background = element_blank(), legend.position = "none",
          axis.ticks.length.x = unit(0, "cm"), axis.ticks.length.y = unit(.2, "cm"),
          axis.text.x = element_blank(), axis.text.y = element_text(size = 20),
          axis.line = element_line(colour="black", linewidth = 1),
          axis.title = element_text(size = 20)) +
    scale_colour_manual(values = c("calm" = "#440154FF", "nki" = "#1F968BFF")) +
    scale_fill_manual(values = c("calm" = "#440154FF", "nki" = "#1F968BFF"))
  ggsave(filename = paste0('visualisation/group.gradients/',modality,'.gradients.box.plot.png'),
         plot = variance_explained_boxplot, width = 5, height = 5, dpi = 700)
  # Assign to output list
  variability_boxplot_list[[modality_idx]] = variance_explained_boxplot
}
### PART 5 - Principal Structural and Functional Gradients Remain Dominant Throughout Development ####
# For each modality, plot the relationship between age and variance explained by
# the principal component, as a function of data set.
scatter_plot_list = vector("list", 2)
for (modality_idx in 1:length(modalities)){
  modality = modalities[modality_idx]
  individual_dme_modality_df = individual_dme_modality_df_list[[modality_idx]]
  # For each data set and modality, extract the ages associated with variance
  # explained.
  individual_dme_modality_df_calm = 
    individual_dme_modality_df[which(individual_dme_modality_df$dataset == "calm"),] %>%
    subset(variable == "G1_var") %>%
    rename(., all_of(c('G1_var' = 'variance'))) %>%
    # Divide percentage by 100 to get a proportion so we can merge with meta-data
    mutate(G1_var = G1_var/100) %>%
    # Merge with the meta-data sheet
    merge(., get(paste0("calm_dme_and_metadata_",modality)), by = "G1_var") %>%
    mutate(id = as.character(id)) %>%
    # Convert back to a percentage
    mutate(G1_var = G1_var*100)
  # Now repeat for NKI! 
  individual_dme_modality_df_nki = 
    individual_dme_modality_df[which(individual_dme_modality_df$dataset == "nki"), ] %>%
    subset(variable == "G1_var") %>%
    rename(., all_of(c('G1_var' = 'variance', 'session' = 'timepoint'))) %>%
    mutate(G1_var = G1_var/100) %>%
    merge(., get(paste0("nki_dme_and_metadata_", modality)), by = c("G1_var", "session")) %>%
    mutate(G1_var = G1_var*100)
  # Merge individual_dme_modality_df data frames for both data sets, and create
  # a data set column.
  individual_dme_modality_df = 
    bind_rows(individual_dme_modality_df_calm, individual_dme_modality_df_nki, .id = "dataset") %>%
    mutate(dataset = factor(ifelse(dataset == 1, 'calm', 'nki')))
  # Specify X axis breaks (age)
  xaxis_breaks = seq(from = round(min(individual_dme_modality_df$scan_age)), 
                     to = max(individual_dme_modality_df$scan_age),
                     by = 3)
  # And the Y axis breaks (variance explained)
  yaxis_breaks = seq(from = min(individual_dme_modality_df$G1_var),
                     to = max(individual_dme_modality_df$G1_var),
                     length.out = 4)
  # Create a scatter plot of age (X axis) against variance explained (Y axis)
  # for the first component, and represent the data sets as different coloured
  # lines.
  age_first_component_variance_scatter = 
    ggplot(data = individual_dme_modality_df, 
           mapping = aes(x = scan_age, y = .data[['G1_var']], group=dataset)) +
      geom_point(aes(colour=dataset), size=3, alpha = 0.5) + 
      geom_smooth(method = lm, se=FALSE, aes(colour=dataset), linewidth = 2, fullrange = TRUE) +
      labs(x = "Age (Years)", y = "Component 1 percentage variance") +
      theme(panel.background = element_blank(), axis.line = element_line(colour="black", linewidth = 1),
            axis.ticks.length = unit(.2, "cm"), axis.text = element_text(size = 20),
            axis.title = element_text(size = 20), legend.position = 'none') +
      scale_colour_manual(values = c("#440154FF", "#1F968BFF")) +
      scale_x_continuous(breaks = xaxis_breaks) +
      scale_y_continuous(breaks = yaxis_breaks, labels = 
                           formatC(yaxis_breaks, digits = 2, format = "f"))
  ggsave(filename = sprintf('visualisation/group.gradients/%s.variance.comp.1.age.scatter.png',modality),
         plot = age_first_component_variance_scatter, height = 5, width = 5, dpi = 700)
  scatter_plot_list[[modality_idx]] = age_first_component_variance_scatter
}
### PART 6 - Composite Figure of Parts 4 and 5! ####
composite_variability_age_and_component_plot = 
  plot_grid(plotlist = c(variability_boxplot_list, scatter_plot_list))
ggsave('visualisation/group.gradients/composite.variability.age.and.component.plot.png',
       composite_variability_age_and_component_plot, height = 10, width = 10, dpi = 700)

### PART 7 - Visualize Grey-Scale 3D Brain for fMRI Visualisation ####
# Generate 50 hues of grey
example_df = data.frame(
  colours_col = rep(c(scales::alpha("grey", seq(from = .25, to = 1, length.out = 50))),4),
  region = schaefer200x7_labels)
example_grey_brain_visualisation = 
  ggseg3d(.data = example_df, atlas = schaefer7_200_3d, hemisphere = "left",
        colour = "colours_col", surface = "inflated", show.legend = F) %>%
  remove_axes() %>% pan_camera("left medial")
scope$transform(example_grey_brain_visualisation, file = 
                  'visualisation/group.gradients/example.grey.scale.brain.png',
                height = 4, width = 4, dpi = 700)
### PART 8 - Visualize Schaefer 200-Node 7-Network Parcellation ####
# We want to visualise each of Yeo's 7 networks in a different colour from the 
# viridis colour palette. 
yeo_7_networks = c("Vis","SomMot","DorsAttn","SalVentAttn","Limbic","Cont","Default")
# Find which regions belong to which regions
yeo_7_networks_assignment = vector(length = 200)
yeo_7_networks_colours = vector(length = 200)
for (idx in 1:200){
  yeo_7_networks_assignment[idx] = strsplit(schaefer200x7_labels[idx],'_')[[1]][3]
  # Assign colours to each network...
  if (yeo_7_networks_assignment[idx] == "Vis"){
    yeo_7_networks_colours[idx] = alpha("#440154FF", 0.5)
  } else if (yeo_7_networks_assignment[idx] == "SomMot"){
    yeo_7_networks_colours[idx] = alpha("#39568CFF", 0.5)
  } else if (yeo_7_networks_assignment[idx] == "DorsAttn"){
    yeo_7_networks_colours[idx] = alpha("#1F968BFF", 0.5)
  } else if (yeo_7_networks_assignment[idx] == "SalVentAttn"){
    yeo_7_networks_colours[idx] = alpha("#73D055FF", 0.5)
  } else if (yeo_7_networks_assignment[idx] == "Limbic"){
    yeo_7_networks_colours[idx] = alpha("#453781FF", 0.5)
  } else if (yeo_7_networks_assignment[idx] == "Cont"){
    yeo_7_networks_colours[idx] = alpha("#287D8EFF", 0.5)
  } else{
    yeo_7_networks_colours[idx] = alpha("#3CBB75FF",0.5)
  }
}
yeo_7_network_schaefer_ggseg3d_df = 
  data.frame(region = schaefer200x7_labels,
             assignment = factor(yeo_7_networks_assignment),
             colours = yeo_7_networks_colours) 
yeo_7_network_schaefer_ggseg3d_left_lateral = 
  ggseg3d(.data = yeo_7_network_schaefer_ggseg3d_df, atlas = schaefer7_200_3d, 
          hemisphere = "left", colour = "colours", 
          surface = "inflated", show.legend = F) %>%
  remove_axes() %>% pan_camera("left lateral")
scope$transform(yeo_7_network_schaefer_ggseg3d_left_lateral,
                file = "visualisation/group.gradients/schaefer200_7_yeo_left_lateral.png",
                width = 5, height = 5)



### PART 9 - Visualization of Group-Level Gradients and Individual Variation ####
# Plot group-representative gradients on the cortical surface, for each data 
# set and modality. We use the ggseg3d package to plot 2 views (medial and 
# lateral) for each hemisphere in each modality, for the principal 3 gradients. 
view_list = c("left medial","right lateral", "left lateral", "right medial")
modality_list = c("structural","functional")
# Load the group-level gradients.
group_gradients = read_mat('data/group.gradients/calm.nki.group.gradients.mat')
# Extract group-gradients for each data set. Note, we select the referred group
# in CALM.
calm_group_gradients = group_gradients[["calm"]][[1]]
nki_group_gradients = group_gradients[["nki"]]
all_eigenvectors = list(calm_group_gradients, nki_group_gradients)
# To ensure valid comparisons between data sets and modalities, flip the left
# hemisphere for the CALM structural gradients, and flip both hemispheres for 
# the first functional gradient in CALM.
all_eigenvectors[[1]][1,1:100,1:2] = -1*all_eigenvectors[[1]][1,1:100,1:2]
all_eigenvectors[[1]][2,,1] = -1 * all_eigenvectors[[1]][2,,1]
# Create the section breaks for the colour palette from viridis
eigenvector_quantiles = quantile(unlist(all_eigenvectors), probs = seq(.1,.9, by = .1))
colour_map = c("#440154FF" = min(unlist(all_eigenvectors)),
               "#482173FF" = eigenvector_quantiles[[1]],
               "#414487FF" = eigenvector_quantiles[[2]],
               "#35608DFF" = eigenvector_quantiles[[3]],
               "#2A788EFF" = eigenvector_quantiles[[4]],
               "#21908CFF" = eigenvector_quantiles[[5]],
               "#22A884FF" = eigenvector_quantiles[[6]],
               "#43BF71FF" = eigenvector_quantiles[[7]],
               "#7AD151FF" = eigenvector_quantiles[[8]],
               "#BBDF27FF" = eigenvector_quantiles[[9]],
               "#FDE725FF" = max(unlist(all_eigenvectors)))
view_list = c("left medial","right lateral", "left lateral", "right medial")
# Now loop through each data set, plotting manifolds for each modality.
for (dataset_idx in 1:length(datasets)){
  dataset = datasets[dataset_idx]
  eigenvectors = all_eigenvectors[[dataset_idx]]
  for (modality_idx in 1:length(modalities)){
    for (component in 1:3){
      gradient_df = tibble(eigenvectors = eigenvectors[modality_idx,,component],
                           region = schaefer200x7_labels)
      for (view_idx in 1:length(view_list)){
        scope <- kaleido()
        p = ggseg3d(.data=gradient_df,atlas=schaefer7_200_3d,
                    colour="eigenvectors",surface="inflated",
                    # This is from the viridis palette 
                    palette=colour_map, show.legend = T, 
                    hemisphere = strsplit(view_list[view_idx], " ")[[1]][1]) %>%
          pan_camera(view_list[view_idx]) %>%
          remove_axes()
        # Set the file name to save the plot!
        save_filename = paste0(
          "visualisation/group.gradients/",dataset,".group.",
          modalities[modality_idx],".manifold.", component,".",
          view_list[view_idx],".png")
        scope$transform(p,save_filename,width=4,height=4)
        rm(p,save_filename)
        print(paste("Plotted the",view_list[view_idx], "surface for", 
                    modalities[modality_idx], "component", component, 
                    "in the", dataset, "data set."))
        rm(scope);gc()
      }
    } 
  } 
}
