# This script compiles diffusion map embedding outputs with meta-data for CALM
# and NKI, across time points and modalities. For open-access, we provide the
# code for creating the DME summary data frames. However, due to data access 
# restrictions, we cannot provide the meta-data sheets for NKI, as this requires
# a data usage agreement. Further, with CALM, we include minimal meta-data, 
# including derived meta-data such as head motion. 
### PART 1 - Setting up the Work Space ####
rm(list = ls())
library(raveio)
library(dplyr)
setwd('U:/gradients_open_access')
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
### PART 2A - Load CALM Meta-Data ####
calm_sc_metadata = read.csv('data/calm/calm_sc_covariates.csv')
calm_fc_metadata = read.csv('data/calm/calm_fc_covariates.csv')
# For each of the meta-data sheets...
calm_metadata = list(calm_sc_metadata,calm_fc_metadata)
for (modality_idx in 1:length(modalities)){
  # Format categorical variables into factors
  calm_metadata[[modality_idx]][,c("BIDS","ID","Sex","scannerID","timepoint")] =
    lapply(calm_metadata[[modality_idx]][,c("BIDS","ID","Sex","scannerID","timepoint")],factor)
  # Convert all column names to lowercase
  colnames(calm_metadata[[modality_idx]]) = tolower(colnames(calm_metadata[[modality_idx]]))
  # Re-code the time point column to a character factor
  calm_metadata[[modality_idx]]$timepoint = recode_factor(calm_metadata[[modality_idx]]$timepoint, `0` = "baseline", `1` = "followup")
  # And repeat for scanner ID
  calm_metadata[[modality_idx]]$scannerid = recode_factor(calm_metadata[[modality_idx]]$scannerid, `0` = "before", `1` = "after")
}

### PART 2B - Format CALM Meta-Data with DME Outputs ####
# Load the DME outputs (eccentricities, variance explained etc.) for each time
# point, and for our main analysis parcellation (Schaefer 200-node 7-network).
calm_baseline_dme = read_mat('data/calm/dme/referred_schaefer200x7_baseline.mat')
calm_followup_dme = read_mat('data/calm/dme/referred_schaefer200x7_followup.mat')
for (modality_idx in 1:length(modalities)){
  modality = modalities[modality_idx]
  # Select the appropriate meta data sheet. Note that the DME outputs may 
  # contain fewer participants than those in the meta data sheet because some 
  # participants may've been excluded in QC. 
  metadata = calm_metadata[[modality_idx]]
  # Initialize a list to hold DME data frames for each time point
  dme_df_list_modality = vector("list",length(calm_timepoints))
  for (timepoint_idx in 1:length(calm_timepoints)){
    # Subset meta-data by timepoint 
    metadata_subset = subset(metadata,timepoint == calm_timepoints[timepoint_idx])
    # Extract the participants with DME data 
    sub_list = c(get(paste0("calm_",calm_timepoints[timepoint_idx],"_dme"))[[paste0(modality,".sub.list")]])
    nsub = length(sub_list)
    # Remove space at end of each subject ID
    sub_list = gsub(" ","",sub_list)
    # For each time point for this modality, extract the manifold eccentricity 
    # per node and subject lists.
    nodal_df = data.frame(sub_list, get(paste0("calm_",calm_timepoints[timepoint_idx],"_dme"))[[paste0(modality,".manifold.eccentricity")]])
    colnames(nodal_df) = c("bids",schaefer200x7_labels)
    # Average manifold eccentricity across nodes for each participant
    nodal_df$global_manifold_eccentricity = rowMeans(nodal_df %>% dplyr::select(starts_with("7Networks")))
    # Extract the variance explained by each component, and average across 
    # hemispheres for each participant
    variance_explained = get(paste0("calm_",calm_timepoints[timepoint_idx],"_dme"))[[paste0(modality,".variance.explained")]]
    variance_explained_array = array(NA,dim=c(nsub,ncomp))
    for (sub_idx in 1:nsub){
      variance_explained_array[sub_idx,] = colMeans(variance_explained[sub_idx,,])
    }
    # Add sub_list and make into a data frame
    variance_explained_df = data.frame(sub_list,variance_explained_array)
    colnames(variance_explained_df) = c("bids","G1_var","G2_var","G3_var")
    # Merge variance_explained_df with nodal_df, attach associated meta-data, 
    # and assign to output list
    dme_df_list_modality[[timepoint_idx]] = merge(merge(variance_explained_df, nodal_df),metadata_subset)
  }
  # Bind rows for dme_df_list_modality, and assign to output list
  dme_df_list_modality_bound = bind_rows(dme_df_list_modality)
  # After collecting data for both time points, assign to modality output list.
  if (modality == "structural"){
    calm_dme_and_metadata_structural = dme_df_list_modality_bound
  } else{
    calm_dme_and_metadata_functional = dme_df_list_modality_bound
  }
}

### PART 2C - Load NKI Meta-Data ####
# Load DME outputs for each time point
nki_bas1_dme = read_mat('data/nki/dme/schaefer200x7_bas1.mat')
nki_flu1_dme = read_mat('data/nki/dme/schaefer200x7_flu1.mat')
nki_flu2_dme = read_mat('data/nki/dme/schaefer200x7_flu2.mat')
# Load motion estimates (mean frame-wise displacement/FWD). For each session, 
# create a new data frame with subject, session, and mean FWD for SC and FC.
nki_motion_estimates = read_mat('data/nki/motion_estimates.mat')[['fd.metadata']]
nki_ses_names = c("bas1","flu1","flu2")
# For each modality, create a new data frame with id, session, and mean FWD.
for (modality_idx in 1:length(modalities)){
  nki_modality_motion_estimates = nki_motion_estimates[[modality_idx]]
  motion_estimate_df = data.frame(
    meanfwd = unlist(nki_modality_motion_estimates[[1]]),
    id = unlist(nki_modality_motion_estimates[[2]]),
    session = c(rep("bas1",length(unlist(nki_modality_motion_estimates[[2]][[1]]))),
                rep("flu1",length(unlist(nki_modality_motion_estimates[[2]][[2]]))),
                rep("flu2",length(unlist(nki_modality_motion_estimates[[2]][[3]])))))
  # Name appropriately!
  if (modalities[modality_idx] == "structural"){
    nki_motion_estimates_sc = motion_estimate_df
  } else{
    nki_motion_estimates_fc = motion_estimate_df
  }
  rm(motion_estimate_df)
}
# Load the neuroimaging participants data sheet (minimal meta-data). For data 
# privacy reasons, we cannot provide this as open-access.
nkir_neuroimaging_participants = 
  read.table("//cbsu/data/imaging/projects/external/nkir/data/enhanced_nkir/participants.tsv",header = T)
# Retain ID, session, age, and gender from the neuro-imaging data sheet
nkir_neuroimaging_participants = nkir_neuroimaging_participants %>%
  dplyr::select(c("subject","session","gender","age")) %>%
  distinct() %>%
  # Add a 'sub-' prefix for subjects for easier merging with other data frames
  mutate(subject = paste0('sub-',subject)) %>%
  # Convert categorical variables to factors
  mutate_at(c("subject","session","gender"),as.factor) %>%
  # And convert the session variables to lower case
  mutate(session = tolower(session)) %>%
  # Rename subject to id, and gender to sex, in line with other data frames
  rename("id"="subject","sex"="gender")

### PART 2D - Format NKI Meta-Data and DME Outputs ####
# For NKI, find how many participants have one, two, or three sessions. To 
# prepare for the models assessing age effects on DME outputs, create meta-data
# data frames for each number of sessions with ID, sex, mean FWD, and a list 
# with the nodal eccentricities.
for (modality_idx in 1:length(modalities)){
  modality = modalities[modality_idx]
  # For this modality, and across all 3 sessions, load the list of participants
  # with DME data.
  nki_dme_pps_list = list(nki_bas1_dme[[paste0(modalities[modality_idx],".sub.list")]],
                          nki_flu1_dme[[paste0(modalities[modality_idx],".sub.list")]],
                          nki_flu2_dme[[paste0(modalities[modality_idx],".sub.list")]])
  # Select the appropriate motion estimates data frame
  if (modalities[modality_idx] == "structural"){
    motion_estimates_df = nki_motion_estimates_sc
  } else{
    motion_estimates_df = nki_motion_estimates_fc
  }
  # For each session and each participant, create a data frame with nodal DME 
  # outputs, variance explained by each component, and associated meta-data.
  dme_and_metadata_list = vector("list",length(nki_ses_names))
  for (nki_ses_idx in 1:length(nki_ses_names)){
    nki_session = nki_ses_names[nki_ses_idx]
    # Create a data frame with nodal DME values and ID
    nodal_dme_df = data.frame(nki_dme_pps_list[[nki_ses_idx]], 
                              get(paste0('nki_',nki_session,'_dme'))[[paste0(modalities[modality_idx],".manifold.eccentricity")]])
    # Format column names with ID and Schaefer 200-node parcellation labels
    colnames(nodal_dme_df) = c("id",schaefer200x7_labels)
    # Average manifold eccentricity across nodes to get a global measure
    nodal_dme_df$global_manifold_eccentricity = rowMeans(nodal_dme_df %>% dplyr::select(starts_with("7Networks")))
    # Load the variance explained by each DME component for each participant.
    variance_explained = get(paste0('nki_',nki_session,'_dme'))[[paste0(modalities[modality_idx],".variance.explained")]]
    # Average across hemispheres for each component and participant
    nsub = dim(variance_explained)[1]
    variance_explained_array = array(NA,dim=c(nsub,ncomp))
    for (sub_idx in 1:nsub){
      variance_explained_array[sub_idx,] = colMeans(variance_explained[sub_idx,,])
    }
    # Add subject list, and format column names
    variance_explained_df = data.frame(nki_dme_pps_list[[nki_ses_idx]],variance_explained_array)
    colnames(variance_explained_df) = c("id","G1_var","G2_var","G3_var")
    # Merge with nodal_dme_df
    variance_and_nodal_dme = merge(nodal_dme_df, variance_explained_df)
    # Now move onto the meta-data. First, subset the neuroimaging meta-data by
    # this session, and merge with motion estimates for participants with high-
    # quality neuroimaging data in this session. We then add the DME out-puts,
    # and assign to the output list.
    dme_and_metadata_list[[nki_ses_idx]] = nkir_neuroimaging_participants %>% 
      subset(session == nki_session) %>% merge(motion_estimates_df[which(motion_estimates_df$session == nki_session),],all.y = T) %>%
      # Now merge with variance_and_nodal_dme
      merge(variance_and_nodal_dme, by = "id", all.y = T) %>%
      # And convert session to a factor
      mutate(session = factor(session)) %>%
      # And change 'age' to 'scan_age', as we retrieved this from the neuro-
      # imaging data sheet.
      rename(scan_age = age)
  }
  # After looping through all sessions for this modality, bind the rows of the 
  # dme and meta-data lists to get one complete data frame, and assign to the 
  # relevant output.
  if (modality == "structural"){
    nki_dme_and_metadata_structural = bind_rows(dme_and_metadata_list)
  } else{
    nki_dme_and_metadata_functional = bind_rows(dme_and_metadata_list)
  }
}

### PART 2E - Saving Formatted Meta-Data and DME Outputs ####
calm_dme_and_metadata_structural$id = as.character(calm_dme_and_metadata_structural$id)
names(calm_dme_and_metadata_structural) = gsub('7Networks_*', replacement = '', x = names(calm_dme_and_metadata_structural))
write.csv(calm_dme_and_metadata_structural, file = "data/calm/dme/dme.and.metadata.structural.connectivity.csv",
          quote = FALSE)
calm_dme_and_metadata_functional$id = as.character(calm_dme_and_metadata_functional$id)
names(calm_dme_and_metadata_functional) = gsub('7Networks_*', replacement = '', x = names(calm_dme_and_metadata_functional))
write.csv(calm_dme_and_metadata_functional, file = "data/calm/dme/dme.and.metadata.functional.connectivity.csv",
          quote = FALSE)
names(nki_dme_and_metadata_structural) = gsub('7Networks_*', replacement = '', x = names(nki_dme_and_metadata_structural))
write.csv(nki_dme_and_metadata_structural, file = "data/nki/dme/dme.and.metadata.structural.connectivity.csv",
          quote = FALSE)
names(nki_dme_and_metadata_functional) = gsub('7Networks_*', replacement = '', x = names(nki_dme_and_metadata_functional))
write.csv(nki_dme_and_metadata_functional, file = "data/nki/dme/dme.and.metadata.functional.connectivity.csv",
          quote = FALSE)


### PART 3 - Descriptive Statistics for Paper ####
### PART 3A - Main Analysis ####
### NKI ####
for (modality in modalities){
  graph_theory_metrics = read.csv(sprintf('data/nki/connectomes/participant.%s.graph.theory.metrics.csv',modality))
  metadata = get(paste0('nki_dme_and_metadata_',modality))
  names(graph_theory_metrics)[names(graph_theory_metrics) == 'Timepoint'] <- 'session'
  names(graph_theory_metrics)[names(graph_theory_metrics) == 'subid'] <- 'id'
  merged_df = merge(graph_theory_metrics, metadata, by = c('session', 'id'), all.x = TRUE)
  merged_df$session = factor(merged_df$session)
  age_ranges_by_session = tapply(merged_df$scan_age, merged_df$session, range)
  age_mean_by_session = tapply(merged_df$scan_age, merged_df$session, mean)
  print(sprintf('For %s analyses in NKI, we examined %.1f children across %.1f timepoints, 
        spanning %.2f to %.2f years old [%.2f percent male]. %.1f scans were at baseline [Mean age of %.2f, range of %.2f and %.2f],
        %.1f second-scans [Mean age of %.2f, range of %.2f and %.2f years],
        and %.1f with three scans [Mean age of %.2f, range of %.2f and %.2f years].', modality, nrow(merged_df), 
          nlevels(merged_df$session), min(merged_df$scan_age), max(merged_df$scan_age),
          (summary(merged_df$sex)[['M']]/nrow(merged_df))*100,
          summary(merged_df$session)[[1]], age_mean_by_session[[1]], age_ranges_by_session[[1]][1], age_ranges_by_session[[1]][2],
          summary(merged_df$session)[[2]], age_mean_by_session[[2]], age_ranges_by_session[[2]][1], age_ranges_by_session[[2]][2],
          summary(merged_df$session)[[3]], age_mean_by_session[[3]], age_ranges_by_session[[3]][1], age_ranges_by_session[[3]][2]))
}
### CALM ####
for (modality in modalities){
  graph_theory_metrics = read.csv(sprintf('data/calm/connectomes/participant.%s.graph.theory.metrics.csv',modality))
  metadata = get(paste0('calm_dme_and_metadata_',modality))
  names(graph_theory_metrics)[names(graph_theory_metrics) == 'Timepoint'] <- 'timepoint'
  names(graph_theory_metrics)[names(graph_theory_metrics) == 'subid'] <- 'bids'
  merged_df = merge(graph_theory_metrics, metadata, by = c('timepoint', 'bids'), all.x = TRUE)
  merged_df$timepoint = factor(merged_df$timepoint)
  age_ranges_by_session = tapply(merged_df$scan_age, merged_df$timepoint, range)
  age_mean_by_session = tapply(merged_df$scan_age, merged_df$timepoint, mean)
  print(sprintf('For %s analyses in CALM, we examined %.1f children across %.1f timepoints, 
        spanning %.2f to %.2f years old [%.2f percent male]. %.1f scans were at baseline [Mean age of %.2f, range of %.2f and %.2f],
        %.1f second-scans [Mean age of %.2f, range of %.2f and %.2f years].', modality, nrow(merged_df), 
                nlevels(merged_df$timepoint), min(merged_df$scan_age), max(merged_df$scan_age),
                (summary(merged_df$sex)[['M']]/nrow(merged_df))*100,
                summary(merged_df$timepoint)[[1]], age_mean_by_session[[1]], age_ranges_by_session[[1]][1], age_ranges_by_session[[1]][2],
                summary(merged_df$timepoint)[[2]], age_mean_by_session[[2]], age_ranges_by_session[[2]][1], age_ranges_by_session[[2]][2]))

}

