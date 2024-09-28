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
library(readxl)
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
### PART 2 - Formatting CALM Meta-Data with DME Outputs ####
# Load the DME outputs (eccentricities, variance explained etc.) for each time
# point, and for our main analysis parcellation (Schaefer 200-node 7-network).
calm_baseline_dme = read_mat('data/calm/dme/referred_schaefer200x7_baseline.mat')
calm_followup_dme = read_mat('data/calm/dme/referred_schaefer200x7_followup.mat')
# Get the CALM meta-data
calm_metadata_sc = read_xlsx('//cbsu/data/imaging/projects/external/nkir/analyses/Alicja_calm_structural_master.xlsx')
calm_metadata_fc = read_xlsx('//cbsu/data/imaging/projects/external/nkir/analyses/Alicja_calm_functional_master.xlsx')
for (modality_idx in 1:length(modalities)){
  modality = modalities[modality_idx]
  # Get the appropriate meta-data
  if (modality == "structural"){
    metadata = calm_metadata_sc
  } else{
    metadata = calm_metadata_fc
  }
  # Rename column ...22 to conners_peer_relations_t
  colnames(metadata)[which(names(metadata) == "...22")] <- "conners_peer_relations_t"
  # Recode time points!
  metadata$timepoint = recode_factor(metadata$timepoint, `0` = "baseline", `1` = "followup")
  # Initialize a list to hold the DME outputs for each timepoint
  dme_df_list_modality = vector("list",length(calm_timepoints))
  for (timepoint_idx in 1:length(calm_timepoints)){
    timepoint = calm_timepoints[timepoint_idx]
    # Subset meta-data by timepoint. 0 is baseline, 1 is followup.
    metadata_subset = metadata[which(metadata$timepoint == timepoint), ]
    # Make the BIDS column lower-case
    colnames(metadata_subset)[which(names(metadata_subset) == "BIDS")] <- "bids"
    # Find the participants with DME data
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
    # Merge variance_explained_df with nodal_df, 
    variance_and_nodal_df = merge(variance_explained_df, nodal_df)
    # Attach associated meta-data and assign to output list
    dme_df_list_modality[[timepoint_idx]] = merge(variance_and_nodal_df, metadata_subset, by = "bids", all.x = TRUE)
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
# Add a categorical factor for sex
calm_dme_and_metadata_structural$sex = factor(ifelse(calm_dme_and_metadata_structural$Sex == 0, 'M', 'F'))
calm_dme_and_metadata_functional$sex = factor(ifelse(calm_dme_and_metadata_functional$Sex == 0, 'M', 'F'))
# And rename the meanFWD column to meanfwd
names(calm_dme_and_metadata_structural)[names(calm_dme_and_metadata_structural)=="meanFWD"] <- "meanfwd"
names(calm_dme_and_metadata_functional)[names(calm_dme_and_metadata_functional)=="meanFWD"] <- "meanfwd"
### PART 3 - Load NKI Meta-Data ####
# Load DME outputs for each time point
nki_bas1_dme = read_mat('data/nki/dme/schaefer200x7_bas1.mat')
nki_flu1_dme = read_mat('data/nki/dme/schaefer200x7_flu1.mat')
nki_flu2_dme = read_mat('data/nki/dme/schaefer200x7_flu2.mat')
# Load motion estimates (mean frame-wise displacement/FWD). For each session, 
# create a new data frame with subject, session, and mean FWD for SC and FC.
# motion_estimates_v2 was created by concatenating FWD data from NKI connectomes
# on 26th March 2024 from thresholded_structural_and_functional_connectomes.mat
nki_motion_estimates = read_mat('data/nki/motion_estimates_v2.mat')[['fd.metadata']][[1]]
nki_ses_names = c("bas1","flu1","flu2")
# For each modality, create a new data frame with id, session, and mean FWD.
for (modality_idx in 1:length(modalities)){
  nki_modality_motion_estimates = nki_motion_estimates[[modality_idx]]
  motion_estimate_df = data.frame(
    meanfwd = c(unlist(nki_modality_motion_estimates[[1]][[2]]), 
                unlist(nki_modality_motion_estimates[[2]][[2]]),
                unlist(nki_modality_motion_estimates[[3]][[2]])),
    id = c(unlist(nki_modality_motion_estimates[[1]][[1]]), 
           unlist(nki_modality_motion_estimates[[2]][[1]]),
           unlist(nki_modality_motion_estimates[[3]][[1]])),
    session = c(rep("bas1",length(unlist(nki_modality_motion_estimates[[1]][[2]]))),
                rep("flu1",length(unlist(nki_modality_motion_estimates[[2]][[2]]))),
                rep("flu2",length(unlist(nki_modality_motion_estimates[[3]][[2]])))))
  # Correct a small mistake - for follow-up sessions, the data for meanFWD and 
  # ID are swapped. 
  if (modalities[modality_idx] == "functional"){
    idx_to_change = which(motion_estimate_df$session=="flu1"| motion_estimate_df$session=="flu2")
    new_meanfwd = as.numeric(motion_estimate_df$id[idx_to_change])
    new_id = as.character(motion_estimate_df$meanfwd[idx_to_change])
    motion_estimate_df$meanfwd[idx_to_change] = new_meanfwd
    motion_estimate_df$id[idx_to_change] = new_id
  }
  # Name appropriately!
  if (modalities[modality_idx] == "structural"){
    nki_motion_estimates_sc = motion_estimate_df
  } else{
    nki_motion_estimates_fc = motion_estimate_df
  }
  rm(motion_estimate_df)
}
# Load the participants data sheet (minimal meta-data). For data privacy 
# reasons, we cannot provide this as open-access.
nkir_neuroimaging_participants = 
  read_excel("//cbsu/data/imaging/projects/external/nkir/data/enhanced_nkir/data-2023-04-09T22_13_49.162Z.xlsx",na = c(".", "MD"))
# Retain ID, session, age, and gender from the neuro-imaging data sheet
nkir_neuroimaging_participants = nkir_neuroimaging_participants %>%
  dplyr::select(c("id","session","dem_002","age_04")) %>%
  distinct() %>%
  # Add a 'sub-' prefix for subjects for easier merging with other data frames
  mutate(id = paste0('sub-',id)) %>%
  # Convert categorical variables to factors
  mutate_at(c("id","session","dem_002"),as.factor) %>%
  # And convert the session variables to lower case
  mutate(session = factor(tolower(session))) %>%
  # Rename dem_002 to sex, in line with other data frames
  rename("sex"="dem_002", "scan_age" = "age_04") %>%
  # Convert sex to M and F
  mutate(sex = factor(sex, labels = c("M", "F")))


### PART 4 - Format NKI Meta-Data and DME Outputs ####
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
      subset(session == nki_session) %>% merge(motion_estimates_df[which(motion_estimates_df$session == nki_session),],all.y = T, by = c("id", "session")) %>%
      # Now merge with variance_and_nodal_dme
      merge(variance_and_nodal_dme, by = "id", all.y = T) %>%
      # And convert session to a factor
      mutate(session = factor(session))
      # And change 'age' to 'scan_age', as we retrieved this from the neuro-
      # imaging data sheet.
      # rename(scan_age = age)
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
# Remove any rows with missing values
nki_dme_and_metadata_structural = na.omit(nki_dme_and_metadata_structural)
nki_dme_and_metadata_functional = na.omit(nki_dme_and_metadata_functional)

### PART 5 - Saving Formatted Meta-Data and DME Outputs ####
calm_dme_and_metadata_structural$id = as.character(calm_dme_and_metadata_structural$ID)
names(calm_dme_and_metadata_structural) = gsub('7Networks_*', replacement = '', x = names(calm_dme_and_metadata_structural))
write.csv(calm_dme_and_metadata_structural, file = "data/calm/dme/dme.and.metadata.structural.connectivity.csv",
          quote = FALSE)
calm_dme_and_metadata_functional$id = as.character(calm_dme_and_metadata_functional$ID)
names(calm_dme_and_metadata_functional) = gsub('7Networks_*', replacement = '', x = names(calm_dme_and_metadata_functional))
write.csv(calm_dme_and_metadata_functional, file = "data/calm/dme/dme.and.metadata.functional.connectivity.csv",
          quote = FALSE)
names(nki_dme_and_metadata_structural) = gsub('7Networks_*', replacement = '', x = names(nki_dme_and_metadata_structural))
write.csv(nki_dme_and_metadata_structural, file = "data/nki/dme/dme.and.metadata.structural.connectivity.csv",
          quote = FALSE)
names(nki_dme_and_metadata_functional) = gsub('7Networks_*', replacement = '', x = names(nki_dme_and_metadata_functional))
write.csv(nki_dme_and_metadata_functional, file = "data/nki/dme/dme.and.metadata.functional.connectivity.csv",
          quote = FALSE)


### PART 6 - Descriptive Statistics for Paper ####
### PART 6A - Main Analysis ####
### NKI ####
for (modality in modalities){
  graph_theory_metrics = read.csv(sprintf('data/nki/connectomes/participant.%s.graph.theory.metrics.csv',modality))
  metadata = get(paste0('nki_dme_and_metadata_',modality))
  names(graph_theory_metrics)[names(graph_theory_metrics) == 'Timepoint'] <- 'session'
  names(graph_theory_metrics)[names(graph_theory_metrics) == 'subid'] <- 'id'
  merged_df = merge(graph_theory_metrics, metadata, by = c('session', 'id'), all.x = TRUE)
  # Drop missing values (if any)
  merged_df = na.omit(merged_df)
  merged_df$session = factor(merged_df$session)
  age_ranges_by_session = tapply(merged_df$scan_age, merged_df$session, range)
  age_mean_by_session = tapply(merged_df$scan_age, merged_df$session, mean)
  print(sprintf('For %s analyses in NKI, we examined %.1f children across %.1f timepoints, 
        spanning %.2f to %.2f years old [%.2f percent male]. %.1f scans were at baseline [Mean age of %.2f, range of %.2f and %.2f],
        %.1f second-scans [Mean age of %.2f, range of %.2f and %.2f years],
        and %.1f with three scans [Mean age of %.2f, range of %.2f and %.2f years].', modality, nrow(merged_df), 
          nlevels(merged_df$session), min(merged_df$scan_age), max(merged_df$scan_age),
          (summary(merged_df$gender)[['M']]/nrow(merged_df))*100,
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

