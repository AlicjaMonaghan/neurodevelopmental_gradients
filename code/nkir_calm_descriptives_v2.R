# This script details imputation and deriving summary statistics for two data 
# sets: the center for attention, learning and memory (CALM), and the Nathan
# Rockland Kline Institute Longitudinal Study of Brain Discovery (NKI). CALM has
# two time points, and NKI has three. 

#### PART 1 - Set up the work space ####
rm(list = ls())
library(tidyverse)
library(readxl)
library(mice)
library(R.matlab)
library(rhdf5)
setwd('//cbsu/data/imaging/projects/external/nkir/')
# Load the CALM data sheets for all time points
baseline_calm_800_datasheet = read_xlsx('data/baseline_calm/Project 72_CALM 800datasheet_160421.xlsx')
baseline_calm_200_datasheet = read_xlsx('data/baseline_calm/Project 72_CALM 200datasheet_171221.xlsx')
age_at_scan_calm_200 = read_xlsx('data/baseline_calm/age_at_scan_calm_200.xlsx')
missing_metadata_da = read_xlsx('data/baseline_calm/Meta_data_duncan.xlsx')
missing_calm_bids_da = read_xlsx('data/baseline_calm/baseline_bids_missing.xlsx')
baseline_bids_missing_metadata = read_xlsx('data/baseline_calm/located_and_converted_bids_ids_DA.xlsx')
longitudinal_calm_datasheet = read_xlsx('data/longitudinal_calm/Project 78_CALM II_datasheet27072021.xlsx')
calm_mri_ids = read_xlsx('data/longitudinal_calm/Updated_CALM_MRI_IDs.xlsx')
unprocessed_bids = read_xlsx('data/baseline_calm/located_and_converted_bids_ids.xlsx')
calm_ethnicity = read.csv('/imaging/astle/users/jj02/CALM/info/ethnicityinfo_080222.csv')
calm_800_200_processed_log = read_xlsx('data/baseline_calm/mri-calm-800-200-processed-log.xlsx')
# And the same for NKI...
nkir_demos = read.csv('data/enhanced_nkir/assessment_data/8100_Demos_20230530.csv')

### PART 2 - Format CALM Neuroimaging Data Sheets ####
### Format the baseline BIDS which were previously unprocessed! ####
# Find the meta-data for those who were previously not identified
baseline_bids_missing_metadata = baseline_bids_missing_metadata[,-c(20,21)]
# Move 'age_in_months' (scanning) to the Age-Scan1 column
baseline_bids_missing_metadata$...26[which(baseline_bids_missing_metadata$...26=="AGE at scan (Months)")] <- NA
baseline_bids_missing_metadata$`Age-Scan1`[which(!is.na(baseline_bids_missing_metadata$...26))] <- 
  baseline_bids_missing_metadata$...26[which(!is.na(baseline_bids_missing_metadata$...26))]
# Subset by participants who have scans
baseline_bids_missing_metadata = subset(baseline_bids_missing_metadata,`...22`=="Y")
# Remove the 5 participants who have repeat scans 
baseline_bids_missing_metadata = baseline_bids_missing_metadata[-c(which(baseline_bids_missing_metadata$...27=="REPEAT SCAN")),]
# Remove the redundant columns
baseline_bids_missing_metadata = baseline_bids_missing_metadata[,c(1:18)]
# Format column names
baseline_bids_missing_metadata = 
  baseline_bids_missing_metadata %>% rename_with(~str_remove(.x, "\\.y$")) %>%
  select(-c("MRI-ID-T1.x")) %>%
  rename_with(~str_remove(.x, "\\.x$")) 
# Convert age columns to numeric 
baseline_bids_missing_metadata[c("Age-T1","Age-Scan1")] = sapply(baseline_bids_missing_metadata[c("Age-T1","Age-Scan1")],as.numeric)
# Replace all NA strings with missing values
for (col_idx in 1:length(colnames(baseline_bids_missing_metadata))){
  baseline_bids_missing_metadata[,col_idx][which(baseline_bids_missing_metadata[,col_idx]=="NA"),] <- NA
}
# Remove white space in the CBU IDs
baseline_bids_missing_metadata$`MRI-ID-T1` = gsub(" ","",baseline_bids_missing_metadata$`MRI-ID-T1`)
# Find the BIDS IDs of these participants
baseline_bids_missing_metadata = 
  merge(baseline_bids_missing_metadata,unprocessed_bids,by=c("MRI-ID-T1")) %>%
  rename_with(~str_remove(.x, "\\.y$")) %>%
  select(-c("BIDS-ID-T1.x"))
### Merge the newly BIDS-ed baseline participants with the larger data sheet ####
# First, ensure that the column classes are the same
numeric_cols = c("Age-T2","Time-Between-Tests","Scan-T1","Scan-T2","Time-Between-Scans","Age-Scan2","T1-months-beh-to-scan","T2-months-beh-to-scan")
baseline_bids_missing_metadata[numeric_cols] <- sapply(baseline_bids_missing_metadata[numeric_cols],as.numeric)
calm_mri_ids[numeric_cols] <- sapply(calm_mri_ids[numeric_cols],as.numeric)
baseline_bids_missing_metadata$`Date-Scan-T1` = as.POSIXct(baseline_bids_missing_metadata$`Date-Scan-T1`,format="%Y-%m-%d")
baseline_bids_missing_metadata$`Date-Scan-T2` = as.POSIXct(baseline_bids_missing_metadata$`Date-Scan-T2`,format="%Y-%m-%d")
# Convert age from months to years
baseline_bids_missing_metadata$`Age-Scan1` = baseline_bids_missing_metadata$`Age-Scan1`/12
# Now merge! Start with the two IDs which are common across both data frames
calm_mri_ids = rows_patch(calm_mri_ids,baseline_bids_missing_metadata[1:2,],by="ID",unmatched="ignore")
# And add the unique IDs
calm_mri_ids = bind_rows(calm_mri_ids,baseline_bids_missing_metadata[-c(1:2),])
# Add the age at scan for the CALM 200 participants whose data was missed out
# from the baseline_calm_200_datasheet
colnames(age_at_scan_calm_200) = c("ID","MRI-ID-T1","Age-Scan1")
calm_mri_ids = merge(calm_mri_ids, age_at_scan_calm_200, by = "ID", all = T)
### Merge the missing metadata data frame with the calm_mri_ids data frame ####
# Check that all these participants are from baseline only
missing_metadata_da$ID %in% longitudinal_calm_datasheet$`ID No.`[which(!is.na(longitudinal_calm_datasheet$Age_in_months_T2))]
# Select columns from missing_metadata_da 
missing_metadata_da = 
  missing_metadata_da[c("ID","MRI-ID-T1.y","BIDS-ID-T1.y","test_age","scan_age")]
colnames(missing_metadata_da) = c("ID","MRI-ID-T1","BIDS-ID-T1","Age-T1","Age-Scan1")
# Bind rows with the calm_mri_ids data frame. Use bind_rows because we're adding
# unique IDs
calm_mri_ids = bind_rows(missing_metadata_da,calm_mri_ids)
# To remove duplicate columns, find missing values in MRI-ID-T1.x and replace 
# with corresponding values from MRI-ID-T1.y
calm_mri_ids$`MRI-ID-T1.x`[which(is.na(calm_mri_ids$`MRI-ID-T1.x`))] = calm_mri_ids$`MRI-ID-T1.y`[which(is.na(calm_mri_ids$`MRI-ID-T1.x`))]
# And remove duplicate columns 
calm_mri_ids = select(calm_mri_ids,-c("MRI-ID-T1","MRI-ID-T1.y"))

### PART 3 - Format CALM Behavioural and Cognitive Data ####
# Here, we'll first combine the two baseline cohorts (CALM-800 and CALM-200), 
colnames(baseline_calm_800_datasheet)[1] = "ID"
colnames(baseline_calm_200_datasheet)[1] = "ID"
colnames(longitudinal_calm_datasheet)[1] = "ID"
calm_ethnicity = calm_ethnicity[1:2]
colnames(calm_ethnicity)[1] = "ID"
colnames(calm_800_200_processed_log)[1:3] = c("ID","MRI-ID-T1","BIDS-ID-T1")
calm_800_200_processed_log$`CALM_200/800` = factor(calm_800_200_processed_log$`CALM_200/800`)
# Keep distinct entries
calm_800_200_processed_log = calm_800_200_processed_log %>% distinct()
# Bind both baseline groups, but ensure that the column names match
baseline_calm_200_datasheet = rename(baseline_calm_200_datasheet, "Possible_ADHD" = "Possible_ ADHD")
baseline_calm_200_datasheet = rename(baseline_calm_200_datasheet, "Under_SLT" = "Under SLT")
calm_baseline = rbind(baseline_calm_800_datasheet,baseline_calm_200_datasheet)
# And add ethnicity data. Keep all x values because not all participants have 
# ethnicity data.
calm_baseline = merge(calm_baseline, calm_ethnicity,all.x=T)
calm_baseline$Ethnicity_Category = as.factor(calm_baseline$Ethnicity_Category)
calm_ethnicity_complete = calm_ethnicity[-which(is.na(calm_ethnicity$Ethnicity_Category)),]
calm_ethnicity_complete$Ethnicity_Category = factor(calm_ethnicity_complete$Ethnicity_Category)
levels(calm_ethnicity_complete$Ethnicity_Category) = 
  c("White","Mixed/Multiple Ethnic Groups","Asian/Asian British",
    "Black/African/Carribean/Black British","Other Ethnic Group")
summary(calm_ethnicity_complete$Ethnicity_Category)
sprintf('Ethnicity data was provided for %.2f of the referred sample, %.2f of whom 
        were White, %.2f mixed, %.2f asian/asian british, and %.2f black.',
        nrow(calm_ethnicity_complete)/nrow(calm_baseline)*100,
        summary(calm_ethnicity_complete$Ethnicity_Category)[[1]]/nrow(calm_ethnicity_complete)*100,
        summary(calm_ethnicity_complete$Ethnicity_Category)[[2]]/nrow(calm_ethnicity_complete)*100,
        summary(calm_ethnicity_complete$Ethnicity_Category)[[3]]/nrow(calm_ethnicity_complete)*100,
        summary(calm_ethnicity_complete$Ethnicity_Category)[[4]]/nrow(calm_ethnicity_complete)*100)

# Merge the two MRI information data frames together. Subset by those with 
# neuroimaging data!
calm_updated_baseline_longitudinal_mri = 
  calm_mri_ids %>% full_join(calm_800_200_processed_log,by="ID",suffix=c("",".y")) %>%
  # Remove duplicate columns
  select(-ends_with(".y")) %>%
  # And remove duplicate rows
  distinct(.) %>%
  # Subset by those with neuroimaging data
  filter(., !is.na(`BIDS-ID-T1`) | !is.na(`BIDS-ID-T2`))
# To ensure that all ages are in years, not months, convert any ages greater than
# 18 in the Age-T1 column to their year equivalents
calm_updated_baseline_longitudinal_mri$`Age-T1`[which(calm_updated_baseline_longitudinal_mri$`Age-T1`>18)] = 
  calm_updated_baseline_longitudinal_mri$`Age-T1`[which(calm_updated_baseline_longitudinal_mri$`Age-T1`>18)]/12
# Find how many CALM 800 and CALM 200 children have neuroimaging data!
sprintf("%s and %s CALM-800 and CALM-200 children, respectively, have neuroimaging data.",
        sum(!is.na(calm_updated_baseline_longitudinal_mri$`BIDS-ID-T1`) & calm_updated_baseline_longitudinal_mri$`CALM_200/800`=="800"),
        sum(!is.na(calm_updated_baseline_longitudinal_mri$`BIDS-ID-T1`) & calm_updated_baseline_longitudinal_mri$`CALM_200/800`=="200"))


### PART 4A - Descriptive statistics for baseline CALM - whole sample! ####
# For the baseline and longitudinal sample, we detail the number of participants
# recruited, age range, sex distribution, and the subset with neuroimaging data.
calm_baseline$Gender = factor(calm_baseline$Gender)
levels(calm_baseline$Gender) = c("M","F")
calm_baseline$Primary_Reason = factor(calm_baseline$Primary_Reason)
levels(calm_baseline$Primary_Reason) = 
  c("Attention", "Literacy", "Maths", "Language", 
    "Poor School Progress (Literacy and Maths)", "Memory")
calm_baseline$Diagnosis = factor(calm_baseline$Diagnosis)
# Create a factor level to differentiate between CALM 800 and CALM 200
calm_baseline$'800/200' = c(rep('800',nrow(baseline_calm_800_datasheet)),rep('200',nrow(baseline_calm_200_datasheet)))
calm_baseline$`800/200` = factor(calm_baseline$`800/200`)
# Correct spelling mistakes in column names
names(calm_baseline)[names(calm_baseline) == "PhAb_Allieration_Standard_Score_For_Analysis"] <- "PhAb_Alliteration_Standard_Score_For_Analysis"
# Find which, if any, participants have all empty data! we find that 7 
# participants were removed from baseline CALM, and 1 should be removed from 
# longitudinal CALM.
calm_baseline = calm_baseline[-which(is.na(calm_baseline$Age_in_months)),]
### PART 4B - Find which baseline CALM participants have neuroimaging data ####
# Find which baseline participants have T1w anatomical data and dwi data. Then 
# add in which have T2w and BOLD fMRI data.
calm_baseline_bids_dir = '//cbsu/data/Imaging/projects/cbu/calm/CALM-1000_BIDS/'
calm_baseline_neuroimaging_pps = dir(calm_baseline_bids_dir,pattern="sub-")
# Initialise output vectors for each type of file
calm_baseline_t1w_pps = c()
calm_baseline_t2w_pps = c()
calm_baseline_functional_pps = c()
calm_baseline_dwi_pps = c()
# For each participant, find whether they have T1w, T2w, DWI and rsfMRI BOLD NIFTI files.
for (sub_idx in 1:length(calm_baseline_neuroimaging_pps)){
  sub_id = calm_baseline_neuroimaging_pps[sub_idx]
  anatomical_files = list.files(paste0(calm_baseline_bids_dir,sub_id,'/anat/'))
  functional_files = list.files(paste0(calm_baseline_bids_dir,sub_id,'/func/'))
  dwi_files = list.files(paste0(calm_baseline_bids_dir,sub_id,'/dwi/'))
  # Find matches!
  if (paste0(sub_id,"_T1w.nii.gz")%in%anatomical_files & paste0(sub_id,"_T1w.json")%in%anatomical_files){
    calm_baseline_t1w_pps = append(calm_baseline_t1w_pps,sub_id)
  } else{
    print(paste0(sub_id," does not have T1w files."))
  }
  if (paste0(sub_id,"_T2w.nii.gz")%in%anatomical_files & paste0(sub_id,"_T2w.json")%in%anatomical_files){
    calm_baseline_t2w_pps = append(calm_baseline_t2w_pps,sub_id)
  } else{
    print(paste0(sub_id," does not have T2w files."))
  }
  if (paste0(sub_id,"_task-rest_bold.nii.gz")%in%functional_files & paste0(sub_id,"_task-rest_bold.json")%in%functional_files){
    calm_baseline_functional_pps = append(calm_baseline_functional_pps,sub_id)
  } else{
    print(paste0(sub_id," does not have functional files."))
  }
  if (paste0(sub_id,"_dwi.bval")%in%dwi_files & paste0(sub_id,"_dwi.bvec")%in%dwi_files & paste0(sub_id,"_dwi.json")%in%dwi_files & paste0(sub_id,"_dwi.nii.gz")%in%dwi_files){
    calm_baseline_dwi_pps = append(calm_baseline_dwi_pps,sub_id)
  } else{
    print(paste0(sub_id," does not have dwi files."))
  }
  
}
### PART 4C - Descriptive statistics for longitudinal CALM - whole sample ####
# Find the longitudinal children
longitudinal_calm_datasheet = longitudinal_calm_datasheet[-which(is.na(longitudinal_calm_datasheet$Age_in_months_T2)),]
longitudinal_calm_datasheet$Gender = factor(longitudinal_calm_datasheet$Gender)
longitudinal_calm_datasheet$Diagnosis_T2 = factor(longitudinal_calm_datasheet$Diagnosis_T2)
### PART 4D - Find which longitudinal CALM participants have neuroimaging data ####
calm_longitudinal_bids_dir = '//cbsu/data/Imaging/projects/cbu/calm/CALM-II_BIDS/'
calm_longitudinal_neuroimaging_pps = dir(calm_longitudinal_bids_dir,pattern="sub-")
# Initialise output vectors for each type of file. Note that no T2w scans were
# collected.
calm_longitudinal_t1w_pps = c()
calm_longitudinal_functional_pps = c()
calm_longitudinal_dwi_pps = c()
# For each participant, find whether they have T1w, DWI and rsfMRI BOLD NIFTI files.
for (sub_idx in 1:length(calm_longitudinal_neuroimaging_pps)){
  sub_id = calm_baseline_neuroimaging_pps[sub_idx]
  anatomical_files = list.files(paste0(calm_longitudinal_bids_dir,sub_id,'/anat/'))
  functional_files = list.files(paste0(calm_longitudinal_bids_dir,sub_id,'/func/'))
  dwi_files = list.files(paste0(calm_longitudinal_bids_dir,sub_id,'/dwi/'))
  # Find matches!
  if (paste0(sub_id,"_T1w.nii.gz")%in%anatomical_files & paste0(sub_id,"_T1w.json")%in%anatomical_files){
    calm_longitudinal_t1w_pps = append(calm_longitudinal_t1w_pps,sub_id)
  } else{
    print(paste0(sub_id," does not have T1w files."))
  }
  if (paste0(sub_id,"_task-rest_bold.nii.gz")%in%functional_files & paste0(sub_id,"_task-rest_bold.json")%in%functional_files){
    calm_longitudinal_functional_pps = append(calm_longitudinal_functional_pps,sub_id)
  } else{
    print(paste0(sub_id," does not have functional files."))
  }
  if (paste0(sub_id,"_dwi.bval")%in%dwi_files & paste0(sub_id,"_dwi.bvec")%in%dwi_files & paste0(sub_id,"_dwi.json")%in%dwi_files & paste0(sub_id,"_dwi.nii.gz")%in%dwi_files){
    calm_longitudinal_dwi_pps = append(calm_longitudinal_dwi_pps,sub_id)
  } else{
    print(paste0(sub_id," does not have dwi files."))
  }
  
}
### PART 4E - Format baseline and longitudinal CALM meta-data for imputation ####
# Merge baseline and longitudinal CALM data frames after merging with MRI
baseline_calm_with_mri = merge(calm_baseline,calm_mri_ids,all = T,by="ID")
longitudinal_calm_with_mri = merge(longitudinal_calm_datasheet,calm_mri_ids,all=T)
baseline_longitudinal_calm = 
  baseline_calm_with_mri %>% left_join(longitudinal_calm_with_mri,by="ID",suffix = c("", ".y")) %>% 
  # This removes the duplicate columns!
  select(-ends_with(".y"))
# Add neuroimaging information i.e. whether the participant has T1w, T2w, DWI, 
# and rsfMRI files. Start with anatomical files...
baseline_longitudinal_calm$'T1w_T1' <- 
  factor(ifelse(baseline_longitudinal_calm$ID %in% str_remove(gsub("sub-","",calm_baseline_t1w_pps),"^0+"),"1","0"))
baseline_longitudinal_calm$'T2w_T1' <- 
  factor(ifelse(baseline_longitudinal_calm$ID %in% str_remove(gsub("sub-","",calm_baseline_t2w_pps),"^0+"),"1","0"))
baseline_longitudinal_calm$'T1w_T2' <- 
  factor(ifelse(baseline_longitudinal_calm$ID %in% str_remove(gsub("sub-","",calm_longitudinal_t1w_pps),"^0+"),"1","0"))
# Now DWI files...
baseline_longitudinal_calm$'DWI_T1' <- 
  factor(ifelse(baseline_longitudinal_calm$ID %in% str_remove(gsub("sub-","",calm_baseline_dwi_pps),"^0+"),"1","0"))
baseline_longitudinal_calm$'DWI_T2' <- 
  factor(ifelse(baseline_longitudinal_calm$ID %in% str_remove(gsub("sub-","",calm_longitudinal_dwi_pps),"^0+"),"1","0"))
# And rsfMRI files...
baseline_longitudinal_calm$'rsfMRI_T1' <- 
  factor(ifelse(baseline_longitudinal_calm$ID %in% str_remove(gsub("sub-","",calm_baseline_functional_pps),"^0+"),"1","0"))
baseline_longitudinal_calm$'rsfMRI_T2' <- 
  factor(ifelse(baseline_longitudinal_calm$ID %in% str_remove(gsub("sub-","",calm_longitudinal_functional_pps),"^0+"),"1","0"))
write.table(baseline_longitudinal_calm,
            "/imaging/projects/external/nkir/data/baseline_calm/baseline_longitudinal_calm_cognitive_neuroimaging.csv",
            row.names = F, sep=",")

# Select the CALM columns to retain. We'll use these to harmonise connectome data 
# before and after the scanner change in March 2021. 
myvars = c("ID","Gender",
           "Matrix_Reasoning_T_Score_for_analysis","Matrix_Reasoning_T_Score_for_analysis_T2",
           "MRI-ID-T1.x","MRI-ID-T2",
           "BIDS-ID-T1","BIDS-ID-T2",
           "Date-Scan-T1","Date-Scan-T2",
           "Age-Scan1","Age-Scan2",
           "PhAb_Alliteration_Standard_Score_For_Analysis","PhAb_Allieration_Standard_Score_For_Analysis_T2",
           "AWMA_Digit_Recall_Standard","AWMA_Digit_Recall_Standard_T2",
           "PhAB_Object_RAN_RT_Standard_Score_for_analysis","PhAB_Object_RAN_RT_Standard_Score_for_analysis_T2",
           "AWMA_Backward_Digit_Standard","AWMA_Backward_Digit__Standard_T2",
           "AWMA_Dot_Matrix_Standard", "AWMA_Dot_Matrix_Standard_T2")
calm_scanner_metadata_for_haromisation = baseline_longitudinal_calm[myvars]
# Make all column names have underscores instead of '-' for MICE imputation
colnames(calm_scanner_metadata_for_haromisation) = str_replace_all(colnames(calm_scanner_metadata_for_haromisation),"-","_")
# Subset by those with neuroimaging data
calm_scanner_metadata_for_haromisation = 
  calm_scanner_metadata_for_haromisation %>% filter(!is.na(BIDS_ID_T1) | !is.na(BIDS_ID_T2))
# And now separate by baseline and longitudinal
calm_scanner_metadata_for_haromisation_baseline = 
  calm_scanner_metadata_for_haromisation %>% filter(!is.na(BIDS_ID_T1)) %>%
  select(!ends_with("2")) %>%
  distinct()
sprintf('For the baseline CALM with %d participants,%d were missing gender, %d missing matrix reasoning t scores,
        %d missing alliteration, %d missing object counting, %d missing digit recall, %d missing dot matrix,
        %d missing backward digit span, and %d missing age at scan',
        nrow(calm_scanner_metadata_for_haromisation_baseline), sum(is.na(calm_scanner_metadata_for_haromisation_baseline$Gender)), 
        sum(is.na(calm_scanner_metadata_for_haromisation_baseline$Matr_Reasoning_T_Score_for_analysis)),
        sum(is.na(calm_scanner_metadata_for_haromisation_baseline$PhAb_Allieration_Standard_Score_For_Analysis)), 
        sum(is.na(calm_scanner_metadata_for_haromisation_baseline$PhAB_Object_RAN_RT_Standard_Score_for_analysis)), 
        sum(is.na(calm_scanner_metadata_for_haromisation_baseline$AWMA_Digit_Recall_Standard)),
        sum(is.na(calm_scanner_metadata_for_haromisation_baseline$AWMA_Dot_Matr_Standard)), 
        sum(is.na(calm_scanner_metadata_for_haromisation_baseline$AWMA_Backward_Digit_Standard)),
        sum(is.na(calm_scanner_metadata_for_haromisation_baseline$`Age_Scan1`) & 
              !is.na(calm_scanner_metadata_for_haromisation_baseline$`BIDS_ID_T1`)))

calm_scanner_metadata_for_haromisation_longitudinal = 
  calm_scanner_metadata_for_haromisation %>% filter(!is.na(BIDS_ID_T2)) %>%
  select(!ends_with("1")) %>%
  select(!ends_with("for_analysis")) %>%
  select(!ends_with("Standard"))
sprintf('For longitudinal CALM with %d participants,%d were missing gender, %d missing matrix reasoning t scores,
        %d missing alliteration, %d missing object counting, %d missing digit recall, %d missing dot matrix,
        %d missing backward digit span, and %d missing age at scan',
        nrow(calm_scanner_metadata_for_haromisation_longitudinal), sum(is.na(calm_scanner_metadata_for_haromisation_longitudinal$Gender)), 
        sum(is.na(calm_scanner_metadata_for_haromisation_longitudinal$Matr_Reasoning_T_Score_for_analysis)),
        sum(is.na(calm_scanner_metadata_for_haromisation_longitudinal$PhAb_Allieration_Standard_Score_For_Analysis)), 
        sum(is.na(calm_scanner_metadata_for_haromisation_longitudinal$PhAB_Object_RAN_RT_Standard_Score_for_analysis)), 
        sum(is.na(calm_scanner_metadata_for_haromisation_longitudinal$AWMA_Digit_Recall_Standard)),
        sum(is.na(calm_scanner_metadata_for_haromisation_longitudinal$AWMA_Dot_Matr_Standard)), 
        sum(is.na(calm_scanner_metadata_for_haromisation_longitudinal$AWMA_Backward_Digit_Standard)),
        sum(is.na(calm_scanner_metadata_for_haromisation_longitudinal$`Age_Scan1`) & 
              !is.na(calm_scanner_metadata_for_haromisation_longitudinal$`BIDS_ID_T1`)))

# Impute missing data and save complete data frames!
calm_scanner_metadata_for_haromisation_baseline_impute = mice(calm_scanner_metadata_for_haromisation_baseline[,c(2:3,7:12)],seed = 1)
calm_scanner_metadata_for_haromisation_baseline_imputed = 
  cbind(complete(calm_scanner_metadata_for_haromisation_baseline_impute,1),calm_scanner_metadata_for_haromisation_baseline[,c(1,4:6)])
write.table(calm_scanner_metadata_for_haromisation_baseline_imputed,
            "analyses/connectomes/calm/harmonisation/calm_scanner_metadata_for_haromisation_baseline_imputed.csv",
            row.names = F,sep=",")

calm_scanner_metadata_for_haromisation_longitudinal_impute = mice(calm_scanner_metadata_for_haromisation_longitudinal[,c(2:3,8:13)],seed = 1)
calm_scanner_metadata_for_haromisation_longitudinal_imputed = 
  cbind(complete(calm_scanner_metadata_for_haromisation_longitudinal_impute,1),calm_scanner_metadata_for_haromisation_longitudinal[,c(1,4:7)])
# Correct column name spellings
names(calm_scanner_metadata_for_haromisation_longitudinal_imputed)[names(calm_scanner_metadata_for_haromisation_longitudinal_imputed) == "PhAb_Allieration_Standard_Score_For_Analysis_T2"] <- "PhAb_Alliteration_Standard_Score_for_Analysis_T2"
write.table(calm_scanner_metadata_for_haromisation_longitudinal_imputed,
            "analyses/connectomes/calm/harmonisation/calm_scanner_metadata_for_haromisation_longitudinal_imputed.csv",
            row.names = F,sep=",")
### PART 4F - Summary statistics for final CALM sample ####
# Load in the imputed baseline and longitudinal meta-data
calm_final_sample_baseline_imputed = read.csv("analyses/connectomes/calm/harmonisation/calm_scanner_metadata_for_haromisation_baseline_imputed.csv")
calm_final_sample_longitudinal_imputed = read.csv("analyses/connectomes/calm/harmonisation/calm_scanner_metadata_for_haromisation_longitudinal_imputed.csv")
# Format both data frames with factors as required
calm_final_sample_baseline_imputed = calm_final_sample_baseline_imputed %>%
  mutate_if(., is.character, as.factor)
calm_final_sample_longitudinal_imputed = calm_final_sample_longitudinal_imputed %>%
  mutate_if(., is.character, as.factor)
# Find proportion of baseline participants who were referred or non-referred
calm_final_sample_baseline_imputed$'CALM_200/800' = 
  calm_updated_baseline_longitudinal_mri$`CALM_200/800`[
    calm_updated_baseline_longitudinal_mri$ID %in% calm_final_sample_baseline_imputed$ID]
# For the referred sample, find ethnicity!
calm_final_sample_baseline_imputed_referred = calm_final_sample_baseline_imputed %>%
  subset(., `CALM_200/800`=="800") %>%
  left_join(., calm_ethnicity_complete, by = "ID") %>%
  # Convert all character columns to factors
  mutate_if(., is.character, as.factor)
# Create summary statistics for these participants
sprintf('At baseline, %.1f participants had neuroimaging data available, aged between %.2f and %.2f years, 
        with mean age %.2f and standard deviation of %.2f. %.2f percent were male. %.1f children were referred, of which
        %.2f had ethnicity data, where %.2f percent were white, %.2f percent mixed, 
        %.2f percent asian, %.2f percent black, and %.2f percent other.',
        nrow(calm_final_sample_baseline_imputed), min(calm_final_sample_baseline_imputed$Age_Scan1),
        max(calm_final_sample_baseline_imputed$Age_Scan1), mean(calm_final_sample_baseline_imputed$Age_Scan1),
        sd(calm_final_sample_baseline_imputed$Age_Scan1), 
        summary(calm_final_sample_baseline_imputed$Gender)[[2]]/nrow(calm_final_sample_baseline_imputed)*100,
        summary(calm_final_sample_baseline_imputed$`CALM_200/800`)[[2]],
        sum(!is.na(calm_final_sample_baseline_imputed_referred$Ethnicity_Category))/nrow(calm_final_sample_baseline_imputed_referred)*100,
        summary(calm_final_sample_baseline_imputed_referred$Ethnicity_Category)[[1]]/nrow(calm_final_sample_baseline_imputed_referred)*100,
        summary(calm_final_sample_baseline_imputed_referred$Ethnicity_Category)[[2]]/nrow(calm_final_sample_baseline_imputed_referred)*100,
        summary(calm_final_sample_baseline_imputed_referred$Ethnicity_Category)[[3]]/nrow(calm_final_sample_baseline_imputed_referred)*100,
        summary(calm_final_sample_baseline_imputed_referred$Ethnicity_Category)[[4]]/nrow(calm_final_sample_baseline_imputed_referred)*100,
        summary(calm_final_sample_baseline_imputed_referred$Ethnicity_Category)[[5]]/nrow(calm_final_sample_baseline_imputed_referred)*100)
# Repeat the same process for longitudinal participants. 
calm_final_sample_longitudinal_imputed = calm_final_sample_longitudinal_imputed %>%
  left_join(., calm_ethnicity_complete, by = "ID") %>%
  # Add time between scans for participants with both baseline and longitudinal scans!
  left_join(., calm_mri_ids[c("ID","Time-Between-Scans")], by = "ID")
  mutate_if(., is.character, as.factor)
sprintf('At follow-up, %.1f participants had neuroimaging data available, aged between %.2f and %.2f years,
        with mean age %.2f and standard deviation %.2f, with %.2f percent male. %.1f of these 
        participants had accompanying baseline scans, who returned an average of %.2f years later. 
        All children were referred. Ethnicity data was provided for %.2f percent of the sample, of whom %.2f were 
        white, %.2f mixed, %.2f asian, %.2f black, and %.2f other.',
        nrow(calm_final_sample_longitudinal_imputed), min(calm_final_sample_longitudinal_imputed$Age_Scan2),
        max(calm_final_sample_longitudinal_imputed$Age_Scan2), mean(calm_final_sample_longitudinal_imputed$Age_Scan2),
        sd(calm_final_sample_longitudinal_imputed$Age_Scan2), summary(calm_final_sample_longitudinal_imputed$Gender)[[2]]/
          nrow(calm_final_sample_longitudinal_imputed)*100, 
        length(intersect(calm_scanner_metadata_for_haromisation_baseline_imputed$ID,
                         calm_scanner_metadata_for_haromisation_longitudinal_imputed$ID)),
        mean(calm_final_sample_longitudinal_imputed$`Time-Between-Scans`, na.rm=T),
        sum(!is.na(calm_final_sample_longitudinal_imputed$Ethnicity_Category))/nrow(calm_final_sample_longitudinal_imputed)*100,
        summary(calm_final_sample_longitudinal_imputed$Ethnicity_Category)[[1]]/nrow(calm_final_sample_longitudinal_imputed)*100,
        summary(calm_final_sample_longitudinal_imputed$Ethnicity_Category)[[2]]/nrow(calm_final_sample_longitudinal_imputed)*100,
        summary(calm_final_sample_longitudinal_imputed$Ethnicity_Category)[[3]]/nrow(calm_final_sample_longitudinal_imputed)*100,
        summary(calm_final_sample_longitudinal_imputed$Ethnicity_Category)[[4]]/nrow(calm_final_sample_longitudinal_imputed)*100,
        summary(calm_final_sample_longitudinal_imputed$Ethnicity_Category)[[5]]/nrow(calm_final_sample_longitudinal_imputed)*100)
# Find the average interval between cognitive testing and neuroimaging.
baseline_longitudinal_calm_updated = baseline_longitudinal_calm %>%
  distinct() %>%
  # Find cognitive-neuroimaging testing interval for baseline participants.
  mutate(cognitive_neuro_interval_t1 = difftime(`Date-Scan-T1`,DOT,units='days')) %>%
  # Repeat for longitudinal participants
  mutate(cognitive_neuro_interval_t2 = difftime(`Date-Scan-T2`,DOT_T2,units='days'))
sprintf('For the %.1f participants with BIDS data at baseline, they completed the neuroimaging scans 
        an average of %.2f days after cognitive testing. For the %.1f participants with BIDS data at 
        follow-up, they completed the neuroimaging scans an average of %.2f days after cognitive.',
        sum(!is.na(baseline_longitudinal_calm_updated$`BIDS-ID-T1`)), 
        mean(baseline_longitudinal_calm_updated$cognitive_neuro_interval_t1,na.rm=T),
        sum(!is.na(baseline_longitudinal_calm_updated$`BIDS-ID-T2`)),
        mean(baseline_longitudinal_calm_updated$cognitive_neuro_interval_t2,na.rm=T))



### PART 5 - Format NKI Data Frame ####
# Format column names for the demographic data sheet
colnames(nkir_demos) = nkir_demos[1,]
nkir_demos = nkir_demos[-1,]
nkir_participants = nkir_demos %>%
  # Find participants part of the longitudinal child study
  subset(., SUB_STUDY == "Long_child") %>%
  # Convert columns to numeric
  mutate_at(c("SUB_STUDY_DAY_LAG","DAY_LAG","AGE","DEM_001","DEM_002",
              "DEM_003","DEM_004","DEM_006"),as.numeric) %>%
  # Format column names
  rename(., sex = DEM_002, ethnicity = DEM_003, race = DEM_004) %>%
  # And convert categorical variables to factors 
  mutate(sex = factor(sex, labels = c("M","F")),
         ethnicity = factor(ethnicity, labels = c("Non-Hispanic Latino/Spanish","Hispanic Latino/Spanish")),
         race = factor(race, labels = c("American Indian/Native Alaskan","Asian","Black/African American",
                                        "Native Hawiian or Other Pacific Islander","White","Other Race")))

# Count number of times each participant appears in the spreadsheet, to tell us 
# how many sessions they have 
nkir_participants_sessions = as.data.frame(table(nkir_participants$ID))
colnames(nkir_participants_sessions) = c("ID","session")
nkir_participants_sessions$session = factor(nkir_participants_sessions$session)
# Add information about race
nkir_participants_session = nkir_participants %>%
  select(c("ID","AGE","sex","race")) %>%
  distinct(ID,.keep_all = T) %>%
  right_join(., nkir_participants_sessions, by = "ID")
sprintf('%.1f NKI children were recruited, aged between %.1f and %.1f years old, 
        with %.2f percent males. %.f had a single session, %f had one follow-up session,
        and %f had two follow-up sessions.',
        nrow(nkir_participants_sessions),min(nkir_participants$AGE),
        max(nkir_participants$AGE),summary(nkir_participants$sex)[[1]]/nrow(nkir_participants)*100,
        summary(nkir_participants_sessions$session)[[1]], summary(nkir_participants_sessions$session)[[2]],
        summary(nkir_participants_sessions$session)[[3]])

### PART 3B - Add information about neuro-imaging! ####
# Load the neuroimaging participants data sheet (minimal meta-data)
nkir_neuroimaging_participants = read.table("/imaging/projects/external/nkir/data/enhanced_nkir/participants.tsv",header = T)
# Remove filepath column and keep unique entries
nkir_neuroimaging_participants = 
  subset(nkir_neuroimaging_participants,select=-c(filepath)) %>% distinct(.)
# Format column names
names(nkir_neuroimaging_participants)[names(nkir_neuroimaging_participants) == "subject"] = "ID"
names(nkir_neuroimaging_participants)[names(nkir_neuroimaging_participants) == "gender"] = "sex"
names(nkir_neuroimaging_participants)[names(nkir_neuroimaging_participants) == "age"] = "AGE"
# Find which participants have which types of neuroimaging data
nkir_bids_dir = "/imaging/projects/external/nkir/data/enhanced_nkir/"
nkir_neuroimaging_pps = dir(nkir_bids_dir,pattern="sub-*")
# Initialize output lists of vectors for each data type, corresponding to session
# number i.e. ses-BAS1 = 1, ses-FLU1 = 2, ses-FLU2 = 3
nkir_t1w_pps = list(c(),c(),c())
nkir_t2w_pps = list(c(),c(),c())
nkir_rsfmri_pps = list(c(),c(),c())
nkir_dwi_pps = list(c(),c(),c())
# Loop across participants and extract matches!
for (sub_idx in 1:length(nkir_neuroimaging_pps)){
  sub_id = nkir_neuroimaging_pps[sub_idx]
  # Find at which time points this participant has neuro-imaging data for
  sessions = list.files(paste0(nkir_bids_dir,sub_id))
  for (session_number in 1:length(sessions)){
    session = sessions[session_number]
    # Set the correct index for this session, so that we can index the output
    # lists of vectors appropriately
    if (session == "ses-BAS1"){
      session_idx = 1
    } else if (session == "ses-FLU1"){
      session_idx = 2
    } else {
      session_idx = 3
    }
    # List the files in this session
    anatomical_files = list.files(paste0(nkir_bids_dir,sub_id,'/',session,'/anat/'))
    functional_files = list.files(paste0(nkir_bids_dir,sub_id,'/',session,'/func/'))
    diffusion_files = list.files(paste0(nkir_bids_dir,sub_id,'/',session,'/dwi/'))
    # Find matches and append sub_id to relevant lists
    if (paste0(sub_id,"_",session,"_T1w.nii.gz")%in%anatomical_files & paste0(sub_id,"_",session,"_T1w.json")%in%anatomical_files){
      nkir_t1w_pps[[session_idx]] = append(nkir_t1w_pps[[session_idx]],sub_id)
    } else{
      print(paste(sub_id,"in",session, "does not have T1w files."))
    }
    if (paste0(sub_id,"_",session,"_T2w.nii.gz")%in%anatomical_files & paste0(sub_id,"_",session,"_T2w.json")%in%anatomical_files){
      nkir_t2w_pps[[session_idx]] = append(nkir_t2w_pps[[session_idx]],sub_id)
    } else{
      print(paste(sub_id,"in",session, "does not have T2w files."))
    }
    if (paste0(sub_id,"_",session,"_dwi.bval")%in%diffusion_files & paste0(sub_id,"_",session,"_dwi.bvec")%in%diffusion_files & 
        paste0(sub_id,"_",session,"_dwi.json")%in%diffusion_files & paste0(sub_id,"_",session,"_dwi.nii.gz")%in%diffusion_files){
      nkir_dwi_pps[[session_idx]] = append(nkir_dwi_pps[[session_idx]],sub_id)
    } else{
      print(paste(sub_id,"in",session, "does not have dwi files."))
    }
    if (paste0(sub_id,"_",session,"_task-rest_acq-1400_bold.json")%in%functional_files & paste0(sub_id,"_",session,"_task-rest_acq-1400_bold.nii.gz")%in%functional_files){
      nkir_rsfmri_pps[[session_idx]] = append(nkir_rsfmri_pps[[session_idx]],sub_id)
    } else{
      print(paste(sub_id,"in",session, "does not have rsfmri files."))
    }
  }
}
# Format nkir_participants data frame 
nkir_participants$ID = factor(nkir_participants$ID)
# For each list, find participants who feature once, twice, or three times.
for (total_number_of_scans in 1:3){
  # Unlist the lists of participants
  t1w_participants_unlisted = data.frame("ID" = unlist(nkir_t1w_pps))
  rsfmri_participants_unlisted = data.frame("ID" = unlist(nkir_rsfmri_pps))
  dwi_participants_unlisted = data.frame("ID" = unlist(nkir_dwi_pps))
  for (modality_idx in 1:2){
    # Specify whether we're searching for participants with T1w and either DWI
    # or fMRI scans.
    if (modality_idx == 1){
      second_modality = dwi_participants_unlisted
      second_modality_name = "DWI"
    } else{
      second_modality = rsfmri_participants_unlisted
      second_modality_name = "rsfMRI"
    }
    
    # Find how many participants have X number of T1w and second modality scans
    t1w_participants_sliced = t1w_participants_unlisted %>% group_by(ID) %>% slice(total_number_of_scans)
    second_modality_participants_sliced = second_modality %>% group_by(ID) %>% slice(total_number_of_scans)
    third_modality_participants_sliced = rsfmri_participants_unlisted %>% group_by(ID) %>% slice(total_number_of_scans)
    t1w_and_second_modality_participants_sliced = intersect(t1w_participants_sliced,second_modality_participants_sliced)
    t1w_and_third_modality_participants_sliced = intersect(t1w_participants_sliced,third_modality_participants_sliced)
    # Remove sub- prefix for easier indexing into other data sheets
    t1w_and_second_modality_participants_sliced$ID = str_replace(t1w_and_second_modality_participants_sliced$ID,"sub-","")
    
    # Extract demographic information for these participants by slicing the 
    # demographic data frame according to which instance (scan number) we're 
    # examining
    sliced_demo_df = 
      nkir_participants[nkir_participants$ID %in% t1w_and_second_modality_participants_sliced$ID,] %>% 
      group_by(ID) %>% slice(total_number_of_scans)
    # If there are some participants without demographic matches in the nkir_participants data frame,
    # then use the neuroimaging data frame!
    missing_demo_participants = setdiff(t1w_and_second_modality_participants_sliced$ID, sliced_demo_df$ID)
    if (length(missing_demo_participants) > 0){
      sliced_demo_df = nkir_neuroimaging_participants[nkir_neuroimaging_participants$ID %in% missing_demo_participants,] %>%
        group_by(ID) %>% slice(total_number_of_scans) %>% bind_rows(., sliced_demo_df) %>% mutate_if(., is.character, as.factor)
      }
    summary_statistics_to_user = sprintf('%.1f participants had %.1f T1w and %s scan(s). %.2f percent were male, with
            mean and sd ages of %.2f and %.2f respectively, with range of %.2f and %.2f. 
            %.2f participants had race information available, of which %.2f 
            percent were American Indian, %.2f percent Asian, %.2f percent Black, 
            %.2f percent Pacific Islander, %.2f percent White, and %.2f Other. ',
            nrow(t1w_and_second_modality_participants_sliced), total_number_of_scans, second_modality_name,
            summary(sliced_demo_df$sex)[[1]]/nrow(sliced_demo_df)*100,
            mean(sliced_demo_df$AGE), sd(sliced_demo_df$AGE), min(sliced_demo_df$AGE),
            max(sliced_demo_df$AGE), sum(!is.na(sliced_demo_df$race)), 
            summary(sliced_demo_df$race)[[1]]/nrow(sliced_demo_df)*100,
            summary(sliced_demo_df$race)[[2]]/nrow(sliced_demo_df)*100,
            summary(sliced_demo_df$race)[[3]]/nrow(sliced_demo_df)*100,
            summary(sliced_demo_df$race)[[4]]/nrow(sliced_demo_df)*100,
            summary(sliced_demo_df$race)[[5]]/nrow(sliced_demo_df)*100,
            summary(sliced_demo_df$race)[[6]]/nrow(sliced_demo_df)*100)
    print(summary_statistics_to_user)
  }
}
# For each list, find participants who feature once, twice, or three times.
for (total_number_of_scans in 1:3){
  # Unlist the lists of participants
  t1w_participants_unlisted = data.frame("ID" = unlist(nkir_t1w_pps))
  rsfmri_participants_unlisted = data.frame("ID" = unlist(nkir_rsfmri_pps))
  dwi_participants_unlisted = data.frame("ID" = unlist(nkir_dwi_pps))
  # Find how many participants have X number of the different scan types
  t1w_participants_sliced = t1w_participants_unlisted %>% group_by(ID) %>% slice(total_number_of_scans)
  rsfmri_participants_sliced = rsfmri_participants_unlisted %>% group_by(ID) %>% slice(total_number_of_scans)
  dwi_participants_sliced = dwi_participants_unlisted %>% group_by(ID) %>% slice(total_number_of_scans)
  # Find participants who have T1w scans and either DWI or rsfMRI scans
  t1w_and_rsfmri_sliced = intersect(t1w_participants_sliced,rsfmri_participants_sliced)
  t1w_and_dwi_sliced = intersect(t1w_participants_sliced,dwi_participants_sliced)
  # Take the bigger of the two lists
  if (length(t1w_and_dwi_sliced$ID) > length(t1w_and_rsfmri_sliced$ID)){
    more_common_modality = t1w_and_dwi_sliced$ID
  } else {
    more_common_modality = t1w_and_rsfmri_sliced$ID
  }
  # Remove sub- prefix for easier indexing into other data sheets
  more_common_modality = str_replace(more_common_modality,"sub-","")
  # Extract demographic information for these participants by slicing the 
  # demographic data frame according to which instance (scan number) we're 
  # examining
  sliced_demo_df = 
    nkir_participants[nkir_participants$ID %in% more_common_modality,] %>% 
    group_by(ID) %>% slice(total_number_of_scans)
  # If there are some participants without demographic matches in the nkir_participants data frame,
  # then use the neuroimaging data frame!
  missing_demo_participants = setdiff(more_common_modality, sliced_demo_df$ID)
  if (length(missing_demo_participants) > 0){
    sliced_demo_df = nkir_neuroimaging_participants[nkir_neuroimaging_participants$ID %in% missing_demo_participants,] %>%
      group_by(ID) %>% slice(total_number_of_scans) %>% bind_rows(., sliced_demo_df) %>% mutate_if(., is.character, as.factor)
  }
  summary_statistics_to_user = sprintf(
  '%.1f participants had %.1f T1w and at least one of rsfMRI or DWI scan(s). 
  %.2f percent were male, with mean and sd ages of %.2f and %.2f respectively, 
 with range of %.2f and %.2f. %.2f percent had race information available, 
 of which %.2f percent were American Indian, %.2f percent Asian, %.2f percent Black, 
 %.2f percent Pacific Islander, %.2f percent White, and %.2f Other.',
   length(more_common_modality), total_number_of_scans,
   summary(sliced_demo_df$sex)[[1]]/nrow(sliced_demo_df)*100,
   mean(sliced_demo_df$AGE), sd(sliced_demo_df$AGE), min(sliced_demo_df$AGE),
   max(sliced_demo_df$AGE), sum(!is.na(sliced_demo_df$race))/nrow(sliced_demo_df)*100, 
   summary(sliced_demo_df$race)[[1]]/nrow(sliced_demo_df)*100,
   summary(sliced_demo_df$race)[[2]]/nrow(sliced_demo_df)*100,
   summary(sliced_demo_df$race)[[3]]/nrow(sliced_demo_df)*100,
   summary(sliced_demo_df$race)[[4]]/nrow(sliced_demo_df)*100,
   summary(sliced_demo_df$race)[[5]]/nrow(sliced_demo_df)*100,
   summary(sliced_demo_df$race)[[6]]/nrow(sliced_demo_df)*100)
  print(summary_statistics_to_user)
}
