# This script pulls cognitive and behavioural (Conners) data from CALM and NKI, across all time points.
import os
import pandas as pd
import numpy as np

os.chdir('/home/am10/gradients_open_access/')
# Load the DME results for both modalities
calm_dme_metadata_sc = pd.read_csv('data/calm/dme/dme.and.metadata.structural.connectivity.csv')
calm_dme_metadata_fc = pd.read_csv('data/calm/dme/dme.and.metadata.functional.connectivity.csv')
# Merge...


# Subset by the columns we want!
calm_columns_to_subset = ['BIDS', 'scan_age', 'Sex'. 'timepoint']
# Concatenate
calm_metadata_concat = pd.merge(calm_struc_master, calm_func_master, how='inner', on=['BIDS', 'timepoint', 'Sex', 'scan_age'])
# Load the meta- and DME-data for structural and functional connectivity in NKI
nki_dme_metadata_sc = pd.read_csv(
    '/imaging/astle/am10/diffusion.map.embedding/data/nki/dme/dme.and.metadata.structural.connectivity.csv')
nki_dme_metadata_fc = pd.read_csv(
    '/imaging/astle/am10/diffusion.map.embedding/data/nki/dme/dme.and.metadata.functional.connectivity.csv')
# Merge the two data sheets to find common entries
nki_metadata_concat = pd.merge(nki_dme_metadata_sc, nki_dme_metadata_fc, how='inner',
                               on=['id', 'session', 'sex', 'scan_age'])
# Subset by the identification columns
ID_columns = ['id', 'session', 'sex', 'scan_age']
nki_metadata_concat = nki_metadata_concat.loc[:, ID_columns]
# At this point, we have 347 participants in NKI, with no missing values.
# Load the cognitive data sheet, and select the columns we need i.e. raw forward and backward digit span, and the D-KEFS
# Tower sub-scale.
nki_cognitive_loris = pd.read_excel(
    '/imaging/projects/external/nkir/data/enhanced_nkir/data-2023-04-09T22_13_49.162Z.xlsx', na_values='.').loc[:,
                      ['id', 'session', 'age_04', 'dspan_01', 'dspan_04', 'tow_46']]
# Add a 'sub-' prefix to ID and then make session lower-case
nki_cognitive_loris['id'] = 'sub-' + nki_cognitive_loris['id'].astype(str)
nki_cognitive_loris['session'] = nki_cognitive_loris['session'].str.lower()
# Merge with the concatenated meta-data
nki_metadata_concat = nki_metadata_concat.merge(nki_cognitive_loris, on=['id', 'session'], how='left')
missing_count = nki_metadata_concat.loc[:, ID_columns]
# At this point, we still have 347 participants, with no missing demographic data.
# Load the trails measure for NKI
nki_trails = pd.read_csv(
    '/imaging/projects/external/nkir/data/enhanced_nkir/assessment_data/8100_DKEFS_Trails_20230530.csv').loc[1:,
             ['Anonymized ID', 'Visit', 'Sub Study Label', 'Visual Scanning', 'Number-Letter Switching', 'Motor Speed']]
# Subset by the longitudinal child discovery study
# nki_trails = nki_trails[nki_trails['Sub Study Label'] == 'Long_child']
# Recode VA and V1 as BAS1, V4 as FLU1, and V5 as FLU2.
timepoint_mapping = {'VA': 'BAS1', 'V1': 'BAS1', 'V2': 'BAS1', 'V4': 'FLU1', 'V5': 'FLU2'}
nki_trails = nki_trails.assign(session=nki_trails["Visit"].map(timepoint_mapping))
# Rename the Anonymized ID variable in trails for easier merging with the rest of the cognitive data
nki_trails = nki_trails.rename(columns={'Anonymized ID': 'id'})
# Add the 'sub-' prefix to ID and rename session to lower-case variables
nki_trails['id'] = 'sub-' + nki_trails['id'].astype(str)
nki_trails['session'] = nki_trails['session'].str.lower()
# Remove rows with 2 or more missing values
# nki_trails = nki_trails[nki_trails.isnull().sum(axis=1) < 2]
# Merge with the concatenated meta-data
nki_metadata_concat = nki_metadata_concat.merge(nki_trails.drop_duplicates(subset=['id', 'session']), how='left')
missing_count = nki_metadata_concat.loc[:, ID_columns].isnull().sum()
# For NKI, load the two Conners data sheets. Remove the first row, and bind together.
nki_conners_3ps = pd.read_csv(
    '/imaging/projects/external/nkir/data/enhanced_nkir/assessment_data/8100_Conners_3-P(S)_20230530.csv').loc[1:, ]
# Subset by the columns that we want!
nki_conners = nki_conners_3ps.loc[:,
              ["Anonymized ID", "Visit", "INATTENTION T-SCORE", "HYPERACTIVITY/IMPULSIVITY T-SCORE",
               "LEARNING PROBLEMS T-SCORE", "EXECUTIVE FUNCTIONING T-SCORE", "AGGRESSION T-SCORE",
               "PEER RELATIONS T-SCORE", "Sub Study Label"]]
# Rename the Anonymized ID variable for easier merging with the meta-data
nki_conners = nki_conners.rename(columns={'Anonymized ID': 'id'})
# Add the 'sub-' prefix
nki_conners['id'] = 'sub-' + nki_conners['id']
# Recode timepoints
nki_conners = nki_conners.assign(session=nki_conners["Visit"].map(timepoint_mapping))
# And make session lower-case
nki_conners['session'] = nki_conners['session'].str.lower()
# Find the participants in the longitudinal child development study
# nki_conners = nki_conners[nki_conners['Sub Study Label'] == 'Long_child']
# Merge with the cognitive and meta-data
nki_metadata_concat = nki_conners.drop_duplicates(
    subset=['id', 'session']).merge(nki_metadata_concat, how='right', on=['id','session'])
missing_count = nki_metadata_concat.loc[:, ID_columns].isnull().sum()
# Add a data set variable to both the CALM and NKI data frames
nki_metadata_concat['dataset'] = 'nki'
calm_func_master['dataset'] = 'calm'
# Subset the columns we want from the NKI data set
nki_cognitive_conners_complete_with_metadata = \
    nki_metadata_concat.loc[:, ['id', 'session', 'scan_age', 'dspan_01', 'dspan_04', 'tow_46',
                                'Visual Scanning', 'Number-Letter Switching', 'Motor Speed',
                                'INATTENTION T-SCORE', 'HYPERACTIVITY/IMPULSIVITY T-SCORE',
                                'LEARNING PROBLEMS T-SCORE', 'EXECUTIVE FUNCTIONING T-SCORE',
                                'AGGRESSION T-SCORE', 'PEER RELATIONS T-SCORE', 'dataset']]
# Rename the columns so that they're the same as CALM
nki_cognitive_conners_complete_with_metadata.rename(
    columns={'session': 'timepoint', 'scan_age': 'age_in_months', 'dspan_01': 'awma_digit_recall_raw', 'dspan_04':
        'awma_backward_digit_raw', 'Visual Scanning': 'trails_visual_scanning_raw', 'Number-Letter Switching':
                 'trails_number_letter_switching_raw', 'Motor Speed': 'trails_motor_speed_raw',
             'tow_46': 'tower_total_achievement_raw', 'INATTENTION T-SCORE': 'conners_inattention_t',
             'HYPERACTIVITY/IMPULSIVITY T-SCORE': 'conners_hyperactivity_impulsivity_t', 'LEARNING PROBLEMS T-SCORE':
                 'conners_learning_problems_t', 'EXECUTIVE FUNCTIONING T-SCORE': 'conners_executive_function_t',
             'AGGRESSION T-SCORE': 'conners_aggression_t', 'PEER RELATIONS T-SCORE': 'conners_peer_relations_t'},
    inplace=True)
# For NKI, make sure that the age is in months!
nki_cognitive_conners_complete_with_metadata['age_in_months'] = \
    nki_cognitive_conners_complete_with_metadata['age_in_months']*12
# Merge the two data sets together and save
nki_calm_combined_cognitive_conners = pd.concat([nki_cognitive_conners_complete_with_metadata,
                                                 calm_cognitive_conners_metadata])
nki_calm_combined_cognitive_conners.to_csv('phenotypic/nki.calm.combined.cognitive.conners.csv')
