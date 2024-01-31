# This script pulls cognitive and behavioural (Conners) data from CALM and NKI, across all time points.
import os
import pandas as pd
import numpy as np

os.chdir('/imaging/astle/am10/diffusion.map.embedding/data/')
# Load the cognitive scales for CALM
calm_baseline = pd.read_csv('/imaging/projects/external/nkir/data/calm_cognition_behavioural/project78_calm800.csv')
# And the Conners scale
calm_baseline_conners = \
    pd.read_excel('/imaging/projects/external/nkir/data/calm_cognition_behavioural/calm800_conners_updated.xlsx')
# And for follow-up
calm_followup = pd.read_csv('/imaging/projects/external/nkir/data/calm_cognition_behavioural/'
                            'project78_calmII_masterVesrsion27102023.csv')
# Specify the columns we want. Note that the naming criteria differs for baseline and follow-up in CALM.
calm_cognitive_baseline_columns = ['id', 'age_in_months', 'awma_digit_recall_raw', 'awma_backward_digit_raw',
                                   'trails_visual_scanning_raw', 'trails_number_letter_switching_raw',
                                   'trails_motor_speed_raw', 'tower_total_achievement_raw']
calm_conners_columns = \
    ['conners_hyperactivity_impulsivity_t', 'conners_learning_problems_t', 'conners_executive_function_t',
     'conners_aggression_t', 'conners_peer_relations_t', 'conners_inattention_t', 'id']
calm_followup.rename(columns={
    'ID No.': 'id', 'Age_in_months_T2': 'age_in_months', 'Matrix_Reasoning_Raw_T2': 'matrix_reasoning_raw',
    'AWMA_Digit_Recall_Raw_T2': 'awma_digit_recall_raw', 'AWMA_Backward_Digit__Raw_T2': 'awma_backward_digit_raw',
    'TRAILS_Visual_Scanning_raw_T2': 'trails_visual_scanning_raw', 'TRAILS_Motor_Speed_Raw_T2':
        'trails_motor_speed_raw', 'TRAILS_Number_Letter_Switching_Raw_T2': 'trails_number_letter_switching_raw',
    'Tower_total_achievement_raw_T2': 'tower_total_achievement_raw', 'Conners_hyperactivity_impulsivity_T_T2':
        'conners_hyperactivity_impulsivity_t', 'Conners_learning_problems_T_T2': 'conners_learning_problems_t',
    'Conners_ExecutiveFunction_T_T2': 'conners_executive_function_t', 'Conners_aggression_T_T2': 'conners_aggression_t',
    'Conners_PeerRelations_T_T2': 'conners_peer_relations_t', 'Conners_Inattention_T_T2': 'conners_inattention_t'},
    inplace=True)
# Add a column to each CALM data frame to indicate time point
calm_baseline['timepoint'] = 'baseline'
calm_baseline_conners['timepoint'] = 'baseline'
calm_followup['timepoint'] = 'followup'
# Add time point to the column list
calm_cognitive_baseline_columns.append('timepoint')
calm_conners_columns.append('timepoint')
# Combine the features from CALM baseline and followup into a single data frame
calm_cognitive = pd.concat([calm_baseline.loc[:, calm_cognitive_baseline_columns],
                            calm_followup.loc[:, calm_cognitive_baseline_columns]]).drop_duplicates()
calm_conners = pd.concat([calm_baseline_conners.loc[:, calm_conners_columns],
                          calm_followup.loc[:, calm_conners_columns]]).drop_duplicates()
calm_cognitive_conners = calm_cognitive.merge(calm_conners, on=['id', 'timepoint'])
# Load the meta-data for CALM for structural and functional connectivity, and select relevant columns
calm_metadata_sc = pd.read_csv('calm/dme/dme.and.metadata.structural.connectivity.csv').loc[:, ['id', 'timepoint']]
calm_metadata_fc = pd.read_csv('calm/dme/dme.and.metadata.functional.connectivity.csv').loc[:, ['id', 'timepoint']]
# Merge the two to find common entries
calm_metadata_concat = pd.merge(calm_metadata_sc, calm_metadata_fc, how='inner', on=['id', 'timepoint'])
# Merge cognitive and Conners data in CALM with the meta-data for participants with structural and functional gradients.
calm_cognitive_conners_metadata = calm_metadata_concat.merge(calm_cognitive_conners, on=['id', 'timepoint'], how='left')
calm_cognitive_conners_metadata.to_csv('calm/calm_cognitive_conners_metadata.csv')
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
nki_trails = nki_trails[nki_trails['Sub Study Label'] == 'Long_child']
# Recode VA and V1 as BAS1, V4 as FLU1, and V5 as FLU2.
timepoint_mapping = {'VA': 'BAS1', 'V1': 'BAS1', 'V2': 'BAS1', 'V4': 'FLU1', 'V5': 'FLU2'}
nki_trails = nki_trails.assign(session=nki_trails["Visit"].map(timepoint_mapping))
# Rename the Anonymized ID variable in trails for easier merging with the rest of the cognitive data
nki_trails = nki_trails.rename(columns={'Anonymized ID': 'id'})
# Add the 'sub-' prefix to ID and rename session to lower-case variables
nki_trails['id'] = 'sub-' + nki_trails['id'].astype(str)
nki_trails['session'] = nki_trails['session'].str.lower()
# Remove rows with 2 or more missing values
nki_trails = nki_trails[nki_trails.isnull().sum(axis=1) < 2]
# Merge with the concatenated meta-data
nki_metadata_concat = nki_metadata_concat.merge(nki_trails.drop_duplicates(subset=['id', 'session']), how='left')
missing_count = nki_metadata_concat.loc[:, ID_columns].isnull().sum()
# At this point, we still have 347 data points, with no missing demographic data.
nki_metadata_concat.to_csv('nki/nki_metadata_cognitive.csv')
# For NKI, load the two Conners data sheets. Remove the first row, and bind together.
nki_conners_3ps = pd.read_csv(
    '/imaging/projects/external/nkir/data/enhanced_nkir/assessment_data/8100_Conners_3-P(S)_20230530.csv').loc[1:, ]
nki_conners_3srs = pd.read_csv(
    '/imaging/projects/external/nkir/data/enhanced_nkir/assessment_data/8100_Conners_3-SR(S)_20230530.csv').loc[1:, ]
"""
Wednesday 31st January 2024: Noticed a small mistake in the next line. The columns for the parent short-form (3ps) and
youth self-report (3srs) for the Conners scale are not the same. Therefore, there's no point in concatenating them!
3 lines down, we select T-scores for the Conners scale. This corresponds to the parent-report measures, just like CALM.
"""
nki_conners = pd.concat([nki_conners_3ps, nki_conners_3srs])
# Subset by the columns that we want!
nki_conners = nki_conners.loc[:, ["Anonymized ID", "Visit", "INATTENTION T-SCORE", "HYPERACTIVITY/IMPULSIVITY T-SCORE",
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
nki_conners = nki_conners[nki_conners['Sub Study Label'] == 'Long_child']
# Merge with the cognitive and meta-data
nki_metadata_concat = nki_metadata_concat.merge(nki_conners.drop_duplicates(subset=['id', 'session']), how='left')
missing_count = nki_metadata_concat.loc[:, ID_columns].isnull().sum()
# We still have 347 data points, with no missing demographic data!
nki_metadata_concat.to_csv('nki/nki_cognitive_conners_complete_with_metadata.csv')
# Add a data set variable to both the CALM and NKI data frames
nki_metadata_concat['dataset'] = 'nki'
calm_cognitive_conners_metadata['dataset'] = 'calm'
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
# Merge the two data sets together and save
nki_calm_combined_cognitive_conners = pd.concat([nki_cognitive_conners_complete_with_metadata,
                                                 calm_cognitive_conners_metadata])
nki_calm_combined_cognitive_conners.to_csv('phenotypic/nki.calm.combined.cognitive.conners.csv')
