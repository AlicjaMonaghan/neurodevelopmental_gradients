# This python script calculates structure-function coupling for CALM and NKI. This is the Spearman correlation between
# the Euclidean distance of each node relative to all others in structural and functional manifold space. Updated in
# July-August 2025 by Alicja Monaghan, MRC Cognition and Brain Sciences Unit.

from scipy.stats import spearmanr, shapiro, pearsonr
import pandas as pd
import os
import numpy as np
from scipy.io import loadmat
import glob
import re
import statsmodels.formula.api as smf
from factor_analyzer import FactorAnalyzer
import h5py
from scipy.spatial.distance import cdist
from scipy.stats import mannwhitneyu, kstest, ttest_ind
from pingouin import compute_effsize
from sklearn.impute import KNNImputer
from sklearn.pipeline import Pipeline
from sklearn.metrics import root_mean_squared_error

os.chdir('/Users/alicjamonaghan/Desktop/neurodevelopmental_gradients')


# This function calculates structure-function coupling for the user-specified data set, as a spearman-rank correlation
# between the euclidean distance between each node and all others, for each modality.
def run_structure_function_euclidean(dataset, parcellation, return_within_modality_distance_matrices):
    individual_gradients_dir = os.getcwd() + '/data/' + dataset + '/dme/'
    # Load the labels for the parcellation scheme.
    parcellation_metadata = loadmat(glob.glob(os.getcwd() + '/data/' + parcellation + '*.mat')[0], squeeze_me=True)
    # Extract the last key i.e. that which holds the labels and meta-data
    parcellation_metadata = parcellation_metadata[list(parcellation_metadata)[-1]]
    parcellation_labels = parcellation_metadata['name'].item()
    # If we're using the Schaefer parcellations, remove the '7Networks_' suffix from each label, as this produces an X
    # prefix otherwise when exporting to a CSV.
    if 'schaefer' in parcellation:
        for idx, label in enumerate(parcellation_labels):
            parcellation_labels[idx] = np.char.split(label, sep='_', maxsplit=1).item()[1]
    # For each time point in the directory, load the individual-level structural and functional gradients. If we
    # use CALM, make sure we're selecting the referred participants.
    if dataset == 'calm':
        individual_gradients_dir = individual_gradients_dir + 'referred_'
    # Create a list to hold the output data frames
    output_df_list = []
    within_manifold_distance_matrices_list = []
    for file in glob.glob(individual_gradients_dir + parcellation + '_*.mat'):
        individual_gradients = loadmat(file)
        # Extract the time point
        timepoint = file.partition('schaefer200x7_')[2].split('.mat')[0]
        # Find which participants have both structural and functional gradients, and return the indices.
        idx_both_modalities = \
            np.stack([(xi, xp) for (xi, x) in enumerate(list(individual_gradients['structural_sub_list']))
                      for (xp, y) in enumerate(list(individual_gradients['functional_sub_list'])) if x == y])
        print('For {} {}, {} participants had structural gradients, {} functional, and {} with both.'.
              format(re.search('schaefer200x7_(.*).mat', file).group(1), dataset,
                     len(individual_gradients['structural_sub_list']),
                     len(individual_gradients['functional_sub_list']), len(idx_both_modalities)))
        # Create a new array with corresponding structural and functional gradients for the relevant participants,
        # and append to output array
        individual_gradients_array = np.stack(
            (individual_gradients['structural_eigenvectors'][:, :, idx_both_modalities[:, 0]],
             individual_gradients['functional_eigenvectors'][:, :, idx_both_modalities[:, 1]]))
        # Within each modality, calculate the euclidean distance between each node and all others
        nroi, nsub = individual_gradients_array.shape[1], individual_gradients_array.shape[3]
        within_manifold_nodal_euclidean = np.zeros((nsub, nroi, nroi, 2))
        between_manifold_nodal_euclidean = np.zeros((nsub, nroi))
        for sub_idx in range(0, nsub):
            for modality_idx in range(0, 2):
                within_manifold_nodal_euclidean[sub_idx, :, :, modality_idx] = cdist(
                    individual_gradients_array[modality_idx, :, :, sub_idx],
                    individual_gradients_array[modality_idx, :, :, sub_idx])
                # Append to the output list
                within_manifold_distance_matrices_list.append(within_manifold_nodal_euclidean)
        # After we've calculated the relative position of each node within its manifold, calculate the spearman-rank
        # correlation coefficient between each set of euclidean distances for each node.
        for sub_idx in range(0, nsub):
            for roi in range(0, nroi):
                between_manifold_nodal_euclidean[sub_idx, roi] = \
                    spearmanr(within_manifold_nodal_euclidean[sub_idx, roi, :, 0],
                              within_manifold_nodal_euclidean[sub_idx, roi, :, 1])[0]
        # After looping through all participants and modalities for this time point, create a summary data frame
        between_manifold_nodal_euclidean_df = pd.DataFrame(between_manifold_nodal_euclidean)
        # Name the columns with the parcellated regions
        between_manifold_nodal_euclidean_df = between_manifold_nodal_euclidean_df.set_axis([parcellation_labels],
                                                                                           axis=1)
        # Add subject ID
        if dataset == 'nki':
            # For NKI, SUB IDs are universal across time points, meaning that sub-A0001 at baseline is the same as sub-
            # A0001 at the first and second follow-up. This is not the case for CALM - there is one integer ID to find
            # the same participant across time points, but this is almost always different from their sub ID.
            between_manifold_nodal_euclidean_df['id'] = list(
                individual_gradients['structural_sub_list'][idx_both_modalities[:, 0]])
        else:
            # Therefore, for CALM, load the structural meta-data.
            calm_structural_metadata = pd.read_csv('data/calm/dme/dme.and.metadata.structural.connectivity.csv')
            # Subset meta-data by this timepoint
            calm_structural_metadata = calm_structural_metadata[calm_structural_metadata['timepoint'] == timepoint]
            # Find the universal IDs associated with the BIDS (sub-) identifiers
            BIDS = pd.DataFrame({'bids': individual_gradients['structural_sub_list'][idx_both_modalities[:, 0]]})
            merged_bids_id_df = BIDS.merge(calm_structural_metadata.loc[:, ['bids', 'id']], on='bids')
            between_manifold_nodal_euclidean_df['id'] = merged_bids_id_df['id'].values
        # Append to the output list
        output_df_list.append(between_manifold_nodal_euclidean_df)
    if return_within_modality_distance_matrices is True:
        return output_df_list, within_manifold_distance_matrices_list
    else:
        return output_df_list


# This function imputes missing cognitive/behavioural data, and then conducts a principal component analysis. In the
# case of raw cognitive scores, we regress out age to standardise them. This produces two outputs - one data frame
# showing the loadings of each measure onto each component, and a second data frame showing the loadings of each
# participant onto each component.
def principal_component_analysis(domain, ncomp):
    # Load the data sheet with cognitive and psychopathology measures across cohorts.
    cognitive_psychopathology_data_sheet = pd.read_csv(
        'data/phenotypic/august25.updated.nki.calm.combined.cognitive.conners.csv')
    # We will always include the following columns
    id_columns = ['id', 'age_in_months', 'timepoint', 'dataset']
    # 'domain' takes one of two options: cognitive or psychopathology.
    if domain == 'cognitive':
        measure_columns = ['awma_digit_recall_raw', 'awma_backward_digit_raw', 'tower_total_achievement_raw',
                           'trails_visual_scanning_raw', 'trails_number_letter_switching_raw', 'trails_motor_speed_raw']
    else:
        measure_columns = ['conners_inattention_t', 'conners_hyperactivity_impulsivity_t', 'conners_peer_relations_t',
                           'conners_learning_problems_t', 'conners_executive_function_t', 'conners_aggression_t']
    # Set the number of data points we have and initialise the imputer
    num_data_points = cognitive_psychopathology_data_sheet.shape[0]
    pipeline = Pipeline([('scaler', StandardScaler()),('knn_imputer', KNNImputer(n_neighbors=5))])
    # Loop through each measure, and calculate the percentage of missing data. If there is missing data, impute using
    # k-nearest neighbours.
    for measure in measure_columns:
        x = cognitive_psychopathology_data_sheet[measure].values.reshape(num_data_points, 1)
        print("%.2f percent of %s data is missing." % (round((np.isnan(x).sum() / num_data_points) * 100, 2), measure))
        cognitive_psychopathology_data_sheet[measure] = pipeline.fit_transform(x)
    # Save the imputed datasets!
    cognitive_psychopathology_data_sheet.to_csv('data/phenotypic/calm_nki_imputed_' + domain + '.csv', index=False)
    # If we're processing cognitive data, we need to derive age-standardised scores by regressing out age.
    if domain is 'cognitive':
        age = cognitive_psychopathology_data_sheet['age_in_months'].values.reshape(num_data_points, 1)
        for cognitive_measure in measure_columns:
            # Add a column in the data sheet for age squared
            cognitive_psychopathology_data_sheet['age_squared'] = age ** 2
            age_regression = smf.ols(formula=cognitive_measure + ' ~ age + age_squared',
                                     data=cognitive_psychopathology_data_sheet)
            res = age_regression.fit()
            # Extract the residuals and assign to the data frame!
            cognitive_psychopathology_data_sheet[cognitive_measure + '_std'] = res.resid
        # After looping through all measures, update the measure columns
        measure_columns = [i + '_std' for i in measure_columns]
    # Initialise an array to hold the Z-scored version of the scales
    z_scored_array = np.zeros(shape=(num_data_points, len(measure_columns)))
    for column_idx, column in enumerate(measure_columns):
        if 'std' in column:
            z_score_measure = (cognitive_psychopathology_data_sheet[column].values - 100) / 15
        elif '_t' in column:
            z_score_measure = (cognitive_psychopathology_data_sheet[column].values - 50) / 10
        if np.isnan(z_score_measure).sum() > 0:
            print("Missing values in {}.".format(column))
        # Allocate the z-scored measure to the new array
        z_scored_array[:, column_idx] = z_score_measure
    # Initialise and fit the PCA with varimax rotation
    fa = FactorAnalyzer(rotation="varimax", n_factors=ncomp, method='principal')
    fa.fit(z_scored_array)
    # Get the rotated factor pattern
    loadings = pd.DataFrame(fa.loadings_, index=measure_columns, columns=[f"Factor{i + 1}" for i in range(ncomp)])
    # Get the variance explained by each component and cumulatively
    variance_explained = fa.get_factor_variance()
    print("For the %s scales, the cumulative variance explained is %.2f per cent" % (
        domain, round(variance_explained[2][-1] * 100, 2)))
    for comp in range(0, ncomp):
        print("Component %.1f variance explained: %.2f percent." % (
            comp + 1, round(variance_explained[1][comp] * 100, 2)))
    # Create a data frame to hold the loadings for each participant, with their ID, time point, data set, and age.
    individual_loadings = pd.DataFrame(fa.transform(z_scored_array), columns=[f"Factor{i + 1}" for i in range(ncomp)])
    individual_loadings_pd = individual_loadings.join(cognitive_psychopathology_data_sheet[id_columns])
    # Return two data frames - one with the loadings of the measures onto each factor, and a second with the loadings
    # of each participant onto each factor.
    return loadings, individual_loadings_pd

########################################################
# PART 1 - Structure-Function Relationships #
########################################################
# For each data set and participant, calculate the Euclidean distance for each node in relation to others within the
# same manifold, and then the distance between these measures in the structural and functional manifolds.
datasets = ['calm', 'nki']
parcellation = 'schaefer200x7'
for dataset in datasets:
    # This function will calculate between-manifold euclidean distances and produce a summary data frame for each time
    # point. Additionally, for NKI, collect the within-manifold euclidean distances for visualisation.
    if dataset is 'calm':
        euclidean_correlation_df = run_structure_function_euclidean(dataset, parcellation='schaefer200x7',
                                                                    return_within_modality_distance_matrices=False)
        timepoints = ['baseline', 'followup']
    else:
        euclidean_correlation_df, within_manifold_euclid = run_structure_function_euclidean(
            dataset, parcellation='schaefer200x7', return_within_modality_distance_matrices=True)
        timepoints = ['bas1', 'flu1', 'flu2']
    for timepoint_idx, timepoint in enumerate(timepoints):
        # Add time point to the data frame
        euclidean_correlation_df[timepoint_idx]['timepoint'] = timepoint
        # And save...
        euclidean_correlation_df[timepoint_idx].to_csv(
            'data/' + dataset + '/structure.function/euclidean.distance.spearman.' + timepoint + '.' +
            parcellation + '.csv', header=True, index=False)
    # Additionally, save the within-manifold euclidean distances for NKI
    if dataset is 'nki':
        for timepoint_idx, timepoint in enumerate(timepoints):
            hf = h5py.File(os.getcwd() + '/data/' + dataset + '/structure.function/within.manifold.euclidean.distance.'
                           + timepoint + '.h5', 'w')
            hf.create_dataset('within.manifold.euclidean.distances', data=within_manifold_euclid[timepoint_idx])
            hf.close()

########################################################
# PART 2 - Updating Meta-Data Sheet #
########################################################
domains = ['cognitive', 'psychopathology']
# Load the NKI and CALM psychopathology and cognitive data
nki_calm_psych_cog = pd.read_csv('data/phenotypic/nki.calm.combined.cognitive.conners.csv', index_col=False).iloc[:, 1:]
# Find the number of CALM participants we have
calm_nsub = nki_calm_psych_cog['dataset'].value_counts()['calm']
# We need to replace the ages for CALM with the correct ages. Load up the master sheet with the correct ages.
calm_func_master = pd.read_excel('data/calm/Alicja_calm_fc_master.xlsx')
# Convert scan age to age in months
# calm_func_master['updated_age_in_months'] = calm_func_master['scan_age'] * 12
# Recode time point from numbers to character labels
calm_func_master['timepoint'] = np.where(calm_func_master['timepoint'] == 0, 'baseline', 'followup')
# Replace ID with id...
calm_func_master = calm_func_master.rename(columns={'ID': 'id'})
# Add a dataset column
calm_func_master['dataset'] = 'calm'
# Convert id and timepoint to objects, for easier merging...
calm_func_master = calm_func_master.astype({'id': 'object', 'timepoint': 'object', 'dataset': 'object'})
# Subset by the columns we want
calm_func_master = calm_func_master.loc[:, ['id', 'timepoint', 'Age_at_test(years)', 'scan_age']]
calm_func_master['id'] = calm_func_master['id'].astype(str)
calm_func_master['timepoint'] = calm_func_master['timepoint'].astype(str)
# Subset the NKI and CALM psychopathology and cognitive data by CALM
calm_psych_cog = nki_calm_psych_cog[nki_calm_psych_cog['dataset'] == 'calm']
# Reset the index
calm_psych_cog.reset_index(drop=True, inplace=True)
# Merge with the correct age data!
calm_psych_cog = calm_psych_cog.merge(calm_func_master, on=['id', 'timepoint'], how='left')
# Remove old age_in_months column
# calm_psych_cog.drop(columns='age_in_months', inplace=True)
# And rename updated_age_in_months to age_in_months, in line with the original data sheet
calm_psych_cog.rename(columns={'updated_age_in_months': 'age_in_months'}, inplace=True)
# Extract measures of psychopathology and cognition for NKI
nki_psych_cog = nki_calm_psych_cog[nki_calm_psych_cog['dataset'] == 'nki']
# We have structure-function coupling values for 346 data points, but 347 values for psychopathology and cognition.
# Load the SF-coupling data frame, and find the participant psychopathology and cognition data that we need to remove:
# this is sub-A00037126 BAS1 scan from NKI.
coupling_data_df = pd.read_csv('data/structure.function/coupling.data.df.csv')
nki_coupling_data_df = coupling_data_df[coupling_data_df['dataset'] == 'nki']
participant_to_remove = pd.concat(
    [nki_coupling_data_df.loc[:, ['id', 'dataset', 'timepoint']],
     nki_psych_cog.loc[:, ['id', 'timepoint', 'dataset']]]).drop_duplicates(keep=False)
# Remove this participant from nki_psych_cog. Note that we only index ID because this participant has one time point.
idx_to_remove = nki_psych_cog[nki_psych_cog['id'] == participant_to_remove[['id']].values[0].item()].index[0]
nki_psych_cog.drop(index=idx_to_remove, inplace=True)
# Load the FC DME meta-data sheet, as we need to correct the ages for NKI - at the moment, the age_in_months is
# equivalent to the real age in years rounded down!
nki_fc_dme = pd.read_csv('data/nki/dme/dme.and.metadata.functional.connectivity.csv').rename(
    columns={'bids': 'id', 'session': 'timepoint'})
nki_fc_dme['age_in_months'] = nki_fc_dme['scan_age'] * 12
nki_fc_dme['dataset'] = 'nki'
# Merge the two data frames together!
nki_psych_cog_merged = nki_psych_cog.merge(nki_fc_dme[['id', 'timepoint', 'age_in_months', 'dataset']],
                                           on=['id', 'timepoint'], how='left', suffixes=('', '_new'))
# Drop the helper/merge columns!
nki_psych_cog_merged.drop(columns=['age_in_months', 'dataset_new'], inplace=True)
# Re-name age_in_months_new as age_in_months, to match the CALM dataset
nki_psych_cog_merged.rename(columns={'age_in_months_new': 'age_in_months'}, inplace=True)
# Now merge with the CALM cognitive and psychopathology variables
updated_psych_cog = pd.concat([nki_psych_cog_merged, calm_psych_cog], ignore_index=True)
# And save...
updated_psych_cog.to_csv('data/phenotypic/august25.updated.nki.calm.combined.cognitive.conners.csv')

########################################################
# PART 3 - Dimensions of Cognition and Psychopathology
########################################################
# We use principal component analysis to extract dimensions of cognition using 6 scales common across NKI and CALM, and
# dimensions of psychopathology using 6 Conners scales.
for domain in domains:
    # The function below produces two outputs: one data frame showing the loadings of each measure onto each component,
    # and a second data frame showing the loadings of each participant onto each component.
    if domain == 'cognitive':
        ncomp = 2
    else:
        ncomp = 3
    measure_loadings, participant_loadings = principal_component_analysis(domain=domain, ncomp=ncomp)
    # Save each of these data frames...
    measure_loadings.to_csv('data/phenotypic/' + domain + '.measure.loadings.august25.csv')
    participant_loadings.to_csv('data/phenotypic/' + domain + '.participant.loadings.august25.csv', index=False)
# For each domain, conduct a Mann-Whitney U test to test the null hypothesis that the mean participant loadings for each
# data set are equal.
for domain in domains:
    df = pd.read_csv('data/phenotypic/' + domain + '.participant.loadings.csv')
    # Find how many factors there are
    factors = df.filter(like='Factor').columns
    # For each factor, test for a significant difference in loadings between CALM and NKI.
    for factor in factors:
        statistic, pval = mannwhitneyu(x=df[df['dataset'] == 'calm'][factor], y=df[df['dataset'] == 'nki'][factor])
        print('Testing for difference in %s %s loadings: Mann-Whitney U statistic of %.3f, p-value of %.3f.' %
              (domain, factor, round(statistic, 3), round(pval, 3)))

########################################################
# PART 4 - Assessing impact of missingness
########################################################
# Assess whether the distribution of continuous participant characteristics, such as age, and head motion, for
# participants with and without imputed data. We will follow this tutorial:
# https://medium.com/analytics-vidhya/effectiveness-of-knn-imputation-part-i-the-iris-dataset-e784f157275a
cognitive_psychopathology_data_sheet = pd.read_csv(
        'data/phenotypic/august25.updated.nki.calm.combined.cognitive.conners.csv')
# Add a column to indicate whether the participant has any missing cognitive or psychopathology data.
cognitive_psychopathology_data_sheet['missing_data'] = (
    cognitive_psychopathology_data_sheet.loc[:,'awma_digit_recall_raw': 'conners_peer_relations_t'].isnull().any(axis=1))
n_missing = cognitive_psychopathology_data_sheet.missing_data.sum()
print(f'{n_missing} participants have one or missing cognitive or psychopathology measures, equal to '
      f'{np.round((n_missing/cognitive_psychopathology_data_sheet.shape[0])*100, 2)} percent.')
# Load the structure-function coupling data frame which has the head motion estimates, and merge
coupling_df = pd.read_csv('data/structure.function/coupling.data.df.csv')
cognitive_psychopathology_df = cognitive_psychopathology_data_sheet.merge(coupling_df, on=['timepoint', 'dataset', 'id'])
missing_data_participants = cognitive_psychopathology_df[cognitive_psychopathology_df['missing_data'] == True]
no_missing_data_participants = cognitive_psychopathology_df[cognitive_psychopathology_df['missing_data'] == False]
# We will assess whether participants with missing data have specific age or head motion characteristics
missingness_age_ttest = ttest_ind(no_missing_data_participants.scan_age_y, missing_data_participants.scan_age_y, equal_var=False)
missingness_age_effect = compute_effsize(no_missing_data_participants.scan_age_y, missing_data_participants.scan_age_y)
# We will now assess the effectiveness of the KNN-imputation using root mean squared error.
cognitive_psychopathology_subset = cognitive_psychopathology_data_sheet.iloc[:, 3:15]
# Find how many values are missing of the total nsubject x nmeasure, where 27.50% of rows have missing data
num_val = cognitive_psychopathology_subset.shape[0] * cognitive_psychopathology_subset.shape[1]
percent_nan = (cognitive_psychopathology_subset.isna().sum().sum()/num_val)*100
print(f'{np.round(percent_nan, 3)} percent of the values are missing!')
# Set seed for reproducibility
np.random.seed(123)
# Select rows in cognitive_psychopathology_subset with no missing values
no_missing_subset = cognitive_psychopathology_subset.dropna(axis=0, how='any')
# Create a copy of this into which we will add missing values at random
missing_subset = no_missing_subset.copy()
# Since around 7% of the data is missing, we shall introduce around 10% missing (with replacement)
for feature in cognitive_psychopathology_subset.columns:
    missing_subset.loc[no_missing_subset.sample(frac=0.1, replace=True).index, feature] = np.nan
# Impute the missing data with uniform weights and 5 neighbours
impute = KNNImputer()
missing_subset_knn = impute.fit_transform(missing_subset)
# For each psychopathological and cognitive measure, calculate the normalized root mean square error and express this
# as a percentage of the observable range.
column_wise_rmse = np.zeros((len(cognitive_psychopathology_subset.columns)))
for feature_idx in range(len(cognitive_psychopathology_subset.columns)):
    feature_data_complete = no_missing_subset.iloc[:, feature_idx]
    feature_range = feature_data_complete.max() - feature_data_complete.min()
    rmse = (root_mean_squared_error(feature_data_complete, missing_subset_knn[:, feature_idx])/feature_range)*100
    column_wise_rmse[feature_idx] = rmse
# Create into a data frame and save!
RMSE_accuracy = pd.DataFrame(data=column_wise_rmse, index=cognitive_psychopathology_subset.columns, columns=['RMSE'])
RMSE_accuracy.to_csv('data/phenotypic/psychopathology_cognition_imputation_accuracy.csv')

