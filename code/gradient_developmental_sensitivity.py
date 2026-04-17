"""
This script covers eLife sensitivity analyses relating to the construction of the affinity matrix, use of Procrustes
rotation, the stability of the structural and functional gradients across development (including head motion), and the
utility of manifold eccentricity beyond traditional graph-theory metrics. Written by Alicja Monaghan in June-September
2025, at the MRC Cognition and Brain Sciences Unit (University of Cambridge).
"""

import os
from mat73 import loadmat
import scipy.io as sio
from brainspace.gradient.kernels import compute_affinity
from brainspace.gradient.embedding import diffusion_mapping
from brainspace.gradient.alignment import procrustes
import numpy as np
from scipy.stats import spearmanr, mode, skew, false_discovery_control, zscore, ttest_rel
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
from nibabel.freesurfer.io import read_annot
import re
from pingouin import mediation_analysis, compute_effsize_from_t
from bct.algorithms import density_und
from netneurotools.metrics import communicability_wei


os.chdir('/Users/alicjamonaghan/Desktop/neurodevelopmental_gradients')


def run_diffusion_map_embedding(input_connectome):
    # This function derives group-level gradients using the same methodology in the paper i.e. for each hemisphere
    # separately, and then align the left hemisphere to the right using a Procrustes rotation. Default parameters are
    # used for the diffusion map embedding algorithm itself. We extract three components in the Schaefer 200-node
    # parcellation. Start by initialising an array to hold the variance explained across each hemisphere, as well as the
    # eigenvectors/values.
    nroi = input_connectome.shape[0]
    hemi_variance_explained = np.zeros((2, 3))
    hemi_eigenvectors = np.zeros((nroi, 3))
    for hemisphere in range(2):
        if hemisphere == 0:
            nroi_range = range(0, int(nroi / 2))
        else:
            nroi_range = range(int(nroi / 2), nroi)
        affinity = compute_affinity(input_connectome[np.ix_(nroi_range, nroi_range)], kernel='normalized_angle')
        diffusion_eigenvectors, diffusion_eigenvalues = diffusion_mapping(affinity, n_components=3, random_state=1)
        # Assign hemisphere-specific eigenvectors to output array
        hemi_eigenvectors[np.ix_(nroi_range), :] = diffusion_eigenvectors
        # Calculate the variance explained by each component in this hemisphere
        hemi_variance_explained[hemisphere, :] = np.squeeze(
            [x / sum(diffusion_eigenvalues) for x in diffusion_eigenvalues])
    for component in range(3):
        # Align the left hemisphere to the right using a Procrustes rotation, and replace in the original array
        rotated_eigenvectors = procrustes(
            hemi_eigenvectors[0:int(nroi / 2), component].reshape(int(nroi / 2), 1),
            hemi_eigenvectors[int(nroi / 2):, component].reshape(int(nroi / 2), 1))
        hemi_eigenvectors[0:int(nroi / 2), component] = np.squeeze(rotated_eigenvectors)
        # Once you've looped through both hemispheres, report how much variance is explained by each component
        mean_variance_explained = np.round(np.mean(hemi_variance_explained[:, component]) * 100, 2)
        print(f'Component {component} explains {mean_variance_explained} percent of variance')
    return hemi_eigenvectors, hemi_variance_explained


def extract_glm_summary(GLM_object):
    # This function extracts the estimates, 95% upper/lower confidence intervals, and p-values for parameters in a GLM.
    estimates = GLM_object.params
    conf_int = GLM_object.conf_int(alpha=0.05)
    p_values = GLM_object.pvalues
    standard_errors = GLM_object.bse
    # Remove the first row, as this corresponds to the intercept (which is significant in all our models)
    summary_df = pd.DataFrame(
        {'estimates': estimates, 'CI_Lower': conf_int[0], 'CI_Upper': conf_int[1],
         'raw_pval': p_values, 'std_err': standard_errors}).iloc[1:]
    summary_df['significant'] = summary_df['raw_pval'].apply(lambda p: 'Yes' if p < 0.05 else 'No')
    return summary_df


def row_wise_thresholding(connectome, threshold):
    # This function thresholds the connectome at the specified thresholded between 0 and 1 e.g. threshold of 0.10 will
    # retain the top 10% connections.
    nroi = connectome.shape[0]
    thresholded_connectome = np.zeros((nroi, nroi))
    n_to_keep = int(np.floor(threshold * nroi))
    for row in range(nroi):
        # Sort descendingly...
        idx = np.argsort((connectome[row, :]))[::-1]
        indices_to_keep = idx[0:n_to_keep]
        thresholded_connectome[row, indices_to_keep] = connectome[row, indices_to_keep]
        # Check if any of the values are negative
        possible_negative_values = np.where(connectome[row, indices_to_keep] < 0)
        if possible_negative_values[0].size != 0:
            print(f'Negative values in row {row}')
    return thresholded_connectome


def gradient_variance_explained(connectome, group_gradients):
    # Calculate the Frobenius norm, which is the size or magnitude of a matrix
    fro_norm_sq = np.linalg.norm(connectome, ord='fro') ** 2
    # Specify the number of regions in the connectome and the number of group-level gradients we're examining
    nroi = connectome.shape[0]
    n_gradients = group_gradients.shape[1]
    # Find how much variance in the connectome is explained by the direction defined by each group gradient
    var_explained = np.zeros(n_gradients)
    for i in range(n_gradients):
        scalar_projection = group_gradients[:, i].T @ connectome @ group_gradients[:, i]
        squared_projection = scalar_projection ** 2
        var_explained[i] = squared_projection / fro_norm_sq
    return var_explained

# Initialise a seed for permutations
np.random.seed(123)
# Load the CALM consensus connectomes
referred_calm_consensus = sio.loadmat('data/calm/connectomes/consensus.mat')['consensus']['sc'].item()['referred']
# Load the group-level gradients
calm_nki_group_gradients = sio.loadmat('data/calm.nki.group.gradients.mat')
# Load the individual thresholded communicability and functional connectivity matrices for NKI and CALM
nki_connectomes = loadmat('data/nki/connectomes/thresholded_structural_and_functional_connectomes.mat')['thresholded']
# Format the NKI meta-data appropriately
nki_metadata = pd.read_excel('data/phenotypic/data-2023-04-09T22_13_49.162Z.xlsx')
nki_metadata['id'] = 'sub-' + nki_metadata.id
calm_connectomes = loadmat('data/calm/connectomes/thresholded_structural_and_functional_connectomes.mat')['thresholded']
# Load the meta-data sheet for all baseline and longitudinal CALM participants: use SC, as this has more participants!
calm_metadata = pd.read_excel('data/calm/Alicja_calm_sc_master.xlsx')
# Load the centroid coordinates for the Schaefer 200-node 7-network parcellation, in MNI space
centroid_coordinates = pd.read_csv('data/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv')
schaefer200_mat = sio.loadmat('data/schaefer200x7_1mm_info.mat', squeeze_me=True)['schaefer200x7_1mm_info']
# Decode bytes to strings
schaefer200_yeo7_labels = np.squeeze([re.findall(r'H_(.*?)_', s) for s in schaefer200_mat['name'].item()])
# Save the network assignments to use in R for visualisation
schaefer200_yeo7_labels_pd = pd.DataFrame(schaefer200_yeo7_labels, columns=['Labels'])
schaefer200_yeo7_labels_pd.to_csv('data/schaefer200_yeo7_labels.csv', index=None)
unique_schaefer200_yeo7_labels = schaefer200_yeo7_labels_pd.Labels.unique()
# Load the DME outputs for the referred CALM cohort
referred_calm_sc_dme = pd.read_csv('data/calm/dme/dme.and.metadata.structural.connectivity.csv')
referred_calm_fc_dme = pd.read_csv('data/calm/dme/dme.and.metadata.functional.connectivity.csv')
# Load the master files containing demographic data for all CALM participants
calm_sc_master = pd.read_excel('data/calm/Alicja_calm_sc_master.xlsx')
calm_fc_master = pd.read_excel('data/calm/Alicja_calm_fc_master.xlsx')
# Load the DME outputs for the NKI cohort
nki_sc_dme = pd.read_csv('data/nki/dme/dme.and.metadata.structural.connectivity.csv')
nki_fc_dme = pd.read_csv('data/nki/dme/dme.and.metadata.functional.connectivity.csv')

#################################
# Parcellate von Economo/Koskinas
#################################
vonEconomo_fsa5 = np.loadtxt('data/economo_koskinas_fsa5.csv')
# Split the von Economo assignments into the left and right hemisphere
vonEconomo_hemisphere_fsa5 = np.zeros((10242, 2))
vonEconomo_hemisphere_fsa5[:, 0] = vonEconomo_fsa5[0:10242]
vonEconomo_hemisphere_fsa5[:, 1] = vonEconomo_fsa5[10242:]
lh_annot = read_annot('data/atl-Schaefer2018_space-fsaverage5_hemi-L_desc-200Parcels7Networks_deterministic.annot')
rh_annot = read_annot('data/atl-Schaefer2018_space-fsaverage5_hemi-R_desc-200Parcels7Networks_deterministic.annot')
annot_files = [lh_annot, rh_annot]
# Initialise an array to hold the von Economo assignments for each hemisphere (LH, then RH)
vonEconomo_assignments = np.zeros((100, 2))
# Each hemisphere has 10242 vertices and 101 unique parcels (0 to 101): 0 indicates the background/medial wall
for hemisphere_idx in range(2):
    for parcel in range(1, 101):
        # Find the indices in the annotation files matching this unique parcel
        idx = np.where(annot_files[hemisphere_idx][0] == parcel)[0]
        # Extract the associated von Economo assignments and remove any 0's
        parcel_ve = vonEconomo_hemisphere_fsa5[idx, hemisphere_idx]
        parcel_ve = parcel_ve[parcel_ve != 0]
        # Assign the mode to the output - we subtract 1 to account for zero-indexing
        vonEconomo_assignments[parcel - 1, hemisphere_idx] = mode(parcel_ve)[0]
# Concatenate vonEconomo_assignments into a single 1D array
vonEconomo_assignments_formatted = np.squeeze(vonEconomo_assignments.reshape(-1, 1))
ve_label_map = {1: 'Agranular', 2: 'Frontal', 3: 'Parietal', 4: 'Polar', 5: 'Granular'}
vonEconomo_assignments_formatted = [ve_label_map[val] for val in vonEconomo_assignments_formatted]
ve_schaefer200_pd = pd.DataFrame(vonEconomo_assignments_formatted, columns=['Labels'])
ve_schaefer200_pd.to_csv('data/vonEconomo_Koskinas_schaefer200.csv', index=False)

#################################
# Impact of raw streamline counts
#################################
# Starting with the referred subset of CALM, log-transform the raw SC matrix, calculate the affinity and then embed.
referred_calm_logsc = np.log(referred_calm_consensus.item()['weighted'].item()['schaefer200x7'].item() + 1)
calm_logsc_gradient = run_diffusion_map_embedding(referred_calm_logsc)
# Correlate the original referred CALM gradients with those derived from log SC
for component in range(3):
    alignment = spearmanr(
        calm_logsc_gradient[0][:, component], calm_nki_group_gradients['calm']['referred'].item()[0, :, component])
    print(f'Alignment for gradient {component} in CALM is {np.round(alignment[0], 3)}')
# Now load the NKI group-level structural connectome
nki_sc_group = sio.loadmat('data/nki/connectomes/nki_group_structural_connectome_schaefer200x7.mat')
# Log-transform the raw SC matrix, calculate the affinity and then embed
nki_logsc = np.log(nki_sc_group['group_distance_sc'] + 1)
nki_logsc_gradient = run_diffusion_map_embedding(nki_logsc)
# Correlate the original NKI gradients with those derived from log SC
for component in range(3):
    alignment = spearmanr(nki_logsc_gradient[0][:, component], calm_nki_group_gradients['nki'][0, :, component])
    print(f'Alignment for gradient {component} in NKI is {np.round(alignment[0], 3)}')

#################################
# Effect of grouping
#################################
# To attempt to replicate the findings produced by Dong and colleagues in their 2021 PNAS paper ("Shifting gradients of
# macro-scale cortical organization mark the transition from childhood to adolescence"), we will divide participants
# based on their age: children are those aged 6 to 12 years old, whilst adolescents are 12 and above. Again, for ease,
# we restrict this to baseline scans.
age_groups = ['child', 'adolescent']
datasets = ['calm', 'nki']
modalities = ['sc', 'fc']
# Store the Spearman correlation coefficients of child- and adult-derived principal and secondary connectivity gradients
grouping_corr_array = np.zeros((len(datasets), len(modalities), 2))
for dataset_idx, dataset in enumerate(datasets):
    # Select the appropriate meta-data sheets
    if dataset == "calm":
        baseline_metadata = calm_metadata[
            calm_metadata.timepoint == 0][['BIDS', 'scan_age']].rename(columns={'BIDS': 'id'})
        baseline_session = "baseline"
    else:
        baseline_metadata = nki_metadata[
            nki_metadata.session == "BAS1"][['id', 'age_04']].rename(columns={'age_04': 'scan_age'})
        baseline_session = "bas1"
    for modality_idx, modality in enumerate(modalities):
        if dataset == "calm":
            baseline_connectome_data = calm_connectomes[modality]['harmonised'][baseline_session]['referred']
        else:
            baseline_connectome_data = nki_connectomes[modality][baseline_session]
        if modality_idx == 0 and dataset == "calm":
            individual_connectomes = baseline_connectome_data['schaefer200x7']['individual']
        elif modality_idx == 0 and dataset == "nki":  # Extract communicability matrices, not streamline counts...
            individual_connectomes = baseline_connectome_data['schaefer200x7']['individual'][:, :, :, 1]
        else:
            individual_connectomes = baseline_connectome_data['schaefer200x7']['individual']
        # Get the correct subject list and create a data frame
        subject_list = np.squeeze(baseline_connectome_data['sub'])
        subject_list_pd = pd.DataFrame(subject_list, columns=['id'])
        # Merge with the meta-data sheet that has participant ages
        metadata_pd = baseline_metadata.merge(subject_list_pd, on='id')
        # Classify participants as either children or adolescents
        metadata_pd['group'] = ['child' if x < 12 else 'adolescent' for x in metadata_pd['scan_age']]
        # Initialize a list to hold the group gradients from the two age-groups
        group_gradient_list = []
        for age_group in age_groups:
            print(f'{modality}: Out of {len(metadata_pd)} participants in {dataset},'
                  f'{metadata_pd.group.value_counts()[age_group]} were {age_group}')
            # Find the indices of connectomes corresponding to this age group
            group_idx = np.where(metadata_pd.group == age_group)[0]
            group_representative_connectome = np.mean(individual_connectomes[group_idx, :, :], axis=0)
            # Derive the first 3 gradients without Procrustes rotation (except for aligning hemispheres)
            group_gradients = run_diffusion_map_embedding(group_representative_connectome)[0]
            # Create a data frame and save
            group_gradient_df = pd.DataFrame(data=group_gradients, columns=[f'G{i}' for i in range(3)])
            group_gradient_df.to_csv('data/sensitivity/' + dataset + '_' + age_group + '_' + modality + '.csv')
            group_gradient_list.append(group_gradient_df)
        # For this specific modality and dataset, extract the Spearman rank correlation coefficient for the concordance
        # between child-derived and adolescent-derived principal gradients
        for i in range(2):
            grouping_corr_array[dataset_idx, modality_idx, i] = spearmanr(
                group_gradient_list[0].iloc[:,i], group_gradient_list[1].iloc[:,i])[0]
d0, d1, d2 = grouping_corr_array.shape
# Create a multi-index of all coordinates
coords = [(i, j, k) for i in range(d0) for j in range(d1) for k in range(d2)]
grouping_corr_df = pd.DataFrame(coords, columns=['dataset', 'modality', 'gradient'])
grouping_corr_df['corr'] = grouping_corr_array.flatten()
grouping_corr_df.to_csv('data/sensitivity/child_adolescent_gradient_comparison.csv', index=False)

########################################################
# Effect of Procrustes alignment on spatial similarity
########################################################
# For each modality (SC and FC), assess the spatial similarity of individual-level gradients with group-level gradients.
for modality_idx, modality in enumerate(modalities):
    individual_connectomes = calm_connectomes[modality]['harmonised']['baseline']['referred']['schaefer200x7'][
        'individual']
    subject_list = np.squeeze(calm_connectomes[modality]['harmonised']['baseline']['referred']['sub'])
    # Create an output array for the spatial alignment
    calm_spatial_alignment = np.zeros((len(subject_list)))
    for sub in range(len(subject_list)):
        # Generate the gradient without Procrustes alignment
        gradients = run_diffusion_map_embedding(individual_connectomes[sub, :, :])
        # Calculate the similarity between the first individual-level and group-level gradient
        calm_spatial_alignment[sub] = abs(
            spearmanr(gradients[0][:, 0], np.squeeze(calm_nki_group_gradients['calm'])
                                          ['referred'].item()[modality_idx, :, 0])[0])
    # After looping through all participants, create a summary data frame
    alignment_summary_pd = pd.DataFrame(calm_spatial_alignment, columns=['principal_alignment'])
    alignment_summary_pd['BIDS'] = subject_list
    # Merge with the relevant meta-data frame, and subset by baseline
    metadata = pd.read_excel('data/calm/Alicja_calm_' + modality + '_master.xlsx')
    metadata_baseline = metadata[metadata['timepoint'] == 0]
    alignment_summary_pd = alignment_summary_pd.merge(metadata_baseline[['BIDS', 'scan_age', 'Sex', 'meanFWD']],
                                                      on='BIDS')
    # Save the alignment data frame so that we can visualise in R!
    alignment_summary_pd.to_csv('data/calm/procrustes_alignment_df_' + modality + '.csv', index=None)
    # Z-score the numeric columns!
    numeric_cols = alignment_summary_pd.select_dtypes(include='number').columns
    alignment_summary_pd[numeric_cols] = zscore(alignment_summary_pd[numeric_cols])
    # Report the Spearman correlation between scan age and principal alignment before accounting for head motion
    alignment_age_without_motion = spearmanr(alignment_summary_pd.scan_age, alignment_summary_pd.principal_alignment)
    print(f'Alignment-age correlation is {np.round(alignment_age_without_motion[0], 3)} '
          f'with p-value of {np.round(alignment_age_without_motion[1], 3)}')
    alignment_glm = smf.glm(formula='principal_alignment ~ scan_age + Sex + meanFWD', data=alignment_summary_pd).fit()
    alignment_glm_summary = extract_glm_summary(alignment_glm)
    alignment_glm_summary.to_csv('data/sensitivity/calm_principal_' + modality + '_alignment_with_age.csv')
    # Report correlation between age and mean FWD
    age_motion_association = spearmanr(alignment_summary_pd.scan_age, alignment_summary_pd.meanFWD)
    print(f'Correlation between age and {modality} motion is {np.round(age_motion_association[0], 3)},'
          f' with p-value of {np.round(age_motion_association[1], 3)}')
# Repeat this pipeline with the NKI baseline data!
for modality_idx, modality in enumerate(modalities):
    if modality is "sc":
        individual_connectomes = nki_connectomes[modality]['bas1']['schaefer200x7']['individual'][:, :, :, 1]
        meanfwd = nki_connectomes[modality]['bas1']['meanfwd']
    else:
        individual_connectomes = nki_connectomes[modality]['bas1']['schaefer200x7']['individual']
        meanfwd = np.squeeze(nki_connectomes[modality]['mean_fwd'][0])
    subject_list = np.squeeze(nki_connectomes[modality]['bas1']['sub'])
    nki_spatial_alignment = np.zeros((len(subject_list)))
    for sub in range(len(subject_list)):
        gradients = run_diffusion_map_embedding(individual_connectomes[sub, :, :])
        nki_spatial_alignment[sub] = abs(
            spearmanr(gradients[0][:, 0], calm_nki_group_gradients['nki'][modality_idx, :, 0])[0])
    alignment_summary_pd = pd.DataFrame(nki_spatial_alignment, columns=['principal_alignment'])
    alignment_summary_pd['id'] = subject_list
    # Add motion estimates!
    alignment_summary_pd['meanFWD'] = meanfwd
    metadata = pd.read_excel('data/phenotypic/data-2023-04-09T22_13_49.162Z.xlsx')
    metadata_baseline = metadata[metadata.session == "BAS1"]
    metadata_baseline.id = 'sub-' + metadata_baseline.id
    alignment_summary_pd = alignment_summary_pd.merge(metadata_baseline[['id', 'age_04', 'dem_002']], on='id')
    # Convert age_04 to a float
    alignment_summary_pd['age_04'] = alignment_summary_pd['age_04'].astype(float)
    # Rename age_04 to scan_age and dem_002 to Sex to match the CALM dataset
    alignment_summary_pd.rename(columns={'age_04': 'scan_age', 'dem_002': 'Sex'}, inplace=True)
    numeric_cols = alignment_summary_pd.select_dtypes(include='number').columns
    alignment_summary_pd[numeric_cols] = zscore(alignment_summary_pd[numeric_cols])
    alignment_age_without_motion = spearmanr(alignment_summary_pd.scan_age, alignment_summary_pd.principal_alignment)
    print(f'Alignment-age correlation is {np.round(alignment_age_without_motion[0], 3)} '
          f'with p-value of {np.round(alignment_age_without_motion[1], 3)}')
    alignment_glm = smf.glm(formula='principal_alignment ~ scan_age + Sex + meanFWD', data=alignment_summary_pd).fit()
    alignment_summary = extract_glm_summary(alignment_glm)
    alignment_summary.to_csv('data/sensitivity/nki_principal_' + modality + '_alignment_with_age.csv')
    # Report correlation between head motion and age!
    age_motion_association = spearmanr(alignment_summary_pd.scan_age, alignment_summary_pd.meanFWD)
    print(f'Correlation between age and {modality} motion is {np.round(age_motion_association[0], 3)},'
          f' with p-value of {np.round(age_motion_association[1], 3)}')
    # Save the alignment data frame for NKI FC as we will test whether head motion mediates the relationship between age
    # and alignment.
    if modality is "fc":
        alignment_summary_pd.to_csv('data/sensitivity/nki_fc_alignment_with_age_df.csv')
# Since head motion for functional scans was negatively correlated with age in NKI, and alignment was significantly
# associated with age, we tested whether the apparent relationship between alignment and age may be explained by motion.
nki_fc_alignment_df = pd.read_csv('data/sensitivity/nki_fc_alignment_with_age_df.csv', index_col=0)
nki_mediation = mediation_analysis(nki_fc_alignment_df, x='scan_age', y='principal_alignment', m='meanFWD', covar='Sex')
mediation_effect_size = (nki_mediation.coef[0] * nki_mediation.coef[1]) / (
        (nki_mediation.coef[0] * nki_mediation.coef[1]) + nki_mediation.coef[2])

#################################
# Methodological thresholding
#################################
# Load the group-level functional connectome for NKI baseline participants
nki_group_fc_unthresholded = sio.loadmat('data/nki/connectomes/NKI_group_functional_connectome_unthresholded.mat')
# Take the strongest 10% of connections as the thresholded example
fc_thresholded = row_wise_thresholding(nki_group_fc_unthresholded['group_functional_connectome'], .10)
# And create a mask where we extract all positive connections, and report the density
nki_fc_positives = nki_group_fc_unthresholded['group_functional_connectome'].copy()
nki_fc_positives[nki_fc_positives < 0] = 0
print(f'Density of all-positive FC is {np.round(density_und(nki_fc_positives)[0] * 100, 2)} %')
# Compare the gradients from the thresholded and non-negative connectomes
fc_thresholded_gradient = run_diffusion_map_embedding(fc_thresholded)
fc_positives_gradient = run_diffusion_map_embedding(nki_fc_positives)
# Now load the group-level structural connectome
nki_group_sc_unthresholded = sio.loadmat('data/nki/connectomes/NKI_group_structural_connectome_unthresholded.mat')
unthresholded_sc_density = density_und(nki_group_sc_unthresholded['group_structural_connectome'])
print(f'Density of un-thresholded SC is {np.round(unthresholded_sc_density[0] * 100, 2)} %')
# Take the strongest 10% of connections and convert to weighted communicability
sc_thresholded = row_wise_thresholding(nki_group_sc_unthresholded['group_structural_connectome'], .10)
sc_thresholded_communicability = communicability_wei(sc_thresholded)
# And calculate the communicability with all streamline counts
sc_unthresholded_communicability = communicability_wei(nki_group_sc_unthresholded['group_structural_connectome'])
# Calculate the communicability gradients for the unthresholded and thresholded structural connectomes
sc_thresholded_gradient = run_diffusion_map_embedding(sc_thresholded_communicability)
sc_unthresholded_gradient = run_diffusion_map_embedding(sc_unthresholded_communicability)

########################################################
# Variance explained in individual FC matrices by group
########################################################
# We assess how much of the variance in individual FC matrices is explained by each of the group-level gradients (G1,
# G2, and G3, for both FC and SC). Start with CALM...
calm_fc = calm_connectomes['fc']['harmonised']['baseline']['referred']['schaefer200x7']['individual']
calm_fc_subject_list = np.squeeze(calm_connectomes['fc']['harmonised']['baseline']['referred']['sub'])
calm_fc_gradients = calm_nki_group_gradients['calm']['referred'].item()[1, :, :]
calm_fc_var_explained = np.zeros((len(calm_fc_subject_list), 3))
for sub in range(len(calm_fc_subject_list)):
    # Multiply by 100 to get the percentage of variance explained
    calm_fc_var_explained[sub, :] = gradient_variance_explained(calm_fc[sub, :, :], calm_fc_gradients) * 100
# Merge with the CALM FC meta-data at baseline
calm_fc_metadata = pd.read_excel('data/calm/Alicja_calm_fc_master.xlsx')
calm_fc_metadata_baseline = calm_fc_metadata[calm_fc_metadata.timepoint == 0]
calm_var_explained_pd = pd.DataFrame(calm_fc_var_explained, columns=['G1', 'G2', 'G3'])
calm_var_explained_pd['BIDS'] = calm_fc_subject_list
calm_var_explained_pd = calm_var_explained_pd.merge(calm_fc_metadata_baseline[['BIDS', 'scan_age', 'Sex', 'meanFWD']], on='BIDS')
# Save this data frame to visualise the developmental change in R!
calm_var_explained_pd.to_csv('data/calm/individual_variance_explained_by_group_CALM_fc.csv', index=None)
# Standardise the numeric columns so that we can extract normalised beta coefficients
numeric_cols = calm_var_explained_pd.select_dtypes(include='number', exclude='category').columns
calm_var_explained_pd[numeric_cols] = calm_var_explained_pd[numeric_cols].apply(zscore)
# Model variance explained by the first group-level gradient using age, sex, and head motion
calm_fc_var_explained_glm = smf.glm('G1 ~ scan_age + Sex + meanFWD', data=calm_var_explained_pd).fit()
calm_fc_var_explained_glm_summary = extract_glm_summary(calm_fc_var_explained_glm)
calm_fc_var_explained_glm_summary.to_csv('data/sensitivity/calm_individual_fc_var_explained.csv', index=None)
# Calculate Cohen's D as an effect size measure using the t-statistic for age
scan_age_d = compute_effsize_from_t(calm_fc_var_explained_glm.tvalues['scan_age'], N = len(calm_fc_subject_list))
# Make sure to report the non-significant associations between age and G2/G3! Then, move onto NKI...
nki_fc = nki_connectomes['fc']['bas1']['schaefer200x7']['individual']
nki_fc_subject_list = np.squeeze(nki_connectomes['fc']['bas1']['sub'])
nki_fc_var_explained = np.zeros((len(nki_fc_subject_list), 3))
for sub in range(len(nki_fc_var_explained)):
    nki_fc_var_explained[sub, :] = gradient_variance_explained(
        nki_fc[sub, :, :], calm_nki_group_gradients['nki'][1, :, :]) * 100
# Merge with the NKI FC meta-data at baseline - note, we have an additional participant in the subject list that is not
# in the NKI FC meta-data
nki_fc_metadata = pd.read_csv('data/nki/dme/dme.and.metadata.functional.connectivity.csv')
nki_fc_metadata_baseline = nki_fc_metadata[nki_fc_metadata.session == "bas1"]
var_explained_pd = pd.DataFrame(nki_fc_var_explained, columns=['G1', 'G2', 'G3'])
var_explained_pd['id'] = nki_fc_subject_list
var_explained_pd = var_explained_pd.merge(nki_fc_metadata_baseline[['id', 'scan_age', 'sex', 'meanfwd']], on='id')
# Model variance explained by the first group-level gradient using age, sex, and head motion
nki_fc_var_explained_glm = smf.glm('G1 ~ scan_age + sex + meanfwd', data=var_explained_pd).fit()
nki_fc_var_explained_glm_summary = extract_glm_summary(nki_fc_var_explained_glm)
# Here, we see that the variance explained by the principal component has a non-significant relationship with age after
# we account for head motion!
nki_fc_var_explained_glm_summary.to_csv('data/sensitivity/nki_individual_fc_var_explained.csv', index=None)

########################################################
# Variance explained in individual SC matrices by group
########################################################
calm_sc_var_explained = np.zeros((len(calm_connectomes['sc']['harmonised']['baseline']['referred']['sub'])))
for sub in range(len(calm_sc_var_explained)):
    # Only extract the variance explained for the first component, as they're all very small! (<2%)
    communicability_connectome = calm_connectomes['sc']['harmonised']['baseline']['referred']['schaefer200x7']['individual'][sub,:,:]
    calm_sc_var_explained[sub] = gradient_variance_explained(
        communicability_connectome, calm_nki_group_gradients['calm']['referred'].item()[0, :, :])[0]*100
# Could be that compression is more effective for FC than SC, possibly because of FC's highly modular structure, whereas
# SC is a lot smaller? Calculate the variability of connectomes and gradients for SC and FC in CALM.

########################################################
# Developmental change in the gradient variance
########################################################
# Load the individual-level gradients for the referred portion of CALM
neurodivergent_calm_dme = sio.loadmat('data/calm/dme/referred_schaefer200x7_baseline.mat')
# Calculate the standard deviation of functional gradients, z-score, and merge with the CALM variance explained data
# frame - the ordering of the gradients is the same!
calm_var_explained_pd['std_raw'] = np.std(neurodivergent_calm_dme['functional_eigenvectors'][:,0,:], axis=0)
calm_var_explained_pd['std_z'] = zscore(calm_var_explained_pd['std_raw'])
# Add the raw scan age (not z-scored) so that we can add a little panel to Figure 1
calm_var_explained_pd = calm_var_explained_pd.merge(calm_fc_metadata_baseline[['BIDS', 'scan_age']], on='BIDS')
calm_var_explained_pd.rename(columns={'scan_age_x': 'scan_age_z', 'scan_age_y': 'scan_age'}, inplace=True)
calm_var_explained_pd.to_csv('data/calm/fc_std_development.csv', index=None)
calm_fc_std_glm = smf.glm('std_z ~ Sex + meanFWD + scan_age_z', data=calm_var_explained_pd).fit()
# Extract Cohen's D as an effect size measure for age
fc_std_age_cohen = compute_effsize_from_t(calm_fc_std_glm.tvalues['scan_age_z'], N = len(calm_fc_subject_list))
# Repeat for the communicability principal gradient
calm_std_sc_raw = np.std(neurodivergent_calm_dme['structural_eigenvectors'][:,0,:], axis=0)
calm_sc_std_pd = pd.DataFrame(
    {'std_raw': calm_std_sc_raw, 'BIDS': neurodivergent_calm_dme['structural_sub_list'],
     'std_z': zscore(calm_std_sc_raw)})
calm_sc_master_baseline = calm_sc_master[calm_sc_master.timepoint == 0]
calm_sc_std_pd = calm_sc_std_pd.merge(calm_sc_master_baseline[['BIDS', 'meanFWD', 'Sex', 'scan_age']], on='BIDS')
# Z-score head motion and scan age so we can extract standardised beta coefficients
calm_sc_std_pd[['meanFWD']] = zscore(calm_sc_std_pd[['meanFWD']])
calm_sc_std_pd['scan_age_z'] = zscore(calm_sc_std_pd.scan_age)
calm_sc_std_pd.to_csv('data/calm/sc_std_development.csv', index=None)
calm_sc_std_glm = smf.glm('std_z ~ Sex + meanFWD + scan_age_z', data=calm_sc_std_pd).fit()
sc_std_age_cohen = compute_effsize_from_t(calm_sc_std_glm.tvalues['scan_age_z'], N = len(calm_sc_std_pd))
# Find common baseline CALM participants with structural and functional gradients.
sc_fc_std_pd = pd.merge(calm_sc_std_pd, calm_var_explained_pd, how='inner', on='BIDS')
sc_fc_std_comparison = ttest_rel(sc_fc_std_pd.std_raw_x, sc_fc_std_pd.std_raw_y)

########################################################
# Psychopathology coverage within CALM and NKI
########################################################
# We excluded high-motion participants, which may have systematically removed those with high psychopathology scores.
# Therefore, report the skewness of each of the 6 Conners sub-scales, and the proportion of youth who have an elevated
# t-score (i.e. a value of 70+ indicates a very elevated score, where the child exhibits many more concerns than are
# typically reported).
calm_nki_imputed_conners = pd.read_csv('data/phenotypic/calm_nki_imputed_psychopathology.csv')
conners_columns = ['conners_inattention_t', 'conners_hyperactivity_impulsivity_t', 'conners_learning_problems_t',
                   'conners_executive_function_t', 'conners_aggression_t', 'conners_peer_relations_t']
for measure in conners_columns:
    # Calculate the Fisher-Pearson coefficient of skewness (g1)
    measure_skewness = skew(calm_nki_imputed_conners[measure])
    print(f'{measure} has a skewness coefficient of {np.round(measure_skewness, 2)}')
    # Calculate the proportion of children with a t-statistic including or exceeding 65, indicating an elevated or very
    # elevated score...
    psychopathology_prop = (len(np.where(calm_nki_imputed_conners[measure] > 65)[0])/len(calm_nki_imputed_conners))*100
    print(f'{np.round(psychopathology_prop, 2)} of children had an elevated or very elevated {measure} score.')

########################################################
# Assessing effect of referral type on gradient values
########################################################
modalities_full = ['structural', 'functional']
calm_master_sheets = [calm_sc_master, calm_fc_master]
for modality_idx in range(2):
    # Create a data frame representing the neurodivergent and neurotypical participants, including a referral column
    neurodivergent_calm_pd = pd.DataFrame(
        {'BIDS': neurodivergent_calm_dme[modalities_full[modality_idx] + '_sub_list'], 'Referred': 1, 'timepoint': 0})
    neurotypical_calm_pd = pd.DataFrame(
        {'BIDS': pd.Series(neurotypical_calm_dme[modalities_full[modality_idx] + '_sub_list']).str.strip(),
         'Referred': 0, 'timepoint': 0})
    gradient_pd = pd.concat((neurodivergent_calm_pd, neurotypical_calm_pd)).merge(calm_master_sheets[modality_idx])
    nsub_referred = gradient_pd.Referred.value_counts()[1]
    nsub_control = gradient_pd.Referred.value_counts()[0]
    print(f'{nsub_referred} referred and {nsub_control} controls for {modalities_full[modality_idx]}')
    gradient_pd = sm.add_constant(gradient_pd)
    gradient_nodal_values = np.concatenate((
        neurodivergent_calm_dme[modalities_full[modality_idx] + '_eigenvectors'][:,0,:],
        neurotypical_calm_dme[modalities_full[modality_idx] + '_eigenvectors'][:,0,:]), axis=1)
    # Z-score each column so that we can extract standardized beta coefficients
    gradient_nodal_values_df = zscore(gradient_nodal_values)
    # Z-score meanFWD as a continuous covariate
    gradient_pd['meanFWD'] = zscore(gradient_pd['meanFWD'])
    # Initialise an array to hold the t-statistic, raw p-value, and FDR-corrected p-values
    nodal_model_statistics = np.zeros((200, 4))
    for node in range(200):
        model = sm.GLM(gradient_nodal_values[node,:], gradient_pd[['const', 'Referred', 'Sex', 'meanFWD', 'scan_age']]).fit()
        nodal_model_statistics[node, 0] = model.tvalues['Referred']
        nodal_model_statistics[node, 1] = model.pvalues['Referred']
        # Compute the effect size from the referral group t-statistic
        nodal_model_statistics[node, 2] = compute_effsize_from_t(model.tvalues['Referred'], nsub_referred, nsub_control)
    # After looping through all nodes, apply FDR corrections, and convert into a data frame
    nodal_model_statistics[:, 3] = false_discovery_control(nodal_model_statistics[:, 1])
    nodal_model_statistics_pd = pd.DataFrame(nodal_model_statistics, columns=['t', 'p', 'cohen', 'FDR_p'])
    nodal_model_statistics_pd['region'] = schaefer200_mat['name']
    nodal_model_statistics_pd.to_csv('data/nodal_value_comparison_' + modalities_full[modality_idx] + '.csv')

########################################################
# Assessing group effects in variance explained
########################################################
# Specify the columns we need for the linear mixed effects modelling!
columns = ['id', 'sex', 'scan_age', 'meanfwd', 'G1_var', 'G2_var', 'G3_var']
for modality_idx, modality in enumerate(modalities):
    calm_dme = pd.read_csv(f'data/calm/dme/dme.and.metadata.{modalities_full[modality_idx]}.connectivity.csv')[columns]
    calm_dme['dataset'] = 'calm'
    nki_dme = globals()[f'nki_{modality}_dme'][columns]
    nki_dme['dataset'] = 'nki'
    # Bind the two data frames together, and make sure categorical variables are encoded as categories!
    dme_summary = pd.concat([calm_dme, nki_dme])
    dme_summary = dme_summary.astype({col: 'category' for col in dme_summary.select_dtypes(include='object').columns})
    dme_summary[dme_summary.select_dtypes(include='number').columns] = dme_summary.select_dtypes(include='number').apply(zscore)
    # Convert all IDs to strings!
    dme_summary['id'] = dme_summary['id'].astype(str).astype('category')
    # Report how many participants have each number of observations!
    print(f"{len(dme_summary.id.unique())} unique {modality} participants.")
    print(dme_summary.id.value_counts().value_counts())
    # For each principal gradient, create a linear mixed effects model with ID as the grouping variable!
    for gradient in range(3):
        mixedlm_formula = f"G{gradient+1}_var ~ meanfwd + sex + scan_age + dataset"
        md = smf.mixedlm(mixedlm_formula, dme_summary, groups=dme_summary['id']).fit()
        # Extract the summary for the dataset effect!
        dataset_effect = extract_glm_summary(md).loc['dataset[T.nki]']
        print(f"Effect of dataset on {modality} G{gradient+1}_var: Beta coefficient of {dataset_effect.estimates:.4},"
              f" with 95% lower CI of {dataset_effect.CI_Lower:.4}, 95% upper CI of {dataset_effect.CI_Upper:.4},"
              f"and p-value of {dataset_effect.raw_pval:.4}")
