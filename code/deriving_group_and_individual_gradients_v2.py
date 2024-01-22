# This python script details conducting diffusion map embedding (DME) on structural (weighted communicability) and
# functional connectomes from 3 time points in the nathan kline institute rockland (nki) sample longitudinal discovery
# of brain development trajectories sub-study, and 2 time points from the centre for attention, learning, and memory
# (calm). First, log onto a cluster node, activate neuroconda v2 and load python3. Written by Alicja Monaghan - MRC CBU,
# University of Cambridge. All correspondence to Alicja.Monaghan@mrc-cbu.cam.ac.uk

# PART 1 - Importing Packages and Setting Up Work Space #
import scipy
from scipy.io import savemat
from scipy.stats import spearmanr, shapiro, pearsonr
import os
import numpy as np
import rpy2.robjects as robjects
from brainspace.gradient import diffusion_mapping, compute_affinity, laplacian_eigenmaps, GradientMaps
from brainspace.gradient.alignment import procrustes
from scipy.io import loadmat
import mat73

os.chdir('/home/am10/gradients_open_access/')


# PART 2 - Specify functions!
# This function conducts diffusion map embedding (DME), using the Brain Space library, across different data sets,
# DME parameters (the anisotropic diffusion parameter alpha, diffusion time t, and kernels), modalities (streamline
# counts from probabilistic tractography converted into communicability matrices, alongside resting-state functional
# magnetic resonance imaging connectomes).
def run_group_dimensionality_reduction(dataset, modality, parcellation, ncomp, user_kernel, user_alpha, user_t,
                                       embedding_approach, subset):
    # The data set argument has two options - CALM or NKI. Load the thresholded connectomes for the given modality
    # (sc or fc).
    thresholded_connectomes = mat73.loadmat(
        os.getcwd() + '/data/' + dataset +
        '/connectomes/thresholded_structural_and_functional_connectomes.mat')['thresholded'][modality]
    # Format modality names for output
    if modality == 'sc':
        modality_name = 'structural'
    else:
        modality_name = 'functional'
    # If the data set is calm, extract the harmonised connectomes. We didn't do any harmonisation for NKI as all neuro-
    # imaging data was collected on the same scanner without any software changes.
    if dataset == 'calm':
        thresholded_connectomes = thresholded_connectomes['harmonised']
    group_thresholded_connectome = thresholded_connectomes['group']
    # Load the parcellated group connectomes, using one of 3 options: schaefer100x7, schaefer200x7, or brainnetome246.
    # If using CALM, specify whether we want the referred or non-referred subset using the 'subset' argument.
    if dataset == 'calm':
        group_thresholded_connectome = group_thresholded_connectome[subset][parcellation]
    else:
        group_thresholded_connectome = group_thresholded_connectome[parcellation]
    # Find the number of regions of interest and initialise output arrays for eigenvectors for each hemisphere (left
    # collected first, then right). The user specifies how many components to extract using ncomp.
    nroi = group_thresholded_connectome.shape[0]
    eigenvectors = np.zeros([nroi, ncomp])
    eigenvalues = np.zeros([ncomp, 2])
    # Conduct dimensionality reduction separately for each hemisphere as, otherwise, a left-right hemisphere split
    # would be detected as the principal eigenvector.
    hemi_variance_explained_array = np.zeros(shape=(2, ncomp))
    for hemi_id in range(0, 2):
        if hemi_id == 0:
            hemi_name = "left"
            # Set the region limits
            nroi_range = range(0, int(nroi / 2))
        else:
            hemi_name = "right"
            nroi_range = range(int(nroi / 2), nroi)
        # Construct the affinity matrix using the brainspace module! user_kernel describes the kernel for the affinity.
        affinity = compute_affinity(group_thresholded_connectome[np.ix_(nroi_range, nroi_range)],
                                    kernel=user_kernel, sparsity=None)
        # Check that there are no missing values in the affinity matrix, and raise a warning if there is!
        if np.isnan(affinity).any():
            print("Missing values in {} {} hemisphere affinity in the {} data set.".
                  format(modality_name, hemi_name, dataset))
        # Conduct the dimensionality reduction technique specified in the approach argument: this is diffusion-map
        # embedding ('dm'), laplacian eigenmaps ('le'), or principal components analysis ('pca').
        if embedding_approach == 'dm':
            # For DME, use the specified anisotropic diffusion (user_alpha) and diffusion times (user_t).
            eigenvectors[nroi_range, :], eigenvalues[:, hemi_id] = \
                diffusion_mapping(affinity, n_components=ncomp, alpha=user_alpha,
                                  diffusion_time=user_t, random_state=10)
        elif embedding_approach == 'le':
            # Using default parameters i.e. normalised laplacian
            eigenvectors[nroi_range, :], eigenvalues[:, hemi_id] = \
                laplacian_eigenmaps(affinity, n_components=ncomp, random_state=10)
        else:
            # Also use default parameters for PCA!
            gm_object_fitted = GradientMaps(n_components=ncomp, approach='pca', random_state=10).fit(affinity)
            eigenvectors[nroi_range, :], eigenvalues[:, hemi_id] = \
                gm_object_fitted.gradients_, gm_object_fitted.lambdas_
        # Find how much variance each component explains for each hemisphere!
        for comp in range(0, ncomp):
            hemi_variance_explained_array[hemi_id, comp] = eigenvalues[comp, hemi_id] / sum(eigenvalues[:, hemi_id])
            print("For the {} hemisphere, component {} explains {} percent of variance.".
                  format(hemi_name, comp, hemi_variance_explained_array[hemi_id, comp] * 100))
    # And find the variance explained averaged across both hemispheres
    for comp in range(0, ncomp):
        print("Across hemispheres, component {} explains {} percent of variance.".format(
            comp + 1, round(np.mean(hemi_variance_explained_array[:, comp]) * 100, 2)))
    # When processing structural gradients in CALM, swap the first and second gradient around in the left hemisphere.
    # This is because in CALM, the first and second principal gradients explain a very similar amount of variance,
    # meaning that their order can be switched. We know this because the left and right hemispheres should be mirror
    # images of each other. We don't see this effect in the non-referred subset.
    if dataset == 'calm' and subset == 'referred' and modality == 'sc':
        # Align the first component in the left hemisphere to the second component in the right
        lh_component_1 = procrustes(eigenvectors[range(0, int(nroi / 2)), 0].reshape(100, 1),
                                    eigenvectors[range(int(nroi / 2), nroi), 1].reshape(100, 1))
        lh_component_2 = procrustes(eigenvectors[range(0, int(nroi / 2)), 1].reshape(100, 1),
                                    eigenvectors[range(int(nroi / 2), nroi), 0].reshape(100, 1))
        lh_component_3 = procrustes(eigenvectors[range(0, int(nroi / 2)), 2].reshape(100, 1),
                                    eigenvectors[range(int(nroi / 2), nroi), 2].reshape(100, 1))
        rotated_eigenvectors = np.hstack((lh_component_1, lh_component_2, lh_component_3))
        # Append the rotated left hemisphere eigenvectors to the right hemisphere.
        rotated_eigenvectors = np.concatenate(
            [rotated_eigenvectors, eigenvectors[range(int(nroi / 2), nroi), :]])
    else:
        # Align the left hemisphere to the right using a Procrustes rotation.
        rotated_eigenvectors = \
            procrustes(eigenvectors[range(0, int(nroi / 2)), :], eigenvectors[range(int(nroi / 2), nroi), :])
        # Append the rotated left hemisphere eigenvectors to the right hemisphere.
        rotated_eigenvectors = np.concatenate([rotated_eigenvectors, eigenvectors[range(int(nroi / 2), nroi), :]])
    if dataset == 'calm':
        print("For the {} subset...".format(subset))
    if embedding_approach == 'dm':
        print("Conducted {} {} in the {} data set in the {} parcellation, with alpha of {}, diffusion time of {}, "
              "and {} kernel.".format(modality_name, embedding_approach, dataset, parcellation, user_alpha, user_t,
                                      user_kernel))
    else:
        print("Conducted {} {} in the {} data set in the {} parcellation.".
              format(modality_name, embedding_approach, dataset, parcellation))
    # Now find which regions anchor the gradients! Start by loading the parcellation meta-data
    if parcellation is 'schaefer100x7' or 'schaefer200x7':
        parcellation_metadata = scipy.io.loadmat(
            '/imaging/astle/users/da04/PhD/qsiprep_data/data/' + parcellation + '_1mm_info.mat',
            simplify_cells=True)[parcellation + '_1mm_info']['name']
    else:
        parcellation_metadata = scipy.io.loadmat(
            '/imaging/astle/users/da04/PhD/qsiprep_data/data/' + parcellation + '_info.mat',
            simplify_cells=True)[parcellation]['name']
    # Loop through each component. Find the top 3 regions with the largest positive and negative eigenvectors,
    # respectively.
    nroi = len(parcellation_metadata)
    # The parcellation names are sorted by increasing eigenvector values i.e. largest negative eigenvectors are first,
    # and largest positive eigenvectors are last.
    for comp in range(0, ncomp):
        sorted_parcellation_metadata = parcellation_metadata[np.argsort(rotated_eigenvectors[:, comp])]
        print('{} gradient for component {}: anchored at one end by {}, and at the other by {}.'.format(
            modality_name, comp + 1, ','.join(sorted_parcellation_metadata[0:5]),
            ','.join(sorted_parcellation_metadata[nroi - 5:nroi])))
    del thresholded_connectomes, group_thresholded_connectome, hemi_variance_explained_array
    return rotated_eigenvectors, eigenvalues


# This is a Python translation of a spin-test for parcellated brain data, initially developed by Dr. Frantisek Vasa in
# R (see 'Adolescent Turning of Association Cortex in Human Structural Brain Networks' in Cerebral Cortex, 2018). For
# the toolbox in R, see https://github.com/frantisekvasa/rotate_parcellation
def perm_sphere(x, y, corr_type):
    # perm_id describes the array of permutations from regions to themselves on the sphere, and was generated using the
    # rotate.parcellation function in R. we've saved these mappings as a .csv file, so load this now! Subtract 1 from
    # each perm_id value to account for zero-indexing in Python
    perm_id = np.array(robjects.r.readRDS('data/rotated.parcellation.rds')).astype(int) - 1
    # Get the number of regions of interest (nroi) and permutations
    nroi = perm_id.shape[0]
    nperm = perm_id.shape[1]
    # Empirical correlation between the two parcellated cortical maps X and Y, using inputted correlation method
    if corr_type == 'spearman':
        corr_function = spearmanr
    else:
        corr_function = pearsonr
    empirical_correlation = corr_function(x, y)[0]
    # Initialise array to hold the permuted x values
    x_perm = np.zeros((nroi, nperm))
    y_perm = np.zeros((nroi, nperm))
    # Loop across permutations. For each, loop across each region of interest, and calculate the correlation between the
    # randomly-shuffled region
    for r in range(0, nperm):
        for i in range(0, nroi):
            x_perm[i, r] = x[perm_id[i, r]]
            y_perm[i, r] = y[perm_id[i, r]]
    # Correlation to un-permuted measures
    rho_null_xy = np.zeros(nperm)
    rho_null_yx = np.zeros(nperm)
    for r in range(0, nperm):
        rho_null_xy[r] = corr_function(x_perm[:, r], y)[0]
        rho_null_yx[r] = corr_function(y_perm[:, r], x)[0]
    # Depending on the sign of the empirical correlation, find the permuted p-values for x mapped onto y and vice versa
    if empirical_correlation > 0:
        p_perm_xy = sum(rho_null_xy > empirical_correlation) / nperm
        p_perm_yx = sum(rho_null_yx > empirical_correlation) / nperm
    else:
        p_perm_xy = sum(rho_null_xy < empirical_correlation) / nperm
        p_perm_yx = sum(rho_null_yx < empirical_correlation) / nperm
    # Return average p-value and empirical correlation coefficient
    return list([empirical_correlation, (p_perm_xy + p_perm_yx) / 2])


# This function conducts diffusion map embedding (DME), using the Brain Space library, for each participant for two
# modalities (structural connectivity measured by communicability, and functional connectivity measured by resting-
# state functional magnetic resonance imaging), in the user-define data set (CALM or NKI).
def run_individual_diffusion_map_embedding(dataset, modality, parcellation, ncomp, timepoint, subset, return_affinity):
    # Set the modality index for the input modality ('sc' or 'fc')
    if modality == "sc":
        modality_idx = 0
        formatted_modality_name = "structural"
    else:
        modality_idx = 1
        formatted_modality_name = "functional"
    # Load the group-level diffusion-map embeddings for this data set.
    group_level_dme = loadmat(os.getcwd() + '/data/calm.nki.group.gradients.mat')[dataset]
    # Extract the correct modality, If working with CALM, load the appropriate subset first
    if dataset == 'calm':
        group_level_dme = group_level_dme[subset].item()[modality_idx, :, :]
    else:
        group_level_dme = group_level_dme[modality_idx, :, :]
    # Load the individual-level thresholded connectomes for this data set and modality
    filename = os.getcwd() + '/data/' + dataset + '/connectomes/thresholded_structural_and_functional_connectomes.mat'
    individual_connectomes = mat73.loadmat(filename)['thresholded'][modality]
    # If we're processing CALM, ensure we select the harmonised individual-level connectomes. Further, at baseline only,
    # select the correct subset (referred or non-referred).
    if dataset == 'calm':
        individual_connectomes = individual_connectomes['harmonised']
    if dataset == 'calm' and timepoint == 'baseline':
        participant_list = individual_connectomes[timepoint][subset]['sub']
        print('For the {} subset of CALM...'.format(subset))
    else:
        participant_list = individual_connectomes[timepoint]['sub']
    print('{} participants have {} connectomes for the {} time point of {}.'.
          format(len(participant_list), formatted_modality_name, timepoint, dataset))
    # Initialise output arrays for the eigenvectors and variance explained for the parcellation inputted (
    # 'schaefer100x7', 'schaefer200x7', or 'brainnetome246') and components (ncomp)
    nroi = individual_connectomes[timepoint][parcellation]['individual'].shape[1]
    individual_eigenvectors = np.zeros([nroi, ncomp, len(participant_list)])
    variance_explained = np.zeros([len(participant_list), 2, ncomp])
    affinity_array = np.zeros([len(participant_list), 2, int(nroi / 2), int(nroi / 2)])
    # Note that we specify two hemispheres when collecting the eigen-values because we rotate the eigen vectors to
    # ensure that both hemispheres are aligned, but cannot do this for the eigen values (which are scalar and direction-
    # invariant anyway). Loop across participants...
    for participant_idx, participant_id in enumerate(participant_list):
        # Extract connectome for this participant, modality, time point, and data set. If we're extracting structural
        # connectomes from CALM, then we must select the first index (not noughth) of the fourth dimension, as these
        # are the communicability matrices! Further, select the correct subset for CALM (referred vs non-referred)
        if dataset == "calm" and timepoint == 'baseline' and modality == "sc":
            participant_connectome = \
                individual_connectomes[timepoint][subset][parcellation]['individual'][participant_idx, :, :, 1]
        elif dataset == 'calm' and timepoint == 'followup' and modality == 'sc':
            participant_connectome = \
                individual_connectomes[timepoint][parcellation]['individual'][participant_idx, :, :, 1]
        elif dataset == "calm" and timepoint == 'baseline' and modality == "fc":
            participant_connectome = \
                individual_connectomes[timepoint][subset][parcellation]['individual'][participant_idx, :, :]
        elif dataset == "nki" and modality == "sc":
            participant_connectome = individual_connectomes[timepoint][parcellation]['individual'][participant_idx, :,
                                     :, 1]
        else:
            participant_connectome = \
                individual_connectomes[timepoint][parcellation]['individual'][participant_idx, :, :]
        # Initialise output array for rotated eigenvectors
        rotated_eigenvectors = np.zeros([nroi, ncomp])
        for hemi_id in range(0, 2):
            if hemi_id == 0:
                nroi_range = range(0, int(nroi / 2))
            else:
                nroi_range = range(int(nroi / 2), nroi)
            # Construct the affinity matrix using a normalised angle kernel (in line with Park et al., 2021, eLife)
            affinity = compute_affinity(participant_connectome[np.ix_(nroi_range, nroi_range)],
                                        kernel='normalized_angle', sparsity=None)
            # Append the affinity to the output array
            affinity_array[participant_idx, hemi_id, :, :] = affinity
            # Now apply DME (using default parameters) and assign outputs to relevant arrays
            individual_eigenvectors[nroi_range, :, participant_idx], individual_eigenvalues = \
                diffusion_mapping(affinity, n_components=ncomp, random_state=None)
            # Calculate the variance explained by each component and assign to output
            variance_explained[participant_idx, hemi_id, :] = \
                [i / individual_eigenvalues.sum() for i in individual_eigenvalues]
            # Align the individual's hemispheric eigen-vectors to the corresponding group-level hemisphere DME.
            rotated_eigenvectors[nroi_range, :] = procrustes(individual_eigenvectors[nroi_range, :, participant_idx],
                                                             group_level_dme[nroi_range, :])
        # After processing both hemispheres, check that the rotation and DME has worked!
        if np.isnan(rotated_eigenvectors).any():
            print("{} {} connectivity DME failed, for {} {}.".format(
                participant_id[0], formatted_modality_name, timepoint, dataset))
        # If there is no warning, then assign the rotated eigenvectors to the group output array
        individual_eigenvectors[:, :, participant_idx] = rotated_eigenvectors
        print("{} {} connectivity DME complete, for {} {}.".format(
            participant_id[0], formatted_modality_name, timepoint, dataset))
    # After looping through all participants, return the eigen-vectors, variance explained, and participant lists.
    # Return the affinity matrices too if requested.
    if return_affinity:
        return individual_eigenvectors, variance_explained, participant_list, affinity_array
    else:
        return individual_eigenvectors, variance_explained, participant_list


# PART 3 - Deriving Group Gradients! #
"""
OPEN-ACCESS NOTE: Since both CALM and NKI are managed-access, we cannot provide thresholded connectomes. Therefore, the
following section's code is provided for transparency only. Further, variations in DME parameters are coded in case 
reviewers ask to see the effect of DME parameters or further sensitivity analyses. 
"""
# The function run_group_dimensionality_reduction takes the following arguments:
# * dataset = calm or nki
datasets = ['nki', 'calm']
# * modality = sc or fc. sc denotes structural connectivity, comprised of streamline counts from probabilistic
# tractography, thresholded to retain the strongest 10% connections in each row, then transformed into weighted
# communicability matrices. To derive group-representative structural connectomes, we used Betzel's (2019) distance-
# dependent consensus thresholding, for each data set. fc denotes functional connectivity, namely pair-wise
# z-transformed pearson time-series correlations, with the top 10% connections in each row retained.
# Group-representative fc are simply individual-level functional connectomes averaged across participants.
modalities = ['sc', 'fc']
# * parcellation = We focus our analysis on the Schaefer 200-node 7-network parcellation, in line with prior work
# (Park et al., 2021). For sensitivity, we derive group gradients in a structural parcellation of equivalent spatial
# resolution (Brainnetome 246-node), and a coarser functional parcellation (Schaefer 200-node 7-network).
parcellations = ['schaefer100x7', 'schaefer200x7', 'brainnetome246']
# * ncomp = Number of components to extract. We use the default 10.
# * user_kernel = We assess differences in kernels used to construct the affinity matrix based on recent evidence
# (Watson and Andrews, 2023) that different kernels, particularly cosine, normalized angle, and correlation-based
# kernels, can produce qualitatively different components.
kernels = ['pearson', 'spearman', 'cosine', 'gaussian', 'normalized_angle']
# * user_alpha = Anisotropic diffusion for diffusion-map embedding (DME).
alpha_vector = [0, 0.25, 0.5, 0.75, 1]
# * user_t = Diffusion time for DME.
t_vector = [0, 1, 2, 3]
# * embedding_approach = We consider 3 embedding approaches, namely diffusion-map embedding ('dm'), laplacian eigen-maps
# ('le'), and principal components analysis ('pca').
embedding_approaches = ['dm', 'le', 'pca']
# Using default DME parameters, we'll derive 3 sets of group-level gradients: NKI, referred CALM, and non-referred CALM.
all_subsets = ['nki', 'calm referred', 'calm nonreferred']
default_dme_eigenvectors = np.zeros([len(all_subsets), len(modalities), 200, 3])
for subset_idx, dataset_and_subset in enumerate(all_subsets):
    # If dataset_and_subset contains the string 'calm', split into two and extract the subset name.
    if 'calm' in dataset_and_subset.split():
        dataset, subset = dataset_and_subset.split()
    else:
        dataset = dataset_and_subset.split()[0]
        subset = []
    # Now loop across modalities
    for modality_idx, modality in enumerate(modalities):
        default_dme_eigenvectors[subset_idx, modality_idx, :, :] = \
            run_group_dimensionality_reduction(dataset, subset=subset, modality=modality, parcellation='schaefer200x7',
                                               ncomp=3, user_kernel='normalized_angle', user_alpha=0.5, user_t=0,
                                               embedding_approach='dm')[0]
# Save the group gradients as a MATLAB file for easier importing into R for visualisation.
group_gradients = {"nki": default_dme_eigenvectors[0, :, :, :],
                   "calm": {"referred": default_dme_eigenvectors[1, :, :, :],
                            "non-referred": default_dme_eigenvectors[2, :, :, :]}}
scipy.io.savemat("data/calm.nki.group.gradients.mat", group_gradients)

# PART 4 - Comparing Group-Level CALM and NKI Gradients #
# Load the labels for the schaefer 200-node 7-network parcellation
schaefer200x7_metadata = scipy.io.loadmat('data/schaefer200x7_1mm_info.mat', simplify_cells=True)
schaefer200x7_labels = schaefer200x7_metadata['schaefer200x7_1mm_info']['name']
# Load the group-level gradients, and conduct a spearman-rank correlation for corresponding CALM and NKI gradients.
group_gradients = scipy.io.loadmat('data/calm.nki.group.gradients.mat')
for modality_idx, modality in enumerate(modalities):
    for comp in range(0, 3):
        calm_gradient = group_gradients['calm']['referred'].item()[modality_idx, :, comp]
        nki_gradient = group_gradients['nki'][modality_idx, :, comp]
        rho, pval = perm_sphere(calm_gradient, nki_gradient, corr_type='spearman')
        # Note that the absolute value of the correlation is considered - the sign of eigenvectors is arbitrary
        print("Spearman-rank correlation of {:}, with p-value of {:}, for {:} gradient {:} between CALM and NKI.".
              format(round(rho, 3), round(pval, 3), modality, comp + 1))
        # Find the regions with the smallest and largest absolute dissimilarity in eigenvectors
        diff = abs(calm_gradient - nki_gradient)
        most_similar_regions = schaefer200x7_labels[np.argsort(diff)[0:5:]]
        print("Regions with the largest similarity between CALM and NKI included {}.".format(
            ', '.join(most_similar_regions)))
        most_dissimilar_regions = schaefer200x7_labels[np.argsort(diff)[-5:]]
        print("Regions with the largest dissimilarity between CALM and NKI included {}.".format(
            ', '.join(most_dissimilar_regions)))

# PART 5 - Deriving Individual Gradients and Calculating Manifold Eccentricity #
parcellation = 'schaefer200x7'
# We shall calculate manifold eccentricity for NKI and the referred CALM subset.
for dataset_idx, dataset in enumerate(datasets):
    # Specify time points for each data set, and extract the group gradient. Specify the relevant sub sets, leaving a
    # blank subset variable for NKI.
    if dataset == "calm":
        timepoints = ["baseline", "followup"]
        subset = 'referred'
        dataset_group_gradient = group_gradients[dataset][subset].item()
    else:
        timepoints = ["bas1", "flu1", "flu2"]
        dataset_group_gradient = group_gradients[dataset]
        subset = []
    # Loop across time points
    for timepoint_idx, timepoint in enumerate(timepoints):
        mdic_list = []
        # Specify the file name to which we'll save the DME outputs
        if subset:
            save_filename = 'data/' + dataset + '/dme/' + subset + '_' + parcellation + '_' + timepoint + '.mat'
        else:
            save_filename = 'data/' + dataset + '/dme/' + parcellation + '_' + timepoint + '.mat'
        # Now loop across modalities
        for modality_idx, modality in enumerate(modalities):
            # Format the modality name nicely for saving
            if modality == "sc":
                formatted_modality_name = "structural"
            else:
                formatted_modality_name = "functional"
            # Calculate the group-manifold origin
            group_manifold_origin = np.mean(dataset_group_gradient[modality_idx, :, :], axis=0)
            # Conduct DME for each participant. For CALM, we use the referred subset. Entering this for NKI will not
            # affect results.
            individual_eigenvectors, variance_explained, participants, affinity = \
                run_individual_diffusion_map_embedding(
                    dataset=dataset, modality=modality, parcellation=parcellation, ncomp=3, timepoint=timepoint,
                    subset=subset, return_affinity=True)
            # For each participant's node, calculate the euclidean distance with the group manifold origin i.e. the
            # manifold eccentricity
            nroi = individual_eigenvectors.shape[0]
            nsub = individual_eigenvectors.shape[2]
            manifold_eccentricity = np.zeros([nsub, nroi])
            for sub_idx in range(0, nsub):
                for roi in range(0, nroi):
                    manifold_eccentricity[sub_idx, roi] = np.linalg.norm(
                        group_manifold_origin - individual_eigenvectors[roi, :, sub_idx])
            print("Individual-level {} DME and eccentricities complete for {} {}.".format(modality, timepoint, dataset))
            # After processing all eigenvectors for this modality and time point, create a dictionary
            mdic_modality = {formatted_modality_name + '_eigenvectors': individual_eigenvectors,
                             formatted_modality_name + '_manifold_eccentricity': manifold_eccentricity,
                             formatted_modality_name + '_sub_list': sum(participants, []),
                             formatted_modality_name + '_variance_explained': variance_explained,
                             formatted_modality_name + '_affinity': affinity}
            # Assign this dictionary to the list and delete other variables
            mdic_list.append(mdic_modality)
            print("Appended output to list.")
            del individual_eigenvectors, manifold_eccentricity, mdic_modality, variance_explained, affinity
        # After processing both modalities, concatenate the two dictionaries, and save as a MATLAB file
        savemat(save_filename, mdict={**mdic_list[0], **mdic_list[1]})
