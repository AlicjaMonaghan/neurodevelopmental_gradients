# This script calculates graph theory metrics for individual functional and structural connectomes from CALM and NKI.
# We then examine the relationship between mean manifold eccentricity and graph theory metrics, as well as the stability
# of variance explained by individual gradients across time.
# PART 1 - Importing Packages and Setting Up Work Space #
from scipy.stats import spearmanr, shapiro, pearsonr
import pandas as pd
import os
import numpy as np
import rpy2.robjects as robjects
from scipy.io import loadmat
import mat73
from itertools import chain
from bct import modularity_und, transitivity_wu, charpath, clustering_coef_wu, participation_coef, module_degree_zscore
import statsmodels.api as sm
from statsmodels.formula.api import ols
import scikit_posthocs as sp
os.chdir('/home/am10/gradients_open_access/')


# This functions extracts a set of graph theory metrics for thresholded structural (streamline counts) and functional
# connectomes, using a Python adaptation of the MATLAB-based Brain Connectivity Toolbox. Note, whilst we constructed
# structural gradients using communicability, we cannot meaningfully conduct graph theory analyses on those matrices.
def retrieve_graph_theory_metrics_and_manifold_eccentricity(dataset, modality):
    # Load the Schaefer 200-node 7-network parcellation labels.
    schaefer200x7_labels = loadmat('data/schaefer200x7_1mm_info.mat', simplify_cells=True)['schaefer200x7_1mm_info'][
        'name']
    # Assign each node to one of Yeo's (2011) 7 intrinsic functional connectivity networks.
    icn_names = ['Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default']
    icn_vector = np.zeros(200)
    for icn_name_idx, icn_name in enumerate(icn_names):
        matching_index = [idx for idx, s in enumerate(schaefer200x7_labels) if icn_name in s]
        icn_vector[matching_index] = icn_name_idx
    # Extract the thresholded connectomes for the user-specified data set ('calm' or 'nki') and modality ('sc' or 'fc').
    thresholded_connectomes = mat73.loadmat(
        os.getcwd() + '/data/' + dataset +
        '/connectomes/thresholded_structural_and_functional_connectomes.mat')['thresholded'][modality]
    # If CALM, extract harmonised connectomes.
    if dataset == "calm":
        thresholded_connectomes = thresholded_connectomes['harmonised']
    # Extract the time points for this data set, so that we can concatenate individual connectomes
    timepoints = [x for x in list(thresholded_connectomes.keys()) if x not in {'meanfwd', 'mean_fwd', 'sub', 'group'}]
    # Create lists for the local and global graph theory metrics across all time points
    local_graph_theory_metrics_list = []
    global_graph_theory_metrics_list = []
    subids_list = []
    timepoint_list = []
    manifold_eccentricity_list = []
    # Extract connectomes for this time point, then extract graph theory metrics.
    for timepoint in timepoints:
        # For the baseline time point of CALM, select the referred participants.
        if dataset == 'calm' and timepoint == 'baseline' and modality == 'sc':
            connectomes = thresholded_connectomes[timepoint]['referred']['schaefer200x7']['individual'][:, :, :, 0]
        elif dataset == 'calm' and timepoint == 'baseline' and modality == 'fc':
            connectomes = thresholded_connectomes[timepoint]['referred']['schaefer200x7']['individual']
        elif timepoint is not 'baseline' and modality == 'sc':
            connectomes = thresholded_connectomes[timepoint]['schaefer200x7']['individual'][:, :, :, 0]
        else:
            connectomes = thresholded_connectomes[timepoint]['schaefer200x7']['individual']
        # And now extract the participant lists, and format the modality names for later indexing!
        if dataset == 'calm' and timepoint == 'baseline':
            subids = list(chain(*thresholded_connectomes[timepoint]['referred']['sub']))
        else:
            subids = list(chain(*thresholded_connectomes[timepoint]['sub']))
        # Append subids to a larger list
        subids_list.append(subids)
        if modality == 'sc':
            formatted_modality_name = 'structural'
        else:
            formatted_modality_name = 'functional'
        # Extract manifold eccentricity for this data set, timepoint, and modality, and append to corresponding list.
        # Note that we only include the referred subset of CALM baseline. All follow-up participants were referred.
        # However, there are no subset in NKI, therefore select accordingly:
        if dataset == 'calm':
            subset = 'referred_'
        else:
            subset = ''
        manifold_eccentricity_list.append(
            loadmat(os.getcwd() + '/data/' + dataset + '/dme/' + subset + 'schaefer200x7_' +
                    timepoint + '.mat')[formatted_modality_name + '_manifold_eccentricity'])
        # Find out how many participants we have, and initialise empty arrays for local and global graph theory metrics.
        nsub = connectomes.shape[0]
        print('Extracting metrics for {} participants at {} for {}.'.format(nsub, timepoint, dataset))
        global_graph_theory_metrics = np.zeros((nsub, 3))
        local_graph_theory_metrics = np.zeros((nsub, 200, 4))
        # Create a vector where the specific timepoint is repeated nsub times! We need this when creating the summary
        # data frames.
        timepoint_list.append(np.repeat(np.array([timepoint]), nsub, axis=0))
        for sub_idx in range(0, nsub):
            sub_connectome = connectomes[sub_idx, :, :]
            # Modularity coefficient quantifies the extent to which the network can be divided into distinct
            # subnetworks, using Newman's spectral community detection algorithm.
            local_graph_theory_metrics[sub_idx, :, 0], global_graph_theory_metrics[sub_idx, 0] = \
                modularity_und(sub_connectome)
            # Transitivity is the proportion of triangles to triplets - higher transitivity suggests that the network
            # tends to have modules which are highly interconnected.
            global_graph_theory_metrics[sub_idx, 1] = transitivity_wu(sub_connectome)
            # Characteristic path length is the average shortest path length in the network. The smaller this is, the
            # more densely interconnected the network is.
            global_graph_theory_metrics[sub_idx, 2] = charpath(sub_connectome)[0]
            # The higher the nodal clustering coefficient, the greater the proportion of its neighbouring nodes that are
            # neighbours of each other. The corresponding global measure is modularity.
            local_graph_theory_metrics[sub_idx, :, 1] = clustering_coef_wu(sub_connectome)
            # Participation coefficient measures the diversity of the modules which the node's connections belongs to.
            # We specify community structure using Yeo's 7 intrinsic functional connectivity networks (2011), in line
            # with prior structure-function coupling work (Baum and colleagues, 2020, PNAS).
            local_graph_theory_metrics[sub_idx, :, 2] = participation_coef(sub_connectome, icn_vector)
            # Within-module z-score degree centrality measures the number of connections each node within a community
            # has.
            local_graph_theory_metrics[sub_idx, :, 3] = module_degree_zscore(sub_connectome, icn_vector)
        # After processing all participants at this time point, append the local and global metrics to lists
        local_graph_theory_metrics_list.append(local_graph_theory_metrics)
        global_graph_theory_metrics_list.append(global_graph_theory_metrics)
        print('Retrieved {} graph theory metrics for {} {}.'.format(modality, timepoint, dataset))
    # After processing all time points, create two output data frames. The first will have subject IDs, global metrics,
    # and local metrics averaged across nodes. The second will have local metrics averaged across participants. First,
    # concatenate each list to create one large array for local metrics, and another for global.
    local_graph_theory_metrics_long = np.concatenate(local_graph_theory_metrics_list)
    global_graph_theory_metrics_long = np.concatenate(global_graph_theory_metrics_list)
    manifold_eccentricity_long = np.concatenate(manifold_eccentricity_list)
    # Create the data frame with nodal measures
    local_measures_nodal = np.mean(local_graph_theory_metrics_long, axis=0)
    nodal_measures_pd = pd.DataFrame(
        {'Local.Modularity': local_measures_nodal[:, 0], 'Clustering.Coef': local_measures_nodal[:, 1],
         'Participation.Coef': local_measures_nodal[:, 2], 'Within.Module.Z.Degree': local_measures_nodal[:, 3],
         'Mean.Manifold.Eccentricity': np.mean(manifold_eccentricity_long, axis=0)})
    nodal_measures_pd.to_csv(os.getcwd() + '/data/' + dataset + '/connectomes/nodal.' + formatted_modality_name +
                             '.graph.theory.metrics.csv')
    print('Saved local nodal measures averaged across participants.')
    local_measures_participant = np.mean(local_graph_theory_metrics_long, axis=1)
    participant_measures_pd = pd.DataFrame(
        {'Local.Modularity': local_measures_participant[:, 0],
         'Clustering.Coef': local_measures_participant[:, 1], 'Participation.Coef': local_measures_participant[:, 2],
         'Within.Module.Z.Degree': local_measures_participant[:, 3], 'Participant.Manifold.Eccentricity':
             np.mean(manifold_eccentricity_long, axis=1), 'Global.Modularity': global_graph_theory_metrics_long[:, 0],
         'Transitivity': global_graph_theory_metrics_long[:, 1],
         'Char.Path.Length': global_graph_theory_metrics_long[:, 2], 'Timepoint': np.concatenate(timepoint_list),
         'subid': np.concatenate(subids_list)})
    participant_measures_pd.to_csv(os.getcwd() + '/data/' + dataset + '/connectomes/participant.' +
                                   formatted_modality_name + '.graph.theory.metrics.csv')
    print('Saved global and local measures averaged across nodes.')
    return nodal_measures_pd, participant_measures_pd


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


# PART 1 - Calculating Graph Theory Metrics for Individual Connectomes #
# This function will calculate 4 local graph theory metrics (local modularity, clustering coefficient, participation
# coefficient, and within-module z-scored degree) and 3 global metrics (global modularity coefficient, transitivity, and
# characteristic path length). The function saves and outputs two data frames, the first with local graph theory metrics
# averaged across participants (with manifold eccentricity), and the second with local metrics averaged across nodes
# (with global manifold eccentricity). Local metrics are given as participant-averages and nodal-averages.
datasets = ['calm', 'nki']
modalities = ['sc', 'fc']
for dataset in datasets:
    for modality in modalities:
        nodal_measures_pd, participant_measures_pd = retrieve_graph_theory_metrics_and_manifold_eccentricity(
            dataset=dataset, modality=modality)
        # Conduct a permutation spin-test, preserving spatial auto-correlation, for the relationship between manifold
        # eccentricity and local nodal graph theory metrics averaged across participants. We visualise the relationship
        # between brain-wide manifold eccentricity and global graph theory metrics, as well as the local measures, in R.
        local_measures_to_evaluate_nodal_average = \
            list(set(nodal_measures_pd.keys()) - {'ID', 'Mean.Manifold.Eccentricity'})
        print('For {} metrics across all time points in {}...'.format(modality, dataset))
        for measure in local_measures_to_evaluate_nodal_average:
            rho, pval = perm_sphere(nodal_measures_pd['Mean.Manifold.Eccentricity'], nodal_measures_pd[measure],
                                    corr_type='spearman')
            print('Spearman rank correlation of {:.3}, with p as {:.3}, for relationship between mean manifold '
                  'eccentricity and {}.'.format(rho, pval, measure))

# PART 2 - Variance explained by individual-level gradients and their stability #
# For each component, modality, and dataset, report the median variance explained
for dataset in datasets:
    for comp in range(0, 3):
        for modality in modalities:
            if modality is 'sc':
                modality_full = 'structural'
            else:
                modality_full = 'functional'
            df = pd.read_csv('data/' + dataset + '/dme/dme.and.metadata.' + modality_full + '.connectivity.csv')
            print('Median variance explained by component {} in {} connectivity in {} is {} percent, with a range '
                  'between {} and {} percent, respectively.'.format(
                comp + 1, modality_full, dataset, round(df['G' + str(comp + 1) + '_var'].median() * 100, 2),
                round(df['G' + str(comp + 1) + '_var'].min() * 100, 2),
                round(df['G' + str(comp + 1) + '_var'].max() * 100), 2))

# PART 3 - Examining coefficient of variation in manifolds
"""
Open-access note: Make sure to run the R script for Figure 2 so that the coefficient of variation data frame is updated!
"""
# Load the data frame holding the coefficient of variation across nodes, modalities, and data sets.
cv_df = pd.read_csv('data/coefficient.of.variation.csv')
test_of_normality = shapiro(cv_df['eccentricity.cv'])
print('The Shapiro-Wilk test statistic is %.2f, with a p-value of %.2f.' % (
    round(test_of_normality[0], 2), round(test_of_normality[1], 2)))
# Since the data is not normally distributed, we'll conduct a non-parametric equivalent of a two-way ANOVA, called the
# Scheirer-Ray-Hare test. First, rank all input elements.
cv_df['rank'] = cv_df['eccentricity.cv'].rank()
# Recode numerical factors into categorical
cv_df['dataset'] = np.where(cv_df['dataset'] == 1, "calm", "nki")
cv_df['modality'] = np.where(cv_df['modality'] == 1, "structural", "functional")
# Ensure that eccentricity.cv has an underscore
cv_df.rename(columns={'eccentricity.cv': 'eccentricity_cv'}, inplace=True)
# Run a two-way ANOVA using the ranked data
model = ols('eccentricity_cv ~ C(dataset) + C(modality) + C(dataset):C(modality)', data=cv_df).fit()
sm.stats.anova_lm(model, typ=2)
# Apply Dunn's Test to determine which means are different from the rest.
sp.posthoc_dunn(cv_df, p_adjust='fdr_bh', group_col='Group', val_col='eccentricity_cv')
# Get the mean manifold eccentricity CV for each group
cv_df.groupby(['Group'])['eccentricity_cv'].mean()
