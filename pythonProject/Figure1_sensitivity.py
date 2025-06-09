"""
This script covers eLife sensitivity analyses relating to the construction of the affinity matrix, use of Procrustes
rotation, and the stability of the structural and functional gradients across development. Written by Alicja Monaghan in
June 2025.
"""

import os
from mat73 import loadmat
import scipy.io as sio
from brainspace.gradient.kernels import compute_affinity
from brainspace.gradient.embedding import diffusion_mapping
from brainspace.gradient.alignment import procrustes
import numpy as np
from scipy.stats import spearmanr

os.chdir('/Users/alicjamonaghan/Desktop/neurodevelopmental_gradients')


def run_diffusion_map_embedding(input_connectome):
    # This function derives group-level gradients using the same methodology in the paper i.e. for each hemisphere
    # separately, and then align the left hemisphere to the right using a Procrustes rotation. Default parameters are
    # used for the diffusion map embedding algorithm itself. We extract three components in the Schaefer 200-node
    # parcellation. Start by initialising an array to hold the variance explained across each hemisphere, as well as the
    # eigenvectors/values.
    hemi_variance_explained = np.zeros((2, 3))
    hemi_eigenvectors = np.zeros((200, 3))
    for hemisphere in range(2):
        if hemisphere == 0:
            nroi_range = range(0, 100)
        else:
            nroi_range = range(100, 200)
        affinity = compute_affinity(input_connectome[np.ix_(nroi_range, nroi_range)], kernel='normalized_angle')
        diffusion_eigenvectors, diffusion_eigenvalues = diffusion_mapping(affinity, n_components=3)
        # Assign hemisphere-specific eigenvectors to output array
        hemi_eigenvectors[np.ix_(nroi_range), :] = diffusion_eigenvectors
        # Calculate the variance explained by each component in this hemisphere
        hemi_variance_explained[hemisphere, :] = np.squeeze(
            [x / sum(diffusion_eigenvalues) for x in diffusion_eigenvalues])
    for component in range(3):
        # Align the left hemisphere to the right using a Procrustes rotation, and replace in the original array
        rotated_eigenvectors = procrustes(
            hemi_eigenvectors[0:100, component].reshape(100, 1), hemi_eigenvectors[100:, component].reshape(100, 1))
        hemi_eigenvectors[0:100, component] = np.squeeze(rotated_eigenvectors)
        # Once you've looped through both hemispheres, report how much variance is explained by each component
        mean_variance_explained = np.round(np.mean(hemi_variance_explained[:, component]) * 100, 2)
        print(f'Component {component} explains {mean_variance_explained} percent of variance')
    return hemi_eigenvectors


# Load the CALM consensus connectomes
referred_calm_consensus = sio.loadmat('data/calm/connectomes/consensus.mat')['consensus']['sc'].item()['referred']
# Load the group-level gradients
calm_nki_group_gradients = sio.loadmat('data/calm.nki.group.gradients.mat')
# Starting with the referred subset of CALM, log-transform the raw SC matrix, calculate the affinity and then embed.
referred_calm_logsc = np.log(referred_calm_consensus.item()['weighted'].item()['schaefer200x7'].item() + 1)
calm_logsc_gradient = run_diffusion_map_embedding(referred_calm_logsc)
# Correlate the original referred CALM gradients with those derived from log SC
for component in range(3):
    alignment = spearmanr(
        calm_logsc_gradient[:,component], calm_nki_group_gradients['calm']['referred'].item()[0,:,component])
    print(f'Alignment for gradient {component} in CALM is {np.round(alignment[0], 3)}')
# Now load the