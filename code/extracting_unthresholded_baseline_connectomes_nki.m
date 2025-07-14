% This script creates group-representative unthresholded structural and
% functional connectomes for the baseline session of NKI. Written by Alicja
% Monaghan in July 2025, MRC Cognition and Brain Sciences Unit. 

clear;clc;
% Set working directory
cd('//cbsu/data/imaging/projects/external/nkir/analyses/');
% Add path to distance-dependent consensus group thresholding toolbox
% (Betzel et al., 2019)
addpath('/imaging/astle/am10/toolboxes/distanceDependent/');
% Load the structural and functional connectomes
nki_sc = load('/imaging/astle/qsiprep_analysis_da04/qsiprep_analysis_processed/nkir_qsiprep_csd_ACT.mat');
nki_fc = load('connectomes/nki/cleaned_and_parcellated_functional_connectomes.mat');
% Load the DME summary files for structural and functional connectivity
structural_dme = readtable('U:/gradients_open_access/data/nki/dme/dme.and.metadata.structural.connectivity.csv');
functional_dme = readtable('U:/gradients_open_access/data/nki/dme/dme.and.metadata.functional.connectivity.csv');
% Find the IDs matching high-quality low-motion structural(N = 222) and
% functional (N = 213) connectomes at baseline. 
structural_ids = structural_dme(strcmp(structural_dme.session, 'bas1'), 'id');
functional_ids = functional_dme(strcmp(functional_dme.session, 'bas1'), 'id');
% Find the indices of the high-quality structural and functional
% connectomes in the structure holding the connectomes
[~, structuralidx] = ismember(structural_ids.id, nki_sc.nkir_qsiprep_csd_ACT.sample.id.sub);
[~, functionalidx] = ismember(functional_ids.id, nki_fc.nki_fmriprep.sample.sub.bas1);
% Extract the associated connectomes!
structural_connectomes = squeeze(nki_sc.nkir_qsiprep_csd_ACT.schaefer200x7.connectivity(1, structuralidx, 3, :, :));
functional_connectomes = nki_fc.nki_fmriprep.schaefer200x7.connectivity.bas1(functionalidx, :, :);
% Create the group-representative functional connectome (averaging)
group_functional_connectome = squeeze(mean(functional_connectomes, 1));
save('NKI_group_functional_connectome_unthresholded.mat', "group_functional_connectome");
% For the group-representative structural connectome, we'll use the
% distance-dependent consensus thresholding method developed by Betzel and
% colleagues in their 2019 Network Neuroscience paper. 
schaefer200_metadata = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer200x7_1mm_info.mat');
coordinates = [schaefer200_metadata.schaefer200x7_1mm_info.x_mni, schaefer200_metadata.schaefer200x7_1mm_info.y_mni,...
    schaefer200_metadata.schaefer200x7_1mm_info.z_mni];
% Calculate the Euclidean distance between nodes 
D = squareform(pdist(coordinates));
% Set the hemisphere ids (1 = left hemisphere, 2 = right)
hemiid = zeros(200,1);
hemiid(1:200/2,1) = 1;
hemiid(200/2:end,1) = 2;
% Set the dimensions of the structural connectomes as nroi x nroi x nsub
structural_connectomes = permute(structural_connectomes, [3, 2, 1]);
% Create the binary group-representative structural mask
sc_mask = fcn_group_bins(structural_connectomes, D, hemiid, 100);
% Extract the mean nodal streamline count and apply the mask
group_structural_connectome = mean(structural_connectomes, 3).*sc_mask;
save('NKI_group_structural_connectome_unthresholded.mat', "group_structural_connectome");