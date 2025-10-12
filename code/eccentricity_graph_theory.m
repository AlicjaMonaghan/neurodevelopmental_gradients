% This script calculates graph theory metrics for group-level structural
% and functional connectivity matrices, and compares them with manifold
% eccentricity. Originally written by Dr. Danyal Akarca in March 2024, and
% updated by Alicja Monaghan in August 2025 at the MRC Cognition and Brain
% Sciences Unit, University of Cambridge.

%% Part 1 - Setting up the workspace
clear; clc;
cd('/Users/alicjamonaghan/Desktop/neurodevelopmental_gradients');
% Add paths to the Brain Connectivity Toolbox (R2019)
addpath('2019_03_03_BCT');
% We want to conduct dominance analysis - we'll use the MATLAB toolbox
% developed by Shiyong Yu (168806-dominance-analysis) based from the
% methodology developed by David Budescu (Psychological Bulletin)
addpath('code/');

% Load the gradients and parcellation information
load('data/schaefer200x7_1mm_info.mat');
grad = load('data/calm.nki.group.gradients.mat');
calm = load('data/calm/connectomes/thresholded_structural_and_functional_connectomes.mat');
nki = load('data/nki/connectomes/thresholded_structural_and_functional_connectomes.mat');
% Set parcellation coordinates and name
coordinates = [schaefer200x7_1mm_info.x_mni,schaefer200x7_1mm_info.y_mni,schaefer200x7_1mm_info.z_mni];
name = schaefer200x7_1mm_info.name;

%% Part 2 - Setting up gradient and connectivity arrays
calm_grad = grad.calm.referred;
nki_grad = grad.nki;
% Set the connectivity
calm_sc = calm.thresholded.sc.harmonised.group.referred.schaefer200x7;
calm_fc = calm.thresholded.fc.harmonised.group.referred.schaefer200x7;
nki_sc = nki.thresholded.sc.group.schaefer200x7;
nki_fc = nki.thresholded.fc.group.schaefer200x7;
% Update nki_fc to have the same zero elements as calm
nki_fc(calm_fc==0) = 0;
% Gather hyperparameters
nsamp = 2;
sample_labels = {'CALM', 'NKI'};
nmode = 2;
modality_labels = {'SC', 'FC'};
nnode = 200;
ncomp = 3;
% Form full data structures
gradients = zeros(nsamp,nmode,nnode,ncomp);
gradients(1,:,:,:) = calm_grad;
gradients(2,:,:,:) = nki_grad;
connectivity = zeros(nsamp,nmode,nnode,nnode);
connectivity(1,1,:,:) = calm_sc; connectivity(1,2,:,:) = calm_fc;
connectivity(2,1,:,:) = nki_sc; connectivity(2,2,:,:) = nki_fc;

%% Part 3 - Computing graph theory metrics
nstat = 4;
connectivity_measures = zeros(nsamp,nmode,nnode,nstat);
for samp = 1:nsamp
    for mode = 1:nmode
        w = squeeze(connectivity(samp,mode,:,:));
        % Weighted clustering coefficient
        connectivity_measures(samp,mode,:,1) = clustering_coef_wu(w);
        [Ci,q] = modularity_und(w);
        % Within-module Z-score 
        connectivity_measures(samp,mode,:,2) = module_degree_zscore(w,Ci);
        % Participation coefficient
        connectivity_measures(samp,mode,:,3) = participation_coef(w,Ci); 
        % Calculate nodal strength, and divide by the total number of
        % connections possible (i.e. 200 for the Schaefer 200 parcellation)
        connectivity_measures(samp,mode,:,4) = strengths_und(w)/(200-1);
    end
end
stat_label = {'Weighted clustering coefficient', 'Within-module degree Z score',...
    'Participation coefficient', 'Normalized degree centrality'};
% save('data/graph_theory_measures.mat', 'connectivity_measures');
%% Part 4 - Correlate graph theory metrics with manifold eccentricity
gradients_norm = zeros(nsamp, nmode, nnode);
% Initialise an array to hold the correlation coefficients and p-values
corr_stat = zeros(nsamp, nmode, nstat);
corr_pval = zeros(nsamp, nmode, nstat);
for samp = 1:nsamp
    for mode = 1:nmode
        for node = 1:nnode
            v = squeeze(gradients(samp,mode,node,:));
            gradients_norm(samp,mode,node) = norm(v);
        end
    end
end
% Correlate normalized manifold eccentricity with graph theory measures
for samp = 1:nsamp
    for mode = 1:nmode
        for stat = 1:nstat
            x = squeeze(connectivity_measures(samp, mode, :, stat));
            y = squeeze(gradients_norm(samp, mode, :));
            [corr_stat(samp, mode, stat), corr_pval(samp, mode, stat)] = corr(x, y);
        end
    end
end
eccentricity_corr = struct('stat', corr_stat, 'pval', corr_pval, 'stat_label', stat_label);
save('data/eccentricity_graph_theory_correlations.mat', 'eccentricity_corr');
% For the referred subset of CALM, conduct a dominance analysis to assess
% the relative importance of each of these metrics in accounting for
% variance in structural and functional manifold eccentricity. Start by
% getting the collective variance explained.
X = squeeze(connectivity_measures(1, 1, :, :));
Y = squeeze(gradients_norm(1, 1, :));
collective_Rsq_sc = fitlm(X, Y).Rsquared.Ordinary;
calm_sc_dominance = dominance(X, Y);
writetable(calm_sc_dominance, 'data/calm/calm_sc_dominance.csv');
% Repeat for functional connectivity in CALM
X = squeeze(connectivity_measures(1, 2, :, :));
Y = squeeze(gradients_norm(1, 2, :));
collective_Rsq_fc = fitlm(X, Y).Rsquared.Ordinary;
calm_fc_dominance = dominance(X, Y);
writetable(calm_fc_dominance, 'data/calm/calm_fc_dominance.csv');