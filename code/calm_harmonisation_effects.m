% this script checks for scanner effects in CALM data (centre for
% attention, learning, and memory). written by alicja monaghan,
% alicja.monaghan@mrc-cbu.cam.ac.uk, 20/03/2024

%% part 1 - setting up work space %%
clear;clc;
% set working directory
cd('//cbsu/data/imaging/projects/external/nkir/analyses/');
% add path to brain connectivity toolbox (Rubinov and Sporns, 2010)
addpath('/imaging/astle/am10/toolboxes/2019_03_03_BCT/'); 
% load the harmonised and thresholded connectomes and metadata
thresholded_connectomes = load('/imaging/projects/external/nkir/analyses/duncan_thresholded_structural_and_functional_connectomes.mat');
thresholded_connectomes = thresholded_connectomes.thresholded;
metadata = load('/imaging/projects/external/nkir/analyses/duncan_meta.mat');
metadata = metadata.Meta;
% and the group-level consensus networks
consensus = load('/imaging/projects/external/nkir/analyses/consensus.mat');
% set parameters
modalities = {'sc', 'fc'};
harmonised_status_cell = {'harmonised', 'unharmonised'};
graph_theory_metrics_list = ["efficiency","assortativity","density","transitivity","modularity"];
groups = ["referred", "control"];
% initialise cell array of modalities x harmonisation status
graph_theory_metrics_all_modalities = cell(2,2);
%% part 2 - calculating graph theory metrics %%
for modality_idx = 1:length(modalities)
    modality = modalities{modality_idx};
    for harmonisation_idx = 1:length(harmonised_status_cell)
        harmonisation_status = harmonised_status_cell{harmonisation_idx};
        % concatenate data from baseline, followup (all referred) and
        % control groups
        concat_connectomes = cat(1, thresholded_connectomes.(modality).(harmonisation_status).baseline.schaefer200x7.individual, ...
            thresholded_connectomes.(modality).(harmonisation_status).followup.schaefer200x7.individual, ...
            thresholded_connectomes.(modality).(harmonisation_status).control.schaefer200x7.individual);
        % initialise array to hold graph theory metrics
        graph_theory_metrics = zeros(size(concat_connectomes, 1), length(graph_theory_metrics_list));
        for sub = 1:size(concat_connectomes, 1)
            a = squeeze(concat_connectomes(sub, :, :));
            graph_theory_metrics(sub,1) = efficiency_wei(a,0);
            % assortativity
            graph_theory_metrics(sub,2) = assortativity_wei(a,1);
            % density
            graph_theory_metrics(sub,3) = density_und(a);
            % transitivity
            graph_theory_metrics(sub,4) = transitivity_wu(a);
            % maximised modularity coefficient
            [~, graph_theory_metrics(sub,5)] = modularity_und(a);
            fprintf('Statistics for participant number %.1f using %s %s.\n', sub, harmonisation_status, modality);
        end
        % assign to output cell
        graph_theory_metrics_all_modalities{modality_idx, harmonisation_idx} = graph_theory_metrics; 
    end
end
save('harmonisation_graph_theory_metrics.mat', 'graph_theory_metrics_all_modalities');
%% part 3 - conduct general linear models to assess harmonisation %%
% load graph theory metrics
graph_theory_metrics_all_modalities = load('harmonisation_graph_theory_metrics.mat');
graph_theory_metrics_all_modalities = graph_theory_metrics_all_modalities.graph_theory_metrics_all_modalities;
for modality_idx = 1:length(modalities)
    modality = modalities{modality_idx};
    % concatenate the meta-data for this modality, across time points
    concat_metadata = cat(1, metadata.baseline.(modality), metadata.followup.(modality), metadata.control.(modality));
    for harmonisation_idx = 1:length(harmonised_status_cell)
        harmonisation_status = harmonised_status_cell{harmonisation_idx};
        for metric_idx = 1:length(graph_theory_metrics_list)
            metric = graph_theory_metrics_list(metric_idx);
            % extract the specific graph theory measure across timepoints
            metrics_array = graph_theory_metrics_all_modalities{modality_idx, harmonisation_idx}(:, metric_idx);
            [b,dev,stats]= glmfit(table2array(concat_metadata(:, 2:7)), metrics_array);
            fprintf('%s %s: scanner id had a beta value of %d, with p value of %d, when predicting %s\n', ...
                modality, harmonisation_status, stats.beta(5), stats.p(5), metric);
        end
    end
end
%% part 4 - create group-representative connectomes %
for modality_idx = 1:length(modalities)
    modality = modalities{modality_idx};
    for group_idx = 1:length(groups)
        group = groups(group_idx);
        if group == "referred"
            all_modality_matrices = thresholded_connectomes.(modality).harmonised.baseline.schaefer100x7.individual;
        end
    end
end
%% part 5 - calculate communicability of structural connectomes 
% initialise structure
thresholded = struct();
thresholded.sc.harmonised = struct();
% initialise empty array for baseline individuals
thresholded.sc.harmonised.baseline.referred.schaefer200x7.individual = zeros(size(thresholded_connectomes.sc.harmonised.baseline.schaefer200x7.individual));
% and create communicability matrices
addpath('U:/gradients_open_access/code');
for connectome_idx = 1:length(thresholded.sc.harmonised.baseline.referred.schaefer200x7.individual)
    a = squeeze(thresholded_connectomes.sc.harmonised.baseline.schaefer200x7.individual(connectome_idx, :, :));
    thresholded.sc.harmonised.baseline.referred.schaefer200x7.individual(connectome_idx, :, :) = communicability(a);
end
% repeat for referred individuals
thresholded.sc.harmonised.baseline.nonreferred.schaefer200x7.individual = zeros(size(thresholded_connectomes.sc.harmonised.control.schaefer200x7.individual));
for connectome_idx = 1:size(thresholded.sc.harmonised.baseline.nonreferred.schaefer200x7.individual, 1)
    a = squeeze(thresholded_connectomes.sc.harmonised.control.schaefer200x7.individual(connectome_idx,:,:));
    thresholded.sc.harmonised.nonreferred.schaefer200x7.individual(connectome_idx,:,:) = communicability(a);
end
% add the subject IDs
thresholded.sc.harmonised.baseline.referred.sub = metadata.baseline.sc.BIDS;
thresholded.sc.harmonised.baseline.nonreferred.sub = metadata.control.sc.BIDS;
% and add the meta-data
thresholded.sc.harmonised.baseline.referred.sub_data = metadata.baseline.sc;
thresholded.sc.harmonised.baseline.nonreferred.sub_data = metadata.control.sc;
% repeat for follow-up participants
thresholded.sc.harmonised.followup.schaefer200x7.individual = zeros(size(thresholded_connectomes.sc.harmonised.followup.schaefer200x7.individual));
for connectome_idx = 1:size(thresholded.sc.harmonised.followup.schaefer200x7.individual, 1)
    a = squeeze(thresholded_connectomes.sc.harmonised.followup.schaefer200x7.individual(connectome_idx, :,:));
    thresholded.sc.harmonised.followup.schaefer200x7.individual(connectome_idx,:,:) = communicability(a);
end
% add the subject IDs
thresholded.sc.harmonised.followup.sub = metadata.followup.sc.BIDS;
% and the associated meta-data
thresholded.sc.harmonised.followup.sub_data = metadata.followup.sc;
% now calculate communicability for the group-representative connectome for
% the referred subset
thresholded.sc.harmonised.group.referred.schaefer200x7 = communicability(consensus.consensus.sc.referred.weighted.schaefer200x7);
thresholded.sc.harmonised.group.nonreferred.schaefer200x7 = communicability(consensus.consensus.sc.control.weighted.schaefer200x7);
%% part 6 - format functional connectomes in same way as nki %
thresholded.fc = struct();
% note, we're only populating the structure with harmonised connectomes
thresholded.fc.harmonised = struct();
% repeat, but with functional connectomes!
thresholded.fc.harmonised.baseline.nonreferred.schaefer200x7.individual = thresholded_connectomes.fc.harmonised.control.schaefer200x7.individual;
thresholded.fc.harmonised.baseline.referred.schaefer200x7.individual = thresholded_connectomes.fc.harmonised.baseline.schaefer200x7.individual;
thresholded.fc.harmonised.baseline.referred.sub = metadata.baseline.fc.BIDS;
thresholded.fc.harmonised.baseline.nonreferred.sub = metadata.control.fc.BIDS;
% and add the subject meta-data
thresholded.fc.harmonised.baseline.referred.sub_data = metadata.baseline.fc;
thresholded.fc.harmonised.baseline.nonreferred.sub_data = metadata.control.fc;
thresholded.fc.harmonised.followup.schaefer200x7.individual = thresholded_connectomes.fc.harmonised.followup.schaefer200x7.individual;
thresholded.fc.harmonised.followup.sub = metadata.followup.fc.BIDS;
% and the subject meta-data
thresholded.fc.harmonised.followup.sub_data = metadata.followup.fc;
thresholded.fc.harmonised.group.referred.schaefer200x7 = consensus.consensus.fc.referred.weighted.schaefer200x7;
thresholded.fc.harmonised.group.nonreferred.schaefer200x7 = consensus.consensus.fc.control.weighted.schaefer200x7;
% and save!
save('U:/gradients_open_access/data/calm/connectomes/thresholded_structural_and_functional_connectomes.mat', 'thresholded', '-v7.3');


