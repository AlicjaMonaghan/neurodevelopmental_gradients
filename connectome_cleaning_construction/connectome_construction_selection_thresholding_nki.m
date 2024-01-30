% this script details selecting high-quality structural and functional 
% connectomes for participants across 3 time points of the longitudinal 
% nathan-kline institute (nki) study. 

%% part 1 - setting up the work space! %%
clear;clc;
% set working directory
cd('/imaging/projects/external/nkir/analyses/');
% add path to brain connectivity toolbox (Rubinov and Sporns, 2010)
addpath('/imaging/astle/am10/toolboxes/2019_03_03_BCT/');
% add path to distance-dependent consensus group thresholding toolbox
% (Betzel et al., 2019)
addpath('/imaging/astle/am10/toolboxes/distanceDependent/');
% load the cleaned and parcellated functional connectomes
nki_fc = load('connectomes/nki/cleaned_and_parcellated_functional_connectomes.mat');
timepoints_list = ["bas1","flu1","flu2"];
parcellations = ["schaefer100x7","schaefer200x7","brainnetome246"];
%% part 2 - selecting participants with high-quality functional connectomes %%
% for each participant and at each session, calculate mean framewise
% displacement, and find participants for whom over 20% of their spikes are
% high-motion (i.e. framewise displacement > .5). 
for timepoint_idx = 1:length(timepoints_list)
    % find the original number of functional connectomes for this timepoint
    sub_list = nki_fc.nki_fmriprep.sample.sub.(timepoints_list(timepoint_idx));
    fprintf('for ses-%s of the longitudinal nki child brain study, %d participants had functional connectomes.\n', ...
        timepoints_list{timepoint_idx},length(sub_list));
    % initialise output array for mean framewise displacement and
    % proportion of high-motion spikes out of all spikes
    fd_qc = zeros(length(sub_list),2);
    for sub_idx = 1:length(sub_list)
        fd = nki_fc.nki_fmriprep.framewisedisplacement.(timepoints_list(timepoint_idx)){sub_idx};
        % calculate mean framewise displacement
        fd_qc(sub_idx,1) = mean(fd);
    end
    % find participants with high mean framewise displacement
    high_mean_fd = find(fd_qc(:,1)>.5);
    % remove these participants from the parcellated connectomes
    for parcellation_idx = 1:length(parcellations)
        nki_fc.nki_fmriprep.(parcellations(parcellation_idx)).connectivity.(timepoints_list(timepoint_idx))(high_mean_fd,:,:) = [];
    end
    % and remove corresponding cells in the framewise displacement array
    % and subject list
    nki_fc.nki_fmriprep.framewisedisplacement.(timepoints_list(timepoint_idx))(high_mean_fd) = [];
    sub_list(high_mean_fd) = [];
    % now find participants with a high proportion of high-motion spikes
    for sub_idx = 1:length(sub_list)
        fd = nki_fc.nki_fmriprep.framewisedisplacement.(timepoints_list(timepoint_idx)){sub_idx};
        fd_qc(sub_idx,2) = sum(fd>.5)/length(fd);
    end
    % and participants with a high proportion of high-motion spikes
    high_motion_spikes_participants = find(fd_qc(:,2)>.2);
    % and remove from parcellated connectomes
    for parcellation_idx = 1:length(parcellations)
        nki_fc.nki_fmriprep.(parcellations(parcellation_idx)).connectivity.(timepoints_list(timepoint_idx))(high_motion_spikes_participants,:,:) = [];
    end
    % update subject list and framewise displacement 
    sub_list(high_motion_spikes_participants) = [];
    nki_fc.nki_fmriprep.framewisedisplacement.(timepoints_list(timepoint_idx))(high_motion_spikes_participants) = [];
    fprintf(['for ses-%s of nki, we removed %d participants with high mean fd, and a' ...
        ' further %d participants with a high proportion of high-motion spikes.\n'], ...
        timepoints_list{timepoint_idx},length(high_mean_fd),length(high_motion_spikes_participants))
    % re-calculate mean framewise displacement 
    fd_qc(high_mean_fd,:) = [];
    fd_qc(high_motion_spikes_participants,:) = [];
    fprintf(['for ses-%s of nki, %d participants had high-quality functional connectomes, ' ...
        'with mean framewise displacement of %f., and standard deviation of %f.\n'], ...
        timepoints_list{timepoint_idx},length(sub_list),mean(fd_qc(:,1)),std((fd_qc(:,1))));
    % remove the participants from the sub in the structure
    nki_fc.nki_fmriprep.sample.sub.(timepoints_list(timepoint_idx))(high_mean_fd) = [];
    nki_fc.nki_fmriprep.sample.sub.(timepoints_list(timepoint_idx))(high_motion_spikes_participants) = [];
end
%% part 3 - thresholding functional connectomes %%
% in line with margulies and colleagues (2016), for each row, retain the
% top 10% absolute strongest connections. first, create a new structure to 
% hold the thresholded connectomes.
thresholded = struct();
thresholded.fc = struct();
thresholded.fc.sub = cell(3,1);
thresholded.fc.mean_fwd = cell(3,1);
for timepoint_idx = 1:length(timepoints_list)
    % initialise cell to hold IDs of participants whose thresholding was
    % unsuccessful at this time point, for any parcellation i.e. not 
    % enough non-zero values to take the top 10%.
    failed_fc_thresholding = cell(length(parcellations),1);
    thresholded_connectomes_cell = cell(length(parcellations),1);
    for parcellation_idx = 1:length(parcellations)
        unthresholded_connectomes = nki_fc.nki_fmriprep.(parcellations(parcellation_idx)).connectivity.(timepoints_list(timepoint_idx));
        nroi = size(unthresholded_connectomes,2);
        nsub = size(unthresholded_connectomes,1);
        mean_fwd = zeros(nsub,1);
        thresholded_connectomes = zeros(nsub,nroi,nroi);
        for sub_idx = 1:nsub
            sub_connectivity = abs(squeeze(unthresholded_connectomes(sub_idx,:,:)));
            % extract the framewise displacement for this participant and
            % calculate the mean.
            mean_fwd(sub_idx,1) = mean(nki_fc.nki_fmriprep.framewisedisplacement.(timepoints_list(timepoint_idx)){sub_idx});
            % check the density of these connectomes. for fc, these should
            % be fully connected. 
            if density_und(sub_connectivity) ~= 1
                sprintf('functional connectome for %s in the %s parcellation at %s has a density of %.3f percent.', ...
                    nki_fc.nki_fmriprep.sample.sub.(timepoints_list{timepoint_idx}){sub_idx}, parcellations(parcellation_idx), ...
                    timepoints_list(timepoint_idx), density_und(sub_connectivity));
            end
            % for each row, find the 10% absolute strongest
            % connections, and set all others to 0
            for row = 1:nroi
                [sorted_values,I] = sort(sub_connectivity(row,:),'descend');
                % if the number of regions we have does not produce an
                % integer when dividing by 10 (i.e. 90th percentile),
                % round to 2 significant figures
                indices_to_keep = I(1:floor(round(nroi*.10)));
                thresholded_connectomes(sub_idx,row,indices_to_keep) = sorted_values(1:floor(round(nroi*.10)));
                if nnz(squeeze(thresholded_connectomes(sub_idx,row,:))) ~= round(nroi*.10)
                    fprintf("row %d for sub idx of %d does not have enough non-zero values.\n", ...
                        row, sub_idx);
                    % Add this participant to the failed_fc_thresholding
                    % cell array!
                    failed_fc_thresholding{parcellation_idx,1} = cat(1,nki_fc.nki_fmriprep.sample.sub.(timepoints_list{timepoint_idx}){sub_idx}, ...
                        failed_fc_thresholding{parcellation_idx});
                end
            end
        end
        fprintf('thresholded ses-%s functional connectomes in the %s parcellation.\n', ...
            timepoints_list{timepoint_idx}, parcellations(parcellation_idx));
        % assign the thresholded connectomes to the output cell
        thresholded_connectomes_cell{parcellation_idx,1} = thresholded_connectomes;
    end
    % after processing all parcellations for this time point, see if there
    % are any participants who failed the thresholding procuedre. if so,
    % remove them.
    if isempty(failed_fc_thresholding{timepoint_idx}) == 0
        failed_participants = unique(string(failed_fc_thresholding{timepoint_idx}));
        failed_participant_idx = find(strcmp(string(nki_fc.nki_fmriprep.sample.sub.(timepoints_list(timepoint_idx))), failed_participants));
        for parcellation_idx = 1:length(parcellations)
            thresholded_connectomes_cell{parcellation_idx,1}(failed_participant_idx,:,:) = [];
        end
        % And from the mean FWD cell
        mean_fwd(failed_participant_idx,:) = [];
        % And from the participant list in the nki_fc structure!
        nki_fc.nki_fmriprep.sample.sub.(timepoints_list(timepoint_idx))(failed_participant_idx) = [];
        fprintf('Removed %s from %s in the %s parcellation because their connectome was not fully connected.\n', ...
            failed_participants,timepoints_list(timepoint_idx),parcellations(parcellation_idx));
    end
    % assign the thresholded connectomes to the output
    for parcellation_idx = 1:length(parcellations)
        thresholded.fc.(timepoints_list{timepoint_idx}).(parcellations(parcellation_idx)).individual = thresholded_connectomes_cell{parcellation_idx,1};
    end
    clear thresholded_connectomes unthresholded_connectomes
    % add sub ids and mean framewise displacement
    thresholded.fc.(timepoints_list{timepoint_idx}).sub = nki_fc.nki_fmriprep.sample.sub.(timepoints_list(timepoint_idx));
    thresholded.fc.mean_fwd{timepoint_idx,:} = mean_fwd;
end
% across parcellations, concatenate across time points. 
for parcellation_idx = 1:length(parcellations)
    concatenated_fc = cat(1, thresholded.fc.bas1.(parcellations(parcellation_idx)).individual, ...
        thresholded.fc.flu1.(parcellations(parcellation_idx)).individual,...
        thresholded.fc.flu2.(parcellations(parcellation_idx)).individual);
    % check that there are no missing values in the array!
    if isempty(find(isnan(concatenated_fc), 1)) == 1
        fprintf("No missing values for thresholded functional connectomes " + ...
            "across all time points in the %s parcellation.\n", ...
            parcellations(parcellation_idx));
    end
    % average across participants to create the group-representative
    % functional connectome
    thresholded.fc.group.(parcellations(parcellation_idx)) = squeeze(mean(concatenated_fc));
    clear concatenated_fc
end
clear nki_fc
%% part 4 - select participants with high-quality structural connectomes %%
nki_sc = load('/imaging/astle/qsiprep_analysis_da04/qsiprep_analysis_processed/nkir_qsiprep_csd_ACT.mat');
high_quality_sc_nkir.sub = cell(3,1);
high_quality_sc_nkir.mean_fwd = cell(3,1);
for timepoint_idx = 1:length(timepoints_list)
    high_quality_sc_nkir.(timepoints_list(timepoint_idx)) = struct();
    % find the original number of structural connectomes for this session.
    % we're using the sift_invnodevol_radius2_count_connectivity field as 
    % our structural connectivity measure (index of 3). 
    sub_indices_into_struct = find(nki_sc.nkir_qsiprep_csd_ACT.tracker.tracker(:,timepoint_idx)==1);
    sub_list = nki_sc.nkir_qsiprep_csd_ACT.sample.id.sub(sub_indices_into_struct)';
    fprintf('for ses-%s of longitudinal nki, %d participants had structural connectomes.\n', ...
        timepoints_list{timepoint_idx},length(sub_list));
    % initialise output array for mean framewise displacement
    fd_qc = zeros(length(sub_list),1);
    for sub_idx = 1:length(sub_list)
        fd_participants_timepoint_data = squeeze(nki_sc.nkir_qsiprep_csd_ACT.qc.summary(timepoint_idx,sub_indices_into_struct,:));
        % assign mean framewise displacement
        fd_qc(sub_idx,1) = fd_participants_timepoint_data(sub_idx,1);
    end
    % find participants with high mean framewise displacement. note that we
    % use a threshold of 3mm here, rather than .5mm in the case of
    % functional connectomes. we're using the
    % sift_invnodevol_radius2_count_connectivity field as our structural
    % connectivity measure (index of 3).
    high_mean_fd = find(fd_qc(:,1)>3);
    % create an array to index participants whose structural connectome
    % density is more than 3 standard deviations lower than the mean
    low_density_participants = cell(length(parcellations),1);
    connectome_density = zeros(length(parcellations),nsub,1);
    for parcellation_idx = 1:length(parcellations)
        high_quality_sc_nkir.(timepoints_list(timepoint_idx)).(parcellations(parcellation_idx)) = struct();
        timepoint_all_matrices = squeeze(nki_sc.nkir_qsiprep_csd_ACT.(parcellations(parcellation_idx)).connectivity(timepoint_idx,sub_indices_into_struct,3,:,:));
        % first, remove participants with high head motion
        timepoint_all_matrices(high_mean_fd,:,:) = [];
        % for each participant, find their structural connectome density.
        % the density_und function from the brain connectivity toolbox
        % occasionally provides density values greater than 1. therefore,
        % we shall calculate density ourselves, as the number of
        % off-diagonal non-zero elements divided by the number of
        % connections (excluding self-connections).
        nsub = size(timepoint_all_matrices,1);
        nroi = size(timepoint_all_matrices,2);
        for sub_idx = 1:nsub
            sub_connectivity = squeeze(timepoint_all_matrices(sub_idx,:,:));
            % Set diagonal to zero
            sub_connectivity = sub_connectivity - diag(diag(sub_connectivity));
            connectome_density(parcellation_idx,sub_idx,1) = nnz(sub_connectivity) / ((nroi*nroi)-nroi);
        end
        % find which participants have a structural connectome density less
        % than 3 standard deviations than the mean, and assign to the
        % output cell array.
        low_density_participants{parcellation_idx,1} = find(connectome_density(parcellation_idx,:,:) < ...
            mean(connectome_density(parcellation_idx,:,:)) - 3*std(connectome_density(parcellation_idx,:,:)));
    end
    % find the participants with low-density connectomes across one or more
    % parcellations
    low_density_participants_unique = unique(cell2mat(low_density_participants));
    fprintf(['removing %d participants at the %s time point whose structural connectome ' ...
        'density was more than 3 sd less than the mean.\n'], ...
        length(low_density_participants_unique),timepoints_list(timepoint_idx));
    % remove high-motion participants from the framewise displacement array
    % and subject list 
    fd_qc(high_mean_fd,:) = [];
    sub_list(high_mean_fd,:) = [];
    % and remove low density participants
    fd_qc(low_density_participants_unique,:) = [];
    sub_list(low_density_participants_unique,:) = [];
    % across each parcellation, remove these participants' connectomes, and
    % report the new framewise displacements
    for parcellation_idx = 1:length(parcellations)
        % remove these participants from the structural connectome matrices
        high_quality_sc_nkir.(timepoints_list(timepoint_idx)).(parcellations(parcellation_idx)) = struct();
        timepoint_all_matrices = squeeze(nki_sc.nkir_qsiprep_csd_ACT.(parcellations(parcellation_idx)).connectivity(timepoint_idx,sub_indices_into_struct,3,:,:));
        % first, remove participants with high head motion
        timepoint_all_matrices(high_mean_fd,:,:) = [];
        % and now those with low-density connectomes
        timepoint_all_matrices(low_density_participants_unique,:,:) = [];
        % check that there are no missing values!
        if isempty(find(isnan(timepoint_all_matrices), 1)) == 1
            fprintf("No missing values for unthresholded structural " + ...
                "connectomes across all time points in the %s parcellation.\n", ...
                parcellations(parcellation_idx));
        end
        % and assign to a new structure
        high_quality_sc_nkir.(timepoints_list(timepoint_idx)).(parcellations(parcellation_idx)).connectivity = timepoint_all_matrices;
        clear timepoint_all_matrices
    end
    % re-calculate mean framewise displacement 
    fprintf('for ses-%s of nki, %d participants had high-quality structural connectomes, with mean framewise displacement of %f, and sd of %f.\n', ...
        timepoints_list{timepoint_idx},length(sub_list),mean(fd_qc(:,1)),std(fd_qc(:,1)));
    % add sub list and framewise displacement to the new sc structure
    high_quality_sc_nkir.sub{timepoint_idx,1} = sub_list;
    high_quality_sc_nkir.mean_fwd{timepoint_idx,1} = fd_qc;
    clear sub_list fd_qc
end
% also save the mean framewise displacement for these participants, to be
% added as a covariate in generalised additive mixed models to assess the
% relationship between age and manifold eccentricity.
fd_metadata = struct();
fd_metadata.sc.meanfwd = high_quality_sc_nkir.mean_fwd;
fd_metadata.fc.meanfwd = thresholded.fc.mean_fwd;
fd_metadata.sc.sub = high_quality_sc_nkir.sub;
% Ensure that all sub IDs are cell arrays!
for ses_idx = 1:3
    fd_metadata.sc.sub{ses_idx} = cellstr(fd_metadata.sc.sub{ses_idx});
end
fd_metadata.fc.sub = thresholded.fc.sub;
clear nki_sc
%% part 5 - thresholding structural connectomes %%
% to allow a similar density to the functional connectomes, whilst ensuring
% a connected network, for each row, retain the top 10% strongest
% connections. 
thresholded.sc = struct();
thresholded.sc.sub = cell(3,1);
thresholded.sc.mean_fwd = cell(3,1);
for timepoint_idx = 1:length(timepoints_list)
    % initialise an array to hold the IDs of participants who have problems
    % with their thresholding/communicability e.g. missing values, which
    % makes diffusion-map embedding not possible
    problematic_thresholding_subs = cell(length(parcellations),1);
    % and create cell array for thresholded connectomes across all
    % parcellations for this time point
    thresholded_connectomes_cell = cell(length(parcellations),1);
    % loop over parcellations...
    for parcellation_idx = 1:length(parcellations)
        % extract unthresholded connectomes for this time point and
        % parcellation
        unthresholded_connectomes = high_quality_sc_nkir.(timepoints_list(timepoint_idx)).(parcellations(parcellation_idx)).connectivity;
        fprintf('thresholding ses-%s structural connectomes in the %s parcellation.\n', ...
            timepoints_list{timepoint_idx},parcellations(parcellation_idx));
        % set the number of regions of interest (nroi)
        nroi = size(unthresholded_connectomes,2);
        % find the number of participants (nsub)
        nsub = size(unthresholded_connectomes,1);
        % new arrays for thresholded connectomes (first slice of 4th
        % dimension) and communicability (2nd slice of 4th dimension).
        thresholded_connectomes_array = zeros(nsub,nroi,nroi);
        for sub_idx = 1:nsub
            sub_connectivity = squeeze(unthresholded_connectomes(sub_idx,:,:));
            % for each row, find the 10% absolute strongest
            % connections, and set all others to 0
            for row = 1:nroi
                [~,I] = sort(sub_connectivity(row,:),'descend');
                I = I';
                indices_to_keep = I(1:floor(nroi*.10));
                thresholded_connectomes_array(sub_idx,row,indices_to_keep,1)  = sub_connectivity(row,indices_to_keep);
            end
            % now convert to communicability matrices. first find the
            % strength matrix.
            a = squeeze(thresholded_connectomes_array(sub_idx,:,:,1));
            s = diag(sum(a,2));
            pow = (s^-.5)*a*(s^-.5);
            c = expm(pow);
            % assign c to output
            thresholded_connectomes_array(sub_idx,:,:,2) = c;
            % if there are no missing values in c, assign to output.
            if isempty(find(isnan(c), 1)) == 0
                fprintf(['Missing values for weighted communicability for ' ...
                    'participant with sub index of %d, at %s in the %s parcellation.\n'], ...
                    sub_idx,timepoints_list(timepoint_idx),parcellations(parcellation_idx));
                subject_id = high_quality_sc_nkir.sub{timepoint_idx,1}(sub_idx);
                problematic_thresholding_subs{parcellation_idx,1} = cat(2,problematic_thresholding_subs{parcellation_idx,1}, subject_id);
            end
            clear a s pow c
        end
        % assign the thresholded_connectomes_array to corresponding cell
        thresholded_connectomes_cell{parcellation_idx,1} = thresholded_connectomes_array;
        fprintf('thresholding for ses-%s in the %s parcellation complete.\n', ...
            timepoints_list(timepoint_idx), parcellations(parcellation_idx));
    end
    % now that we've thresholded all structural connectomes and calculated
    % communicability across all parcellations for this time point, remove 
    % participants with missing communicability values across one or
    % more parcellations at this time point
    if isempty(problematic_thresholding_subs) == 0
        % remove empty cell array contents 
        problematic_thresholding_subs = problematic_thresholding_subs(~cellfun('isempty',problematic_thresholding_subs));
        % find unique participant ids
        problematic_thresholding_subids = unique([problematic_thresholding_subs{:}]);
        % find the indices of these in the larger sub list and remove
        problematic_thresholding_subidx = find(ismember(high_quality_sc_nkir.sub{timepoint_idx,1},problematic_thresholding_subids));
        for parcellation_idx = 1:length(parcellations)
            thresholded_connectomes_cell{parcellation_idx,1}(problematic_thresholding_subidx,:,:,:) = [];
        end
        high_quality_sc_nkir.sub{timepoint_idx,1}(problematic_thresholding_subidx) = [];
        % now remove associated framewise displacement estimates
        fd_metadata.sc.meanfwd{timepoint_idx,1}(problematic_thresholding_subidx) = [];
        fprintf('removed %d participants with missing communicability matrices at ses-%s.\n', ...
            length(problematic_thresholding_subidx),timepoints_list(timepoint_idx));
    end
    % once we've processed all connectomes for this time point, assign to
    % the output structure with the appropriate parcellation names
    for parcellation_idx = 1:length(parcellations)
        thresholded.sc.(timepoints_list(timepoint_idx)).(parcellations(parcellation_idx)).individual = thresholded_connectomes_cell{parcellation_idx,1};
        thresholded.sc.(timepoints_list(timepoint_idx)).sub = cellstr(high_quality_sc_nkir.sub{timepoint_idx,1});
    end
    % finally, add mean frame-wise displacement estimates
    thresholded.sc.(timepoints_list(timepoint_idx)).meanfwd = fd_metadata.sc.meanfwd{timepoint_idx,1};
    clear thresholded_connectomes_array thresholded_connectomes_cell
    fprintf('thresholding complete.\n');
end
% conduct distance-dependent thresholding to generate a
% group-representative connectome on the thresholded connectomes
for parcellation_idx = 1:length(parcellations)
    % load coordinates for the parcellation
    if parcellation_idx == 1
        schaefer100_metadata = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer100x7_1mm_info.mat');
        coordinates = [schaefer100_metadata.schaefer100x7_1mm_info.x_mni, schaefer100_metadata.schaefer100x7_1mm_info.y_mni,...
            schaefer100_metadata.schaefer100x7_1mm_info.z_mni];
    elseif parcellation_idx == 2
        schaefer200_metadata = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer200x7_1mm_info.mat');
        coordinates = [schaefer200_metadata.schaefer200x7_1mm_info.x_mni, schaefer200_metadata.schaefer200x7_1mm_info.y_mni,...
            schaefer200_metadata.schaefer200x7_1mm_info.z_mni];
    else
        brainnetome246_metadata = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/brainnetome246_info.mat');
        brainnetome246_metadata = brainnetome246_metadata.brainnetome246;
        coordinates = [brainnetome246_metadata.x_mni, brainnetome246_metadata.y_mni, brainnetome246_metadata.z_mni];
    end
    % euclidean distance between the coordinates.
    D = squareform(pdist(coordinates));
    nroi = size(coordinates,1);
    % set the hemisphere ids (1 = left hemisphere, 2 = right)
    hemiid = zeros(nroi,1);
    hemiid(1:nroi/2,1) = 1;
    hemiid(nroi/2:end,1) = 2;
    % concatenate the thresholded connectomes across timepoints
    concatenated_thresholded_sc = cat(1,thresholded.sc.bas1.(parcellations(parcellation_idx)).individual(:,:,:,1),...
        thresholded.sc.flu1.(parcellations(parcellation_idx)).individual(:,:,:,1),...
        thresholded.sc.flu2.(parcellations(parcellation_idx)).individual(:,:,:,1));
    % transpose the connectomes concatenated across time to nroi x nroi
    % x nsub
    concatenated_thresholded_sc = permute(concatenated_thresholded_sc,[2 3 1 ]);
    % calculate the group distance-dependent connectome
    [group_distance_sc, ~] = fcn_group_bins(concatenated_thresholded_sc,D,hemiid,100);
    % find the average connectivity
    average_connectivity = squeeze(mean(concatenated_thresholded_sc,3));
    % replace '1' values in the binary distance-dependent connectome with
    % the corresponding average connectivity values
    idx = find(group_distance_sc);
    group_distance_sc(idx) = average_connectivity(idx);
    % now convert to weighted communicability. 
    % find the strength matrix
    s = diag(sum(group_distance_sc,2));
    % normalise matrix by strength
    pow = (s^-.5)*group_distance_sc*(s^-.5);
    % take the normalised weighted communicability
    c = expm(pow);
    % assign normalised weighted communicability array to output
    thresholded.sc.group.(parcellations(parcellation_idx)) = c;
    fprintf('calculated group-representative communicability matrix for %s structural connectomes.\n', ...
        parcellations(parcellation_idx));
end
save('/imaging/astle/am10/diffusion.map.embedding/data/nki/connectomes/thresholded_structural_and_functional_connectomes.mat','thresholded','-v7.3');
