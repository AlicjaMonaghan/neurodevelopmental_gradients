% this script details selecting high-quality structural and functional 
% connectomes for participants from the baseline and follow-up studies of 
% the centre for attention, learning, and memory (calm). we then harmonise
% across scanner software, to correct for the software upgrade occuring in
% the longitudinal subset of calm. 

%% part 1 - setting up the work space %%
clear;clc;
% set working directory
cd('/imaging/projects/external/nkir/analyses/');
% add path to brain connectivity toolbox (Rubinov and Sporns, 2010)
addpath('/imaging/astle/am10/toolboxes/2019_03_03_BCT/');
% add path to distance-dependent consensus group thresholding toolbox
% (Betzel et al., 2019)
addpath('/imaging/astle/am10/toolboxes/distanceDependent/');
% and path to the combat harmonisation tool!
addpath('/imaging/astle/am10/toolboxes/ComBat/');
% load the structural connectomes from the baseline
calm_baseline_qsiprep = load('/imaging/projects/cbu/calm/baseline_calm_csd_ACT_qsiprep.mat');
calm_baseline_qsiprep = calm_baseline_qsiprep.complete_baseline_calm;
% remove unnecessary parcellations to save space
fields_to_remove = {'aicha384','gordon333','power264','schaefer100x17',...
    'schaefer200x17','schaefer400x17','aal116','schaefer400x7'};
calm_baseline_qsiprep = rmfield(calm_baseline_qsiprep,fields_to_remove);
% now load structural connectomes from the follow-up
calm_followup_qsiprep = load('//cbsu/data/imaging/projects/cbu/calm/calm_2_csd_ACT_qsiprep.mat');
calm_followup_qsiprep = calm_followup_qsiprep.calm_qsiprep;
% and remove unnecessary fields
calm_followup_qsiprep = rmfield(calm_followup_qsiprep,fields_to_remove);
% concatenate the baseline and followup
calm_qsiprep = struct();
calm_qsiprep.baseline = calm_baseline_qsiprep;
calm_qsiprep.followup = calm_followup_qsiprep;
clear calm_baseline_qsiprep calm_followup_qsiprep
% and load the functional connectomes
calm_baseline_fmriprep = load('connectomes/baseline_calm_cleaned_functional_connectomes.mat');
calm_baseline_fmriprep = calm_baseline_fmriprep.calm_fmriprep;
calm_followup_fmriprep = load('connectomes/longitudinal_calm_cleaned_functional_connectomes.mat');
calm_followup_fmriprep = calm_followup_fmriprep.calm_fmriprep;
% create a fc structure
calm_fmriprep = struct();
calm_fmriprep.baseline = calm_baseline_fmriprep;
calm_fmriprep.followup = calm_followup_fmriprep;
clear calm_baseline_fmriprep calm_followup_fmriprep
% set the parcellations
parcellations = ["schaefer100x7","schaefer200x7","brainnetome246"];
% and the number of time points
timepoints_list = {'baseline','followup'};
modalities = {'sc','fc'};
% read the participant meta-data across time points
calm_metadata_baseline = readtable('connectomes/calm/harmonisation/calm_scanner_metadata_for_haromisation_baseline_imputed.csv');
calm_metadata_longitudinal = readtable('connectomes/calm/harmonisation/calm_scanner_metadata_for_haromisation_longitudinal_imputed.csv');
% and read in information about whether the children were referred or not
project72_calm200 = readtable('//cbsu/data/Imaging/projects/external/nkir/data/baseline_calm/Project 72_CALM 200datasheet_171221.xlsx');
project72_calm800 = readtable('//cbsu/data/Imaging/projects/external/nkir/data/baseline_calm/Project 72_CALM 800datasheet_160421.xlsx');
% for baseline participants, add a column to identify which were referred.
% note that all longitudinal participants were referred.
calm_metadata_baseline.referred = ismember(calm_metadata_baseline.ID, project72_calm800.IDNo_);
fprintf('%d out of %d baseline participants with neuroimaging data were referred.\n', ...
    length(find(calm_metadata_baseline.referred==1)),height(calm_metadata_baseline));

%% part 2 - removing repeat scans from baseline calm %%
repeat_scans = ["sub-1019","sub-1013","sub-1011","sub-1007","sub-1009"];
[~,repeat_scans_idx] = intersect(calm_fmriprep.baseline.sample.sub,repeat_scans);
calm_fmriprep.baseline.schaefer100x7.connectivity(repeat_scans_idx,:,:) = [];
calm_fmriprep.baseline.schaefer200x7.connectivity(repeat_scans_idx,:,:) = [];
calm_fmriprep.baseline.brainnetome246.connectivity(repeat_scans_idx,:,:) = [];
calm_fmriprep.baseline.framewise_displacement(repeat_scans_idx,:) = [];
calm_fmriprep.baseline.sample.sub(repeat_scans_idx,:) = [];
% repeat for structural connectomes!
[~,repeat_scans_idx] = intersect(calm_qsiprep.baseline.sample.id.sub,repeat_scans);
calm_qsiprep.baseline.schaefer100x7.connectivity(repeat_scans_idx,:,:,:) = [];
calm_qsiprep.baseline.schaefer200x7.connectivity(repeat_scans_idx,:,:,:) = [];
calm_qsiprep.baseline.brainnetome246.connectivity(repeat_scans_idx,:,:,:) = [];
calm_qsiprep.baseline.qc.framedisplacement(repeat_scans_idx,:) = [];
calm_qsiprep.baseline.sample.id.sub(repeat_scans_idx,:) = [];

%% part 3 - selecting participants with high-quality functional connectomes %%
% for each participant and at each session, calculate mean framewise
% displacement, and find participants for whom over 20% of their spikes are
% high-motion (i.e. framewise displacement > .5). 
for timepoint_idx = 1:length(timepoints_list)
    % find the original number of functional connectomes for this timepoint
    sub_list = calm_fmriprep.(timepoints_list{timepoint_idx}).sample.sub;
    fprintf('for the %s time point of calm, %d participants had functional connectomes.\n', ...
        timepoints_list{timepoint_idx},length(sub_list));
    % initialise output array for mean framewise displacement and
    % proportion of high-motion spikes out of all spikes
    fd_qc = zeros(length(sub_list),2);
    for sub_idx = 1:length(sub_list)
        fd = calm_fmriprep.(timepoints_list{timepoint_idx}).framewise_displacement{sub_idx};
        % calculate mean framewise displacement
        fd_qc(sub_idx,1) = mean(fd);
    end
    % find participants with high mean framewise displacement
    high_mean_fd = find(fd_qc(:,1)>.5);
    % remove these participants from the parcellated connectomes
    for parcellation_idx = 1:length(parcellations)
        calm_fmriprep.(timepoints_list{timepoint_idx}).(parcellations{parcellation_idx}).connectivity(high_mean_fd,:,:) = [];
    end
    % and remove corresponding cells in the framewise displacement array
    % and subject list
    calm_fmriprep.(timepoints_list{timepoint_idx}).framewise_displacement(high_mean_fd) = [];
    sub_list(high_mean_fd) = [];
    % now find participants with a high proportion of high-motion spikes
    for sub_idx = 1:length(sub_list)
        fd = calm_fmriprep.(timepoints_list{timepoint_idx}).framewise_displacement{sub_idx};
        fd_qc(sub_idx,2) = sum(fd>.5)/length(fd);
    end
    % and participants with a high proportion of high-motion spikes
    high_motion_spikes_participants = find(fd_qc(:,2)>.2);
    % and remove from parcellated connectomes
    for parcellation_idx = 1:length(parcellations)
        calm_fmriprep.(timepoints_list{timepoint_idx}).(parcellations{parcellation_idx}).connectivity(high_motion_spikes_participants,:,:) = [];
    end
    % update subject list and framewise displacement 
    sub_list(high_motion_spikes_participants) = [];
    calm_fmriprep.(timepoints_list{timepoint_idx}).framewise_displacement(high_motion_spikes_participants) = [];
    fprintf(['for the %s time point of calm, we removed %d participants with high mean fd, and a' ...
        ' further %d participants with a high proportion of high-motion spikes.\n'],timepoints_list{timepoint_idx},length(high_mean_fd),length(high_motion_spikes_participants))
    % re-calculate mean framewise displacement 
    fd_qc(high_mean_fd,:) = [];
    fd_qc(high_motion_spikes_participants,:) = [];
    fprintf('for the %s time point of calm, %d participants had high-quality functional connectomes, with mean framewise displacement of %f., and standard deviation of %f.\n', ...
        timepoints_list{timepoint_idx},length(sub_list),mean(fd_qc(:,1)),std((fd_qc(:,1))));
    % remove the participants from the sub in the structure
    calm_fmriprep.(timepoints_list{timepoint_idx}).sample.sub(high_mean_fd) = [];
    calm_fmriprep.(timepoints_list{timepoint_idx}).sample.sub(high_motion_spikes_participants) = [];
end
clear baseline_calm_fmriprep_qced followup_calm_fmriprep_qced
%% part 4 - select participants with high-quality structural connectomes %%
for timepoint_idx = 1:length(timepoints_list)
    % find the original number of structural connectomes for this session
    sub_list = calm_qsiprep.(timepoints_list{timepoint_idx}).sample.id.sub;
    fprintf('for the %s time point of calm, %d participants had structural connectomes.\n', ...
        timepoints_list{timepoint_idx},length(sub_list));
    % initialise output array for mean framewise displacement and
    % proportion of high-motion spikes out of all spikes
    fd_qc = zeros(length(sub_list),2);
    for sub_idx = 1:length(sub_list)
        % assign mean framewise displacement
        fd_qc(sub_idx,1) = calm_qsiprep.(timepoints_list{timepoint_idx}).qc.summary(sub_idx,1);
    end
    % find participants with high mean framewise displacement. note that we
    % use a threshold of 3mm here, rather than .5mm in the case of
    % functional connectomes. we're using the
    % sift_invnodevol_radius2_count_connectivity field as our structural
    % connectivity measure (index of 3).
    high_mean_fd = find(fd_qc(:,1)>3);
    % remove these participants from the parcellated connectomes
    for parcellation_idx = 1:length(parcellations)
        calm_qsiprep.(timepoints_list{timepoint_idx}).(parcellations{parcellation_idx}).connectivity(high_mean_fd,:,:,:) = [];
    end
    % and remove corresponding cells in the framewise displacement array
    % and subject list
    calm_qsiprep.(timepoints_list{timepoint_idx}).qc.framedisplacement(high_mean_fd,:) = [];
    sub_list(high_mean_fd) = [];
    fprintf('for the %s time point of calm, we removed %d participants with high mean fd.\n',timepoints_list{timepoint_idx},length(high_mean_fd));
    % re-calculate mean framewise displacement 
    fd_qc(high_mean_fd,:) = [];
    fprintf('for the %s time point of calm, %d participants had high-quality structural connectomes, with mean framewise displacement of %f, and sd of %f.\n', ...
        timepoints_list{timepoint_idx},length(sub_list),mean(fd_qc(:,1)),std(fd_qc(:,1)));
    % remove the participants from the sub in the structure
    calm_qsiprep.(timepoints_list{timepoint_idx}).sample.id.sub(high_mean_fd) = []; 
end
%% part 5 - setting up combat %%
% for baseline and longitudinal scans, extract age at scan, gender, matrix 
% reasoning, phonological reasoning (alliteration and object counting), and
% working memory (forward and backward digit span, and dot matrix). for all
% cognitive measures, we use age-adjusted scores. do this for structural
% and functional connectivity. 
[baseline_participants_indexed_functional,idx_baseline_fc] = intersect(calm_metadata_baseline.BIDS_ID_T1,calm_fmriprep.baseline.sample.sub);
[baseline_participants_indexed_structural,idx_baseline_sc] = intersect(calm_metadata_baseline.BIDS_ID_T1,calm_qsiprep.baseline.sample.id.sub);
demo_vars_baseline = ["Age_Scan1","Matrix_Reasoning_T_Score_for_analysis","PhAb_Alliteration_Standard_Score_For_Analysis",...
    "PhAB_Object_RAN_RT_Standard_Score_for_analysis","AWMA_Digit_Recall_Standard","AWMA_Dot_Matrix_Standard",...
    "AWMA_Backward_Digit_Standard"];
baseline_participants_metadata_fc = table2array(calm_metadata_baseline(idx_baseline_fc,cellstr(demo_vars_baseline)));
baseline_participants_metadata_sc = table2array(calm_metadata_baseline(idx_baseline_sc,cellstr(demo_vars_baseline)));
% add gender separately because we need to dummy code this categorical variable.
baseline_participants_metadata_fc_gender = dummyvar(categorical(calm_metadata_baseline.Gender(idx_baseline_fc)));
baseline_participants_metadata_sc_gender = dummyvar(categorical(calm_metadata_baseline.Gender(idx_baseline_sc)));
% repeat this process for the longitudinal subset.
[longitudinal_participants_indexed_functional,idx_long_fc] = intersect(calm_metadata_longitudinal.BIDS_ID_T2,calm_fmriprep.followup.sample.sub);
[longitudinal_participants_indexed_structural,idx_long_sc] = intersect(calm_metadata_longitudinal.BIDS_ID_T2,calm_qsiprep.followup.sample.id.sub);
demo_vars_longitudinal = ["Age_Scan2","Matrix_Reasoning_T_Score_for_analysis_T2","PhAb_Alliteration_Standard_Score_for_Analysis_T2",...
    "PhAB_Object_RAN_RT_Standard_Score_for_analysis_T2","AWMA_Digit_Recall_Standard_T2","AWMA_Dot_Matrix_Standard_T2",...
    "AWMA_Backward_Digit__Standard_T2"];
longitudinal_participants_metadata_fc = table2array(calm_metadata_longitudinal(idx_long_fc,cellstr(demo_vars_longitudinal)));
longitudinal_participants_metadata_sc = table2array(calm_metadata_longitudinal(idx_long_sc,cellstr(demo_vars_longitudinal)));
longitudinal_participants_metadata_fc_gender = dummyvar(categorical(calm_metadata_longitudinal.Gender(idx_long_fc)));
longitudinal_participants_metadata_sc_gender = dummyvar(categorical(calm_metadata_longitudinal.Gender(idx_long_sc)));
% the major system fault in the scanner occurred on 20th April 2021,
% therefore any data collected after this needs to be treated as if it were
% collected in a different scanning site. therefore, find which
% participants had data collected after this! this only affects
% longitudinal participants.
idx_longitudinal_dates_fc = ("2021-04-20" < calm_metadata_longitudinal.Date_Scan_T2(idx_long_fc));
idx_longitudinal_dates_sc = ("2021-04-20" < calm_metadata_longitudinal.Date_Scan_T2(idx_long_sc));
% find proportion of longitudinal fc and sc scans affected
prop_longitudinal_fc_affected = sum(idx_longitudinal_dates_fc(:)==0)/length(idx_long_fc)*100;
prop_longitudinal_sc_affected = sum(idx_longitudinal_dates_sc(:)==0)/length(idx_long_sc)*100;
% create an array which shows which participants were affected across both
% time points
scanner_id_fc = cat(1,repelem(0,length(idx_baseline_fc))',idx_longitudinal_dates_fc);
scanner_id_sc = cat(1,repelem(0,length(idx_baseline_sc))',idx_longitudinal_dates_sc);
% create a table which will hold the variables required for non-linear
% mixed effect models of age onto manifold eccentricity. these covariates
% will be age, sex, mean framewise displacement, and scanner id. 
meanfwd_baseline_fc = zeros(length(baseline_participants_indexed_functional),1);
for sub = 1:length(baseline_participants_indexed_functional)
    meanfwd_baseline_fc(sub,1) = mean(cell2mat(calm_fmriprep.baseline.framewise_displacement(sub,1)));
end
meanfwd_longitudinal_fc = zeros(length(longitudinal_participants_indexed_functional),1);
for sub = 1:length(longitudinal_participants_indexed_functional)
    meanfwd_longitudinal_fc(sub,1) = mean(cell2mat(calm_fmriprep.followup.framewise_displacement(sub,1)));
end
% add the universal IDs across participants and sessions
[~,baseline_idx] = intersect(calm_metadata_baseline.BIDS_ID_T1,baseline_participants_indexed_functional);
[~,longitudinal_idx] = intersect(calm_metadata_longitudinal.BIDS_ID_T2,longitudinal_participants_indexed_functional);
calm_fc_metadata = table(cat(1,baseline_participants_indexed_functional,longitudinal_participants_indexed_functional), ...
    cat(1,table2array(calm_metadata_baseline(idx_baseline_fc,cellstr("Age_Scan1"))),table2array(calm_metadata_longitudinal(idx_long_fc,cellstr("Age_Scan2")))),...
    cat(1,table2array(calm_metadata_baseline(idx_baseline_fc,cellstr("Gender"))),table2array(calm_metadata_longitudinal(idx_long_fc,cellstr("Gender")))),...
    cat(1,meanfwd_baseline_fc,meanfwd_longitudinal_fc),scanner_id_fc, ...
    cat(1,repelem(0,length(idx_baseline_fc))',repelem(1,length(idx_long_fc))'), ...
    cat(1,calm_metadata_baseline.ID(baseline_idx),calm_metadata_longitudinal.ID(longitudinal_idx)));
calm_fc_metadata.Properties.VariableNames = cellstr(["BIDS","scan_age","Sex","meanFWD","scannerID","timepoint","ID"]);
% repeat for sc!
[~,idx] = intersect(calm_qsiprep.baseline.sample.id.sub,baseline_participants_indexed_structural);
meanfwd_baseline_sc = calm_qsiprep.baseline.qc.summary(idx,1);
[~,idx] = intersect(calm_qsiprep.followup.sample.id.sub,longitudinal_participants_indexed_structural);
meanfwd_longitudinal_sc = calm_qsiprep.followup.qc.summary(idx,1);
[~,baseline_idx] = intersect(calm_metadata_baseline.BIDS_ID_T1,baseline_participants_indexed_structural);
[~,longitudinal_idx] = intersect(calm_metadata_longitudinal.BIDS_ID_T2,longitudinal_participants_indexed_structural);
calm_sc_metadata = table(cat(1,baseline_participants_indexed_structural,longitudinal_participants_indexed_structural), ...
    cat(1,table2array(calm_metadata_baseline(idx_baseline_sc,cellstr("Age_Scan1"))),table2array(calm_metadata_longitudinal(idx_long_sc,cellstr("Age_Scan2")))),...
    cat(1,table2array(calm_metadata_baseline(idx_baseline_sc,cellstr("Gender"))),table2array(calm_metadata_longitudinal(idx_long_sc,cellstr("Gender")))),...
    cat(1,meanfwd_baseline_sc,meanfwd_longitudinal_sc),scanner_id_sc,...
    cat(1,repelem(0,length(idx_baseline_sc))',repelem(1,length(idx_long_sc))'), ...
    cat(1,calm_metadata_baseline.ID(baseline_idx),calm_metadata_longitudinal.ID(longitudinal_idx)));
calm_sc_metadata.Properties.VariableNames = cellstr(["BIDS","scan_age","Sex","meanFWD","scannerID","timepoint","ID"]);
%% part 6 - combat harmonisation! %%
% to load the harmonised connectomes, uncomment the following code:
harmonized_calm_connectomes = load('connectomes/calm/harmonisation/harmonized_connectomes.mat');
harmonized_calm_connectomes = harmonized_calm_connectomes.harmonized_calm_connectomes;
% create a structure to hold the harmonized connectomes for each
% parcellation
harmonized_calm_connectomes = struct();
harmonized_calm_connectomes.baseline.sc = struct();
harmonized_calm_connectomes.baseline.fc = struct();
harmonized_calm_connectomes.followup.sc = struct();
harmonized_calm_connectomes.followup.fc = struct();
for parcellation_idx = 1:length(parcellations)
    parcellation = parcellations(parcellation_idx);
    fprintf('Harmonising for the %s parcellation.\n',parcellation);
    % extract the number of regions of interest (nroi)
    nroi_extracted = extract(parcellation,digitsPattern);
    nroi = str2double(nroi_extracted(1));
    % stack the fc baseline and longitudinal matrices
    [~,idx_baseline_fc] =  intersect(calm_fmriprep.baseline.sample.sub,calm_metadata_baseline.BIDS_ID_T1);
    [~,idx_longitudinal_fc] =  intersect(calm_fmriprep.followup.sample.sub,calm_metadata_longitudinal.BIDS_ID_T2);
    stacked_fc = cat(1,calm_fmriprep.baseline.(parcellation).connectivity(idx_baseline_fc,:,:),calm_fmriprep.followup.(parcellation).connectivity(idx_longitudinal_fc,:,:));
    % check that the number of participants is as expected
    if length(stacked_fc) ~= length(calm_fmriprep.baseline.sample.sub) + length(calm_fmriprep.followup.sample.sub)
        fprintf('incorrect number of participants for fc harmonisation.\n');
    end
    % repeat for sc - we do this for the
    % 'sift_invnodevol_radius2_count_connectivity' sc measure,
    % corresponding to the third index
    [~,idx_baseline_sc] =  intersect(calm_qsiprep.baseline.sample.id.sub,calm_metadata_baseline.BIDS_ID_T1);
    [~,idx_longitudinal_sc] =  intersect(calm_qsiprep.followup.sample.id.sub,calm_metadata_longitudinal.BIDS_ID_T2);
    stacked_sc = cat(1,squeeze(calm_qsiprep.baseline.(parcellation).connectivity(idx_baseline_sc,3,:,:)),...
        squeeze(calm_qsiprep.followup.(parcellation).connectivity(idx_longitudinal_sc,3,:,:)));
    % again, check we have the correct number of participants
    if length(stacked_sc) ~= length(calm_qsiprep.baseline.sample.id.sub) + length(calm_qsiprep.followup.sample.id.sub)
        fprintf('incorrect number of participants for sc harmonisation.\n');
    end
    %% step 6b - turn the matrices into vectors %%
    vectorconnect_fc = zeros(size(stacked_fc,1),(nroi*(nroi-1))/2);
    % make bottom half of matrix into a vector
    for i=1:size(stacked_fc,1)
        values = tril(squeeze(stacked_fc(i,:,:)));
        vectorconnect_fc(i,:)= squareform((values-diag(diag(values)).'));
        clear values
    end
    clear i
    % find indices of columns that are all zero
    zerocols_fc = all(vectorconnect_fc == 0);
    vectorconnect_fc_nnz = vectorconnect_fc;
    vectorconnect_fc_nnz(:,zerocols_fc)=[];
    % repeat for sc
    vectorconnect_sc = zeros(size(stacked_sc,1),(nroi*(nroi-1))/2);
    % make bottom half of matrix into a vector
    for i=1:size(stacked_sc,1)
        values = tril(squeeze(stacked_sc(i,:,:)));
        vectorconnect_sc(i,:)= squareform((values-diag(diag(values)).'));
        clear values
    end
    clear i
    zerocols_sc = all(vectorconnect_sc == 0);
    vectorconnect_sc_nnz = vectorconnect_sc;
    vectorconnect_sc_nnz(:,zerocols_sc)=[];
    %% step 6c - harmonize using combat! %%
    % horizontally concatenate the baseline and longitudinal gender
    % dummy variable
    gender_baseline_longitudinal_fc = cat(1,baseline_participants_metadata_fc_gender,longitudinal_participants_metadata_fc_gender);
    gender_baseline_longitudinal_sc = cat(1,baseline_participants_metadata_sc_gender,longitudinal_participants_metadata_sc_gender);
    % and then the remaining covariates
    formatted_participant_metadata_fc = cat(1,baseline_participants_metadata_fc,longitudinal_participants_metadata_fc);
    formatted_participant_metadata_sc = cat(1,baseline_participants_metadata_sc,longitudinal_participants_metadata_sc);
    % harmonise! omit one category from the categorical variable (gender)
    harmonized_fc = combat(vectorconnect_fc_nnz',scanner_id_fc,[formatted_participant_metadata_fc gender_baseline_longitudinal_fc(:,2)],1);
    fprintf('harmonized functional connectomes for the %s parcellation.\n',parcellations(parcellation_idx));
    harmonized_sc = combat(vectorconnect_sc_nnz',scanner_id_sc,[formatted_participant_metadata_sc gender_baseline_longitudinal_sc(:,2)],1);  
    fprintf('harmonized structural connectomes for the %s parcellation.\n',parcellations(parcellation_idx));
    %% step 6d - convert harmonized vectors to matrices %%
    % insert zeros where they belong. start with fc...
    harmonized_fc_vector = zeros(size(vectorconnect_fc));
    harmonized_fc_vector_data = harmonized_fc;
    for n = 1:length(harmonized_fc_vector)
        if zerocols_fc(n) == 0
            harmonized_fc_vector(:,n) = harmonized_fc_vector_data(1,:);
            harmonized_fc_vector_data = harmonized_fc_vector_data(2:end,:);
        end
    end
    fprintf('replaced all-zero columns for harmonized functional connectomes.\n');
    % repeat for sc...
    harmonized_sc_vector = zeros(size(vectorconnect_sc));
    harmonized_sc_vector_data = harmonized_sc;
    for n = 1:length(harmonized_sc_vector)
        if zerocols_sc(n) == 0
            harmonized_sc_vector(:,n) = harmonized_sc_vector_data(1,:);
            harmonized_sc_vector_data = harmonized_sc_vector_data(2:end,:);
        end
    end
    fprintf('replaced all-zero columns for harmonized structural connectomes.\n');
    % now convert the non-zero harmonised connectomes to matrices, starting
    % with fc...
    harmonized_fc_connectomes = zeros(size(stacked_fc));
    for n = 1:size(harmonized_fc_connectomes,1)
        vectordata = harmonized_fc_vector(n,:);
        matrix = tril(ones(nroi),-1);
        matrix(matrix > 0) = vectordata;
        harmonized_fc_connectomes(n,:,:) = matrix + matrix';
        clear vectordata matrix
    end
    % extract baseline vs longitudinal outputs and assign to output
    % structure
    harmonized_calm_connectomes.baseline.fc.(parcellation) = squeeze(harmonized_fc_connectomes(1:length(baseline_participants_indexed_functional),:,:));
    harmonized_calm_connectomes.followup.fc.(parcellation) = squeeze(harmonized_fc_connectomes(1:length(longitudinal_participants_indexed_functional),:,:));
    % repeat for sc!
    harmonized_sc_connectomes = zeros(size(stacked_sc));
    for n = 1:size(harmonized_sc_connectomes,1)
        vectordata = harmonized_sc_vector(n,:);
        matrix = tril(ones(nroi),-1);
        matrix(matrix > 0) = vectordata;
        harmonized_sc_connectomes(n,:,:) = matrix + matrix';
        clear vectordata matrix
    end
    harmonized_calm_connectomes.baseline.sc.(parcellations(parcellation_idx)) = squeeze(harmonized_sc_connectomes(1:length(baseline_participants_indexed_structural),:,:));
    harmonized_calm_connectomes.followup.sc.(parcellations(parcellation_idx)) = squeeze(harmonized_sc_connectomes(1:length(longitudinal_participants_indexed_structural),:,:));
end
% save the harmonised connectomes
save('connectomes/calm/harmonisation/harmonized_connectomes.mat','harmonized_calm_connectomes','-v7.3');
%% part 7 - checking densities of structural and functional connectomes 
for parcellation_idx = 1:length(parcellations)
    for timepoint_idx = 1:length(timepoints_list)
        for modality_idx = 1:length(modalities)
            harmonized_connectomes = harmonized_calm_connectomes.(timepoints_list{timepoint_idx}).(modalities{modality_idx}).(parcellations(parcellation_idx));
            nsub = size(harmonized_connectomes,1);
            connectome_density = zeros(nsub,1);
            for sub_idx = 1:nsub
                connectome_density(sub_idx,1) = density_und(squeeze(harmonized_connectomes(sub_idx,:,:)));
            end
            % find any connectomes whose density is more than 3 standard
            % deviations less than the mean!
            participants_to_remove = find(connectome_density < mean(connectome_density) - 3*std(connectome_density));
            fprintf('%s connectomes in the %s parcellation at %s: %d have low density.\n', ...
                modalities{modality_idx}, parcellations{parcellation_idx}, timepoints_list{timepoint_idx}, length(participants_to_remove));
        end
    end
end

%% part 8 - thresholding functional connectomes %%
% in line with margulies and colleagues (2016), for each row, retain the
% top 10% absolute strongest connections. we shall do this for both the
% harmonised and unharmonised connectomes. first, create a new structure to 
% hold the thresholded connectomes.
thresholded = struct();
thresholded.fc = struct();
harmonisation_list = ["unharmonised","harmonised"];
for harmonisation_idx = 1:length(harmonisation_list)
    for timepoint_idx = 1:length(timepoints_list)
        % extract subject list
        sub_list = calm_fmriprep.(timepoints_list{timepoint_idx}).sample.sub;
        for parcellation_idx = 1:length(parcellations)
            % create cell to hold the IDs of participants for whom thresholding
            % was unsuccessful i.e. not enough non-zero values
            failed_thresholding = {};
            % extract the unharmonised or harmonised unthresholded 
            % functional connectomes for this time point and parcellation.
            if harmonisation_idx == 1
                unthresholded_connectomes = calm_fmriprep.(timepoints_list{timepoint_idx}).(parcellations(parcellation_idx)).connectivity;
            else
                unthresholded_connectomes = harmonized_calm_connectomes.(timepoints_list{timepoint_idx}).fc.(parcellations(parcellation_idx));
            end            
            % set the number of regions of interest (nroi)
            nroi = size(unthresholded_connectomes,2);
            % find the number of participants (nsub)
            nsub = size(unthresholded_connectomes,1);
            % output arrays for mean framewise displacement and thresholded
            % connectomes
            mean_fwd = zeros(nsub,1);
            thresholded_connectomes = zeros(nsub,nroi,nroi);
            for sub_idx = 1:length(sub_list)
                sub_connectivity = squeeze(unthresholded_connectomes(sub_idx,:,:));
                % for each row, find the 10% absolute strongest
                % connections, and set all others to 0
                for row = 1:nroi
                    [sorted_values,I] = sort(abs(sub_connectivity(row,:)),'descend');
                    % if the number of regions we have does not produce an
                    % integer when dividing by 10 (i.e. 90th percentile),
                    % round to 2 significant figures
                    indices_to_keep = I(1:floor(round(nroi*.10)));
                    thresholded_connectomes(sub_idx,row,indices_to_keep)  = sorted_values(1:floor(round(nroi*.10)));
                    if nnz(squeeze(thresholded_connectomes(sub_idx,row,:))) ~= round(nroi*.10)
                        fprintf("row %d for %s does not have enough non-zero values.\n", ...
                            row, sub_list{sub_idx});
                        % add this ID to the failed_thresholding cell
                        failed_thresholding = [failed_thresholding; sub_list{sub_idx,1}];
                    end
                end
            end
            fprintf('thresholded %s %s functional connectomes in the %s parcellation.\n', ...
                timepoints_list{timepoint_idx},harmonisation_list(harmonisation_idx), ...
                parcellations(parcellation_idx));
        % assign the thresholded connectomes to the output
        thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).(parcellations(parcellation_idx)).individual = thresholded_connectomes;
        fprintf('%d %s participants in the %s parcellation at %s failed functional thresholding.\n', ...
            length(unique(cellstr(failed_thresholding))), harmonisation_list{harmonisation_idx}, parcellations{parcellation_idx}, timepoints_list{timepoint_idx});
        % add the subject IDs. note that the sub IDs for the harmonised
        % and unharmonised data are the same, so it doesn't matter
        % where we pull from 
        thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub = cellstr(calm_fmriprep.(timepoints_list{timepoint_idx}).sample.sub);
        clear thresholded_connectomes
        end
        % to validate the normative range of structural and functional gradient
        % development from nki, split the calm data set into referred and non-
        % referred subsets. note that we only do this for baseline
        % participants, as all longitudinal participants were referred.
        if timepoint_idx == 1
            % find the BIDS IDs of functional connectivity for all referred
            % participants, regardless of QC
            referred_bids = calm_fc_metadata.BIDS(ismember(table2array(calm_fc_metadata(calm_fc_metadata.timepoint == 0,cellstr("ID"))),project72_calm800.IDNo_)==1);
            % find the indices of these BIDS IDs in the thresholded structure
            referred_bids_idx_into_struct = find(ismember(thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub,referred_bids)==1);
            nonreferred_bids_idx_into_struct = find(ismember(thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub,referred_bids)==0);
            % create a field for referred and non-referred participants at
            % baseline. now populate with the functional connectomes across
            % different parcellations.
            for parcellation_idx = 1:length(parcellations)
                thresholded.fc.(harmonisation_list(harmonisation_idx)).baseline.referred.(parcellations(parcellation_idx)).individual = ...
                    thresholded.fc.(harmonisation_list(harmonisation_idx)).baseline.(parcellations(parcellation_idx)).individual(referred_bids_idx_into_struct,:,:);  
                thresholded.fc.(harmonisation_list(harmonisation_idx)).baseline.nonreferred.(parcellations(parcellation_idx)).individual = ...
                    thresholded.fc.(harmonisation_list(harmonisation_idx)).baseline.(parcellations(parcellation_idx)).individual(nonreferred_bids_idx_into_struct,:,:);  
            end
            % and add subject IDs!
            thresholded.fc.(harmonisation_list(harmonisation_idx)).baseline.referred.sub = thresholded.fc.(harmonisation_list(harmonisation_idx)).baseline.sub(referred_bids_idx_into_struct);
            thresholded.fc.(harmonisation_list(harmonisation_idx)).baseline.nonreferred.sub = thresholded.fc.(harmonisation_list(harmonisation_idx)).baseline.sub(nonreferred_bids_idx_into_struct);
            fprintf('%d referred baseline participants had high-quality thresholded functional connectomes.\n',length(referred_bids_idx_into_struct));
        end
    end
    % average fc across participants across time-points to generate a
    % group-representative connectome. for consistency, do this separately
    % for referred and non-referred participants.
    for parcellation_idx = 1:length(parcellations)
        % concatenate the individual (referred) connectomes across time
        % points
        concatenated_referred = cat(1,thresholded.fc.(harmonisation_list(harmonisation_idx)).baseline.referred.(parcellations(parcellation_idx)).individual,...
            thresholded.fc.(harmonisation_list(harmonisation_idx)).followup.(parcellations(parcellation_idx)).individual);
        % average across participants and time points to find the
        % group-representative connectome
        thresholded.fc.(harmonisation_list(harmonisation_idx)).group.referred.(parcellations(parcellation_idx)) = squeeze(mean(concatenated_referred));
        % average across non-referred participants to find the
        % group-representative connectome (again, only at baseline, where
        % we have non-referred participants).
        thresholded.fc.(harmonisation_list(harmonisation_idx)).group.nonreferred.(parcellations(parcellation_idx)) = squeeze(mean(thresholded.fc.(harmonisation_list(harmonisation_idx)).baseline.nonreferred.(parcellations(parcellation_idx)).individual));
    end
end 

%% part 9 - thresholding structural connectomes %%
% to allow a similar density to the functional connectomes, whilst ensuring
% a connected network, for each row, retain the top 10% strongest
% connections. 
thresholded.sc = struct();
referral_status = ["referred","nonreferred"];
for harmonisation_idx = 1:length(harmonisation_list)
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
        for timepoint_idx = 1:length(timepoints_list)
            if harmonisation_idx == 1
                % extract the 'sift_invnodevol_radius2_count_connectivity' 
                % sc measure, corresponding to the third index
                unthresholded_connectomes = squeeze(calm_qsiprep.(timepoints_list{timepoint_idx}).(parcellations(parcellation_idx)).connectivity(:,3,:,:));
            else
                unthresholded_connectomes = harmonized_calm_connectomes.(timepoints_list{timepoint_idx}).sc.(parcellations(parcellation_idx));
            end
            fprintf('thresholding %s %s structural connectomes in the %s parcellation.\n',harmonisation_list(harmonisation_idx), ...
                timepoints_list{timepoint_idx},parcellations(parcellation_idx));
            % set the number of regions of interest (nroi)
            nroi = size(unthresholded_connectomes,2);
            % find the number of participants (nsub)
            nsub = size(unthresholded_connectomes,1);
            % new arrays for thresholded connectomes (first slice of 4th
            % dimension) and communicability (2nd slice of 4th dimension).
            % we retain the thresholded connectomes to test for the
            % effects of harmonisation.
            thresholded_connectomes = zeros(nsub,nroi,nroi,2);
            for sub_idx = 1:nsub
                sub_connectivity = squeeze(unthresholded_connectomes(sub_idx,:,:));
                % for each row, find the 10% absolute strongest
                % connections, and set all others to 0
                for row = 1:nroi
                    [~,I] = sort(sub_connectivity(row,:),'descend');
                    I = I';
                    indices_to_keep = I(1:floor(nroi*.10));
                    thresholded_connectomes(sub_idx,row,indices_to_keep,1)  = sub_connectivity(row,indices_to_keep);
                end
                % now convert to communicability matrices. first find the
                % strength matrix.
                a = squeeze(thresholded_connectomes(sub_idx,:,:,1));
                s = diag(sum(a,2));
                pow = (s^-.5)*a*(s^-.5);
                c = expm(pow);
                % raise warning if there are any missing values
                if isempty(find(ismissing(c), 1)) == 0
                    fprintf('missing communicability values for %s structural connectomes for subject with index of %d in the %s parcellation at %s.\n', ...
                        harmonisation_list(harmonisation_idx), sub_idx, parcellations{parcellation_idx}, timepoints_list{timepoint_idx});
                end
                % assign to output
                thresholded_connectomes(sub_idx,:,:,2) = c;
                clear a s pow c
            end
            % assign to output
            thresholded.sc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).(parcellations(parcellation_idx)).individual = thresholded_connectomes;
            clear thresholded_connectomes
            fprintf('thresholding complete.\n');
            % add sub IDs, either from the harmonised or unharmonised
            % connectomes (both will be the same)
            thresholded.sc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub = cellstr(calm_qsiprep.(timepoints_list{timepoint_idx}).sample.id.sub);
            % when thresholding participants at baseline, separate them
            % into referred and non-referred. this is because non-referred
            % participants were only recruited at baseline. first, find the
            % indices of BIDs IDs for structural connectivity for all
            % referred participants, regardless of QC
            if timepoint_idx == 1
                referred_bids = calm_sc_metadata.BIDS(ismember(table2array(calm_sc_metadata(calm_sc_metadata.timepoint == 0, cellstr("ID"))),project72_calm800.IDNo_) == 1);
                % find the indices of these BIDS IDs in the thresholded structure
                referred_bids_idx_into_struct = find(ismember(thresholded.sc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub,referred_bids)==1);
                nonreferred_bids_idx_into_struct = find(ismember(thresholded.sc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub,referred_bids)==0);
                % create a field for referred and non-referred participants at
                % baseline. now populate with the functional connectomes across
                % different parcellations.
                thresholded.sc.(harmonisation_list(harmonisation_idx)).baseline.referred.(parcellations(parcellation_idx)).individual = ...
                    thresholded.sc.(harmonisation_list(harmonisation_idx)).baseline.(parcellations(parcellation_idx)).individual(referred_bids_idx_into_struct,:,:,:);  
                thresholded.sc.(harmonisation_list(harmonisation_idx)).baseline.nonreferred.(parcellations(parcellation_idx)).individual = ...
                    thresholded.sc.(harmonisation_list(harmonisation_idx)).baseline.(parcellations(parcellation_idx)).individual(nonreferred_bids_idx_into_struct,:,:,:);  
            end
            % and add subject IDs!
            thresholded.sc.(harmonisation_list(harmonisation_idx)).baseline.referred.sub = thresholded.sc.(harmonisation_list(harmonisation_idx)).baseline.sub(referred_bids_idx_into_struct);
            thresholded.sc.(harmonisation_list(harmonisation_idx)).baseline.nonreferred.sub = thresholded.sc.(harmonisation_list(harmonisation_idx)).baseline.sub(nonreferred_bids_idx_into_struct);
        end
        % conduct distance-dependent thresholding to generate a
        % group-representative connectome on the thresholded connectomes
        % euclidean distance between the coordinates.
        D = squareform(pdist(coordinates));
        nroi = size(coordinates,1);
        % set the hemisphere ids (1 = left hemisphere, 2 = right)
        hemiid = zeros(nroi,1);
        hemiid(1:nroi/2,1) = 1;
        hemiid(nroi/2:end,1) = 2;
        % concatenate the thresholded connectomes across timepoints. do
        % this separately for referred (baseline + longitudinal) and
        % non-referred (baseline only) participants
        for referral_status_idx = 1:length(referral_status)
            if referral_status_idx == 1
                concatenated_thresholded_sc = vertcat(squeeze(thresholded.sc.(harmonisation_list(harmonisation_idx)).baseline.referred.(parcellations(parcellation_idx)).individual(:,:,:,1)), ...
                    squeeze(thresholded.sc.(harmonisation_list(harmonisation_idx)).followup.(parcellations(parcellation_idx)).individual(:,:,:,1)));
            else
                concatenated_thresholded_sc = squeeze(thresholded.sc.(harmonisation_list(harmonisation_idx)).baseline.nonreferred.(parcellations(parcellation_idx)).individual(:,:,:,1));
            end
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
            thresholded.sc.(harmonisation_list(harmonisation_idx)).group.(referral_status(referral_status_idx)).(parcellations(parcellation_idx)) = c;
            fprintf('calculated group-representative communicability matrix for %s %s structural connectomes in the %s parcellation.\n', ...
                harmonisation_list(harmonisation_idx),referral_status(referral_status_idx),parcellations(parcellation_idx));
        end
    end
end
% no missing communicability values for harmonised connectomes!
save('/imaging/astle/am10/diffusion.map.embedding/data/calm/connectomes/thresholded_structural_and_functional_connectomes.mat','thresholded','-v7.3');

%% part 8 - test effectiveness of harmonisation
% for each modality, load the unharmonized and harmonized thresholded 
% connectomes across both time points (baseline and longitudinal), and 
% calculate global efficiency, assortativity, density, transitivity, and
% maximised modularity coefficient. 
graph_theory_metrics_list = ["global efficiency","assortativity","density","transitivity","maximised modularity coefficient"];
% hold the graph theory metrics in an nmodality (structural, functional) *
% nparcellation cell array
graph_theory_metrics_all_modalities = cell(2,length(parcellations));
meta_data_all_modalities = cell(2,1);
modality_list = ["sc","fc"];
harmonisation_list = ["unharmonised","harmonised"];
for parcellation_idx = 1:length(parcellations)
    parcellation = parcellations(parcellation_idx);
    for modality_idx = 1:length(modality_list)
        tic;
        if modality_idx == 1
            % extract the unharmonised connectomes for this modality
            [~,idx_baseline] = intersect(calm_qsiprep.baseline.sample.id.sub,calm_metadata_baseline.BIDS_ID_T1);
            [~,idx_longitudinal] = intersect(calm_qsiprep.followup.sample.id.sub,calm_metadata_longitudinal.BIDS_ID_T2);
            formatted_modality_name = "structural";
        else
            [~,idx_baseline] = intersect(calm_fmriprep.baseline.sample.sub,calm_metadata_baseline.BIDS_ID_T1);
            [~,idx_longitudinal] = intersect(calm_fmriprep.followup.sample.sub,calm_metadata_longitudinal.BIDS_ID_T2);
            formatted_modality_name = "functional";
        end
        stacked_unharmonized_connectomes = cat(1,thresholded.(modality_list(modality_idx)).unharmonised.baseline.(parcellation).individual(idx_baseline,:,:,1),...
            thresholded.(modality_list(modality_idx)).unharmonised.followup.(parcellation).individual(idx_longitudinal,:,:,1));
        % extract the harmonised connectomes 
        stacked_harmonized_connectomes = cat(1,thresholded.(modality_list(modality_idx)).harmonised.baseline.(parcellation).individual(idx_baseline,:,:,1),...
            thresholded.(modality_list(modality_idx)).harmonised.followup.(parcellation).individual(idx_longitudinal,:,:,1));
        loading_time = toc;
        fprintf('loaded the harmonised and unharmonised %s connectomes for the %s parcellation in %d seconds.\n',formatted_modality_name,parcellation,loading_time);
        % extract the participant metadata 
        if modality_idx == 1
            metadata_without_gender = cat(1,baseline_participants_metadata_sc,longitudinal_participants_metadata_sc);
            gender = cat(1,baseline_participants_metadata_sc_gender,longitudinal_participants_metadata_sc_gender);
            % merge gender with the remaining meta data. we select the first
            % column of gender to provide a reference. also add scanner id. 
            metadata = cat(2,metadata_without_gender,gender(:,1),scanner_id_sc);
        else
            metadata_without_gender = cat(1,baseline_participants_metadata_fc,longitudinal_participants_metadata_fc);
            gender = cat(1,baseline_participants_metadata_fc_gender,longitudinal_participants_metadata_fc_gender);
            metadata = cat(2,metadata_without_gender,gender(:,1),scanner_id_fc);
        end
        % add meta-data to the output cell
        meta_data_all_modalities{modality_idx,1} = metadata;
        clear metadata_without_gender gender loading_time
        % find the number of participants
        nsub = length(metadata);
        % initialise an output array for the graph theory metrics for the
        % harmonised and unharmonised connectomes 
        graph_theory_metrics = zeros(nsub,length(graph_theory_metrics_list),2);
        for harmonisation_idx = 1:length(harmonisation_list)
            tic;
            if harmonisation_idx == 1
                connectomes = stacked_unharmonized_connectomes;
            else
                connectomes = stacked_harmonized_connectomes;
            end
            for sub = 1:nsub
                % for each participant, calculate global efficiency
                a = squeeze(connectomes(sub,:,:));
                graph_theory_metrics(sub,1,harmonisation_idx) = efficiency_wei(a,0);
                % assortativity
                graph_theory_metrics(sub,2,harmonisation_idx) = assortativity_wei(a,1);
                % density
                graph_theory_metrics(sub,3,harmonisation_idx) = density_und(a);
                % transitivity
                graph_theory_metrics(sub,4,harmonisation_idx) = transitivity_wu(a);
                % maximised modularity coefficient
                [~, graph_theory_metrics(sub,5,harmonisation_idx)] = modularity_und(a);
            end
            graph_theory_time = toc;
            fprintf('calculated all graph theory metrics for %s %s connectomes in the %s parcellation in %d seconds.\n', ...
                harmonisation_list(harmonisation_idx),formatted_modality_name,parcellation,graph_theory_time);
            clear tic graph_theory time 
        end
        % add to the output cell array
        graph_theory_metrics_all_modalities{modality_idx,parcellation_idx} = real(graph_theory_metrics);
        clear graph_theory_metrics
    end
end
save('connectomes/calm/harmonisation/graph_theory_metrics_before_and_after_harmonisation.mat','graph_theory_metrics_all_modalities');
save('connectomes/calm/harmonisation/meta_data_all_modalities_harmonisation.mat','meta_data_all_modalities');
% for each modality, conduct a general linear model with 8 covariates and
% scanner id as the predictors of graph theory metrics. if the
% harmonisation is successful, we'd expect the scanner id variable to be
% significant in the unharmonised data, but not significant in the
% harmonised data. 
graph_theory_metrics_all_modalities = load('connectomes/calm/harmonisation/graph_theory_metrics_before_and_after_harmonisation.mat');
graph_theory_metrics_all_modalities = graph_theory_metrics_all_modalities.graph_theory_metrics_all_modalities;
meta_data_all_modalities = load('connectomes/calm/harmonisation/meta_data_all_modalities_harmonisation.mat');
meta_data_all_modalities = meta_data_all_modalities.meta_data_all_modalities;
% for speed, only assess harmonisation effectiveness for the parcellation
% used for the main analysis (schaefer 200-node 7-network)
schaefer200_parcellation_idx = strmatch('schaefer200x7',parcellations);
for modality_idx = 1:length(modality_list)
    for metric = 1:length(graph_theory_metrics_list)
        % create the glm, setting sex and scanner id (columns 8 and 9) as
        % categorical.
        unharmonised_mod = fitglm(meta_data_all_modalities{modality_idx,1},graph_theory_metrics_all_modalities{modality_idx,schaefer200_parcellation_idx}(:,metric,1),'CategoricalVars',[8,9]);
        harmonised_mod = fitglm(meta_data_all_modalities{modality_idx,1},graph_theory_metrics_all_modalities{modality_idx,schaefer200_parcellation_idx}(:,metric,2),'CategoricalVars',[8,9]);
        fprintf(['Before harmonisation for %s connectomes,scanner id had beta of %d and p-value of %d in a model predicting %s, ' ...
            'compared to beta of %d and p-value of %d after harmonisation.\n'], modality_list(modality_idx), ...
            round(unharmonised_mod.Coefficients.Estimate(10),3,"significant"),... 
            round(unharmonised_mod.Coefficients.pValue(10),3,"significant"), ...
            graph_theory_metrics_list(metric),...
            round(harmonised_mod.Coefficients.Estimate(10),3,"significant"), ...
            round(harmonised_mod.Coefficients.pValue(10),3,"significant"));
    end
end