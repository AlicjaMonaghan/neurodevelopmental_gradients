% this script pulls cleaned functional connectomes from baseline and
% follow-up sessions from two data sets: the centre for attention,
% learning, and memory (calm), and the nathan kline institute longitudinal
% sample discovery of brain development trajectories. the outputs are from
% the python function clean_and_parcellate_fmri_calm and
% clean_and_parcellate_fmri_nki, respectively. we pull calm connectomes 
% first, and then nki.written by alicja monaghan, 04/2023.

%% part 1 - setting up the work space to pull calm connectomes 
clear;clc;
calm_datasets = {'baseline_calm','longitudinal_calm'};
%% part 2 - load and organise external data to be saved
% parcellation metadata
aal116 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/aal116_info.mat');
aal116 = aal116.aal116;
brainnetome246 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/brainnetome246_info.mat');
brainnetome246 = brainnetome246.brainnetome246;
schaefer100x7 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer100x7_1mm_info.mat');
schaefer100x7 = schaefer100x7.schaefer100x7_1mm_info;
schaefer200x7 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer200x7_1mm_info.mat');
schaefer200x7 = schaefer200x7.schaefer200x7_1mm_info;
schaefer400x7 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer400x7_1mm_info.mat');
schaefer400x7 = schaefer400x7.schaefer400x7_1mm_info;
% organise data for later allocation
% aal116
aal116_coordinates = [aal116.x_mni aal116.y_mni aal116.z_mni];
aal116_hemi = string(aal116.hemi);
aal116_lobe = string(aal116.lobe);
aal116_region_labels = string(aal116.name);
% brainnetome246
brainnetome246_coordinates = [brainnetome246.x_mni brainnetome246.y_mni brainnetome246.z_mni];
brainnetome246_hemi = string(brainnetome246.hemi);
brainnetome246_lobe = string(brainnetome246.lobe);
brainnetome246_region_labels = string(brainnetome246.name);
% schaefer100x7
schaefer100x7_coordinates = [schaefer100x7.x_mni schaefer100x7.y_mni schaefer100x7.z_mni];
schaefer100x7_region_labels = string(schaefer100x7.name);
% schaefer200x7
schaefer200x7_coordinates = [schaefer200x7.x_mni schaefer200x7.y_mni schaefer200x7.z_mni];
schaefer200x7_region_labels = string(schaefer200x7.name);
% schaefer400x7
schaefer400x7_coordinates = [schaefer400x7.x_mni schaefer400x7.y_mni schaefer400x7.z_mni];
schaefer400x7_region_labels = string(schaefer400x7.name);
%% part 3 - collect calm functional connectomes and framewise displacement!
for dataset_idx = 1:length(calm_datasets)
    % set the fmriprep output path for the calm data set we are processing
    if dataset_idx == 1
        fmriprep_dir = '/imaging/projects/external/nkir/analyses/calm-I/fmriprep/outdir/';
    else
        fmriprep_dir = '/imaging/projects/external/nkir/analyses/calm-II/fmriprep/outdir/';
    end
    % find the processed participants
    files = dir(strcat(fmriprep_dir,'final_output/cleaned_connectomes/sub-*_cleaned_functional_connectomes.h5'));
    processed_sub = cell(length(files),1);
    nsub = length(files);
    % initialise an output cell for framewise displacement across all
    % participants
    framewise_displacement_cell = cell(nsub,1);
    % initialise output arrays for functional connectomes for the different
    % parcellations
    schaefer100x7 = zeros(nsub,100,100);
    aal116 = zeros(nsub,116,116);
    schaefer200x7 = zeros(nsub,200,200);
    brainnetome246 = zeros(nsub,246,246);
    schaefer400x7 = zeros(nsub,400,400);
    % loop through each participant 
    for i = 1:length(files)
        % get the sub id
        sub = extractBefore(files(i).name,'_cleaned');
        % save the sub id 
        processed_sub{i,1} = sub;
        % load the confound file
        confounds = readtable(sprintf('%sfinal_output/%s/%s_task-rest_desc-confounds_timeseries.tsv', ...
            fmriprep_dir,sub,sub),'FileType','delimitedtext');
        % extract framewise displacement and convert nan to 0
        framewise_displacement = confounds.framewise_displacement;
        framewise_displacement(isnan(framewise_displacement)) = [];
        % assign to output cell
        framewise_displacement_cell{i,1} = framewise_displacement;
        % load the cleaned functional connectomes, and assign to output
        % arrays
        schaefer100x7(i,:,:) = h5read(sprintf('%sfinal_output/cleaned_connectomes/%s_cleaned_functional_connectomes.h5',fmriprep_dir,sub),'/schaefer100x7');
        aal116(i,:,:) = h5read(sprintf('%sfinal_output/cleaned_connectomes/%s_cleaned_functional_connectomes.h5',fmriprep_dir,sub),'/aal116');
        schaefer200x7(i,:,:) = h5read(sprintf('%sfinal_output/cleaned_connectomes/%s_cleaned_functional_connectomes.h5',fmriprep_dir,sub),'/schaefer200x7');
        brainnetome246(i,:,:) = h5read(sprintf('%sfinal_output/cleaned_connectomes/%s_cleaned_functional_connectomes.h5',fmriprep_dir,sub),'/brainnetome246');
        schaefer400x7(i,:,:) = h5read(sprintf('%sfinal_output/cleaned_connectomes/%s_cleaned_functional_connectomes.h5',fmriprep_dir,sub),'/schaefer400x7');
        fprintf('loaded data for %s from %s.\n',sub,calm_datasets{dataset_idx});
    end
    % assign the outputs to a structure
    calm_fmriprep = struct();
    calm_fmriprep.sample.sub = processed_sub;
    % quality control
    calm_fmriprep.framewise_displacement = framewise_displacement_cell;
    % functional connectomes
    % aal116
    calm_fmriprep.aal116.coordinates = aal116_coordinates;
    calm_fmriprep.aal116.connectivity = aal116;
    calm_fmriprep.aal116.hemi = aal116_hemi;
    calm_fmriprep.aal116.lobe = aal116_lobe;
    calm_fmriprep.aal116.regionlabels = aal116_region_labels;
    % brainnetome246
    calm_fmriprep.brainnetome246.coordinates = brainnetome246_coordinates;
    calm_fmriprep.brainnetome246.connectivity = brainnetome246;
    calm_fmriprep.brainnetome246.hemi = brainnetome246_hemi;
    calm_fmriprep.brainnetome246.lobe = brainnetome246_lobe;
    calm_fmriprep.brainnetome246.regionlabels = brainnetome246_region_labels;
    % schaefer100x7
    calm_fmriprep.schaefer100x7.coordinates = schaefer100x7_coordinates;
    calm_fmriprep.schaefer100x7.connectivity = schaefer100x7;
    calm_fmriprep.schaefer100x7.regionlabels = schaefer100x7_region_labels;
    % schaefer200x7
    calm_fmriprep.schaefer200x7.coordinates = schaefer200x7_coordinates;
    calm_fmriprep.schaefer200x7.connectivity = schaefer200x7;
    calm_fmriprep.schaefer200x7.regionlabels = schaefer200x7_region_labels;
    % schaefer400x7
    calm_fmriprep.schaefer400x7.coordinates = schaefer400x7_coordinates;
    calm_fmriprep.schaefer400x7.connectivity = schaefer400x7;
    calm_fmriprep.schaefer400x7.regionlabels = schaefer400x7_region_labels;
    % and save the structure
    save(sprintf('/imaging/projects/external/nkir/analyses/connectomes/%s_cleaned_functional_connectomes.mat',calm_datasets{dataset_idx}),'calm_fmriprep','-v7.3');
end
%% part 4 - initialise outputs for nki functional connectomes and framewise displacement
% enhanced nki has 3 possible sessions: ses-BAS1, ses-FLU1, and ses-FLU2.
% some participants have all sessions, or only some. therefore, we need to
% find which participants have which processed sessions.
files = dir('/imaging/projects/external/nkir/analyses/enhanced_nkir/fmriprep/outdir/final_output/*_cleaned_functional_connectomes.h5');
nfiles = length(files);
% initialise arrays to hold the subject IDs for each session type.
ses_bas1_sub = {};
ses_flu1_sub = {};
ses_flu2_sub = {};
% loop through each file, extract the session type, and append that
% participant's ID to the appropriate list.
all_participant_ids = {};
for file_idx = 1:nfiles
    % extract the sub id and session id 
    sub_id = files(file_idx).name(1:13);
    ses_id = files(file_idx).name(15:22);
    % append to the appropriate lists
    if matches(ses_id,'ses-BAS1') == 1
        ses_bas1_sub{end+1} = sub_id;
    elseif matches(ses_id, 'ses-FLU1') == 1
        ses_flu1_sub{end+1} = sub_id;
    else 
        ses_flu2_sub{end+1} = sub_id;
    end
    all_participant_ids{end+1} = sub_id;
end
% find unique sub ids
unique_sub_ids = unique(all_participant_ids)';
% find how many participants have how many processed sessions.
bas1_nsub = length(ses_bas1_sub);
flu1_nsub = length(ses_flu1_sub);
flu2_nsub = length(ses_flu2_sub);

% initialise output arrays and cells - one for each session.
framedisplace = struct();
framedisplace.bas1 = cell(bas1_nsub,1);
framedisplace.flu1 = cell(flu1_nsub,1);
framedisplace.flu2 = cell(flu2_nsub,1);
fns = fieldnames(framedisplace);
% schaefer100x7
schaefer100x7 = struct();
schaefer100x7.bas1 = zeros(bas1_nsub,100,100);
schaefer100x7.flu1 = zeros(flu1_nsub,100,100);
schaefer100x7.flu2 = zeros(flu2_nsub,100,100);
% aal116
aal116 = struct();
aal116.bas1 = zeros(bas1_nsub,116,116);
aal116.flu1 = zeros(flu1_nsub,116,116);
aal116.flu2 = zeros(flu2_nsub,116,116);
% schaefer200x7
schaefer200x7 = struct();
schaefer200x7.bas1 = zeros(bas1_nsub,200,200);
schaefer200x7.flu1 = zeros(flu1_nsub,200,200);
schaefer200x7.flu2 = zeros(flu2_nsub,200,200);
% brainectome246
brainnetome246 = struct();
brainnetome246.bas1 = zeros(bas1_nsub,246,246);
brainnetome246.flu1 = zeros(flu1_nsub,246,246);
brainnetome246.flu2 = zeros(flu2_nsub,246,246);
% schaefer400x7
schaefer400x7 = struct();
schaefer400x7.bas1 = zeros(bas1_nsub,400,400);
schaefer400x7.flu1 = zeros(flu1_nsub,400,400);
schaefer400x7.flu2 = zeros(flu2_nsub,400,400);
%% part 5 - collect nki functional connectomes across sessions!
% set the new fmriprep_dir
fmriprep_dir = '/imaging/projects/external/nkir/analyses/enhanced_nkir/fmriprep/outdir/final_output/';
for i = 1:nfiles
    tic;
    % find the subject ID
    sub = files(i).name(1:13);
    % extract the session
    session = files(i).name(15:22);
    % find the index of this participant in the relevant session cell array
    if session == 'ses-BAS1'
        subi = find(strcmp(ses_bas1_sub,sub));
        sessioni = 1;
    elseif session == 'ses-FLU1'
        subi = find(strcmp(ses_flu1_sub,sub));
        sessioni = 2;
    else
        subi = find(strcmp(ses_flu2_sub,sub));
        sessioni = 3;
    end
    % display
    fprintf('loading %s data from %s...\n',sub,session);
    % load the confound file
    confounds = readtable(sprintf('%s%s/%s/%s_%s_task-rest_acq-1400_desc-confounds_timeseries.tsv', ...
        fmriprep_dir,sub,session,sub,session),'FileType','delimitedtext');
    % extract framewise displacement and convert nan to 0
    framewise_displacement = confounds.framewise_displacement;
    framewise_displacement(isnan(framewise_displacement)) = [];
    % strip session to just include the session name, without ses- prefix.
    session = lower(extractAfter(session,'ses-'));
    % assign to output cell
    framedisplace.(fns{sessioni}){subi,1} = framewise_displacement;
    % load the cleaned functional connectomes, and assign to output
    % arrays
    fc_file_name = sprintf('%s%s_ses-%s_cleaned_functional_connectomes.h5',fmriprep_dir,sub,upper(session));
    schaefer100x7.(session)(subi,:,:) = h5read(fc_file_name,'/schaefer100x7');
    aal116.(session)(subi,:,:) = h5read(fc_file_name,'/aal116');
    schaefer200x7.(session)(subi,:,:) = h5read(fc_file_name,'/schaefer200x7');
    brainnetome246.(session)(subi,:,:) = h5read(fc_file_name,'/brainnetome246');
    schaefer400x7.(session)(subi,:,:) = h5read(fc_file_name,'/schaefer400x7');
    fprintf('loaded data for %s from ses-%s.\n',sub,session);
end
% assign to a struct
nki_fmriprep = struct();
% sub ids
nki_fmriprep.sample.sub.bas1 = ses_bas1_sub';
nki_fmriprep.sample.sub.flu1 = ses_flu1_sub';
nki_fmriprep.sample.sub.flu2 = ses_flu2_sub';
% framewise displacement
nki_fmriprep.framewisedisplacement.bas1 = framedisplace.bas1;
nki_fmriprep.framewisedisplacement.flu1 = framedisplace.flu1;
nki_fmriprep.framewisedisplacement.flu2 = framedisplace.flu2;
% functional connectomes
% schaefer100x7
nki_fmriprep.schaefer100x7.coordinates = schaefer100x7_coordinates;
nki_fmriprep.schaefer100x7.connectivity.bas1 = schaefer100x7.bas1;
nki_fmriprep.schaefer100x7.connectivity.flu1 = schaefer100x7.flu1;
nki_fmriprep.schaefer100x7.connectivity.flu2 = schaefer100x7.flu2;
nki_fmriprep.schaefer100x7.regionlabels = schaefer100x7_region_labels;
% aal116
nki_fmriprep.aal116.coordinates = aal116_coordinates;
nki_fmriprep.aal116.connectivity.bas1 = aal116.bas1;
nki_fmriprep.aal116.connectivity.flu1 = aal116.flu1;
nki_fmriprep.aal116.connectivity.flu2 = aal116.flu2;
nki_fmriprep.aal116.regionlabels = aal116_region_labels;
nki_fmriprep.aal116.hemi = aal116_hemi;
nki_fmriprep.aal116.lobe = aal116_lobe;
% schaefer200x7
nki_fmriprep.schaefer200x7.coordinates = schaefer200x7_coordinates;
nki_fmriprep.schaefer200x7.connectivity.bas1 = schaefer200x7.bas1;
nki_fmriprep.schaefer200x7.connectivity.flu1 = schaefer200x7.flu1;
nki_fmriprep.schaefer200x7.connectivity.flu2 = schaefer200x7.flu2;
nki_fmriprep.schaefer200x7.regionlabels = schaefer200x7_region_labels;
% brainnetome246
nki_fmriprep.brainnetome246.coordinates = brainnetome246_coordinates;
nki_fmriprep.brainnetome246.connectivity.bas1 = brainnetome246.bas1;
nki_fmriprep.brainnetome246.connectivity.flu1 = brainnetome246.flu1;
nki_fmriprep.brainnetome246.connectivity.flu2 = brainnetome246.flu2;
nki_fmriprep.brainnetome246.regionlabels = brainnetome246_region_labels;
nki_fmriprep.brainnetome246.hemi = brainnetome246_hemi;
nki_fmriprep.brainnetome246.lobe = brainnetome246_lobe;
% schaefer200x7
nki_fmriprep.schaefer400x7.coordinates = schaefer400x7_coordinates;
nki_fmriprep.schaefer400x7.connectivity.bas1 = schaefer400x7.bas1;
nki_fmriprep.schaefer400x7.connectivity.flu1 = schaefer400x7.flu1;
nki_fmriprep.schaefer400x7.connectivity.flu2 = schaefer400x7.flu2;
nki_fmriprep.schaefer400x7.regionlabels = schaefer400x7_region_labels;
% save the structure!
save('/imaging/projects/external/nkir/analyses/connectomes/longitudinal_nki/cleaned_and_parcellated_functional_connectomes.mat','nki_fmriprep','-v7.3');

    
