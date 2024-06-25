% this script details selecting high-quality structural and functional 
% connectomes for participants from the baseline and follow-up studies of 
% the centre for attention, learning, and memory (calm). we then harmonise
% across scanner software, to correct for the software upgrade occuring in
% the longitudinal subset of calm. 

%% part 1 - setting up the work space %%
clear;clc;
% set working directory
cd('//cbsu/data/imaging/projects/external/nkir/analyses/');
% add path to brain connectivity toolbox (Rubinov and Sporns, 2010)
%addpath('/imaging/astle/am10/toolboxes/2019_03_03_BCT/'); %I don't have permission to add this
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT/');
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
%calm_followup_qsiprep = load('calm_2_csd_ACT_qsiprep.mat');
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

%% part 5 - combat harmonisation! %%
% to load the harmonised connectomes, uncomment the following code:
%harmonized_calm_connectomes = load('connectomes/calm/harmonisation/harmonized_connectomes.mat');
%harmonized_calm_connectomes = harmonized_calm_connectomes.harmonized_calm_connectomes;
% create a structure to hold the harmonized connectomes for each
% parcellation
harmonized_calm_connectomes = struct(); %this creates a slightly adjusted structure for the connectomes
harmonized_calm_connectomes.baseline.sc = struct();
harmonized_calm_connectomes.control.sc = struct();
harmonized_calm_connectomes.followup.sc = struct();
harmonized_calm_connectomes.baseline.fc = struct();
harmonized_calm_connectomes.control.fc = struct();
harmonized_calm_connectomes.followup.fc = struct();

%Generate a list of bids IDs, to double check we have the order right
all_sc_ids=[calm_qsiprep.baseline.sample.id.sub;calm_qsiprep.followup.sample.id.sub];
all_fc_ids=[calm_fmriprep.baseline.sample.sub;calm_fmriprep.followup.sample.sub];

%bring the mega/cog/behavioural data in here
Meta_SC=readtable('Alicja_calm_structural_master.xlsx');
Meta_FC=readtable('Alicja_calm_functional_master.xlsx');
%now define somethings we'll need for the harmonisation
ScannerID_SC=squeeze(table2array(Meta_SC(:,5)));
Covariates_SC=table2array(Meta_SC(:,[2,3,6,7,9,11:16]));%So this should take age at scan, sex, timepoint, referred, age at test, and the cog data
ScannerID_FC=squeeze(table2array(Meta_FC(:,5)));
Covariates_FC=table2array(Meta_FC(:,[2,3,6,7,9,11:16]));
Group_SC=squeeze(table2array(Meta_SC(:,6)))+squeeze(table2array(Meta_SC(:,7))); 
Group_FC=squeeze(table2array(Meta_FC(:,6)))+squeeze(table2array(Meta_FC(:,7))); %makes a variable in which 0=comparison, 1=referred baseline, 2=referred followup

%The first time I ran it, I did a tiny bit of imputation, because
%covariates need to be compelte for combat. These imputed data are now
%included in the master sheets
Covariates_SC=knnimpute(Covariates_SC',10);
Covariates_SC=Covariates_SC';
Covariates_FC=knnimpute(Covariates_FC',10);
Covariates_FC=Covariates_FC';





for parcellation_idx = 1:length(parcellations)
    parcellation = parcellations(parcellation_idx);
    fprintf('Harmonising for the %s parcellation.\n',parcellation);
    % extract the number of regions of interest (nroi)
    %nroi_extracted = extract(parcellation,digitsPattern);
    nroi = size(squeeze(calm_qsiprep.baseline.(parcellation).connectivity(:,3,:,:)),3);
    %nroi = str2double(nroi_extracted(1));
    
        %% step 6b - turn the matrices into vectors %%
    %%%bring in the connectomes here
    stacked_sc = cat(1,squeeze(calm_qsiprep.baseline.(parcellation).connectivity(:,3,:,:)),...
        squeeze(calm_qsiprep.followup.(parcellation).connectivity(:,3,:,:)));
    stacked_fc = cat(1,squeeze(calm_fmriprep.baseline.(parcellation).connectivity(:,:,:)),...
        squeeze(calm_fmriprep.followup.(parcellation).connectivity(:,:,:)));
    
    
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
    
    harmonized_fc = combat(vectorconnect_fc_nnz',ScannerID_FC,Covariates_FC,1);
    fprintf('harmonized functional connectomes for the %s parcellation.\n',parcellations(parcellation_idx));
    harmonized_sc = combat(vectorconnect_sc_nnz',ScannerID_SC,Covariates_SC,1);  
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
    % extract baseline vs control vs longitudinal outputs and assign to output
    % structure
    harmonized_calm_connectomes.control.fc.(parcellation) = squeeze(harmonized_fc_connectomes(Group_FC==0,:,:));
    harmonized_calm_connectomes.baseline.fc.(parcellation) = squeeze(harmonized_fc_connectomes(Group_FC==1,:,:));
    harmonized_calm_connectomes.followup.fc.(parcellation) = squeeze(harmonized_fc_connectomes(Group_FC==2,:,:));
    % repeat for sc!
    harmonized_sc_connectomes = zeros(size(stacked_sc));
    for n = 1:size(harmonized_sc_connectomes,1)
        vectordata = harmonized_sc_vector(n,:);
        matrix = tril(ones(nroi),-1);
        matrix(matrix > 0) = vectordata;
        harmonized_sc_connectomes(n,:,:) = matrix + matrix';
        clear vectordata matrix
    end
    harmonized_calm_connectomes.control.sc.(parcellation) = squeeze(harmonized_sc_connectomes(Group_SC==0,:,:));
    harmonized_calm_connectomes.baseline.sc.(parcellation) = squeeze(harmonized_sc_connectomes(Group_SC==1,:,:));
    harmonized_calm_connectomes.followup.sc.(parcellation) = squeeze(harmonized_sc_connectomes(Group_SC==2,:,:));
end

% save the harmonised connectomes - Alicja - haven't run this yet because I
% am going to save them later on, once checked
%save('connectomes/calm/harmonisation/harmonized_connectomes.mat','harmonized_calm_connectomes','-v7.3');

%% Part 6b - split out the behavioural / meta data so it matches the sections 
Meta.control.sc=Meta_SC(Group_SC==0,:);
Meta.baseline.sc=Meta_SC(Group_SC==1,:);
Meta.followup.sc=Meta_SC(Group_SC==2,:);
Meta.control.fc=Meta_FC(Group_FC==0,:);
Meta.baseline.fc=Meta_FC(Group_FC==1,:);
Meta.followup.fc=Meta_FC(Group_FC==2,:);

%make a bit of code for checking the lengths of connectomes and associated
%data
timepoints_list = {'baseline','followup','control'}; 

for timepoints_idx = 1:length(timepoints_list)
    if size(Meta.(timepoints_list{timepoints_idx}).sc,1)==size(harmonized_calm_connectomes.(timepoints_list{timepoints_idx}).sc.(parcellation),1)
     fprintf('the structural connectomes and the meta data are the correct length for the %s data.\n',timepoints_list{timepoints_idx});
    else 
     fprintf('ERROR the structural connectomes and the meta data are *NOT* the correct length for the %s data.\n',timepoints_list{timepoints_idx});
    end
end
%same for the functional data
for timepoints_idx = 1:length(timepoints_list)
    if size(Meta.(timepoints_list{timepoints_idx}).fc,1)==size(harmonized_calm_connectomes.(timepoints_list{timepoints_idx}).fc.(parcellation),1)
     fprintf('the functional connectomes and the meta data are the correct length for the %s data.\n',timepoints_list{timepoints_idx});
    else 
     fprintf('ERROR the functional connectomes and the meta data are *NOT* the correct length for the %s data.\n',timepoints_list{timepoints_idx});
    end
end

%%we need to get the unharmonized connectomes into exactly the same format
%%as their counterparts

for parcellation_idx = 1:length(parcellations)
    parcellation = parcellations(parcellation_idx);
   
    stacked_sc = cat(1,squeeze(calm_qsiprep.baseline.(parcellation).connectivity(:,3,:,:)),...
        squeeze(calm_qsiprep.followup.(parcellation).connectivity(:,3,:,:)));
    stacked_fc = cat(1,squeeze(calm_fmriprep.baseline.(parcellation).connectivity(:,:,:)),...
        squeeze(calm_fmriprep.followup.(parcellation).connectivity(:,:,:)));
    %structural connectomes
    UNharmonized_calm_connectomes.control.sc.(parcellation) = squeeze(stacked_sc(Group_SC==0,:,:));
    UNharmonized_calm_connectomes.baseline.sc.(parcellation) = squeeze(stacked_sc(Group_SC==1,:,:));
    UNharmonized_calm_connectomes.followup.sc.(parcellation) = squeeze(stacked_sc(Group_SC==2,:,:));
    %functional connectomes
    UNharmonized_calm_connectomes.control.fc.(parcellation) = squeeze(stacked_fc(Group_FC==0,:,:));
    UNharmonized_calm_connectomes.baseline.fc.(parcellation) = squeeze(stacked_fc(Group_FC==1,:,:));
    UNharmonized_calm_connectomes.followup.fc.(parcellation) = squeeze(stacked_fc(Group_FC==2,:,:));
    
end

%%now let's check the size of those against the behavioural / meta data -
%%should match

for timepoints_idx = 1:length(timepoints_list)
    if size(Meta.(timepoints_list{timepoints_idx}).sc,1)==size(UNharmonized_calm_connectomes.(timepoints_list{timepoints_idx}).sc.(parcellation),1)
     fprintf('the structural connectomes and the meta data are the correct length for the %s data.\n',timepoints_list{timepoints_idx});
    else 
     fprintf('ERROR the structural connectomes and the meta data are *NOT* the correct length for the %s data.\n',timepoints_list{timepoints_idx});
    end
end
%same for the functional data
for timepoints_idx = 1:length(timepoints_list)
    if size(Meta.(timepoints_list{timepoints_idx}).fc,1)==size(UNharmonized_calm_connectomes.(timepoints_list{timepoints_idx}).fc.(parcellation),1)
     fprintf('the functional connectomes and the meta data are the correct length for the %s data.\n',timepoints_list{timepoints_idx});
    else 
     fprintf('ERROR the functional connectomes and the meta data are *NOT* the correct length for the %s data.\n',timepoints_list{timepoints_idx});
    end
end


%% part 7 - checking densities of structural and functional connectomes 

%we need to adjust this list, because now we've separated out the control
%participants (i.e. the no-referred, so let's add a third category
timepoints_list = {'baseline','followup','control'}; 

%now I think this next section should run
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
        %sub_list = calm_fmriprep.(timepoints_list{timepoint_idx}).sample.sub;
        if timepoint_idx==1, sub_list=all_fc_ids(Group_FC==1);
        end
        if timepoint_idx==2, sub_list=all_fc_ids(Group_FC==2);
        end        
        if timepoint_idx==3, sub_list=all_fc_ids(Group_FC==0);
        end             
            
        for parcellation_idx = 1:length(parcellations)
            % create cell to hold the IDs of participants for whom thresholding
            % was unsuccessful i.e. not enough non-zero values
            failed_thresholding = {};
            % extract the unharmonised or harmonised unthresholded 
            % functional connectomes for this time point and parcellation.
            if harmonisation_idx == 1
                %unthresholded_connectomes = calm_fmriprep.(timepoints_list{timepoint_idx}).(parcellations(parcellation_idx)).connectivity;
                unthresholded_connectomes = UNharmonized_calm_connectomes.(timepoints_list{timepoint_idx}).fc.(parcellations(parcellation_idx));
            else
                unthresholded_connectomes = harmonized_calm_connectomes.(timepoints_list{timepoint_idx}).fc.(parcellations(parcellation_idx));
            end            
            % set the number of regions of interest (nroi)
            nroi = size(unthresholded_connectomes,2);
            % find the number of participants (nsub)
            nsub = size(unthresholded_connectomes,1);
            % output arrays for mean framewise displacement and thresholded
            % connectomes
           % mean_fwd = zeros(nsub,1); %shouldn't need this now as in matched meta data
            thresholded_connectomes = zeros(nsub,nroi,nroi);
            for sub_idx = 1:nsub
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
        if timepoint_idx==1; %add the baseline ids
        %thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub = cellstr(calm_fmriprep.(timepoints_list{timepoint_idx}).sample.sub);
        thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub = all_fc_ids(Group_FC==1);
        end
        if timepoint_idx==2; %add the follow-up ids
        %thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub = cellstr(calm_fmriprep.(timepoints_list{timepoint_idx}).sample.sub);
        thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub = all_fc_ids(Group_FC==2);
        end
        if timepoint_idx==3; %add the Control ids
        %thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub = cellstr(calm_fmriprep.(timepoints_list{timepoint_idx}).sample.sub);
        thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub = all_fc_ids(Group_FC==0);
        end        
        
        clear thresholded_connectomes
        end
       
      
    end
    
    %% I have commented this bit out because I think easier to calc at the end 
    % average fc across participants across time-points to generate a
    % group-representative connectome. for consistency, do this separately
    % for referred and non-referred participants.
    %for parcellation_idx = calm_fmriprep1:length(parcellations)
        % concatenate the individual (referred) connectomes across time
        % points
    %    concatenated_referrcalm_fmripreped = cat(1,thresholded.fc.(harmonisation_list(harmonisation_idx)).baseline.referred.(parcellations(parcellation_idx)).individual,...
    %        thresholded.fc.(harmonisation_list(harmonisation_idx)).followup.(parcellations(parcellation_idx)).individual);
        % average across participants and time points to find the
        % group-representatcalm_fmriprepive connectome
    %    thresholded.fc.(harmonisation_list(harmonisation_idx)).group.referred.(parcellations(parcellation_idx)) = squeeze(mean(concatenated_referred));
        % average across non-referred participants to find the
        % group-representatcalm_fmriprepive connectome (again, only at baseline, where
        % we have non-referred participants).
    %    thresholded.fc.(harmonisation_list(harmonisation_idx)).group.nonreferred.(parcellations(parcellation_idx)) = squeeze(mean(thresholded.fc.(harmonisation_list(harmonisation_idx)).baseline.nonreferred.(parcellations(parcellation_idx)).individual));
   % end
end 

%% part 8b - thresholding functional connectomes %%
% in line with margulies and colleagues (2016), for each row, retain the
% top 10% absolute strongest connections. we shall do this for both the
% harmonised and unharmonised connectomes. first, create a new structure to 
% hold the thresholded connectomes.
thresholded.sc = struct();
harmonisation_list = ["unharmonised","harmonised"];
for harmonisation_idx = 1:length(harmonisation_list)
    for timepoint_idx = 1:length(timepoints_list)
        % extract subject list
        %sub_list = calm_fmriprep.(timepoints_list{timepoint_idx}).sample.sub;
        if timepoint_idx==1, sub_list=all_sc_ids(Group_SC==1);
        end
        if timepoint_idx==2, sub_list=all_sc_ids(Group_SC==2);
        end        
        if timepoint_idx==3, sub_list=all_sc_ids(Group_SC==0);
        end             
            
        for parcellation_idx = 1:length(parcellations)
            % create cell to hold the IDs of participants for whom thresholding
            % was unsuccessful i.e. not enough non-zero values
            failed_thresholding = {};
            % extract the unharmonised or harmonised unthresholded 
            % functional connectomes for this time point and parcellation.
            if harmonisation_idx == 1
                %unthresholded_connectomes = calm_fmriprep.(timepoints_list{timepoint_idx}).(parcellations(parcellation_idx)).connectivity;
                unthresholded_connectomes = UNharmonized_calm_connectomes.(timepoints_list{timepoint_idx}).sc.(parcellations(parcellation_idx));
            else
                unthresholded_connectomes = harmonized_calm_connectomes.(timepoints_list{timepoint_idx}).sc.(parcellations(parcellation_idx));
            end            
            % set the number of regions of interest (nroi)
            nroi = size(unthresholded_connectomes,2);
            % find the number of participants (nsub)
            nsub = size(unthresholded_connectomes,1);
            % output arrays for mean framewise displacement and thresholded
            % connectomes
           % mean_fwd = zeros(nsub,1); %shouldn't need this now as in matched meta data
            thresholded_connectomes = zeros(nsub,nroi,nroi);
            for sub_idx = 1:nsub
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
        thresholded.sc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).(parcellations(parcellation_idx)).individual = thresholded_connectomes;
        fprintf('%d %s participants in the %s parcellation at %s failed functional thresholding.\n', ...
            length(unique(cellstr(failed_thresholding))), harmonisation_list{harmonisation_idx}, parcellations{parcellation_idx}, timepoints_list{timepoint_idx});
        % add the subject IDs. note that the sub IDs for the harmonised
        % and unharmonised data are the same, so it doesn't matter
        % where we pull from
        if timepoint_idx==1; %add the baseline ids
        %thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub = cellstr(calm_fmriprep.(timepoints_list{timepoint_idx}).sample.sub);
        thresholded.sc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub = all_sc_ids(Group_SC==1);
        end
        if timepoint_idx==2; %add the follow-up ids
        %thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub = cellstr(calm_fmriprep.(timepoints_list{timepoint_idx}).sample.sub);
        thresholded.sc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub = all_sc_ids(Group_SC==2);
        end
        if timepoint_idx==3; %add the Control ids
        %thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub = cellstr(calm_fmriprep.(timepoints_list{timepoint_idx}).sample.sub);
        thresholded.sc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoint_idx}).sub = all_sc_ids(Group_SC==0);
        end        
        
        clear thresholded_connectomes
        end
        
    end
    
    %% I have commented this bit out because I think easier to calc at the end 
    % average fc across participants across time-points to generate a
    % group-representative connectome. for consistency, do this separately
    % for referred and non-referred participants.
    %for parcellation_idx = calm_fmriprep1:length(parcellations)
        % concatenate the individual (referred) connectomes across time
        % points
    %    concatenated_referrcalm_fmripreped = cat(1,thresholded.fc.(harmonisation_list(harmonisation_idx)).baseline.referred.(parcellations(parcellation_idx)).individual,...
    %        thresholded.fc.(harmonisation_list(harmonisation_idx)).followup.(parcellations(parcellation_idx)).individual);
        % average across participants and time points to find the
        % group-representatcalm_fmriprepive connectome
    %    thresholded.fc.(harmonisation_list(harmonisation_idx)).group.referred.(parcellations(parcellation_idx)) = squeeze(mean(concatenated_referred));
        % average across non-referred participants to find the
        % group-representatcalm_fmriprepive connectome (again, only at baseline, where
        % we have non-referred participants).
    %    thresholded.fc.(harmonisation_list(harmonisation_idx)).group.nonreferred.(parcellations(parcellation_idx)) = squeeze(mean(thresholded.fc.(harmonisation_list(harmonisation_idx)).baseline.nonreferred.(parcellations(parcellation_idx)).individual));
   % end
end 

%%%Part 8c - now let's repeat our checks that everything still aligns

%make a bit of code for checking the lengths of connectomes and associated
%data
timepoints_list = {'baseline','followup','control'}; 

for timepoints_idx = 1:length(timepoints_list)
    if size(Meta.(timepoints_list{timepoints_idx}).sc,1)==size(thresholded.sc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoints_idx}).(parcellation).individual,1)
     fprintf('the structural connectomes and the meta data are the correct length for the %s data.\n',timepoints_list{timepoints_idx});
    else 
     fprintf('ERROR the structural connectomes and the meta data are *NOT* the correct length for the %s data.\n',timepoints_list{timepoints_idx});
    end
end
%same for the functional data
for timepoints_idx = 1:length(timepoints_list)
    if size(Meta.(timepoints_list{timepoints_idx}).fc,1)==size(thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoints_idx}).(parcellation).individual,1)
     fprintf('the functional connectomes and the meta data are the correct length for the %s data.\n',timepoints_list{timepoints_idx});
    else 
     fprintf('ERROR the functional connectomes and the meta data are *NOT* the correct length for the %s data.\n',timepoints_list{timepoints_idx});
    end
end

timepoints_idx=3 %specify which partition of the data you are checking.
%%But let's also check that the ordering is identical
size(Meta.(timepoints_list{timepoints_idx}).sc(1,1))
size(thresholded.sc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoints_idx}).sub(1))

clear same
for i = 1:size(thresholded.sc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoints_idx}).sub,1)
same(i)=mean(cell2mat(thresholded.sc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoints_idx}).sub(i)) == cell2mat(table2array(Meta.(timepoints_list{timepoints_idx}).sc(i,1))))==1
end

timepoints_idx=3 %specify which partition of the data you are checking.
clear same
for i = 1:size(thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoints_idx}).sub,1)
same(i)=mean(cell2mat(thresholded.fc.(harmonisation_list(harmonisation_idx)).(timepoints_list{timepoints_idx}).sub(i)) == cell2mat(table2array(Meta.(timepoints_list{timepoints_idx}).fc(i,1))))==1
end



%%%Part 9 - I just want to check that there are no scanner effects in the
%%%harmonised data

%structural data first
harmonised_list = {'harmonised','unharmonised'};
timepoints_list = {'baseline','followup','control'}; 
for harmonised_idx = 1:length(harmonised_list);
for timepoints_idx = 1:length(timepoints_list);
for i= 1:size(thresholded.sc.unharmonised.(timepoints_list{timepoints_idx}).schaefer200x7.individual,1)
sc.(harmonised_list{harmonised_idx}).efficiency.(timepoints_list{timepoints_idx})(i)=efficiency_wei(squeeze(thresholded.sc.(harmonised_list{harmonised_idx}).baseline.schaefer200x7.individual(i,:,:)));
sc.(harmonised_list{harmonised_idx}).density.(timepoints_list{timepoints_idx})(i)=density_und(squeeze(thresholded.sc.(harmonised_list{harmonised_idx}).baseline.schaefer200x7.individual(i,:,:)));
[~,sc.(harmonised_list{harmonised_idx}).modularity.(timepoints_list{timepoints_idx})(i)]=modularity_und(squeeze(thresholded.sc.(harmonised_list{harmonised_idx}).baseline.schaefer200x7.individual(i,:,:)));
sc.(harmonised_list{harmonised_idx}).transitivity.(timepoints_list{timepoints_idx})(i)=transitivity_wu(squeeze(thresholded.sc.(harmonised_list{harmonised_idx}).baseline.schaefer200x7.individual(i,:,:)));
sc.(harmonised_list{harmonised_idx}).assortativity.(timepoints_list{timepoints_idx})(i)=assortativity_wei(squeeze(thresholded.sc.(harmonised_list{harmonised_idx}).baseline.schaefer200x7.individual(i,:,:)),1);
end
end
end

%now do some simple GLMs to test if there are scanner effects in the
%harmonised data
graph_theory_metrics_list = ["efficiency","assortativity","density","transitivity","modularity"];
k=5; %select the graph metric you want
harmonised_idx=1; %select if you want to look at harmonised data or not 1=harmonised 2=unharmonised
y=transpose([sc.(harmonised_list{harmonised_idx}).(graph_theory_metrics_list{k}).baseline sc.(harmonised_list{harmonised_idx}).(graph_theory_metrics_list{k}).followup sc.(harmonised_list{harmonised_idx}).(graph_theory_metrics_list{k}).control])
X=vertcat(table2array(Meta.baseline.sc(:,[2:7])),table2array(Meta.followup.sc(:,[2:7])),table2array(Meta.control.sc(:,[2:7])))
[b,dev,stats]=glmfit(X,y)
stats.p%scanner effects are the 5th p value. But also encouraging to see others (1st = constant, 2nd = age, 3rd = sex, 4th = mean FWD, 6th time point, 7th referred)

%functional data next
harmonised_list = {'harmonised','unharmonised'};
timepoints_list = {'baseline','followup','control'}; 
for harmonised_idx = 1:length(harmonised_list);
for timepoints_idx = 1:length(timepoints_list);
for i= 1:size(thresholded.fc.unharmonised.(timepoints_list{timepoints_idx}).schaefer200x7.individual,1)
fc.(harmonised_list{harmonised_idx}).efficiency.(timepoints_list{timepoints_idx})(i)=efficiency_wei(squeeze(thresholded.fc.(harmonised_list{harmonised_idx}).baseline.schaefer200x7.individual(i,:,:)));
fc.(harmonised_list{harmonised_idx}).density.(timepoints_list{timepoints_idx})(i)=density_und(squeeze(thresholded.fc.(harmonised_list{harmonised_idx}).baseline.schaefer200x7.individual(i,:,:)));
[~,fc.(harmonised_list{harmonised_idx}).modularity.(timepoints_list{timepoints_idx})(i)]=modularity_und(squeeze(thresholded.fc.(harmonised_list{harmonised_idx}).baseline.schaefer200x7.individual(i,:,:)));
fc.(harmonised_list{harmonised_idx}).transitivity.(timepoints_list{timepoints_idx})(i)=transitivity_wu(squeeze(thresholded.fc.(harmonised_list{harmonised_idx}).baseline.schaefer200x7.individual(i,:,:)));
fc.(harmonised_list{harmonised_idx}).assortativity.(timepoints_list{timepoints_idx})(i)=assortativity_wei(squeeze(thresholded.fc.(harmonised_list{harmonised_idx}).baseline.schaefer200x7.individual(i,:,:)),1);
end
end
end

%now do some simple GLMs to test if there are scanner effects in the
%harmonised data
graph_theory_metrics_list = ["efficiency","assortativity","density","transitivity","modularity"];
k=5; %select the graph metric you want
harmonised_idx=1; %select if you want to look at harmonised data or not 1=harmonised 2=unharmonised
y=transpose([fc.(harmonised_list{harmonised_idx}).(graph_theory_metrics_list{k}).baseline fc.(harmonised_list{harmonised_idx}).(graph_theory_metrics_list{k}).followup fc.(harmonised_list{harmonised_idx}).(graph_theory_metrics_list{k}).control])
X=vertcat(table2array(Meta.baseline.fc(:,[2:7])),table2array(Meta.followup.fc(:,[2:7])),table2array(Meta.control.fc(:,[2:7])))
[b,dev,stats]=glmfit(X,y)
stats.p%scanner effects are the 5th p value. But also encouraging to see others (1st = constant, 2nd = age, 3rd = sex, 4th = mean FWD, 6th time point, 7th referred)


%%%Part 10 - SAVE THE STRUCTURES :)
% no missing communicability values for harmonised connectomes!
save('/imaging/projects/external/nkir/analyses/duncan_thresholded_structural_and_functional_connectomes.mat','thresholded','-v7.3');
save('/imaging/projects/external/nkir/analyses/duncan_meta.mat','Meta','-v7.3');


