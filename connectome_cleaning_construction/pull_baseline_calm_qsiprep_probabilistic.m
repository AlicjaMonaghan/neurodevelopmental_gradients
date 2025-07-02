%% pull qsi data from qsiprep and qsirecon
% written by danyal akarca and adapted by alicja monaghan. in this script,
% we pull qsiprep outputs for ~20 baseline CALM participants whose
% diffusion-weighted scans were previously unprocessed.
clear; clc;
% set qsiprep directory
qsiprepDir = '/imaging/projects/external/nkir/analyses/calm-I/qsiprep/outdir/final_output/';
% change to current directory
cd(qsiprepDir);
% add path to read required files
addpath('/imaging/astle/am10/Toolboxes/');
%% set sample hyperparameters
sub_id_structure = dir(qsiprepDir);
% remove the first two rows of the structure which show decimal points.
sub_id_structure(1:2) = [];
nsub = length(sub_id_structure);
%% set labels
% parcellation labels
parcellations = {'aal116','aicha384','brainnetome246',...
    'gordon333','power264','schaefer100x17','schaefer100x7',...
    'schaefer200x17','schaefer200x7','schaefer400x17','schaefer400x7'};
nparcellations = length(parcellations);
% types of data
datatypes = {'radius2_count_connectivity','radius2_meanlength_connectivity', ...
    'sift_invnodevol_radius2_count_connectivity','sift_radius2_count_connectivity'};
ndatatypes = length(datatypes);
% quality control labels
qclabels = string(['mean_framewise_displacement','max_framewise_displacement',...
    'max_rotation','max_translation','max_rel_rotation',...
    'max_rel_translation','t1_dice_distance']);
%% load and organise external data to be saved
% parcellation metadata
aal116 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/aal116_info.mat');
aal116 = aal116.aal116;
aicha384 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/aicha384_info.mat');
aicha384 = aicha384.aicha384_info;
brainnetome246 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/brainnetome246_info.mat');
brainnetome246 = brainnetome246.brainnetome246;
gordon333 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/gordon333_info.mat');
gordon333 = gordon333.gordon333;
power264 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/power264_info.mat');
power264 = power264.power264;
schaefer100x17 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer100x17_1mm_info.mat');
schaefer100x17 = schaefer100x17.schaefer100x17_1mm_info;
schaefer100x7 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer100x7_1mm_info.mat');
schaefer100x7 = schaefer100x7.schaefer100x7_1mm_info;
schaefer200x17 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer200x17_1mm_info.mat');
schaefer200x17 = schaefer200x17.schaefer200x17_1mm_info;
schaefer200x7 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer200x7_1mm_info.mat');
schaefer200x7 = schaefer200x7.schaefer200x7_1mm_info;
schaefer400x17 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer400x17_1mm_info.mat');
schaefer400x17 = schaefer400x17.schaefer400x17_1mm_info;
schaefer400x7 = load('/imaging/astle/users/da04/PhD/qsiprep_data/data/schaefer400x7_1mm_info.mat');
schaefer400x7 = schaefer400x7.schaefer400x7_1mm_info;
% organise data for later allocation
% aal116
aal116_coordinates = [aal116.x_mni aal116.y_mni aal116.z_mni];
aal116_hemi = string(aal116.hemi);
aal116_lobe = string(aal116.lobe);
% aicha384
aicha384_coordinates = [aicha384.x_mni aicha384.y_mni aicha384.z_mni];
aicha384_volume = aicha384.volume;
% brainnetome246
brainnetome246_coordinates = [brainnetome246.x_mni brainnetome246.y_mni brainnetome246.z_mni];
brainnetome246_hemi = string(brainnetome246.hemi);
brainnetome246_lobe = string(brainnetome246.lobe);
% gordon333
gordon333_coordinates = [gordon333.x_mni gordon333.y_mni gordon333.z_mni];
gordon333_hemi = string(gordon333.hemi);
gordon333_lobe = string(gordon333.lobe);
% power264
power264_coordinates = [power264.x_mni power264.y_mni power264.z_mni];
power264_hemi = string(power264.hemi);
power264_lobe = string(power264.lobe);
% schaefer100x17
schaefer100x17_coordinates = [schaefer100x17.x_mni schaefer100x17.y_mni schaefer100x17.z_mni];
% schaefer100x7
schaefer100x7_coordinates = [schaefer100x7.x_mni schaefer100x7.y_mni schaefer100x7.z_mni];
% schaefer200x17
schaefer200x17_coordinates = [schaefer200x17.x_mni schaefer200x17.y_mni schaefer200x17.z_mni];
% schaefer200x7
schaefer200x7_coordinates = [schaefer200x7.x_mni schaefer200x7.y_mni schaefer200x7.z_mni];
% schaefer400x17
schaefer400x17_coordinates = [schaefer400x17.x_mni schaefer400x17.y_mni schaefer400x17.z_mni];
% schaefer400x7
schaefer400x7_coordinates = [schaefer400x7.x_mni schaefer400x7.y_mni schaefer400x7.z_mni];
% partipant data
participants = readtable('/imaging/projects/cbu/calm/CALM-II_BIDS/Updated_CALM_MRI_IDs.xlsx', ...
    'ReadVariableNames',true,'VariableNamingRule','preserve','Format','auto');
participants.participant_id = string(participants.ID);
participants.CBU_ID = string(participants.("MRI-ID-T2"));
%% initialise ready to collect sample data
% initialise the sample data
% qc
framedisplace = cell(nsub,1);
qcsummary = zeros(nsub,7);
% aal116
aal116 = zeros(nsub,ndatatypes,116,116);
% aicha384
aicha384 = zeros(nsub,ndatatypes,384,384);
% brainectome246
brainnetome246 = zeros(nsub,ndatatypes,246,246);
% gordon333
gordon333 = zeros(nsub,ndatatypes,333,333);
% power264
power264 = zeros(nsub,ndatatypes,264,264);
% schaefer100x17
schaefer100x17 = zeros(nsub,ndatatypes,100,100);
% schaefer100x7
schaefer100x7  = zeros(nsub,ndatatypes,100,100);
% schaefer200x17
schaefer200x17 = zeros(nsub,ndatatypes,200,200);
% schaefer200x7
schaefer200x7 = zeros(nsub,ndatatypes,200,200);
% schaefer400x17
schaefer400x17 = zeros(nsub,ndatatypes,400,400);
% schaefer400x7
schaefer400x7 = zeros(nsub,ndatatypes,400,400);
%% index subject files from the qsiprep and qsirecon analysis
for i = 1:nsub
    tic;
    % set the subject from the list
    sub = sub_id_structure(i).name;
    % display
    fprintf('loading %s data...\n',sub);
    % set confounds directory
    confoundsDir = sprintf('%s%s/%s_confounds.tsv', qsiprepDir,sub,sub);
    % set qc directory
    qcDir = sprintf('%s%s/%s_desc-ImageQC_dwi.csv',qsiprepDir,sub,sub);
    % set network directory
    connectivityDir = sprintf('%s%s/combined_connectivity.mat',qsiprepDir,sub);
    % load movement data
    % load confounds (note, load to no variable to see the names)
    confounds = readtable(confoundsDir,'FileType','delimitedtext');
    % load quality control summary
    qc = readtable(qcDir);
    qc = [qc.mean_fd qc.max_fd qc.max_rotation qc.max_translation qc.max_rel_rotation qc.max_rel_translation qc.t1_dice_distance];
    % compute framewise displacement for a single subject
    fd = confounds.framewise_displacement';
    fd(isnan(fd)) = [];
    % load connectivity data
    load(connectivityDir);
    % get the correct index
    subi = i;
    % assign connectivity measures    
    % framedisplace
    framedisplace{subi,1} = fd;
    % qcsummary
    qcsummary(subi,:) = qc;
    % aal116
    aal116(subi,1,:,:)        = aal116_radius2_count_connectivity;
    aal116(subi,2,:,:)        = aal116_radius2_meanlength_connectivity;
    aal116(subi,3,:,:)        = aal116_sift_invnodevol_radius2_count_connectivity;
    aal116(subi,4,:,:)        = aal116_sift_radius2_count_connectivity;
    % aicha384
    aicha384(subi,1,:,:)        = aicha384_radius2_count_connectivity;
    aicha384(subi,2,:,:)        = aicha384_radius2_meanlength_connectivity;
    aicha384(subi,3,:,:)        = aicha384_sift_invnodevol_radius2_count_connectivity;
    aicha384(subi,4,:,:)        = aicha384_sift_radius2_count_connectivity;
    % brainnectome246
    brainnetome246(subi,1,:,:)        = brainnetome246_radius2_count_connectivity;
    brainnetome246(subi,2,:,:)        = brainnetome246_radius2_meanlength_connectivity;
    brainnetome246(subi,3,:,:)        = brainnetome246_sift_invnodevol_radius2_count_connectivity;
    brainnetome246(subi,4,:,:)        = brainnetome246_sift_radius2_count_connectivity;
    % gordon333
    gordon333(subi,1,:,:)        = gordon333_radius2_count_connectivity;
    gordon333(subi,2,:,:)        = gordon333_radius2_meanlength_connectivity;
    gordon333(subi,3,:,:)        = gordon333_sift_invnodevol_radius2_count_connectivity;
    gordon333(subi,4,:,:)        = gordon333_sift_radius2_count_connectivity;
    % power264
    power264(subi,1,:,:)        = power264_radius2_count_connectivity;
    power264(subi,2,:,:)        = power264_radius2_meanlength_connectivity;
    power264(subi,3,:,:)        = power264_sift_invnodevol_radius2_count_connectivity;
    power264(subi,4,:,:)        = power264_sift_radius2_count_connectivity;
    % schaefer100x17
    schaefer100x17(subi,1,:,:)        = schaefer100x17_radius2_count_connectivity;
    schaefer100x17(subi,2,:,:)        = schaefer100x17_radius2_meanlength_connectivity;
    schaefer100x17(subi,3,:,:)        = schaefer100x17_sift_invnodevol_radius2_count_connectivity;
    schaefer100x17(subi,4,:,:)        = schaefer100x17_sift_radius2_count_connectivity;
    % schaefer100x7
    schaefer100x7(subi,1,:,:)        = schaefer100x7_radius2_count_connectivity;
    schaefer100x7(subi,2,:,:)        = schaefer100x7_radius2_meanlength_connectivity;
    schaefer100x7(subi,3,:,:)        = schaefer100x7_sift_invnodevol_radius2_count_connectivity;
    schaefer100x7(subi,4,:,:)        = schaefer100x7_sift_radius2_count_connectivity;
    % schaefer200x17
    schaefer200x17(subi,1,:,:)        = schaefer200x17_radius2_count_connectivity;
    schaefer200x17(subi,2,:,:)        = schaefer200x17_radius2_meanlength_connectivity;
    schaefer200x17(subi,3,:,:)        = schaefer200x17_sift_invnodevol_radius2_count_connectivity;
    schaefer200x17(subi,4,:,:)        = schaefer200x17_sift_radius2_count_connectivity;
    % schaefer200x7
    schaefer200x7(subi,1,:,:)        = schaefer200x7_radius2_count_connectivity;
    schaefer200x7(subi,2,:,:)        = schaefer200x7_radius2_meanlength_connectivity;
    schaefer200x7(subi,3,:,:)        = schaefer200x7_sift_invnodevol_radius2_count_connectivity;
    schaefer200x7(subi,4,:,:)        = schaefer200x7_sift_radius2_count_connectivity;
    % schaefer400x17
    schaefer400x17(subi,1,:,:)        = schaefer400x17_radius2_count_connectivity;
    schaefer400x17(subi,2,:,:)        = schaefer400x17_radius2_meanlength_connectivity;
    schaefer400x17(subi,3,:,:)        = schaefer400x17_sift_invnodevol_radius2_count_connectivity;
    schaefer400x17(subi,4,:,:)        = schaefer400x17_sift_radius2_count_connectivity;
    % schaefer400x7
    schaefer400x7(subi,1,:,:)        = schaefer400x7_radius2_count_connectivity;
    schaefer400x7(subi,2,:,:)        = schaefer400x7_radius2_meanlength_connectivity;
    schaefer400x7(subi,3,:,:)        = schaefer400x7_sift_invnodevol_radius2_count_connectivity;
    schaefer400x7(subi,4,:,:)        = schaefer400x7_sift_radius2_count_connectivity;
    % display end
    t = toc;
    fprintf('loaded %s data (%.3g seconds)\n',sub,t);
end
%% assign to a struct
% form a empty struct
additional_calm_qsiprep = struct();
% sample
additional_calm_qsiprep.sample.id.sub = {sub_id_structure.name}';
% info
additional_calm_qsiprep.info.datatypes = datatypes';
additional_calm_qsiprep.info.meta = 'all parcellated data is in the order: subject, data type, data';
additional_calm_qsiprep.info.parcellations = parcellations';
% quality control
additional_calm_qsiprep.qc.framedisplacement = framedisplace;
additional_calm_qsiprep.qc.summary = qcsummary;
additional_calm_qsiprep.qc.summarylabels = qclabels';
% aal 116
additional_calm_qsiprep.aal116.connectivity = aal116;
% aicha384
additional_calm_qsiprep.aicha384.connectivity = aicha384;
% brainnetome246
additional_calm_qsiprep.brainnetome246.connectivity = brainnetome246;
% gordon333
additional_calm_qsiprep.gordon333.connectivity = gordon333;
% power264
additional_calm_qsiprep.power264.connectivity = power264;
% schaefer100x17
additional_calm_qsiprep.schaefer100x17.connectivity = schaefer100x17;
% schaefer100x7
additional_calm_qsiprep.schaefer100x7.connectivity = schaefer100x7;
% schaefer200x17
additional_calm_qsiprep.schaefer200x17.connectivity = schaefer200x17;
% schaefer200x7
additional_calm_qsiprep.schaefer200x7.connectivity = schaefer200x7;
% schaefer400x17
additional_calm_qsiprep.schaefer400x17.connectivity = schaefer400x17;
% schaefer400x7
additional_calm_qsiprep.schaefer400x7.connectivity = schaefer400x7;
%% reconcile with existing calm qsiprep probabilistic connectomes!
% load those processed by danyal akarca
calm_qsiprep = load('/imaging/astle/users/da04/Postdoc/qsiprep_optimise/data/calm_qsiprep_csd_ACT.mat');
calm_qsiprep = calm_qsiprep.calm_qsiprep_csd_ACT;
% remove 'prep' field from calm_qsiprep for easier concatenation
calm_qsiprep = rmfield(calm_qsiprep,'prep');
combined = [calm_qsiprep,additional_calm_qsiprep];
% format nicely
% dti
% aal116
complete_baseline_calm = struct();
complete_baseline_calm.aal116.connectivity = cat(1,combined(1).aal116.connectivity,combined(2).aal116.connectivity);
complete_baseline_calm.aal116.coordinates = aal116_coordinates;
complete_baseline_calm.aal116.hemi = aal116_hemi;
complete_baseline_calm.aal116.lobe = aal116_lobe;
complete_baseline_calm.aal116.regionlabels = string(aal116_region_labels);
% aicha384
complete_baseline_calm.aicha384.connectivity = cat(1,combined(1).aicha384.connectivity,combined(2).aicha384.connectivity);
complete_baseline_calm.aicha384.coordinates = aicha384_coordinates;
complete_baseline_calm.aicha384.regionlabels = string(aicha384_region_labels);
complete_baseline_calm.aicha384.volume = aicha384_volume;
% brainnetome246
complete_baseline_calm.brainnetome246.connectivity = cat(1,combined(1).brainnetome246.connectivity,combined(2).brainnetome246.connectivity);
complete_baseline_calm.brainnetome246.hemi = brainnetome246_hemi;
complete_baseline_calm.brainnetome246.lobe = brainnetome246_lobe;
complete_baseline_calm.brainnetome246.regionlabels = string(brainnetome246_region_labels);
% gordon333
complete_baseline_calm.gordon333.connectivity = cat(1,combined(1).gordon333.connectivity,combined(2).gordon333.connectivity);
complete_baseline_calm.gordon333.coordinates = gordon333_coordinates;
complete_baseline_calm.gordon333.hemi = gordon333_hemi;
complete_baseline_calm.gordon333.lobe = gordon333_lobe;
complete_baseline_calm.gordon333.regionlabels = string(gordon333_region_labels);
% power264
complete_baseline_calm.power264.connectivity = cat(1,combined(1).power264.connectivity,combined(2).power264.connectivity);
complete_baseline_calm.power264.coordinates = power264_coordinates;
complete_baseline_calm.power264.hemi = power264_hemi;
complete_baseline_calm.power264.lobe = power264_lobe;
complete_baseline_calm.power264.regionlabels = string(power264_region_labels);
% schaefer100x7
complete_baseline_calm.schaefer100x7.connectivity = cat(1,combined(1).schaefer100x7.connectivity,combined(2).schaefer100x7.connectivity);
complete_baseline_calm.schaefer100x7.coordinates = schaefer100x7_coordinates;
complete_baseline_calm.schaefer100x7.regionlabels = string(schaefer100x7_region_labels);
% schaefer100x17
complete_baseline_calm.schaefer100x17.connectivity = cat(1,combined(1).schaefer100x17.connectivity,combined(2).schaefer100x17.connectivity);
complete_baseline_calm.schaefer100x17.coordinates = schaefer100x17_coordinates;
complete_baseline_calm.schaefer100x17.regionlabels = string(schaefer100x17_region_labels);
% schaefer200x7
complete_baseline_calm.schaefer200x7.connectivity = cat(1,combined(1).schaefer200x7.connectivity,combined(2).schaefer200x7.connectivity);
complete_baseline_calm.schaefer200x7.coordinates = schaefer200x7_coordinates;
complete_baseline_calm.schaefer200x7.regionlabels = string(schaefer200x7_region_labels);
% schaefer200x17
complete_baseline_calm.schaefer200x17.connectivity = cat(1,combined(1).schaefer200x17.connectivity,combined(2).schaefer200x17.connectivity);
complete_baseline_calm.schaefer200x17.coordinates = schaefer200x17_coordinates;
complete_baseline_calm.schaefer200x17.regionlabels = string(schaefer200x17_region_labels);
% schaefer400x7
complete_baseline_calm.schaefer400x7.connectivity = cat(1,combined(1).schaefer400x7.connectivity,combined(2).schaefer400x7.connectivity);
complete_baseline_calm.schaefer400x7.coordinates = schaefer400x7_coordinates;
complete_baseline_calm.schaefer400x7.regionlabels = string(schaefer400x7_region_labels);
% schaefer400x17
complete_baseline_calm.schaefer400x17.connectivity = cat(1,combined(1).schaefer400x17.connectivity,combined(2).schaefer400x17.connectivity);
complete_baseline_calm.schaefer400x17.coordinates = schaefer400x17_coordinates;
complete_baseline_calm.schaefer400x17.regionlabels = string(schaefer400x17_region_labels);
% sample
complete_baseline_calm.sample.id.sub = vertcat(cellstr(combined(1).sample.id.sub)',combined(2).sample.id.sub);
% info
complete_baseline_calm.info.datatypes = datatypes';
complete_baseline_calm.info.meta = 'all parcellated data is in the order: subject, data type, data';
complete_baseline_calm.info.parcellations = parcellations';
% quality control
complete_baseline_calm.qc.framedisplacement = vertcat(combined(1).qc.framedisplacement,cell2mat(combined(2).qc.framedisplacement));
complete_baseline_calm.qc.summary = vertcat(combined(1).qc.summary,combined(2).qc.summary);
complete_baseline_calm.qc.summarylabels = qclabels';
% prep
complete_baseline_calm.prep = "qsiprep and qsirecon run by danyal akarca and alicja monaghan, 05/23";
% and save!
save('/imaging/projects/cbu/calm/baseline_calm_csd_ACT_qsiprep','complete_baseline_calm','-v7.3');

