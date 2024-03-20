% (C) Copyright Jacek Matuszewski 2024
%% Set directories
clear;
clc;

this_dir = '/Volumes/Slim_Reaper/Projects/Language_MVPA/code/src/bidspm-winner_takes_all';

%addpath(fullfile(this_dir, '..', '..', 'lib', 'bidspm'));

%bidspm();

opt.dir.root = fullfile(this_dir, '..', '..', '..');
opt.dir.derivatives = fullfile(opt.dir.root, 'outputs', 'derivatives');
opt.dir.stats = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-stats');
opt.dir.GroupStats = fullfile(opt.dir.stats, 'derivatives', 'bidspm-groupStats');
opt.dir.wta = fullfile(opt.dir.stats, 'derivatives', 'bidspm-groupStats','WinnerTakesAll');
opt.dir.rois = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-roi');

%% %% Pick your subjects

%List the sub numbers
subNum = {'01', '02','03','04','05',...
    '06','07','08','09','10','11','12',...
    '13','14','15','16','17','18','19','20'};

% subNum = {'01'};
%List the groups
opt.group = {'blind', 'sighted'};

%paste numbers with group to create full subjects list
opt.subjects = {};
for i = 1:length(opt.group)
    opt.subjects = [opt.subjects, strcat(opt.group{i}, subNum)];
end

%TESTING ON 1 sub

%opt.subjects = {'blind01'};

%% Set Design & Task & Regions of interest

opt.modality = {'Read', 'Speech'};

opt.taskName = 'MultimodalReadSpeech';

opt.dir.statsTask = 'task-MultimodalReadSpeech_space-IXI549Space_FWHM-6_node-univSummary';

opt.roi.name = {'hemi-bilateral_space-MNI_atlas-JuBrain_label-V1_mask.nii', ...
    'hemi-L_space-MNI_atlas-JuBrain_label-FG2_mask.nii', ...
    'hemi-L_space-MNI_atlas-JuBrain_label-FG4_mask.nii', ...
    'hemi-L_space-MNI_atlas-JuBrain_label-MTG_mask.nii', ...
    'hemi-L_space-MNI_atlas-JuBrain_label-Broca_mask.nii'};

opt.roi.label = {'V1', 'FG2', 'FG4', 'MTG', 'Broca'};



%% READING WTA

%READING CONTRAST NUMBERS
% 49 : words
% 50 : pseudo
% 51 : control


for iSub = 1:numel(opt.subjects)
    
    for iROI = 1:numel(opt.roi.name)
        
        %Define location for ROI mask in current subject
        mask = fullfile(opt.dir.rois, ['sub-', opt.subjects{iSub}], 'roi',...
            ['sub-', opt.subjects{iSub},'_',opt.roi.name{iROI}]);
        
        %Output file name
        output_name = ['sub-', opt.subjects{iSub},'_atlas-JuBrain_desc-WTA6runsRead_', 'label-',opt.roi.label{iROI},'_dseg.nii'];
        
        %Files to read (contrast images vs betas?)
        stat_dir = fullfile(opt.dir.stats, ['sub-', opt.subjects{iSub}], opt.dir.statsTask);
        
        files = [fullfile(stat_dir, 'con_0049.nii');
            fullfile(stat_dir, 'con_0050.nii');
            fullfile(stat_dir, 'con_0051.nii')];
        
        % OUTPUT IS THE stat_dir!
        %wta_read_hdr = winner_take_all(files, mask, output_name);
        
        winner_take_all(files, mask, output_name);
        % THIS RESULTS IN a NAN MAP!
        %view_results(files, wta_read_hdr);
        
    end
    
end

%% SPEECH WTA

%SPEECH CONTRAST NUMBERS
% 52 : words
% 53 : pseudo
% 54 : control

for iSub = 1:numel(opt.subjects)
    
    for iROI = 1:numel(opt.roi.name)
        
        mask = fullfile(opt.dir.rois, ['sub-', opt.subjects{iSub}], 'roi',...
            ['sub-', opt.subjects{iSub},'_',opt.roi.name{iROI}]);
        
        %Output file name
        output_name = ['sub-', opt.subjects{iSub},'_atlas-JuBrain_desc-WTA6runsSpeech_', 'label-',opt.roi.label{iROI},'_dseg.nii'];
        
        %Files to read (contrast images vs betas?)
        stat_dir = fullfile(opt.dir.stats, ['sub-', opt.subjects{iSub}], opt.dir.statsTask);
        
        files = [fullfile(stat_dir, 'con_0052.nii');
            fullfile(stat_dir, 'con_0053.nii');
            fullfile(stat_dir, 'con_0054.nii')];
        
        % OUTPUT IS THE stat_dir!
        winner_take_all(files, mask, output_name);
        
        %if you would like to see the results with check_reg (for testing)
        %wta_read_hdr = winner_take_all(files, mask, output_name);
        %view_results(files, wta_read_hdr);
    end
end


