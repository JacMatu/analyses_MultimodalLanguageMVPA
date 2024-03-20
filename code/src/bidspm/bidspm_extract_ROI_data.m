% (C) Copyright 2019 bidspm developers

%%% Script for preprocessing pipeline with bidspm()
%%% Due to SliceTiming errors in DCM files slice timing is skipped
%%% to run, needs bidspm() installed and path added and saved. 

clear;
clc;

%% Paths after moving to external hard drive! 

%% initialize bidspm() for this matlab session
%bidspm in code/lib/bidspm
addpath(fullfile(pwd, '..', '..','lib', 'bidspm'));

bidspm();

%% set up BIDS folders and paths
this_dir = fileparts(mfilename('fullpath'));
root_dir = fullfile(this_dir, '..', '..', '..');

opt.dir.bids_dir = fullfile(root_dir,'inputs','RAWDATA');
opt.dir.deriv = fullfile(root_dir, 'outputs','derivatives');
opt.dir.gr_stats = fullfile(opt.dir.deriv, 'bidspm-stats','derivatives','bidspm-groupStats');
opt.dir.roi = fullfile(opt.dir.deriv, 'bidspm-roi');
%% define subjects and tasks to preprocess
%subject_label = {''}; % ALL SUBJECTS FROM RAW
%subject_label = {'blind'}; % ONLY BLIND GROUP
opt.sub.group = {'blind','sighted'}; % ONLY SIGHTED GROUP

opt.space = 'space-IXI549Space';

opt.task_label = 'MultimodalReadSpeech'; %1-back task with spoken and written words, pseudowords and sensory control 
opt.fwhm.preproc = '6';
opt.fwhm.contrast = '0';

%List ROI files
%opt.roi.list = cellstr(spm_select('List',fullfile(opt.dir.roi,'group'), '.*JuBrain.*.nii'));
opt.roi.list = cellstr(spm_select('List',fullfile(opt.dir.roi,'group'), '.*VWFA.*.nii'));

%List contrasts that define GLMs to extract the data from
opt.stats.contrasts = {'ReadControl', 'ReadWord', 'ReadPseudoword', 'SpeechControl', 'SpeechWord', 'SpeechPseudoword'};
%Add 'block' if necessary
%opt.stats.contrasts = strcat('block', opt.stats.contrasts);

% This is how model DIRS are called: 
% sub-blind_task-MultimodalReadSpeech_space-IXI549Space_FWHM-6_conFWHM-0_node-withinGroup_contrast-blockReadControl

% This is how model should be addresed with opt parameters
% ['sub-', opt.sub.group, '_task-', opt.task_label, '_', opt.space,...
% '_FWHM-',opt.fwhm.preproc,'_conFWHM-',opt.fwhm.contrast,...
%     '_node-withinGroup_contrast-',opt.stats.contrasts];

%% 
% Testing
%iGr = 1; iRoi=1;iCon = 1;

%This is your output csv file for all the ROIs, groups and conditions
%f=fopen(fullfile(opt.dir.gr_stats,'JuBrain_stats.csv'),'w');
f=fopen(fullfile(opt.dir.gr_stats,'cVWFA_stats.csv'),'w');
 
 
for iGr = 1:numel(opt.sub.group)
    for iRoi = 1:numel(opt.roi.list)
        for iCon = 1:numel(opt.stats.contrasts)
            
            % Get a proper GLM name and dir
            glm_name =  ['sub-', opt.sub.group{iGr}, '_task-', ...
                opt.task_label, '_', opt.space,'_FWHM-',...
                opt.fwhm.preproc,'_conFWHM-',opt.fwhm.contrast,...
                '_node-withinGroup_contrast-','block',opt.stats.contrasts{iCon}];
            
            glm_file = fullfile(opt.dir.gr_stats, glm_name,'SPM.mat');
            
            % Get ROI mask dir
            roi_msk = fullfile(opt.dir.roi,'group', opt.roi.list{iRoi});
            
            % Load GLM and get the input file list
            load(glm_file);
            
            for i=1:size(SPM.xY.VY,1)
                files{i}=SPM.xY.VY(i).fname; 
            end
            
            % Extract the data from each ROI and file
            %roi_vals = [];
            
            %output file
           % f=fopen(fullfile(opt.dir.gr_stats,'stats.csv'),'w');

            for i=1:size(SPM.xY.VY,1)
                roi_mean=mean(spm_summarise(files{i},roi_msk),2);
                %fprintf(f,'%s,%g\n',files{i},roi_mean);
                  
                fprintf(f,'%s,%s,%s,%s,%g\n',...,
                    files{i},...
                    opt.sub.group{iGr},...
                    opt.roi.list{iRoi},...
                    opt.stats.contrasts{iCon},...
                    roi_mean);
                    
            end
            
           % fclose(f);

            % Rename the output file? 
            %fname = ['group-',opt.sub.group{iGr},'_',...
            %    opt.roi.list{iRoi}(1:end-4),'_',...
            %    opt.stats.contrasts{iCon},'_stats.csv'];
            
            %movefile(fullfile(opt.dir.gr_stats,'stats.csv'), ...
            %    fullfile(opt.dir.gr_stats, fname));
            
        end
        
    end
end

 fclose(f);
