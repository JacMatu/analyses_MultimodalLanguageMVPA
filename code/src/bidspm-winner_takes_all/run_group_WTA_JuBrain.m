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

%% Set Design & Task & Regions of interest
opt.group = {'blind', 'sighted'};

opt.modality = {'Read', 'Speech'};

opt.taskName = 'MultimodalReadSpeech';

opt.dir.statsTask = 'task-MultimodalReadSpeech_space-IXI549Space_FWHM-6_node-univSummary';


%WATCH OUT THESE ARE "NEW" RESLICED ROIS IN WTA DIR 
opt.roi.name = {'hemi-bilateral_space-MNI_atlas-JuBrain_label-V1_desc-ReslicedToStats_mask.nii', ...
    'hemi-L_space-MNI_atlas-JuBrain_label-FG2_desc-ReslicedToStats_mask.nii', ...
    'hemi-L_space-MNI_atlas-JuBrain_label-FG4_desc-ReslicedToStats_mask.nii', ...
    'hemi-L_space-MNI_atlas-JuBrain_label-MTG_desc-ReslicedToStats_mask.nii', ...
    'hemi-L_space-MNI_atlas-JuBrain_label-Broca_desc-ReslicedToStats_mask.nii'};

opt.roi.label = {'V1', 'FG2', 'FG4', 'MTG', 'Broca'};

%do you need to reslice your ROIs to group stats? (typically only once)
opt.roi.reslice = false; 
%% ADD GROUP MAPS FOR VISUALIZATION PURPOSSES
% https://elifesciences.org/articles/50732#s4 FIG 1B description:
% average BETAS across groups (and modalities) and run classification on
% singular group maps


if ~exist(opt.dir.wta, 'dir')
    mkdir(opt.dir.wta)
end

if ~exist(fullfile(opt.dir.wta, 'group_vis'), 'dir')
    mkdir(fullfile(opt.dir.wta, 'group_vis'))
end

if ~exist(fullfile(opt.dir.wta, 'stats'), 'dir')
    mkdir(fullfile(opt.dir.wta, 'stats'))
end


%desc-ReslicedToStats_

%% RESLICE ROIS to STATS if needed 

if opt.roi.reslice
    
    opt.roi.list = cellstr(spm_select('FPList', fullfile(opt.dir.rois, 'group'), '.*atlas-JuBrain.*_mask.nii'));
    opt.roi.names = cellstr(spm_select('List', fullfile(opt.dir.rois, 'group'), '.*atlas-JuBrain.*_mask.nii'));
    
    
    betaReference = fullfile(opt.dir.GroupStats,'sub-blind_task-MultimodalReadSpeech_space-IXI549Space_FWHM-6_conFWHM-0_node-withinGroup_contrast-blockReadWord', 'con_0001.nii');
    outputPath = fullfile(opt.dir.wta,'roi');
    
    if ~exist(outputPath, 'dir')
        mkdir(outputPath)
    end
    
    for iROI = 1:length(opt.roi.list)

        this_roi = opt.roi.list{iROI};
        
        reslicedImages = resliceRoiImages(betaReference, this_roi);
        bidslikeName = strrep(opt.roi.names{iROI}, '_mask.nii', '_desc-ReslicedToStats_mask.nii');
        movefile(reslicedImages, fullfile(outputPath,bidslikeName),'f');
        
        
        
    end
    
    
end


%% RUN GROUP WTA'S FOR VISUALIZATIONS !
%loop that across all groupxmodalityxcondition to get group maps

for iGr = 1:numel(opt.group)
    
    for iMod = 1:numel(opt.modality)
        
        for iROI = 1:numel(opt.roi.name)
            
            %Find the base OneSampleT name to be finished with conditions
            %below
            stat_basename = ['sub-',opt.group{iGr},'_task-MultimodalReadSpeech_space-IXI549Space_FWHM-6_conFWHM-0_node-withinGroup_contrast-block'];
            
            %stat_dir = fullfile(opt.dir.stats.group, stat_name);
            
            %Grab a GROUP roi mask
            %RESLICED TO STATS, NOW IN opt.dir.wta, 'roi'!
            mask = fullfile(opt.dir.wta, 'roi', opt.roi.name{iROI});
            
            %Pick a name for output mask
            output_name = ['group-',opt.group{iGr},'_atlas-JuBrain_desc-WTA6', opt.modality{iMod},'_label-',opt.roi.label{iROI},'_dseg.nii'];
            
            %Grab files from current groupxmodality across all condition
            % 1 = word
            % 2 = pseudoword
            % 3 = control
             files = [cellstr(fullfile(opt.dir.GroupStats, [stat_basename,opt.modality{iMod},'Word'], 'con_0001.nii'));
                 cellstr(fullfile(opt.dir.GroupStats, [stat_basename,opt.modality{iMod},'Pseudoword'], 'con_0001.nii'));
                 cellstr(fullfile(opt.dir.GroupStats, [stat_basename,opt.modality{iMod},'Control'], 'con_0001.nii'))];

             
            winner_take_all(files, mask, output_name);

            
            %Move the output 
            movefile(fullfile(opt.dir.GroupStats, [stat_basename,opt.modality{iMod},'Word'],output_name), ...
                fullfile(opt.dir.wta, 'group_vis', output_name));
            
        end
    end
    
end

%% RUN SOME STATS ON THIS! 

%[RHO,PVAL] = corr(x,y,'Type','Spearman');

%STEFANIA'S PAPER: 

%CORRELATION COEFFICIENTS
%Finally, to compare how similar are the topographical selectivity maps in the three groups we followed, 
%for each pair of groups (i.e. 1.SCv-EBa; 2.SCv-SCa; 3.EBa-SCa) these steps: (1) We computed the Spearman’s
%correlation between the topographical selectivity map of each subject from Group one with the averaged 
%selectivity map of Group two and we compute the mean of these values. (2) We computed the Spearman’s 
%correlation between the topographical selectivity map of each subject from Group two with the averaged 
%selectivity map of Group one and we computed the mean of these values. (3) We averaged the two mean 
%values obtained from step 1 and step 2, in order to have one mean value for each group comparison 
%(see the section ‘Statistical analyses’ for details about the assessment of statistical differences).


%P VALUES
% Therefore, to test statistical differences we used a permutation test (10.000 iterations): 
% (4) We randomly permuted the conditions of the vector of each subject from Group 1 and of the mean vector 
% of Group 2 and we computed the correlation (as in Step 1). (5) We randomly permuted the conditions of the vector 
% of each subject from Group 2 and of the mean vector of Group 1 and we computed the correlation (as in Step 2). 
% (6) We averaged the 2 mean values obtained from step 4 and step 5. (7) We repeated these steps 10.000 times to 
% obtain a distribution of correlations simulating the null hypothesis that the two vectors are unrelated 
% (Kriegeskorte et al., 2008b). If the actual correlation falls within the top α ×100% 
% of the simulated null distribution of correlations, the null hypothesis of unrelated vectors can be 
% rejected with a false-positives rate of α.

%PSEUDOCODE


for each modality
for each ROI mask
for each group
for each subject
    
correlate that subject's map with the other group's map
add to temp corr coefficients

end

(after all subjects are done)
average temp corr coefficients

end

(after both groups are done)
average averaged corr coeff between groups
save coefficient in struct

end
end

             
            


