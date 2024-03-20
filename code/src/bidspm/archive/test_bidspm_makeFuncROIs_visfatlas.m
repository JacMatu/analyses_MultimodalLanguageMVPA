% THis script creates ROIs for MVPA analyses based on the visfatlas
% (Rosenke et al., 2021), additionally intersected with a group mask from
% all subjects in the project. Thus, it results in atlas ROIs that share
% space in every subject. Additionally, results are written in a BIDS
% format, in which each subject has it's own ROI (even though they are
% IDENTICAL). This is sub-optimal, but expected by MVPA scripts, so for now
% it is what it is. 
% This script TRIES to follow BIDSPM workflow, but sometimes crude
% solutions are just easier than trying to reinvent the wheel. But hey, it
% works.
% Workflow: 
% - create a group mask (intersection of all 'mask.nii' from GLMs)
% - extract atlas ROIs from chosen atlas
% - intersect each ROI with the group mask
% - reslice each ROI with subject's beta image and save in BIDS-like
% fashion 
% Jacek Matuszewski, CPP-lab

%% Create / initialize the basics
clear;
clc;



this_dir = fileparts(mfilename('fullpath'));

opt.dir.root = fullfile(this_dir,'..', '..', '..');
bids_dir = fullfile(opt.dir.root,'inputs','RAWDATA');
opt.dir.output = fullfile(opt.dir.root, 'outputs');
opt.dir.preproc = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-preproc');
opt.dir.stats = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-stats');
opt.dir.rois = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-roi');

%bidspm in code/lib/bidspm
addpath(fullfile(pwd, '..', '..','lib', 'bidspm'));

bidspm();


% get options
% this one only for tact
%opt.subjects = {'blind01'};

%opt.subjects = {'blind01', 'sighted01'}; %testing one for each group first? 
%opt.subjects = {'blind01', 'blind02','blind03','blind03','blind04','blind05',...
%    'blind06','blind07','blind08','blind09','blind10','blind11','blind12',...
%    'blind13','blind14','blind15','blind16','blind17','blind18','blind19','blind20'}; 
 
%opt.subjects = {'sighted01', 'sighted02','sighted03','sighted03','sighted04','sighted05',...
%    'sighted06','sighted07','sighted08','sighted09','sighted10','sighted11','sighted12',...
%    'sighted13','sighted14','sighted15','sighted16','sighted17','sighted18','sighted19','sighted20'}; 

opt.subjects = {'blind01', 'blind02','blind03','blind04','blind05',...
    'blind06','blind07','blind08','blind09','blind10','blind11','blind12',...
    'blind13','blind14','blind15','blind16','blind17','blind18','blind19','blind20', ...
    'sighted01', 'sighted02','sighted03','sighted04','sighted05',...
    'sighted06','sighted07','sighted08','sighted09','sighted10','sighted11','sighted12',...
    'sighted13','sighted14','sighted15','sighted16','sighted17','sighted18','sighted19','sighted20'};


%% Create a group mask to mask the atlas ROIs

% Pick a GLM 
task = 'task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas';
% Paste paths to mask.nii from all subjects 
inputs_masks = fullfile(opt.dir.stats, strcat('sub-',opt.subjects), task, 'mask.nii');
%gr_mask_name = 'group_task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas_mask.nii';
gr_mask_name = 'group_task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas_mask.nii';
gr_overlap_name = 'group_task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas-desc_sumoverlap_mask.nii';


    matlabbatch{1}.spm.util.imcalc.input = inputs_masks';
    matlabbatch{1}.spm.util.imcalc.output = gr_mask_name;
    matlabbatch{1}.spm.util.imcalc.outdir = cellstr(fullfile(opt.dir.rois, 'group'));
    matlabbatch{1}.spm.util.imcalc.expression = 'sum(X)==40'; % this has to be equal to the N of subjects 
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 0;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;


    spm('defaults', 'FMRI');
    spm_jobman('run',matlabbatch);
    
    clear('matlabbatch');

    % Overlap between subjects (%)
    matlabbatch{1}.spm.util.imcalc.input = inputs_masks';
    matlabbatch{1}.spm.util.imcalc.output = gr_overlap_name;
    matlabbatch{1}.spm.util.imcalc.outdir = cellstr(fullfile(opt.dir.rois, 'group'));
    %matlabbatch{1}.spm.util.imcalc.expression = '(sum(X)./40).*100'; % If you want to have % Overlap 
    matlabbatch{1}.spm.util.imcalc.expression = 'sum(X)'; % If you want to have N subject/voxel
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;


    spm('defaults', 'FMRI');
    spm_jobman('run',matlabbatch);
    
    clear('matlabbatch');
%% Create atlas ROIs

%THIS CREATES ONLY A GROUP LEVEL MASK

    % Build character-selective vOTC (L)
     bidspm(bids_dir, opt.dir.output, 'subject', ...
         'participant_label', opt.subjects, ...
         'action', 'create_roi', ...
        'preproc_dir', opt.dir.preproc, ...
        'roi_atlas', 'visfatlas', ...
        'hemisphere', {'L'},...
        'roi_name', {'IOS', 'pOTS'}, ...
        'space', {'IXI549Space'});

    % Build a V1 mask (L + R, ventral + dorsal)
     bidspm(bids_dir, opt.dir.output, 'subject', ...
         'participant_label', opt.subjects, ...
         'action', 'create_roi', ...
        'preproc_dir', opt.dir.preproc, ...
        'roi_atlas', 'visfatlas', ...
        'hemisphere', {'L', 'R'},...
        'roi_name', {'v1v', 'v1d'}, ...
        'space', {'IXI549Space'});

    % paste the V1 mask from sub-parts 
    output_name = 'hemi-bilateral_space-MNI_atlas-visfatlas_label-v1combined_mask.nii';
    
    %Grab all V1 masks, no time for BIDS querry, use oldschool spm_select
    %inputs_v1 = cellstr(spm_select('FPList', fullfile(opt.dir.rois, 'group'), 'hemi.*label-v1.*_mask.nii'));
    inputs_v1 = [cellstr(spm_select('FPList', fullfile(opt.dir.rois, 'group'), 'hemi.*label-v1d.*_mask.nii'));...
        cellstr(spm_select('FPList', fullfile(opt.dir.rois, 'group'), 'hemi.*label-v1v.*_mask.nii'))];
    
    % Run batch creating a binarized sum of all masks 
    matlabbatch{1}.spm.util.imcalc.input = inputs_v1;
    matlabbatch{1}.spm.util.imcalc.output = output_name;
    matlabbatch{1}.spm.util.imcalc.outdir = cellstr(fullfile(opt.dir.rois, 'group'));
    matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3+i4) > 0';  
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;


    spm('defaults', 'FMRI');
    spm_jobman('run',matlabbatch);
    
    clear('matlabbatch');
    %clear sub-parts of V1 to avoid confusion in future ROI-grabbing. 
    %Most likely you'll never split it anyway.
    for i = 1:length(inputs_v1)
      delete(inputs_v1{i})
      delete(strrep(inputs_v1{i}, '.nii','.json'))
    end
    

%% Intersect the ROIs with the group mask to ensure common space for all
%get names of ROIs created with atlas
opt.roi.list = cellstr(spm_select('FPList', fullfile(opt.dir.rois, 'group'), '.*atlas-visfatlas.*_mask.nii'));
opt.roi.names = cellstr(spm_select('List', fullfile(opt.dir.rois, 'group'), '.*atlas-visfatlas.*_mask.nii'));
%run loop to intersect
for i = 1:length(opt.roi.list)
   
    masked_roi_name = [opt.roi.names{i}(1:end-8),'desc-groupmasked_mask.nii'];
    
    matlabbatch{1}.spm.util.imcalc.input = {opt.roi.list{i}; fullfile(opt.dir.rois, 'group',gr_mask_name)};
    matlabbatch{1}.spm.util.imcalc.output = masked_roi_name;
    matlabbatch{1}.spm.util.imcalc.outdir = cellstr(fullfile(opt.dir.rois, 'group'));
    matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';  
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;


    spm('defaults', 'FMRI');
    spm_jobman('run',matlabbatch);
    
    clear('matlabbatch');
    
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  

%% Reslice to betas and save to proper dirs 

% Try to use bidspm for that 
% reslicedImages = resliceRoiImages(referenceImage, imagesToCheck, dryRun)

opt.roi.list = cellstr(spm_select('FPList', fullfile(opt.dir.rois, 'group'), '.*atlas-visfatlas.*desc-groupmasked_mask.nii'));
opt.roi.names = cellstr(spm_select('List', fullfile(opt.dir.rois, 'group'), '.*atlas-visfatlas.*desc-groupmasked_mask.nii'));

for iSub = 1:length(opt.subjects)
   for iROI = 1:length(opt.roi.list)
        % Get subject name in BIDS and the beta
        subName = ['sub-', num2str(opt.subjects{iSub})];
        betaReference = fullfile(opt.dir.stats, subName, task,'beta_0001.nii');


        reslicedImages = resliceRoiImages(betaReference, opt.roi.list{iROI});
        bidslikeName = [subName, '_',opt.roi.names{iROI}];
        movefile(reslicedImages, fullfile(opt.dir.rois, subName,'roi',bidslikeName),'f');

   end
end
