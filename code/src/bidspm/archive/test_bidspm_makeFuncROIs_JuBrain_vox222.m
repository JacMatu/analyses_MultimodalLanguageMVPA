%Assumes you already have JuBrain ROIs built with Anatomy toolbox and
%renamed to bids-like fashion! 

% Workflow: 
% - querry atlas ROIs from chosen atlas
% - reslice each ROI with subject's beta, intersects with subject's brain mask from GLM and save in BIDS-like
% fashion 
% Jacek Matuszewski, CPP-lab

%% Create / initialize the basics
clear;
clc;



this_dir = fileparts(mfilename('fullpath'));

opt.dir.root = fullfile(this_dir,'..', '..', '..');
bids_dir = fullfile(opt.dir.root,'inputs','RAWDATA');
opt.dir.output = fullfile(opt.dir.root, 'outputs');
opt.dir.preproc = fullfile(opt.dir.root, 'outputs', 'derivatives','test_vox222', 'derivatives','bidspm-preproc');
opt.dir.stats = fullfile(opt.dir.root, 'outputs', 'derivatives', 'test_vox222','derivatives','bidspm-stats');
opt.dir.rois = fullfile(opt.dir.root, 'outputs', 'derivatives', 'test_vox222','derivatives','bidspm-roi');

%bidspm in code/lib/bidspm
addpath(fullfile(pwd, '..', '..','lib', 'bidspm'));

bidspm();


% get options

task = 'task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas';

%subNum = {'01', '02','03','04','05',...
%    '06','07','08','09','10','11','12',...
%    '13','14','15','16','17','18','19','20'}; 

subNum = {'04'};
group = {'sighted'};
opt.subjects = {};
%paste numbers with group to create full subjects list
for i = 1:length(group)
    opt.subjects = [opt.subjects, strcat(group{i}, subNum)];
end
 
%% Reslice to betas and save to proper dirs 

% Try to use bidspm for that 
% reslicedImages = resliceRoiImages(referenceImage, imagesToCheck, dryRun)

opt.roi.list = cellstr(spm_select('FPList', fullfile(opt.dir.rois, 'group'), '.*atlas-JuBrain.*_mask.nii'));
opt.roi.names = cellstr(spm_select('List', fullfile(opt.dir.rois, 'group'), '.*atlas-JuBrain.*_mask.nii'));

for iSub = 1:length(opt.subjects)
   for iROI = 1:length(opt.roi.list)
        % Get subject name in BIDS and the beta
        subName = ['sub-', num2str(opt.subjects{iSub})];
        
        %make sub dir in rois if necessary
        if ~exist(fullfile(opt.dir.rois, subName,'roi'), 'dir')
            mkdir(fullfile(opt.dir.rois, subName,'roi'));
        end
        
        
        betaReference = fullfile(opt.dir.stats, subName, task,'beta_0001.nii');
        outputPath = fullfile(opt.dir.rois, subName , 'roi');
        
        %copyfile(opt.roi.list{iROI}, fullfile(opt.dir.rois, subName,'roi',opt.roi.names{iROI}), 'f');
        
        %this_roi = fullfile(opt.dir.rois, subName,'roi',opt.roi.names{iROI});
        this_roi = opt.roi.list{iROI};
        
        reslicedImages = resliceRoiImages(betaReference, this_roi);
        bidslikeName = [subName, '_',opt.roi.names{iROI}];
        movefile(reslicedImages, fullfile(opt.dir.rois, subName,'roi',bidslikeName),'f');
        
        %Try to mask with subject's 'mask.nii'
        %Try to also intersect with subject's 'mask.nii'?
        brainMask = fullfile(opt.dir.stats, subName, task, 'mask.nii'); 
        
        
        matlabbatch{1}.spm.util.imcalc.input = [cellstr(brainMask); cellstr(fullfile(opt.dir.rois, subName,'roi',bidslikeName))];
        matlabbatch{1}.spm.util.imcalc.output = bidslikeName;
        matlabbatch{1}.spm.util.imcalc.outdir = cellstr(outputPath);
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';  
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 0;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;


        spm('defaults', 'FMRI');
        spm_jobman('run',matlabbatch);
    
        clear('matlabbatch');

   end
end

disp('ROIs done!');