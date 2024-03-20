%% Create / initialize the basics
clear;
clc;



this_dir = fileparts(mfilename('fullpath'));

opt.dir.root = fullfile(this_dir,'..', '..', '..');
opt.dir.stats = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-stats');
opt.dir.rois = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-roi');

%bidspm in code/lib/bidspm
addpath(fullfile(pwd, '..', '..','lib', 'bidspm'));

bidspm();


% get options
opt.saveROI = true;

% suject to run in each group

% these has beta ref for aud vis
% opt.subjects = {'01','02', '03', '04', ...
%                 '08', '09', '10', ...
%                 '11', '12', '13', '14',  ...
%                 '15', '17', '18', '19' };
            
% this one only for tact
%opt.subjects = {'blind01'};

%opt.subjects = {'blind01', 'sighted01'}; %testing one for each group first? 
%opt.subjects = {'blind01', 'blind02','blind03','blind03','blind04','blind05',...
%    'blind06','blind07','blind08','blind09','blind10','blind11','blind12',...
%    'blind13','blind14','blind15','blind16','blind17','blind18','blind19','blind20'}; 
 
opt.subjects = {'sighted01', 'sighted02','sighted03','sighted03','sighted04','sighted05',...
    'sighted06','sighted07','sighted08','sighted09','sighted10','sighted11','sighted12',...
    'sighted13','sighted14','sighted15','sighted16','sighted17','sighted18','sighted19','sighted20'}; 


%opt.roiList = {'lexVWFA', 'perVWFA'}; 
%opt.radius = 7; %mm
%opt.roiPeaks = {
%     [-39 -71 -8],...
%     [-42 -58 -10]};
% Radius of the sphere around the peak



% get individual coordinates
%individual = roi_getMNICoords(opt.subjects);


opt.roiList = {'cVWFA'};
opt.radius = 12;
opt.roiPeaks = {[-41 -57 -16]};

% Get the ROIs (actually just the spheres)
for iSub = 1:length(opt.subjects)

    % Get subject number
    subName = ['sub-', num2str(opt.subjects{iSub})];
    
    %mkdir(fullfile(opt.dir.rois, subName, 'roi'))

    % for each region this subject has
    for iReg = 1:numel(opt.roiList)

        % if the region is defined (vwfa-br for all, but in some is not present)
        if not(isnan(opt.roiList{iReg}))
            % Get the center
            ROI_center = opt.roiPeaks{iReg};

            % Get the name of the roi for filename
            switch iReg
                case 1, regName = 'cVWFA';
                %case 2, regName = 'perVWFA';
            end

            for rad = opt.radius
                % Set up bids-like name
                bidslikeName = [subName,'_','space-IXI549Space','_', ...
                    'label-',regName,'_','radius-',num2str(rad),'mm','_mask'];

                betaReference = fullfile(opt.dir.stats, subName, '/task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas','beta_0001.nii');

                % anatomical mask is brain mask, make sure that the sphere
                % doesn't excedd the limits of what makes sense
                brainMask = fullfile(opt.dir.stats, subName, 'task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas', 'mask.nii'); 

                % specify the sphere characteristics for each of them
                sphereParams = struct;
                sphereParams.location = ROI_center;
                sphereParams.radius = rad;

                % % specify the object to pass: mask + sphere
                specification = struct('mask1', brainMask, ...
                                       'mask2', sphereParams);

                % specify the path for each subject
                outputPath = fullfile(opt.dir.rois, subName , 'roi');

                sphereMask = createRoi('intersection', specification, betaReference, outputPath, opt.saveROI);

                
                path = opt.dir.rois;
                nativeName = ['label-',sphereMask.label,'_mask'];
                movefile(fullfile(path, subName, 'roi', [nativeName,'.nii']), fullfile(path,  subName, 'roi', [bidslikeName,'.nii']),'f')
                movefile(fullfile(path, subName, 'roi', [nativeName,'.json']), fullfile(path,  subName, 'roi', [bidslikeName,'.json']),'f')

                % reslice
                roiPath = fullfile(opt.dir.rois, subName, 'roi', [bidslikeName, '.nii']);
                sphereMask = resliceRoiImages(betaReference, roiPath);
                
                movefile(fullfile(path, subName, 'roi', ['r' bidslikeName,'.nii']), fullfile(path,  subName, 'roi', [bidslikeName,'.nii']),'f')

            end
        end
    end
end

