% (C) Copyright 2019 bidspm developers

%%% Script for statistical univariate analysis with bidspm()
%%% to run, needs bidspm() installed and path added and saved. 
clear;
clc;

%% initialize bidspm() for this matlab session

%bidspm in code/lib/bidspm
addpath(fullfile(pwd, '..', '..','lib', 'bidspm'));

bidspm();

%% set up BIDS folders and path
this_dir = fileparts(mfilename('fullpath'));


root_dir = fullfile(this_dir, '..', '..', '..');
opt.dir.bids = fullfile(root_dir, 'inputs','RAWDATA');
opt.dir.roi = fullfile(root_dir, 'outputs','derivatives', 'bidspm-roi');
opt.dir.roi_group = fullfile(root_dir, 'outputs','derivatives', 'bidspm-roi','group');
opt.dir.stats = fullfile(root_dir, 'outputs','derivatives', 'bidspm-stats');
opt.dir.preproc = fullfile(root_dir, 'outputs','derivatives', 'bidspm-preproc');

%% define & create ROIs with visfatlas

% opt.dir.roi = output_dir;
% opt.roi.atlas = 'visfatlas';
% opt.roi.hemi = {'L'};
% opt.roi.name = {'IOS', 'pOTS'};
% opt.roi.space = {'IXI549Space'};
% opt.verbosity = 0;
% 
% bidsCreateROI(opt);

%% define & create sphere ROIs for VWFA: aVWFA, cVWFA, pVWFA
%clear('opt');

% vwfa.locations = {[-45 -51 -12],...
%     [-45 -57 -12],...
%     [-45 -72 -10],...
%     [-39 -71 -8],...
%     [-42 -58 -10]};
% 
% vwfa.names = {'aVWFA','cVWFA','pVWFA','lexVWFA','perVWFA'};
% 
% dataImage = '/Volumes/Slim_Reaper/Projects/Language_MVPA/outputs/derivatives/bidspm-stats/sub-blind01/task-MultimodalReadSpeech_space-IXI549Space_FWHM-6_node-univSummary/beta_0001.nii';
% 
% opt.sphere.radius = 7;
% 
% 
% for iROI = 1:numel(vwfa.names) 
%     
%     
%     
%     opt.sphere.location = vwfa.locations{iROI};
%     
%     opt.verbosity = 0;
%     opt.saveImg = true;
%     
%         
%     
%     createRoi('sphere', ...
%         opt.sphere,...
%         dataImage,...
%         opt.dir.roi_group,...
%         opt.saveImg);
% 
%     
%     mni_x = num2str(vwfa.locations{iROI}(1));
%     mni_y = num2str(vwfa.locations{iROI}(2));
%     mni_z = num2str(vwfa.locations{iROI}(3));
%     
%     %set up some names for ROIs
%     nativeName = ['label-sphere',num2str(opt.sphere.radius),'xMinus',mni_x(2:end) ,'yMinus', mni_y(2:end),'zMinus',mni_z(2:end),'_mask'];
%     bidslikeName= ['group_','space-IXI549Space','_', ...
%                     'label-',vwfa.names{iROI},'_','radius-',num2str(opt.sphere.radius),'mm','_mask'];
%     
%     %rename ROIs to BIDS-like
%     movefile(fullfile(opt.dir.roi_group, [nativeName,'.nii']), fullfile(opt.dir.roi_group, [bidslikeName,'.nii']),'f')
%     movefile(fullfile(opt.dir.roi_group, [nativeName,'.json']), fullfile(opt.dir.roi_group, [bidslikeName,'.json']),'f')
%     
% end

%% try Subject-specific ROIs to overcome BIDS issues (caution: all the ROIs are identical anyway!)

% vwfa.locations = {[-45 -51 -12],...
%     [-45 -57 -12],...
%     [-45 -72 -10],...
%     [-39 -71 -8],...
%     [-42 -58 -10]};
% 
% vwfa.names = {'aVWFA','cVWFA','pVWFA','lexVWFA','perVWFA'};

vwfa.locations = {
     [-39 -71 -8],...
     [-42 -58 -10]};
 
vwfa.names = {'lexVWFA','perVWFA'};

opt.subjects = {'blind01', 'sighted01'};

opt.sphere.radius = 7;


for iSub = 1:numel(opt.subjects) 
    
    dataImage = ['/Volumes/Slim_Reaper/Projects/Language_MVPA/outputs/derivatives/bidspm-stats/sub-',opt.subjects{iSub},'/task-MultimodalReadSpeech_space-IXI549Space_FWHM-6_node-univSummary/beta_0001.nii'];
    
    for iROI = 1:numel(vwfa.names) 



        opt.sphere.location = vwfa.locations{iROI};

        opt.verbosity = 0;
        opt.saveImg = true;



        createRoi('sphere', ...
            opt.sphere,...
            dataImage,...
            opt.dir.roi_group,...
            opt.saveImg);


        mni_x = num2str(vwfa.locations{iROI}(1));
        mni_y = num2str(vwfa.locations{iROI}(2));
        mni_z = num2str(vwfa.locations{iROI}(3));

        %set up some names for ROIs
        nativeName = ['label-sphere',num2str(opt.sphere.radius),'xMinus',mni_x(2:end) ,'yMinus', mni_y(2:end),'zMinus',mni_z(2:end),'_mask'];
        bidslikeName= ['sub-',opt.subjects{iSub},'_','space-IXI549Space','_', ...
                        'label-',vwfa.names{iROI},'_','radius-',num2str(opt.sphere.radius),'mm','_mask'];

        mkdir(fullfile(opt.dir.roi, ['sub-',opt.subjects{iSub}],'roi'));            
                    
        %rename ROIs to BIDS-like
        movefile(fullfile(opt.dir.roi_group, [nativeName,'.nii']), fullfile(opt.dir.roi,['sub-',opt.subjects{iSub}],'roi', [bidslikeName,'.nii']),'f')
        movefile(fullfile(opt.dir.roi_group, [nativeName,'.json']), fullfile(opt.dir.roi,['sub-',opt.subjects{iSub}],'roi', [bidslikeName,'.json']),'f')

    end
    
end