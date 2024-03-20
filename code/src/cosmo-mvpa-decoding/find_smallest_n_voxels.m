%% Find optimal number of voxels per ROI for MVPA feature selection
% For MVPAs equal number of voxels between subjects is more important than
% common space (identical ROI mask). Here, without the independent
% functional localizers I'm using atlas-based ROIs. However, given the
% individual differences in anatomy and normalization, final 'mask.nii' are
% slightly different across subjects. Using an uniformed group mask (based
% on the conjunction of all subject-level masks) is too strict. Therefore,
% one has to "hard-code" the number of voxels per ROI to keep it constant
% across subjects (rather than using opt.cosmomvpa.ratioToKeep 0-1
% (proportion), since it will be different for each subject. 

% An alternative solution is to run cosmo with all the voxels (ratioToKeep
% = 1) in each subject, then look for the smallest number of voxels in each
% ROI (= "the worst subject") and then re-runing the MVPA with these values
% in each ROI. This is what this script will attempt to do. 

%Pre-requisits: cosmomvpa ran with ratioToKeep = 1 for each subject. 

clc;
clear;

%% set which file, condition and roi label are we testing

% load the .mat file with decoding accuracies
tic

%% Try to use similar opt structure as in regular MVPA call? 
this_dir = '/Volumes/Slim_Reaper/Projects/Language_MVPA/code/src/cosmo-mpva';

addpath(fullfile(this_dir, '..', '..', 'lib', 'bidspm'));

bidspm();

opt.dir.root = fullfile(this_dir, '..', '..', '..');
opt.dir.derivatives = fullfile(opt.dir.root, 'outputs', 'derivatives');
opt.dir.stats = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-stats');
opt.dir.rois = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-roi'); 
opt.dir.mvpa = fullfile(opt.dir.derivatives, 'cosmo-mvpa');


opt.cosmomvpa.funcFWHM = 2; %changed from 0
%opt.cosmomvpa.ROIlabel = 'cVWFA';
%opt.cosmomvpa.ROIlist = {};
opt.cosmomvpa.space = 'MNI';
opt.cosmomvpa.ROIlabel = 'cVWFA'; % 'JuBrain' 'visfatlas' 'hcpex' 'cVWFA'
    switch opt.cosmomvpa.ROIlabel
        case 'visfatlas'
            opt.cosmomvpa.ROIlist = {'IOS','pOTS','v1combined'}; %based on visfatlas
        case 'hcpex'
            opt.cosmomvpa.ROIlist = {'Broca','STSa','STSp'}; %based on visfatlas
        case 'JuBrain'
            opt.cosmomvpa.ROIlist = {'Broca', 'FG2', 'FG4', 'MTG', 'V1'};
        case 'cVWFA'
            opt.cosmomvpa.ROIlist = {'cVWFA'};
    end
opt.cosmomvpa.modalities = {'reading', 'speech'}; %does this have to be compatibile with BIDS models/conditions?
opt.cosmomvpa.conditions = {'WordPseudoword','WordControl', 'PseudowordControl'};
opt.cosmomvpa.mvpa = 'crossmodal'; %unimodal / crossmodal
opt.taskName = 'MultimodalReadSpeech';
opt.dir.statsTask = 'task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas';
opt.cosmomvpa.pathOutput = fullfile(opt.dir.derivatives, 'cosmo-mvpa', opt.dir.statsTask, [opt.cosmomvpa.ROIlabel,'_all_voxels'], opt.cosmomvpa.mvpa);



%subject numbers to be pasted with 'sub-' and chosen group, i.e. 'sub-blind01'
subNum = {'01', '02','03','04','05',...
    '06','07','08','09','10','11','12',...
    '13','14','15','16','17','18','19','20'}; 
 
%paste numbers with group
subList = [strcat('blind', subNum), strcat('sighted',subNum)];


%nbRowsInAccu = 12; % HARDCODED: NB ROI * NB MODALITIES * CONDITIONS
nbRowsInAccu = length(opt.cosmomvpa.modalities) * length(opt.cosmomvpa.conditions) * length(opt.cosmomvpa.ROIlist);


%% Load acuu



%should grab 40 files, 20 blind 20 sighted
%accuFile = dir(['/Volumes/Slim_Reaper/Projects/Language_MVPA/outputs/derivatives/cosmo-mvpa/task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas/',decodTitle,'/',mvpa,'/permutations/*.mat']);
accuFile = dir(fullfile(opt.cosmomvpa.pathOutput, 'permutations', '*.mat'));


count = 1;

    
accuGroup = struct( ...
        'sub', [], ...
        'roiArea', [], ...
        'roiDimension', [], ...
        'roiNbVoxels', [], ...
        'ffxResults', [], ...
        'conditions', [], ...
        'modality', [], ...
        'accuracy', [], ...
        'permutation', [], ...
        'pValue', []);
        
    
    for iFile = 1:length(accuFile)
        
        load(fullfile(accuFile(iFile).folder, accuFile(iFile).name));
        
        for iRow = 1:nbRowsInAccu
            % store results
            accuGroup(count).sub = accu(iRow).sub;
            accuGroup(count).roiArea = accu(iRow).roiArea;
            accuGroup(count).roiDimension = accu(iRow).roiDimension;
            accuGroup(count).roiNbVoxels = accu(iRow).roiNbVoxels;
            accuGroup(count).ffxResults = accu(iRow).ffxResults;
            accuGroup(count).conditions = accu(iRow).conditions;
        %accuGroup(count).modality = accu(iRow).modality;
            accuGroup(count).accuracy = accu(iRow).accuracy;
            accuGroup(count).permutation = accu(iRow).permutation;
            accuGroup(count).pValue = accu(iRow).pValue;

                count = count + 1;
                
        end
    end

accu = accuGroup;


