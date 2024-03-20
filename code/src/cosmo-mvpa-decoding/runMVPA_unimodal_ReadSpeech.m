% Adjustting Marco's script to my data

clear;
clc;

this_dir = '/Volumes/Slim_Reaper/Projects/Language_MVPA/code/src/cosmo-mvpa';

addpath(fullfile(this_dir, '..', '..', 'lib', 'bidspm'));

bidspm();

opt.dir.root = fullfile(this_dir, '..', '..', '..');
opt.dir.derivatives = fullfile(opt.dir.root, 'outputs', 'derivatives');
opt.dir.stats = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-stats');
opt.dir.rois = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-roi'); 
opt.dir.mvpa = fullfile(opt.dir.derivatives, 'cosmo-mvpa');

if ~exist(opt.dir.mvpa, 'dir')
    mkdir(opt.dir.mvpa)
end

%% MVPA settings

opt.cosmomvpa.funcFWHM = 2; %changed from 0
opt.cosmomvpa.space = 'IXI549Space'; %preproc/normalisation space
opt.cosmomvpa.ROInature = 'atlas'; %'atlas' or 'sphere'
%opt.cosmomvpa.ROIlabel = 'cVWFA'; %changed from 'MT' other labels to test: aVWFA, cVWFA, pVWFA, perVWFA
%opt.cosmomvpa.ROIlabel = 'LexPerVWFA';
%opt.cosmomvpa.ROIlabel = 'visfatlas';
%opt.cosmomvpa.ROIlabel = 'hcpex';
opt.cosmomvpa.ROIlabel = 'JuBrain';

% set which ROI sphere dimension (mm) to use
%opt.cosmomvpa.roiDimension = [ 7 ]; % specify different sphere radiuses if you have them, e.g. [7,10,15...]
opt.cosmomvpa.roiDimension = [ 10 ];

opt.cosmomvpa.ratioToKeep = 1; %watch out, with this ROI not all people have the same n of voxels
%opt.cosmomvpa.ratioToKeep = [100];

opt.cosmomvpa.normalization = 'zscore';

opt.cosmomvpa.child_classifier = @cosmo_classify_libsvm;
%opt.cosmomvpa.child_classifier = @cosmo_classify_lda;

opt.cosmomvpa.feature_selector = @cosmo_anova_feature_selector;


% design info
opt.cosmomvpa.nbTrialRepetition = 1; % ??? 

opt.cosmomvpa.nbRun = 6; %changed from 5

% decoding info (modalities, conditions to decode, labels)
opt.cosmomvpa.modalities = {'reading', 'speech'}; %does this have to be compatibile with BIDS models/conditions?

% define conditions to test based on labels
opt.cosmomvpa.conditions = {'WordPseudoword','WordControl', 'PseudowordControl'};
%opt.cosmomvpa.conditions = {'WordPseudowordControl'}; %start slowly with one multiclass

opt.cosmomvpa.labels = {'reading word', 'reading pseudoword', 'reading control', ...
    'speech word', 'speech pseudoword', 'speech control'};

%% Setup output dirs for MVPA

opt.taskName = 'MultimodalReadSpeech';
opt.dir.statsTask = 'task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas';

opt.cosmomvpa.pathOutput = fullfile(opt.dir.derivatives, 'cosmo-mvpa', opt.dir.statsTask,opt.cosmomvpa.ROIlabel,'unimodal');

if ~exist(opt.cosmomvpa.pathOutput, 'dir')
    mkdir(opt.cosmomvpa.pathOutput)
end

if ~exist(fullfile(opt.cosmomvpa.pathOutput, 'accuracy'), 'dir')
    mkdir(fullfile(opt.cosmomvpa.pathOutput, 'accuracy'))
end

if ~exist(fullfile(opt.cosmomvpa.pathOutput, 'permutations'), 'dir')
    mkdir(fullfile(opt.cosmomvpa.pathOutput, 'permutations'))
end

if ~exist(fullfile(opt.cosmomvpa.pathOutput, 'confusion_matrices'), 'dir')
    mkdir(fullfile(opt.cosmomvpa.pathOutput, 'confusion_matrices'))
end


%% Pick your subjects

%List the sub numbers
subNum = {'01', '02','03','04','05',...
    '06','07','08','09','10','11','12',...
    '13','14','15','16','17','18','19','20'}; 
%List the groups
group = {'blind','sighted'};

%paste numbers with group to create full subjects list
opt.subjects = {};
for i = 1:length(group)
    opt.subjects = [opt.subjects, strcat(group{i}, subNum)];
end

%% RUN mvpa analyses

% set which type of ffx results you want to use
% Initial tests say beta is better, but maybe plot it against T MAPS also? 
opt.cosmomvpa.ffxResults = {'beta'};%{'beta'}; 

cosmomvpaRoiCrossValidation_ReadSpeech(opt); % ADJUSTED?

% set which type of ffx results you want to use
%opt.cosmomvpa.ffxResults = {'tmap'};
%cosmomvpaRoiCrossValidation_ReadSpeech(opt) % ADJUSTED?


disp('Unimodal MVPA done!')

%% After you are done, don't forget to calculate P values with nonParametric 