%% RUN SEARCHLIGHT STEPS IN COSMOMVPA
% Script from Stefania/Ceren/Mohamed, then Iqra/Alice and then Jacek and Marco

clear;
clc;

%this_dir = fileparts(mfilename('fullpath'));
this_dir = '/Volumes/Slim_Reaper/Projects/Language_MVPA/code/src/cosmo-mvpa-searchlight';

addpath(fullfile(this_dir, '..', '..', 'lib', 'bidspm'));

bidspm();

opt.dir.root = fullfile(this_dir, '..', '..', '..');
opt.dir.derivatives = fullfile(opt.dir.root, 'outputs', 'derivatives');
opt.dir.stats = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-stats');
opt.dir.rois = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-roi');
opt.dir.mvpa = fullfile(opt.dir.derivatives, 'cosmo-mvpa');

% in Iqra's is `opt.pathOutput`
opt.dir.mvpaSearchlight = fullfile(opt.dir.derivatives, 'cosmo-mvpa-searchlight');

if ~exist(opt.dir.mvpaSearchlight, 'dir')
    mkdir(opt.dir.mvpaSearchlight)
end

% check that svm works properly
cosmo_check_external('libsvm');

%% SET OPTIONS   

% subjects
%subject_label = {'01', '02', '03', '04', '08', '09', '10' ...
%'11', '12', '13', '14', '15', '17', '18', '19'}; %#ok<*NASGU>

%TEST WITH 2
subject_label = {'blind01', 'sighted01'}; %#ok<*NASGU>

opt.taskName = 'MultimodalReadSpeech';

opt.cosmomvpa.modalities = {'reading', 'speech'};

% define the 4D maps to be used
opt.cosmomvpa.funcFWHM = 2;

% define the exact path for the results from stats
opt.dir.statsTask = ['task-MultimodalReadSpeech_space-IXI549Space_FWHM-',num2str(opt.cosmomvpa.funcFWHM),'_node-mvpa6betas'];


% define a neighborhood with approximately nb voxels or radius size in each searchlight.
opt.cosmomvpa.searchlightVoxelNb = 3; % 100 150 'count', or 3 - 5 with 'radius'
opt.cosmomvpa.sphereType = 'radius'; % 'radius' or 'count'

% set which type of ffx results you want to use
opt.cosmomvpa.ffxResults = {'beta'}; % in Iqra's is opt.mvpa.map4D

% set whole brain or another mask
opt.cosmomvpa.roiSource = 'wholeBrain';

% design info
opt.cosmomvpa.nbTrialRepetition = 1;
opt.cosmomvpa.nbRun = 6;

% cosmo options
opt.cosmomvpa.tool = 'cosmo';

% Use the cosmo_cross_validation_measure and set its parameters
% (classifier and partitions) in a measure_args struct.
opt.cosmomvpa.measure = @cosmo_crossvalidation_measure;

% Define which classifier to use, using a function handle.
% Alternatives are @cosmo_classify_{svm,matlabsvm,libsvm,nn,naive_bayes, lda}
opt.cosmomvpa.classifier = @cosmo_classify_libsvm;
opt.cosmomvpa.classifierName = 'libsvm';



%% STEP 1: RUN SEARCHLIGHT

parfor iSub = 1:length(subject_label)
    
    % perform the searchlight
   % searchlight_audVis(opt, subject_label{iSub});
   searchlight_ReadSpeech(opt, subject_label{iSub});

end

%% STEP 2: SMOOTH SEARCHLIGHT RESULTS

funcFWHM2Level = 4;

parfor iSub = 1:length(subject_label)
    
    % smooth the searchlight results
    smoothSearchlight(opt, subject_label{iSub}, funcFWHM2Level);

end

funcFWHM2Level = 6;

parfor iSub = 1:length(subject_label)
    
    % smooth the searchlight results
    smoothSearchlight(opt, subject_label{iSub}, funcFWHM2Level);

end

funcFWHM2Level = 8;

parfor iSub = 1:length(subject_label)
    
    % smooth the searchlight results
    smoothSearchlight(opt, subject_label{iSub}, funcFWHM2Level);

end

%% STEP 3: CREATE SL RESULTS MAPS

% set a specific condition and modality to create results with, these are used to filter files such as
% `s8_sub-01_ffxResult-beta_condition-lowVsHighFreq_modality-visual_radius-3_date-202401072330_mvpaSearchlight.nii`

funcFWHM2Level = [ 4 6 8 ];

condition = 'lowVsHighFreq'; 


for iFWHM = 1:length(funcFWHM2Level)
    for iMod = 1:length(opt.cosmomvpa.modalities)

        modality = opt.cosmomvpa.modalities{iMod};

        createSearchlightResultsMaps(opt, subject_label, condition, modality, funcFWHM2Level(iFWHM));
    
    end
end

%% STEP 4: CREATE CONTRAST MAPS

% we do binary decoding, so the chance is 50%
opt.cosmomvpa.chanceLevel = 50;

condition = 'lowVsHighFreq'; 

funcFWHM2Level = [ 4 6 8 ];

condition = 'lowVsHighFreq'; 


for iFWHM = 1:length(funcFWHM2Level)
    

    modality = 'visual';

    createContrastMaps(opt, subject_label, condition, modality, funcFWHM2Level)

    modality = 'auditory';

    createContrastMaps(opt, subject_label, condition, modality, funcFWHM2Level)

    modality = 'fullMotionAudVis';

    createContrastMaps(opt, subject_label, condition, modality, funcFWHM2Level)

end

createContrastMaps(opt, subject_label, condition, modality, funcFWHM2Level)

%% STEP 5: RUN RFX

% runRFX(maps, funcFWHM2Level)
