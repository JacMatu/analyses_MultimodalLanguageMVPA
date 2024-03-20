%% statistical significance in mvpa
% non-parametric technique by combining permutations and bootstrapping

% step1:
% For each subject, the labels of the different conditions (eg. motion_vertical and motion_horizontal) were permuted,
% and the same decoding analysis was performed.
% The previous step was repeated 100 times for each subject.

% DONE in our decoding scripts

% step2:
% A bootstrap procedure was applied in order to obtain a group-level null distribution
% that was representative of the whole group.
% From each subject’s null distribution, one value was randomly chosen (with replacement)
% and averaged across all participants.
% This step was repeated 100,000 times resulting in a group-level null distribution of 100,000 values.

% step3:
% The statistical significance of the MVPA results was estimated by comparing the observed result
% to the group-level null distribution. This was done by calculating the proportion of observations
% in the null distribution that had a classification accuracy higher than the one obtained in the real test.

% step4:
% To account for the multiple comparisons, all p values were corrected using false discovery rate (FDR) correction


% JM comment: try to reuse opt from previous scripts to loop across classifications, ROIs and modalities? 
% also, classification name (number of classes) is stored in 'accu' mat files (should be loaded by this script!)

clc;
clear;

%% build opt file mimicking MVPAs
clear;
clc;

this_dir = '/Volumes/Slim_Reaper/Projects/Language_MVPA/code/src/cosmo-mpva';

addpath(fullfile(this_dir, '..', '..', 'lib', 'bidspm'));

bidspm();

opt.dir.root = fullfile(this_dir, '..', '..', '..');
opt.dir.derivatives = fullfile(opt.dir.root, 'outputs', 'derivatives');
opt.dir.stats = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-stats');
opt.dir.rois = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-roi'); %changed from cpp_roi
opt.dir.mvpa = fullfile(opt.dir.derivatives, 'cosmo-mvpa');

% MVPA SETTINGS

opt.cosmomvpa.funcFWHM = 2; %changed from 0
opt.cosmomvpa.space = 'IXI549Space'; %preproc/normalisation space
opt.cosmomvpa.ROIlabel = {'lexVWFA','perVWFA'}; %changed from 'MT' other labels to test: aVWFA, cVWFA, pVWFA, perVWFA

opt.cosmomvpa.ratioToKeep = 0.8;

opt.cosmomvpa.normalization = 'zscore';

opt.cosmomvpa.child_classifier = @cosmo_classify_libsvm;
%opt.cosmomvpa.child_classifier = @cosmo_classify_lda;

opt.cosmomvpa.feature_selector = @cosmo_anova_feature_selector;

% set which ROI sphere dimension (mm) to use
opt.cosmomvpa.roiDimension = [7]; % specify different sphere radiuses if you have them, e.g. [7,10,15...]

% design info
opt.cosmomvpa.nbTrialRepetition = 1; % ??? 

opt.cosmomvpa.nbRun = 6; %changed from 5

% TASK SETTINGS
opt.dir.statsTask = 'task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas';

opt.cosmomvpa.pathOutput = fullfile(opt.dir.derivatives, 'cosmo-mvpa', opt.dir.statsTask);

opt.taskName = 'MultimodalReadSpeech';

%opt.subjects = {'blind01', 'sighted01'}; %testing one for each group first? 
opt.subjects = {'blind01'}; %testing one for each group first? 
            
opt.cosmomvpa.modalities = {'reading', 'speech'}; %does this have to be compatibile with BIDS models/conditions?

% define conditions to test based on labels
opt.cosmomvpa.conditions = {'WordPseudoword','WordPseudowordControl'};

opt.cosmomvpa.labels = {'reading word', 'reading pseudoword', 'reading control', ...
    'speech word', 'speech pseudoword', 'speech control'};

% set which type of ffx results you want to use
opt.cosmomvpa.ffxResults = {'beta'}; 


%% set which file, condition and roi label are we testing

% load the .mat file with decoding accuracies

% Consider looping across: 
% - ROI labels
% - ROI dimensions
% - Modalities
% - Decoding conditions

tic

decodTitle = 'spatialFrequencies';

decodingModality = 'visual'; %visual or auditory

roiList = {'leftMTp', 'rightMTp', 'leftMTa', 'rightMTa'};

subList = {'01', '02', '03', '04', '08', '09', '10' ...
    '11', '12', '13', '14', '15', '17', '18', '19'};

nbRowsInAccu = 12; % HARDCODED: NB ROI * NB MODALITIES OR CONDITIONS

im = 'beta'; %'tmap', 'beta'

smooth = 0;

featureRatio = 200;

voxNb = '100';

% number of iterations for group level null distribution
nbIter = 100000;

% Load acuu

accuFile = dir(fullfile(opt.cosmomvpa.pathOutput,'permutations','*.mat'));

%accuFile = dir(['/Volumes/MICK/analisys_high-re_multiSensSpatFreq/outputs/derivatives/cosmo-mvpa/task-audVisMotSpatialFreq_space-MNI152NLin2009cAsym_FWHM-' num2str(smooth) '_node-mvpaBlockAverage/*.mat']);

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
        
        for iRow = 1:12
            % store results
            accuGroup(count).sub = accu(iRow).sub;
            accuGroup(count).roiArea = accu(iRow).roiArea;
            accuGroup(count).roiDimension = accu(iRow).roiDimension;
            accuGroup(count).roiNbVoxels = accu(iRow).roiNbVoxels;
            accuGroup(count).ffxResults = accu(iRow).ffxResults;
            accuGroup(count).conditions = accu(iRow).conditions;
            accuGroup(count).modality = accu(iRow).modality;
            accuGroup(count).accuracy = accu(iRow).accuracy;
            accuGroup(count).permutation = accu(iRow).permutation;
            accuGroup(count).pValue = accu(iRow).pValue;

                count = count + 1;
                
        end
    end

accu = accuGroup;


for iRoi = 1:length(roiList)
    
    roiLabel = roiList(iRoi);
    
    fprintf('roi: %s \n\n', roiLabel{1})
    
    %% STEP 1: DONE
    
    %% STEP 2: create group null distribution
    timeStart = datestr(now,'HH:MM');
    
    subSamp = zeros(length(subList), nbIter);
    
%     L = length(subList);
    
    for iIter = 1:nbIter
        
        for iAccu = 1:length(accu)
            
            for iSub = 1:length(subList)
                
                subID = subList(iSub);
                                
                if strcmp(char({accu(iAccu).sub}.'), char(subID)) == 1
                    
                    if strcmp(char({accu(iAccu).ffxResults}.'), im) == 1 
                    
                    %check if all the parameters and conditions match
                        
                        if strcmp(string({accu(iAccu).modality}.'), decodingModality)==1
                            
                            if strcmp(string({accu(iAccu).roiArea}.'),roiLabel) == 1
                                
                                %read the subject level permutations = accu(iAccu).permutation;
                                %pick one decoding accuracy randomly with replacement
                                
                                if iIter == 1
                                    fprintf('sapmling from: %s sub-%s %s %s \n\n', ...
                                        accu(iAccu).roiArea, ...
                                        accu(iAccu).sub, ...
                                        accu(iAccu).ffxResults, ...
                                        accu(iAccu).modality)
                                end
                                
                                
                                subSamp(iSub, iIter) = datasample(accu(iAccu).permutation, 1);
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
    timeEnd = datestr(now,'HH:MM');
    
    groupNullDistr = mean(subSamp);
    
    disp('step 2 done')
    
    %% STEP 3: check where does the avg accu of the group falls in the group level null ditribution
    % calculate the proportion of values in the group null ditribution which are above the actual decoding
    % accuracy for a one-tailed test. accordingly change for two-tailed test.
    % p = sum(accuracy < acc0) / nbIter; %from Ceren
    % pValue = (sum(Pooled_Perm>accuracy)+1)/(1+NrPermutations); % from Mohamed
    
    subAccu=zeros(length(subList),1);
    
    % subObsPVal=zeros(length(subList),1);
    
    for iAccu=1:length(accu)
        
        for iSub=1:length(subList)
            
            subID=subList(iSub);
            
            if strcmp(char({accu(iAccu).sub}.'),char(subID)) == 1
                
                %check if all the parameters and conditions match
                    
                    if strcmp(string({accu(iAccu).modality}.'),decodingModality)==1
                        
                        if strcmp(string({accu(iAccu).roiArea}.'),roiLabel)==1
                            
                            %read the actual decoding accuracy
                            subAccu(iSub)=[accu(iAccu).accuracy].';
                            
                        end
                    end
                
            end
        end
    end
    
    subObsPVal(iRoi) = sum(mean(subAccu)<groupNullDistr)/nbIter;
    
    disp('step 3 done')

    
end


%% STEP 4: correct the obtained p-value

% function mafdr([vector of pvalues], BHFDR, 'true') % from Stefania
fdrCorrPVal = mafdr(subObsPVal, 'BHFDR', 'true');


%% save the outout
% save the
% decoding type
% decoding condition
% roi List
% group Null distribution
% subObsPVal
% FDR corrected values

% set output folder/name
pathOutput='/Volumes/MICK/analisys_high-re_multiSensSpatFreq/outputs/derivatives/cosmo-mvpa';

savefileMat = fullfile(pathOutput, ...
    ['stats', '_',  decodingModality , '_featureRatio-', num2str(featureRatio), '_smooth-', num2str(smooth), '_ffx-', im, '_', datestr(now, 'yyyymmddHHMM'), '.mat']);

% set structure array for keeping the results
% mvpaStats = struct( ...
%             'decodTitle', [], ...
%             'decodCondition', [], ...
%             'roiList', [], ...
%             'groupNullDis', [], ...
%             'obsPVal', [], ...
%             'fdrCorPVal', []);

%store output
mvpaStats.decodTitle = decodTitle;
mvpaStats.decodCondition = decodingModality;
mvpaStats.roiList = roiList; % this tells the order of corresponding p-values
mvpaStats.groupNullDistr = groupNullDistr;
mvpaStats.obsPVal = subObsPVal;
mvpaStats.fdrCorPVal = fdrCorrPVal;

% mat file
save(savefileMat, 'mvpaStats');

toc