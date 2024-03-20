clear ;
%% set paths
% spm
warning('off');
addpath(genpath('/Users/shahzad/Documents/MATLAB/spm12'));

% cosmo
cosmo = '~/Documents/MATLAB/CoSMoMVPA';
addpath(genpath(cosmo));
cosmo_warning('once');

% libsvm
libsvm = '~/Documents/MATLAB/libsvm';
addpath(genpath(libsvm));
% verify it worked.
cosmo_check_external('libsvm'); % should not give an error

% add cpp repo
% run ../../visTacMotion_fMRI_analysis/src/../lib/CPP_SPM/initCppSpm.m;

% bidspmDir=fullfile(fileparts(mfilename('fullpath')),'..', 'lib','bidspm');
% addpath(fullfile(fileparts(mfilename('fullpath')),'..', 'lib','bidspm'))
% bidspm();

%%% ATTENTION %%%
% With bidspm the same piece of code was not running 
% Therefore, Inititalise the cpp_spm v1.1.4 which is here 
% '/Users/shahzad/Files/fMRI/visTacMotionFoR/code/visTacMotion_ROI/src' and then run this code

cd ('/Users/shahzad/Files/fMRI/visTacMotionFoR/code/visTacMotion_ROI/src')
run '../lib/CPP_SPM/initCppSpm.m'
cd ('/Users/shahzad/Files/fMRI/visTac/fMRI_analysis/code/calcSearchlight')

%% run mvpa 
% load your options
opt = getOptionSearchlight();

% perform the searchlight
info = step1RunSearchlight_noMask(opt);


info = step1RunSearchlight(opt);  %visDir
info = step1RunSearchlightAcrossHandAnat(opt) ; %across hand anat
info = step1RunSearchlightAcrossHandExt(opt) ; %across hand ext
info = step1RunCrossModalSearchlight(opt); %crossmodalExt-both
info = step1RunSearchlightCrossModalAnatBoth(opt); %crossModalAnat-Both
  
% smoothing
funcFWHM2Level = 8;
maps = 'beta';%beta
% condition = 'BodyParts5';
step2SmoothSearchlight(maps, funcFWHM2Level);

condition = 'AcrossHand_Ext';%conditionNameForStep3
step3CreateSLResultsMaps(condition,maps, funcFWHM2Level)

step4CreateContrastMaps(maps, funcFWHM2Level)

step5RunRFX(maps, funcFWHM2Level)
