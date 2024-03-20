% (C) Copyright 2019 bidspm developers

%%% Script for preprocessing pipeline with bidspm()
%%% Due to SliceTiming errors in DCM files slice timing is skipped
%%% to run, needs bidspm() installed and path added and saved. 

clear;
clc;

%% Paths after moving to external hard drive! 

%% initialize bidspm() for this matlab session
%bidspm in code/lib/bidspm
addpath(fullfile(pwd, '..', '..','lib', 'bidspm'));

bidspm();

%% set up BIDS folders and paths
this_dir = fileparts(mfilename('fullpath'));
root_dir = fullfile(this_dir, '..', '..', '..');

bids_dir = fullfile(root_dir,'inputs','RAWDATA');
output_dir = fullfile(root_dir, 'outputs','derivatives');

%% define subjects and tasks to preprocess
%subject_label = {''}; % ALL SUBJECTS FROM RAW
%subject_label = {'blind'}; % ONLY BLIND GROUP
subject_label = 'sighted'; % ONLY SIGHTED GROUP


task_label = {'MultimodalReadSpeech'}; %1-back task with spoken and written words, pseudowords and sensory control 


%% does preproc without Slice Timing Correction

subject_n = {'01','02','03','04','05','06','07','08','09','10',...
    '11','12','13','14','15','16','17','18','19','20'};

parfor i = 1:length(subject_n)
    bidspm(bids_dir, output_dir, 'subject', ...
        'participant_label', cellstr([subject_label, subject_n{i}]), ...
        'action', 'preprocess', ...
        'task', task_label, ...
        'verbosity', 0, ...
        'ignore', {'slicetiming'}, ... 
        'space', {'individual', 'IXI549Space'}, ...
        'fwhm', 2);
end
  
  %% smooth with 6 for univariate summaries
parfor i = 1:length(subject_n)
     bidspm(bids_dir, output_dir, 'subject', ...
           'participant_label', [subject_label,subject_n(i)], ...
           'action', 'smooth', ...
           'task', task_label, ...
          'space', {'IXI549Space'},...
          'fwhm', 6);
end
