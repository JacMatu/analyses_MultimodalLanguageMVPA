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
subject_label = {''}; % ALL SUBJECTS FROM RAW
%subject_label = {'blind'}; % ONLY BLIND GROUP
%subject_label = {'sighted'}; % ONLY SIGHTED GROUP


task_label = {'MultimodalReadSpeech'}; %1-back task with spoken and written words, pseudowords and sensory control 


  %% smooth with 6 for univariate summaries
  
 bidspm(bids_dir, output_dir, 'subject', ...
        'participant_label', subject_label, ...
        'action', 'smooth', ...
        'task', task_label, ...
        'space', {'individual', 'IXI549Space'},...
        'ignore',{'qa'}, ...
        'fwhm', 6);
