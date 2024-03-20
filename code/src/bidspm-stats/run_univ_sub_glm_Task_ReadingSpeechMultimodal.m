% (C) Copyright 2019 bidspm developers

%%% Script for statistical univariate analysis with bidspm()
%%% to run, needs bidspm() installed and path added and saved. 
clear;
clc;

%% initialize bidspm() for this matlab session

%bidspm in code/lib/bidspm
addpath(fullfile(pwd, '..','..','lib', 'bidspm'));

bidspm();

%% set up BIDS folders and path
this_dir = fileparts(mfilename('fullpath'));

root_dir = fullfile(this_dir, '..', '..', '..');
bids_dir = fullfile(root_dir, 'inputs','RAWDATA');
output_dir = fullfile(root_dir, 'outputs','derivatives');
preproc_dir = fullfile(root_dir, 'outputs','derivatives', 'bidspm-preproc');

%% define sub and number of betas

subject_label = {''}; % or 'sighted' for full sighted group or '' for EVERYONE
%group = ;
model_file = fullfile(this_dir,'models', 'model-UNIV_group_blocks_ReadSpeech.json');

 
%% subject level (use bold with 6 FWHM)

  bidspm(bids_dir, output_dir, 'subject', ...
           'participant_label', subject_label, ...
           'action', 'stats', ...
           'preproc_dir',preproc_dir,...
           'model_file',model_file,...
           'verbosity', 3, ...
           'concatenate', false, ...
           'space', {'IXI549Space'},...
           'fwhm', 6);  
