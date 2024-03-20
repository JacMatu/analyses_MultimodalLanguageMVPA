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

%  bidspm(bids_dir, output_dir, 'subject', ...
%           'participant_label', subject_label, ...
%         'action', 'stats', ...
%           'preproc_dir',preproc_dir,...
%           'model_file',model_file,...
%           'verbosity', 3, ...
%           'concatenate', false, ...
%           'space', {'IXI549Space'},...
%           'fwhm', 6);  


 %% dataset level (use bold with 6 FWHM)
 
 
 %group level summaries: one sample t test per each contrast
   bidspm(bids_dir, output_dir, 'dataset', ...
         'participant_label', subject_label, ...
         'action', 'stats', ...
         'preproc_dir', preproc_dir, ...
         'model_file', model_file, ...,
         'verbosity',1,...
         'space', {'IXI549Space'},...
         'fwhm', 6,...,
         'concatenate', false)

 % 

 %% RESULTS p = 000.1 unc
 
 
 
 % define the results to be saved as output
%results = defaultResultsStructure();% do we need this and what does it do?
%results.nodeName = 'within_group'; 
% results.name = {'motion','static','motion_gt_static','static_gt_motion'};
%results.name = {'blockReadControl', 'blockReadWord', 'blockReadPseudoword',...
%    'blockSpeechControl', 'blockSpeechWord', 'blockSpeechPseudoword',...
%    'readWordPseudoGtControl', 'speechWordPseudoGtControl'};
%results.MC =  'none';
%results.p = 0.001;
%results.k = 0;
%results.png = true();
%results.csv = true();
%results.binary = false; 
%results.montage.do = false();
%results.threshSpm = true();
%results.nidm = false();
%results.montage.slices = -12:4:60;
%results.montage.orientation = 'axial';
%results.montage.background = struct('suffix', 'T1w', ...
%                                     'desc', 'preproc', ...
%                                     'modality', 'anat');                                                             
% results.montage.background =[];
%opt.results = results;

%bidspm(bids_dir, output_dir, 'dataset', ...
%       'participant_label', subject_label, ...
%       'action', 'results', ...
%       'preproc_dir', preproc_dir, ...
%       'model_file', model_file,...
%       'space', {'IXI549Space'}, ...
%       'fwhm', 6,...
%       'options', opt);
   
%% RESULTS p = 0.05 FWE
 

 % define the results to be saved as output
%results = defaultResultsStructure();% do we need this and what does it do?
%results.nodeName = 'within_group'; 
% results.name can't be empty, empty =/= takes all specified
%results.name = {'blockReadControl', 'blockReadWord', 'blockReadPseudoword',...
%    'blockSpeechControl', 'blockSpeechWord', 'blockSpeechPseudoword',...
%    'readWordPseudoGtControl', 'speechWordPseudoGtControl'};
%results.MC =  'FWE';
%results.p = 0.05;
%results.k = 10;
%results.png = false();
%results.csv = false();
%results.binary = false; 
%results.montage.do = false();
%results.threshSpm = true();
%results.nidm = false();
%results.montage.slices = -12:4:60;
%results.montage.orientation = 'axial';
%results.montage.background = struct('suffix', 'T1w', ...
%                                     'desc', 'preproc', ...
%                                     'modality', 'anat');                                                             
% results.montage.background =[];
%opt.results = results;

%bidspm(bids_dir, output_dir, 'dataset', ...
%       'participant_label', subject_label, ...
%       'action', 'results', ...
%       'preproc_dir', preproc_dir, ...
%       'model_file', model_file,...
%       'space', {'IXI549Space'}, ...
%       'fwhm', 6,...
%       'options', opt);