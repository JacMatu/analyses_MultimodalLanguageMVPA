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

subject_label = 'sighted'; % or 'sighted' for full sighted group or '' for EVERYONE
%n_betas = {'6betas', '18betas'}; %6 (folds by RUN) or 18 (folds by BLOCK)
n_betas = {'6betas'};

 
%% run bidspm 

%opt. pour changer des options
%par exemple si on veut ajouter des options pour les r�sultats :
%opt.results.xxxx (see documentation)

%+ 'options', opt, .... dans bidspm command 

subject_n = {'01','02','03','04','05','06','07','08','09','10',...
    '11','12','13','14','15','16','17','18','19','20'};

parfor i = 1:length(subject_n)

    for iBeta = 1:numel(n_betas)

      model_file = fullfile(this_dir,'models', ['model-MVPA_',n_betas{iBeta},'_blocks_ReadSpeech.json']);

      bidspm(bids_dir, output_dir, 'subject', ...
               'participant_label', cellstr([subject_label, subject_n{i}]), ...
               'action', 'stats', ...
               'preproc_dir',preproc_dir,...
               'model_file',model_file,...
               'verbosity', 0, ...
               'concatenate', true, ...
               'space', {'IXI549Space'},...
               'fwhm', 2);   

    end

end
   
%%   
fprintf('Done');

%%
%'action', 'contrasts', ... %to run only the contrasts and not the
%statistical values

 %% dataset level
 % 
 % opt.results = struct('nodeName',  'dataset_level', ...
 %                      'name', {{'VisMot_gt_VisStat'}}, ...
 %                      'Mask', false, ...
 %                      'MC', 'none', ...
 %                      'p', 0.05, ...
 %                      'k', 10, ...
 %                      'NIDM', true);
 % 
 

