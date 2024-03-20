% (C) Copyright 2019 CPP BIDS SPM-pipeline developpers

function opt = getOptionSearchlight()
  % returns a structure that contains the options chosen by the user to run
  % decoding with cosmo-mvpa.

  if nargin < 1
    opt = [];
  end

  % suject to run in each group
  opt.subjects = {'001'};
%   %{'001','002','003','004','005','006','007','008',...
%              '009','010','011','014','015','016','017',...
%              'pil001','pil002','pil004','013','pil005'};%,
         
  % Uncomment the lines below to run preprocessing
  % - don't use realign and unwarp
  opt.realign.useUnwarp = true;

  % we stay in native space (that of the T1)
  opt.space ='MNI';

  % The directory where the data are located
  
  opt.pathOutput = fullfile(fileparts(mfilename('fullpath')),'..', '..','outputs','derivatives',...
      'searchlightOutput','searchlightAcrossHandExt');
  
  % multivariate
  opt.model.file = fullfile(fileparts(mfilename('fullpath')), '..', ...
                            'model', 'model-FoR2mvpa_smdl.json');

  % task to analyze
  opt.taskName = 'FoR2';

  opt.parallelize.do = false;
  opt.parallelize.nbWorkers = 1;
  opt.parallelize.killOnExit = true;

  %% DO NOT TOUCH
  opt = checkOptions(opt);
  saveOptions(opt);
  % we cannot save opt with opt.mvpa, it crashes

  %% mvpa options

  % define the 4D maps to be used
  opt.funcFWHM = 2;

% Define a neighborhood with approximately 100 voxels in each searchlight.
  opt.mvpa.searchlightVoxelNb = 3; % 100 150 'count', or 3 - 5 with 'radius'
  opt.mvpa.sphereType = 'radius'; % 'radius' or 'count'
  
  % set which type of ffx results you want to use
  opt.mvpa.map4D = {'beta'};%,'tmap'};%
  
% whole brain or another mask?
  opt.mvpa.roiSource = 'wholeBrain';
  
  % design info
  opt.mvpa.nbRun = 10;
  opt.mvpa.nbTrialRepetition = 1;%

  % cosmo options
  opt.mvpa.tool = 'cosmo';
  
  % Use the cosmo_cross_validation_measure and set its parameters
  % (classifier and partitions) in a measure_args struct.
  opt.mvpa.measure = @cosmo_crossvalidation_measure;

  % Define which classifier to use, using a function handle.
  % Alternatives are @cosmo_classify_{svm,matlabsvm,libsvm,nn,naive_bayes, lda}
  opt.mvpa.classifier = @cosmo_classify_libsvm;
  opt.mvpa.className = 'libsvm';
  
end
