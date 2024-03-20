function info = searchlight_ReadSpeech(opt, subID)

% be very verbose please
cosmo_warning('off')

% get the smoothing parameter for 4D map
funcFWHM = opt.cosmomvpa.funcFWHM;   

opt.cosmomvpa.pathOutput = fullfile(opt.dir.derivatives, 'cosmo-mvpa-searchlight', opt.dir.statsTask, ['sub-' subID]);

if ~exist(opt.cosmomvpa.pathOutput, 'dir')
    mkdir(opt.cosmomvpa.pathOutput)
end

% set cosmo mvpa structure
condLabelNb = [1 2 3 ...
                4 5 6];

%condLabelName = {'visual_sf_high', 'visual_sf_low', 'visual_full_motion', ...
%                  'auditory_sf_high', 'auditory_sf_low', 'auditory_full_motion'};

condLabelName = opt.cosmomvpa.labels;

decodingConditionList = opt.cosmomvpa.conditions; % there is a loop below so it can have multiple conditions! 

decodingConditionLabel = opt.cosmomvpa.conditions; % what is it? 

%% let's get going!

% this can be coded better (eg it is good to have the exact nb of runs for each subject even 
% though it is the same nb repeated)
nbRun = opt.cosmomvpa.nbRun;


% count = 1;

% for iSub = 1:numel(opt.subjects)

  % get FFX path
  % subID = iSub;

  % set structure array for keeping the results
  info = struct( ...
  'subID', [], ...
  'maskPath', [], ...
  'maskVoxNb', [], ...
  'searchlightVoxelNb', [], ...
  'image', [], ...
  'ffxSmooth', [], ...
  'roiSource', [], ...
  'decodingConditions', [], ...
  'imagePath', []);

  fprintf(['\n For subj: ' subID '\n\n']);

  opt.cosmomvpa.pathData = fullfile(opt.dir.stats, ...
                                    ['sub-' subID], ...
                                    opt.dir.statsTask);

  pathData = opt.cosmomvpa.pathData;

  for iFfxResult = 1:length(opt.cosmomvpa.ffxResults)

    fprintf(['\n For ffx result: ' opt.cosmomvpa.ffxResults{iFfxResult} '\n\n']);

    % get stats image to run the searchlight on
    resultsImage = fullfile(pathData, ['sub-' subID ...
                                      '_task-' opt.taskName ...
                                        '_space-' opt.cosmomvpa.space ... 
                                      '_desc-4D_' opt.cosmomvpa.ffxResults{iFfxResult} ...
                                      '.nii']);
        
    maskName = 'mask.nii';
    
    mask = fullfile(pathData, ...
                    maskName);
        
    for iDecodingCondition = 1:length(decodingConditionList)

      decodingCondition = decodingConditionList{iDecodingCondition};

      fprintf(['\n For condition: ' decodingCondition '\n\n']);

      % loop through modalities
      for iModality = 1:numel(opt.cosmomvpa.modalities)

        fprintf(['\n For modality: ' opt.cosmomvpa.modalities{iModality} '\n\n']);
          
        % load cosmo input
        ds = cosmo_fmri_dataset(resultsImage, 'mask', mask);

        % Getting rid off zeros
        zeroMask = all(ds.samples == 0, 1);
        ds = cosmo_slice(ds, ~zeroMask, 2);
          
        % calculate the mask size
        maskVoxel = size(ds.samples, 2);

        % set cosmo structure
        % ds = setCosmoStructure(opt, nbRun, ds, condLabelNb, condLabelName);
        ds = setTargetsChunksLabels(opt, ds, condLabelNb, condLabelName, nbRun);

        % slice the ds according to your targers (COPIED FROM DECODING SCRIPT!)
        switch opt.cosmomvpa.modalities{iModality}
                                                        
          case 'reading'
              
              switch decodingCondition
                  
              %Multiclass (3)
                  case 'WordPseudowordControl'
                      
                      ds = cosmo_slice(ds, ds.sa.targets == 1 | ds.sa.targets == 2 | ...
                          ds.sa.targets == 3);
              %Pairwise (2)        
                  case 'WordPseudoword'
              
                      ds = cosmo_slice(ds, ds.sa.targets == 1 | ds.sa.targets == 2);
                  case 'WordControl' 
                      
                     ds = cosmo_slice(ds, ds.sa.targets == 1 | ds.sa.targets == 3);
                      
                  case 'PseudowordControl'
                     ds = cosmo_slice(ds, ds.sa.targets == 2 | ds.sa.targets == 3);
              end
              
          case 'speech'
              
              switch decodingCondition
                  
              %Multiclass (3)
                  case 'WordPseudowordControl'
                      
                      ds = cosmo_slice(ds, ds.sa.targets == 4 | ds.sa.targets == 5 | ...
                          ds.sa.targets == 6);
              %Pairwise (2)        
                  case 'WordPseudoword'
              
                      ds = cosmo_slice(ds, ds.sa.targets == 4 | ds.sa.targets == 5);
                  case 'WordControl' 
                      
                     ds = cosmo_slice(ds, ds.sa.targets == 4 | ds.sa.targets == 6);
                      
                  case 'PseudowordControl'
                     ds = cosmo_slice(ds, ds.sa.targets == 5 | ds.sa.targets == 6);
              end
              

              
      end
          
        % remove constant features
        ds = cosmo_remove_useless_data(ds);

        % partitioning, for test and training : cross validation
        opt.cosmomvpa.partitions = cosmo_nfold_partitioner(ds);
          
        % define a neightborhood
        nbrhood = cosmo_spherical_neighborhood(ds, ...
                                               opt.cosmomvpa.sphereType, ...
                                               opt.cosmomvpa.searchlightVoxelNb);          

        % cosmo_disp(nbrhood);

        % Run the searchlight
        svm_results = cosmo_searchlight(ds, ...
                                        nbrhood, ...
                                        opt.cosmomvpa.measure, ...
                                        opt.cosmomvpa);

        % Store results to disc
        nameOutputNiftiSL = fullfile(char(opt.cosmomvpa.pathOutput), ...
                           [['sub-', subID], ...
                           '_ffxResult-', opt.cosmomvpa.ffxResults{iFfxResult}, ...
                            '_condition-', decodingCondition, ...
                           '_modality-', opt.cosmomvpa.modalities{iModality}, ...
                           '_', opt.cosmomvpa.sphereType, '-', num2str(opt.cosmomvpa.searchlightVoxelNb), ...
                           '_date-', datestr(now, 'yyyymmddHHMM'), '_mvpaSearchlight.nii']);
                    
        % cosmo_disp(svm_results);
        % cosmo_plot_slices(svm_results)
          
        cosmo_map2fmri(svm_results, nameOutputNiftiSL);      
        
      end

    end

  end

  % store the relevant info
  info.subID = subID;
  info.maskPath = mask;
  info.maskVoxNb = maskVoxel;
  info.decodingConditions = decodingConditionLabel;
  info.modality = opt.cosmomvpa.modalities{iModality};
  info.searchlightVoxelNb = opt.cosmomvpa.searchlightVoxelNb;
  info.image = opt.cosmomvpa.ffxResults{iFfxResult};
  info.ffxSmooth = funcFWHM;
  % info.roiSource = opt.cosmomvpa.roiSource;
  info.imagePath = resultsImage;
    
  % increase the counter and allons y!
  % count = count + 1;

  fprintf(['sub-'  subID  ' done']);
  
  %% save output

  % mat file
  % set output filename
  nameMatFile = fullfile(opt.cosmomvpa.pathOutput, ...
                         [['sub-', subID], ...
                          '_ffxResult-', opt.cosmomvpa.ffxResults{iFfxResult}, ...
                          '_condition-', decodingCondition, ...
                          '_modality-', opt.cosmomvpa.modalities{iModality}, ...
                          '_', opt.cosmomvpa.sphereType, '-', num2str(opt.cosmomvpa.searchlightVoxelNb), ...
                          '_date-', datestr(now, 'yyyymmddHHMM'), '_mvpaSearchlight.mat']);

  save(nameMatFile, 'info');


% end


end

function ds = setTargetsChunksLabels(opt, ds, stim, labels, nbRun)

  % set chunks (runs by trial_type), target (stimulation type - modality),
  % names (stimulation type name)
  trialsPerRun = length(stim) * opt.cosmomvpa.nbTrialRepetition;

  chunks = repmat((1:nbRun)', 1, opt.cosmomvpa.nbTrialRepetition)';
  chunks = chunks(:);
  chunks = repmat(chunks, trialsPerRun, 1);

  targets = repmat(stim', 1, nbRun)';
  targets = targets(:);

  labels = repmat(labels', 1, nbRun)';
  labels = labels(:);

  ds.sa.targets = targets;
  ds.sa.chunks = chunks;
  ds.sa.labels = labels;
  
  % horzcat(ds.sa.chunks, ds.sa.targets); %, char(ds.sa.labels))

end

% function ds = setCosmoStructure(opt,nbRun, ds, condLabelNb, condLabelName)
%   % sets up the target, chunk, labels by stimuli condition labels, runs,
%   % number labels.

%   % design info from opt
% %   nbRun = opt.cosmomvpa.nbRun;
%   betasPerCondition = opt.cosmomvpa.nbTrialRepetition;

%   % chunk (runs), target (condition), labels (condition names)
%   conditionPerRun = length(condLabelNb);
%   betasPerRun = betasPerCondition * conditionPerRun;

%   chunks = repmat((1:nbRun)', 1, betasPerRun);
%   chunks = chunks(:);

%   targets = repmat(condLabelNb', 1, nbRun)';
%   targets = targets(:);
%   targets = repmat(targets, betasPerCondition, 1);

%   condLabelName = repmat(condLabelName', 1, nbRun)';
%   condLabelName = condLabelName(:);
%   condLabelName = repmat(condLabelName, betasPerCondition, 1);

%   % assign our 4D image design into cosmo ds git
%   ds.sa.targets = targets;
%   ds.sa.chunks = chunks;
%   ds.sa.labels = condLabelName;

% end
