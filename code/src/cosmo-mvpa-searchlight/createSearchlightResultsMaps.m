function createSearchlightResultsMaps(opt, subject_label, condition, modality, funcFWHM2Level)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Step 2 was smoothing the SL result maps
  % Here Step 3 is producing mean classification accuracy
  % and standard deviation map
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % maps = 'beta'% 't_maps';
  % funcFWHM = 0 %  2;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  maps = opt.cosmomvpa.ffxResults{1};
  
  % make the prefix for SL output files
  prefixSmooth = [spm_get_defaults('smooth.prefix'), num2str(funcFWHM2Level), '_'];
  
  % set prefix to nothing if 
  if funcFWHM2Level == 0
  
    prefixSmooth = [];
  
  end
  
  % get the folder name to pick files from
  opt.cosmomvpa.pathOutput = fullfile(opt.dir.derivatives, 'cosmo-mvpa-searchlight', opt.dir.statsTask);
  
  resultFolder = opt.cosmomvpa.pathOutput;                   
  
  % s8_sub-01_ffxResult-beta_condition-lowVsHighFreq_modality-visual_radius-3_date-202401072330_mvpaSearchlight
  
  midFilePattern = ['ffxResult-', maps, ...
                    '_condition-', condition, ...
                    '_modality-', modality, ...
                    '_', opt.cosmomvpa.sphereType, '-', num2str(opt.cosmomvpa.searchlightVoxelNb), ...
                    '*mvpaSearchlight.nii'];
    
  % read specific smoothed files only
  slNiiFile = dir(fullfile(resultFolder, 'sub-*', [prefixSmooth, 'sub*_', midFilePattern]));
  slNiiFile([slNiiFile.isdir]) = [];
  
  for iSub = 1:numel(subject_label)
  
    % prepare subject full file for spm input
    subName = slNiiFile(iSub).name;
    subFullPath = fullfile(resultFolder, ['sub-' subject_label{iSub}], subName);
  
    % save the path +nii files for spm smoothing
    subjFullfile{iSub} = subFullPath;   %#ok<AGROW>
  
  end
  
  disp('Working on these files: ');
  for iFile = 1:length(slNiiFile)
    disp([slNiiFile(iFile).name]);
  end
  disp('\n');
  
  numSubjects = numel(subject_label);
  fprintf('NumSubjects: %i \n', numSubjects);
  
  % Empty matrix of 4 dimensions (first 3 dimensions are the brain image,
  % the fourth dimention is the subject number)
  z = [];

  % first subject Number
  subjectNb = 1;
  
  % loop for each subject
  while subjectNb <= numSubjects

    % load the searchlight map
    % fun fact: load_nii does not like the ", so char() is needed... 
    temp = load_nii(char(subjFullfile{subjectNb}));

    fprintf('Loading of Map %.0f finished. \n', subjectNb);

    % multiply by 100 to get number in percent
    k = temp.img * 100;

    % concatenate each subject to the 4th dimension
    z = cat(4, z, k);

    % increase the counter
    subjectNb = subjectNb + 1;
  end
  
    %% Mean Accuracy
    % calcuate mean accuracy maps
    % copy structure of one map
    meanMap = temp;
  
    % erase the values in the image
    meanMap.img = [];
  
    means = [];
  
    % Calculate mean of each voxel across subjects (4th dimension)
    means = mean(z, 4);
  
    meanMap.img = means;
  
    % save mean accuracy map as nifti map
    save_nii(meanMap, ...
             char(fullfile(resultFolder, ...
                      ['AverageAcc_', ...
                       prefixSmooth, ...
                       midFilePattern(1:end - 20), ...
                       '_subNb-', num2str(numSubjects), ...
                       '_mvpaSearchlight.nii'])));
  
    %% Standard deviation maps
    % copy structure of one map
    stdMap = temp;
  
    % Calculate standard deviation of each voxel across subjects
    stds = std(z, [], 4);
  
    stdMap.img = stds;
  
    % save standard deviation map as nifti map
    save_nii(stdMap, char(fullfile(resultFolder, ...
                              ['StdAcc_', ...
                                prefixSmooth, ...
                                midFilePattern(1:end - 20), ...
                                '_subNb-', num2str(numSubjects), ...
                                '_mvpaSearchlight.nii'])));
  
  end