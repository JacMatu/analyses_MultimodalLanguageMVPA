function step3CreateSLResultsMaps(condition,maps, funcFWHM2Level)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Step 2 was smoothing the SL result maps
  % Here Step 3 is producing mean classification accuracy
  % and standard deviation map
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % maps = 'beta'% 't_maps';
  % funcFWHM = 0 %  2;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  opt = getOptionSearchlight();

  if nargin == 0

    maps = opt.mvpa.map4D{2};

    funcFWHM2Level = 8;

  end

  % get the smoothing of 4D images
  funcFWHM = opt.funcFWHM;

  % make the prefix for SL output files
  prefixSmooth = [spm_get_defaults('smooth.prefix'), num2str(funcFWHM2Level), '_'];

  % dummy call for ffxDir
%   ffxDir = getFFXdir(opt.subjects{1}, funcFWHM, opt);  
%   [~, folderName] = fileparts(ffxDir);

  % get the folder name to pick files from
%   resultFolder = fullfile(opt.pathInput,...
%                            [folderName,  '_',opt.mvpa.roiSource, ...
%                            '_', opt.mvpa.sphereType, ...
%                            '-', num2str(opt.mvpa.searchlightVoxelNb), ...
%                            '_classifier-', opt.mvpa.className]);
  resultFolder = opt.pathOutput ;                     
                       
  midFilePattern = ['4D-', maps, ...
                    '_', condition, '_', ...
                    opt.mvpa.sphereType, '-', num2str(opt.mvpa.searchlightVoxelNb),...
                    '*.nii'];
  
  
  % read smoothed files only
  slNiiFile = dir(fullfile(resultFolder, [prefixSmooth, '*_', midFilePattern]));
  slNiiFile([slNiiFile.isdir]) = [];

  for iSub = 1:numel(opt.subjects)

    % prepare subject full file for spm input
    subName = [slNiiFile(iSub).name];
    subFullPath = fullfile(resultFolder, subName);

    % save the path +nii files for spm smoothing
    subjFullfile{iSub} = subFullPath;   %#ok<AGROW>

  end

  disp('Files:');
  for iFile = 1:length(slNiiFile)
    disp([slNiiFile(iFile).name]);
  end

  numSubjects = numel(opt.subjects);
  fprintf('NumSubjects: %i \n', numSubjects);

  % Empty matrix of 4 dimensions (first 3 dimensions are the brain image,
  % the fourth dimention is the subject number)
  z = [];

  % first subject Number
  subjectNb = 1;

  % loop for each subject
  while subjectNb <= numSubjects

    % load the searchlight map
    temp = load_nii(subjFullfile{subjectNb});
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
           fullfile(resultFolder, ...
                    ['AverageAcc_', ...
                     prefixSmooth, ...
                     midFilePattern(1:end - 5), ...
                     '_subNb-', num2str(numSubjects), ...
                     '.nii']));

  %% Standard deviation maps
  % copy structure of one map
  stdMap = temp;

  % Calculate standard deviation of each voxel across subjects
  stds = std(z, [], 4);

  stdMap.img = stds;

  % save standard deviation map as nifti map
  save_nii(stdMap, fullfile(resultFolder, ...
                            ['StdAcc_', ...
                             prefixSmooth, ...
                             midFilePattern(1:end - 5), ...
                             '_subNb-', num2str(numSubjects), ...
                             '.nii']));

end