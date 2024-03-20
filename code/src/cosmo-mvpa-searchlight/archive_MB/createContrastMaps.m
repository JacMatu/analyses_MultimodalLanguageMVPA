function createContrastMaps(opt, subject_label, condition, modality, funcFWHM2Level)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3 was producing mean classification accuracy
% and standard deviation map
% Here Step 4 is omitting the change level from the maps and make
% contrast maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maps = 'beta'% 't_maps';
% funcFWHM = 0 %  2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maps = opt.cosmomvpa.ffxResults{1};

chanceLevel = opt.cosmomvpa.chanceLevel;

% make the prefix for SL output files
prefixSmooth = [spm_get_defaults('smooth.prefix'), num2str(funcFWHM2Level), '_'];
  
% get the folder name to pick files from
opt.cosmomvpa.pathOutput = fullfile(opt.dir.derivatives, 'cosmo-mvpa-searchlight', opt.dir.statsTask);

% define where the searchlight files are
resultFolder = opt.cosmomvpa.pathOutput;                     
                       
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

disp('Files:');
for iFile = 1:length(slNiiFile)
  disp([slNiiFile(iFile).name]);
end

% extract total subject number
totalSubNb = numel(subject_label);
fprintf('NumSubjects: %i \n', totalSubNb);

% load individual accuracy maps from each subjecy
iSub = 1;

% loop for each subject
while iSub <= totalSubNb

  % get subject label
  subLabel = subject_label{iSub};

  % load the searchlight map
  temp = load_nii(char(subjFullfile{iSub}));

  % calculate accuracies (max 100)
  accuracies = temp.img * 100;

  % extract the chance level
  accuracies(accuracies > 0) = accuracies(accuracies > 0) - chanceLevel;

  % omit the zeros
  accuracies(accuracies == 0) = nan;

  % assign accuracies into new map
  con = accuracies;

  % init a dummy .nii skeleton
  con_map = temp;

  % delete the .img values
  con_map.img = [];

  % assign con values into .img
  con_map.img = con;

  % save mean accuracy map as nifti map
  save_nii(con_map, ...
            char(fullfile(resultFolder, ...
                          ['sub-' subject_label{iSub}], ...
                          ['conMap_', ...
                          prefixSmooth, ...
                          ['sub-', subLabel, '_'], ...
                          midFilePattern(1:end - 20), ...
                          '_mvpaSearchlight.nii'])));

  % increase subject the counter
  iSub = iSub + 1;

end

end