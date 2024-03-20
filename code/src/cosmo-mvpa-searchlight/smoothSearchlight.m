function smoothSearchlight(opt, subID, funcFWHM2Level)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1 was running searchlight analysis
% Here is Step 2, smoothing the SL result maps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maps = 'beta'% 't_maps';
% funcFWHM = 0 %  2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maps = opt.cosmomvpa.ffxResults{1};

% make the prefix for SL output files
prefixSmooth = [spm_get_defaults('smooth.prefix'), num2str(funcFWHM2Level), '_'];

% get the folder name to pick files from
opt.cosmomvpa.pathOutput = fullfile(opt.dir.derivatives, 'cosmo-mvpa-searchlight', opt.dir.statsTask, ['sub-' subID]);

slResultFolder = opt.cosmomvpa.pathOutput;

midFilePattern = ['ffxResult-', maps,'*mvpaSearchlight.nii'];
  
% define where the sl files are
slNiiFile = dir(fullfile(slResultFolder, ['sub*_', midFilePattern]));

searchLightNiftis = cell(size(slNiiFile, 1), 1);

for iFile = 1:size(slNiiFile, 1)

  searchLightNiftis{iFile} = fullfile(slResultFolder, slNiiFile(iFile).name);

end

% get rid of the nans to avoid that outer brain is smoothed together with the brain
fprintf('Converting zero values to nans\n\n');
  
for iFile = 1:size(searchLightNiftis, 1)

  % convert 0 to nan in .nii files
  tmp = load_nii(searchLightNiftis{iFile});
  tmp.img(tmp.img == 0) = nan;
  save_nii(tmp, searchLightNiftis{iFile});

end

%% SMOOTH
% spm batch for smoothing
matlabbatch = [];
matlabbatch{1}.spm.spatial.smooth.data = searchLightNiftis(:);
matlabbatch{1}.spm.spatial.smooth.fwhm = [funcFWHM2Level funcFWHM2Level funcFWHM2Level];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 1;
matlabbatch{1}.spm.spatial.smooth.prefix = prefixSmooth;

spm_jobman('run', matlabbatch);

end
