function step5RunRFX(maps, funcFWHM2Level)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % maps = 'beta'% 't_maps';
  % funcFWHM = 0 %  2;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
%   initEnv();

  opt = getOptionSearchlight();

  if nargin == 0

    maps = opt.mvpa.map4D{2};

    funcFWHM2Level = 8;

  end

  % get the smoothing of 4D images
  funcFWHM = opt.funcFWHM;

  % make the prefix for SL output files
  prefixSmooth = [spm_get_defaults('smooth.prefix'), num2str(funcFWHM2Level), '_'];

  % choose which mask to be used for RFX
%   brainMask = [fullfile(opt.derivativesDir, 'group', 'meanMask.nii'), ',1'];
  brainMask = '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/BinariseMeanAnat.nii,1';
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
                       
  midFilePattern = ['4D-', char(maps), ...
                    '_', 'AcrossHand_Ext', '_', ...
                    opt.mvpa.sphereType, '-', num2str(opt.mvpa.searchlightVoxelNb),...
                    '*.nii'];                     
%   midFilePattern = ['4D-', maps, '_', opt.mvpa.sphereType, ...
%                     '-', num2str(opt.mvpa.searchlightVoxelNb), '*.nii'];
  
  % make folder to save rfx output
  RFXDir = [resultFolder, ['/RFX_SecondLevelSmooth_', ...
                        prefixSmooth(1:end - 1)], '/', char(maps), '_', num2str(funcFWHM)];

  if ~exist(RFXDir, 'dir')
    mkdir(RFXDir);
  end

  %% get the .nii files
  % read conMap files only
  slNiiFile = dir(fullfile(resultFolder, ['conMap_', prefixSmooth, '*_', midFilePattern]));
  slNiiFile([slNiiFile.isdir]) = [];

  % display the con map files to be used in RFX
  disp('Files:');
  for iFile = 1:length(slNiiFile)
    disp([slNiiFile(iFile).name]);
  end

  % load the files to create contrast maps
  for iSub = 1:numel(opt.subjects)

    % prepare subject full file for spm input
    subName = [slNiiFile(iSub).name, ',1'];
    subFullPath = fullfile(resultFolder, subName);

    % save the path +nii files for creating spm contrast maps
    subjFullfile{iSub, 1} = subFullPath;   %#ok<AGROW>

  end

  %% run RFX

  % prepare the matlab batch
  % factorial design specification
  matlabbatch = [];
  matlabbatch{1}.spm.stats.factorial_design.dir = {RFXDir};
  matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = subjFullfile;
  matlabbatch{1}.spm.stats.factorial_design.cov = ...
      struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
  matlabbatch{1}.spm.stats.factorial_design.multi_cov = ...
      struct('files', {}, 'iCFI', {}, 'iCC', {});
  matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
  matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
  matlabbatch{1}.spm.stats.factorial_design.masking.em = {brainMask};
  matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
  matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
  matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

  % model estimation
  matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = ...
      cfg_dep('Factorial design specification: SPM.mat File', ...
              substruct('.', 'val', '{}', {1}, ...
                        '.', 'val', '{}', {1}, ...
                        '.', 'val', '{}', {1}), ...
              substruct('.', 'spmmat'));
  matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
  matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

  % contrast manager
  matlabbatch{3}.spm.stats.con.spmmat(1) = ...
      cfg_dep('Model estimation: SPM.mat File', ...
              substruct('.', 'val', '{}', {2}, ...
                        '.', 'val', '{}', {1}, ...
                        '.', 'val', '{}', {1}), ...
              substruct('.', 'spmmat'));

  matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'SC';
  matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = 1;
  matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

  matlabbatch{3}.spm.stats.con.delete = 1; % 0

  save([RFXDir '/SecondLevel_matlabbatch.mat'], 'matlabbatch');
  spm_jobman('run', matlabbatch);
  matlabbatch = {};

end