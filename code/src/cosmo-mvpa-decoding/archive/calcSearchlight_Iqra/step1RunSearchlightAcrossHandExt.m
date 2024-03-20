function info = step1RunSearchlightAcrossHandExt(opt)

  % Step 1 : running searchlight (SL) 

  % get the smoothing parameter for 4D map
  funcFWHM = opt.funcFWHM;   

  % set output folder/name
  savefileMat = fullfile(opt.pathOutputTemp, ...
                         [opt.taskName, ...
                         '_SL_AcrossHand_Ext',...
                          '_', datestr(now, 'yyyymmddHHMM'), '.mat']);

  %% MVPA options

  % set cosmo mvpa structure
  condLabelName = {'handDown_pinkyThumb', 'handDown_fingerWrist', 'handUp_pinkyThumb', 'handUp_fingerWrist', 'visual_horizontal','visual_vertical' };
  decodingConditionList = {'HDPT_HUPT_vs_HDFW_HUFW','HUPT_HDFW_vs_HDPT_HUFW','HDPT_HDFW_vs_HUPT_HUFW'};


  %% let's get going!

  % set structure array for keeping the results
  
  info = struct( ...
    'subID', [], ...
    'maskPath', [], ...
    'maskVoxNb', [], ...
    'searchlightVoxelNb', [], ...
    'image', [], ...
    'ffxSmooth', [], ...
    'roiSource', [], ...
    'decodingCondition', [], ...
    'imagePath', []);

  count = 1;

  for iSub = 1:numel(opt.subjects)

    % get FFX path
    subID = opt.subjects{iSub};
    ffxDir = getFFXdir(subID, funcFWHM, opt);
    if strcmp(subID,'pil001')==1 || strcmp(subID,'pil002')==1 || strcmp(subID,'pil005')==1
        nbRun=9;
    else
        nbRun=10;
    end

    for iImage = 1:length(opt.mvpa.map4D)
        % 4D image
        imageName = ['sub-',num2str(subID),'_task-FoR2_space-IXI549Space_FWHM-', num2str(funcFWHM),'_desc-4D_', opt.mvpa.map4D{iImage},'.nii'];
        image=fullfile(fileparts(mfilename('fullpath')),'..', '..','outputs','derivatives','mvpaVol',strcat('sub-',subID),imageName);
        
        maskName = 'mask.nii';
        mask = fullfile(fileparts(mfilename('fullpath')),'..', '..','outputs','derivatives','SLmasks',...
            strcat('r','sub-',subID,'_ses-001_space-IXI549Space_res-r1pt0_desc-brain_',maskName));

        for iDecodingCondition=2%:length(decodingConditionList) %see the types in decoding conditionlist and the filename
            decodingCondition=decodingConditionList(iDecodingCondition);
            
            if strcmp (decodingCondition, decodingConditionList(1))==1
            % set cosmo mvpa structure
                condLabelNb = [1 2 1 2 5 6 ]; 
                targets_of_interest = [1,2];
            elseif strcmp (decodingCondition, decodingConditionList(2))==1
            % set cosmo mvpa structure
                condLabelNb = [2 1 1 2 5 6 ]; 
                targets_of_interest = [1,2];
            elseif strcmp (decodingCondition, decodingConditionList(3))==1
            % set cosmo mvpa structure
                condLabelNb = [1 1 2 2 5 6 ]; 
                targets_of_interest = [1,2];
            end
            
            % load cosmo input
            ds = cosmo_fmri_dataset(image, 'mask', mask);

            % Getting rid off zeros
            zeroMask = all(ds.samples == 0, 1);
            ds = cosmo_slice(ds, ~zeroMask, 2);
            
             % calculate the mask size
            maskVoxel = size(ds.samples, 2);

            % set cosmo structure
            ds = setCosmoStructure(opt,nbRun, ds, condLabelNb, condLabelName);
            
            % slice the ds according to your targers (choose your
            % train-test conditions
%             ds = cosmo_slice(ds, ds.sa.targets == condLabelNb(1) | ds.sa.targets == condLabelNb(2));
            ds = cosmo_slice(ds, ds.sa.targets == targets_of_interest(1) | ds.sa.targets == targets_of_interest(2));
            
            % remove constant features
            ds = cosmo_remove_useless_data(ds);

            % partitioning, for test and training : cross validation
            opt.mvpa.partitions = cosmo_nfold_partitioner(ds);
            
            % define a neightborhood
            nbrhood = cosmo_spherical_neighborhood(ds, ...
                                                   opt.mvpa.sphereType, ...
                                                   opt.mvpa.searchlightVoxelNb);                                    
            %cosmo_disp(nbrhood);
            
            % Run the searchlight
            svm_results = cosmo_searchlight(ds, ...
                                            nbrhood, ...
                                            opt.mvpa.measure, ...
                                            opt.mvpa);
            % store the relevant info
            info(count).subID = subID;
            info(count).maskPath = mask;
            info(count).maskVoxNb = maskVoxel;
            info(count).decodingConditions = decodingCondition;
            info(count).searchlightVoxelNb = opt.mvpa.searchlightVoxelNb;
            info(count).image = opt.mvpa.map4D{iImage};
            info(count).ffxSmooth = funcFWHM;
%             info(count).roiSource = opt.mvpa.roiSource;
            info(count).imagePath = image;
            
            % increase the counter and allons y!
            count = count + 1;

            fprintf(['Sub'  subID  'done']);
            
           % Store results to disc
            savingResultFile = fullfile(char(opt.pathOutputTemp), ...
                [['sub-', subID], ...
                '_4D-', opt.mvpa.map4D{iImage}, ...
                '_', 'AcrossHand_Ext', ...
                '_', opt.mvpa.sphereType, ...
                '-', num2str(opt.mvpa.searchlightVoxelNb), ...
                '_', datestr(now, 'yyyymmddHHMM'), '.nii']);
                      
%             cosmo_disp(svm_results);
%             cosmo_plot_slices(svm_results)
            
            cosmo_map2fmri(svm_results, savingResultFile);      
        end
    end
  end
  %% save output

  % mat file
  save(savefileMat, 'info');

end

function ds = setCosmoStructure(opt,nbRun, ds, condLabelNb, condLabelName)
  % sets up the target, chunk, labels by stimuli condition labels, runs,
  % number labels.

  % design info from opt
%   nbRun = opt.mvpa.nbRun;
  betasPerCondition = opt.mvpa.nbTrialRepetition;

  % chunk (runs), target (condition), labels (condition names)
  conditionPerRun = length(condLabelNb);
  betasPerRun = betasPerCondition * conditionPerRun;

  chunks = repmat((1:nbRun)', 1, betasPerRun);
  chunks = chunks(:);

  targets = repmat(condLabelNb', 1, nbRun)';
  targets = targets(:);
  targets = repmat(targets, betasPerCondition, 1);

  condLabelName = repmat(condLabelName', 1, nbRun)';
  condLabelName = condLabelName(:);
  condLabelName = repmat(condLabelName, betasPerCondition, 1);

  % assign our 4D image design into cosmo ds git
  ds.sa.targets = targets;
  ds.sa.chunks = chunks;
  ds.sa.labels = condLabelName;

end
