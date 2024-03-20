function info = step1RunCrossModalSearchlight(opt)

  % Step 1 : running searchlight (SL) 

  % get the smoothing parameter for 4D map
  funcFWHM = opt.funcFWHM;   

  % set output folder/name
  savefileMat = fullfile(opt.pathOutput, ...
                         [opt.taskName, ...
                         '_SL_CrossModalExt',...
                          '_', datestr(now, 'yyyymmddHHMM'), '.mat']);

  %% MVPA options

  % set cosmo mvpa structure
%   condLabelNb = [1 2 3 4 5 6 ];
%   condLabelName = {'handDown_pinkyThumb', 'handDown_fingerWrist', 'handUp_pinkyThumb', 'handUp_fingerWrist', 'visual_horizontal','visual_vertical' };
%   decodingConditionList = {'handDown_pinkyThumb_vs_handDown_fingerWrist', 'handUp_pinkyThumb_vs_handUp_fingerWrist',...
%       'visual_vertical_vs_visual_horizontal'};

  % set cosmo mvpa structure
  condNb = [ 1 2 3 4 5 6];
  condName = {'handDown_pinkyThumb', 'handDown_fingerWrist', 'handUp_pinkyThumb', 'handUp_fingerWrist', 'visual_horizontal','visual_vertical' };
  directionNb = [1 2 2 1 2 1 ];
  directionName = {'vertical', 'horizontal','horizontal','vertical','horizontal','vertical', };
  modalityNb = [1 1 1 1 2 2];
  modalityName = {'tactile','tactile','tactile','tactile','visual','visual'};
  decodingConditionList = {'trainVisual_testTactile','trainTactile_testVision','both'};


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
%         mask = fullfile(fileparts(mfilename('fullpath')),'..', '..','outputs','derivatives','bidspm-preproc',...
%             strcat('sub-',subID),'ses-001','anat',strcat('r','sub-',subID,'_ses-001_space-IXI549Space_res-r1pt0_desc-brainp00pt50_',maskName));
        mask = fullfile(fileparts(mfilename('fullpath')),'..', '..','outputs','derivatives','SLmasks',...
            strcat('r','sub-',subID,'_ses-001_space-IXI549Space_res-r1pt0_desc-brain_',maskName));
%         mask = fullfile(fileparts(mfilename('fullpath')),'..', '..','outputs','derivatives','bidspm-stats',strcat('sub-',subID),'task-handDown_space-IXI549Space_FWHM-2',maskName);
%        mask='/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-001_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii';
        for iCrossDecodType=3%1:3%length(decodingConditionList) 
            decodingCondition=decodingConditionList(iCrossDecodType);
            
            if iCrossDecodType==1
                test=1;
                decodingCondition=decodingConditionList(1);
            elseif iCrossDecodType==2
                test = 2;
                decodingCondition=decodingConditionList(2);
            elseif iCrossDecodType==3
                test = [];
                decodingCondition=decodingConditionList(3);
            end
            
            % load cosmo input
            ds = cosmo_fmri_dataset(image, 'mask', mask);
            
            ds = cosmo_remove_useless_data(ds);

            % Getting rid off zeros
%             zeroMask = all(ds.samples == 0, 1);
%             ds = cosmo_slice(ds, ~zeroMask, 2);
            

            % set cosmo structure
            ds = setCosmoStructure(opt, nbRun,ds, modalityNb, directionNb );
            
            % Demean  every pattern to remove univariate effect differences
            meanPattern = mean(ds.samples,2);  % get the mean for every pattern
            meanPattern = repmat(meanPattern,1,size(ds.samples,2)); % make a matrix with repmat
            ds.samples  = ds.samples - meanPattern; % remove the mean from every every point in each pattern
            
            % Slice the dataset accroding to modality
            modIdx = (ds.sa.modality==1) | (ds.sa.modality==2) ;
            ds = cosmo_slice(ds,modIdx) ;
            
            % Slice the dataset accroding to motion direction  
            %slice according to directions
            ds = cosmo_slice(ds,(ds.sa.targets==1) | (ds.sa.targets==2)) ;
%             ds = cosmo_slice(ds,strcmp(ds.sa.targets,'vertical') | strcmp(ds.sa.labels,'horizontal')) ;
           
%             % slice the ds according to your targers (choose your
%             % train-test conditions
%             if strcmp (decodingCondition, decodingConditionList(1))==1
%                     ds = cosmo_slice(ds, ds.sa.targets == condLabelNb(1) | ds.sa.targets == condLabelNb(2));
%             elseif strcmp (decodingCondition, decodingConditionList(2)) ==1
%                     ds = cosmo_slice(ds, ds.sa.targets == condLabelNb(3) | ds.sa.targets == condLabelNb(4));
%             elseif strcmp (decodingCondition, decodingConditionList(3)) ==1
%                     ds = cosmo_slice(ds, ds.sa.targets == condLabelNb(5) | ds.sa.targets == condLabelNb(6));
%             end
            
            % remove constant features
%             ds = cosmo_remove_useless_data(ds);

            % partitioning, for test and training : cross validation
%             opt.mvpa.partitions = cosmo_nfold_partitioner(ds);
            opt.mvpa.partitions = cosmo_nchoosek_partitioner(ds, 1, 'modality',test);
            
             % calculate the mask size
            maskVoxel = size(ds.samples, 2);
            
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
            savingResultFile = fullfile(char(opt.pathOutput), ...
                [['sub-', subID], ...
                '_4D-', opt.mvpa.map4D{iImage}, ...
                '_', char(decodingCondition), ...
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

function ds = setCosmoStructure(opt, nbRun,ds, modalityNb, directionNb )
  % sets up the target, chunk, labels by stimuli condition labels, runs,
  % number labels.

  % design info from opt
%   nbRun = opt.mvpa.nbRun;
  betasPerCondition = opt.mvpa.nbTrialRepetition;

  chunks =  repmat(1 :(nbRun*betasPerCondition),1,6);
  chunks = chunks(:);
  
  modalityNb =  repmat(modalityNb,(nbRun*betasPerCondition),1);
  modalityNb = modalityNb(:);
  
  directionNb =  repmat(directionNb,(nbRun*betasPerCondition),1);
  directionNb = directionNb(:);
  
  % assign our 4D image design into cosmo ds git
  ds.sa.chunks = chunks;
  ds.sa.modality = modalityNb;
  ds.sa.targets = directionNb;

end