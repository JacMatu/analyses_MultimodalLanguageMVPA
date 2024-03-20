% cosmo MVPA

function cosmomvpaRoiCrossValidation_ReadSpeech_crossmodal(opt)

% be very verbose please
cosmo_warning('off')

%% Loop through subjects

for iSub = 1:numel(opt.subjects)
    
    opt.cosmomvpa.pathData = fullfile(opt.dir.stats, ...
                                      ['sub-' opt.subjects{iSub}], ...
                                      opt.dir.statsTask);
    BIDS = bids.layout(opt.dir.rois, ... 
                       'use_schema', false);

   % querry ROIs based on atlas or spheres                
   if strcmp(opt.cosmomvpa.ROInature, 'atlas')
        opt.query = [];
        opt.query.sub = opt.subjects{iSub};
        opt.query.modality = 'roi';
        opt.query.suffix = 'mask';
        opt.query.atlas = opt.cosmomvpa.ROIlabel;
        opt.query.space = 'MNI';
        opt.query.desc = ''; %take only ROIs without desc-groupmasked!
   elseif strcmp(opt.cosmomvpa.ROInature, 'sphere') 
        opt.query = [];
        opt.query.sub = opt.subjects{iSub};
        opt.query.modality = 'roi';
        opt.query.suffix = 'mask';
        opt.query.radius = [num2str(opt.cosmomvpa.roiDimension),'mm'];
        opt.query.label = opt.cosmomvpa.ROIlabel;
        opt.query.space = opt.cosmomvpa.space;
    end
    opt.cosmomvpa.roiFileNames = bids.query(BIDS, 'data', ...
                                            opt.query);
    
                                        
    opt.cosmomvpa.csvFileName = fullfile(opt.cosmomvpa.pathOutput, ...
        'accuracy', ...
        ['sub-', opt.subjects{iSub}, ...
        '_task-', opt.taskName, ...
        '_label-', opt.cosmomvpa.ROIlabel, ...
        '_desc-smth', num2str(opt.cosmomvpa.funcFWHM), ...
        '_ffx-', opt.cosmomvpa.ffxResults{1}, ...
        '_featureRatio-', 'WorstSubjectPerROI', ...
        '_date-', datestr(now, 'yyyymmddHHMM'), ...
        '_mvpa.csv' ]);
    
    matName = fullfile(opt.cosmomvpa.pathOutput, ...
        'permutations', ...
        ['sub-', opt.subjects{iSub}, ...
        '_task-', opt.taskName, ...
        '_label-', opt.cosmomvpa.ROIlabel, ...
        '_desc-smth', num2str(opt.cosmomvpa.funcFWHM), ...
        '_ffx-', opt.cosmomvpa.ffxResults{1}, ...
        '_featureRatio-', 'WorstSubjectPerROI', ...
        '_date-', datestr(now, 'yyyymmddHHMM'), ...
        '_mvpa.mat' ]);
    
    % set structure array for keeping the results
    accu = struct( ...
        'sub', [], ...
        'roiArea', [], ...
        'roiDimension', [], ...
        'roiNbVoxels', [], ...
        'mvpaFeatures', [], ...
        'ffxResults', [], ...
        'conditions', [], ...
        'accuracy', []);
    
    %   'modality', [], ...
    
    % get FFX path
    subID = opt.subjects{iSub};
    
    % this can be coded better (eg it is good to have the exact nb of runs for each subject even 
    % though it is the same nb repeated)
    nbRun = opt.cosmomvpa.nbRun(1);
    
    pathData = opt.cosmomvpa.pathData;
    
    fprintf(['\n For subj: ' subID '\n\n']);
    
    % loop through ffxResults (beta and tmaps)
    
    count = 1;
    
    for iFfxResult = 1:length(opt.cosmomvpa.ffxResults)
        
        fprintf(['\n For ffx result: ' opt.cosmomvpa.ffxResults{iFfxResult} '\n\n']);
        
%         % set the 4D result image
%         resultsImage = fullfile(pathData, ['sub-' opt.subjects{iSub} ...
%                                           '_task-' 'audVisTactMotSpatialFreq' ...
%                                           '_space-MNI152NLin2009cAsym' ... 
%                                           '_desc-4D_' opt.cosmomvpa.ffxResults{iFfxResult} ...
%                                           '.nii']);
                                      
        resultsImage = fullfile(pathData, ['sub-' opt.subjects{iSub} ...
                                          '_task-' opt.taskName ...
                                          '_space-' opt.cosmomvpa.space ... 
                                          '_desc-4D_' opt.cosmomvpa.ffxResults{iFfxResult} ...
                                          '.nii']);
        
        % loop through ROI dimension and area
        
        for iRoiDimension = 1:length(opt.cosmomvpa.roiDimension)
            
            %                 fprintf(['\n For ROI dimension: ' num2str(opt.cosmomvpa.roiDimension(iRoiDimension)) 'mm\n\n']);
            
            % set ratio to keep depending on the ROI dimension
            %opt.cosmomvpa.feature_selection_ratio_to_keep = opt.cosmomvpa.ratioToKeep(iRoiDimension);
            
            %                 % load ROI masks for a specific dimension
            %                 roiPattern = ['*' num2str(opt.cosmomvpa.roiDimension(iRoiDimension)) 'mm.nii'];
            %
            %                 roiFileNames = dir(fullfile(pathROI, roiPattern));
            
            roiFileNames = opt.cosmomvpa.roiFileNames;
            
            for i = 1:length(roiFileNames)
                
                roiMasks{i} = roiFileNames{i};  %#ok<*AGROW>
                
            end
            
            %% ROI MVPA analyses
            
            for iMask = 1:length(roiMasks)
                
                % define brain mask
                mask = roiMasks{iMask};
                
                roiName = roiMasks{iMask};
                
                [~, roiName, ~] = fileparts(roiName);
                
                bidspart = split(roiName, "_");
                
                 % atlas and spheres have LABEL in different part of name
                if strcmp(opt.cosmomvpa.ROInature, 'atlas')
                    label = split(bidspart(5), '-');
                else
                    label = split(bidspart(3), '-');
                end
                
                roiName = label{2};
                
                fprintf(['\nMask: ' roiName '\n\n']);
                
                
                % Set the number of voxels for this ROI!
                % take into account whether it varies from ROI to ROI or is
                % fixed single number (either ratio or number of voxels!)
                if numel(opt.cosmomvpa.ratioToKeep) > 1
                    
                    opt.cosmomvpa.feature_selection_ratio_to_keep = opt.cosmomvpa.ratioToKeep(iMask);
                    
                elseif numel(opt.cosmomvpa.ratioToKeep) == 1
                    
                    opt.cosmomvpa.feature_selection_ratio_to_keep = opt.cosmomvpa.ratioToKeep;
                    
                end
                
                featureRatio = opt.cosmomvpa.feature_selection_ratio_to_keep;
                
                %% 3 conditions (Words, Pseudowords, Control) by 2 modalities
                stim = [1 2 3  ...
                        1 2 3];
                
                modalityNb = [1 1 1 ...
                              2 2 2];

                conditionsToTest = {'WordPseudoword','WordControl', 'PseudowordControl'};
              
                TrainTestSplit = {'trainRead_testSpeech', ... 
                                    'trainSpeech_testRead'};
   
                                
                for iConditionToTest = 1:length(conditionsToTest)
                    
                    % loop through modalities
                     for iModality = 1:numel(TrainTestSplit)
                        
                        fprintf([' Training Modality: ' opt.cosmomvpa.modalities{iModality} '\n\n']);
                        
                        TrainTest = TrainTestSplit{iModality};
                        
                        % define the data structure
                        ds = cosmo_fmri_dataset(resultsImage, 'mask', mask);
                        
                        % getting rid off zeros
                        zero_msk = all(ds.samples == 0, 1);
                        
                        ds = cosmo_slice(ds, ~zero_msk, 2);
                        
                        mask_size = size(ds.samples, 2);
                        
                        % remove constant features
                        ds = cosmo_remove_useless_data(ds);
                        
                        % Demean  every pattern to remove univariate effect differences
                        meanPattern = mean(ds.samples, 2);  % get the mean for every pattern
                        meanPattern = repmat(meanPattern, 1, size(ds.samples, 2)); % make a matrix with repmat
                        ds.samples  = ds.samples - meanPattern; % remove the mean from every every point in each pattern
                        
                        
                        % set chunks, targets and labels
                        
                        ds = setTargetsChunksLabels(opt, ds, stim, modalityNb, nbRun);
                        
                        %always train on speech test on reading? 
                        %train_modality = 2;
                        %test_modality = 1;
                       
                         switch TrainTestSplit{iModality}
                                                        
                            case 'trainRead_testSpeech'
                                
                                train_modality = 1;
                                
                                test_modality = 2;
                                
                            case 'trainSpeech_testRead'
                                
                                train_modality = 2;
                                
                                test_modality = 1;
   
                        end                                
   
                       
                        
                        % get only the two modalities of interest
                        msk = cosmo_match(ds.sa.modality, [1 2]);
                        ds_sel = cosmo_slice(ds, msk);
                        
                        % avoid condition 3 full motion or the moment
                        msk = cosmo_match(ds_sel.sa.targets, [1 2 3]);
                        ds_sel = cosmo_slice(ds_sel, msk);
                        
                        
                         switch conditionsToTest{iConditionToTest}
                                    
                                %Pairwise (2)        
                                    case 'WordPseudoword'
                                
                                        ds_sel = cosmo_slice(ds_sel, ds_sel.sa.targets == 1 | ds.sa.targets == 2);
                                    case 'WordControl' 
                                        
                                       ds_sel = cosmo_slice(ds_sel, ds_sel.sa.targets == 1 | ds.sa.targets == 3);
                                       
                                    case 'PseudowordControl'
                                        
                                        ds_sel = cosmo_slice(ds_sel, ds_sel.sa.targets == 2 | ds.sa.targets == 3);
                          end     
                                
                        % partitioning, for test and training : cross validation
                        partitions = cosmo_nchoosek_partitioner(ds_sel, 1, 'modality', test_modality);
                        
                        % ROI mvpa analysis
                        
                        % I use @cosmo_meta_feature_selection_classifier
                        % Iqra uses @cosmo_classify_meta_feature_selection
                        % Donnow the difference

                        [pred, accuracy] = cosmo_crossvalidate(ds_sel, ...
                                                               @cosmo_classify_meta_feature_selection, ...
                                                               partitions, opt.cosmomvpa);
                                                           

                        
                    %% Confusion matrices    
                        ConfMat = sum(cosmo_confusion_matrix(ds_sel.sa.targets,pred), 3);
                        
                        ConfMatName = fullfile(opt.cosmomvpa.pathOutput, ...
                            'confusion_matrices',...
                            ['sub-', opt.subjects{iSub}, ...
                            '_task-', opt.taskName, ...
                            '_label-', roiName, ...
                            '_desc-smth', num2str(opt.cosmomvpa.funcFWHM), ...
                            '_ffx-', opt.cosmomvpa.ffxResults{1}, ...
                            '_featureRatio-', 'WorstSubjectPerROI', ...
                            '_date-', datestr(now, 'yyyymmddHHMM'), ...
                            '_desc-', [TrainTestSplit{iModality},'_',conditionsToTest{iConditionToTest}],...
                            '_mvpa_ConfMat.mat' ]);
                        
                        
                        save(ConfMatName, 'ConfMat');
                        %% PERMUTATION PART - ASK ALICE FOR SOLUTION!
                        
%                         if opt.mvpa.permutate  == 1
                            % number of iterations
                            nbIter = 100;
                            
                            % allocate space for permuted accuracies
                            acc0 = zeros(nbIter, 1);
                            
                            % make a copy of the dataset
                            ds0 = ds_sel;
                            
                            % for _niter_ iterations, reshuffle the labels and compute accuracy
                            % Use the helper function cosmo_randomize_targets
                            for k = 1:nbIter
                                
                                % This randomization doesn't work for
                                % crossmodal decoding
                                % ds0.sa.targets = cosmo_randomize_targets(ds_sel);
                                
                                % Use this approach instead: randomize chunk by chunk and separating modalities
                                for iChunk=1:max(ds0.sa.chunks)
                                    for iTestModality = 1:max(ds0.sa.modality)
                                        
                                        ds0.sa.targets(ds0.sa.chunks==iChunk & ds0.sa.modality==iTestModality) = Shuffle(ds0.sa.targets(ds0.sa.chunks==iChunk & ds0.sa.modality==iTestModality));
                                    
                                    end
                                end
                                
                                [~, acc0(k)] = cosmo_crossvalidate(ds0, ...
                                    @cosmo_meta_feature_selection_classifier, ...
                                    partitions, opt.cosmomvpa);
                                
                            end
                            
                            
                            
                            
                            p = sum(accuracy < acc0) / nbIter;     
                       
                            
                            
                                                % store results to be saved
                        accu = storeResults(accu, count, ...
                            subID, ...
                            roiName, ...
                            opt.cosmomvpa.roiDimension(iRoiDimension), ...
                            mask_size, ...
                            featureRatio, ...
                            opt.cosmomvpa.ffxResults{iFfxResult}, ...
                            conditionsToTest{iConditionToTest}, ...
                            TrainTest, ...
                            accuracy, ...
                            acc0, ...
                            p);    
 %%                

                        count = count + 1;
                        
                        fprintf(['  - condition: ' conditionsToTest{iConditionToTest} ', accuracy: ' num2str(accuracy) '\n\n\n']);
                        
                     end
                    
                end
                
            end
            
        end
        
    end
    
    savefileCsv = opt.cosmomvpa.csvFileName;
    
    writetable(struct2table(accu), savefileCsv)
    
    save(matName, 'accu');
   
end

end


function ds = setTargetsChunksLabels(opt, ds, stim, modalityNb, nbRun)

    % set chunks (runs by trial_type), target (stimulation type - modality),
    % names (stimulation type name)
    trialsPerRun = length(stim) * opt.cosmomvpa.nbTrialRepetition;

    chunks = repmat((1:nbRun)', 1, opt.cosmomvpa.nbTrialRepetition)';
    chunks = chunks(:);
    chunks = repmat(chunks, trialsPerRun, 1);

    targets = repmat(stim', 1, nbRun)';
    targets = targets(:);


    modalityNb =  repmat(modalityNb,(nbRun*opt.cosmomvpa.nbTrialRepetition),1);
    modalityNb = modalityNb(:);

    ds.sa.targets = targets;
    ds.sa.chunks = chunks;
    ds.sa.modality = modalityNb;
    
%     horzcat(ds.sa.chunks, ds.sa.targets, ds.sa.modality)


end


function  accu = storeResults(accu, count, subID, roiName, ...
    roiDimension, mask_size,  featureRatio, ffxResult, conditionName,TrainTest,accuracy, acc0, p)
 

    % store results
    accu(count).sub = subID;
    accu(count).roiArea = roiName;
    accu(count).roiDimension = roiDimension;
    accu(count).roiNbVoxels = mask_size;
    accu(count).mvpaFeatures = featureRatio;
    accu(count).ffxResults = ffxResult;
    accu(count).conditions = conditionName;
    accu(count).TrainTest = TrainTest;
    accu(count).accuracy = accuracy;
    accu(count).permutation = acc0';
    accu(count).pValue = p;
    %accu(count).predictors = pred;

end
