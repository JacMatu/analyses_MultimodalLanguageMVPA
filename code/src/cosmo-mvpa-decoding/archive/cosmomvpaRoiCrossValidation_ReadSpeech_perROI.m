% cosmo MVPA

function cosmomvpaRoiCrossValidation_ReadSpeech(opt)

% be very verbose please
cosmo_warning('off')

%% Loop through subjects

for iSub = 1:numel(opt.subjects)
    
    opt.cosmomvpa.pathData = fullfile(opt.dir.stats, ...
                                      ['sub-' opt.subjects{iSub}], ...
                                      opt.dir.statsTask);
    
                                  
    %Grab the ROIs based on opt.cosmomvpa.ROIlabel!                              
    BIDS = bids.layout(opt.dir.rois, ... 
                       'use_schema', false);

    opt.query = [];
    opt.query.sub = opt.subjects{iSub};
    opt.query.modality = 'roi';
    opt.query.suffix = 'mask';
    opt.query.label = opt.cosmomvpa.ROIlabel;
    opt.query.space = opt.cosmomvpa.space;
    
    opt.cosmomvpa.roiFileNames = bids.query(BIDS, 'data', ...
                                            opt.query);
    
    %LOOP ACROSS ROIs should start here? Now there is ROI label specified in                                    
    
    roiFileNames = opt.cosmomvpa.roiFileNames;
            
            for i = 1:length(roiFileNames)
                
                roiMasks{i} = roiFileNames{i};  %#ok<*AGROW>
                
            end
    
    
    for iMask = 1:length(roiMasks)
                
                % define brain mask
                mask = roiMasks{iMask};
                
                roiName = roiMasks{iMask};
                
                [~, roiName, ~] = fileparts(roiName);
                
                bidspart = split(roiName, "_");
                
                label = split(bidspart(3), '-');
                
                roiName = label{2};
                
                fprintf(['\nMask: ' roiName '\n\n']);
    
    
    opt.cosmomvpa.csvFileName = fullfile(opt.cosmomvpa.pathOutput, ...
        'accuracy', ...
        ['sub-', opt.subjects{iSub}, ...
        '_task-', opt.taskName, ...
        '_label-', opt.cosmomvpa.ROIlabel{iMask}, ...
        '_desc-smth', num2str(opt.cosmomvpa.funcFWHM), ...
        '_ffx-', opt.cosmomvpa.ffxResults{1}, ...
        '_featureRatio-', num2str(opt.cosmomvpa.ratioToKeep), ...
        '_date-', datestr(now, 'yyyymmddHHMM'), ...
        '_mvpa.csv' ]);
    
    matName = fullfile(opt.cosmomvpa.pathOutput, ...
        'permutations', ...
        ['sub-', opt.subjects{iSub}, ...
        '_task-', opt.taskName, ...
        '_label-', opt.cosmomvpa.ROIlabel{iMask}, ...
        '_desc-smth', num2str(opt.cosmomvpa.funcFWHM), ...
        '_ffx-', opt.cosmomvpa.ffxResults{1}, ...
        '_featureRatio-', num2str(opt.cosmomvpa.ratioToKeep), ...
        '_date-', datestr(now, 'yyyymmddHHMM'), ...
        '_mvpa.mat' ]);
    
    
    % set structure array for keeping the results
    accu = struct( ...
        'sub', [], ...
        'roiArea', [], ...
        'roiDimension', [], ...
        'roiNbVoxels', [], ...
        'ffxResults', [], ...
        'conditions', [], ...
        'modality', [], ...
        'accuracy', []);
    
%  ...
%  'predictors', []);
    
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
        
        % set the 4D result image
        resultsImage = fullfile(pathData, ['sub-' opt.subjects{iSub} ...
                                          '_task-' opt.taskName ...
                                          '_space-' opt.cosmomvpa.space ... 
                                          '_desc-4D_' opt.cosmomvpa.ffxResults{iFfxResult} ...
                                          '.nii']);
        
        % loop through ROI dimension and area
        
        for iRoiDimension = 1:length(opt.cosmomvpa.roiDimension)
            
            %                 fprintf(['\n For ROI dimension: ' num2str(opt.cosmomvpa.roiDimension(iRoiDimension)) 'mm\n\n']);
            
            % set ratio to keep depending on the ROI dimension
            opt.cosmomvpa.feature_selection_ratio_to_keep = opt.cosmomvpa.ratioToKeep(iRoiDimension);
           

            
            %% ROI MVPA analyses
            
            
                
                %% 3 conditions (words, pseudo, control) x 2 modalities (reading, speech) in both groups (blind, sighted)
                stim = [1 2 3  ...
                        4 5 6];
                
                labels = opt.cosmomvpa.labels;
                
                conditionsToTest = opt.cosmomvpa.conditions;
                
                for iConditionToTest = 1:length(conditionsToTest)
                    
                    % loop through modalities
                    for iModality = 1:numel(opt.cosmomvpa.modalities)
                        
                        fprintf([' Modlaity: ' opt.cosmomvpa.modalities{iModality} '\n\n']);
                        
                        % define the data structure
                        ds = cosmo_fmri_dataset(resultsImage, 'mask', mask);
                        
                        % getting rid off zeros
                        zero_msk = all(ds.samples == 0, 1);
                        
                        ds = cosmo_slice(ds, ~zero_msk, 2);
                        
                        mask_size = size(ds.samples, 2);
                        
                        % set chunks, targets and labels
                        
                        ds = setTargetsChunksLabels(opt, ds, stim, labels, nbRun);
                        
                        switch opt.cosmomvpa.modalities{iModality}
                                                        
                            case 'reading'
                                
                                switch conditionsToTest{iConditionToTest}
                                    
                                    case 'WordPseudoword'
                                
                                        ds = cosmo_slice(ds, ds.sa.targets == 1 | ds.sa.targets == 2);
                                
                                    case 'WordPseudowordControl'
                                        
                                        ds = cosmo_slice(ds, ds.sa.targets == 1 | ds.sa.targets == 2 | ...
                                            ds.sa.targets == 3);
                                end
                                
                            case 'speech'
                                
                                switch conditionsToTest{iConditionToTest}
                                    
                                    case 'WordPseudoword'
                                
                                        ds = cosmo_slice(ds, ds.sa.targets == 4 | ds.sa.targets == 5);
                                
                                    case 'WordPseudowordControl'
                                        
                                        ds = cosmo_slice(ds, ds.sa.targets == 4 | ds.sa.targets == 5 | ...
                                            ds.sa.targets == 6);
                                end
                                

                                
                        end
                        
                        % remove constant features
                        ds = cosmo_remove_useless_data(ds);
                        
                        % partitioning, for test and training : cross validation
                        partitions = cosmo_nfold_partitioner(ds);
                        
                        % ROI mvpa analysis
                        [pred, accuracy] = cosmo_crossvalidate(ds, ...
                            @cosmo_meta_feature_selection_classifier, ...
                            partitions, opt.cosmomvpa);
                        
                        % Create and save confusion matrices for each subject
                        %ConfMat_folds = cosmo_confusion_matrix(ds.sa.targets,pred);
                        %ConfMat = sum(ConfMat_folds, 3);
                        ConfMat = sum(cosmo_confusion_matrix(ds.sa.targets,pred), 3);
                        
                        ConfMatName = fullfile(opt.cosmomvpa.pathOutput, ...
                            'confusion_matrices',...
                            ['sub-', opt.subjects{iSub}, ...
                            '_task-', opt.taskName, ...
                            '_label-', opt.cosmomvpa.ROIlabel{iMask}, ...
                            '_desc-smth', num2str(opt.cosmomvpa.funcFWHM), ...
                            '_ffx-', opt.cosmomvpa.ffxResults{1}, ...
                            '_featureRatio-', num2str(opt.cosmomvpa.ratioToKeep), ...
                            '_date-', datestr(now, 'yyyymmddHHMM'), ...
                            '_desc-', [opt.cosmomvpa.modalities{iModality},conditionsToTest{iConditionToTest}],...
                            '_mvpaConfMat.mat' ]);
                        
                        
                        save(ConfMatName, 'ConfMat');
                        %% PERMUTATION PART
                        
%                         if opt.mvpa.permutate  == 1
                            % number of iterations
                            nbIter = 100;
                            
                            % allocate space for permuted accuracies
                            acc0 = zeros(nbIter, 1);
                            
                            % make a copy of the dataset
                            ds0 = ds;
                            
                            % for _niter_ iterations, reshuffle the labels and compute accuracy
                            % Use the helper function cosmo_randomize_targets
                            for k = 1:nbIter
                                ds0.sa.targets = cosmo_randomize_targets(ds);
                                
                                [~, acc0(k)] = cosmo_crossvalidate(ds0, ...
                                    @cosmo_meta_feature_selection_classifier, ...
                                    partitions, opt.cosmomvpa);
                                
                            end
                            
                            p = sum(accuracy < acc0) / nbIter;
                            
                            

%                         end
                        
                        %%
                        
                        % store results to be saved
                        accu = storeResults(accu, count, ...
                            subID, ...
                            roiName, ...
                            opt.cosmomvpa.roiDimension(iRoiDimension), ...
                            mask_size, ...
                            opt.cosmomvpa.ffxResults{iFfxResult}, ...
                            conditionsToTest{iConditionToTest}, ...
                            opt.cosmomvpa.modalities{iModality}, ...
                            accuracy, ...
                            acc0, ...
                            p);
                        
                        count = count + 1;
                        
                        fprintf(['  - condition: ' conditionsToTest{iConditionToTest} ', accuracy: ' num2str(accuracy) '\n']);
                        
                        fprintf('   - %d permutations: accuracy=%.3f, p=%.4f\n\n\n', nbIter, accuracy, p);

                        
                    end
                    
                end
                
            
            
            end
        
        end
    
    savefileCsv = opt.cosmomvpa.csvFileName;
    
    writetable(struct2table(accu), savefileCsv)

    save(matName, 'accu');
    end
end

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

end


function  accu = storeResults(accu, count, subID, roiName, ...
    roiDimension, mask_size, ffxResult, conditionName, modality, accuracy, acc0, p)

    % store results
    accu(count).sub = subID;
    accu(count).roiArea = roiName;
    accu(count).roiDimension = roiDimension;
    accu(count).roiNbVoxels = mask_size;
    accu(count).ffxResults = ffxResult;
    accu(count).conditions = conditionName;
    accu(count).modality = modality;
    accu(count).accuracy = accuracy;
    %     accu(count).predictors = pred;
    
    % save permuted accuracies
    accu(count).permutation = acc0';
    accu(count).pValue = p;
    

end
