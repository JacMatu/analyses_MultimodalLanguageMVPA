%% statistical significance in mvpa
% non-parametric technique by combining permutations and bootstrapping

% step1:
% For each subject, the labels of the different conditions (eg. motion_vertical and motion_horizontal) were permuted,
% and the same decoding analysis was performed.
% The previous step was repeated 100 times for each subject.

% DONE in our decoding scripts

% step2:
% A bootstrap procedure was applied in order to obtain a group-level null distribution
% that was representative of the whole group.
% From each subject’s null distribution, one value was randomly chosen (with replacement)
% and averaged across all participants.
% This step was repeated 100,000 times resulting in a group-level null distribution of 100,000 values.

% step3:
% The statistical significance of the MVPA results was estimated by comparing the observed result
% to the group-level null distribution. This was done by calculating the proportion of observations
% in the null distribution that had a classification accuracy higher than the one obtained in the real test.

% step4:
% To account for the multiple comparisons, all p values were corrected using false discovery rate (FDR) correction

clc;
clear;

%% STEP 1: SET THINGS UP

% load the .mat file with decoding accuracies
tic



decodTitle = 'ReadSpeech_JuBrain';

decodingCondition = {'WordPseudoword', 'WordControl', 'PseudowordControl'};

decodingModality = {'reading', 'speech'};

roiList = {'Broca', 'FG2', 'FG4', 'MTG', 'V1'};

subGroup = {'blind', 'sighted'};


nbRowsInAccu = numel(decodingCondition) * numel(decodingModality) * numel(roiList); 

im = 'beta'; %'tmap', 'beta'

smooth = 2;


% number of iterations for group level null distribution
%nbIter = 10000;
nbIter = 10;

% Setup the structure for p values! - MAYBE THAT CAN BE OMMITED, IT's JUST
% PREALOCATION

% for g = 1:numel(subGroup)
%     for m = 1:numel(decodingModality)
%         for r = 1:numel(roiList)
%             for c = 1:numel(decodingCondition)
%                 
%                 Pvalues.(subGroup{g}).(decodingModality{m}).(roiList{r}).(decodingCondition{c}) = [];
%                 
%             end
%         end
%     end
% end



% GROUP LOOP SHOULD START HERE? 
for iGr = 1:numel(subGroup)
    
     %List the sub numbers
    subNum = {'01', '02','03','04','05',...
        '06','07','08','09','10','11','12',...
        '13','14','15','16','17','18','19','20'}; 
    %paste numbers with group to create full subjects list

    subList = strcat(subGroup{iGr}, subNum);
    
    % Load acuu from one group only!
    accuFile = dir(['/Volumes/Slim_Reaper/Projects/Language_MVPA/outputs/derivatives/cosmo-mvpa/task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas/JuBrain/unimodal/permutations/','sub-',subGroup{iGr},'*.mat']);
    
    count = 1;
    
    
    accuGroup = struct( ...
        'sub', [], ...
        'roiArea', [], ...
        'roiDimension', [], ...
        'roiNbVoxels', [], ...
        'mvpaFeatures',[], ...
        'ffxResults', [], ...
        'conditions', [], ...
        'modality', [], ...
        'accuracy', [], ...
        'permutation', [], ...
        'pValue', []);
    
    %Go through all the loaded file and join their content to a group structure
    % Expected size =
    
    for iFile = 1:length(accuFile)
        
        load(fullfile(accuFile(iFile).folder, accuFile(iFile).name));
        
        for iRow = 1:nbRowsInAccu
            % store results
            accuGroup(count).sub = accu(iRow).sub;
            accuGroup(count).roiArea = accu(iRow).roiArea;
            accuGroup(count).roiDimension = accu(iRow).roiDimension;
            accuGroup(count).roiNbVoxels = accu(iRow).roiNbVoxels;
            accuGroup(count).mvpaFeatures = accu(iRow).mvpaFeatures;
            accuGroup(count).ffxResults = accu(iRow).ffxResults;
            accuGroup(count).conditions = accu(iRow).conditions;
            accuGroup(count).modality = accu(iRow).modality;
            accuGroup(count).accuracy = accu(iRow).accuracy;
            accuGroup(count).permutation = accu(iRow).permutation;
            accuGroup(count).pValue = accu(iRow).pValue;
            
            count = count + 1;
            
        end
    end
    
    accu = accuGroup;
    
    
    %MODALITY LOOP SHOULD START HERE?
    
    for iMod = 1:numel(decodingModality)
        
        for iRoi = 1:length(roiList)
            
            for iCond = 1:numel(decodingCondition)
                
                roiLabel = roiList(iRoi);
                
                fprintf('roi: %s \n\n', roiLabel{1})
                
                
                
                %% STEP 2: create group null distribution
                timeStart = datestr(now,'HH:MM');
                
                subSamp = zeros(length(subList), nbIter);
                
                % Now go into the group accu and randomly draw nbIter values from
                % permutations for each subject in that given modality/decoding condition
                
                for iIter = 1:nbIter
                    
                    for iAccu = 1:length(accu)
                        
                        for iSub = 1:length(subList)
                            
                            subID = subList(iSub);
                            
                            if strcmp(char({accu(iAccu).sub}.'), char(subID)) == 1
                                
                                if strcmp(char({accu(iAccu).ffxResults}.'), im) == 1
                                    
                                    %check if all the parameters and conditions match
                                    
                                    if strcmp(string({accu(iAccu).modality}.'), decodingModality{iMod})==1 %CHANGED THIS TO MODALITY
                                        
                                        if strcmp(string({accu(iAccu).conditions}.'), decodingCondition{iCond})==1
                                            
                                            if strcmp(string({accu(iAccu).roiArea}.'),roiLabel) == 1
                                                
                                                %read the subject level permutations = accu(iAccu).permutation;
                                                %pick one decoding accuracy randomly with replacement
                                                
                                                if iIter == 1
                                                    fprintf('sapmling from: %s sub-%s %s %s \n\n', ...
                                                        accu(iAccu).roiArea, ...
                                                        accu(iAccu).sub, ...
                                                        accu(iAccu).ffxResults, ...
                                                        accu(iAccu).modality)
                                                end
                                                
                                                
                                                subSamp(iSub, iIter) = datasample(accu(iAccu).permutation, 1);
                                                
                                            end
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                %subSamp should now have a size of N subjects x N Iterations
                %(e.g. 20 x 10000)
                
                timeEnd = datestr(now,'HH:MM');
                
                % This averages all subjects PER iteration, resulting in a 1 x NbIter
                % distribution
                groupNullDistr = mean(subSamp);
                
                disp('step 2 done')
                
                %% STEP 3: check where does the avg accu of the group falls in the group level null ditribution
                % calculate the proportion of values in the group null ditribution which are above the actual decoding
                % accuracy for a one-tailed test. accordingly change for two-tailed test.
                % p = sum(accuracy < acc0) / nbIter; %from Ceren
                % pValue = (sum(Pooled_Perm>accuracy)+1)/(1+NrPermutations); % from Mohamed
                
                subAccu=zeros(length(subList),1);
                
                %subObsPVal=zeros(length(subList),1);
                
                for iAccu=1:length(accu)
                    
                    for iSub=1:length(subList)
                        
                        subID=subList(iSub);
                        
                        if strcmp(char({accu(iAccu).sub}.'),char(subID)) == 1
                            
                            %check if all the parameters and conditions match
                            
                            if strcmp(string({accu(iAccu).modality}.'),decodingModality{iMod})==1
                                
                                if strcmp(string({accu(iAccu).conditions}.'),decodingCondition{iCond})==1
                                    
                                    if strcmp(string({accu(iAccu).roiArea}.'),roiLabel)==1
                                        
                                        %read the actual decoding accuracy
                                        subAccu(iSub)=[accu(iAccu).accuracy].';
                                        
                                    end
                                end
                            end
                        end
                    end
                end
                
                % Check how many times group average Accuracy (OBSERVED) is lower than
                % the random one?
                
                %subObsPVal(iRoi) = sum(mean(subAccu)<groupNullDistr)/nbIter;
                
                % If you want a detailed struct with each condition bound
                % to a single p value
                %subObsPVal = sum(mean(subAccu)<groupNullDistr)/nbIter;
                %Pvalues.(subGroup{iGr}).(decodingModality{iMod}).(roiList{iRoi}).(decodingCondition{iCond}) = subObsPVal;
                
                subObsPVal(iCond) = sum(mean(subAccu)<groupNullDistr)/nbIter;
                Pvalues.(subGroup{iGr}).(decodingModality{iMod}).(roiList{iRoi}) = subObsPVal;
                
                disp('step 3 done')
                
                
            end
        end
    end
    
end

%% STEP 4: correct the obtained p-values
    %Typically it's easy to use mafrd function, but it requires a
    %bioinformatics matlab toolbox licence: 
    
    % function mafdr([vector of pvalues], BHFDR, 'true') % from Stefania
    %fdrCorrPVal = mafdr(subObsPVal, 'BHFDR', 'true');
    
    % "Free" solution
    % try using fdr_bh function from matlab file exchange
    % https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh
    % added to code/lib!
    % Usage:
    %  >> [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
   
    % OR IF YOU JUST WANT ADJUSTED P VALUES WITH DEFAULTS
    %[~, ~, ~,  p_adj] = fdr_bh(p_vector_or_matrix);
    % Add FDR function from MATLAB FORUM, no mafdr, toolbox licence missing :(
    addpath('/Volumes/Slim_Reaper/Projects/Language_MVPA/code/lib/fdr_bh');
    
    

    


%% STEP 5 save the outout
    
    %TBA: maybe some basic metadata about the decoding? 

    % set output folder/name
    pathOutput='/Volumes/Slim_Reaper/Projects/Language_MVPA/outputs/derivatives/cosmo-mvpa/task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas/JuBrain/unimodal';
    
    %Uncorr p values
    save(fullfile(pathOutput, ['Cosmo_stats_unimodal_',decodTitle,'_', '_uncorr']), 'Pvalues');
    
    %FDR-corr p values
    

    
    

toc