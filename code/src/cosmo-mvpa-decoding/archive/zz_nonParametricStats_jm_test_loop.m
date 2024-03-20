




clc;
clear;

%% set which file, condition and roi label are we testing

% load the .mat file with decoding accuracies
tic


decodTitle = 'visfatlas';

%decodingModality = {'reading','speech'}; %sensory modalities

%decodingConditions = {'WordPseudoword', 'WordControl','PseudowordControl'};

%roiList = {'pOTS', 'IOS','v1combined'}; %based on visfatlas

% TEST PHASE: JUST PICK ONE
decodingModality = {'reading'};
decodingConditions = {'WordPseudoword'};
roiList = {'IOS'};

group = 'sighted'; %blind sighted
mvpa = 'unimodal';
 
%subject numbers to be pasted with 'sub-' and chosen group, i.e. 'sub-blind01'
subList = {'01', '02','03','04','05',...
    '06','07','08','09','10','11','12',...
    '13','14','15','16','17','18','19','20'}; 
 
%paste numbers with group

subList = strcat(group, subList);

%nbRowsInAccu = 12; % HARDCODED: NB ROI * NB MODALITIES * CONDITIONS
nbRowsInAccu = length(decodingModality) * length(decodingConditions) * length(roiList);


im = 'beta'; %'tmap', 'beta'

smooth = 2;

featureRatio = 0.8;

%voxNb = '100';
% number of iterations for group level null distribution
nbIter = 100000;
%nbIter = 100; %just test if p values are different from 0

% Load acuu



%should grab 40 files, 20 blind 20 sighted
accuFile = dir(['/Volumes/Slim_Reaper/Projects/Language_MVPA/outputs/derivatives/cosmo-mvpa/task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas/',decodTitle,'/',mvpa,'/permutations/*.mat']);


%This is "hardcoded", assumes that only 1 MVPA output per subject is in the file, grabs all and slices it in half into blind or sighted, 
% could be coded better with BIDS query/spm_select maybe? 
if strcmp(group,'blind')
    
    accuFile = accuFile(1:20);
elseif strcmp(group,'sighted')
    
    accuFile = accuFile(21:40);

end
count = 1;

    
accuGroup = struct( ...
        'sub', [], ...
        'roiArea', [], ...
        'roiDimension', [], ...
        'roiNbVoxels', [], ...
        'ffxResults', [], ...
        'conditions', [], ...
        'modality', [], ...
        'accuracy', [], ...
        'permutation', [], ...
        'pValue', []);
        
    
    for iFile = 1:length(accuFile)
        
        load(fullfile(accuFile(iFile).folder, accuFile(iFile).name));
        
        for iRow = 1:12
            % store results
            accuGroup(count).sub = accu(iRow).sub;
            accuGroup(count).roiArea = accu(iRow).roiArea;
            accuGroup(count).roiDimension = accu(iRow).roiDimension;
            accuGroup(count).roiNbVoxels = accu(iRow).roiNbVoxels;
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

%Add a modality loop? 

for iMod = 1:length(decodingModality)

    %subObsPVal = zeros(length(roiList), length(decodingConditions)); %prealocate?

    for iRoi = 1:length(roiList)
        
        roiLabel = roiList(iRoi);
        
        fprintf('roi: %s \n\n', roiLabel{1})
        
        for iCond = 1:length(decodingConditions)
            
            fprintf('condition: %s \n\n', decodingConditions{iCond})
            %% STEP 1: DONE

            %% STEP 2: create group null distribution
                timeStart = datestr(now,'HH:MM');

                subSamp = zeros(length(subList), nbIter);

        %     L = length(subList);

            for iIter = 1:nbIter

                for iAccu = 1:length(accu)

                    for iSub = 1:length(subList)

                        subID = subList(iSub);

                        if strcmp(char({accu(iAccu).sub}.'), char(subID)) == 1

                            if strcmp(char({accu(iAccu).ffxResults}.'), im) == 1 

                            %check if all the parameters and conditions match

                                if strcmp(string({accu(iAccu).modality}.'), decodingModality{iMod})==1

                                    if strcmp(string({accu(iAccu).roiArea}.'),roiLabel) == 1
                                        
                                        if strcmp(string({accu(iAccu).conditions}.'),decodingConditions{iCond}) == 1

                                        %read the subject level permutations = accu(iAccu).permutation;
                                        %pick one decoding accuracy randomly with replacement

                                        if iIter == 1
                                            fprintf('sapmling from: %s sub-%s %s %s \n\n', ...
                                                accu(iAccu).roiArea, ...
                                                accu(iAccu).sub, ...
                                                accu(iAccu).ffxResults, ...
                                                accu(iAccu).modality, ...
                                                accu(iAccu).conditions)
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

            timeEnd = datestr(now,'HH:MM');

            groupNullDistr = mean(subSamp);

            disp('step 2 done')

            %% STEP 3: check where does the avg accu of the group falls in the group level null ditribution
            % calculate the proportion of values in the group null ditribution which are above the actual decoding
            % accuracy for a one-tailed test. accordingly change for two-tailed test.
            % p = sum(accuracy < acc0) / nbIter; %from Ceren
            % pValue = (sum(Pooled_Perm>accuracy)+1)/(1+NrPermutations); % from Mohamed

            subAccu=zeros(length(subList),1);

            % subObsPVal=zeros(length(subList),1);

            for iAccu=1:length(accu)

                for iSub=1:length(subList)

                    subID=subList(iSub);

                    if strcmp(char({accu(iAccu).sub}.'),char(subID)) == 1

                        %check if all the parameters and conditions match

                            if strcmp(string({accu(iAccu).modality}.'),decodingModality{iMod})==1

                                if strcmp(string({accu(iAccu).roiArea}.'),roiLabel)==1

                                    %read the actual decoding accuracy
                                    subAccu(iSub)=[accu(iAccu).accuracy].';

                                end
                            end

                    end
                end
            end

           subObsPVal(iRoi,iCond) = sum(mean(subAccu)<groupNullDistr)/nbIter;

            disp('step 3 done')

        end 
    end


    % STEP 4: correct the obtained p-value

    % function mafdr([vector of pvalues], BHFDR, 'true') % from Stefania
    %vectorize subObsPVal

    % P VALUES ARE STORED BY ROI, THEN BY CONDITIONS, 
    % e.g. with 2 ROI 3 CONDITIONS it's a vector of 
    % [ROI1_cond1, ROI1_cond2, ROI1_cond3, ROI2_cond1, ROI2_cond2, ROI2_cond3]
    subObsPValVec = [subObsPVal(1,:), subObsPVal(2,:)];

    fdrCorrPVal = mafdr(subObsPValVec, 'BHFDR', 'true');

    %% save the outout
    % save the
    % decoding type
    % decoding condition
    % roi List
    % group Null distribution
    % subObsPVal
    % FDR corrected values

    % set output folder/name
    %pathOutput='/Volumes/MICK/analisys_high-re_multiSensSpatFreq/outputs/derivatives/cosmo-mvpa';
    pathOutput = ['/Volumes/Slim_Reaper/Projects/Language_MVPA/outputs/derivatives/cosmo-mvpa/task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas/',mvpa,'/stats'];

    if ~exist(pathOutput, 'dir')
        mkdir(pathOutput)
    end

    savefileMat = fullfile(pathOutput, ...
        ['stats', '_', 'group_' group, decodingModality{iMod} , '_featureRatio-', num2str(featureRatio), '_smooth-', num2str(smooth), '_ffx-', im, '_', datestr(now, 'yyyymmddHHMM'), '.mat']);

    % set structure array for keeping the results
    % mvpaStats = struct( ...
    %             'decodTitle', [], ...
    %             'decodCondition', [], ...
    %             'roiList', [], ...
    %             'groupNullDis', [], ...
    %             'obsPVal', [], ...
    %             'fdrCorPVal', []);

    %store output
    mvpaStats.decodTitle = decodTitle;
    mvpaStats.decodModality = decodingModality{iMod};
    mvpaStats.decodConditions = decodingConditions;
    mvpaStats.roiList = roiList; % this tells the order of corresponding p-values
    mvpaStats.groupNullDistr = groupNullDistr;
    mvpaStats.obsPVal = subObsPValVec;
    mvpaStats.fdrCorPVal = fdrCorrPVal;

    % mat file
    save(savefileMat, 'mvpaStats');

end
toc