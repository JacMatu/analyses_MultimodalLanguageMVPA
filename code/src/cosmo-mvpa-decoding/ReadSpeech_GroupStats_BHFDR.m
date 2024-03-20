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

    
    
%% SETUP VARIABLES
% Add FDR function from MATLAB FORUM, no mafdr, toolbox licence missing :(
addpath('/Volumes/Slim_Reaper/Projects/Language_MVPA/code/lib/fdr_bh');

decodTitle = 'ReadSpeech_JuBrain';

decodingCondition = {'WordPseudoword', 'WordControl', 'PseudowordControl'};

decodingModality = {'reading', 'speech'};

roiList = {'Broca', 'FG2', 'FG4', 'MTG', 'V1'};

subGroup = {'blind', 'sighted'};


nbRowsInAccu = numel(decodingCondition) * numel(decodingModality) * numel(roiList); 

im = 'beta'; %'tmap', 'beta'

smooth = 2;

% number of iterations for group level null distribution
nbIter = 1000;



%% BH FDR CORRECTION

%Load your uncorrected PValues from the Stelzer script
pathOutput='/Volumes/Slim_Reaper/Projects/Language_MVPA/outputs/derivatives/cosmo-mvpa/task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas/JuBrain/unimodal';
    
nameUncorr = ['Cosmo_stats_unimodal_',decodTitle,'_NIter_',num2str(nbIter),'_uncorr.mat'];
  
load(fullfile(pathOutput, nameUncorr)); 

% Grab all P values to a one big matrix with 1 row per group/modality
    % with 5 rois x 3 conditions per row (roi1 con1 2 3, roi2 con1 2 3 ...)

    puncorr_matrix = [];
    
    for iGr=1:numel(subGroup)
        for iMod=1:numel(decodingModality)
            
            puncorr_matrix = [puncorr_matrix; ...
                struct2array(PvaluesUncorr.(subGroup{iGr}).(decodingModality{iMod}))];
                
        end
    end
    
    
    % Correct the p values with FDR: think about the size of the family :)
    bh_fdr_family_size = 60;
    
    [~, ~, ~,  p_adj60] = fdr_bh(puncorr_matrix);
    
    
    % LOOP TO WRITE FDR CORRECTED VALUES BACK TO A STRUCT SIMILAR TO
    % PValuesUncorr
    
    start_row = 1;
    for iGr=1:numel(subGroup)
        
        for iMod=1:numel(decodingModality)
            
            
            
            for iRoi = 1:numel(roiList)
                
                col = 3* iRoi - 2; %starting column in the vector (previous ROI + 3)
                vec_temp = zeros(1,numel(decodingCondition));
                
                for iCond = 1:numel(decodingCondition)
                    
                    vec_temp(iCond) = p_adj60(start_row, col); %Grab 1 p value/condition for this ROI and add it to a vector
                    col = col+1;
                    
                end
                
                PvaluesBHFDR.(subGroup{iGr}).(decodingModality{iMod}).(roiList{iRoi}) = vec_temp;
                
                %start_col = start_col + 3;
                
            end
            
            start_row = start_row + 1;
            
        end
    end

 %% Save the output   
    %FDR-corr p values
    
    pathOutput='/Volumes/Slim_Reaper/Projects/Language_MVPA/outputs/derivatives/cosmo-mvpa/task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas/JuBrain/unimodal';
    nameOutput = ['Cosmo_stats_unimodal_',decodTitle,'_NIter_',num2str(nbIter),'_uncorr'];
   
    
    save(fullfile(pathOutput, nameOutput), 'PvaluesBHFDR');