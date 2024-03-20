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

decodingCondition = {'WordPseudoword', 'WordControl', 'PseudowordControl'};
decodingModality = {'reading', 'speech'};
roiList = {'Broca', 'FG2', 'FG4', 'MTG', 'V1'};
subGroup = {'blind', 'sighted'};

load('/Volumes/Slim_Reaper/Projects/Language_MVPA/outputs/derivatives/cosmo-mvpa/task-MultimodalReadSpeech_space-IXI549Space_FWHM-2_node-mvpa6betas/JuBrain/unimodal/Cosmo_stats_unimodal_ReadSpeech_JuBrain__uncorr.mat');
    %% Tests
   

    p_uncorr = [1:15; 16:30;31:45; 46:60];
    
    
%Grab conditions by 3s 
start_row = 1;
for iGr=1:numel(subGroup)
    
    for iMod=1:numel(decodingModality)
        
        start_col = 1;
        
        for iRoi = 1:numel(roiList)
            
            Ptest.(subGroup{iGr}).(decodingModality{iMod}).(roiList{iRoi}) = p_uncorr(start_row,start_col:(start_col+2));
            
            start_col = start_col + 3;
           
        end
        
        start_row = start_row + 1;
        
    end
end
    

% try to grab conditions in a loop
start_row = 1;
for iGr=1:numel(subGroup)
    
    for iMod=1:numel(decodingModality)
        
        
        
        for iRoi = 1:numel(roiList)
            
            col = 3* iRoi - 2; %starting column in the vector (previous ROI + 3)
            vec_temp = zeros(1,numel(decodingCondition));
            
            for iCond = 1:numel(decodingCondition)
               
                vec_temp(iCond) = p_uncorr(start_row, col); %Grab 1 p value/condition for this ROI and add it to a vector
                col = col+1;
                
            end
            
            PvaluesBHFDR.(subGroup{iGr}).(decodingModality{iMod}).(roiList{iRoi}) = vec_temp;
            
            %start_col = start_col + 3;
           
        end
        
        start_row = start_row + 1;
        
    end
end
    
    
    % Correct them with the most hardcore FDR possible (family of 60)
   % [~, ~, ~,  p_adj] = fdr_bh(p_uncorr);