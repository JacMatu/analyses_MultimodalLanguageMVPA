% % % step 0 : to change the dimension of the mask
% % %first initiate CPP_SPM in batchSearchlight.m
% % 
% % % suject to run in each group
% % opt.subjects = {'001','002','003','004','005','006','007','008',...
% %              '009','010','011','014','015','016','017',...
% %              'pil001','pil002','pil004','013','pil005'};%,
% %          
% % for iSub=1%:length(opt.subjects)
% %     subID=char(opt.subjects(iSub));
% %     
% %     %choose the mask for searchlight and reslice it to match the dimensions
% % %     imageToCheck = fullfile(fileparts(mfilename('fullpath')),'..', '..','outputs','derivatives','bidspm-preproc',...
% % %     strcat('sub-',subID),'ses-001','anat',strcat('sub-',subID,'_ses-001_space-IXI549Space_res-r1pt0_desc-brain_','mask.nii'));
% %     imageToCheck = fullfile(fileparts(mfilename('fullpath')),'..', '..','outputs','derivatives','SLmasks',...
% %     strcat('sub-',subID,'_ses-001_space-IXI549Space_res-r1pt0_desc-brain_','mask.nii'));
% %     
% %     % reference image to which the dimensions should be matched
% %     referenceImage='sub-00_ses-001_task-handDown_run-001_space-IXI549Space_desc-preproc_bold.nii';
% % %     referenceImage=fullfile(fileparts(mfilename('fullpath')),'..', '..','outputs','derivatives','bidspm-preproc',...
% % %     strcat('sub-',subID),'ses-001','func',strcat('sub-',subID,'_ses-001_task-handDown_run-001_space-IXI549Space_desc-preproc_bold.nii'));
% %     
% %     mask_fn = resliceRoiImages(referenceImage, imageToCheck);
% % end
% % 
% % for iSub=1%:length(opt.subjects)
% %     subID=char(opt.subjects(iSub));
% %     
% %     inputImage = fullfile(fileparts(mfilename('fullpath')),'..', '..','outputs','derivatives','SLmasks',...
% %     strcat('r','sub-',subID,'_ses-001_space-IXI549Space_res-r1pt0_desc-brain_','mask.nii'));
% %     
% % %     inputImage ='rsub-001_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii';
% %     inputImage = '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/meanAnat.nii';
% %     threshold=0;
% %     outputImage = thresholdToMask(inputImage, threshold);
% % %     outputImage = thresholdToMaskThisFolder(inputImage, threshold);
% % 
% % end    

%% to reslice manually
%do it manually as the functions are not giving nice results
%http://rfmri.org/node/88

%% to binarise
% use imcalc function and do it manually 
%threshold at >=1
