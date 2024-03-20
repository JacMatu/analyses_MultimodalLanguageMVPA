% motion_winner_take_all run cross validation

subLabel = 'pilot001';

derivativesDir = '/Users/barilari/data/V5_high-res_pilot_bimodalMotion_boldClassic/outputs/derivatives';

statsDirs = fullfile(derivativesDir, 'cpp_spm-stats');

dir = fullfile(statsDirs, ...
                'sub-pilot001/task-gratingBimodalMotion_space-individual_FWHM-0_desc-gratingBimodalMotionBOLDSplitConditionsMVPALatest');

MT_mask = fullfile(derivativesDir, ...
                    'cpp_spm-roi/sub-pilot001/ses-001/roi/rsub-pilot001_ses-001_acq-r0p75_hemi-R_space-EPIT1w_label-rightMTP_mask.nii');

PT_mask = fullfile(derivativesDir, ...
                    'cpp_spm-roi/sub-pilot001/ses-001/roi/rsub-pilot001_ses-001_acq-r0p75_hemi-R_space-EPIT1w_label-rightPT_mask.nii');

% VISUAL

% to be changed/checked for each subject

% run 1 - hor = 13 | vert = 31
% run 2 - hor = 14 | vert = 32
% run 3 - hor = 15 | vert = 33
% run 4 - hor = 16 | vert = 34
% run 5 - hor = 17 | vert = 35
% run 6 - hor = 18 | vert = 36

fprintf('VISUAL:\n')

visual_hor = 13:18;
visual_vert = 31:36;

for iRun = 1:length(visual_hor)

    files = [fullfile(dir, ['con_00' num2str(visual_hor(iRun), '%02.f') '.nii']); ...
             fullfile(dir, ['con_00' num2str(visual_vert(iRun), '%02.f') '.nii']) ...
             ];
         
    fprintf (" run %i - working with hor beta00%s and vert beta00%s \n", ...
             iRun, num2str(visual_hor(iRun), '%02.f'), num2str(visual_vert(iRun), '%02.f'));

    output_name = ['run-' num2str(iRun, '%03.f') '_desc-wta_label-visualMT_dseg.nii'];

    winner_take_all(files, MT_mask, output_name);     

    output_name = ['run-' num2str(iRun, '%03.f') '_desc-wta_label-visualPT_dseg.nii'];

    winner_take_all(files, PT_mask, output_name);  

end

% keep voxels with same labels

output_name = 'desc-wtaCommonRuns_label-visualMT_dseg.nii';

files = spm_select('FPList', dir, '^run-.*label-visualMT_dseg.nii');
keep_common_labels(files, output_name);

output_name = 'desc-wtaCommonRuns_label-visualPT_dseg.nii';

files = spm_select('FPList', dir, '^run-.*label-visualPT_dseg.nii');
keep_common_labels(files, output_name);

% AUDITORY

% to be changed/checked for each subject

% run 1 - hor = 7  | vert = 25
% run 2 - hor = 8  | vert = 26
% run 3 - hor = 9  | vert = 27
% run 4 - hor = 10 | vert = 28
% run 5 - hor = 11 | vert = 29
% run 6 - hor = 12 | vert = 30

fprintf('AUDITORY:\n')

auditory_hor = 7:12;
auditory_vert = 25:30;

for iRun = 1:length(auditory_hor)

    files = [fullfile(dir, ['con_00' num2str(auditory_hor(iRun), '%02.f') '.nii']); ...
             fullfile(dir, ['con_00' num2str(auditory_vert(iRun), '%02.f') '.nii']) ...
             ];
         
    fprintf (" run %i - working with hor beta00%s and vert beta00%s \n", ... 
             iRun, num2str(auditory_hor(iRun), '%02.f'), num2str(auditory_vert(iRun), '%02.f'));

    output_name = ['run-' num2str(iRun, '%03.f') '_desc-wta_label-auditoryMT_dseg.nii'];

    winner_take_all(files, MT_mask, output_name);     

    output_name = ['run-' num2str(iRun, '%03.f') '_desc-wta_label-auditoryPT_dseg.nii'];

    winner_take_all(files, PT_mask, output_name);  

end 

% keep voxels with same labels

output_name = 'desc-wtaCommonRuns_label-auditoryMT_dseg.nii';

files = spm_select('FPList', dir, '^run-.*label-auditoryMT_dseg.nii');
keep_common_labels(files, output_name);

output_name = 'desc-wtaCommonRuns_label-auditoryPT_dseg.nii';

files = spm_select('FPList', dir, '^run-.*label-auditoryPT_dseg.nii');
keep_common_labels(files, output_name);

% BIMODAL

% to be changed/checked for each subject

% run 1 - hor = 1 | vert = 19
% run 2 - hor = 2 | vert = 20
% run 3 - hor = 3 | vert = 21
% run 4 - hor = 4 | vert = 22
% run 5 - hor = 5 | vert = 23
% run 6 - hor = 6 | vert = 24

fprintf('BIMODAL:\n')

bimodal_hor = 1:6;
bimodal_vert = 19:24;

for iRun = 1:length(bimodal_hor)

    files = [fullfile(dir, ['con_00' num2str(bimodal_hor(iRun), '%02.f') '.nii']); ...
             fullfile(dir, ['con_00' num2str(bimodal_vert(iRun), '%02.f') '.nii']) ...
             ];
         
    fprintf (" run %i - working with hor beta00%s and vert beta00%s \n", ... 
             iRun, num2str(bimodal_hor(iRun), '%02.f'), num2str(bimodal_vert(iRun), '%02.f'));
         
    output_name = ['run-' num2str(iRun, '%03.f') '_desc-wta_label-bimodalMT_dseg.nii'];

    winner_take_all(files, MT_mask, output_name);     

    output_name = ['run-' num2str(iRun, '%03.f') '_desc-wta_label-bimodalPT_dseg.nii'];

    winner_take_all(files, PT_mask, output_name);  

end 

% keep voxels with same labels

output_name = 'desc-wtaCommonRuns_label-bimodalMT_dseg.nii';

files = spm_select('FPList', dir, '^run-.*label-bimodalMT_dseg.nii');
keep_common_labels(files, output_name);

output_name = 'desc-wtaCommonRuns_label-bimodalPT_dseg.nii';

files = spm_select('FPList', dir, '^run-.*bimodalPT_dseg.nii');
keep_common_labels(files, output_name);