% motion_winner_take_all

subLabel = 'pilot004';

derivativesDir = '/Users/barilari/data/V5_high-res_pilot_bimodalMotion_boldClassic/outputs/derivatives';

statsDirs = fullfile(derivativesDir, 'cpp_spm-stats');

dir = fullfile(statsDirs, ...
                [ 'sub-' subLabel '/task-gratingBimodalMotion_space-individual_FWHM-0_desc-gratingBimodalMotionBOLDSplitConditionsMVPALatest' ]);

MT_mask = fullfile(derivativesDir, ...
                    [ 'cpp_spm-roi/sub-' subLabel '/ses-001/roi/rsub-' subLabel '_ses-001_acq-r0p75_hemi-R_space-EPIT1w_label-rightMTP_mask.nii' ]);

PT_mask = fullfile(derivativesDir, ...
                    [ 'cpp_spm-roi/sub-' subLabel '/ses-001/roi/rsub-' subLabel '_ses-001_acq-r0p75_hemi-R_space-EPIT1w_label-rightPT_mask.nii' ]);

% VISUAL

% 39 : visual hor
% 42 : visual vert

files = [fullfile(dir, 'con_0045.nii'); ...
         fullfile(dir, 'con_0048.nii')];
     
output_name = 'desc-wtaMean_label-visualMT_dseg.nii';

winner_take_all(files, MT_mask, output_name);     

output_name = 'desc-wtaMean_label-visualPT_dseg.nii';

winner_take_all(files, PT_mask, output_name);  
     
% AUDITORY

% 38 : auditory hor
% 41 : auditory vert

files = [fullfile(dir, 'con_0044.nii'); ...
         fullfile(dir, 'con_0047.nii')];

output_name = 'desc-wtaMean_label-auditoryMT_dseg.nii';

winner_take_all(files, MT_mask, output_name);     

output_name = 'desc-wtaMean_label-auditoryPT_dseg.nii';

winner_take_all(files, PT_mask, output_name);  

% BIMODAL

% 37 : bimodal hor
% 40 : bimodal vert

files = [fullfile(dir, 'con_0043.nii'); ...
         fullfile(dir, 'con_0046.nii')];

output_name = 'desc-wtaMean_label-bimodalMT_dseg.nii';

winner_take_all(files, MT_mask, output_name);     

output_name = 'desc-wtaMean_label-bimodalPT_dseg.nii';

winner_take_all(files, PT_mask, output_name);  