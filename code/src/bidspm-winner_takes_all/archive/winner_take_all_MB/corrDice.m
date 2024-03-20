% motion_winner_take_all

subLabel = 'pilot005';

derivativesDir = '/Users/barilari/data/V5_high-res_pilot_bimodalMotion_boldClassic/outputs/derivatives';

statsDirs = fullfile(derivativesDir, 'cpp_spm-stats');

dir = fullfile(statsDirs, ...
                [ 'sub-' subLabel '/task-gratingBimodalMotion_space-individual_FWHM-0_desc-gratingBimodalMotionBOLDSplitConditionsMVPALatest' ]);

MT_mask = fullfile(derivativesDir, ...
                    [ 'cpp_spm-roi/sub-' subLabel '/ses-001/roi/rsub-' subLabel '_ses-001_acq-r0p75_hemi-R_space-EPIT1w_label-rightMTP_mask.nii' ]);

PT_mask = fullfile(derivativesDir, ...
                    [ 'cpp_spm-roi/sub-' subLabel '/ses-001/roi/rsub-' subLabel '_ses-001_acq-r0p75_hemi-R_space-EPIT1w_label-rightPT_mask.nii' ]);

                
% MT
visualMapMT = fullfile(dir, 'desc-wtaMean_label-visualMT_dseg.nii'); 

auditoryMapMT = fullfile(dir, 'desc-wtaMean_label-auditoryMT_dseg.nii'); 

visualMT = spm_summarise(visualMapMT, MT_mask);

auditoryMT = spm_summarise(auditoryMapMT, MT_mask);

sum(visualMT(:) == 1);
sum(auditoryMT(:) == 1);help dice

dice(visualMT, auditoryMT)
corr2(visualMT, auditoryMT)

% sub-001 
%   dice 0.5332 0.4550
%   corr -0.019

% sub-002
%   dice 0.5314 0.3766
%   corr -0.0901

% sub-003
%   dice 0.6062 0.5877
%   corr 0.1941

% PT
visualMapPT = fullfile(dir, 'desc-wtaMean_label-visualPT_dseg.nii'); 

auditoryMapPT = fullfile(dir, 'desc-wtaMean_label-auditoryPT_dseg.nii'); 

visualPT = spm_summarise(visualMapPT, PT_mask);

auditoryPT = spm_summarise(auditoryMapPT, PT_mask);

sum(visualPT(:) == 1)
sum(auditoryPT(:) == 1)

dice(visualPT, auditoryPT)

corr2(visualPT, auditoryPT)

% sub-001 
%   dice 0.571 0.5279
%   corr 0.1

% sub-002
%   dice 0.5615 0.4048
%   corr -0.0259

% sub-003
%   dice 0.6291 0.5860
%   corr 0.2158