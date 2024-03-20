% test winner_take_all

addpath(fullfile(pwd, '..'));

nb_conditions = 4;
nb_runs = 2;

generate_data(nb_conditions, nb_runs);

files = spm_select('FPList', pwd, '^sbeta.*.nii');

%% basic
winner_hdr = winner_take_all(files);
view_results(files, winner_hdr);

%% with mask
winner_hdr = winner_take_all(files, fullfile(pwd, 'mask.nii'));
view_results(files, winner_hdr);

%% with ROI mask
ROI_bouding_box = [10, 20; ...
                   20, 30; ...
                   15, 25];
create_mask(dim, ROI_bouding_box);
winner_hdr = winner_take_all(files, fullfile(pwd, 'mask.nii'));
view_results(files, winner_hdr);

%% winner take all for each run
% beta 1 to 4 come from run 1
files = spm_select('FPList', pwd, '^sbeta_000[1-4].nii');
winner_hdr = winner_take_all(files);
winner(1).fname = spm_file(winner_hdr.fname, 'prefix', 'run-1_');
movefile(winner_hdr.fname, winner(1).fname);

% beta 5 to 8 come from run 1
files = spm_select('FPList', pwd, '^sbeta_000[5-8].nii');
winner_hdr = winner_take_all(files);
winner(2).fname = spm_file(winner_hdr.fname, 'prefix', 'run-2_');
movefile(winner_hdr.fname, winner(2).fname);

files = spm_select('FPList', pwd, '^sbeta.*.nii');
view_results(files, winner);

%% keep voxels with same labels
files = spm_select('FPList', pwd, '^run-.*_dseg.nii');
common_hdr = keep_common_labels(files);

spm_check_registration(spm_select('FPList', pwd, '.*_dseg.nii'));
