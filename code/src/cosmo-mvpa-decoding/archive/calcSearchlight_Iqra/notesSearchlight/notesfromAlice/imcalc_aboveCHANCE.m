
clear all
spm fmri




%% Demo: fMRI searchlights LDA classifier
%
%To run MVPA with LDA classifier, crossvalidation method
%
% #   For CoSMoMVPA's copyright information and license terms,   #
% #   see the COPYING file distributed with CoSMoMVPA.           #


clear all;

%PILOT 
sub(1).id      =    'sub-01';
sub(1).Nrun    =    '0';
sub(1).folder  =    'spmT_MVPA';
%sub(1).notes  =    ONLY LOC

sub(2).id      =    'sub-02' ;
sub(2).Nrun    =    '0';
sub(2).folder  =    'spmT_MVPA';
%sub(2).notes  =    ONLY LOC

sub(3).id      =    'sub-03';
sub(3).Nrun    =    '0';
sub(3).folder  =    'spmT_MVPA';
%sub(3).notes  =    PBM SES02 - removed from analysis

sub(4).id      =    'sub-04' ;
sub(4).Nrun    =    '19';
sub(4).folder  =    'spmT_MVPA';
%sub(4).notes  =     NA

sub(5).id      =    'sub-05' ;
sub(5).Nrun    =    '20';
sub(5).folder  =    'spmT_MVPA';
%sub(5).notes  =    NA

sub(6).id      =    'sub-06' ;
sub(6).Nrun    =    '20';
sub(6).folder  =    'spmT_MVPA';
%sub(6).notes  =    NA

sub(7).id      =    'sub-07' ;
sub(7).Nrun    =    '20';
sub(7).folder  =    'spmT_MVPA';
%sub(7).notes  =    some movements especially in MVPA runs of ses-01
%(A3, V3, A5, A8)

sub(8).id      =    'sub-08' ;
sub(8).Nrun    =    '20';
sub(8).folder  =    'spmT_MVPA';
%sub(8).notes  =    was tired/fell asleep in ses-02 from run-05 to run-10 ?

sub(9).id      =    'sub-09' ;
sub(9).Nrun    =    '20';
sub(9).folder  =    'spmT_MVPA';
%sub(9).notes  =    phonoloc : did not hear enough ?? keypress weird ! 

sub(10).id      =    'sub-10' ;
sub(10).Nrun    =    '20';
sub(10).folder  =    'spmT_MVPA';
%sub(10).notes  =    mvmt in ses-01_A1 and V8 and ses-02_V6

sub(11).id      =    'sub-11' ;
sub(11).Nrun    =    '20';
sub(11).folder  =    'spmT_MVPA';
%sub(11).notes  =    LEFT HANDED !!! double check localizers! 
%sub(11).notes  =    fell asleep in ses-02_A3 ? 

sub(12).id      =    'sub-12' ;
sub(12).Nrun    =    '20';
sub(12).folder  =    'spmT_MVPA';
%sub(12).notes  =    mvmt in ses-02_V10

sub(13).id      =    'sub-13' ;
sub(13).Nrun    =    '20';
sub(13).folder  =    'spmT_MVPA';
%sub(13).notes  =    NA



no_sub=[4:13];
task = 'Aud'; %Vis or Aud
feature = 'Vowels'; %'Cons', 'Vowels' or 'Speak'
current_directory=cd;


for isub=no_sub
%imcalc
matlabbatch{1}.spm.util.imcalc.input = {strcat(current_directory,'/',feature,'/',sub(isub).id,'_searchlight_100vx_svm_MVPA',task,'_',feature,'_beta.nii,1')};
matlabbatch{1}.spm.util.imcalc.output = strcat('AbvChance_', sub(isub).id,'_searchlight_100vx_svm_MVPA',task,'_',feature,'_beta.nii');
matlabbatch{1}.spm.util.imcalc.outdir = {fullfile(current_directory,feature)};
matlabbatch{1}.spm.util.imcalc.expression = 'i1-0.333'; %%i1-chance level
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

%smoothing DOES NOT WORK - output not saved
matlabbatch{2}.spm.spatial.smooth.data = {strcat(current_directory,'/',feature,'/', 'AbvChance_', sub(isub).id,'_searchlight_100vx_svm_MVPA',task,'_',feature,'_beta.nii')};
matlabbatch{2}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's8';
matlabbatch{2}.spm.spatial.smooth.outdir = {strcat(current_directory,'/res_aboveChance_smoothed')};

spm_jobman('run',matlabbatch);
end
