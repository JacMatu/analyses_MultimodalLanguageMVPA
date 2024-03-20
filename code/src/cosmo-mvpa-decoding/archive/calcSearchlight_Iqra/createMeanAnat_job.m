%-----------------------------------------------------------------------
% Job saved on 08-Oct-2023 22:49:58 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%%
matlabbatch{1}.spm.util.imcalc.input = {
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-001_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-002_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-003_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-004_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-005_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-006_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-007_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-008_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-009_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-010_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-011_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-013_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-014_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-015_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-016_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-017_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-pil001_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-pil002_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-pil004_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        '/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks/rsub-pil005_ses-001_space-IXI549Space_res-r1pt0_desc-brain_mask.nii,1'
                                        };
%%
matlabbatch{1}.spm.util.imcalc.output = 'meanAnat';
matlabbatch{1}.spm.util.imcalc.outdir = {'/Volumes/IqraMacFmri/visTac/fMRI_analysis/outputs/derivatives/SLmasks'};
matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+i16+i17+18+i19+i20)/20';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
