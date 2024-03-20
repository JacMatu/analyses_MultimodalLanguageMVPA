%% Demo: fMRI searchlights LDA classifier
%
%To run MVPA with LDA classifier, crossvalidation method
%
% #   For CoSMoMVPA's copyright information and license terms,   #
% #   see the COPYING file distributed with CoSMoMVPA.           #


clear all;


%subject definition (see function in same dir)
global sub
sub=sub_data;
no_sub= 4:9; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data and brain mask to run searchlight within 
%config=cosmo_config();

%% Define dir&paths
mvpa_dir = fileparts(mfilename('fullpath'));
mask_dir = fullfile(mvpa_dir, 'ROIs','Masks');
output_dir = fullfile(mvpa_dir, 'Searchlight');

lipspeech_dir = fullfile(mvpa_dir, '..'); %root dir
bids_dir = fullfile(lipspeech_dir, 'Lip_BIDS');
stats_dir = fullfile(bids_dir, 'derivatives', 'bidspm-stats');

%% Define subjects, ROIs, tasks to work on

sub_all=sub_data; %call function in same dir named sub_data.m - with all info on each subject
nsub_all = length(sub_all);

algo = 'svm'; %'lda' or 'svm'
val = 'beta'; %or 'beta' ?

%%Define paths 
%%%%%study_path= cd;%one folder back 

%tmap to use as data 
task_list = {'Aud', 'Vis'}; %Vis or Aud
model_list = {'Vowels', 'Cons', 'Speak'}; %'Cons', 'Vowels' or 'Speak'

for t = 1:length(task_list)
    task = task_list{t};
    
    for f = 1:length(model_list)
        model = model_list{f};
       
        data_img = strcat('_task-MVPA',task, '_space-IXI549Space_desc-4D_tmap.nii');

            % reset citation list
            %cosmo_check_external('-tic');


        for s = no_sub

            %exctract the subname to select the correct data 4D file 
            sub_name=(sub(s).id);
            
            working_on=strcat(sub_name,'_searchlight_100vx_', algo, '_MVPA', task, '_', model, '_', val);
            disp(working_on)

            mask_fn=fullfile(mask_dir,'fullbrainMASK_mni.nii');% whole brain mask excluding cerebellum

            %find the data of this subject
            data_img=fullfile(stats_dir, sub_name, ...
                strcat('/task-MVPA', task, '_space-IXI549Space_FWHM-2_node-MVPA', task, model), ...
                strcat(sub_name, '_task-MVPA', task, '_space-IXI549Space_desc-4D_', val, '.nii'));
            
            %% LDA classifier searchlight analysis %%
            % This analysis identified brain regions where the categories can be
            % distinguished using a a Linear Discriminant Analysis (LDA) classifier.

            %set targets
            targets=repmat(1:3,1,str2num(sub(s).Nrun))'; %there are 3 consonants (mean value for the 9 iteration of the cons), and each is repeated in each run. 
            targets= sort(targets); %comment this line if want to see a "random" decoding - labels will not correspond anymore

            %set chunks, corresponds to the number of runs
            chunks = repmat(1:str2num(sub(s).Nrun), 1, 3)';

            ds_per_run = cosmo_fmri_dataset(data_img, 'mask', mask_fn,...
                                            'targets',targets,'chunks',chunks);

            %remove constant features (due to liberal masking)
            ds_per_run=cosmo_remove_useless_data(ds_per_run);

            % print dataset
            fprintf('Dataset input:\n');
            cosmo_disp(ds_per_run);


            % Use the cosmo_cross_validation_measure and set its parameters
            % (classifier and partitions) in a measure_args struct.
            measure = @cosmo_crossvalidation_measure;
            measure_args = struct();

            % Define which classifier to use, using a function handle.
            % Alternatives are @cosmo_classify_{svm,matlabsvm,libsvm,nn,naive_bayes}
            if strcmp(algo, 'lda')
                measure_args.classifier = @cosmo_classify_lda;
            elseif strcmp(algo, 'svm')
                measure_args.classifier = @cosmo_classify_svm;
            end 

            %%  Set partition scheme. odd_even is fast; for publication-quality analysis nfold_partitioner is recommended.
            %I use nfold partitions here

            % Other alternatives are:
            % - cosmo_nfold_partitioner    (take-one-chunk-out crossvalidation)
            % - cosmo_nchoosek_partitioner (take-K-chunks-out  "             ").
            % measure_args.partitions = cosmo_oddeven_partitioner(ds_per_run);
            measure_args.partitions = cosmo_nfold_partitioner(ds_per_run);
            measure_args.normalization = 'zscore';
            % print measure and arguments
            fprintf('Searchlight measure:\n');
            cosmo_disp(measure);
            fprintf('Searchlight measure arguments:\n');
            cosmo_disp(measure_args);

            %% Define the size of the voxels' sphere you want to take for each searchlight

            % Here I define a neighborhood with approximately 100 voxels in each searchlight.
            nvoxels_per_searchlight=100;
            nbrhood=cosmo_spherical_neighborhood(ds_per_run,...
                                    'count',nvoxels_per_searchlight);


            % Run the searchlight
            sl_results = cosmo_searchlight(ds_per_run,nbrhood,measure,measure_args);

            % print output dataset
            fprintf('Dataset output:\n');
            cosmo_disp(sl_results);

            % Plot the output
            cosmo_plot_slices(sl_results);

            % Define the name of the file to be saved in the output location
            output_fn = fullfile(output_dir, strcat(working_on, '.nii'));

            % Store results to disc
            cosmo_map2fmri(sl_results, output_fn);

            % Show citation information
            cosmo_check_external('-cite');
        end
    end 
end 