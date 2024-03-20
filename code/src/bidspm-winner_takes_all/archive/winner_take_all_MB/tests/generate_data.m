function generate_data(nb_conditions, nb_runs)
    %
    % generates a set of 3D beta maps for ``nb_conditions`` conditions
    % with one activation in a given location for each
    %
    % across runs only the background noise changes:
    % activations location, size and amplitude stays the same
    %
    % to resemble a typical beta map of SPM, this also:
    %
    % - adds a margin of NaN values on the edge of the beta images
    % - smooths the beta maps
    %

    if nargin < 1
        nb_conditions = 4;
    end

    if nargin < 2
        nb_runs = 2;
    end

    % start from scratch
    delete *.nii;

    % nb of voxels in x, y and z dimensions
    dim = 40;

    % in mm
    vox_size = 1;

    % nb of voxels used as margin of the side of the image
    margin = 5;

    %%
    mask = create_mask(dim, margin, vox_size);

    %% generate data

    % define activations for each condition
    %
    % thay are all cube activations with
    % - size
    % - location: defines the (x1, y1, z1) corner of the activation
    % - amplitude: randomly distributed
    %
    sz = randi([3 7], nb_conditions, 1);
    amp = randn(nb_conditions, 1) * 5;
    for cdt = 1:nb_conditions
        loc(cdt, 1:3) = randi([margin + sz(cdt), dim - margin - sz(cdt)], 1, 3);
    end

    % prepare image headers
    hdr_betas(1:(nb_conditions * nb_runs)) = deal(hdr_beta(dim, vox_size));

    beta_idx = 1;

    for run = 1:nb_runs

        for cdt = 1:nb_conditions

            hdr_betas(beta_idx).fname   = sprintf('beta_%04d.nii', beta_idx);
            hdr_betas(beta_idx).descrip = sprintf('beta (%04d) - condition %i - run %i', ...
                                                  beta_idx, ....
                                                  cdt, ...
                                                  run);

            vol = nan(hdr_betas(cdt).dim);

            % add normally distributed noise
            vol(mask) = randn(repmat(numel(margin:dim - margin), 1, 3));

            vol(loc(cdt, 1):loc(cdt, 1) + sz(cdt) - 1, ...
                loc(cdt, 2):loc(cdt, 2) + sz(cdt) - 1, ...
                loc(cdt, 3):loc(cdt, 3) + sz(cdt) - 1) = amp(cdt) + ...
                                                           randn(sz(cdt), sz(cdt), sz(cdt));

            spm_write_vol(hdr_betas(beta_idx), vol);

            beta_idx = beta_idx + 1;

        end

    end

    %% smooth betas

    files = spm_select('FPList', pwd, '^beta.*.nii');

    matlabbatch{1}.spm.spatial.smooth.data = cellstr(files);
    matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{1}.spm.spatial.smooth.im = 1;

    spm_jobman('run', matlabbatch);

end

function hdr = hdr_beta(dim, vox_size)

    hdr = struct('fname',   [], ...
                 'dim',     [dim, dim, dim], ...
                 'dt',      [spm_type('float32') spm_platform('bigend')], ...
                 'mat',     transformation_matrix(vox_size), ...
                 'pinfo',   [1 0 0]', ...
                 'descrip', 'spm_spm:beta');

end
