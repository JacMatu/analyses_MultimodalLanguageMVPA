function mask = create_mask(dim, location, vox_size)

    if nargin < 1
        dim = 40;
    end

    if nargin < 2
        location = 1;
    end

    if nargin < 3
        vox_size = 1;
    end

    hdr_mask = struct('fname',   'mask.nii', ...
                      'dim',     [dim, dim, dim], ...
                      'dt',      [spm_type('uint8') spm_platform('bigend')], ...
                      'mat',     transformation_matrix(vox_size), ...
                      'pinfo',   [1 0 0]', ...
                      'descrip', 'resultant analysis mask');

    % generate mask
    mask = false(hdr_mask.dim);

    if numel(location) == 1

        location = [location, dim - location; ...
                    location, dim - location; ...
                    location, dim - location];

    elseif all(size(location) == [3, 2])

    else
        error(['location must either be:\n', ...
               '- a scalar that defines the margin width\n', ...
               '- 3 X 2 defining the mask position in x, y, z.']);
    end

    mask(location(1, 1):location(1, 2), ...
         location(2, 1):location(2, 2), ...
         location(3, 1):location(3, 2)) = true;

    spm_write_vol(hdr_mask, mask);
end
