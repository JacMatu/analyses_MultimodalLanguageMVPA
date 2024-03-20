function winner_hdr = winner_take_all(beta_images, mask_image, file_name)
    %
    % USAGE::
    %
    %     winner_hdr = winner_take_all(beta_images)
    %     winner_hdr = winner_take_all(beta_images, mask_image)
    %
    % :param beta_images: list of all images to compare
    % :type  beta_images: cellstr
    %
    % :param mask_image: Path to a mask image.
    %                    Must be of same resolution as the beta_images.
    % :type  mask_image: char
    %

    if nargin < 2
        mask_image = '';
    end

    betas_hdr = spm_vol(beta_images);
    betas = spm_read_vols(betas_hdr);

    [~, max_indices] = max(betas, [], 4);

    winner_hdr = betas_hdr(1);
    winner_hdr.descrip = 'winner take all labels';
    winner_hdr.fname = spm_file(winner_hdr.fname, 'filename', file_name);

    winner = nan(winner_hdr.dim);

    if ~isempty(mask_image)
        mask = spm_read_vols(spm_vol(mask_image));
        winner(mask == 1) = max_indices(mask == 1);

    else
        winner = max_indices;

    end

    spm_write_vol(winner_hdr, winner);

end
