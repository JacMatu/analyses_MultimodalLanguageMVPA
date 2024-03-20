function common_hdr = keep_common_labels(images_list, file_name)

    hdr = spm_vol(images_list);
    vols = spm_read_vols(hdr);

    % if all indices are the same for all volumes
    % then the sum across 4th dimension
    % should be equal to nb_vol X value of of 1rst volume
    sum_across_vols = sum(vols, 4, 'omitnan');
    first_vol = vols(:, :, :, 1);
    n_times_first_vol = first_vol * size(vols, 4);

    common = sum_across_vols == n_times_first_vol;

    common_labels = zeros(size(common));
    common_labels(common) = first_vol(common);

    common_hdr = hdr(1);
    common_hdr.descrip = 'common labels';
    common_hdr.fname = spm_file(common_hdr.fname, 'filename', file_name);

    spm_write_vol(common_hdr, common_labels);

end
