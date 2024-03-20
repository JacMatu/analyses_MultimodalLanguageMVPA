function view_results(files, winner_hdr)
    to_view = cellstr(files);
    for i = 1:numel(winner_hdr)
        to_view{end + 1} = winner_hdr(i).fname;
    end
    spm_check_registration(char(to_view));
end
