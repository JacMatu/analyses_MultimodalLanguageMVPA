function M = transformation_matrix(vox_size)

    if nargin < 1
        vox_size = 1;
    end

    % voxel to world transformation matrix
    M  = [-vox_size   0           0           84; ...
          0           vox_size    0           -116; ...
          0           0           vox_size    -56; ...
          0           0           0           1];
end
