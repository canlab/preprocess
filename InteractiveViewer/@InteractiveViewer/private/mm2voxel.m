function vox = mm2voxel(mm, M)
    vox = M\[mm(:)' 1]';
    vox = vox(1:3);
end