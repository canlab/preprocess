function [svol,vol1] = create_smoothed_pvals(x,y,z,voxmm,S)
% function [svol, vol1] = create_smoothed_pvals(x,y,z,voxmm,S)
% e.g. vol1 = create_smoothed_pvals(64,64,16,4,10);
% Tor Wager

vol1 = rand(x,y,z);

% smoothing
% create image of 4 mm voxels
V.fname = 'temp_vol.img';
V.dim = [x y z 4];
V.mat = [voxmm 0 0 -32*4; 0 voxmm 0 -32*4; 0 0 voxmm -32; 0 0 0 1];
V.pinfo = [1;0;0];
V.descrip = 'volume with 4 mm voxels';
spm_write_vol(V,vol1);

P = 'temp_vol.img';
Q = 'stemp_vol.img';
spm_smooth(P,Q,S);

[svol,hdr] = readim2(Q(1:end-4),'t');

svol = svol .* hdr.SPM_scale;

% imagesc(vol2(:,:,10))

return