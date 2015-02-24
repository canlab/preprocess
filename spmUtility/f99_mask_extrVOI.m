%function spm_mask(P1,P2, thresh)
% Masks Images.
% FORMAT spm_mask(P1,P2, thresh)
% P1     - matrix of input image filenames from which
%          to compute the mask.
% P2     - matrix of input image filenames on which
%          to apply the mask.
% thresh - optional threshold(s) for defining the mask.
% The masked images are prepended with the prefix `m'.
%
% If any voxel in the series of images is zero (for data types without
% a floating point representation) or does not have a finite value (for
% floating point and double precision images), then that voxel is set to
% NaN or zero in all the images.  If a threshold, or vector of
% thresholds is passed, then the masking is mased on voxels whos
% values are above all the thresholds.
%
% Images sampled in different orientations and positions can be passed
% to the routine.  Providing the `.mat' files are correct, then these
% should be handled appropriately.
%_______________________________________________________________________
% @(#)spm_mask.m	2.7 John Ashburner 99/04/01


% extraction of coordinates of VOI
% ------------------------------------------------------------------------------
[TabDat, selXYZ] = spm_VOI(SPM,VOL,[],[],hReg);

P1=spm_get(Inf,'.img','Select your new mask image');

V1=spm_vol(P1);
m1=prod(size(V1));

% Create headers
% ------------------------------------------------------------------------------
VO=V1;
M   = VO.mat;
dim = VO.dim(1:3);


spm_progress_bar('Init',VO(1).dim(3),'Masking','planes completed')
for j=1:dim(3),

	msk = zeros(dim(1:2));
        ff = find(selXYZ(3,:) == j);
        for i = 1:length(ff)
            msk(selXYZ(1,ff(i)),selXYZ(2,ff(i))) = 1;
        end
	Mi  = spm_matrix([0 0 j]);

	% Write the images.
	M1  = M\V1.mat\Mi;
        img = spm_slice_vol(V1,M1,dim(1:2),[1 0]);
        img = msk;
        spm_write_plane(V1,img,j);

	spm_progress_bar('Set',j);
end;
spm_progress_bar('Clear');
return;
