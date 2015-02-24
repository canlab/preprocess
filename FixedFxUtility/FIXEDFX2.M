function [vol,maskingImage,V] = fixedfx2(fnames,cl_size,crit_p)
% fixedfx2(fnames,cl_size,crit_p)
% ----------------------------------------------------------
% function for doing fixed effects analysis on a set of images
% 
% finds clusters of size >= cl_size with voxels significant
% at p < crit_p, uncorrected.
%
% fnames: matrix of image names
% 	  should be an spmT image for fixed effects analysis
%
% 
%
% Tor Wager, 10/26/01

	% ======================================================
	% * File stuff
	% ======================================================

	% get file names
	% ------------------------------------------------------
	check_fnames(fnames)
	[p,f,e] = fileparts(fnames(1,:));


	% get output directory - omit the last (subject) directory
	% ------------------------------------------------------
	outdir = p(1:max(find(p == filesep)));


	% create mean img
	% ------------------------------------------------------
	meanName = [outdir 'mean_' f '.img'];
	tor_spm_mean_ui(fnames,meanName) 



	% ======================================================
	% * Stats stuff - load and threshold image
	% ======================================================


	nimages = size(fnames,1);
	crit_t = tinv(1-crit_p,nimages-1);

	% Read mean image and multiply by sqrt(n) to get fixed effect t map
	% ------------------------------------------------------
	[vol,hdr] = readim2(meanName,'t');
	vol = vol .* hdr.SPM_scale;			% account for the SPM scaling factor!
	% 90 degree rotation to get in canonical format
	%for i = 1:size(vol,3)
	%	vol2(:,:,i) = rot90(vol(:,:,i));
	%end
	%vol = vol2;

	vol = vol .* sqrt(nimages);

	disp(['Max fixed effects t = ' num2str(max(max(max(vol))))])
	

	% Make masking image
	% ------------------------------------------------------
	disp(['Making masking image at p < ' num2str(crit_p) ', t > ' num2str(crit_t) ', cluster size >= ' num2str(cl_size)])
	[maskedImage maskingImage] = maskImg(vol,crit_t,Inf);
	[maskingImage,numClusters,XYZ] = clusterSizeMask(cl_size,maskingImage);

	if sum(sum(sum(maskingImage))) == 0
		disp('No significant results.  No image written - exiting.')
		vol = []; V = [];
		return
	else
		disp(['Clusters found: ' num2str(numClusters)])
	end

	% Threshold image with mask
	% ------------------------------------------------------
	disp('Thresholding image')
	vol = vol .* maskingImage;

	% Get Z output - intensity values for sig. voxels
	% ------------------------------------------------------
	for i = 1:size(XYZ,2)
		% row is y, col is x
		Z(i) = vol(XYZ(1,i),XYZ(2,i),XYZ(3,i));
	end


	% ======================================================
	% * Write out the image and give summary
	% ======================================================
	% get the correct name and directory
	% ------------------------------------------------------
	VN = spm_vol(meanName);
	[p,f,e] = fileparts(VN.fname);
	f = ['fixedfx_' f];

	% use the V structure from the 1st input image
	% ------------------------------------------------------
	V = spm_vol(fnames(1,:));
	V.fname = [p filesep f e];
	V.descrip = ['fixedfx p = ' num2str(crit_p) ' cl_size = ' num2str(cl_size) ' ' V.descrip];


	% write image
	% ------------------------------------------------------
	disp(['Writing image:' V.fname])
	V = spm_write_vol(V,vol);


	% additional things to save with V in output mat file
	% ------------------------------------------------------
	V.maskingImage = maskingImage;
	V.crit_p = crit_p;
	V.cl_size = cl_size;
	V.crit_t = crit_t;
	V.numClusters = numClusters;
	V.XYZ = XYZ;
	V.Z = Z;
	V.XYZmm = voxel2mm(XYZ,V.mat);

	% for compatibility with SPM struct and cluster analysis
	% ------------------------------------------------------
	V.voxSize = [hdr.xsize hdr.ysize hdr.zsize]';
	V.u = crit_t;
	V.k = cl_size;
	V.title = f;

	% save V in mat file
	% ------------------------------------------------------
	[p,f,e] = fileparts(V.fname);
	f = ['V_' f];
	fname = [p filesep f];
	disp(['Saving V information in: ' fname '.mat'])
	eval(['save ' fname ' V'])

	
return







function check_fnames(fnames)
	warnMe = 0;

	for i = 1:size(fnames,1)
		[p{i},f{i},e] = fileparts(fnames(i,:));
	end

	if any(~strcmp(f{1},f(:))),
		warning('File names are not all the same contrast!')
		fnames
	end
return

	