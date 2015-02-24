function xY = tor_get_roi(connum,imnames,threshold)
% function xY = tor_get_roi(connum,imnames,threshold)
%
% this function gets timeseries data from an ROI
% input:
% connum: an SPM contrast number in string format - '0002'
% imnames: a matrix of image names, in spm_list_files output format
% threshold: a t-threshold at which to image activation blobs
%
% xY.Y is the average of voxels in ROI
% 
% run this in the SPM results directory.
% 10/17/01 by Tor Wager

	[a] = readim2(['spmT_' connum);
	a = double(a);
	[mask,maskingimg] = maskImg(a,threshold,Inf);


	[dummy,dummy,h] = readim2(maskingimg,'p');

	voxels = getVoxels2(h,maskingimg);

	close

	O.mask = maskingimg;
	O.coords = voxels;

	% slow
	% ts = timeseries2('volume',imnames,O);

	ts = timeseries2('multi',imnames,O);

	xY.Y = ts.avg;
	xY.y = ts.indiv;
	xY.xyz = voxels;
	xY.dstr = 'not adjusted';

return

