function [h,ph,fh,vol,TOR.multCompCorrect] = ind_subject_outline(h,ph,fh,u,k,wcon,mycol)
% function [h,ph,fh,vol,TOR.multCompCorrect] = ind_subject_outline(h,ph,fh,u,k,wcon,mycol)
%
% u		height threshold
%		corr. method is hard-coded in ind_subject_outline
% k		extent threshold
% wcon		index number of which contrast to show in results
%
% vol		individual subject sig mask, un-outlined
%
% Tor Wager 2/15/02

% -------------------------------------------------------------
% * process for individual subject
% -------------------------------------------------------------

% get the results for this subject
% - - - - - - - - - - - - - - - - - - - - - - - - - -
	TOR.u = u;
	% height threshold - T or p value

	TOR.k = k;
	% extent threshold - number of voxels

	TOR.Ic = wcon;
	% optional: which contrast to show in results (number)

	TOR.resdir = pwd;
	TOR.maskWOtherContrasts = 0;

	TOR.multCompCorrect = 'FWE'; %'FDR'; 
	% correct for multiple comparisons
	% choices are 'FWE|FDR|uncorrected'

	[hReg,SPM,VOL,xX,xCon,xSDM] = tor_spm_results_ui(TOR);


% write results mask
% - - - - - - - - - - - - - - - - - - - - - - - - - - 
	vol = voxel2mask(SPM.XYZ',VOL.DIM');
	V = spm_vol('mask.img');   
	V.fname = 'res_mask.img';
	V = spm_write_vol(V,vol);

	% load thresholded results mask (with cluster size)
	% [vol,hdr] = readim2('res_mask');


% outline the mask
% - - - - - - - - - - - - - - - - - - - - - - - - - -
ovol = outline(vol);


% change to mm coordinates and back
% so that we can map to standard volume
% not used right now, because we have resliced the anatomical
% to be in the same space as the functional (coregistered)
% - - - - - - - - - - - - - - - - - - - - - - - - - -
%vox = mask2voxel(vol);
%vox = vox';
%vox = voxel2mm(vox,V.mat);
%V2 = spm_vol('/usr/private/spm99/canonical/avg_152_T1_biman5.img');
%V2.M = V2.mat;
%vox = mm2voxel(vox,V2);
%vol = voxel2mask(vox,V2.dim(1:3));


% plot the mask in the specified input color
% figure is created automatically if necessary
% - - - - - - - - - - - - - - - - - - - - - - - - - -
% mycol = rand(1,3);

% if ~exist('h'), h = [];, end
% if ~exist('ph'), ph = [];, end

if ~isempty(fh) 
	[h,ph(end+1),fh] = plotmask3(ovol,mycol,h,fh);
else
	[h,ph(end+1),fh] = plotmask3(ovol,mycol);
end

return