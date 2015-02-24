function overlay_mask
% function overlay_mask
%
% You must start in the directory above the individual subject results directories!
% Background overlay image must be in same space as individual subjects!
% Background image is hard-coded in this script!
%
% mysbs		string of individual subject directory codes
% u		height threshold
%		corr. method is hard-coded in ind_subject_outline
% k		extent threshold
% wcon		index number of which contrast to show in results
%
% gvol		group mask, with value in elements indicating number of subjects 
%		activating that voxel.
%
% Tor Wager, 5/03/02
%
% svol is the individual subject mask of results for each subject0
% gvol is the sum of individual subject masks
% gt is the threshold for how many subjects must activate for the voxel to be plotted
%	if gt is a vector, plot will be separated in colors based on multiple levels
%	order: 'g' 'b' 'r'.  maximum of 3 levels

%bgimg = '/usr/private/spm99/templates/scalped_avg152T1.img';

bgimg = spm_get(1,'*.img','Choose background img');
gimg = spm_get(1,'*.img','Choose results mask');

bgimg = reslice_bkg(bgimg,gimg);
gvol = readim2(gimg,'t');

% ---------------------------------------------------------------
% * read and plot the background image
% ---------------------------------------------------------------
[img,hdr,h] = readim2(bgimg,'p');
colormap(gray)
fh = gcf;
gh = [];

% ---------------------------------------------------------------
% * build color map (now done in plotmask3)
% ---------------------------------------------------------------
%cols = zeros(length(mysbs),3);
%cols(:,3) = (1:length(mysbs))' ./ length(mysbs);



% ---------------------------------------------------------------
% * plot the group mask
% ---------------------------------------------------------------
if ~isempty(fh)
	[h,gh(end+1),fh] = plotmask3(gvol,2,h,fh);
else
	[h,gh(end+1),fh] = plotmask3(gvol,2);
end



return


function bgimg = reslice_bkg(bgimg,imgname)

	
	reslice_imgs(imgname,bgimg,0);
	[d,f,e] = fileparts(bgimg);
	bgimg = fullfile(d,['r' f e]);
	bghdr = fullfile(d,['r' f '.hdr']);
	newbg = fullfile(d,['grp_outline_bkg' e]);
	newhdr = fullfile(d,['grp_outline_bkg.hdr']);
	eval(['!mv ' bgimg ' ' newbg])
	eval(['!mv ' bghdr ' ' newhdr])
	bgimg = newbg;
	cd ..
return

