function gvol = group_outline(mysbs,u,k,wcon,gt,varargin)
% function gvol = group_outline(mysbs,u,k,wcon,gt,[opt] plot ind outlines)
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
% Tor Wager, 2/15/02
%
% svol is the individual subject mask of results for each subject0
% gvol is the sum of individual subject masks
% gt is the threshold for how many subjects must activate for the voxel to be plotted
%	if gt is a vector, plot will be separated in colors based on multiple levels
%	order: 'g' 'b' 'r'.  maximum of 3 levels

bgimg = spm_get(1,'*img','Choose overlay image.'); % '/usr/private/spm99/templates/scalped_avg152T1.img';
bgimg = reslice_bkg(bgimg,mysbs);

h = [];
ph = [];
fh = [];
if length(gt) > 3, error('Maximum of 3 elements in gt'),end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% Get the ind. subject activations & grp. vol and plot outlines
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% ---------------------------------------------------------------
% * load the canonical overlay image and plot
% ---------------------------------------------------------------

[img,hdr,h] = readim2(bgimg,'p');
colormap(gray)
fh = gcf;

% ---------------------------------------------------------------
% * load the outline colors
% ---------------------------------------------------------------

a = fileparts(which('outline'));
a = fullfile(a,'outline_colors');
load(a)							% variable is called mycolors

mycolors = rand(24,3);

% ---------------------------------------------------------------
% * loop through subjects, get results, make outline, plot
% ---------------------------------------------------------------

index = 1;
for snum = mysbs
	
	SubjCode = snum{1};
  	if index == 1, cd([SubjCode]), else  cd(['../' SubjCode]), end

	% get/plot the outline of the results mask for this subject
	% --------------------------------------------------------------
	[h,ph,fh,svol] = ind_subject_outline(h,ph,fh,u,k,wcon,mycolors(index,:));
	myleg{index} = SubjCode;

	% add subject to group mask
	% --------------------------------------------------------------
	if index == 1, gvol = svol;, else,  gvol = gvol + svol;, end

	index = index + 1;

end

legend(ph,myleg)

close

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% Threshold and plot the group mask
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% ---------------------------------------------------------------
% * threshold the group mask
% ---------------------------------------------------------------
gtvol = gvol;
gtvol(gtvol < gt(1)) = 0;

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
	[h,gh(end+1),fh] = plotmask3(gtvol,2,h,fh);
else
	[h,gh(end+1),fh] = plotmask3(gtvol,2);
end


% ---------------------------------------------------------------
% * plot additional colors, if specified
% ---------------------------------------------------------------

if length(gt) > 1
	colcodes = [2 3 1];
	for i = 2:length(gt)
		gtvol(gtvol < gt(i)) = 0;
		if any(any(any(gtvol)))
			[h,gh(end+1),fh] = plotmask3(gtvol,colcodes(i),h,fh);
		end
	end
end

return



function bgimg = reslice_bkg(bgimg,snum)

	eval(['cd ' snum{1}])
	% take the 1st image here
	D = dir;
    index = 3;
    while D(index).isdir | ~findstr('.img',D(index).name)
	%while D(index).isdir | ~strcmp(D(index).name(end-2:end),'img')
       % D(index).isdir | ~strcmp(D(index).name(end-2:end),'img')
		index = index + 1;
    end
    index = index + 1;
	imgname = D(index).name;
	reslice_imgs(imgname,bgimg,0);
    cd ..
	[d,f,e] = fileparts(bgimg);
	bgimg = fullfile(d,['r' f e]);
	bghdr = fullfile(d,['r' f '.hdr']);
	newbg = fullfile(pwd,['grp_outline_bkg' e]);
	newhdr = fullfile(pwd,['grp_outline_bkg.hdr']);
	eval(['!mv ' bgimg ' ' newbg])  % linux/unix
	eval(['!mv ' bghdr ' ' newhdr])
    eval(['!move ' bgimg ' ' newbg])    % windows
	eval(['!move ' bghdr ' ' newhdr])
	bgimg = newbg;

return
