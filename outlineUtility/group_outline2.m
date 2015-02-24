function [gvol,clusters] = group_outline2(mysbs,u,k,wcon,gt,varargin)
% function [gvol,clusters] = group_outline2(mysbs,u,k,wcon,gt,[opt] bkground img name, cell array of colors)
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
% gt is the threshold for how many subjects must activate for the voxel to be plotted
%	if gt is a vector, plot will be separated in colors based on multiple levels
%	order: 'g' 'b' 'r'.  maximum of 3 levels recommended but not necessary.
%	elements of gt should be in ASCENDING ORDER!
%
% Optional inputs:
% background image name; try using spm_get to get a filename of an analyze img file
% cell array of colors: e.g., {'g' 'b' 'r' 'y'}
%
% gvol		group mask, with value in elements indicating number of subjects 
%		activating that voxel.
%
% Tor Wager, 2/15/02    Improved version over group_outline, but may be quite a bit slower.
% This version does NOT plot individual outlines, as the previous group_outline does.
%
% svol is the individual subject mask of results for each subject0
% gvol is the sum of individual subject masks
%
%
% Uses montage_clusters.m
% Example:
% Po = spm_get(1,'*img','Choose anatomical overlay image');
% [gvol,clusters] = group_outline2(EXPT.subjects,.05,5,5,[10 15 20],Po);
% [gvol,clusters] = group_outline2(EXPT.subjects,.05,5,5,[10 15 20],Po,{'b' 'r' 'y'});


if any(gt - sort(gt)), error('gt must be in ascending order.'),end

if length(varargin) > 0, 
    bgimg = varargin{1};, 
else
    bgimg = spm_get(1,'*img','Choose anatomical overlay image.'); % '/usr/private/spm99/templates/scalped_avg152T1.img';
end
bgV = spm_vol(bgimg);

if length(varargin) > 1,
    mycol = varargin{2};
else
    mycol = {'g' 'b' 'r' 'y' 'm' 'c' 'k'};
end

% ---------------------------------------------------------------
% * loop through subjects, get results, make outline, plot
% ---------------------------------------------------------------

index = 1;
for snum = mysbs
	
	SubjCode = snum{1};
  	if index == 1, cd([SubjCode]), else  cd(['../' SubjCode]), end

    
    
    % -------------------------------------------------------------
    % * process for individual subject
    % -------------------------------------------------------------

    % get the results for this subject
    % - - - - - - - - - - - - - - - - - - - - - - - - - -
	TOR.u = u;
	TOR.k = k;
	TOR.Ic = wcon;
	TOR.resdir = pwd;
	TOR.maskWOtherContrasts = 0;
	TOR.multCompCorrect = 'uncorrected'; %'FDR'; 
	[hReg,SPM,VOL,xX,xCon,xSDM] = tor_spm_results_ui(TOR);
    svol = voxel2mask(SPM.XYZ',VOL.DIM');
    mytit{index} = SPM.title;
     
    % check VOL mat files - all subs should be the same
    if index > 1,
        if any(VOL.M - oldVOL.M), 
            error([snum{1} ': Subject VOL.M field differs from previous subject!'])
        end
    end
    oldVOL = VOL;
    
	% add subject to group mask	
	% --------------------------------------------------------------
	if index == 1, gvol = svol;, else,  gvol = gvol + svol;, end

	index = index + 1;

end

cd ..

% -------------------------------------------------------------
% * check titles to make sure they're all the same.
% -------------------------------------------------------------
testdt = str2mat(mytit);
if length(strmatch(testdt(1,:),testdt)) < size(testdt,1)
    warning('Contrast titles do not all match!!!')
end

% -------------------------------------------------------------
% * Make clusters out of mask, for compatibility with montage_clusters
% -------------------------------------------------------------
gt(end+1) = Inf;
for i = 1:length(gt)-1
    
    [x,y,z] = ind2sub(size(gvol),find(gvol >= gt(i) & gvol < gt(i+1)));
    XYZ = [x y z]';
    clusters{i}.XYZmm = voxel2mm(XYZ,VOL.M);
    clusters{i}.title = mytit{1};
    clusters{i}.threshold = gt(i);
    
end
  
% -------------------------------------------------------------
% * Write image file of gvol (borrowing hdr info from an existing con img)
% -------------------------------------------------------------

if i < 10, myz = '000';, elseif i < 100, myz = '00'; else, myz = '0';,end
V = spm_vol(fullfile(TOR.resdir,['con_' myz num2str(i) '.img']));
V.fname = ['group_con' num2str(i) '.img'];
V.descrip = ['u=' num2str(TOR.u) ' k=' num2str(TOR.k) ' corr=' TOR.multCompCorrect ' dir=' TOR.resdir];

spm_write_vol(V,gvol);

% -------------------------------------------------------------
% * Image montage with montage_clusters.m
% -------------------------------------------------------------
str = ['montage_clusters(bgimg,clusters{1}'];
for i = 2:length(gt)-1
    if ~isempty(clusters{i}.XYZmm)
    	str = [str ',clusters{' num2str(i) '}'];
    end
end
% 0 is for suppressing overlap plotting, mycol is color specification
str = [str ',0,mycol);'];

if ~isempty(clusters{1}.XYZmm)
	eval(str)
else
	disp(['Con ' num2str(wcon) ': No voxels significant for ' num2str(gt(1)) ' participants.'])
end

return
