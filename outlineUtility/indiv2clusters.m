function indiv2clusters(infile,numsubs,EXPT,varargin)
% indiv2clusters(numsubs,EXPT,[opt] cluster size threshold)
% 
% by Tor Wager
%
% This function turns a gvol.mat file
% made with group_outline.m in outlineUtility
% into a clusters structure.
%
% gvol contains how many individual subjects
% activated at a particular voxel for a contrast.
% 
% Inputs:
% 1 infile is the ind_gvol*.mat file for the contrast
%   * indexes the number of the contrast
% 2 numsubs is the number of subjects required to
%   be included in the clusters structure
% 3 EXPT is the EXPT info structure from get_expt_info.m
%
% images extracted are con images in EXPT.SNPM.P

cl_size = 2;
if length(varargin) > 0, cl_size = varargin{1};,end

eval(['load ' infile])
gvol(gvol < numsubs) = 0;
disp(infile)

% get which contrast it is from the file name

for i = 1:length(infile), 
    if isreal(str2num(infile(i))) & ~isempty(str2num(infile(i))),
        a(i)=1;,
    else,a(i)=0;,
    end, 
end
connum = str2num(infile(find(a)));

% get con*image names
P = EXPT.SNPM.P{find(EXPT.SNPM.connums == connum)};
disp(P(1,:))
disp(P(2,:))
disp('...etc...')

% get mat file info
 V = spm_vol(deblank(P(1,:)));
 [d,f,e] = fileparts(infile);
 V.fname = [f '_masked.img'];
 
 % get CLU structure
 
    % Get XYZ pointlist from mask and cluster size mask if spec
	% ------------------------------------------------------
	if cl_size > 0
		[maskingImage,numClusters,XYZ] = clusterSizeMask(cl_size,gvol);
	else
		XYZ = mask2voxel(gvol);
		XYZ = XYZ';
	end

	% Get Z output - intensity values for sig. voxels
	% ------------------------------------------------------
	for i = 1:size(XYZ,2)
		% row is y, col is x
		Z(i) = gvol(XYZ(1,i),XYZ(2,i),XYZ(3,i));
	end

	% add other fields of CLU
	% ------------------------------------------------------
	CLU.maskingImage = maskingImage;
	CLU.cl_size = cl_size;
	CLU.crit_t = numsubs;
	CLU.numClusters = numClusters;
	CLU.XYZ = XYZ;
	CLU.Z = Z;
	CLU.XYZmm = voxel2mm(XYZ,V.mat);

	% for compatibility with SPM struct and cluster analysis
	% ------------------------------------------------------
	CLU.voxSize = EXPT.voxSize;
	CLU.u = numsubs;
	CLU.k = cl_size;
	CLU.title = V.fname;
	CLU.VOX = EXPT.voxSize;		% compatible with VOL structure
    CLU.M = V.mat;

% extract clusters
clusters = tor_extract_rois(P,CLU,CLU);

eval(['save ' f '_clusters CLU clusters'])
return
