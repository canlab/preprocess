function [CLU,clusters] = spm2struct(cwd,RFX,varargin)
% [CLU,clusters] = spm2struct(directory name, EXPT.RFX, [opt] P, [opt] reslice)
%
% Given the directory of an SPM analysis,
% Creates a structure with relevant info
% Compatible with spm_add_blobs and cluster functions.
% must have res_mask.img -> mask of results.
%
% RFX structure contains info about the spm results to be obtained.
% fields with example values are below:
%
% RFX.maskWOtherContrasts = 0; ...mask with other contrasts?
% RFX.resdir = pwd;	...results directory
% RFX.Ic = 2;	...which contrast to test (2 is for one-sample t-test)
% RFX.multCompCorrect = 'uncorrected';	...or 'FDR', 'corrected'
% RFX.u = .001; ...height threshold (t or p)
% RFX.k = 10;	...extent threshold
%
% To Display:
% spm_image	(and choose anatomical) or spm_image('init','myname.img')
% spm_orthviews('AddBlobs',1,CLU.XYZ,CLU.Z,CLU.V.mat)
% spm_orthviews('AddColouredBlobs',1,CLU.XYZ,CLU.Z,CLU.V.mat,[0 0 1])
%
% optional: reslice flag (enter 0, 1, or 2)
% if selected, this option will allow you to 
% reslice the SPM filtered t into another space
% and choose new images to extract data from. 
% 
% entering 1 will choose new images only,
% and entering 2 will perform the reslice and get new images.
% 2 doesn't work here - for SnPM only.  see snpm2struct.m
% 
% Tor Wager, 2/12/02

if length(varargin) > 0, P2 = varargin{1};, end
if length(varargin) > 1, resl = varargin{2};,else resl = 0;,end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% go to the input directory and check for results files
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

eval(['cd(''' cwd ''')'])

if ~(exist('SPM.mat') == 2), error('Cannot find SPM.mat file.'), end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% get coordinates and height values from SPM
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
RFX.maskWOtherContrasts = 0;
RFX.resdir = pwd;
if isfield(RFX,'Ic')
    if isempty(RFX.Ic)
        RFX.Ic = 2;
    end
else, RFX.Ic = 2;
end
[SPM,VOL,xX,xCon,xSDM] = tor_spm_getSPM(RFX);
CLU.XYZ = SPM.XYZ;
CLU.XYZmm = SPM.XYZmm;
CLU.Z = SPM.Z;
CLU.V = VOL.M;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% change pseudo-T values to Z scores
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
CLU.df = SPM.df(2);
[CLU.Z] = spm_t2z(CLU.Z,CLU.df);


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% fill in other fields
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% adjust crit_t from 0 to 1 for compatibility with fixed_TSU (Talairach)
% needs to be a non-zero value to display on Talairach
% ------------------------------------------------------
CLU.crit_t = SPM.u;

CLU.cl_size = SPM.k;

% for compatibility with SPM struct and cluster analysis
% ------------------------------------------------------
CLU.voxSize = diag(CLU.V)'; 
CLU.voxSize = CLU.voxSize(1:3);
CLU.VOX = CLU.voxSize;		% compatible with VOL structure
CLU.M = CLU.V;

CLU.u = CLU.crit_t;
CLU.k = CLU.cl_size;

CLU.title = ['spm_' SPM.swd(end-6:end)];

if resl, P2 = spm_get(Inf,'*img','Select imgs for data extraction.');
elseif ~exist('P2')             % no img files - do not extract clusters here.
	[clusters] = tor_extract_rois([],CLU,CLU,1);
elseif isempty(P2), 
	P2 = spm_get(Inf,'*img','Select imgs for data extraction.');
else
        [clusters] = tor_extract_rois(P2,CLU,CLU,1);
        for i = 1:length(clusters), clusters(i).P = P2;, end
end


return