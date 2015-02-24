function [clusters,subclusters] = mask2clusters(P,varargin)
%[clusters,subclusters] = mask2clusters(img mask file with voxels,[imgs to extract data from],[df])
%
% tor wager
% extracts clusters and con img data from mask
% use with mask_intersection.m
%
% to get clusters but not extract data, enter only one argument.
% to get clusters and choose extraction imgs with the GUI, enter an empty [] 2nd argument.
%
% DOES *NOT* CONVERT BETWEEN DIFFERENT VOXEL SIZES AND POSITIONS BETWEEN IMNAMES AND SPM/VOL STRUCTS
%
% see also roi_probe
%
% modification 2/27/03
% if no imgs are entered, Z-scores are values from mask
% if df is entered, values in mask img are converted to Z-scores with spm_t2z.m
% if extract img names are empty and df is entered, assume we're using values from mask as t-values
%   and convert to Z-scores
% 

if isempty(P), P = spm_get(1,'*img','Select mask image with voxels to extract');, end

if length(varargin) > 0,
    imP = varargin{1};
    if isempty(imP) & length(varargin) < 2, 
        cwd = spm_get(1,'*','Select images with data to extract');, 
    end
else
    imP = [];
end

if length(varargin) > 1,
    df = varargin{2};
else
    df = [];
end

resl = 0;


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% get coordinates and height values from the mask
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

[CLU.XYZ,CLU.XYZmm,CLU.Z,CLU.V] = img2voxel(P);
CLU.XYZ = CLU.XYZ(1:3,:);

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% change pseudo-T  or T values to Z scores
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if ~isempty(imP)
    if isempty(df) & length(varargin) > 1, df = size(imP,1) - 1;, end     % to get df
    if df == 0, df = [];, end
end

if ~isempty(df)
    disp(['Converting to Z scores based on ' num2str(df) ' df.'])
    [CLU.Z] = spm_t2z(CLU.Z,df);
else
    df = 0; 
    disp('Saving values in mask file in clusters.Z (no z-score conversion)')
    %CLU.Z = ones(1,size(CLU.XYZ,2));
end

if size(CLU.Z,1) > size(CLU.Z,2), CLU.Z = CLU.Z';, end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% fill in other fields
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% adjust crit_t from 0 to 1 for compatibility with fixed_TSU (Talairach)
% needs to be a non-zero value to display on Talairach
% ------------------------------------------------------
CLU.crit_t = 1;

CLU.cl_size = 0;
CLU.df = df;

% for compatibility with SPM struct and cluster analysis
% ------------------------------------------------------
CLU.voxSize = diag(CLU.V.mat)'; 
CLU.voxSize = CLU.voxSize(1:3);
CLU.VOX = CLU.voxSize;		% compatible with VOL structure
CLU.M = CLU.V(1).mat;

CLU.u = CLU.crit_t;
CLU.k = CLU.cl_size;

[a,b,c] = fileparts(CLU.V.fname);
CLU.title = ['snpm_' b];

if resl, imP = spm_get(Inf,'*img','Select imgs for data extraction.');
else, % input imP              % to get P
end

%if ~exist(P(1,:))
%    warning('Image files cannot be found in ind. subject directories: Path has changed?')
    % tries to re-find P - only works for my specific directory 
%    for i = 1:size(P,1)
%        [d,f,e] = fileparts(P(1,:));
%        [dummy,d] = fileparts(d);
%        newP{i,1}=fullfile('..',d,[f e]);
%    end
%    P = cell2mat(newP);
%    disp(['Found: ' P(1,:) ' etc.'])
%end

[clusters] = tor_extract_rois(imP,CLU,CLU,1);
for i = 1:length(clusters), 
    clusters(i).P = P;, 
    clusters(i).imP = imP;, 
    if size(imP,1) == 1 & df == 0,
        %disp('Saving values in mask file in clusters.Z')
        clusters(i).Z = clusters(i).all_data;
    end
end

if nargout > 1, [subclusters] = cluster_princomp(clusters);, end
return