function [tID, tN, P, sigmask, XYZ] = tor_fdr(df,Timgname,varargin)
% function [tID, tN, P, sigmask, XYZ] = tor_fdr(df,Timgname,[opt] mask img name)
% if img names are empty, pick file names from list
% images must have the same voxel dimensions

% get list of files if argument is empty
% -----------------------------------------------------------------------
if isempty(Timgname),
    Timgname = spm_get(Inf,'*.img','Select a t-map img file (or more to get conjunction)');
end

firstName = deblank(Timgname(1,:));
[Timg,hdr] = readim2(firstName,'t');   % read in first image

moreTimg = [];
go = 0;

% calculate mask
% -----------------------------------------------------------------------
if length(varargin) > 0,
    maskimg = varargin{1};
    while ~go
        go = 1; 
        if isempty(maskimg), 
            maskimg = spm_get(1,'*.img','Select an explicit mask image (include voxels > 0)');
        end
        mask = readim2(maskimg,'t');
        if isempty(mask) | sum(sum(sum(mask))) <= 0, 
            go = 0; maskimg = [];
            disp('Mask is empty.')
        end
    end
else
    mask = ones(size(Timg));
end


% finish masking 1st image
% -----------------------------------------------------------------------
Timg(mask <= 0) = NaN;          % explicit mask (if supplied)
Timg(~isfinite(Timg)) = NaN;
Timg(Timg == 0) = NaN;          % mask out voxels outside brain (exactly 0 in 0-mask)


% read in other images, if any, for conjunction
% -----------------------------------------------------------------------
for i = 2:size(Timgname,1)
    nextTimgname = deblank(Timgname(i,:));
    moreTimg{i-1} = readim2(nextTimgname,'t');   % read in image
    moreTimg{i-1}(mask <= 0) = NaN;              % explicit mask (if supplied)
    moreTimg{i-1}(~isfinite(moreTimg{i-1})) = NaN;
    moreTimg{i-1}(moreTimg{i-1} == 0) = NaN;       % mask out voxels outside brain (exactly 0 in 0-mask)
end


% make P map for 1st image
% -----------------------------------------------------------------------
P        = 1-cdf('T',Timg(:),df);
P = reshape(P,hdr.xdim,hdr.ydim,hdr.zdim);

% if there is more than one p map, multiply them to give joint probability
% -----------------------------------------------------------------------
for i = 1:length(moreTimg)
    
    P2        = 1-cdf('T',moreTimg{i}(:),df);
    P2 = reshape(P2,hdr.xdim,hdr.ydim,hdr.zdim);
    P = P .* P2;
end

% get rid of NaN's in P-map  -- P is now a vector
% -----------------------------------------------------------------------
Pvec = P;
Pvec(isnan(Pvec)) = [];
Pvec = Pvec';

% get FDR t-threshold
% -----------------------------------------------------------------------
[pID pN] = FDR(Pvec,0.05);

tID      = icdf('T',1-pID,df);     % T threshold - ok for imaging
tN       = icdf('T',1-pN,df);      % T threshold, no correl. assumptions
      
if isempty(pID)
    disp(['Result: can''t use FDR (returns empty)'])
else
    sigmask = P < pID;
    XYZ = mask2voxel(sigmask);
    cl = spm_clusters(XYZ');
    disp(['Result at tcrit = ' num2str(tID) ': ' num2str(size(cl,2)) ' significant voxels in ' num2str(max(cl)) ' clusters.'])
end

return
      