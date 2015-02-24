function Pw = ind_model_fit(P,S,DX,nsess,varargin)
% single-subject timeseries extraction and model fitting
% saves full-model betas (model fits) and nonlinear fits to betas.
%
% P     = image file names (str matrix)
% S     = filtering matrix (e.g., high-pass)
% DX    = full model to fit, unfiltered
% vb    = [optional] verbose output level: 0 none, 1 some, 2 lots
% mask  = [optional] mask 3-D volume to apply before extracting
% nsess = number of sessions (intercepts): assumes last nsess columns
%         of DX are run intercepts to be removed before trimming
%
% dims  = dimensions of image files in data
% cols  = no. of columns of DX, also no. of beta images to write
% SDX   = smoothed (filtered) full model to fit
% PSDX  = pseudoinverse of filtered full model
% PSDXS = PSDX * S, ready to multiply with data y for pinv(SX) * Sy
% betas = 4-D array of x,y,z,column
%
% ntrimmed = number of outliers trimmed from timeseries at each voxel
% Pw    = string matrix of output file names

% -------------------------------------------------------------------
% * configure input arguments
% -------------------------------------------------------------------

global dographics
dographics = 1;

if dographics, figure;imagesc(DX);,colormap(gray); drawnow; end

vb = 2; if length(varargin) > 0, vb = varargin{1};, end
if vb > 0, t1 = clock;, fprintf(1,'\n\t\tSubject setup: %s',pwd),end

doridge = 0;
if length(varargin) > 0, doridge = varargin{2};, end

V = spm_vol(P(1,:));

global dims
global cols

% ridge regression option: modify diagonals here
SDX = S * DX;
if doridge
    k = .0001;
    PSDX = inv(SDX' * SDX + (eye(size(SDX,2)) * k)) * SDX';   % pinv, adding ridge constant k
else    
    PSDX = pinv(SDX);
end
PSDXS = PSDX * S;

dims = V.dim(1:3);
cols = size(SDX,2);
betas = NaN * zeros([dims cols]);


if vb > 0, fprintf(1,'\n\t\tbeta dims: %3.0f %3.0f %3.0f %3.0f ',dims(1), dims(2), dims(3), cols),end

mask = ones(dims); 
if length(varargin) > 1, 
    mask = varargin{2};, 
    if isstr(mask), Vm = spm_vol(mask);, mask = spm_read_vols(Vm);,end
end

ntrimmed = NaN * zeros(dims);

if vb > 0, fprintf(1,'\n\t\tFinished in %3.0f s',etime(clock,t1)),end

if dographics, close, end

% -------------------------------------------------------------------
% * for each slice...
% -------------------------------------------------------------------

for slicei = 1:dims(3)
    
    if vb == 2, t1 = clock;, fprintf(1,'\n\t\tSlice %3.0f ',slicei),end
    
    [a,b] = process_slice(slicei,P,S,PSDXS,DX,nsess,mask(:,:,slicei));
    try
        betas(:,:,slicei,:) = a;
    catch
        fprintf(1,'...!error updating betas...')
        whos betas
        whos a
    end
    ntrimmed(:,:,slicei) = b;
    
    if vb == 2, fprintf(1,'%6.0f trimmed, %6.0f s total',sum(sum(sum(~isnan(b)))),etime(clock,t1)),end
end

    
% -------------------------------------------------------------------
% * write beta images
% -------------------------------------------------------------------    
 
if vb > 0, t1 = clock;, fprintf(1,'\n\t\tWriting %3.0f beta images.',cols),end

for bi = 1:cols
    Pwb = write_beta(bi,V,betas(:,:,:,bi));
    if bi == 1, Pw = Pwb;, else, Pw = str2mat(Pw,Pwb);, end
end

% number of trimmed outliers image
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
V.fname = ['ntrimmed.img'];
V.descrip = ['Number of outliers trimmed from timeseries at each voxel.'];
spm_write_vol(V,ntrimmed);

if vb > 0, fprintf(1,'%6.0f s to dir: %s ',etime(clock,t1),fileparts(Pw(1,:))),end

return
    
    
    
    
    
    
    
    
    
function [betas,ntrimmed] = process_slice(slicei,P,S,PSDXS,DX,nsess,varargin)

global dographics 

mask = []; if length(varargin) > 0, mask = varargin{1};, end

global dims
global cols
% -------------------------------------------------------------------
% * load the slice
% -------------------------------------------------------------------
O.z = slicei;
et = clock;
if ~isempty(mask) & ~(sum(sum(sum(mask))) > 0)
    betas = NaN * zeros([dims(1:2) cols]);
    ntrimmed = NaN;
    fprintf(1,'...Empty slice...')
    return
else
    %sl = timeseries2('slice',P,O);
    sl = timeseries_extract_slice(V,slicei);
end
fprintf(1,'...loaded in %3.2f s...',etime(clock,et))
if ~isempty(mask), sl(:,:,1) = sl(:,:,1) .* mask;, end

if dographics, 
    clf
    subplot 221; imagesc(mask);,colormap(gray); title('Mask'),colorbar('horiz'),drawnow; 
    subplot 222; imagesc(sl(:,:,2));,colormap(gray); title('Slice timepoint 2'),colorbar('horiz'),drawnow; 
end

% -------------------------------------------------------------------
% * fit the model to each nonzero voxel
% -------------------------------------------------------------------
betas = NaN * zeros([dims(1:2) cols]);
nsize = size(sl); nsize = nsize(1:2);
ntrimmed = NaN * zeros(nsize);
wvox = find(sl(:,:,1) ~= 0 & ~isnan(sl(:,:,1)));
[i,j] = ind2sub(size(sl(:,:,1)),wvox);
fprintf(1,'%3.0f voxels...',length(i))

et = clock;
for k = 1:length(i)
    % b = pinv(S * X) * (S * y) do PSDX * S first
    % equivalent to: betas(i(k),j(k),:) = PSDX * (S * squeeze(sl(i(k),j(k),:)));
    [y,ntrimmed(i(k),j(k))] = trimts(squeeze(sl(i(k),j(k),:)),3,DX(:,end-nsess+1:end));
    betas(i(k),j(k),:) = PSDXS * y;   %for NaN replacement: PSDXS(:,~isnan(y)) * y(~isnan(y));
    if k == 100, fprintf(1,'%3.0f s per 100 vox...',etime(clock,et)), end
end

if dographics, 
    subplot 223; imagesc(ntrimmed);,colormap(gray); title('ntrimmed'),colorbar('horiz'),drawnow; 
    subplot 224; imagesc(betas(:,:,end));,colormap(gray); title('Last beta (intercept)'),colorbar('horiz'),drawnow; 
end

return



function Pw = write_beta(bi,V,betas)

if bi < 10, myz = '000';, elseif bi < 100, myz = '00';, else myz = '000';,end
V.fname = ['dx_beta_' myz num2str(bi) '.img'];
V.descrip = ['Beta image for column ' num2str(bi)];

spm_write_vol(V,betas);
Pw = which(V.fname);

return
