function P = fdr_conjunction(P,df,varargin)
% P = fdr_conjunction(P,df,[cluster size threshold])
% df is degrees of freedom [interest error], e.g. [1 13]
% Tor Wager 1/23/02

% cluster size for cluster struct, if entered
if nargin > 3, 
    bmask = varargin{2};
else
    bmask = [];
end

% cluster size for cluster struct, if entered
if nargin > 2, 
    cl_size = varargin{1};
else
    cl_size = 0;
end

% get list of files if argument is empty
% -----------------------------------------------------------------------
if isempty(P),
    P = spm_get(Inf,'*.img','Select a t-map img file (or more to get conjunction)');
end

n = size(P,1);
Vs = spm_vol(P);


% reslice images if they have different voxel sizes
% -----------------------------------------------------------------------
a = cat(1,Vs.dim);
if size(unique(a,'rows'),1) > 1         % if voxel sizes are different
    
    disp(['Reslicing images into space of 1st image: ' P(1,:)])
    reslice_imgs(deblank(P(1,:)),P(2:end,:),0)
    % adjust filenames
    for i = 1:size(P,1)
        newP{i} = P(i,:);
        if i > 1
            [p,f,e] = fileparts(P(i,:));
            f = ['r' f];
            newP{i} = fullfile(p,[f e]);
        end
    end
    P = str2mat(newP);
    Vs = spm_vol(P);
end


% calculate min img
% -----------------------------------------------------------------------
% minI = spm_imcalc_ui(P,'min.img',calce);  % trouble with more than 2 imgs?
vol = spm_read_vols(Vs(1));
for i = 2:length(Vs)
    vol = min(vol,spm_read_vols(Vs(i)));
end
V = Vs(1);
V.fname = 'min.img';
spm_write_vol(V,vol);


% FDR on min
% -----------------------------------------------------------------------
if isempty(bmask)
    % implicit 0 mask
    disp('Implicit zero masking')
    [u,Ps,Ts] = spm_uc_FDR(.05,df,'T',n,Vs,0);
else
    % explicit brain mask
    disp(['Explicit masking.'])
    %bmask = spm_vol(bmask);
    [u,Ps,Ts] = spm_uc_FDR(.05,df,'T',n,Vs,bmask);
end

% mask conjunction image at height threshold
% -----------------------------------------------------------------------
%vol = readim2('min.img','t');
[maskedImage maskingImage] = maskImg(vol,u,Inf);


% get cluster structure
% -----------------------------------------------------------------------
V = mask2struct('min.img',u,cl_size);


% get number of clusters and report
% -----------------------------------------------------------------------
numVox = size(V.XYZ,2);
if numVox > 0
    cls = spm_clusters(V.XYZ);
    for i = 1:max(cls)
        clsizes(i) = sum(cls == i);
    end
    
    % display results
    % -----------------------------------------------------------------------
    disp('Choose canonical image.')
    spm_image
    spm_orthviews('AddBlobs',1,V.XYZ,V.Z,V.mat)
    
else
    clsizes = 0;
end

[finalmask,numClusters] = clusterSizeMask(cl_size,maskingImage); % returns mask only
disp(['Result at tcrit = ' num2str(u) ': ' num2str(numVox) ' sig. voxels in ' num2str(numClusters) ' clusters of ' num2str(cl_size) ' or greater.'])
disp(['Cluster sizes are ' num2str(clsizes)])

% write filtered image
% -----------------------------------------------------------------------
Vs = spm_vol('min.img');
Vs.fname = 'filtered_min.img';
spm_write_vol(Vs,maskedImage);



return
