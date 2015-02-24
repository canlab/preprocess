function clusters = roi_probe(varargin)
% function clusters = roi_probe([files to extract data from],[clusters or img mask file to define voxels],[behavior],[size threshold])
% Takes a mask defined by one clusters or img file
% And extracts images from img files you pick, or those defined by another clusters file
% Works now only with clusters file as 2nd input.
% 
% also see mask2clusters.m, a simpler version that just extracts clusters from a mask file.
%
% Tor Wager

diary off

size_t = 0;

P = [];
if length(varargin) > 0
    P = varargin{1};
    [dummy,dummy,V] = load_file(P(1,:));
end

if length(varargin) > 3
    size_t = varargin{4};
end

% -----------------------------------------------------------------------------
% * get mask coordinates
% -----------------------------------------------------------------------------
XYZmm = [];

if length(varargin) > 1,
    f1 = varargin{2};
    [XYZmm] = load_file(f1,size_t);
end

while isempty(XYZmm)
    f1 = spm_get(Inf,'*','Choose mask clusters file or img file with mask');
    disp('Size mask ONLY implemented for image file input.')
    size_t = input('Enter cluster size threshold to impose on mask (can be 0): ');
    [XYZmm] = load_file(f1,size_t);

    if isempty(XYZmm), warning('Problem: no coordinates found - empty mask or no clusters.XYZmm in mat file?'), end
end

% -----------------------------------------------------------------------------
% * load data from images, if not input as argument to this function
% -----------------------------------------------------------------------------
while isempty(P)
    f2 = spm_get(Inf,'*','Choose cl file with img names or img files with data to extract');

    if findstr(f2,'.mat')
        [dummy,P,V] = load_file(f2);
    else
        P = f2;
    end

    if isempty(P), warning('Problem with .mat file selected: No files in clusters(1).imnames?');, end
end

% -----------------------------------------------------------------------------
% * extract data from P from mm coordinates in XYZ
% -----------------------------------------------------------------------------
SPM.title = f1;
SPM.u = 0;
vs = diag(V.mat)';
VOL.VOX = vs(1:3);
VOL.M = V.mat;
SPM.Z = ones(1,size(XYZmm,2));

SPM.XYZ = mm2voxel(XYZmm,VOL)';
SPM.XYZmm = voxel2mm(SPM.XYZ,VOL.M);

disp([num2str(max(spm_clusters(SPM.XYZ))) ' clusters found.'])
clusters = tor_extract_rois(P,SPM,VOL);

% -----------------------------------------------------------------------------
% * one-sample T-test on timeseries values (extracted data)
% -----------------------------------------------------------------------------
for i = 1:length(clusters)
    m = nanmean(clusters(i).all_data);
    se = nanstd(clusters(i).all_data) ./ sqrt(size(clusters(i).all_data,1));
    se(se == 0) = NaN;
    tmp = m ./ se; tmp(isnan(tmp)) = Inf;
    clusters(i).Z = spm_t2z(tmp,size(clusters(i).all_data,1)-1);
    clusters(i).Z(isinf(clusters(i).Z)) = -9999;
    m = nanmean(clusters(i).timeseries);
    se = nanstd(clusters(i).timeseries) ./ sqrt(size(clusters(i).timeseries,1));
    se(se == 0) = NaN;
    clt(i) = m ./ se;
end

% -----------------------------------------------------------------------------
% * re-extract with correct Z values and print table
% -----------------------------------------------------------------------------

SPM.Z = cat(2,clusters.Z);
SPM.XYZ = cat(2,clusters.XYZ);
SPM.XYZmm = cat(2,clusters.XYZmm);
clusters = tor_extract_rois(P,SPM,VOL);

if length(varargin) > 2
    [clusters] = simple_correl_clusters(clusters,varargin{3});
    
    for i = 1:length(clusters), 
    [h,p,ci,stats] = ttest(clusters(i).timeseries);
    r = corrcoef(clusters(i).timeseries,varargin{3}); r= r(1,2);
    [rci,sig,rZ,rp] = r2z(r,23,.05);
    
    clusters(i).str = sprintf('Ind: %d, [%3.0f, %3.0f, %3.0f], %d vox   C: u=%3.2f, t=%3.2f, p=%3.4f   R: r=%3.2f, Z=%3.2f, p=%3.4f', ...
        i,clusters(i).mm_center(1), clusters(i).mm_center(2),clusters(i).mm_center(3),clusters(i).numVox,nanmean(clusters(i).timeseries), stats.tstat, p, r, rZ, rp);
    
    end
    str2mat(clusters.str)

end

diary on

fprintf(1,'\n\n%s\n',clusters(1).title)
fprintf(1,'Images(1): %s\n',clusters(1).imnames(1,:))
cluster_table(clusters);
fprintf(1,'\nCluster one-sample t-tests overall t-value\n')
fprintf(1,'* = p < .05, one tailed.\n')
fprintf(1,'Cluster\tt\tp\tSig\n')
clt = [1:length(clt); clt; 1 - tcdf(clt,size(clusters(1).timeseries,1))];
for i = 1:length(clusters)
    fprintf(1,'\t%3.0f\t%3.3f\t%2.6f\t',clt(:,i))
    if clt(3,i) < .05,fprintf(1,'*');,end
    fprintf(1,'\n')
end

diary off
%for i = 1:length(clusters), fprintf(1,'%3.2f\t',clt(i)), end


return





% -----------------------------------------------------------------------------
% * sub-functions
% -----------------------------------------------------------------------------
    
function [XYZmm,P,V] = load_file(f1,varargin)


XYZ = [];
XYZmm = [];
P = [];

if findstr(f1,'.img')
    % it's an image file - load and get coords
    V = spm_vol(f1);
    vol = spm_read_vols(V);
    
    if length(varargin) > 0,
        size_t = varargin{1};
        if size_t > 0, vol = clusterSizeMask(size_t,vol);, end
    end

    [x,y,z] = ind2sub(size(vol),find(vol > 0));
    XYZmm = voxel2mm([x y z]',V.mat);
    
elseif findstr(f1,'.mat')
    eval(['load ' f1])
    
    if length(varargin) > 0,
        size_t = varargin{1};
        a = cat(1,clusters.numVox) < size_t;
        clusters(a) = [];
    end
    
    if isempty(clusters), error('No clusters, clusters is empty, or no clusters survive extent threshold!'), end
    
    if isfield(clusters,'imnames')
        P = clusters(1).imnames;
        V = spm_vol(P(1,:));
    else
        V = [];
    end
    XYZmm = cat(2,clusters.XYZmm);
    
    
else
    error('File must be a .mat file containing CLU or an .img mask file')
end