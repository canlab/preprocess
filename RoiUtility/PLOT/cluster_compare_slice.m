function cluster_compare_slice(varargin)
% cluster_compare_slice(cl1,cl2,etc.)
%
% This function takes a series of clusters structures and overlays
% activations on a series of slices, with the same scale used in each plot
% so that activation z-scores can be compared across the clusters
% structures.  Useful for showing comparisons between different analysis
% methods or different contrasts.
%
% tor wager
% 11/3/03

for i = 1:length(varargin)
    
    sxyz{i} = cat(2,varargin{i}.XYZ);
    try
        sz{i} = cat(2,varargin{i}.Z);
    catch
        sz{i} = cat(1,varargin{i}.Z)';
    end
    if i == 1, VOL.M = varargin{i}(1).M;,end
    
end
   
tmp = cat(2,sxyz{i});

VOL.dim = [1 1 max(tmp(3,:))];

ovlP = which('scalped_single_subj_T1.img');
compare_slice(ovlP,sxyz,sz,VOL);

return