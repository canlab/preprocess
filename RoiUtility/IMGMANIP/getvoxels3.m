function [voxels,outmask] = getvoxels3(h,mask,wslices,varargin)
% [voxels,outmask] = getvoxels3(h,mask,wslices,[opt] 'sagg')
% h is handle of axis images
% mask is binary 1/0 mask
% wslices is image handle index to volume slice mapping
%
% Tor Wager
% 
% press 1st button to select, press 2nd button after last slice choice
% mask should be original mapped volume, bottom right corner is 1st voxel.
%
% optional 4th argument based on saggital slices
% where 0,0,0 is bottom right back corner
% slices are right to left, in order
%
% functions called:
% C:\matlabR12\toolbox\matlab\elmat\ind2sub.m

if nargin > 3, dosagg = 1;, else, dosagg = 0;, end
button = 1;
voxels = [];
outmask = zeros(size(mask));

while button == 1

    if dosagg
        [newv,newmask,button] = voxelcrunch(h,mask,wslices,'sagg');
    else
        [newv,newmask,button] = voxelcrunch(h,mask,wslices);
    end
    voxels = [voxels; newv];
    outmask(find(newmask)) = 1;
    
end

voxels = unique(voxels,'rows');
    

return




% --------------------------------------------------------------------------
% * voxelcrunch: sub function: gets list of voxels from current axis
% --------------------------------------------------------------------------

function [voxels2,mask2,button] = voxelcrunch(h,mask,wslices,varargin)

if nargin > 3, dosagg = 1;, else, dosagg = 0;, end

[x,y,button] = ginput(2);
   x = round(x);
   y = round(y);
   xy = [];
   for i = min(x): max(x)
      for j = min(y):max(y)
         xy = [xy;i j];
      end
   end
   x = xy(:,1);
   y = xy(:,2);

if dosagg
    gx = x; % g denotes graphic input x and y
    gy = y;
    y = gx;
    z = gy;
    z = size(mask,3) - z;       % reverse z
    % define x (slice)
    xrange = 1:length(h);       % zrange is index of plot handles
    x = xrange(gca == h);       % z is which plot index number is currently used
    x = wslices(x);
    x = size(mask,3) - x;       % reverse x
    voxels = [ones(length(y),1) .* x y z];
    
else
    
    % reverse y: bottom right of head is 1st voxel in Analyze format 
    gy = y;
    y = size(mask,2) - y;

    % define z (slice)
    zrange = 1:length(h);       % zrange is index of plot handles
    z = zrange(gca == h);       % z is which plot index number is currently used
    z = wslices(z);

    % make output - [col,row,slice], high is right, high is bottom
    voxels = [x y];
    voxels(:,3) = z;
end

% * mask voxels
% --------------------------------------------------------------------------
[mask2,voxels2] = fast_voxel2mask(voxels,mask);
%selectedmask = mask .* mask2;
%voxels = mask2voxel(selectedmask);


hold on
for i = 1:size(voxels2,1)
    if dosagg
        plot(voxels2(i,2),size(mask,3) - voxels2(i,3),'r.','MarkerSize',3,'MarkerFaceColor','y')
    else
        plot(voxels2(i,1),size(mask,2) - voxels2(i,2),'r.','MarkerSize',3,'MarkerFaceColor','y')
    end
end
drawnow

return



function [mask2,mskdvox] = fast_voxel2mask(voxels,mask)

    mask2 = zeros(size(mask));
    ch = sub2ind(size(mask),voxels(:,1),voxels(:,2),voxels(:,3));
    mask2(ch) = 1;
    mask2 = mask .* mask2;
    [a,b,c] = ind2sub(size(mask2),find(mask2(:) > 0)); 
    mskdvox = [a b c];
    
return
    