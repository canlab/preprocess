function [voxels,h] = getVoxels2(h,mask)
% function [voxels,h] = getVoxels2(h,mask)
%
% simplified.
% 1) must image all slices on montage, in order from 1:end, bottom to top
% 2) all dimensions and voxel sizes must match between overlay images on plot
%
% returns voxel coordinates from boxes drawn on image with handles h.
% voxels returned are masked by mask of 0's and 1's.


% --------------------------------------------------------------------------
% * get voxels in plot
% --------------------------------------------------------------------------
button = 1;
while ~(any(button == 2))
	[morevox,button] = voxelcrunch(h,mask);
	voxels = [voxels;morevox];
end


% --------------------------------------------------------------------------
% * get rid of duplicates
% --------------------------------------------------------------------------
i = 1;	
while i < size(voxels,1)
	j = i+1;
	while j <= size(voxels,1) 
		if sum(voxels(i,:) == voxels(j,:)) == 3
			voxels(j,:) = [];
		else
			j = j + 1;
		end
	end
	i = i+1;
end

disp(['Selected ' num2str(size(voxels,1)) ' voxels TOTAL.'])
return








% --------------------------------------------------------------------------
% * voxelcrunch: sub function: gets list of voxels from current axis
% --------------------------------------------------------------------------

function [voxels,button] = voxelcrunch(h,mask)

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

% define z (slice)
zrange = 1:length(h);       % zrange is index of plot handles
z = zrange(gca == h);       % z is which plot index number is currently used

% make output - [col,row,slice], high is right, high is bottom
voxels = [x y];
voxels(:,3) = z;


% * mask voxels
% --------------------------------------------------------------------------
mask2 = voxel2mask(voxels,size(mask));
selectedmask = mask .* mask2;
voxels = mask2voxel(selectedmask);



% GRAPHICAL OUTPUT STUFF
% fill polygon
hold on
if button == 1
   color = 'ys';

elseif button == 3
   color = 'gs';
elseif button == 2
   color = 'cs';
else color = 'ms'
end

hold on
for i = 1:size(voxels,1)
   plot(voxels(i,1),voxels(i,2),color,'MarkerSize',3,'MarkerFaceColor',color(1))
end
drawnow


% ==== get rid of duplicates =====
i = 1;	
while i < size(voxels,1)
	j = i+1;
	while j <= size(voxels,1) 
		if sum(voxels(i,:) == voxels(j,:)) == 3
			voxels(j,:) = [];
		else
			j = j + 1;
		end
	end
	i = i+1;
end
	
disp(['Selected ' num2str(size(voxels,1)) ' voxels.'])
return


        