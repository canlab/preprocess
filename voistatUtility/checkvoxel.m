function checkvoxel(coord,redT,TMAP,ts,ROI,pointh,varargin);
% function checkvoxel(coord,redT,tmap,ts,ROI,vox.pointh);
%
% This function is for use with voistatgui, and requires variables assigned
% by the montage, getvoxels, and voxelsurf programs.
%
% Coord should be [x,y,z] voxel coordinates in functional images / T-map
% x is columns of img matrix, y is rows.
%
% TMAP is t-map in canonical orientation, whole thing.
% ts is timeseries to check against.
%
% varargin is for batch mode, where you pass it the dir name.

% CHECKING STUFF
% -----------------------------------------------------------------------------------------
% plot a blue circle on the plot to check which point you actually got the timeseries from.
if size(coord,1) == 1 & strcmp(redT.orient,'ax')	% only if single voxel check in axial.
	if pointh > 1, delete(vox.pointh),end
	axes(redT.hand{1}(redT.originalslice{2} == coord(3)))
	hold on
	if strcmp(redT.flipy,'yes')
		plot(1 + redT.hdr{2}.xdim - coord(1),coord(2),'ro');
	else
		plot(coord(1),coord(2),'ro');
	end
end

% make a slice figure
figure;subplot(1,2,1);imagesc(TMAP(:,:,coord(1,3)));colorbar('horiz')
hold on; 
for i = 1:size(coord,1)
	if coord(i,3) == coord(1,3)
		plot(coord(i,1),coord(i,2),'k.','MarkerFaceColor','k');
	end
	if size(coord,1) == 1
		plot([coord(1) coord(1)],[1 size(TMAP,1)],'k')
		plot([1 size(TMAP,2)],[coord(2) coord(2)],'k')	
	end
end
title('T-map')

% check the ts against values you load in.
check(1:5,1) = ts(1:5);
for i = 1:5
   if nargin >6
      disp('	checkvoxel batch mode.'); P = getfiles(varargin{1}); %([ROI.imgname '000' num2str(i)]); 
   else
      P = getfiles(ROI.imgname); %([ROI.imgname '000' num2str(i)]);
   end
	[b,hdr] = readim2([P{i}(1:end-4)]); check(i,2) = b(coord(1,2),coord(1,1),coord(1,3));
end

% plot the functional slice
subplot(1,2,2);imagesc(b(:,:,coord(1,3)));colorbar('horiz')
hold on; 
for i = 1:size(coord,1)
	if coord(i,3) == coord(1,3)
		plot(coord(i,1),coord(i,2),'k.','MarkerFaceColor','k');
	end
	if size(coord,1) == 1
		plot([coord(1) coord(1)],[1 size(TMAP,1)],'k')
		plot([1 size(TMAP,2)],[coord(2) coord(2)],'k')	
	end
end
title('Functional image 0005')
set(gcf,'Position',[298    92   776   420])
disp('your ts   from files')
check
% -----------------------------------------------------------------------------------------