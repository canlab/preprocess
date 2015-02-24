function cluster_subpeaks(clusters)
% cluster_subpeaks(clusters)
% Creates a surface plot of 3-D activation blobs in clusters structure.
% Older function; not much used.

for mycl = 1:length(clusters)
   if clusters(mycl).numVox > 4 
       plc(clusters(mycl),mycl);
   end
end


return



function plc(cl,mycl)
   
% ---------------------------------------
% * find the modal Z value
% ---------------------------------------
%mmres = 5;  % the resolution, in mm, for defining a 'slice'
mmres = cl.voxSize(3);

makemov = 0;

zp = cl.XYZmm(3,:);
x = [min(zp):mmres:max(zp)];[n] = histc(zp,x);
zbest = x(n == max(n));
zbest = zbest(1);

% ---------------------------------------
% * find voxels
% ---------------------------------------
wh = zp > zbest-.5*mmres & zp < zbest+.5*mmres;
XYZ = cl.XYZmm(:,wh);
z = cl.Z(wh)';

x = XYZ(1,:);
y = XYZ(2,:);
[X,Y] = meshgrid(x,y);

% ---------------------------------------
% * define height (Z) surface
% ---------------------------------------
Z = zeros(size(X));

for i = 1:size(XYZ,2)
    
    whz = find(X > XYZ(1,i) - mmres*.5 & X < XYZ(1,i) + mmres*.5 & Y > XYZ(2,i) - mmres*.5 & Y < XYZ(2,i) + mmres*.5);
    Z(whz) = z(i);
    
end

%XYZo = mm2voxel(XYZ,cl,1);
%m = voxel2mask(XYZo,[size(X,1) size(X,1) 1],z);
%Z = m(:,:,unique(XYZo(:,3)));
%figure;imagesc(Z)

% ---------------------------------------
% * create surface plot
% ---------------------------------------
cla
surf(X,Y,Z)
xlabel('X','FontSize',18)
ylabel('Y','FontSize',18)
zlabel('Z-score','FontSize',18)
set(gca,'YTick',unique(y))
set(gca,'XTick',unique(x))
set(gcf,'Color','w')
title(['Cluster ' num2str(mycl)],'FontSize',18)

axis vis3d
for i = 0:1:360,view(i,30),drawnow,end

return