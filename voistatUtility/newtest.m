%axes(h(end)); hold off; axis on;cla
%hold on
%index = 1;
%for i = 1:size(cindex,1)
%   plot([cindex(i,1) cindex(i,1)],[-8 -3],'Color',cmap(i,:))
%	index=index+1;   
%end
%drawnow
%set(gca,'YTick',[])
%text(clim(1),-16,['Anatomical: ' anatimN])
%text(clim(1),-14,['Overlay: ' mapimN])
%text(clim(1),-12,['Range: ' num2str(clim(1)) ' to ' num2str(clim(2))])
%text(clim(1),-10,['Threshold: ' num2str(thresh)])
%daspect([diff/20 1 1])

axes(h(end)); hold off; axis on;axis square;cla
axis auto;
hold on
for i = 1:size(cindex,1)
   plot([cindex(i,1) cindex(i,1)],[-6 -3],'Color',colmp(i,:))  
end
drawnow
set(gca,'YTick',[])
text(clim(1),-14,['Anatomical: ' anatimN])
text(clim(1),-11,['Overlay: ' mapimN])
text(clim(1),-10,['Range: ' num2str(clim(1)) ' to ' num2str(clim(2))])
text(clim(1),-8,['Threshold: ' num2str(thresh)])
axis image
daspect([diff/20 1 1])


