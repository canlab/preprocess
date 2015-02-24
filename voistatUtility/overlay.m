% inputs
slices = 2:3;
clim = [-5 5];
mapim = a;
anatim = a/2;
thresh = 3;


% set image scale factors
xscale = size(anatim,2) / size(mapim,2);
yscale = size(anatim,1) / size(mapim,1);
ptsize = min([xscale yscale]);

figure; cmap = hsv(1000); colormap(gray)





diff = clim(2) - clim(1); if diff < 0, error('clim must be [min max], where max > min'),end
cindex = (clim(1):diff/(size(cmap,1)-1):clim(2))';

for i = slices
   
   slice = mapim(:,:,i);   
   
   
   for k = 1:size(slice,1)
     for m = 1:size(slice,2)
        if (abs(slice(k,m)) > thresh), 
           newslice(k,m) = slice(k,m);
           
           % if value is greater or less than min/max of index, set equal to min/max
           if round(cindex(end,1)*100) < round(slice(k,m)*100)
              slice(k,m) = cindex(end,1);
           elseif round(cindex(1,1)*100) > round(slice(k,m)*100)
              slice(k,m) = cindex(1,1);
           end
           
           plot(slice(k,m),'Color',cmap(round(cindex*100) == round(slice(k,m)*100),:),'MarkerSize',ptsize)
        else newslice(k,m) = 0;
        end
     end
   end


end

% modify readim2 to take slices
