function [tmask,deacttmask] = gettmask(tmapname, threshold)
% function tmask = gettmask(tmapname, threshold)
% takes an .img file - without the extension - 
% or a 3d array with the tmap in it.

%get masking image


	if isstr(tmapname)
      disp('Loading t-map from file...')
      tmap = readim2(tmapname,'p');
      disp(['making t-mask from ' tmapname ' at z >= ' num2str(threshold)])
   else
      disp('Using 3d t-map array variable...')
      tmap = tmapname;
      disp(['making t-mask from wkspace variable t-map at z >= ' num2str(threshold)])
   end
   
   for z = 1:size(tmap,3)
   		for x = 1:size(tmap,2)
           for y = 1:size(tmap,1)
               if tmap(y,x,z) >= threshold, tmask(y,x,z) = 1;,else tmask(y,x,z) = 0;,end
               if tmap(y,x,z) <= -threshold, deacttmask(y,x,z) = 1;,else deacttmask(y,x,z) = 0;,end
           end
      end
   end
   tmask = uint8(tmask);
   deacttmask = uint8(deacttmask);
   return
   