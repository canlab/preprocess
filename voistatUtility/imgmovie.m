function [M,slicemovie,Cube] = imgmovie(filenames,slice,varargin)
% function M,slicearray,Cube = imgmovie(filenames,slice,skip)
%
% filenames = cell array of image file names
% slice = which slice (z) you want to make a movie of
% nimages = number of images in series, total
% skip = optional.  How many images to skip btwn samples in your movie. default = 1;

%whos filenames
if nargin < 3, skip = 1;
else skip = varargin{1};
end

% offset = hdr.xdim * hdr.ydim * slice;

% get the header for the basename, use this info for all images
basename = filenames{1}(1:end-4);
eval(['hdr = read_hdr(''' basename '.hdr'');'])

%% vvvvvv
% figure out the format of the data file from the header.
% only three types are supported right now ...
switch hdr.datatype     
   case 0
      fmt = 'int8';
   case 2
      fmt = 'uint8';
   case 4
      fmt = 'short';
   case 8
      fmt = 'int';
   case 16
      fmt = 'float';
   case 32
      fmt = 'float';
      xdim = hdr.xdim * 2;
      ydim = hdr.ydim * 2;

   otherwise
         errormesg(sprintf('Data Type %d Unsupported. Aborting',hdr.datatype));
         return
end

index = 0;
for JJ = 1:skip:size(filenames,1),index  = index+1;,end
disp(['reading ' num2str(index) ' image files...'])
%skip
%size(filenames,1)


index = 1;
skipbytes = hdr.xdim * hdr.ydim * slice-1; 

for i = 1:skip:size(filenames,1)
   eval(['[pFile,messg] = fopen(''' filenames{i} ''',''r'');']) 
   outnames{index} = filenames{i};
   if pFile == -1
      disp(messg);   
      return;
   end  
   status = fseek(pFile,skipbytes,-1);
	slicemovie(:,:,index) = (fread(pFile,[hdr.xdim, hdr.ydim], fmt)) ;
   fclose(pFile);
   if mod(i,50) == 0, disp(['		done ' num2str(i)]), end
   index = index+1;
end

disp('calculating global cube value...')
Cube = nanmean(nanmean(slicemovie(10:30,10:30,:)));

disp('building movie...')
figure;
%set(gcf,'Position',[46 16 1106 916]);
colormap(gray)
axis off
axis square
jkl = 0;
for i=1:size(slicemovie,3)   imagesc(slicemovie(:,:,i));
		text(5,5,outnames{i});
		M(i) = getframe;
end
disp('Done!  type movie(your movie variable, output from this function) to see the movie.')

%%%%^^^^^   
   
  % return
