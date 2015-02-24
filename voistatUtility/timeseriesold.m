function ts = timeseries(cmd,basename,varargin)
% function ts = timeseries(cmd,basename,nimgs [opt],slice/voxel/mask [opt])
% *****************************************************
% cmd : 'slice' 'voxel' 'volume'
% nimgs : # of img files to get timeseries from, OR row vector of image numbers.
% basename : 'basename of img file, excluding number' OR matrix of images, time = last dim
% nimgs: optional - number of images to load.  enter if reading from .img files.
% varargin: 
% 		if 'slice' - slice number you want to acquire
%		if 'voxel' - [x y z] coordinates of voxel, y = row #, x = col. #
%		if 'volume' - mask to get average area; then returns single timeseries	
%       if 'multi' - finds each individual timeseries and the avg.
%
% 'all' entered as last argument finds the ts for ALL the voxels in the region
% - normally finds the average of voxels in the region.
%
% note: luis flips x and y.  this enters slices as [y,x].
%			this should matter only if x and y are different.
%
% ASSUMES BRAIN IS SIDEWAYS - so that x coord = col of matrix = y axis of actual brain.
% 'voxel' tested 1/29/01 for single_subj_T1.img.
% 'multi' tested same day, with s8t1.img (symmetrical x,y)
format compact;

% Set up arguments and build list of images if reading from files
% ***************************************************************

if isstr(basename)	% for files
   file = 1;
   %disp('input is a .img file - building list of image names')
   switch nargin
   case 3
      if nargin > 2 & size(varargin{1},2) == 1
         nimgs = varargin{1};
         domask = 0;
  		else error('must enter number of images as 3rd argument when reading .img files')
  		end
        
   case 4
       nimgs = varargin{1};
         switch cmd
            case 'slice'
               z = varargin{2};
            case 'voxel'
               coord = varargin{2};
            case 'volume'
               domask = 1;
               mask = varargin{2};
               mask = double(mask);
            case 'multi'
               coords = varargin{2};
            otherwise error('unknown command in first argument')
         end  % end switch
   otherwise, error('3 or 4 arguments only'),
   end
            
   % build cell array of image names
   % disp(['building list of images: ' num2str(nimgs)])
   % for fixed list of image numbers
   if size(nimgs,2) > 1
      disp('		Timeseries: fixed list of images')
      for i = 1:size(nimgs,2)
         if nimgs(1,i) < 10
            str{i} = ([basename '000' num2str(nimgs(i)) '.img']);
         elseif nimgs(1,i) < 100, str{i} = ([basename '00' num2str(nimgs(i)) '.img']);
         elseif nimgs(1,i) < 1000, str{i} = ([basename '0' num2str(nimgs(i)) '.img']);
         elseif nimgs(1,i) < 10000, str{i} = ([basename num2str(nimgs(i)) '.img']);
         else error('image numbers from 1 to 10000 only.')
         end
      end
		nimgs = size(nimgs,2);
   else
            % build list given a single number
      disp(['building list of images from 1 to ' num2str(nimgs)])
      for i = 1:9
   	    str{i} = ([basename '000' num2str(i) '.img']);
	end
   if nimgs > 99
		for i = 10:99
   		    str{i} = ([basename '00' num2str(i) '.img']);
		end
   else
        for i = 10:nimgs
   		    str{i} = ([basename '00' num2str(i) '.img']);
		end
   end
	if nimgs > 999
		for i = 100:999
   		    str{i} = ([basename '0' num2str(i) '.img']);
   	    end
   	    for i = 1000:nimgs
   		    str{i} = ([basename num2str(i) '.img']);
		end
	else
   	    for i = 100:nimgs
   		    str{i} = ([basename '0' num2str(i) '.img']);
		end
   end
   
	end

   imnames = str;, clear str;
   	if nimgs == 1, clear imnames, imnames{1} = [basename '.img'];,hdr=read_hdr([basename '.hdr']);
	else
		eval(['hdr = read_hdr(''' basename '0001.hdr'');'])
	end
   
	% check to see if all SPM scaling factors are 1!!!  And list them in spmscale.
	for i = 1:nimgs
		hdr = read_hdr([imnames{i}(1:end-3) 'hdr']);
		spmscale(i) = hdr.SPM_scale;
	end
	if ~(mean(spmscale) == spmscale(1)),warning('Timeseries: not all spm scale factors are the same!!!'),end
	
	% swap x and y dims!  xdim on image is ydim in matrix (rows)
	% now xdim, ydim, zdim are in terms of the image, which is sideways, not the brain.
	xdim = hdr.ydim;
 	ydim = hdr.xdim;
   	zdim = hdr.zdim;
   	
	
   switch hdr.datatype     
   case 0
      fmt = 'int8';
	  bitsperbyte = 8;
   case 2
      fmt = 'uint8';
	  bitsperbyte = 8;
	  bytespervoxel = 1;
   case 4
      fmt = 'short';
	  bitsperbyte = 8;
	  bytespervoxel = 2;
   case 8
      fmt = 'int';
   case 16
      fmt = 'float';
   		bitsperbyte = 8;
	  bytespervoxel = 4;   
   case 32
      fmt = 'float';
      xdim = hdr.xdim * 2;
      ydim = hdr.ydim * 2
   otherwise
         error(sprintf('Data Type %d Unsupported. Aborting',hdr.datatype))
   end  % end switch
   disp(['Timeseries: format is ' num2str(fmt)])
   onevox = bytespervoxel * bitsperbyte;
   disp(['bits per voxel = ' num2str(onevox) ', bytes per voxel = ' num2str(bytespervoxel) ', x,y dims = ' num2str(xdim) ' ' num2str(ydim)])
   
% Set up variables if matrix of images is entered directly   
% ***************************************************************   
else		% for matrix image input
   disp('input is a matrix of images - should be 4D.')
   file = 0;
   xdim = size(basename,2);
   ydim = size(basename,1);
   zdim = size(basename,3);
   nimgs = size(basename,4);
   switch cmd
            case 'slice'
               z = varargin{1};
            case 'voxel'
               coord = varargin{1};
            case 'volume'
               if nargin > 2
                  domask = 1;
                  mask = varargin{1};
                  mask = double(mask);
               else domask = 0;
               end
            case 'multi'
               coords = varargin{2};
            otherwise error('unknown command in first argument')
   end  % end switch
end	

clear ts;

% Do the computations 
% *************************************************************** 

switch cmd
   
case 'slice'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if nargin > 2
   	shift = bytespervoxel * (z-1) .* xdim .* ydim;
		%shift = onevox .* zshift;
   else error('must enter slice number as argument.'),end
   if file
   	for i = 1:nimgs
   		eval(['fid = fopen(''' imnames{i} ''',''r'',''b'');']) 
   		status = fseek(fid,shift,-1);
   		ts(:,:,1,i) = fread(fid,[ydim,xdim],fmt); 
         fclose(fid);
         	%if mod(i,50) == 0,disp(['done ' num2str(i) '...']),end
      end
  else 
      ts = basename(:,:,z,:);
  end 
   
case 'voxel'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['voxel = ' num2str(coord)])
disp(['dimensions = ' num2str(xdim) ' ' num2str(ydim) ' ' num2str(zdim)])
% tested 1/29/01 with square img matrix.
if nargin > 3,
	% yshift is in relation to the matrix; it's defined w/coord(1) = x bec. x in matrix is y (row) in file.
   zshift = (coord(3)-1) .* xdim .* ydim;
   yshift = (coord(1)-1) .* ydim   ;            % x should index # of columns
   xshift = (coord(2)-1)  ;                     % y should index # of rows
else error('must enter voxel coordinates as argument.'),end
if file
	shift = bytespervoxel * (zshift + yshift + xshift);	% shifts in bytes; 2 bytes per voxel
	for i = 1:nimgs
   		eval(['fid = fopen(''' imnames{i} ''',''r'',''b'');']) 
   		status = fseek(fid,shift,-1);
		if i == 1,disp(['reading bits ' num2str(ftell(fid)+1) ' - ' num2str(shift + onevox)]),end
		if ~(ftell(fid) == shift ),warning('coordinates exceeded end of file!'),end
   		ts(i,1) = fread(fid,1,fmt);
       fclose(fid);
       if mod(i,100) == 0,disp(['done ' num2str(i) '...']),end
    end  
	
else
   ts(:,1) = basename(coord(2),coord(1),coord(3),:);
end

case 'volume'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
if domask
   disp('masking image with mask variable.')
   if nargin == 5
       % 'all' doesn't work - too much memory!  try uint8 or some scalefactor
       disp('function will return timeseries for EACH voxel.')   
       %pack
   else
       disp('function will return timeseries with average of all voxels in mask.')
   end
   nvox = sum(sum(sum(mask)));
   for i = 1:size(mask,3)
      if sum(sum(mask(:,:,i))) > 0,maxslice = i;,end
   end
   for i = size(mask,3):-1:1;
      if sum(sum(mask(:,:,i))) > 0,minslice = i;,end
   end
   zshift = (minslice-1) .* xdim .* ydim;
   shift = 2 .* zshift;
else disp('no mask specified - finding all values')
   minslice = 1;
   maxslice = zdim;
   shift = 0;
end
zeroa = zeros(ydim,xdim);
for i = 1:nimgs
   if file
      eval(['fid = fopen(''' imnames{i} ''',''r'',''b'');'])
      status = fseek(fid,shift,-1);
      for j=minslice:maxslice 
	     a  = (fread(fid,[ydim,xdim], fmt)) ;
         array (:,:,j) = a;
      end
      for j = maxslice+1:zdim
         array(:,:,j) = zeroa;
      end
      fclose(fid);
      if mod(i,50) == 0,disp(['done ' num2str(i) '...']),end
   else
      array = basename(:,:,:,i);
   end
   if domask,array = array .* mask;,end
   if nargin == 5,ts(:,:,:,i) = array;              % ALL ts 
   else ts(i,1) = sum(sum(sum(array))) / nvox;      % AVG ts
   end
   clear array;
end   

case 'multi'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
for k =  1:size(coords,1)
	% checked for accuracy 1/29/01, along with 'voxel'
    coord = coords(k,:);
   	zshift = (coord(3)-1) .* xdim .* ydim;
   	xshift = (coord(1)-1) .* ydim;               % x should index # of columns
   	yshift = (coord(2)-1);                       % y should index # of rows
    if file
	    shift = bytespervoxel .* (zshift + yshift + xshift);
	    for i = 1:nimgs
   		    eval(['fid = fopen(''' imnames{i} ''',''r'',''b'');']) 
   		    status = fseek(fid,shift,-1);
			if ~(ftell(fid) == shift ),warning('coordinates exceeded end of file!'),end
   		    ts.indiv(i,k) = fread(fid,1,fmt);
            fclose(fid);
        end  
    else
        ts.indiv(:,k) = basename(coord(2),coord(1),coord(3),:);
    end
    disp(['done ' num2str(k) ' voxels. coord = ' num2str(coord(k,:))])
end
ts.avg = mean(ts.indiv,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end					% end switch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%