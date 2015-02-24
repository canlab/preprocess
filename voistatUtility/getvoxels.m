function [voxels,mask,h] = getvoxels(fun,nsliceinputs,varargin)
% function voxels = getvoxels(fun,nsliceinputs, EITHER: h(plothandles),tmask OR redT,whichmask)
%
% fun: 'single', 'square', 'poly'
% h is a vector of plothandles, with handles for each plot 1:i
% nsliceinputs is how many different slices / squares you want to get voxels from 
%		- used with 'poly' option only
%
% tmask [opt] is a 3-d mask file, or an spmT.img file 
%     if you specify a file name, 
%     uses default threshold of 3 to make mask
%
% redT, whichmask:
%   if you're imaging all functional slices and plot # = z value, there's no problem.
%   but if you have a restricted plotting range, you've rotated, or anything else with the montage program,
%   the coordinates you selected must be reconverted to be accurate.
%   ::: here redT is the output of montage4 after plotting, and whichmask is the index number of the mask
%       you want to use - 2 to n, where 2 is the first t-map in the T structure.
%   ::: you can leave h blank ([]).  the script reads it from the T structure.
%
% voxels = getvoxels('poly',[],2,redT,2);
%
%	To NOT MASK, leave off final argument; then it'll get ALL voxels in the region.
% voxels = getvoxels('poly',[],2,redT);
%
% output: list of [x y z] coordinates of voxels and
%		  a 3d mask for 'poly'
%
% voxels: [x,y,z]
%   x: x on plot, y in brain    should be 
%   y: y on plot, z in brain
%   z: plot,      x in brain

%get masking image
if nargin < 4, masked = 0;,
else masked = 1;
end
    if ~isstruct(varargin{1})
        rescale = 0; 
        h = varargin{1};
        if nargin < 3,
            tmapname = varargin{2};
        end
        if isstr(tmapname)
            disp('making mask from file, using default threshold of z = 3')
            tmask = gettmask(tmapname,3);
        else
   	        tmask = tmapname;   
        end
    elseif isstruct(varargin{1})
        T = varargin{1}; 
		figure(T.figh)
        if nargin > 3,whichmask = varargin{2};,else whichmask = 2; masked = 0;end
        [h,scalex,scaley,osliceplot,originalslice,whichplot,tmask,ori,flipy,maskdims] = setupT(T,whichmask);
        rescale = 1; 
    else
        error('incorrect argument spec.  1st arg must be either plot handles or T structure from montage.')
    end


%if iscell(h),h = h{1};,end

switch fun
   
case 'poly'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
voxels = [];
for in = 1:nsliceinputs
   if rescale == 1 & masked == 1
        morevox = getvoxels('square',1,T,whichmask);
    elseif rescale == 1 & masked == 0
        morevox = getvoxels('square',1,T);
    elseif masked == 1
   		morevox = getvoxels('square',1,h,tmask);
   else
      morevox = getvoxels('square',1,h);
   end
   voxels = [voxels;morevox];
end

% make output mask
if rescale
    disp('using montage output; will rotate mask to canonical orientation')
    if size(T.im{1}(3)) > size(h,2),warning('# of anatomical images exceeds # of plots! Results may not be accurate.'),end
    out = uint8(zeros(maskdims));
    for i = 1:size(voxels,1)
      out(voxels(i,2),voxels(i,1),voxels(i,3)) = 1;
    end  
    mask = rotateback(out,ori,flipy);
else   
    warning('no rescaling specified - make sure orientation and dimensions are same as original img files!!!')
    out = zeros(size(tmask,1), size(tmask,2), size(tmask,3));
   for i = 1:size(voxels,1)
      out(voxels(i,2),voxels(i,1),voxels(i,3)) = 1;
   end
   mask = out;
end

% Plot the selected voxels on the original orientation of the brain.
%[a,hdr,h] = readim2(rotateback(T.fullmask{1},ori,T.flipy),'ax','p');
%colormap(gray)
%scalex = size(a,2) / size(mask,2);
%scaley = size(a,1) / size(mask,1);
%plotmask(mask,h,'bs',[scalex scaley]);

%convert back from mask to voxel list in canonical orientation
voxels = [];
index = 1;
for i = 1:size(mask,2)
    for j = 1:size(mask,1)
        for k = 1:size(mask,3)
            if mask(j,i,k) == 1, voxels(index,:) = [i j k];,index = index+1;,end
        end
    end
end

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

disp(['Selected ' num2str(size(voxels,1)) ' voxels TOTAL.'])
return

case 'single'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get input from one slice
[x,y,button] = ginput;
x = round(x);
y = round(y);

%locate left vertex
minx = min(x);
minx = minx(1);
whichy = y(x == minx);
whichy = whichy(1);

%locate bottom vertex
miny = min(y);
miny = miny(1);
whichx = x(y == miny);
whichx = whichx(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'square'
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
end   % end switch   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define z (slice)
zrange = 1:size(h,2);       % zrange is index of plot handles
z = zrange(gca == h);       % z is which plot index number is currently used

disp(['plot index is 	' num2str(z)])
disp(['originalslice is ' num2str(originalslice(whichplot == z))])	% orig. slice for each slice with active voxels [20 21]
disp(['osliceplot is	' num2str(osliceplot(z))])					% orig. slice for each plot - [0 20 0 0 0 21]

% make output - [col,row,slice], high is right, high is bottom
voxels = [x y];
voxels(:,3) = z;

% convert voxel coords  - [col of funct.,row of funct.,slice of funct.] 
if rescale
    voxels = convert(voxels,scalex,scaley,osliceplot);
    voxels(voxels(:,3) == 0,:) = [];
end   


if masked == 1
	whos tmask
   mvoxels = [];
   for i = 1:size(voxels,1)
      if tmask(voxels(i,2),voxels(i,1),voxels(i,3)) == 1, mvoxels = [mvoxels; voxels(i,:)];,end
   end
   voxels = mvoxels;
end

% GRAPHICAL OUTPUT STUFF
% fill polygon
hold on
if button == 1
   color = 'yo';
   %fill(x,y,'b')
elseif button == 3
   color = 'go';
   %fill(x,y,'r')
else color = 'mo'
end

hold on
for i = 1:size(voxels,1)
   plot(voxels(i,1)*scalex,voxels(i,2)*scaley,color,'MarkerSize',3)
end
drawnow

if strcmp(fun,'single') & rescale
    % change voxels based on orientation
    % voxels are [x,y,z] = [col,row,slice]
    nslices = T.hdr{2}.zdim;
    switch ori
    case 'sagg'
        ny = T.hdr{2}.ydim;
        % columns = same (y on brain)
        if strcmp(flipy,'yes'),newvoxels(:,1) = 1+ny-voxels(:,1);,end        % flip back y (columns)
        % z on brain, z in canon. image = rows, flip z b/c small is top in image and bottom in brain
        newvoxels(:,3) = 1+nslices-voxels(:,2);
        % x in brain,rows in canon img. = plot, early rows = l side of brain on both imgs.
        newvoxels(:,2) = voxels(:,3);
    case 'ax'
		strcmp(flipy,'yes')
        % flipy reverses L and R.  flip rows of image.
        newvoxels = voxels;
        if strcmp(flipy,'yes'), nx = T.hdr{2}.xdim;,newvoxels(:,1)=1+nx-voxels(:,1); end
    case 'cor'
        newvoxels = voxels;
        warning('Single voxel rescaling not defined yet for coronal slices.  Coords are wrong!')
    otherwise warning('unknown orientation! voxel coordinates are probably wrong.')    
    end
    voxels = newvoxels;
end

	
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


% ======= sub-functions ========

function [h,scalex,scaley,osliceplot,originalslice,whichplot,tmask,ori,flipy,maskdims] = setupT(T,whichmask)
    h = T.hand{1};
    scalex = T.scalex(whichmask);
    scaley = T.scaley(whichmask);
    osliceplot = T.osliceplot{whichmask};
	originalslice = T.originalslice{whichmask};
	whichplot = T.whichplot{whichmask};
    tmask = T.fullmask{whichmask};
    ori = T.orient;
    flipy = T.flipy;
    maskdims = T.fullimdims{whichmask};
return

function coord = convert(coord,scalex,scaley,osliceplot)
    % apply in order: scale, rotate, flipy
    for i = 1:size(coord,1)
        coord(i,1) = round(coord(i,1) / scalex);
        coord(i,2) = round(coord(i,2) / scaley);
        coord(i,3) = osliceplot(coord(i,3));
    end
return

function out = rotateback(mask,ori,flipy)
    switch ori
    case 'sagg'
        disp('rotating mask from saggital to axial')
        mask = flipdim(mask,1);
        out = ipermute(mask,[3 2 1]);
    case 'cor'
        disp('rotating mask from coronal to axial')
        mask = flipdim(mask,1);
        out = ipermute(mask,[3 1 2]);
    otherwise out = mask;
    end
    if strcmp(flipy,'yes')
        disp('flipping back y dimension of img (columns of array).')
        out = flipdim(out,2);
    end
return
        