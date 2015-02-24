function [h,coord] = ortho(in,varargin)
% function [h,coord] = ortho(in,coord [all remaining are opt],h,voxsize,'get')
%
% in is a 3d matrix containing an image file
%
% h is optional figure handles (4) for existing figure
%   - 1st is figure h, 2nd - 4th are axis h's.  row vector.
% scalefactors is aspect ratio to adjust data by
%
% if 5th argument = 'get', get coordinate to use.
% !!!Enter Coordinates as x,y,z (col,row,slice of array)!!!
%

if nargin > 4,
    if strcmp(varargin{4},'get'),
        getcoord = 1;
    else getcoord = 0;
    end
else getcoord = 0;
end
if nargin > 3, voxs = varargin{3};,else voxs = [3.12 3.12 4];,end
if isempty(voxs),voxs = [3.12 3.12 4];,end
if nargin > 2 & not(isempty(varargin{2}))
    h = varargin{2};,
    if size(h,2) > 4,for i = 5:size(h,2),delete(h(i)),end,end
else h = [];,end
if nargin > 1 & not(isempty(varargin{1})), 
    coord = varargin{1};
else coord = [36 36 14];
end
clim = [min(min(min(in))) max(max(max(in)))];
yscale = voxs(2) / voxs(3);
zsize = size(in,3);

% swap coordinates to read row,col,slice (y,x,z) instead of x,y,z 
a = coord(1); coord(1) = coord(2);coord(2) = a;


% create figure if necessary
if isempty(h)
    h(1) = figure; colormap(gray)
    set(gcf,'Position',[233   402   544   195])
    h(2) = subplot(1,3,1);
    h(3) = subplot(1,3,2);
    h(4) = subplot(1,3,3);
elseif size(h,2) >= 4
    figure(h(1))
else error('input vector of 4+ object handles or []')
end

% plot views - axial
axes(h(2))
slice = in(:,:,coord(3));
imagesc(slice,clim)
axis image
axis off
% plot [miny maxy x x] [y y minx maxx]
hold on; h(5) = plot([0 size(in,2)],[coord(1) coord(1)],'r');
h(6) = plot([coord(2) coord(2)],[0 size(in,1)],'r');

axes(h(3))
im = rotim('ax2sagg',in);
slice = im(:,:,coord(1));   % take everything at selected i
imagesc(slice,clim)
asp = ([1 yscale 1]);  % ydim / zdim
daspect(asp);
axis off
% plot [miny maxy z z] [y y minz maxz]
hold on; h(7) = plot([0 size(in,2)],[1+zsize-coord(3) 1+zsize-coord(3)],'r');
h(8) = plot([coord(2) coord(2)],[0 size(in,3)],'r');

axes(h(4))
im = rotim('ax2cor',in);
%coord = rotim('ax2cor',coord)
slice = im(:,:,coord(2));   % everything at selected j
imagesc(slice,clim)
asp = ([1 yscale 1]);  % ydim / zdim
daspect(asp);
axis off
% plot [minx maxx z z] [x x minz maxz]
hold on; h(9) = plot([0 size(in,1)],[1+zsize-coord(3) 1+zsize-coord(3)],'r');
h(10) = plot([coord(1) coord(1)],[0 size(in,3)],'r');

if getcoord, 
    coord = getcoordinate(h,coord);
    [h,coord] = ortho(in,coord,h,voxs);
end

% swap coordinates back to read x,y,z instead of row,col,slice (y,x,z)
a = coord(1); coord(1) = coord(2);coord(2) = a;


return




function coord = getcoordinate(h,coord)
    [x,y] = ginput(1);
    x = round(x);,y = round(y);
    if gca == h(2)  % axial
        coord = [x y coord(3)];
    elseif gca == h(3)
        coord = [coord(1) y x];
    elseif gca == h(4)
        coord = [y coord(2) x];
    else error('unknown axis')
    end
return




function out = rotim(funct,in)
% function out = rotim(funct,in)
% coord is wrong; vol is ok.

if size(size(in),2) == 2,type = 'coord';
elseif size(size(in),2) == 3, type = 'vol';
else error('incompatible array size of input.')
end
    
switch funct
case 'ax2sagg'
    if strcmp(type,'vol')
        out = permute(in,[3 2 1]);
        for i = 1:size(out,3)
            out(:,:,i) = flipud(out(:,:,i));
        end
    elseif strcmp(type,'coord')
        out = [in(:,3) in(:,2) in(:,1)];
    end
    
case 'sagg2ax'
    if strcmp(type,'vol')
        for i = 1:size(in,3)
            out(:,:,i) = flipud(in(:,:,i));
        end
        out = permute(in,[3 2 1]);
    elseif strcmp(type,'coord')
        out = [in(:,3) in(:,2) in(:,1)];
    end
    
case 'ax2cor'
    if strcmp(type,'vol')
        out = permute(in,[3 1 2]);
        for i = 1:size(out,3)
            out(:,:,i) = flipud(out(:,:,i));
        end
    elseif strcmp(type,'coord')
        out = [in(:,3) in(:,1) in(:,2)];
    end
    
case 'cor2ax'
    if strcmp(type,'vol')
        for i = 1:size(in,3)
            out(:,:,i) = flipud(in(:,:,i));
        end
        out = permute(in,[2 3 1]);
    elseif strcmp(type,'coord')
        out = [in(:,2) in(:,3) in(:,1)];
    end
end % switch
return
