function plotmask(mask,varargin)
% function plotmask(mask,h,color,scalefactor,whichplot,ori)
% optional args
%     h - vector of axis handles, for existing plot
%     color - default = 'ws', white square
%     scalefactor [x y] - index on plot = mask dims * scale[x or y]
%     whichplot - vector of subplot numbers to plot each slice of mask onto
%     Other orientations don't work so well without whichplot, bec. there
%       are, for saggital, 256 slices it tries to find plot handles for.

if nargin > 5
    switch varargin{5}
    case 'sagg', mask = permute(mask,[3 2 1]);
    case 'cor', mask = permute(mask,[3 1 2]);
    end
end    
 
whichplot = round(1:28/size(mask,3):28);
scalex = 1;
scaley = 1;
color = 'ws';
   
if nargin > 4,
    if ~isempty(varargin{4})
        whichplot = varargin{4};
    end
end

if nargin > 3
    if ~isempty(varargin{3})
        scalef = varargin{3};
        scalex = scalef(1);
        scaley = scalef(2);
    end
end

if nargin > 2
    if ~isempty(varargin{2}) 
        color = varargin{2};
    end
end

if nargin > 1
   plotH = varargin{1};
else
   figure;
   for i = 1:28
      plotH(i) = subplot(7,4,i);
   end
end

% T-MAPS .....................................................
% plot mask point by point
      for j = 1:size(mask,3) % All slices of tmap mapped onto closest anatomical
         subplot(plotH(whichplot(j)))
         hold on
         % slice = whichplot(j);
         for x = 1:size(mask,2)
            for y = 1:size(mask,1)
               if mask(y,x,j) == 1,plot(x*scalex,y*scaley,color,'MarkerSize',3,'MarkerFaceColor',color(1)),end
            end
         end
      drawnow
      end % loop of slices for this image
