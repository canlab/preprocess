function [h,ph,fh] = plotmask(mask,varargin)
% function [h,ph,fh] = plotmask(mask,color,axis_handles,fig_handle,scalefactor,whichplot,ori)
% optional args
%     h - vector of axis handles, for existing plot
%     color - must be vector of 3 numbers defining color.  default = white.
%     scalefactor [x y] - index on plot = mask dims * scale[x or y]
%     whichplot - vector of subplot numbers to plot each slice of mask onto
%     Other orientations don't work so well without whichplot, bec. there
%       are, for saggital, 256 slices it tries to find plot handles for.
%
% outputs
%	h	vector of axis handles for slice plots
%	ph	handle of last point plotted; useful for legends
%
% Tor Wager, 1/15/02

% defaults
% - - - - - - - - - - - - - - - - - - - - - - - - - -
whichplot = 1:size(mask,3);
scalex = 1;
scaley = 1;
color = [1 1 1];


% set up input arguments
% - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin > 6
    switch varargin{6}
    case 'sagg', mask = permute(mask,[3 2 1]);
    case 'cor', mask = permute(mask,[3 1 2]);
    end
end    
 

if nargin > 5,
    if ~isempty(varargin{5})
        whichplot = varargin{5};
    end
end

if nargin > 4
    if ~isempty(varargin{4})
        scalef = varargin{4};
        scalex = scalef(1);
        scaley = scalef(2);
    end
end

if nargin > 3
    if ~isempty(varargin{3}) 
        fh = varargin{3};
	figure(fh)
    end
end

if nargin > 2
    if ~isempty(varargin{2}) 
        h = varargin{2};
    end
else
   h = [];
end

if nargin > 1
    if ~isempty(varargin{1}) 
        color = varargin{1};
    end
end



% make the figure if necessary
% - - - - - - - - - - - - - - - - - - - - - - - - - -
if isempty(h)
	
	fh = figure; set(gcf,'Color','w')
	for i = 1:size(mask,3)
		h(i) = subplot(ceil(sqrt(size(mask,3))),ceil(sqrt(size(mask,3))),i);
	end
end



% plot mask point by point
% - - - - - - - - - - - - - - - - - - - - - - - - - -
for j = 1:size(mask,3) 
      	subplot(h(whichplot(j)))
        hold on
	 
	[y x] = find(mask(:,:,j) > 0);

        for k = 1:length(y)
              ph = plot(x(k)*scalex,y(k)*scaley,'s','Color',color,'MarkerSize',2,'MarkerFaceColor',color);
	end

end % loop of slices for this image

return
