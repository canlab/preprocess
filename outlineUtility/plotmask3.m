function [h,ph,fh] = plotmask(mask,varargin)
% function [h,ph,fh] = plotmask(mask,color,axis_handles,fig_handle,scalefactor,whichplot,ori)
% optional args
%     h - vector of axis handles, for existing plot
%     color - must be vector of 3 numbers defining color.  default = white.
%		OR, color can be 1, 2, or 3, indicating graded color
%		of either red, green, or blue, depending on the integer value in the mask!
%		Returns an error for non-integer masks.
%
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
	% set up color gradient
	if ~ischar(color) & length(color) == 1
	
		% if it's an integer, then use
		% graded color!
		% ---------------------------------------------------------------
		% * build color map, for graded color on masks with integer values
		% * color input must be integer from 1-3, indicating r g or b
		% ---------------------------------------------------------------
		numcols = max(max(max(mask)));
		cols = zeros(numcols,3);
		% cols(:,color) = .3 + .3 * ((1:numcols)' ./ numcols);
		cols(:,color) = (.5 + .2 * (1:numcols)') ./ max(.5 + .2 * (1:numcols)');
		color = cols;

	else
		% ---------------------------------------------------------------
		% * build color map, of all-same-color equal to input color
		% ---------------------------------------------------------------
		color = color .* ones(max(1,max(max(max(mask)))),3);
	end

    end
end

if ~(any(any(any(mask)))), 
	disp('plotmask3: empty mask.'),
	ph = plot(0,0,'ko','MarkerSize',0.01);, 
	if ~(exist('fh') == 1), fh = [];,end,
	return, 
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
	 
	% could do this, but need different axis handles!
	% [y,x,z] = ind2sub(size(mask),mask > 0);
	
	[y x] = find(mask(:,:,j) > 0);

        for k = 1:length(y)
	      pcol = color(mask(y(k),x(k),j),:);
              ph = plot(x(k)*scalex,y(k)*scaley,'s','Color',pcol,'MarkerSize',2,'MarkerFaceColor',pcol);
	end

end % loop of slices for this image



if exist('cols') == 1 & size(cols,1) > 0
	% ---------------------------------------------------------------
	% * build color bar
	% ---------------------------------------------------------------
	figure; 
	set(gcf,'Position',[256   569   402   123])
	bar(eye(size(cols,1)),1,'stacked')
	colormap(cols)
end

return