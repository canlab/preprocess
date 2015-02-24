function [h,ph,fh] = plotmask4(mask,varargin)
% function [h,ph,fh] = plotmask4(mask,mycolor,axis_handles,fig_handle,scalefactor,whichplot,ori)
% optional args
%     h - vector of axis handles, for existing plot
%     mycolor - must be vector of 3 numbers defining mycolor.  default = white.
%		OR, mycolor can be 1, 2, or 3, indicating graded mycolor
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
% THIS VERSION WORKS WITH IMAGING THE TRANSPOSE OF SLICES, SO THE BRAIN IS VERTICAL,
% AND WILL ROTATE YOUR POINTS 90 DEGREES!! DO NOT USE IMAGESC(SLICE) AND THIS PLOT.
% USE IMAGESC(SLICE')
%
% Tor Wager, 1/12/03

% defaults
% - - - - - - - - - - - - - - - - - - - - - - - - - -
whichplot = 1:size(mask,3);
scalex = 1;
scaley = 1;
mycolor = [1 1 1];

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
        mycolor = varargin{1};
	% set up mycolor gradient
	if ~ischar(mycolor) & length(mycolor) == 1
	
		% if it's an integer, then use
		% graded mycolor!
		% ---------------------------------------------------------------
		% * build mycolor map, for graded mycolor on masks with integer values
		% * mycolor input must be integer from 1-3, indicating r g or b
		% ---------------------------------------------------------------
		numcols = 60; %max(max(max(mask)));
		cols = zeros(numcols,3);
		% cols(:,mycolor) = .3 + .3 * ((1:numcols)' ./ numcols);
		%cols(:,mycolor) = (.5 + .2 * (1:numcols)') ./ max(.5 + .2 * (1:numcols)');
        cols(:,mycolor) = (.5:.5/(60-1):1)';
		mycolor = cols;

	else
		% ---------------------------------------------------------------
		% * build mycolor map, of all-same-mycolor equal to input mycolor
		% ---------------------------------------------------------------
		mycolor = mycolor .* ones(max(1,max(max(max(mask)))),3);
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
	
	fh = figure; set(gcf,'color','w')
	for i = 1:size(mask,3)
		h(i) = subplot(ceil(sqrt(size(mask,3))),ceil(sqrt(size(mask,3))),i);
	end
end

[histo,bins] = hist(mask(:),60);

% plot mask point by point
% - - - - - - - - - - - - - - - - - - - - - - - - - -
for j = 1:size(mask,3) 
      	subplot(h(whichplot(j)))
        hold on
	 
	% could do this, but need different axis handles!
	% [y,x,z] = ind2sub(size(mask),mask > 0);
	
	[y x] = find(mask(:,:,j) > 0);
    
        for k = 1:length(y)
               if ~ischar(mycolor) & size(mycolor,1) > 1
                   whbin= find((bins - mask(y(k),x(k),j)).^2 == min((bins - mask(y(k),x(k),j)).^2));
                    pcol = mycolor(whbin(1),:);
                else
                    pcol = mycolor;
                end
	      %pcol = mycolor(mask(y(k),x(k),j),:);
          % HERE THIS IS REVERSED, FOR IMAGESC(SLICE')
              ph = plot(y(k)*scaley,x(k)*scalex,'s','color',pcol,'MarkerSize',2,'MarkerFaceColor',pcol);
	end

end % loop of slices for this image



if exist('cols') == 1 & size(cols,1) > 0
	% ---------------------------------------------------------------
	% * build mycolor bar
	% ---------------------------------------------------------------
	%figure; 
	%set(gcf,'Position',[256   569   402   123])
	%bar(eye(size(cols,1)),1,'stacked')
	%colormap(cols)
end

return