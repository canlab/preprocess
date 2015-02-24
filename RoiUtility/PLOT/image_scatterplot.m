function image_scatterplot(image1, image2,varargin)
% image_scatterplot(image1, image2, varargin)
% 
% Inputs:
% Two image files (must be in the same space!) in .img or .nii format
% image1 and image2 can also be image_vector objects
%
% Optional inputs:
% Any of the words below, followed by a new value:
% nbins = 20;
% plotpoints = 0;
% plotlines = 1;
% plotcontour = 1;
% pointcolor = [.5 .5 .5];
% contourcolormap = 'jet';
%
% Examples:
% create_figure('weightcorr'); image_scatterplot(hr_intercept, gsr_intercept, 'plotlines', 0);
% create_figure('weightcorr'); image_scatterplot(hr_intercept, gsr_intercept, 'plotlines', 0, 'nbins', 40);
%
% tor wager, august 2010
% 
 
nbins = 20;
plotpoints = 0;
plotlines = 1;
plotcontour = 1;

pointcolor = [.5 .5 .5];
contourcolormap = 'jet';

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case 'nbins', nbins = varargin{i+1};
            case 'plotpoints', plotpoints = varargin{i+1};
            case 'plotlines', plotlines = varargin{i+1};
            case 'plotcontour', plotcontour = varargin{i+1};
                
            case 'pointcolor', pointcolor = varargin{i+1}; varargin{i+1} = [];
            case 'contourcolormap', contourcolormap = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning('varargs:UnknownOption', 'Unknown input string option: %s\n', varargin{i});
        end
    end
end

% does not work for comma-separated volumes: 
% if ~exist(image1, 'file'), error('%s cannot be found!', image1); end
% if ~exist(image2, 'file'), error('%s cannot be found!', image2); end

if isa(image1, 'image_vector')
    dat = image1.dat;
else
    [tmp, dat] = iimg_read_img(image1, 2);
end

if isa(image2, 'image_vector')
    dat2 = image1.dat;
else
    [tmp, dat2] = iimg_read_img(image2, 2);
end

% get rid of exactly 0 or NaN values
whomit = any([dat == 0 | isnan(dat) dat2 == 0 | isnan(dat2)], 2);
dat(whomit) = [];
dat2(whomit) = [];

edges1 = linspace(min(dat)-eps, max(dat)+eps, nbins);
[tmp, bin] = histc(dat, edges1);

edges2 = linspace(min(dat2)-eps, max(dat2)+eps, nbins);
[tmp, bin2] = histc(dat2, edges2);

[means, stds] = deal(zeros(nbins, 2));

for i = 1:nbins
    wh1 = bin == i;
    %wh2 = bin2 == i;
    
%     % for contour plot
%     bincenters(i, 1) = nanmean(dat(wh1));
%     bincenters(i, 2) = nanmean(dat2(wh2));
    
    % average and std of image2 for each decile of image 1
    means(i, 1) = nanmean(dat(wh1));
    means(i, 2) = nanmean(dat2(wh1));
    stds(i, 1) = nanstd(dat(wh1));
    stds(i, 2) = nanstd(dat2(wh1));
    
end

% for contour plot
for i = 1:nbins
    for j = 1:nbins
        wh1 = bin == i;
        wh2 = bin2 == j;
     
        % stuff for color-mapped 2D histogram
        surfarray(j, i) = sum(wh1 & wh2);
        
    end
end


create_figure('hist2d');

% Points
if plotpoints
    pointh = plot(dat,dat2,'k.', 'MarkerSize', 1, 'Color', pointcolor);
end
drawnow

% Lines
if plotlines
    plot(means(:, 1), means(:, 2), 'ks-', 'MarkerSize', 10, 'MarkerFaceColor', [.0 .0 .0]);
    
    plot(means(:, 1), means(:, 2) + stds(:, 2), 'k--', 'MarkerSize', 6, 'MarkerFaceColor', [.0 .0 .0]);
    plot(means(:, 1), means(:, 2) - stds(:, 2), 'k--', 'MarkerSize', 6, 'MarkerFaceColor', [.0 .0 .0]);
    
    %lineh = errorbar_horizontal(means(:, 1), means(:, 2), stds(:, 1));
    %linehx = errorbar(means(:, 1), means(:, 2), stds(:, 2));
    %set([lineh, linehx], 'Color', 'k');
end
drawnow

% Contour
if plotcontour
    hold on;
    
    % log scale
    surfarray = log(surfarray);
    surfarray(isinf(surfarray)) = -1;
    
    [X, Y] = meshgrid(edges1, edges2);
    %hh = surf(X, Y, surfarray);
    [cc, hh] = contourf(X, Y, surfarray, 20);
    %set(hh, 'FaceAlpha', .5);
    set(get(hh, 'children'), 'FaceAlpha', .5, 'EdgeColor', 'none')
    %view(0, 90)
    
    % adjust colormap to white BG
    colormap(contourcolormap);
    cm = colormap(gca);
    cm(1:1, :) = ones(1, 3);
    colormap(cm);
    
    legend(hh, 'Log frequency')
    axis tight
    drawnow
end

end % function

