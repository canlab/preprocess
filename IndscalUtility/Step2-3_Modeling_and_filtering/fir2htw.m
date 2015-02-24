function [h,t,w,w_times,halfh] = fir2htw(b,varargin)
% [h,t,w,w_times,halfh] = fir2htw(b,[hconstraint],[doplot],[colors cell])
% estimates height, time to peak, and width of FIR response
% 
% tor wager, 2/14/05
% minh = min height

if length(varargin) > 0, 
    hconstraint = varargin{1};
else
    hconstraint = length(b) - 1;
end

if hconstraint > length(b)
    warning('fir2htw is trying to use more betas than there are!  limiting.')
    hconstraint = length(b);
end

if length(varargin) > 1,
    doplot = varargin{2};
else
    doplot = 1;
end

% exit if all vals are same
if all(b == mean(b)), h = NaN;, t = h; w = h; w_times = [NaN NaN];, halfh = NaN;, 
    warning('fir2htw:  all values are the same!')
    return, 
end



colors = {'ro-'};
if length(varargin) > 1, colors = varargin{3};,end 

dat = abs(b(1:hconstraint));            % abs value of b in range
h = max(dat);                           % max abs in first n elements
t = find(dat == h);                     % time this occurs

h = b(t);    % the actual value at max

minh = min(dat);                            % min abs in first n elements
mint = find(dat == minh);                   % time this occurs
minh = b(mint);                             % the actual value at min
% to avoid using max to min, make minh - 0

h = h - minh;       % height = highest - lowest point

% calculations for width
halfh = minh + .5 * h;   %(h - minh);       % 1/2 the max to min distance

d = distance(halfh,b(1:hconstraint));

if h > 0
    wh = find(b(1:hconstraint) > halfh);
else
    wh = find(b(1:hconstraint) < halfh);
end

if isempty(wh), wh = NaN;,end

% first half
x2 = max(wh(1),1);    % first above halfh
x1 = max(wh(1)-1,1);  % first above-half - 1
y2 = b(x2); 
y1 = b(x1); 
m = y2 - y1;          % slope, x reduces to 1 (x2 - x1 = 1)
if m == 0
    w_times(1) = x1;    % exact match
else
    w_times(1) = x1 + (halfh - y1) ./ m; % solve y = mx + b for x* given m; y1 = b
end

% second half
x1 = min(wh(end),hconstraint);      % last one above halfh
x2 = min(wh(end)+1,hconstraint);    % + 1
y2 = b(x2);  y1 = b(x1);
m = y2 - y1;                      % slope, x reduces to 1 (x2 - x1 = 1)
if m == 0
    w_times(2) = x1;    % exact match
else
    w_times(2) = x1 + (halfh - y1) ./ m;     % solve y = mx + b for x* given m; y1 = b
end


w = w_times(2) - w_times(1);        % width in elements


if doplot
    
    try
        
    hold on
    h1 = arrow([t minh],[t h+minh],'Length',12);        % height arrow
    h2 = arrow([0 h+minh],[t h+minh],'Length',12);      % delay arrow
    h3 = arrow([w_times(1) halfh],[w_times(2) halfh],'Length',12);
    h4 = arrow([w_times(2) halfh],[w_times(1) halfh],'Length',12);
    
    h5 = text(t + .1 * w, halfh + .4*halfh,'h','FontSize',16,'FontWeight','bold');
    h6 = text(t - .5 * w, h - .2*halfh,'t','FontSize',16,'FontWeight','bold');
    h7 = text(t - .5 * w, halfh - .2*halfh,'w','FontSize',16,'FontWeight','bold');

    set(h1,'Color',colors{1}(1))
    set(h2,'Color',colors{1}(1))
    set(h3,'Color',colors{1}(1))
    set(h4,'Color',colors{1}(1))
    set(h5,'Color',colors{1}(1))
    set(h6,'Color',colors{1}(1))
    set(h7,'Color',colors{1}(1))
    
    catch
        warning('Error drawing arrows!')
    end
end



return


% OLD EXTRA STUFF

range = max(b(1:hconstraint)) - min(b(1:hconstraint));

% find all values within some tolerance of the half-max height
% if not enough points found to determine width, then interpolate until we
% get values.
% it would be better to interpolate!
tolval = .02; wh = []; resampleval = 1;

tol = tolval * range;   % tolerance of tolval % of range
wh = find(d <= tol);    % find all points within tolerance

while length(find(wh>t)) < 1 | length(find(wh<t)) < 1
    resampleval = resampleval + 1;
    
    b2 = resample(b(1:hconstraint),resampleval,1); % interpolate
    
    d = distance(halfh,b2);
    wh = find(d <= tol);    % find all points within tolerance
    wh = wh ./ resampleval;  % convert back to original time units
    
    if resampleval > 10, warning('Could not find width!'), break, end
end
    
% now find nearest elements to peak that are at max height
d_from_t = wh - t;
w_times = [max(d_from_t(d_from_t < 0)) min(d_from_t(d_from_t > 0))];



if length(w_times) > 1
    w = w_times(2) - w_times(1);        % width in elements
    w_times = t + w_times;              % time indices of elements before and after peak at 1/2 max height
else
    w_times = [NaN NaN];
    w = NaN;                            % width is undefined; 1st or last tp was max
end


