function [ROI,avg,y] = trialavg2(y,delta,window,varargin)
% [ROI,avg,y] = trialavg2(y,delta,window,varargin)
%
% Required arguments:
% y                 data vector
% delta             cell array of delta functions for each condition
%			(or matrix with 1 col for each condition)
% window            [a b] 2-element vector of window around each trial onset
%                   0 is timepoint of stim onset, negative numbers ok.
%                   for window = [-5 1], 1st 5 tps before the trial, 6th at start,
%                   7th after start - so 7 timepoints overall.
%
% Optional arguments:
%
% 'trialbaseline', 	subtract from baseline of each trial, default 1st timepoint in trial
% 'basepoints',     followed by vector: time points in trial window to use as baseline
%                   relative to trial window - so if you include neg numbers in window,
%                   you're adding timepoints, so adjust basepoints accordingly.
% 'trimoverall', 	remove outliers from overall timeseries
% 'trimtrial', 	    remove outliers on a trial-by-trial basis
% 'trimste',        followed by number: std deviations of trimming threshold
% 'groups',         vector of group assignment number (integers) for each condition
% 'append',         followed by existing structure to append output to
%
% 'plot', 		(1,0) plot or not
% 'nosteplot'		1 suppresses plotting of standard errors
% 'colors'		plot colors
% 'legend'		cell array of legend text
% 'title'		title of figure
%
% 'options',        followed by structure with the above fields
%                   supercedes other optional input arguments if different
%
% Tor Wager, 11/01/01

% -------------------------------------------------------------------
% * set up input arguments
% -------------------------------------------------------------------

doplot = 0;
dotrialbaseline = 0;
dotrimoverall = 0;
dotrimtrial = 0;
basepoints = 1;
trimstd = 3;
groups = [];
nosteplot = 0;

for i = 1:length(varargin)
	if isstr(varargin{i})
		switch varargin{i}
		case 'plot', 		    doplot = 1;
        case 'nosteplot',       nosteplot = 1;
        case 'colors',          O.colors = varargin{i+1};
		case 'trialbaseline', 	dotrialbaseline = 1;
		case 'trimoverall', 	dotrimoverall = 1;
		case 'trimtrial', 	    dotrimtrial = 1;
        case 'basepoints',      basepoints = varargin{i+1};
        case 'trimste',         trimstd = varargin{i+1};
        case 'groups',          groups = varargin{i+1};
        case 'append',          ROI = varargin{i+1};
		end % end switch
	end
end

tp = window(2) - window(1);

% keep a record in O structure appended to output struct
O.plot = doplot;
O.nosteplot = nosteplot;
O.trialbaseline = dotrialbaseline;
O.trimoverall = dotrimoverall;
O.trimtrial = dotrimtrial;
O.basepoints = basepoints;
O.trimstd = trimstd;
O.groups = groups;
O.window = window;

% find option structure in input here.
for i = 1:length(varargin)
	if isstr(varargin{i})
		switch varargin{i}
		case 'options', 		    
            inO = varargin{i+1};
            if isfield(inO,'plot'),doplot = inO.plot; O.plot = inO.plot;,end
            if isfield(inO,'trialbaseline'),dotrialbaseline = inO.trialbaseline; O.trialbaseline = inO.trialbaseline;,end
            if isfield(inO,'trimoverall'),dotrimoverall = inO.trimoverall; O.trimoverall = inO.trimoverall;,end
            if isfield(inO,'trimtrial'),dotrimtrial = inO.trimtrial; O.trimtrial = inO.trimtrial;,end
            if isfield(inO,'basepoints'),basepoints = inO.basepoints; O.basepoints = inO.basepoints;,end
            if isfield(inO,'trimstd'),trimstd = inO.trimstd; O.trimstd = inO.trimstd;,end
            if isfield(inO,'groups'),groups = inO.groups; O.groups = inO.groups;,end
            % plotting only
            if isfield(inO,'nosteplot'),O.nosteplot = inO.nosteplot;,end
            if isfield(inO,'colors'),O.colors = inO.colors;,end
            if isfield(inO,'title'),O.title = inO.title;,end
            if isfield(inO,'legend'),O.legend = inO.legend;,end
        end
    end
end
ROI.options = O;


% -------------------------------------------------------------------
% * loop through conditions and get selective average data
% -------------------------------------------------------------------

if ~iscell(delta)
	delta = mat2cell(delta,size(delta,1),ones(size(delta,2)));
end

for i = 1:length(delta)
    
    % find the values to get from y
    a = find(delta{i} == 1);
    b = a + window(1);
    c = a + window(2);
    d = [b c];
    
    % avoid zero or negative indices 
    d(d(:,1) <= 0 | d(:,2) <= 0,:) = NaN;
    
    % avoid going off end of y
    d(d(:,1) > length(y) | d(:,2) > length(y),:) = NaN;
    
    % count trials skipped
    skipped(i) = sum(isnan(d(:,1)));
    
    % remove skipped trials
    d(isnan(d(:,1)),:) = [];
    
    % do overall trimming here
    if dotrimoverall
        tpts = abs(y) > trimstd .* std(y);
        y(tpts) = NaN;
    end
    
    % get data from y: num trials rows, tp columns
    for j =1:size(d,1)
        tdat = y(d(j,1):d(j,2));
        
        % do trimming and baseline subtract here.
        if dotrimtrial
            tdat(abs(tdat) > trimstd .* std(tdat)) = NaN; 
        end
        
        if dotrialbaseline
            tdat = tdat - tdat(basepoints);
        end
        
        try,e(j,:) = tdat';catch,whos e,tdat,error('Averaging function problem!'),end
    end
    
    numtrials(i) = size(e,1);
    
    avgdata{i} = e;
    avg{i} = nanmean(e);
    ste{i} = getste(e);
    
end

ROI.avg = avg;
ROI.ste = ste;
ROI.avgdata = avgdata;
ROI.numtrials = numtrials;
ROI.skipped = skipped;

% -------------------------------------------------------------------
% * deal with grouping of conditions, if necessary
% -------------------------------------------------------------------
if ~isempty(groups)
    for i = 1:max(groups)
        gdat = avgdata(groups == i);
        grpdata{i} = cell2mat(gdat');
        grpavg{i} = nanmean(grpdata{i});
        grptrials(i) = size(gdat,1);
        grpste{i} = getste(grpdata{i});
    end
    ROI.grpavg = grpavg;
    ROI.grpste = grpste;
    ROI.grpdata = grpdata;
    ROI.grptrials = grptrials;
end

% -------------------------------------------------------------------
% * plotting function call
% -------------------------------------------------------------------
if doplot
    if isempty(groups)
        tor_plot_avgs(avg,ste,O);
    else
        tor_plot_avgs(grpavg,grpste,O);
    end
end

return







% -------------------------------------------------------------------
% * sub-functions
% -------------------------------------------------------------------

function ste = getste(m)

    for tpoint = 1:size(m,2)
		vals = m(:,tpoint);
		if sum(~isnan(vals)) > 1
			ste(tpoint) = std(vals(~isnan(vals)))/ sqrt(length(vals(~isnan(vals))));
		else
			ste(tpoint) = NaN;
		end
    end
        
return
    