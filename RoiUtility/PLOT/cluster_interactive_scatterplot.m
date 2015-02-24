function [xdat,ydat] = cluster_interactive_scatterplot(cl, varargin)
% [xdat,ydat] = cluster_interactive_scatterplot(cl, [image names], [X data]);
%
% Set up interactive scatterplots in spm_orthviews
%
% This function plots a behavioral variable (for example), called [X data]
% above, against average brain activity measures within each 'blob' of a
% clusters structure.
% It assumes you have one number per subject stored in the cl.timeseries
% field.  These do not have to be time series values; they are actually
% individual difference values in this case.  Yes, this is confusing
% nomenclature--sorry.
%
% If you enter image names, this function will extract clusters from them.
% In this case, you can enter [] for the cl argument.
% If you enter [] as a second argument, the function will assume you've
% entered clusters in cl struct format as the first argument.
% If you enter a 3rd argument, it will assume that's a variable with one
% element per subject that will be plotted against brain data.
%
%
% need to:
% extract data:
% cl = tor_extract_rois(EXPT.SNPM.P{1},cl);
% specify covariate or seed:
%
% [xdat,ydat] = cluster_interactive_scatterplot(cl, [image names], [X data]);
% 
% Example: In robust regression directory:
% load SETUP  % contains image names, covariates
% cluster_interactive_scatterplot([clpos clneg], SETUP.files, SETUP.covariates)


if nargin > 1 && ~isempty(varargin{1})
    cl = tor_extract_rois(varargin{1},cl);
end

if nargin > 2
    for i = 1:length(cl), cl(i).xdat = varargin{2}; end
end

% default
xdat = []; ydat = [];

% check for orthviews fig and create orthviews if necessary
fh = findobj('Tag','Graphics'); 
if ~isempty(fh),
    figure(fh);
else
    disp('Cannot find Graphics window. Creating.');
    cluster_orthviews(cl,{[1 0 0]});
end

% create new figure if necessary
h = findobj('Tag','scatterplot_axis');
if length(h) > 1, delete(h); end
if isempty(h) || ~ishandle(h)
    fh = findobj('Tag','Graphics'); figure(fh);
    
    % shell
    shellh = axes('Position',[.46 .01 .46 .42]);
    set(shellh,'YColor',[1 1 1],'XColor',[1 1 1],'YTickLabel',[],'XTickLabel',[]);
    
    % axes
    h = axes('Position',[.54 .08 .34 .34]);
    set(h,'Tag','scatterplot_axis','FontSize',16);
    
    % set interactive function to run this
    set(fh,'WindowButtonUpFcn', @cluster_interactive_callback)

    %set(gcf,'WindowButtonUpFcn','[xdat,ydat] = cluster_interactive_scatterplot(cl);');
end

% --------------------------------------------
% attach data to figure
% --------------------------------------------
fh = findobj('Tag','Graphics'); figure(fh);
data = guidata(fh);
data.cl = cl;
%if ~isfield(data,'clusterfield') || isempty(data.clusterfield)
    data.clusterfield = 'timeseries';
%end
%if ~isfield(data,'xlab'),
    data.xlab = spm_input(['Enter x-axis label: '],[],'s','Behavioral score');
%end
%if ~isfield(data,'ylab'),
    data.ylab = spm_input(['Enter y-axis label: '],[],'s','Contrast data');
%end
data.axish = h;         % axis handle
guidata(fh,data);                   % store data 

return

    






function cluster_interactive_callback(varargin)

% get data stored in figure
% should contain: clusterfield, xlab, ylab
% if doesn't exist, create as default: 'timeseries'
fh = findobj('Tag','Graphics'); figure(fh);
data = guidata(fh);

N = fieldnames(data);
for i=1:length(N)
    eval([N{i} ' = data.' N{i} ';']);
end

% check for data in specified field of cl
if ~isfield(cl,clusterfield), display_error('nodata',clusterfield); end
if ~isfield(cl,'xdat'), display_error('noxdata',clusterfield); end

spm_orthviews_showposition;

% activate scatterplot axis
axes(axish);

% find closest cluster and return index
% ---------------------------------------------------
pos = spm_orthviews('Pos')';


% check to see if we're in a cluster
wh = 0; 
centers = cat(1,cl.mm_center);

% find closest cluster, based on center
d = distance(pos,centers); wh = find(d == min(d)); wh = wh(1);

% only accept if cursor is w/i 2 mm of a voxel in the cluster
d = distance(pos,cl(wh).XYZmm'); d = min(d); if d > 15, wh = 0; end



if wh
    %cluster_table(cl(wh));
    fprintf(1,'Cl. %3.0f, Voxels: %3.0f, Coords: %3.0f, %3.0f, %3.0f\n',wh,cl(wh).numVox,cl(wh).mm_center(1), cl(wh).mm_center(2),cl(wh).mm_center(3));

    %spm_orthviews('Reposition',cl(wh).mm_center);
    
    % get data from structure
    ydat = cl(wh).(clusterfield);
    xdat = cl(wh).xdat;
    
    % check for data
    if isempty(ydat),display_error('nodata',clusterfield); end
    if isempty(xdat), display_error('noxdata',clusterfield), end
        
    cla;
  
    % get design matrix, with 1st column as regressor of interest
    dorobust = 1;
    if isfield(cl, 'nuisance') && ~isempty(cl(wh).nuisance)
        X = [xdat cl(wh).nuisance];
    else
        X = xdat;
    end
    if no_intercept(X), X(:,end+1) = 1; end
    
    %[tmp,tmp,r] = partialcor(X,ydat,1,1,dorobust);
    
    plot_correlation(X,ydat,'robust','xlabel',data.xlab,'ylabel',data.ylab);
    
%     if isfield(cl, 'nuisance') && ~isempty(cl(wh).nuisance)
%         [r,str,sig,ry,rx,h,rr] = prplot(ydat,[xdat cl(wh).nuisance],1,1,{'ko'});
%     else
%         plot_correlation_samefig(xdat,ydat,[],'ko',0,1);
%     end
    %axis auto;
%     xlabel(data.xlab);
%     ylabel(data.ylab);

drawnow

% Return data to workspace
disp('Assigning braindata with cluster average in workspace.');
assignin('base', 'braindata', ydat);

else
    cla;
    str = sprintf('No nearby cluster. Distance is %3.0f',d);
    fprintf(1,str); pause(.5); fprintf(repmat('\b',1,length(str)));
end

return



function display_error(meth,clusterfield)

switch meth
    case 'nodata'
        disp(['You have asked to get data from cl.' clusterfield ', but there is no data there.']);
        error('For a fix, try cl = tor_extract_rois(my_image_name_matrix_here,cl);');
        
    case 'noxdat'
        error('You must assign the x-axis data to plot to cl(*).xdat);');
end

return


function val = no_intercept(X)

val = all(any(diff(X)));

return

