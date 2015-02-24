function cl = cluster_barplot_plotonly(cl,varargin)
% cl = cluster_barplot_plotonly(cl,['colors'],[colors],['covs'],[covs])
%
% This function complements cluster_barplot.  It only plots data, but it
% plots data next to a visualization of the cluster.  
% A tiff file of each figure is saved.
%
% Optional arguments
% (Empty arguments are ignored)
% 2 - covariates
% 3 - String Arguments
%       'dorob' : do robust IRLS means and correlations
%       'xlab' : followed by x-axis labels
%       'colors' : followed by colors in cell array: {'r' 'g' 'b'}
%       'covs' : followed by btwn-subjects covariates in columns
%       'contrasts' : followed by contrast values
%                   if entered, runs repeated measures anova
%                   repeated_anovan
%                   saves results in cl(i).BARPLOT.anovares
%       'transform' : followed by function to transform each element of
%                   data by before testing, e.g., 'x.^2' or '1/x'
%       'peakdata'  : use individual peak data within ROI rather than
%                   average of ROI.  cl(i).BARPLOT.peakdata.
%
% Examples:
%
% The following will plot and test 3 contrasts using repeated measures
% anova, with robust fitting for the estimates.  Data are the negative of
% data stored in cl.BARPLOT.data.
% contrasts = [1 1 -1 -1; 1 -1 1 -1; 1 -1 -1 1];
% covs = EXPT.cov;
%xlab = {'Pcontrol' 'Pplacebo' 'Wcontrol' 'Wplacebo'};
%funct = '-1 * x';
% cluster_barplot_plotonly(cl(1),'transform',funct,'xlab',xlab, ...
% 'covs',covs,'contrasts',contrasts,'dorob');


% ----------------------------------------------------  
% > Set up input arguments
% ---------------------------------------------------- 
warning off MATLAB:divideByZero    
dorob = 0; xlab = []; colors = {'r'}; covs = []; contrasts = [];f = [];peakdata = 0;

if length(varargin) > 0,
    for i = 1:length(varargin)
        if strcmp(varargin{i},'dorob'), dorob = 1;, end
        if strcmp(varargin{i},'xlab'), xlab = varargin{i+1};, end
        if strcmp(varargin{i},'colors'), colors = varargin{i+1};, end
        if strcmp(varargin{i},'covs'), covs = varargin{i+1};, end
        if strcmp(varargin{i},'contrasts'), contrasts = varargin{i+1};, end
        if strcmp(varargin{i},'transform'), f = inline(varargin{i+1});, end
        if strcmp(varargin{i},'peakdata'), peakdata = 1;, end
    end
end
  
if isempty(colors), colors = {'r'};,end   
if ~iscell(colors), c = colors; colors = [];, colors{1} = c;, clear c, end

diary cluster_barplot_output_report.txt
dorob,contrasts,covs,colors,f,peakdata
for i = 1:length(nargin), inputname(i),end
    
for i = 1:length(cl)
    
    % ----------------------------------------------------  
    % > Define data, colors, and names
    % ---------------------------------------------------- 
    if peakdata,
        dat = cat(2,cl(i).BARPLOT.peakdata{:});
    else,
        dat = cl(i).BARPLOT.data;
    end
    
    if ~isempty(f)          % transform data if specified
        dat = f(dat);
    end
    
    % names of regions
    
   AddOpts.Resize='on';
   AddOpts.WindowStyle='normal';
   AddOpts.Interpreter='tex';
    if ~isfield(cl(i),'shorttitle'), cluster_orthviews(cl(i),{[1 0 0]});
        answer=inputdlg('Enter the short title for this cluster.','I need your help!',1,{['Cluster ' num2str(i)]},AddOpts);
        cl(i).shorttitle =  answer{1};
    end
    if isempty(cl(i).shorttitle), 
        cluster_orthviews(cl(i),{[1 0 0]});
        answer=inputdlg('Enter the short title for this cluster.','I need your help!',1,{['Cluster ' num2str(i)]},AddOpts);
    
        cl(i).shorttitle = answer{1};
        %cl(i).shorttitle = ['Cluster ' num2str(i)];,
    end
    
    figure('Color','w'); 
    set(gcf,'Position',[1513         519        1064         501]);
    
    if length(colors) < i, colors = [colors colors];,end
    
    % ----------------------------------------------------  
    % > Anatomy plot
    % ---------------------------------------------------- 

    montage_clusters_maxslice([],cl(i),colors(i));    
    tmp = get(gca,'Position');
    tmp(3) = tmp(3) .* .5; tmp(1) = tmp(1) .* .5;
    set(gca,'Position',tmp);
    title(cl(i).shorttitle,'FontSize',24);
    
    % ----------------------------------------------------  
    % > Data plot
    % ---------------------------------------------------- 
    mm = sprintf('Centered at (%3.0f, %3.0f, %3.0f) mm',cl(i).mm_center(1),cl(i).mm_center(2),cl(i).mm_center(3));
    h = axes('position',[.5 tmp(2) .4 tmp(4)]);
    warning off;, barplot_columns(dat,mm,covs,'noind','nofig','dorob');, warning on
    if ~isempty(xlab), set(gca,'XTickLabel',xlab);, end
    
    % tmp for opioid
    ylabel('Opioid activity (3 - Bmax/Kd)')
    
    % ----------------------------------------------------  
    % > Repeated measures anova
    % ----------------------------------------------------  
    if ~isempty(contrasts)
        warning off, [str,res]=repeated_anovan(dat,covs,contrasts,dorob);, warning on
        res.contrasts = contrasts;
        cl(i).BARPLOT.anovares = res;
        if i ==1, fprintf(1,'Name\t%s\t\n',cl(i).BARPLOT.anovares.strnames),end
        fprintf(1,'%s\t%s\t\n',cl(i).shorttitle,cl(i).BARPLOT.anovares.str)
    end
    
    saveas(gcf,['cl' num2str(i) '_slice_bar'],'tif')
    
    % tmp for opioid
    tmp = res.p(1:2,2:end);
    if any(tmp(:) < .05), % leave open 
            saveas(gcf,['cl' num2str(i) '_slice_bar'],'fig')
    elseif any(res.p(1:2,1) < .05), 
            % if pain effect is sig
            set(gcf,'Position',[65         544        1064         501])
    else
        close,
    end
        
    
end

diary off

return


