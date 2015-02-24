function  c = plot_cluster_fir(c,cl,i,varargin)
%  c = plot_cluster_fir(c,cl,i,[do HP filter],[do height time wid plot])
%
% plots deconvolution data given a contrast structure:
% c = 
%        numframes: [30 30 30 30]
%            model: [664x121 double]
%             data: [664x14x2 double]   % time x regions x subjects
%         contrast: [1 1 1 1]
%              ons: {1x24 cell}
%            delta: {[664x1 double]  [664x1 double]  [664x1 double]  [664x1 double]}
%            names: {'CA'  'UA'  'CP'  'UP'}
%               TR: 2
%               bf: 30
%            nsess: 4
%       model_desc: 'DX deconv FIR model'
%          px_desc: 'pinv of c.model'
%               px: [121x664 double]
%              con: [1x121 double]
%         con_desc: 'expanded contrast with all basis functions for each trial type'
%    contrast_desc: 'contrast across trial types (collapsing basis functions)'
%
% cl is clusters structure; if data is not in c.data, tries to get it from
% cl.all_data
%
% i is index of which cluster to plot
% optional last argument is HP FILTER ON (1) or FILTER OFF (2)
%
% tor wager
%
%
% Example:
% % make DX FIR deconv model from SPM.mat file
% --------------------------------------------------
% load SPM
% next line for OLD SPM ONLY!!
% N = whos; for i = 1:length(N), eval(['SPM.' N(i).name ' = 'N(i).name]);,end
% c = spm2dx(SPM);
% 
% put data into c structure
% --------------------------------------------------
% c.data = [];
%for i =1:length(cl)
%    c.data(:,i,:) = cl(i).all_data;
%end
%
% define contrast and HP filter
% c.contrast = [1 1 1 1 1];
% c.HP = 100;
%
% run
% plot_cluster_fir(c,cl,1);     % first cluster
%
% to collapse across conditions, i.e., antic and pain
% c.contrast = [1 1 0 0; 0 0 1 1];
% [c.model,c.con] = apply_dx_contrast(c.model,c.numframes,c.contrast);
% c.px = pinv(c.model);
% c.names = {'Anticipation' 'Antic + pain'};



% --------------------------------------------------
% filtering and plotting options
% --------------------------------------------------

dofilt = 1; dohtw = 1;
if length(varargin) > 0, dofilt = varargin{1};,end
if length(varargin) > 1, dohtw = varargin{2};,end

    
colors = {'ro-' 'go--' 'bd-' 'yd--' 'ks-' 'cs--'};
if isfield(c,'colors'), colors = c.colors;,end

if ~isfield(c,'HP'), c.HP = 120;, end
if ~isfield(c,'filtertype'), c.filtertype = 'spm';,end
if ~isfield(c,'doHP'),c.doHP = 1;,end
if ~isfield(c,'firstimg'),c.firstimg = 0;,end

% --------------------------------------------------
% other required fields
% --------------------------------------------------

if ~isfield(c,'px') % pseudoinverse of model
    if ~isfield(c,'model'), error('You must include a model in c.  see, e.g., spm2dx'),end
    c.px = pinv(c.model);
end

if ~isfield(c,'numframes') % number of regressors for each trial type
    error('You must include c.numframes, the number of regressors for each trial type'),
end

if ~isfield(c,'TR') % number of regressors for each trial type
    error('You must include c.TR, the scan repetition time'),
end

if ~isfield(c,'con') % pseudoinverse of model
    if ~isfield(c,'contrast'), error('You must include c.contrast or c.con (expanded); see expand_contrast.m'),end
    c.con = expand_contrast(c.contrast,c.numframes,c.model);
end

if ~isfield(c,'names') % number of regressors for each trial type
    c.names = cell(1,length(c.numframes));
end

if iscell(c.model)
    if size(c.data,1) ~= size(c.model{1},1), error('Data and model have diff numbers of timepoints!'),end
else
    if size(c.data,1) ~= size(c.model,1), error('Data and model have diff numbers of timepoints!'),end
end

% --------------------------------------------------
% get data and save in y
% --------------------------------------------------

% get data in right format, if necessary
%cl(cat(1,cl.numVox)==1) = [];        % optional -- get rid of one-voxel clusters
if ~isfield(c,'data')
    if size(cl(1).all_data,2) == 1   % only one subject
        c.data = cat(2,cl.all_data);
        c.data = c.data(1:size(c.model,1),:);     % limit data to size of model
    else
        % multi-subject
        c.data = [];
        for i =1:length(cl)
            c.data(:,i,:) = cl(i).all_data;
        end
    end
end

if length(size(c.data)) > 2
    % 3-D data, 3rd dim is subjects
    y = squeeze(c.data(:,i,:)); % now subjects are columns

else
    y = c.data(:,i);   
end


% --------------------------------------------------
% data filtering
% --------------------------------------------------

    c.y = y;
    if dofilt, y = filterAdjust(c);, end
    [dat,V] = smooth_timeseries(y,.01);   % small amount of smoothing
    
    if iscell(c.px),
        % DIFFerent models for each subject
        for subj = 1:length(c.px), b(:,subj) = c.px{subj} * dat(:,subj);,end
    else
        % SAME model for all subjects
        b = c.px * dat;       % fit model to this timeseries; cols of b are subjects
    end
    
    if length(size(c.data)) > 2  % 3-D data, 3rd dim is subjects
        sterr = ste(b')';          % ste across subjects
        b = mean(b,2);             % mean beta across subjects
    end
    
    % --------------------------------------------------
    % make figure and plot cluster on slice
    % --------------------------------------------------
    
    tor_fig(1,2);
    subplot(1,2,1); montage_clusters_maxslice([],cl(i),{'r'});
    subplot(1,2,2);

    % color bars to signify conditions
    hold on

    % --------------------------------------------------
    % plot lines where they should go
    % height/width/delay plots, if specified    
    % --------------------------------------------------
    [lineh,legstr] = lineplots(c.con,b,c,colors,c.names,dohtw);

    
    % --------------------------------------------------
    % standard error plots, if requested
    % --------------------------------------------------
    if exist('sterr') == 1
        steplots(c.con,b,c,colors,sterr);
    end
    
    % --------------------------------------------------
    % box below plot showing period(s)
    % --------------------------------------------------
    periodbox(c.con,b,c);
    
        
    
    legend(lineh,legstr);
    xlabel('Time (TRs)')
    ylabel('BOLD activation')
    drawnow;
    
    
    return
    
    
    
    
    % --------------------------------------------------
    % plot lines where they should go
    % --------------------------------------------------    
    function [lineh,legstr] = lineplots(con,b,c,colors,names,varargin)
    
    lineh = []; % handles for plotted lines
    legstr = {};    % legend string for names
    
    xstart = 1; 
    for cind = 1:size(con,1)
        % for each type (contrast)
        b2 = b(find(con(cind,:) > 0));    % take only betas specified by this contrast

        start = 1;          
        xtext = .5; % start for xtext
        for j = find(c.contrast(cind,:) > 0)    % plot this beta series
  
            x = xstart:xstart+c.numframes(j)-1;   % set x coord range for this type
            
            dat = b2(start:start+c.numframes(j)-1);
            lineh(end+1) = plot(x,dat,colors{j},'LineWidth',3);, 
        
            % height time to peak delay, if requested
            if length(varargin) > 0, if varargin{1} > 0, fir2htw(dat,length(dat)-1,1,colors(j)), end,end
            
            start = start + c.numframes(j);
            
            legname = names{j}; legname(find(legname=='_')) = ' ';
            legstr{end+1} = legname;
            
        end
            
        xstart = xstart + c.numframes(j);   % update x values
    end
    
    return
    
    
    % --------------------------------------------------
    % standard error plots, if requested
    % --------------------------------------------------    
    function steplots(con,b,c,colors,sterr)
        % b = betas
        % c = context structure; uses numframes
        % colors = cell array of colors
        % sterr = standard error values
        
    xstart = 1; 
    for cind = 1:size(con,1)
        % for each type (contrast)
        b2 = b(find(con(cind,:) > 0));    % take only betas specified by this contrast
        sterr2 = sterr(find(con(cind,:) > 0));  
        
        start = 1;          
        
        for j = find(c.contrast(cind,:) > 0)
  
            x = xstart:xstart+c.numframes(j)-1;   % set x coord range for this type

            h =fill_around_line(b2(start:start+c.numframes(j)-1), ...
                sterr2(start:start+c.numframes(j)-1) ...
                ,colors{j}(1),x);
        
            start = start + c.numframes(j);
         
        end
        xstart = xstart + c.numframes(j);   % update x values
    end
    
    return
    
    
    % --------------------------------------------------
    % box below plot showing period(s)
    % --------------------------------------------------
    function periodbox(con,b,c,j)
        
        
    b = b(1:end-1);
    m = min(b); rg = .1 *(max(b) - m); % min and 10% of max height for plotting location
    xstart = 1;    
    for cind = 1:size(con,1)
        
        j = find(c.contrast(cind,:) > 0); j = j(1);
        color = [0 0 0] + 1./cind;
        
        h = fill([xstart xstart+c.numframes(j)-.5 xstart+c.numframes(j)-.5 xstart], [m-1.5*rg m-1.5*rg m-.5*rg m-.5*rg],color);
            set(h,'FaceAlpha',.6);
            %text(xstart+xtext,m-rg,names{j},'Color','k','FontSize',16,'FontWeight','b')
            %start = start + c.numframes(j);
            %xtext = xtext + 5;  % shift text over
            
        xstart = xstart + c.numframes(j);   % update x values
    end

    return
    
