function hout = barplot_columns(dat,varargin)
% axishandle = barplot_columns(dat,[title],[covs])
%
% This function makes a barplot of columns of data, with standard error
% bars.  Optional arguments include removing continuous covariates before plotting,
% robust (IRLS) estimation of means and correlations with covariates, and
% within-subject error bars based on the subject x condition interaction
% (overall), which is not quite the standard error contrasts of interest,
% but is the standard error for a 1-way repeated measures ANOVA.
%
% plots circles around points at z >= 1.96
% plots individual points, unless you enter 4th argument
%
% Examples:  Just plot means and SE
% h = barplot_columns(tmp,'Cluster 1',[],1);
%
% Optional arguments
% 1 - Title for figure
% 2 - covariates
% 3 - String Arguments
%       'nofig' : do not make figure
%       'noind' : do not plot individual scores
%       'plotout': circle potential outliers at z>1.96 in red
%       'dorob' : do robust IRLS means and correlations
%       'groupby' : followed by number to make bars in groups of
%       'legend' : followed by legend string, for grouping
% Examples:
% barplot_columns(ctmp,'RT effects by Switch Type',overall_sw,'nofig','dorob')
% barplot_columns(j2avgso,'J2 Means',[],'noind','groupby',4,'legend',{'No Switch' 'Attribute Switch' 'Object Switch' 'Double Switch'}); saveas(gcf,'j2_means','tiff')
%
% Standard Errors ARE NOT Adjusted for covariate, right now.

% ----------------------------------------------------  
% > Set up input arguments
% ----------------------------------------------------  

dofig = 0; doind = 1; plotout = 0; dorob = 0; groupby = 0; legstr = [];
if length(varargin) > 2,
    for i = 3:length(varargin)
        if strcmp(varargin{i},'nofig'), dofig = 0;, end
        if strcmp(varargin{i},'noind'), doind = 0;, end
        if strcmp(varargin{i},'plotout'), plotout = 1;, end
        if strcmp(varargin{i},'dorob'), dorob = 1;, end
        if strcmp(varargin{i},'groupby'), groupby = varargin{i+1};, end
        if strcmp(varargin{i},'legend'), legstr = varargin{i+1};, end
    end
end

        
if length(varargin) > 1, 
    covs = varargin{2};, if ~isempty(covs), covs = scale(covs,1); ,end, 
else, covs = [];, end

% delete nans casewise
%wh = find(any(isnan(dat),2));
%dat(wh,:) = [];
%if ~isempty(covs), covs(wh,:) = [];,end

% replace nans with mean
for i = 1:size(dat,2),
    if any(isnan(dat(:,i)))
        warning('Some NaNs!')
        %dat(find(isnan(dat(:,i))),i) = nanmean(dat(:,i));
    end
end


dat1 = dat;     % save original data for ind subj plot




% ----------------------------------------------------  
% > Get means and standard error of means
% With robust option, if specified, and removing
% covariates, if there are any.
% ---------------------------------------------------- 


if dorob
    % ROBUST FIT
    stderr = [];
    
    for i = 1:size(dat,2)
        
        
        if ~isempty(covs)
            
            tmpx=zscore(covs); tmpy=dat(:,i);
            wh = find(isnan(tmpx) | isnan(tmpy));
            tmpx(wh,:) = []; tmpy(wh) = [];
        
            n(i) = length(tmpy);
        
            [bb,stats]=robustfit(tmpx,tmpy);
            stderr(i) = stats.se(1);
            [b3,stats]=robustfit(tmpx,zscore(tmpy)); 
            robcor(i) = b3(2);
            
        else
            % just get N and tmpy
            tmpy=dat1(:,i);
            wh = find(isnan(tmpy));
            tmpy(wh) = [];
            n(i) = length(tmpy);
            [bb,stats]=robustfit(ones(size(tmpy,1),1),tmpy,'bisquare',[],'off'); 

        end
        dattmp(i) = bb(1);
    end
    if ~isempty(covs), mycor = robcor;,end
    if isempty(stderr), stderr = ste(dat1);,end
else
    % OLS FIT
    dat = nanmean_t(dat);
    stderr = ste(dat1);
    if ~isempty(covs)
        
        for i = 1:size(dat,2)
        
            tmpx=zscore_t(covs); tmpy=dat1(:,i);
            wh = find(isnan(tmpx) | isnan(tmpy));
            tmpx(wh,:) = []; tmpy(wh) = [];
             n(i) = length(tmpy);
             
            tmp = corrcoef([tmpx tmpy]);
            olscor(i) = tmp(1,2);
            X = [ones(size(tmpx,1),1) tmpx];
            r = tmpy - X * pinv(X) * tmpy;
            stderr(i) = ste(r);
            
        end

        mycor = olscor;
    else
        % just get N
        for i = 1:size(dat,2)
            tmpy=dat1(:,i);
            wh = find(isnan(tmpy));
            tmpy(wh) = [];
            n(i) = length(tmpy);
        end
    end
end

if dorob,

    dat = dattmp;
    %disp(['Robust Means']), disp(dat);
    %if ~isempty(covs), disp(['Robust Correlations (1st cov)']), disp(robcor);, end
else
    %disp(['OLS Means']), disp(dat);
    %if ~isempty(covs), disp(['Correlations (1st cov)']), disp(olscor);, end
end
%disp(['Residual Standard Error']), disp(stderr);


% ----------------------------------------------------  
% > Make figure
% ---------------------------------------------------- 

if dofig
    f = figure('Color','w'); hout = gca; set(gca,'FontSize',12); %hold on; grid on;
else
    f = get(gcf); hout = gca; set(gca,'FontSize',12); hold on;
end
   
got_print_matrix=0;
if got_print_matrix;
% ----------------------------------------------------  
% > Report descriptives and stats
% ----------------------------------------------------
fprintf(1,'Means:\t'), print_matrix(dat)
fprintf(1,'STE:\t'), print_matrix(stderr)
t = dat./stderr; 
for i =1:length(t), p(i) = 2.* (1 - tcdf(abs(t(i)),n(i)-1));,end

fprintf(1,'t-scores:\t'), print_matrix(t)
fprintf(1,'p (2 tailed):\t'), print_matrix(p)
fprintf(1,'\n')

if exist('mycor') == 1
    fprintf(1,'Correl:\t'), print_matrix(mycor)

    for i =1:length(t),
        [rci,sig,z,p(i)] = r2z(abs(mycor(i)),n(i),.05);
    end

    fprintf(1,'correl p:\t'), print_matrix(p)
end
end
% ----------------------------------------------------  
% > BARPLOT
% ---------------------------------------------------- 
if ~groupby
    % standard barplot, no grouping
        h = bar(dat); set(h,'FaceColor',[.8 .8 .8])
        tor_bar_steplot(dat,stderr,{'k'});
        
        set(gca,'XLim',[0 size(dat1,2)+1],'XTick',1:size(dat1,2),'XTickLabel',1:size(dat1,2))
        xlabel('Task Condition'), ylabel('BOLD contrast')
        
        %set(gcf,'Position',[464   283   930   827]), drawnow
       
        %s = std(dat1);
  
else
    % do grouping
    
    set(gca,'FontSize',16)
    switch groupby
    case 4
        offset = - 1 + .55; mult = .185;
        
    case 3
        offset = - 1 + .585; mult = .21;
    otherwise
        offset = - 1 + .65; mult = .23;
    end
        
    cm = [1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 0;0 1 0; 1 1 0; 0 0 1];
    cm = colormap(gray);

    %xm = mean(x); 
    %sterr = ste(x);     % .* tinv(.975,size(x,1));  % 95% CI

    xm2 = reshape(dat,[groupby length(dat)./groupby])';
    sterr = reshape(stderr,[groupby length(stderr)./groupby])';

	hh = bar(xm2); colormap(cm)

    for i = 1:size(sterr,1)
        tor_bar_steplot([xm2(i,:)],[sterr(i,:)],{'k'},i + offset,mult)
    end
    
    if ~isempty(legstr), legend(hh,legstr);,end
    
end % grouping
    
