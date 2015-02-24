function [res,raw] = repeated_anovan(Y,Xn,varargin)
% [res,raw] = repeated_anovan(Y,Xbtwn,keywords and additional inputs);
%
% factorial repeated measures anova.  If no contrasts, assumes 2 levels per
% factor
%
% Tor Wager
%
% t, p, and betas are reported from fastest to slowest factors
% or, if specific contrasts are entered, in the order of entry
%
% Tested: Same as SPSS 11.01 repeated measures for 2 x 2 within / 1 between mixed
%         effects model, 6 / 2 / 03.  (without robust MCD option)
%
%         Same as SPSS 11.01 repeated measures for 2 x 2 x 2 within, 11/7/03
%         Same as SPSS 11.04 general linear model, 8/06
%
% [optional]: 
% enter matrix of within-subjects contrasts
% Default is to specify factorial 2 x 2 contrasts (only works for 2 factors)
% Raw is returned empty right now.
%
% Keywords are:
% 'notable', suppress table output
% 'noplot', suppress plots
% 'robust', robust IRLS fitting
%
% 'btwnnames', followed by names of btwn subjects covariates (cell array)
% 'winames', followed by names of btwn subjects covariates (cell array)
% 'colors', followed by names of line colors (cell array)
%
% Examples:
% [str,res]=repeated_anovan(Y,cov,[1 1 -1 -1; 1 -1 1 -1; 1 -1 -1 1; 1 1 1 1]);
%
% [str,res] = repeated_anovan(Y,EXPT.cov,[1 1 -1 -1; 1 -1 1 -1; 1 -1 -1 1], ...
% 'robust','btwnnames',{'Placebo' 'Order'},'winames',{'Heat' 'P-C' 'P-C*H-W'},'ynames',{'HC' 'HA' 'LC' 'LA'});
%
% Enter an n x 4 matrix Y of repeated-measures data, subjects are rows
% cov is a n x k matrix of k between subjects covariates
% contrasts entered as row vectors
%
% res.t reports t-values: Columns are contrasts, rows are covariate effects
% The first row is the effect of the contrast (vs 0)
% Additional rows are contrast x covariate interactions
% res.p and res.b are ordered the same way.
% res.b has effects (regression weights), which are some scaling of the mean
% contrast value, depending on how contrast weights are scaled.
% for 1 -1 contrasts, they should be half the mean difference.
%
% Automatically centers covariates in the robust case.
%


[c,dorob,dotable,doplot,btwnnames,winames,ynames,res] = setup_inputs(varargin);


% build contrasts within subjects if needed
% -----------------------------------------------------
if isempty(c), c = get_contrasts(Y); end


% -------------------------------------------------------------------------
% get contrast data, including 'contrast' of all ones for btwn Ss effects
% Y-side
% -------------------------------------------------------------------------

c = c';         % c is transformed to column vectors

[c, wh_btwn, nwithin] = prepare_btwn_contrast(c);
nbtwn = size(Xn,2);     % without intercept

Yi = Y * c;     % contrast data
res.c = c;     % save as col vectors
res.wh_btwn = wh_btwn;

% save for posterity, and also needed for partial correlations in table
res.X = Xn;
res.Y = Y;
res.Yi = Yi;

[btwnnames, winames, ynames, res] = setup_names(res,btwnnames,winames,ynames, Y, nbtwn,nwithin);




% -------------------------------------------------------------------------
% fit model of btwn-ss effects for each w/i subjects effect
% -------------------------------------------------------------------------

for i = 1:nwithin + 1

    % fit

    if isempty(Xn) % no between subjects covariates

        if dorob
            [bb,stats]=robustfit(ones(length(Yi(:,i)),1),Yi(:,i),'bisquare',[],'off');
            res.dfe = stats.dfe;
        else
            [h,pp,ci,stats] = ttest(Yi(:,i));
            bb = mean(Yi(:,i));
            stats.t = stats.tstat;
            stats.p = pp;
            res.dfe = stats.df;
        end

        res.b(i) = bb;
        res.t(i) = stats.t;
        res.p(i) = stats.p;

    else  % covariates

        if dorob
            [bb,stats]=robustfit(Xn,Yi(:,i));
        else
            [bb,dev,stats]=glmfit(Xn,Yi(:,i));
        end


        res.dfe = stats.dfe;

        res.b(:,i) = bb;
        res.p(:,i) = stats.p;
        res.t(:,i) = stats.t;

        %res.mse(i) = ssef ./ stats.dfe;

    end

    if isfield(stats,'w'), res.w(:,i) = stats.w; else res.w(:,i) = ones(size(Y,1),1); end
    raw.stats(i) = stats;

end



% get strings for results and names (stretched out format
str = [];strnames=[];
for i = 1:size(res.t,1)
    for j = 1:size(res.t,2)

        % t-value
        if res.p(i,j) < .05, str2='*'; else, str2=''; end
        str = [str sprintf('%3.2f%s\t',res.t(i,j),str2)];

        % name
        if i == 1,strnames = [strnames winames{j} sprintf('\t')];
        else, strnames = [strnames winames{j} '*' btwnnames{i} sprintf('\t') ];
        end
    end
end

res.str = str;
res.strnames = strnames;

raw.Y = Y; raw.Yi = Yi; raw.c = c';

% -----------------------------------------------------
% table
% -----------------------------------------------------

if dotable, res = regtable(res); end

% -----------------------------------------------------
% plot
% -----------------------------------------------------

if doplot
    wh_col = 1;
    res.colors = regplot2(Xn,Y,wh_col,dorob,mean(res.w,2),btwnnames,ynames,res.colors); 
end

end



% -----------------------------------------------------
% -----------------------------------------------------
% -----------------------------------------------------


% Sub-functions


% -----------------------------------------------------
% -----------------------------------------------------
% -----------------------------------------------------




% -----------------------------------------------------

% Setup: inputs

% -----------------------------------------------------

function [c,dorob,dotable,doplot,btwnnames,winames,ynames,res] = setup_inputs(inputs)

dotable = 1;
doplot = 1;
dorob = 0;
c = [];
btwnnames = [];
winames = [];
ynames = [];
res.colors = [];

for i = 1:length(inputs)
    if ismatrix(inputs{i})
        c = inputs{i};
    elseif isstr(inputs{i})
        switch inputs{i}
            % reserved keywords
            case 'notable', dotable = 0;
            case 'noplot', doplot = 0;
            case 'robust', dorob = 1;

            case 'btwnnames', btwnnames = inputs{i+1};
            case 'winames', winames = inputs{i+1};
            case 'ynames', ynames = inputs{i+1};
            case 'colors', res.colors = inputs{i+1};
            otherwise, warning(['Unknown input string option:' inputs{i}]);
        end
    end
end

if dorob, res.flag = 'IRLS';  else, res.flag = 'OLS'; end

end


% -------------------------------------------------------------------------

% SETUP: get names

% -------------------------------------------------------------------------
function [btwnnames, winames, ynames, res] = setup_names(res,btwnnames,winames,ynames, Y, nbtwn,nwithin);

if isempty(btwnnames)
    btwnnames = {'Intercept'};
    for i = 1:nbtwn, btwnnames{i+1} = ['Cov' num2str(i)]; end
else
    btwnnames = [{'Intercept'} btwnnames];
end

if isempty(winames)
    winames = {'Btwn'};    % intercept is first; this is for btwn-subjects effects
    for i = 1:nwithin, winames{i+1} = ['WIcon' num2str(i)]; end
else
    winames = [{'Btwn'} winames];
end

if isempty(ynames)
    for i = 1:size(Y,2), ynames{i} = ['y' num2str(i)]; end
end

res.btwnnames = btwnnames;
res.winames = winames;
res.ynames = ynames;

if length(btwnnames) ~= nbtwn + 1, error('Between names are wrong length.'), end
if length(winames) ~= nwithin + 1, error('Within names are wrong length.'), end

end



function c = get_contrasts(Y);

% number of factors
% -----------------------------------------------------

nfact = log2(size(Y,2));
if mod(nfact,1), error('Wrong number of columns for factorial repeated measures!'),end

for i = nfact:-1:1
    O = ones(1,i);                      % build block
    O = repmat([O -O],1,nfact-i+1);     % should always be size(Y,2) columns in O
    c(i,:) = O;                         % string out
end

%interactions
if nfact > 1
    nint = factorial(nfact-1);
    if nfact == 2, c(3,:) = [1 -1 -1 1]; end
end

end

% -----------------------------------------------------

function [c, wh_btwn, nwithin] = prepare_btwn_contrast(c)
% c is column vectors

% prepare between-subjects effects
wh_btwn = find(all(diff(c) == 0));   % a contrast of all ones, for mean effect
if ~isempty(wh_btwn)                % get rid of it if it exists
    c(:,wh_btwn) = [];
end

nwithin = size(c,2);  % number of within/ss effects, excluding intercept
ndat = size(c,1);     % no data cols
c = [(1 ./ ndat) * ones(ndat,1) c] ;  % add to front
wh_btwn = 1;

end


% -----------------------------------------------------

% OUPUT: Regression table

% -----------------------------------------------------

function res = regtable(res)

nwithin  = size(res.c,2);   % c has one btwn, rest within
nbtwn = size(res.b,1);         % including intercept, which is 1st

% get partial correlations
for i=2:nbtwn
    [xx,yy,partr(i,1)] = partialcor(res.X,res.Yi(:,1),i-1,0,strcmp(res.flag,'IRLS'));
    
    for j=2:nwithin
        [xx,yy,partr(i,j)] = partialcor(res.X,res.Yi(:,j),i-1,0,strcmp(res.flag,'IRLS'));
    end
end
res.r_partial = partr;



fprintf(1,'\nMixed General Linear Model    \n')
fprintf(1,'Method is %s, dfe = %3.0f    \n',res.flag,res.dfe)
fprintf(1,'\nBetween-subjects effects    \n')
fprintf(1,'------------------------    \n')
fprintf(1,'Parameter\tb-hat\tt\tp\tr_partial    \n')
for i=1:nbtwn
    fprintf(1,'%s\t%3.2f\t%3.2f\t%3.4f\t%3.2f\t    \n', ...
        res.btwnnames{i},res.b(i,res.wh_btwn),res.t(i,res.wh_btwn),res.p(i,res.wh_btwn),partr(i,1));
end

fprintf(1,'\nWithin-subjects effects    \n')
fprintf(1,'------------------------    \n')
fprintf(1,'Parameter\tb-hat\tt\tp\t    \n')

for i=2:nwithin
    fprintf(1,'%s\t%3.2f\t%3.2f\t%3.4f\t    \n', ...
        res.winames{i},res.b(1,i),res.t(1,i),res.p(1,i));
end

fprintf(1,'\nBetween * Within     \n')
fprintf(1,'------------------------    \n')
for i=2:nbtwn
    for j=2:nwithin
        myname = [res.btwnnames{i} '*' res.winames{j}];

        fprintf(1,'%s\t%3.2f\t%3.2f\t%3.4f\t%3.2f\t    \n', ...
            myname,res.b(i,j),res.t(i,j),res.p(i,j),res.r_partial(i,j));
    end
end
fprintf(1,'\n')

end



% -------------------------------------------------------------------------

% OUTPUT: regression plot

% -------------------------------------------------------------------------

function colors = regplot2(X,Y,wh_col,dorob,w,btwnnames,ynames,colors)

r = size(Y,2);
if dorob, robstr = 'robust'; else robstr = []; end

tor_fig; set(gca,'FontSize',24)

if isempty(colors), colors = {'ro' 'gs' 'cv' 'mh'}; end
while length(colors) < r, colors = [colors colors]; end

for i = 1:r
    h = plot_correlation(X,Y(:,i),'col',wh_col,'colors',colors(i),robstr,'noprint','weights',w);
    han(i) = h{1}(1);
end

han2 = findobj(gcf,'Type','text');
delete(han2)

% range over which to get fits for line
xrange = max(X(:,1)) - min(X(:,1));
minx = min(X(:,1)) - .1*xrange;
maxx = max(X(:,1)) + .1*xrange;

set(gca,'XLim',[minx maxx]);

% Y range
yrange = max(Y(:)) - min(Y(:));
miny = min(Y(:)) - .1*yrange;
maxy = max(Y(:)) + .1*yrange;
set(gca,'YLim',[miny maxy]);

xlabel(btwnnames{wh_col+1});
ylabel('Adjusted data');
%ynames = {'HC' 'HA' 'LC' 'LA'};
if ~isempty(ynames), legend(han,ynames), end

% change yellow to orange
h = findobj(gca,'Color','y');
set(h,'Color',[1 .5 0]);

scn_export_papersetup(400);


end







% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % old model comparison stuff
% % -----------------------------------------------------
% 
% [br,devr,statsr]=glmfit(Xn,y);
% 
% ssef = stats.resid'*stats.resid;
% sser = statsr.resid'*statsr.resid;
% dfb = statsr.dfe - stats.dfe;
% F = ((sser - ssef) ./ dfb) ./ (ssef ./ stats.dfe);
% omnip = 1-fcdf(F,dfb,stats.dfe);
% 
% str = sprintf('Omni F(%3.0f,%3.0f)=\t%3.2f\t,p = \t%3.4f\t',dfb,stats.dfe,F,omnip);
% res.str = str;
% 
% if omnip(i) < .05, str2='*'; else, str2=''; end
% str = sprintf('%3.2f%s\t',F,str2);
% 
% res.b = b; res.t = t; res.p = p; res.F = F; res.omnip = p; res.mse = ssef ./ stats.dfe;
% res.dfb = dfb; res.dfe = stats.dfe;
% 
% for i = 1:nwithin + nint
%     if p(i) < .05, str2='*'; else, str2=''; end
%     str = [str sprintf('%3.2f%s\t',t(i),str2)];
% end
% 
% raw.Xd = Xd; raw.Xi = Xi; raw.stats = stats; raw.statsr = statsr;
% 
% end
% 

% ISMATRIX: Returns 1 if the input matrix is 2+ dimensional, 0 if it is a scalar 
%           or vector.
%
%     Usage ismat = ismatrix(X)
%
% RE Strauss, 5/19/00

function ismat = ismatrix(X)
  [r,c] = size(X);
  if (r>1 && c>1)
    ismat = 1;
  else
    ismat = 0;
  end

end



%str2 = sprintf('Individual parameters\n');
%for i = 1:nfact + nint
%    str2 = sprintf(