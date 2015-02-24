function [str,res,raw] = repeated_anovan(X,Xn,varargin);
% [str,res,raw] = repeated_anovan(X within,covariates between,[contrasts],[IRLS outlier],[plot]);
% factorial repeated measures anova.  If no contrasts, assumes 2 levels per factor 
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
%
% [optional]: contrasts
% Default is to specify factorial 2 x 2 contrasts (only works for 2 factors)
% Raw is returned empty right now.
%
% [optional]: IRLS
% Default is no robust IRLS fitting.
% Enter 2nd varargin argument as flag (1 = do)
%
% [optional]: plot
% Default is 1
%
% Example:
% [str,res]=repeated_anovan(X,cov,[1 1 -1 -1; 1 -1 1 -1; 1 -1 -1 1; 1 1 1 1]);
%
% Enter an n x 4 matrix X of repeated-measures data, subjects are rows
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

if length(varargin) > 1, domcd = varargin{2};,else,domcd = 0;,end
if length(varargin) > 2, doplot = varargin{3};,else,doplot = 1;,end
if domcd, res.flag = 'IRLS';, else, res.flag = 'OLS';,end


% build contrasts within subjects
% -----------------------------------------------------

if length(varargin) > 0, 
    c = varargin{1};,
else

    % number of factors
    % -----------------------------------------------------

    nfact = log2(size(X,2));
    if mod(nfact,1), error('Wrong number of columns for factorial repeated measures!'),end

    for i = nfact:-1:1
        O = ones(1,i);                      % build block
        O = repmat([O -O],1,nfact-i+1);     % should always be size(X,2) columns in O
        c(i,:) = O;                         % string out 
    end
    
    %interactions
    if nfact > 1
        nint = factorial(nfact-1);
        if nfact == 2, c(3,:) = [1 -1 -1 1];,end
    end

end

% effects - Y side
res.c = c;
c = c';
Yi = X * c;

nfact = size(c,2);

for i = 1:nfact
    
    % fit
    
    if isempty(Xn) % no between subjects covariates
        
        if domcd
            [bb,stats]=robustfit(ones(length(Yi(:,i)),1),Yi(:,i),'bisquare',[],'off');   
            b(i) = bb;
            t(i) = stats.t;
            p(i) = stats.p;
            res.dfe = stats.dfe;
        else
            [h,pp,ci,stats] = ttest(Yi(:,i));
            b(i) = mean(Yi(:,i));
            t(i) = stats.tstat;
            p(i) = pp;
            res.dfe = stats.df;
        end
        
    else  % covariates
        
        if domcd
            [bb,stats]=robustfit(Xn,Yi(:,i));   
        else
            [bb,dev,stats]=glmfit(Xn,Yi(:,i));
        end
        
        b(:,i) = bb;
        p(:,i) = stats.p;
        t(:,i) = stats.t;
    
        res.dfe = stats.dfe;
        %res.mse(i) = ssef ./ stats.dfe;
        
        % test plot
        %if i == 3,[r,infostring,sig] = plot_correlation(Xn(:,2),Yi(:,i));,
        %    xlabel('Fact3,cov2')
        %    [r,str,sig,ry,rx,h] = prplot(Yi(:,i),Xn,2);
        %    xlabel('Robust fact3, cov2')
        %    [r,str,sig,ry,rx,h] = prplot(Yi(:,i),Xn,2,1);
        %    xlabel('Non-robust prplot fact3, cov2')
        %end
    end

    raw.stats(i) = stats;

end

res.b = b; res.t = t; res.p = p; 

str = [];strnames=[];
for i = 1:size(res.t,1)
    for j = 1:size(res.t,2)

        if res.p(i,j) < .05, str2='*';,else, str2='';,end
        str = [str sprintf('%3.2f%s\t',t(i,j),str2)];
        
        if i == 1,strnames = [strnames sprintf('F%d\t',j)];,
        else, strnames = [strnames sprintf('F%d x cov%d\t',j,i)];
        end
    end
end

res.str = str;
res.strnames = strnames;

raw.X = X; raw.Yi = Yi; raw.c = c';


return







% old model comparison stuff
% -----------------------------------------------------

[br,devr,statsr]=glmfit(Xn,y);
 
ssef = stats.resid'*stats.resid;
sser = statsr.resid'*statsr.resid;
dfb = statsr.dfe - stats.dfe;
F = ((sser - ssef) ./ dfb) ./ (ssef ./ stats.dfe);
omnip = 1-fcdf(F,dfb,stats.dfe);

str = sprintf('Omni F(%3.0f,%3.0f)=\t%3.2f\t,p = \t%3.4f\t',dfb,stats.dfe,F,omnip);
res.str = str;

if omnip(i) < .05, str2='*';,else, str2='';,end
str = sprintf('%3.2f%s\t',F,str2);

res.b = b; res.t = t; res.p = p; res.F = F; res.omnip = p; res.mse = ssef ./ stats.dfe;
res.dfb = dfb; res.dfe = stats.dfe;

for i = 1:nfact + nint
    if p(i) < .05, str2='*';,else, str2='';,end
    str = [str sprintf('%3.2f%s\t',t(i),str2)];
end

raw.Xd = Xd; raw.Xi = Xi; raw.stats = stats; raw.statsr = statsr;

return

        
        
%str2 = sprintf('Individual parameters\n');
%for i = 1:nfact + nint
%    str2 = sprintf(