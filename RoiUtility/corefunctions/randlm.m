function out = randlm(Type,y,sub,cond,noint,lincheck,varargin)
% function out = randlm(y,sub,cond,'plot' [opt])
%
% Type  :   'simult' or 'sequential', to fit predictors sim. or one at a time
% y     :   unadjusted data from voi
% sub   :   the struct output of getscans3, with eachevt option
% cond  :   the condition number you want to test, 1 or 2 for visual l/r
% noint :   a mtx of regressors of no interest, without intercept, to be removed prior to fitting
% lincheck: a linear model that you want to check against the randfx fit, without intercept
% plot  :   any additional arguments will call the plot commands for this function.
%
% if you DON'T want to adjust for predictors of no interest or fit a linear model, use [] for these args.

if nargin > 6, plot = 1;,else plot = 0;,end

conds = [];
times = [];
alltsince = [];
allnumin = [];
lagToCount = 4000; 		% Lag in seconds that you want to count prior events in

% Build list of time since last and num in 4/8/whatever
for i = 1:8
    eval(['scan = sub.s' num2str(i) ';'])
    eval(['conds = [conds;sub.s' num2str(i) '(:,2)];'])
    eval(['times = [times;sub.s' num2str(i) '(:,2)];'])
    
    timess = scan(scan(:,2) == cond,1);
    tsincel(1:2,1) = 99999;         % remove first two trials
    numin8(1:2,1) = 0;
    for i = 3:size(timess,1)
        tsincel(i,1) = timess(i) - timess(i-1);
        numin8(i,1) = sum(timess < timess(i) & timess > timess(i) - lagToCount,1);
    end
    alltsince = [alltsince;tsincel];
    allnumin = [allnumin;numin8];
    clear scan, clear timess, clear tsincel, clear numin8
end
alltsince = [99999;alltsince]; allnumin = [0;allnumin];      % add extra space for intercept
conds = conds';

% Build the X matrix
out.Xorig = sub.X(:,conds == cond); out.Xorig = [ones(size(out.Xorig,1),1) out.Xorig];  % add intercept
out.Xfilt = filterX(out.Xorig);
out.ctimes = times(conds == cond,:);
out.conds = conds';

% Filter and adjust X and y
[out.nointX,y] = adjusty(noint,y);		% returns filtered and adjusted y, and the X regressed on it to adjust
out.y = y;

% ===== regress linear model, if spec. =======
if ~isempty(lincheck)
    lincheckX = filterX(lincheck); 
    out.lincheckX = [ones(size(lincheckX,1),1) lincheckX];
    out.linbeta = (out.lincheckX \ y)';
    out.linfits = out.lincheckX * out.linbeta';
    out.linseb = sqrt(diag(inv(out.lincheckX'*out.lincheckX)) * var(y - out.linfits))';
    out.lint = out.linbeta ./ out.linseb;
    out.linp = tdist(out.lint,size(out.linfits,1)-size(out.lincheckX,2));
    lin.fits = out.linfits;
    lin.y = y;
    lin.e = y - out.linfits;
end  

disp('Condition number of X matrix is ' num2str(cond(out.Xfilt)))

% ===== rand fx fit =========================
switch Type
case 'simult'
    out.betas = (out.Xfilt \ y);
    out.fits = out.Xfilt*out.betas;
    out.e = out.y - out.fits;
case 'sequential'
    intercpt = ones(size(out.Xfilt,1),1);                       % fit intercept
    out.betas(1) = intercpt \ y;
    y = y - intercpt * out.betas(1);
    for i = 2:size(out.Xfilt,2)                                 % fit each event and adjust
        m = [intercpt out.Xfilt(:,i)]; b = m \ y;
        out.betas(i,1) = b(2);
        y = y - m*b;
    end
    out.fits = out.Xfilt*out.betas;
    out.e = out.y - out.fits;
case 'ridge'
	k = .1;
	out.betas = ridge(out.y,out.Xfilt,k);
	out.fits = out.Xfilt*out.betas;
    out.e = out.y - out.fits;
		% determine bias
	xtx = out.Xfilt' * out.Xfilt;
	out.bias = k^2 * out.betas' * (xtx + k*eye(size(xtx,1)))^-2 * out.betas;
case 'ridgechoosek'
end

out.meanb = mean(out.betas);
out.seb = sqrt(var(out.betas) / size(out.betas,1));
out.t = out.meanb / out.seb;
out.p = tdist(out.t,size(out.betas,1)-1);
out.data = [out.betas/max(out.betas) alltsince allnumin];   % scale betas between -1 and 1
out.data(alltsince == 99999,:) = [];                        % remove first two trials
out.Xorig = sparse(out.Xorig);
out.Xfilt = sparse(out.Xfilt);

% ======== fit times to betas ===========
y = out.data(:,1);
X = [ones(size(out.data,1)) out.data(:,2) out.data(:,2).^2 out.data(:,3) out.data(:,3).^2];
betas = X \ y; betas = betas / betas(1);                    % scale betas as change from intercept
out.nlbetas = betas;


if plot
    out = plotrandfx(out);
    voistat('plotall',lin);
    voistat('plotall',out);
    drawnow
end

return






function X = filterX(X);
    S = load('S'); S = S.S;
    S.KL = full(S.KL);S.KH = full(S.KH);
    X = X - S.KH*(S.KH'*X); 
return

function y = filtery(y);
    S = load('S'); S = S.S;
    S.KL = full(S.KL);S.KH = full(S.KH);
    y = y - S.KH*(S.KH'*y);  
    y = y*100 / mean(y);                    % scale to be percent change
return

function [nointX,y] = adjusty(noint,y);
	if ~isempty(noint)
			% Fit predictors and temporal derivatives of regressors to be modeled out
			% Input: a convolved set of regressors, minus the intercept (noint), timeseries (y)
    	
		nointX = filterX(noint); 							% smooth and filter the X matrix
 		[dummy,tempder] = gradient(nointX);			 
		nointX = [ones(size(nointX,1),1) nointX tempder];	
			% Return the residuals from the regression
		y = filtery(y);										% smooth and filter y
    	y = y - nointX * (nointX \ y);                  	% adjust y, removing effects of no interest
    	nointX = sparse(nointX);
	end

function out = plotrandfx(out)
    figure; clear M;
    subplot(3,2,5); plot(out.data(:,2),out.data(:,1),'bo');grid on; title('time since last')
    subplot(3,2,6); plot(out.data(:,3),out.data(:,1),'bo');grid on; title('num in last 8')
    subplot(3,2,1);hist(out.data(:,1),25); title('betas')
    subplot(3,2,2); hist(out.data(:,2),25); title('time in ms since last one')
    subplot(3,2,3); hist(out.data(:,3),25); title('number in last 8 seconds')
    
    subplot(3,2,4); plot3(out.data(:,2),out.data(:,3),out.data(:,1),'bo'); grid on; title('z = betas')
    for i = 1:4:720
        view(i,30)
        M(i) = getframe(gca);
    end
    %subplot(3,2,4);movie(M);
    %out.film = M;
return