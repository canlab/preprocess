function [out1,out2] = getv(funct,ts,varargin)
% function [out1,out2] = getv(funct,ts,varargin)
%
% funct:
%	'get' 		: get cross-correlations up to n/4, given a timeseries (ts).
%				[xc,adjustedts] = getV('get',ts,cov);
%	'make'		: construct an autocorrelation matrix (V) from the cross-correlations.
%				[V,xc] = getv('make',xc,dim of V matrix);
%	'get&make'	: both of the above.
%				[V,xc] = getv('get&make',ts,cov);
%
%	'test'		: compute large-lag standard error and t-tests on how many parameters
%					in the xc function are not 0.
%				stats = getv('test',ts,xc);
%
% xc: the autocorrelation function, with xc(1) = 1
%
% ts:
%	a timeseries vector with one column.
%	the ts output is the ajusted ts after removal of covariates and their derivatives.
%
% cov:
%	optional covariates to regress out of timeseries before computing autocorrelation.
%	-	should be column vectors.
%	-	fits covariates you entered and their temporal derivative.
%
% output arguments:
% 'get'
%       out1 = r;                 % 1st output = autocorr. function
%       out2.ts = ts;             % timeseries used, after adjustments
%       out2.se = ser;            % large-lag st. error
%       out2.t = tr;              % t-values
%       out2.p = pr;              % p-values
%
% 'make'
%       out1 = V                  % autocorrelation matrix
%       out2 = xc                 % autocorrelation function used, after tapering.
%
% last edit 4/1/01 by Tor Wager to ensure positive definite

switch funct
case 'get'
	if nargin > 2, 
		[covX,ts] = adjusty(varargin{1},ts);	
	end
	N = size(ts,1);
	
	% This uses Matlab's xcorr function, but it overest. the autocorrelation for random numbers
	%[xc,lags] = xcorr(ts,floor(N/4),'coeff');
	%xc = xc(lags >=0,:);
	%plot(xc)
	%V = xc;
	
	% this is the formula from Box and Jenkins, 1970, p. 32
	tsmean = mean(ts);
	for k = 1:floor(N/4)+1
		clear zt
		clear ztk
		clear prod
		for j = 1:N - (k-1)
			zt(j) = ts(j) - tsmean;
			ztk(j) = ts(j+k-1) - tsmean;
			prod(j) = zt(j) .* ztk(j);
		end
		sumprod(k) = sum(prod);
		c(k) = (1/N) * sumprod(k);
		r(k) = c(k) / c(1);
	end
	out1 = r;                 % 1st output = autocorr. function
    out2 = ts;             % timeseries used, after adjustments
	
	
case 'test'
	r = varargin{1};	% autocorr. function = r
	N = size(ts,1);
    % compute large-lag standard error, after Box and Jenkins, 1970, p. 35.
    ser(1) = sqrt( 1/N );
    for i = 2:size(r,2)
        ser(i) = sqrt( (1 + 2 * sum(r(1:i).^2)) / N );
    end
    tr = r ./ ser;            % t-values
    df = N - 1;          % one-sample t-test, repeated for each lag
    pr = tdist(tr',df)';      % p-values
    
	out1.xc = r;          	  % 1st output = autocorr. function
    out1.se = ser;            % large-lag st. error
    out1.t = tr;              % t-values
    out1.p = pr;              % p-values
    
    
	
case 'make'
	if nargin > 2, dim = varargin{1};,else dim = size(ts,2);,end
	if size(ts,2) <= 1,error('xc should be a row vector.'),end
	xc = ts;
	
	% figure;plot(xc)
	
	% taper data; see Chatfield, 'Analysis of time series', p. 148
	% this is not applied, as it is in the book, to the actual timeseries. ?
	sizexc = size(xc,2);
	num = round(sizexc * .05); 			% taper about 5% of values.
	txc = xc;
	for i = sizexc - num : sizexc
		txc(i) = (sizexc - i)*xc(i) / (num+1);
	end
	
	txc(sizexc+1:dim) = 0;
	V = toeplitz(txc);
	
	%figure;plot(xc)
	%hold on; plot(txc,'r')
	
	% Old method without using the toeplitz function, on original xc.
	%myones = ones(dim,1);
	%V = diag(myones);
	%for i = 2:size(xc,2)
	%	d = myones((i:end),1) .* xc(i);
	%	V = V + diag(d,i-1) + diag(d,-(i-1));
	%end
	
	% increase taper if matrix is still not positive definite
   [dummy,test] = chol(V);  % also pos definite of all eig(V) are positive.  
   if test(1) > 0, warning('	V matrix is not positive definite!'),  %should test value be 0?  
		disp('		Increasing taper to stabilize variance.')
		while test(1) > 0
			num = num + 1;
			txc = xc;									% new tapered xc
			for i = sizexc - num : sizexc				% taper
				txc(i) = (sizexc - i)*xc(i) / (num+1);
			end
			txc(sizexc+1:dim) = 0;						% fill in with zeros
			V = toeplitz(txc);							% make V
			[dummy,test] = chol(V); 					% test V
		end
		disp(['		Tapered ' num2str(num) ' values from xc, which equals ' num2str(num*100/sizexc) ' % of the values.'])
	end
    
    out1 = V;
    out2 = txc;
    
	
	
case 'get&make'
	if nargin > 2
		[xc,adjts] = getv('get',ts,varargin{1});   
	else
		xc = getV('get',ts);
	end
	dim = size(ts,1);
	stats = getv('test',ts,xc);
	if exist('adjts') == 1, stats.adjustedts = adjts;,end
	[V,stats.taperedxc] = getv('make',xc,dim);
	
	out1 = V;
    out2 = stats;
    
end		% end switch



function [nointX,y] = adjusty(nointX,y);
	if ~isempty(nointX)
			% Fit predictors and temporal derivatives of regressors to be modeled out
			% Input: a convolved set of regressors, minus the intercept (noint), timeseries (y)
    	if size(nointX,2) > 1
 			[dummy,tempder] = gradient(nointX);	
		else
			[tempder] = gradient(nointX);	
      end
		nointX = [ones(size(nointX,1),1) nointX tempder];	
			% Return the residuals from the regression
    	y = y - (nointX * (nointX \ y));                  	% adjust y, removing effects of no interest
    	nointX = sparse(nointX);
	end
return
	