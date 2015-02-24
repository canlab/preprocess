function gls = genls(X,y,varargin)
% function gls = genls(X,y, [optional: tolerance])
%
% The genls function uses a whitening strategy to remove the autocorrelation
% from a timeseries (and model) before fitting the model.
%
% It iterates until the change in the residual autocorrelation is less than
% a specified tolerance, which is equal to the average % difference per element 
% of the autocorrelation function.
% 	(1 is 100% difference, 0 is no difference between the elements on average.)
% 
% The default tolerance is 1% (.01)
%
% Also, the intercept should be the first column of X, or you'll get a warning.

format compact
% X is predictor matrix, each col. is a predictor
numpreds = size(X,2);
numpoints = size(y,1);
if not(mean(X(:,1) == 1) & max(X(:,1) == 1)), warning('No intercept in X!'), end;

% set up variables
converged = 0;
if nargin > 2, threshold = varargin{1};,else threshold = 0.01;,end
maxiterations = 10;
niterations = 0;
index = 1;

% set up ridge regression k
k = .01; if nargin > 3, k = varargin{2};,end

disp(['computing initial model with ridge regression, kappa = ' num2str(k) ' ...'])
B = ridge(y,X,k);
e = y - X*B;

% get autocorrelation of residuals	(V or sigma^-1)
xc = getv('get',e);
V = getv('make',xc,size(e,1));
% whiten data
%[si,p] = chol(V);	PROBLEM: not positive definite, can't do this.
%si = real(inv(V^.5));
si = inv(chol(V));

sX = si * X;
sy = si * y;

xcold = xc;

disp('iterating generalized least squares...')
while not(converged)
   % fit model: calculate new autocorrelation  
   B = ridge(sy,sX,k);
	e = sy - sX*B;
   	xc = getv('get',e);
	% gradual adjustment of acf
	xc(2:end) = xcold(2:end) - .2 * (xcold(2:end) - xc(2:end));  
	
	% calculate difference of old and new autocorrelation
	diff = mean(abs(xc-xcold))./mean(abs(xcold));		% avg. percent difference
	gls.diff(index) = diff; index = index + 1;
	if diff < threshold | niterations == maxiterations, converged = 1; 
	else
		% get new sx and sy based on new autocorrelation
		V = getv('make',xc,size(e,1));
		si = inv(chol(V));
		%si = real(inv(V^.5));
		sX = si * X;
		sy = si * y;
		xcold = xc;
   end
   niterations = niterations + 1;
   disp(['iteration ' num2str(niterations)])
end

if gls.diff(end) <= threshold
	disp(['converged in ' num2str(niterations) ' iterations'])
else
	disp(['...did not converge'])
end
disp(['at a threshold of ' num2str(threshold*100) ' %'])

% final model parameters
xtxi = inv(X'*X);
betas = B;
fits = X * betas;
df = size(X,1) - size(X,2);
var = (e' * e) / df;			% (small) sigma squared.
varbetas = xtxi .* var;
sebetas = sqrt(diag(varbetas));
tbetas = betas ./ sebetas;
pbetas = tdist(tbetas',df)';

gls.niterations = niterations;
gls.threshold = threshold;
gls.X = X;
gls.y = y;
gls.e = e;
gls.betas = betas';
gls.sebetas = sebetas';
gls.tbetas = tbetas';
gls.pbetas = pbetas';
gls.var = var;					% (small) sigma squared.
gls.df = df;
gls.V = sparse(V);				% var/cov matrix of errors
gls.xc = xc;
gls.fits = fits;
gls.xtxi = xtxi;


return

