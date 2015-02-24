function [k,s] = hazard(t,y,a)

% HAZARD The Epanechnikov hazard function estimator
% from Trisha van Zandt

% column vector 't' is the points for which the hazard function is to be
% estimated
%
% column vector 'y' is the ordered RT data.
%
% optional constant 'a' determines the degree of smoothing.  

% If output argument 's' is specified, an estimated standard error of the
% hazard estimate 'k' (Silverman, 1986) will be returned.

if (nargin==2)
    c=.3;
elseif (nargin==3)
    c=a;
end

n=length(y);
h = c*min(std(y), iqr(y)/1.349)/n^.2;

fhat = mean( epanech( ((ones(length(t),1)*y')-(t*ones(1,n)))'/h ))'/h;
Fhat = mean(Iepanech( (t*ones(1,n)-ones(length(t),1)*y')'/h))';

k = fhat./(1-Fhat);

if (nargout==2)
    s = sqrt(.2683281571*k.^2./fhat/n/h);
end