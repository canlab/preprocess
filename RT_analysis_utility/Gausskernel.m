function fhat = Gausskernel(t,y,a)

% Gaussian kernel estimator

% from Trisha van Zandt's "Analysis of Response Time Distribution" chapter
%
% given column vector 'y' containing the ordered RTs, a column vector 't'
% containing the time points at which the estimate is desired, and a
% smoothing parameter 'h', the Gaussian kernel estimator is computed with
% the function Gausskernel.

% input arguments: 't' are the time points, 'y' is the ordered RTs, and 'a'
% is an optional multiplicative constant for the bandwidth parameter h.

if (nargin==2)
    c=.9;
elseif(nargin==3)
    c=a;
end

h = c*min(std(y),iqr(y)/1.349)/length(y)^.2;
fhat = mean(normpdf((ones(length(t),1) * y' - t*ones(1,length(y)))',0,h))';