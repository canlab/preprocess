function [b,X] = trial_level_beta(y,varargin)
%
%
% [b,X] = trial_level_beta(y,onsets,k,DX) 
%
% future:
% trial_level_beta(y,onsets,hrf)    % take existing HRF
% trial_level_beta(y,onsets,k)      % build DX and estimate
% onsets for multiple types?

onsets = varargin{1};
k = varargin{2};
DX = varargin{3};

% get HRF from data and DX
% ---------------------------------------------------------------------
hrf = pinv(DX) * y;


% make linear model for each trial
% ---------------------------------------------------------------------
[n,c] = size(DX);    % n is number of time points
t = length(onsets);  % t is number of trials

I = DX(:,k+1:end);   % intercepts for each session

X = zeros(n,t);
for i = 1:t, X(onsets(i),i) = 1;, end

X = getPredictors(X,hrf);
X = [X I];


% fit trial-level model to get betas
% ---------------------------------------------------------------------

b = pinv(X) * y;


