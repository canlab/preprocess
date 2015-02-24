function b = beta_scale(b,baseind,meanind)
% function b = beta_scale(b,baseind,meanind)
%
% Tor Wager
% 
% Scales model fit betas to percent change from mean signal
% Input:
%   b, a vector of betas for a particular voxel
%   baseind, a vector of column indices for resting baseline
%   meanind, a vector of column indices for session intercept terms
% 
% Process:
%   b - (mean intercept b + mean baseline b)
%   Assumes baseline and intercepts use effects coding [1's and 0's]
%   so that the magnitude of the b reflects the mean signal intensity.
%   Convolving a baseline reg, for example, will violate this assumption.
%
% Output:
%   Betas scaled to % change from resting baseline value (*mean timeseries)

% note: Tom Nichols says not to use baseind, because it's colinear w. intercept
% just use mean (intercepts) and baseline scans will contribute to that.

bm = mean(b(meanind)); % + mean(b(baseind));

if ~isnan(bm) & ~isinf(bm) & bm ~= 0  
	b = b ./ bm;
else
	warning('Baseline estimate is 0 or Inf')
	b = b ./ NaN;
end

return
