function [params, height, delay] = ts_stats(DX,y,rint,rbase,varargin)
% ts_stats(DX,y,rint,rbase,[opt] gammfunct)
% by Tor Wager
%
% DX        is deconvolution model matrix to fit to timeseries data
% y         is timeseries data for whole experiment
% rint      is vector of indices of DX columns of interest (all params of interest)
% rbase     is vector of indices of DX columns for baseline estimate
% gammfunct is the function to fit nonlinearly; if empty, skips this.
%
% This function does:
% - Fits DX matrix to data to get parameter estimates
% - Subtracts param estimates from baseline estimates, if specified
% - Returns DX param estimates of interest for group ANOVA
% - Fits 2-parameter (height + delay) gamma function to param est. of interest
% - Returns fitted height for each param of interest
% - Returns fitted delay for each param of interest
%
%
% y and DX should be filtered and adjusted as appropriate before running this function.

height = []; delay = []; gammfunct = [];
if length(varargin) > 0, gammfunct = varargin{1};,end

b = pinv(DX) * y;

% done in tor_setup_anova.m
%a = repmat(rint,numframes,1);
%b = repmat((1:numframes)'-1,1,size(rint,2));
%c = a + b; 
%c = c(:);

if isempty(rbase)
    params = b(rint);
else
    params = b(rint) - mean(b(rbase));
end

if ~isempty(gammfunct)
    % then fit the function
    % nlfit...
end

return

