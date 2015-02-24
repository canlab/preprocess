function [msb,msw,m,ssb,ssw,dfb,dfw] = sumsquares(d)
% [msb,msw,m,ssb,ssw] = sumsquares(d)
% d is data, rows are cases, cols are observations
% msb and msw are mean squared error terms, using statistical df (n-1, n-k)
%
% ssb is sums of squares between (explained by column means)
% ssw is error - residuals of x - avg(x)
%
% Centers rows first, to remove main effects of case
% see also sumsqr.m

wh = find(any(isnan(d),2));
if ~isempty(wh), 
    warning('Some NaN values found in data! Eliminating row-wise.')
    d(wh,:) = [];
end

casem = mean(d,2);
casem = casem(:,ones(1,size(d,2)));

d = d - casem;

m = mean(d);
ssb = sum(m .^ 2);
dfb = length(m) - 1;
msb = ssb ./ dfb;

m = m(ones(size(d,1),1),:);

gm = mean(d(:));
dfw = length(d(:)) - size(d,2); % df within, n - k

ssw = (d - m) .^ 2; msw = sum(ssw(:)) ./ dfw;

m = m(1,:);

return

