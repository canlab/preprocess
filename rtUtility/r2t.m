function [t,p,sig] = r2t(r,n)
% function r2t(r,n)
% pearson's r to t-value with n - 2 df
% http://www2.chass.ncsu.edu/garson/pa765/correl.htm
r2tfun = inline('r.*((n-2).^.5) ./ ((1-r.^2).^.5)','r','n');

t = r2tfun(r,n);
p = 2*(1 - tcdf(abs(t),n-2));

if p < .05, sig = 1; else, sig = 0;,end

return
