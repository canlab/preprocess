function [estimate,sesums,z,p] = correl_contrast_Steiger(b,e,c,d,r,n)
% computes Z-test of contrast [b + e - c - d] for correlation matrix r
% b,e,c,d are 2-element vectors of indices of r for each specified
% correlation
% n is sample size
%
% returns: 
% estimate of contrast
% standard error of estimate
% z-value for test
% p-value for test (1-tailed)
%
% example:
% r = [1 .25 .3 .4; .25 1 .25 .6; .3 .25 1 .3; .4 .6 .3 1]
% n = 100;
% b = [3 1]; e = [4 2]; c = [4 1]; d = [3 2];

global R
R = r;

% variance of the contrast across correlation values
% ::
% contrast values ^ 2 (algebraically)
% [1 1 -1 -1]^2 or [b e -c -d] becomes 
% b^2 + e^2 + c^2 + d^2 + 2be + 2cd - 2(other pairs)
%
% replace ^2 with variance and * with covariance
% use PearsonFilon for the estimates of the covariance/variance of two
% correlations
% variance of the contrast across correlation values

varSums = co(b,b) + co(e,e) + co(c,c) + co(d,d) + ...
          2*co(b,e) + 2*co(c,d) ...
          -2*co(b,c) - 2*co(e,c) - 2*co(b,d) - 2*co(e,d);
sesums = sqrt(varSums ./ n);

fprintf(1,'(%3.2f + %3.2f) - (%3.2f + %3.2f), n = %3.0f\t',r(b(1),b(2)),r(e(1),e(2)),r(c(1),c(2)),r(d(1),d(2)));


estimate = r(b(1),b(2)) + r(e(1),e(2)) - r(c(1),c(2)) - r(d(1),d(2));
%estimate = estimate .^ 2;

z = estimate ./ sesums;
p = 2* (min(normcdf(z),1 - normcdf(z)));     % one-tailed

fprintf(1,'est = %3.3f, SE = %3.3f, Z = %3.2f, p = %3.4f\n',estimate,sesums,z,p);

return


function cv = co(x,y)
% x and y are 2-vector indices of matrix to compute covariance for
% cv is covariance
cv = pearsonf(x(1),x(2),y(1),y(2));
return


function pf = pearsonf(i,j,k,h)

% Pearson-Filon: covariance (or var) of element i,j with element k,h
% depending on correlation values

global R

 pf = (1/2) .* R(i,j) .* R(k,h) .* ...
 ( R(i,k).^2 + R(i,h).^2 + R(j,k).^2 + R(j,h).^2 ) + ...
        R(i,k) .* R(j,h) + R(i,h) .* R(j,k) - ...
        R(i,j) .* ( R(j,k) .* R(j,h) + R(i,k) .* R(i,h) ) - ... 
        R(k,h) .* ( R(j,k) .* R(i,k) + R(j,h) .* R(i,h) );
 
        
return



