function [difference,se,dfe,ratio,p] = ttest2_from_variance(mx,my,s2x,s2y,nx,ny)
% [difference,se,df,ratio,p] = function ttest2_from_variance(mx,my,s2x,s2y,nx,ny)
%
% do a 2-sample t-test, given variances and means (not raw data)
% do not pool variances, use Sattherwaite approx. for df
%
% mx, my are mean of x and y, s2 are variances

difference = mx - my;

% var of mean estimates
s2xbar = s2x ./ nx;
s2ybar = s2y ./ ny;

% sattherwaite
dfe = (s2xbar + s2ybar) .^2 ./ (s2xbar.^2 ./ (nx-1) + s2ybar.^2 ./ (ny-1));

% se diff
se = sqrt(s2xbar + s2ybar);

% t-score
ratio = difference ./ se;

% p-value
p = 1-tcdf(abs(ratio),dfe);

fprintf(1,'u1 = %3.2f, u2 = %3.2f, udiff = %3.2f, se = %3.2f, t(%3.1f) = %3.2f, p = %3.4f\n',mx,my,difference,se,dfe,ratio,p);

return