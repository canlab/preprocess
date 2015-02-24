function F = Iexgauss(t,theta)

% IEXGAUSS The ex-Gaussian CDF
% from Trisha van Zandt

% syntax F=Iexguass(t,theta), where t is a vector of times for which the
% CDF is to be computed, and theta is a vector of exgaussian parameters
% (mu, sigma, tau)

mu = theta(1); sigma=theta(2); tau = theta(3);
part1 = -exp(-t./tau + mu./tau + sigma.^2./2./tau.^2);
part2 = normcdf((t-mu-sigma.^2./tau)./sigma);
part3 = normcdf((t-mu)/sigma);
F = part1.*part2 + part3;