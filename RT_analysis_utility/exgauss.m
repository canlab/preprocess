function fit = exgauss(t,theta)

% ex-Gaussian density function
%
% from Trisha van Zandt's chapter appendix
% gives the density function for the ex-Gaussian distribution
%
% input arguments: 't' is a column vector of times for which the desnity is
% to be computed; 'theta' is a vector of three elements, mu, sigma, and tau
% (in that order)
%
% returns column vectors of the same length as 't'


mu = theta(1); sigma=theta(2); tau=theta(3);
part1 = exp(-t./tau + mu./tau + sigma.^2./2./tau.^2);
part2 = normcdf((t-mu-sigma.^2./tau)./sigma)./tau;
fit = part1.*part2;


