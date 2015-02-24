function [yboot h_hat] = bootstrapme(y)

%bootstrapping routine from Trisha van Zandt "Analysis of Response Time
%Distributions" chapter appendix

% y is a column vector containing the ordered RTs
% yboot will contain the same number of observations as y, sampled with
% replacement from y.

yboot = y([fix(rand(length(y),1)*length(y))+1]);

% ------------- for computing SE of hazard function estimate ----------


% columns in matrix h_hat will contain the estimated hazard functions for
% each bootstrapped sample

for i = 1:n_samples,
    yboot = y([fix(rand(length(y),1)*length(y))+1]);
    h_hat(:,i) = hazard(t,yboot);
end