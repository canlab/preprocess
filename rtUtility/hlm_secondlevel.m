function hlm_secondlevel(Zm,varZm)
%hlm_secondlevel(Zm,varZm)
%
% right now, a hierarchical 1-sample t-test

% Zm is subject estimates
% varZm is vector of variances for each subject

nsub = size(Zm,1);

% -----------------------------------------------------------------------
% variance components and subject weights
% -----------------------------------------------------------------------

varbtwn = repmat(max(0,(var(Zm)-mean(varZm))),nsub,1);                      % Estimate of the between subject variance assuming white noise: MSA - mean(MSE)

varwplusb = varZm+varbtwn;
Wm = (1./varwplusb)./repmat(sum(1./varwplusb),nsub,1);        % Calculate weights



% -----------------------------------------------------------------------
% weighted (hierarchical) regression
% -----------------------------------------------------------------------
grpcontrast = [];
% rb, covb are correl and covariance matrices across time for Zpop
% contrast across means (group difference) if grpcontrast is not empty
[Zpop,df,ZPdiff,dfdiff,regstat] = weighted_reg(Zm,Wm,varwplusb,grpcontrast);



return

