function [varargout]=likelihoodRatioFitTest(X,Y,r_params)

% Usage:
% p=likelihoodRatio(X,Y,r_params);
% [p,chisq]=likelihoodRatio(X,Y,r_params);
% [p,chisq,df]=likelihoodRatio(X,Y,r_params);
% [p,chisq,df,Fstats]=likelihoodRatio(X,Y,r_params);
% [p,chisq,df,Fstats,Rstats]=likelihoodRatio(X,Y,r_params);
% 
% X and Y are prepared exactly as they would be for inputs to glmfit. This
% script will call glmfit(X,Y,'binomial') and glmfit(Xr,Y,'binomial') where
% Xr=X; Xr(:,r_params)=[];
% 
% As implied by the above, r_params is a vector indicating which of the
% columns of X to remove for the reduced model.
% 
% Please note that the design matrix X must *NOT* have a column of 1's for
% the intercept parameter, consistent with calls to glmfit.
% 
% Performs a likelihood ratio test to determine whether the parameters
% contained in Fmodel but not in Rmodel account for a significant amount of
% variance in the response variable.
% 
% Outputs in order are:
% 1. p-value for the test
% 2. test statistic for the test
% 3. df for the test
% 4. stats structure output from glmfit for the full model
% 5. stats structure output from glmfit for the reduced model

Xr=X;Xr(:,r_params)=[];
[B,DEV,Fstats]=glmfit(X,Y,'binomial');
[B,DEV,Rstats]=glmfit(Xr,Y,'binomial');
df=Rstats.dfe-Fstats.dfe;

if size(Y,2)==2
%     Y=Y(:,1)./Y(:,2);
    x=[];y=[];
    for k=1:length(Y)
        for j=1:Y(k,2)
            x(end+1,:)=X(k,:);
        end
        y(end+1:end+Y(k,2))=[ones(Y(k,1),1); zeros(Y(k,2)-Y(k,1),1)];
    end
    X=x;
    Xr=x;Xr(:,r_params)=[];
    Y=y';
end

X=[ones(size(X,1),1) X];
Xr=[ones(size(Xr,1),1) Xr];
Lf=exp(Y'*(X*Fstats.beta)-sum(log(1+exp(X*Fstats.beta))));
Lr=exp(Y'*(Xr*Rstats.beta)-sum(log(1+exp(Xr*Rstats.beta))));

chisq=-2*log(Lr/Lf);
p=1-chi2cdf(chisq,df);

varargout{1}=p;varargout{2}=chisq;varargout{3}=df;varargout{4}=Fstats;varargout{5}=Rstats;
