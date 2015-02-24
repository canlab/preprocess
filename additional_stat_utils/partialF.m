function [varargout]=partialF(X,Y,r_params,varargin)

% Usage:
% p=partialF(X,Y,r_params,varargin);
% [p,F]=partialF(X,Y,r_params,varargin);
% [p,F,df]=partialF(X,Y,r_params,varargin);
% [p,F,df,Fmodel]=partialF(X,Y,r_params,varargin);
% [p,F,df,Fmodel,Rmodel]=partialF(X,Y,r_params,varargin);
% 
% X and Y are inputs to mult_reg, which are passed along with varargin.
% 
% r_params are the columns of X to remove from the full model in building
% the reduced model.
% 
% Performs an F-test on partial sums of squares to determine whether the
% parameters contained in Fmodel but not in Rmodel account for a
% significant amount of variance in the model.
% 
% edited 4/8/07 to be parallel with likelihoodRatio.m

Fmodel=mult_reg(X,Y,varargin{:});
x=X;
x(:,r_params)=[];
Rmodel=mult_reg(x,Y,varargin{:});

df(1)=Fmodel.p-Rmodel.p;
df(2)=Fmodel.n-Fmodel.p;
F=((Rmodel.SSE-Fmodel.SSE)/df(1))/(Fmodel.SSE/df(2));
p=ftest(F,df(1),df(2));

varargout{1}=p;varargout{2}=F;varargout{3}=df;varargout{4}=Fmodel;varargout{5}=Rmodel;
