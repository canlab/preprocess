function [invhess,directn]=updhess(xnew,xold,gradxnew,gradxold,invhess,para)
%UPDHESS Performs the Inverse Hessian Update.
%   Returns direction of search for use with 
%   unconstrained optimization problems. 

%   Copyright 1990-2001 The MathWorks, Inc.
%   $Revision: 1.10 $  $Date: 2001/03/27 19:55:14 $
%   Andy Grace 7-9-90.

u=xnew-xold;
v=gradxnew-gradxold;
if para(1,6)==0
% The BFGS Hessian Update formula:
    invhess=invhess + v*v'/(v'*u)  -invhess*u*u'*invhess'/(u'*invhess*u);
    directn=-invhess\gradxnew;

elseif para(1,6)==1

% The DFP formula
    a=u*u'/(u'*v);
    b=-invhess*v*v'*invhess'/(v'*invhess*v);
    invhess=invhess + a + b;
    directn=-invhess*gradxnew;

elseif para(1,6)==3

% A formula given by Gill and Murray
    a = 1/(v'*u);
    invhess=invhess - a*(invhess*v*u'+u*v'*invhess)+a*(1+v'*invhess*v*a)*u*u' ;
    directn=-invhess*gradxnew;
elseif para(1,6)==2
% Steepest Descent
    directn=-gradxnew;
end
       
