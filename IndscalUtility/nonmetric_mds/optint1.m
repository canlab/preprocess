function [xold,fold,invhess,para]=optint(xnew,fnew,para)
%OPTINT Function to initialize FMINU routine.

%   Copyright 1990-2001 The MathWorks, Inc.
%   $Revision: 1.13 $  $Date: 2001/03/27 19:55:13 $
%   Andy Grace 7-9-90.
lenx=length(xnew);
invhess=eye(lenx);  
xold=xnew;
fold=fnew;
para=foptions(para);
if para(14)==0, para(14)=lenx*100;end 
if para(1)>0, para, end
if para(1)>0,
    disp('')
    if isinf(para(1))
      disp('f-COUNT   FUNCTION    STEP-SIZE      GRAD/SD  LINE-SEARCH')
    else
      disp('f-COUNT   FUNCTION    STEP-SIZE      GRAD/SD')
    end
end