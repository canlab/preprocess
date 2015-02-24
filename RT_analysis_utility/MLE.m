function F = MLE(p,func,y)

% p is a vector of starting values for parameters to be estimated
% func is a character variable specifying the theoretical density
% function, e.g. 'exgauss'
% y is a vector of observed RTs;

% call as, e.g., X = fminsearch ('MLE',p,[],func,y);

F = -sum(log(feval(func,y,p)));