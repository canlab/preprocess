function [x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = lsqcurvefit_t(FUN,x,XDATA,YDATA,LB,UB,options,varargin)

% ------------Initialization----------------
defaultopt = struct('Display','final','LargeScale','on', ...
   'TolX',1e-6,'TolFun',1e-6,'DerivativeCheck','off',...
   'Diagnostics','off',...
   'Jacobian','off','JacobMult',[],...% JacobMult set to [] by default
   'JacobPattern','sparse(ones(Jrows,Jcols))',...
   'MaxFunEvals','100*numberOfVariables',...
   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
   'PrecondBandWidth',0,'TypicalX','ones(numberOfVariables,1)',...
   'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
   'TolPCG',0.1,'MaxIter',400,...
   'LineSearchType','quadcubic','LevenbergMarquardt','on'); 
% If just 'defaults' passed in, return the default options in X
if nargin==1 & nargout <= 1 & isequal(FUN,'defaults')
   x = defaultopt;
   return
end

 if nargin < 7, options= [];
   if nargin < 6, UB = [];
      if nargin < 5, LB = [];
         if nargin < 4, 
            error('LSQCURVEFIT requires four input arguments');
         end, end, end, end
if nargout > 5
   computeLambda = 1;
else 
   computeLambda = 0;
end


% XDATA put at the end so it's bundled with varargin in lsqncommon
caller = 'lsqcurvefit';
[x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = ...
   lsqncommon(FUN,x,YDATA,LB,UB,options,defaultopt,caller,...
              computeLambda,XDATA,varargin{:});

