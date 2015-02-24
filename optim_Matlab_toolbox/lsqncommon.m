function [x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = lsqncommon(FUN,x,YDATA,LB,UB,options,defaultopt,caller,computeLambda,varargin)
%LSQNCOMMON Solves non-linear least squares problems.
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN]=...
%      LSQNCOMMON(FUN,X0,YDATA,LB,UB,OPTIONS,DEFAULTOPT,CALLER,COMPUTELAMBDA,XDATA,VARARGIN...) 
%   contains all the setup code common to both LSQNONLIN and LSQCURVEFIT to call either the 
%   large-scale SNLS or the medium-scale NLSQ.

%   Copyright 1990-2001 The MathWorks, Inc. 
%   $Revision: 1.6 $  $Date: 2001/03/27 19:55:16 $

xstart=x(:);
numberOfVariables=length(xstart);
% Note: XDATA is bundled with varargin already for lsqcurvefit 
lenVarIn = length(varargin);

large = 'large-scale';
medium = 'medium-scale';

switch optimget(options,'Display',defaultopt,'fast')
case {'off','none'}
    verbosity = 0;
case 'iter'
    verbosity = 2;
case 'final'
    verbosity = 1;
case 'testing'
    verbosity = Inf;
otherwise
    verbosity = 1;
end

[xstart,l,u,msg] = checkbounds(xstart,LB,UB,numberOfVariables);
if ~isempty(msg)
    EXITFLAG = -1;
    [Resnorm,FVAL,OUTPUT,LAMBDA,JACOB] = deal([]);
    x(:)=xstart;
    if verbosity > 0
        disp(msg)
    end
    return
end
lFinite = ~isinf(l);
uFinite = ~isinf(u);

if min(min(u-xstart),min(xstart-l)) < 0
    xstart = startx(u,l); 
end

diagnostics = isequal(optimget(options,'Diagnostics',defaultopt,'fast'),'on');
gradflag =  strcmp(optimget(options,'Jacobian',defaultopt,'fast'),'on');
line_search = strcmp(optimget(options,'LargeScale',defaultopt,'fast'),'off'); % 0 means large-scale, 1 means medium-scale
mtxmpy = optimget(options,'JacobMult',[]); % use old
if isequal(mtxmpy,'atamult')
    warnstr = sprintf('%s\n%s\n%s\n', ...
        'Potential function name clash with a Toolbox helper function:',...
        'Use a name besides ''atamult'' for your JacobMult function to',...
        'avoid errors or unexpected results.');
    warning(warnstr)
end

% Convert to inline function as needed
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
    [funfcn, msg] = fprefcnchk(FUN,caller,lenVarIn,gradflag);
else
    errmsg = sprintf('%s\n%s', ...
        'FUN must be a function or an inline object;', ...
        ' or, FUN may be a cell array that contains these type of objects.');
    error(errmsg)
end

fuser = [];  
JAC = [];
x(:) = xstart;
try
    switch funfcn{1}
    case 'fun'
        fuser = feval(funfcn{3},x,varargin{:});
    case 'fungrad'
        [fuser,JAC] = feval(funfcn{3},x,varargin{:});
    case 'fun_then_grad'
        fuser = feval(funfcn{3},x,varargin{:}); 
        JAC = feval(funfcn{4},x,varargin{:});
    otherwise
        errmsg = sprintf('%s%s\n','Undefined calltype in ',upper(caller));
        error(errmsg);
    end
catch
    if ~(isequal(funfcn{1},'fun_then_grad')) | isempty(fuser)
        badfunfcn = funfcn{3};
    else % 'fun_then_grad & ~isempty(fuser) (so it error'ed out on the JAC call)
        badfunfcn = funfcn{4};
    end   
    if isa(badfunfcn,'inline')
        errmsg = sprintf(['User supplied %s ==> %s\n' ...
                'failed with the following error:\n\n%s'],...
            'expression or inline function',formula(badfunfcn),lasterr);
    else % function (i.e. string name of), function handle, inline function, or other
        % save last error in case call to "char" errors out
        lasterrmsg = lasterr;
        % Convert to char if possible
        try
            charOfFunction = char(badfunfcn);
        catch
            charOfFunction = '';
        end
        if ~isempty(charOfFunction)
            errmsg = sprintf(['User supplied %s ==> %s\n' ...
                    'failed with the following error:\n\n%s'],...
                'function',charOfFunction,lasterrmsg);
        else
            errmsg = sprintf(['User supplied %s ' ...
                    'failed with the following error:\n\n%s'],...
                'function',lasterrmsg);
        end
    end
    error(errmsg)
end

if isequal(caller,'lsqcurvefit')
    if ~isequal(size(fuser), size(YDATA))
        error('Function value and YDATA sizes are incommensurate.')
    end
    fuser = fuser - YDATA;  % preserve fuser shape until after subtracting YDATA 
end

f = fuser(:);
nfun=length(f);

if gradflag
    % check size of JAC
    [Jrows, Jcols]=size(JAC);
    if isempty(mtxmpy) 
        % Not using 'JacobMult' so Jacobian must be correct size
        if Jrows~=nfun | Jcols ~=numberOfVariables
            errstr = sprintf('%s\n%s%d%s%d\n',...
                'User-defined Jacobian is not the correct size:',...
                '    the Jacobian matrix should be ',nfun,'-by-',numberOfVariables);
            error(errstr);
        end
    end
else
    Jrows = nfun; 
    Jcols = numberOfVariables;   
end

% trustregion and enough equations (as many as variables) 
if ~line_search & nfun >= numberOfVariables 
    OUTPUT.algorithm = large;
    
    % trust region and not enough equations -- switch to line_search
elseif ~line_search & nfun < numberOfVariables 
    warnstr = sprintf('%s\n%s\n', ...
        'Large-scale method requires at least as many equations as variables; ',...
        '   switching to line-search method instead.  Upper and lower bounds will be ignored.');
    warning(warnstr);
    OUTPUT.algorithm = medium;
    
    % line search and no bounds  
elseif line_search & isempty(l(lFinite)) & isempty(u(uFinite))
    OUTPUT.algorithm = medium;
    
    % line search and  bounds  and enough equations, switch to trust region 
elseif line_search & (~isempty(l(lFinite)) | ~isempty(u(uFinite))) & nfun >= numberOfVariables
    warnstr = sprintf('%s\n%s\n', ...
        'Line-search method does not handle bound constraints; ',...
        '   switching to large-scale method instead.');
    warning(warnstr);
    OUTPUT.algorithm = large;
    
    % can't handle this one:   
elseif line_search & (~isempty(l(lFinite)) | ~isempty(u(uFinite)))  & nfun < numberOfVariables
    errstr = sprintf('%s\n%s\n%s\n', ...
        'Line-search method does not handle bound constraints ',...
        '   and large-scale method requires at least as many equations as variables; ',...
        '   aborting.');
    error(errstr);
end

if diagnostics > 0
    % Do diagnostics on information so far
    constflag = 0; gradconstflag = 0; non_eq=0;non_ineq=0;lin_eq=0;lin_ineq=0;
    confcn{1}=[];c=[];ceq=[];cGRAD=[];ceqGRAD=[];
    hessflag = 0; HESS=[];
    msg = diagnose(caller,OUTPUT,gradflag,hessflag,constflag,gradconstflag,...
        line_search,options,defaultopt,xstart,non_eq,...
        non_ineq,lin_eq,lin_ineq,l,u,funfcn,confcn,f,JAC,HESS,c,ceq,cGRAD,ceqGRAD);
end

% Execute algorithm
if isequal(OUTPUT.algorithm,large)
    if ~gradflag % provide sparsity of Jacobian if not provided.
        Jstr = optimget(options,'JacobPattern',[]);
        if isempty(Jstr)  
            % Put this code separate as it might generate OUT OF MEMORY error
            Jstr = sparse(ones(Jrows,Jcols));
        end
        if ischar(Jstr) 
            if isequal(lower(Jstr),'sparse(ones(jrows,jcols))')
                Jstr = sparse(ones(Jrows,Jcols));
            else
                error('Option ''JacobPattern'' must be a matrix if not the default.')
            end
        end
    else
        Jstr = [];
    end
    [x,FVAL,LAMBDA,JACOB,EXITFLAG,OUTPUT,msg]=...
        snls(funfcn,x,l,u,verbosity,options,defaultopt,f,JAC,YDATA,caller,Jstr,computeLambda,varargin{:});
else
    [x,FVAL,JACOB,EXITFLAG,OUTPUT,msg] = ...
        nlsq(funfcn,x,verbosity,options,defaultopt,f,JAC,YDATA,caller,varargin{:});
    LAMBDA.upper=[]; LAMBDA.lower=[];   
end
Resnorm = FVAL'*FVAL;
if verbosity > 0
    disp(msg);
end


% Reset FVAL to shape of the user-function, fuser
FVAL = reshape(FVAL,size(fuser));

%--end of lsqncommon--

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [allfcns,msg] = fprefcnchk(funstr,caller,lenVarIn,gradflag)
%PREFCNCHK Pre- and post-process function expression for FUNCHK.
%   [ALLFCNS,MSG] = PREFUNCHK(FUNSTR,CALLER,lenVarIn,GRADFLAG) takes
%   the (nonempty) expression FUNSTR from CALLER with LenVarIn extra arguments,
%   parses it according to what CALLER is, then returns a string or inline
%   object in ALLFCNS.  If an error occurs, this message is put in MSG.
%
%   ALLFCNS is a cell array: 
%    ALLFCNS{1} contains a flag 
%    that says if the objective and gradients are together in one function 
%    (calltype=='fungrad') or in two functions (calltype='fun_then_grad')
%    or there is no gradient (calltype=='fun'), etc.
%    ALLFCNS{2} contains the string CALLER.
%    ALLFCNS{3}  contains the objective function
%    ALLFCNS{4}  contains the gradient function (transpose of Jacobian).
%  
%    NOTE: we assume FUNSTR is nonempty.
% Initialize
msg='';
allfcns = {};
funfcn = [];
gradfcn = [];

if gradflag
    calltype = 'fungrad';
else
    calltype = 'fun';
end

% {fun}
if isa(funstr, 'cell') & length(funstr)==1
    % take the cellarray apart: we know it is nonempty
    if gradflag
        calltype = 'fungrad';
    end
    [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
    if ~isempty(msg)
        error(msg);
    end
    
    % {fun,[]}      
elseif isa(funstr, 'cell') & length(funstr)==2 & isempty(funstr{2})
    if gradflag
        calltype = 'fungrad';
    end
    [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
    if ~isempty(msg)
        error(msg);
    end  
    
    % {fun, grad}   
elseif isa(funstr, 'cell') & length(funstr)==2 % and ~isempty(funstr{2})
    
    [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
    if ~isempty(msg)
        error(msg);
    end  
    [gradfcn, msg] = fcnchk(funstr{2},lenVarIn);
    if ~isempty(msg)
        error(msg);
    end
    calltype = 'fun_then_grad';
    if ~gradflag
        warnstr = ...
            sprintf('%s\n%s\n%s\n','Jacobian function provided but OPTIONS.Jacobian=''off'';', ...
            '  ignoring Jacobian function and using finite-differencing.', ...
            '  Rerun with OPTIONS.Jacobian=''on'' to use Jacobian function.');
        warning(warnstr);
        calltype = 'fun';
    end   
    
elseif ~isa(funstr, 'cell')  %Not a cell; is a string expression, function name string, function handle, or inline object
    [funfcn, msg] = fcnchk(funstr,lenVarIn);
    if ~isempty(msg)
        error(msg);
    end   
    if gradflag % gradient and function in one function/M-file
        gradfcn = funfcn; % Do this so graderr will print the correct name
    end  
else
    errmsg = sprintf('%s\n%s', ...
        'FUN must be a function object or an inline object;', ...
        ' or, FUN may be a cell array that contains these type of objects.');
    error(errmsg)
end

allfcns{1} = calltype;
allfcns{2} = caller;
allfcns{3} = funfcn;
allfcns{4} = gradfcn;
allfcns{5}=[];

