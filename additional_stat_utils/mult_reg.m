function [varargout]=mult_reg(varargin)

% Usage:
% 
% results=mult_reg(X,Y)
% results=mult_reg(X,Y,names)
% [results,std_res]=mult_reg(...)
% [results,std_res]=mult_reg(...,'r',c)
% results=mult_reg(...,'O')
% results=mult_reg(...,W)
% results=mult_reg(...,W,'k')
% 
% 
% results is a structure containing parameter estimates, hypothesis tests,
% sums of squares, and other output relevant to the regression analysis.
% 
% X is either a matrix of all of the predictor variables with observations
% in rows, or the design matrix with the column of ones in the first
% column.
% 
% Y is the response variable.
% 
% names: is a cell array of strings with the names of all of the variables
% in the model, starting with Y and progressing through the X variables in
% the same order they appear (from left to right) in the X matrix.
% 
% The inclusion of the string character 'O' will result in regression
% through the origin.
% 
% Inclusion of a second output variable will produce an equivalent results
% structure for the standardized regression model.
% 
% Inclusion of the string character 'r', which must be followed by the
% constant c, will result in the std_res output being the results of a
% ridge regression with bias c. Note that the ridge regression results are
% necessarily for the standardized model.
% 
% W is a weight matrix for weighted least squares regression. Inclusion of
% the string character 'k' as an input variable will result in the
% variances being assumed to be known precisely.
% 
% 

origin=0;
known=0;
ridge=0;
ridge_go=0;
dfa=0;
vifoff=0;
for k=3:nargin
    if ischar(varargin{k})
        if strcmp(varargin{k},'O')
            origin=1;
        elseif strcmp(varargin{k},'k')
            known=1;
        elseif strcmp(varargin{k},'r')
            ridge=1;
            c=varargin{k+1};
        elseif strcmp(varargin{k},'R')
            ridge_go=1;
            c=varargin{k+1};
        elseif strcmp(varargin{k},'dfa')
            dfa=1;
        elseif strcmp(varargin{k},'vifoff')
            vifoff=1;
        end
    end
end

X=varargin{1};
if (origin==0)&&~isequal(X(:,1),ones(size(X,1),1))
    X=[ones(size(X,1),1) X];
end
Y=varargin{2};
if nargin>2
    for k=3:nargin
        if iscell(varargin{k})
            names=varargin{k};
        elseif isnumeric(varargin{k})&&isscalar(varargin{k})==0
            W=varargin{k};
        end
    end
end
if exist('names','var')==0
    names{1}='Response';
    for k=2:size(X,2)
        names{k}=['Predictor ' int2str(k-1)];
    end
end

if exist('W','var')==0
    W=eye(size(X,1));
end

Y=(W^(1/2))*Y;
X=(W^(1/2))*X;

if rcond((X'*X)^-1)<=10^-18
    warning('X''*X matrix is barely or not invertible. Aborting. Please fix your design matrix.');
    varargout{1}=struct('n',[],'p',[],'b',[],'e',[],'H',[],'Yhat',[],'df_R',[],'df_E',[],'SSTO',[],'SSE',[],'SSR',[],'MSR',[],'MSE',[],'var_cov_e',[],'var_cov_b',[],'aR2',[],'R2',[],'R',[],'corr_mat',[],'F_model',[],'p_model',[],'param',[],'y_name',[]);
    varargout{2}=varargout{1};
    return
end

n=size(Y,1);
if ~dfa
    p=size(X,2);
else
    p=size(X,2)+1;
end
if ~ridge_go
    b=((X'*X)^-1)*X'*Y;
else
    b=((X'*X+c*eye(size(X,2)))^-1)*X'*Y;
    if ~vifoff
        for k=1:size(X,2)
            if k~=size(X,2)
                [r,res]=mult_reg([X(:,1:k-1) X(:,k+1:end)],X(:,k),'r',c,'vifoff');
            else
                [r,res]=mult_reg(X(:,1:k-1),X(:,k),'r',c,'vifoff');
            end
            results.VIF(k)=(1-res.R2)^-1;
        end
    end
end
e=Y-X*b;
H=X*((X'*X)^-1)*X';
Yhat=H*Y;

df_R=p-1;
df_E=n-p;

SSE=e'*e;
MSE=SSE/df_E;


if ~known
    var_cov_b=MSE*((X'*X)^-1);
    var_cov_e=MSE*(eye(n)-H)'*(W^-1)*(eye(n)-H);
else
    var_cov_b=((X'*X)^-1);
    var_cov_e=(eye(n)-H)*(W^-1)*(eye(n)-H)';
end

if ~isequal(W,eye(size(X,1)))
    v = var_cov_b; v(1,:) = []; v(:,1) = [];
    SSR = b(2:end)'*v^-1*b(2:end);
    SSTO = SSR + SSE;
else
    SSR=b'*X'*Y-(1/n)*Y'*ones(n)*Y;
    SSTO=Y'*(eye(n)-(1/n)*ones(n))*Y;
end


MSR=SSR/df_R;

 
R2=SSR/SSTO;
R=sqrt(R2);
corr_mat=corr([Y X(:,2:end)]);

F_model=MSR/MSE;
p_model=1-fcdf(F_model,df_R,df_E);
for k=1:size(X,2)
    param(k).t=b(k)/sqrt(var_cov_b(k,k));
    if isinf(param(k).t)
        warning('T value is infinite')
        param(k).p=NaN;
    else
        if param(k).t>0
            param(k).p=(1-tcdf(param(k).t,df_E))*2;
        else
            param(k).p=tcdf(param(k).t,df_E)*2;
        end
    end
    if k==1&&~origin
        param(k).name='Intercept';
    elseif origin
        param(k).name=names{k+1};
    else
        param(k).name=names{k};
    end
end

if ~isequal(W,eye(size(X,1)))
    Yhat=(W^(-1/2))*Yhat;
    e=(W^(-1/2))*e;
end

results.n=n;
results.p=p;
results.b=b;
results.e=e;
results.H=H;
results.Yhat=Yhat;
results.df_R=df_R;
results.df_E=df_E;
results.SSTO=SSTO;
results.SSE=SSE;
results.SSR=SSR;
results.MSR=MSR;
results.MSE=MSE;
results.var_cov_e=var_cov_e;
results.var_cov_b=var_cov_b;
results.aR2=1-((n-1)/(n-p))*(SSE/SSTO);
results.R2=R2;
results.R=R;
results.corr_mat=corr_mat;
results.F_model=F_model;
results.p_model=p_model;
results.param=param;
results.y_name=names{1};

varargout{1}=results;

if nargout>1
    Y=(1/sqrt(n-1))*((Y-mean(Y))/std(Y));
    for k=2:size(X,2)
        X(:,k)=(1/sqrt(n-1))*((X(:,k)-mean(X(:,k)))/std(X(:,k)));
    end
    X=X(:,2:end);
    if ridge==0
        std_res=mult_reg(X,Y,names,'dfa');
    elseif ~vifoff
        std_res=mult_reg(X,Y,names,'dfa','R',c);
    else
        std_res=mult_reg(X,Y,names,'dfa','R',c,'vifoff');
    end
    varargout{2}=std_res;
end

