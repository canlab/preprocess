function results=reg_val(varargin)
 
% 
% This function prints output in a tabular format for each of the requested
% criteria. It is suggested that you copy and paste this output into a
% spreadsheet program (e.g. excel) so that it is appropriately arranged
% into columns. Note that for each criterion, models are sorted such that
% models with minimal values of the criterion appear in the first row, and
% models with maximal values of the criterion appear in the last row
% (reversed for R2 and aR2).
% 
% 
% Usage:
% 
% results=reg_val(X,Y,criterion)
% results=reg_val(...,num_models)
% results=reg_val(...,names)
% 
% results is a structure containing...
% 
% X is the design matrix with the column of ones in the first
% column. This matrix must include the full set of candidate parameters for
% inclusion in the model.
% 
% Y is the response variable.
% 
% criterion is a cell array of strings indicating the criteria to be used
% for model selection. Valid criteria are (see page 353 of Kutner,
% Nachtsheim, & Neter, "Applied Linear Regression Models: 4th Ed."):
% 'SSE', 'R2', 'aR2', 'MSE', 'C', 'AIC', 'SBC', 'PRESS'
% 
% criterion MUST be a cell even if it only has one element for the function
% to behave correctly.
% 
% If you wish to use multiple criteria simultaneously, enter a cell array
% of strings containing those criteria you wish to use, or simply input the
% string 'all'. Note that because the actual calculations are not
% intensive (with exception of PRESS), all of the criteria excluding PRESS
% are calculated for all possible models. As a result, all of the criteria
% will be reported. However, with large numbers of parameters, sorting
% models on the basis of all criteria can be relatively intensive (though 
% not compared to the time required to estimate all the models), and as a 
% result models will only be sorted on the basis of the requested criteria.
% 
% Note that the SSE and R2 criteria produce the same results, as do aR2 and
% MSE. Consequently, if either is requested both statistics will be
% reported.
% 
% num_models: the maximum number of models to report on for each number of
% included parameters. If omitted, this value defaults to 2.
% 
% names: is a cell array of strings with the names of all of the variables
% in the model, starting with Y and progressing through the X variables in
% the same order they appear (from left to right) in the X matrix. IT IS
% HIGHLY RECOMMENDED THAT THIS VARIABLE BE USED, OR IT WILL BE CUMBERSOME
% TO DETERMINE WHICH PARAMETER IS WHICH FROM THE OUTPUT.
% 
% Note that this script will generate a divide by zero error when it
% estimates the model containing no parameters (intercept only).
% 

weighted=0;
known=0;
X=varargin{1};
Y=varargin{2};
criterion=varargin{3};
for k=4:nargin
    if iscell(varargin{k})
        if weighted
            known = 1;
        else
            names=varargin{k};
        end
    elseif isscalar(varargin{k})
        num_models = varargin{k};
    elseif isnumeric(varargin{k})
        W = varargin{k};
    else
        disp(['Error: invalid input: variable ' int2str(k)])
        beep
        return
    end
end

if ~exist('W','var')
    W=eye(size(X,1));
end

if ~isequal(X(:,1),ones(size(X,1),1))
    X=[ones(size(X,1),1) X];
end

if exist('num_models','var')==0
    num_models=2;
end

p=size(X,2);
model=zeros(2^(p-1)-1,p-1);
for k=0:2^(p-1)-1
    num=dec2bin(k);
    for j=1:size(num(:),1)
        model(k+1,end+1-j)=str2num(num(end+1-j));
    end
end

if exist('names','var')
    for k=1:size(model,1);
        params=find(model(k,:))+1;
        if ~known
            res(k)=mult_reg(X(:,[1 params]),Y,W,names);
        else
            res(k)=mult_reg(X(:,[1 params]),Y,W,'k',names);
        end
    end
else
    for k=1:size(model,1)
        params=find(model(k,:))+1;
        if ~known
            res(k)=mult_reg(X(:,[1 params]),Y,W);
        else
            res(k)=mult_reg(X(:,[1 params]),Y,W,'k');
        end
    end
    for k=1:size(model,1)
        res(k).model=model(k,:);
    end
end

k=1;
while k<=length(res)
    if isempty(res(k).n)
        res(k)=[];
        k=k-1;
    end
    k=k+1;
end

if strcmp(criterion{1},'all')
    criterion{1}='SSE';
    criterion{2}='aR2';
    criterion{3}='C';
    criterion{4}='AIC';
    criterion{5}='SBC';
    criterion{6}='PRESS';
end


for k=1:length(res)
    SSE(k)=res(k).SSE;
    R2(k)=res(k).R2;
    MSE(k)=res(k).MSE;
    aR2(k)=res(k).aR2;
    C(k)=res(k).SSE/res(end).MSE-(res(k).n-res(k).p);
    AIC(k)=res(k).n*log(res(k).SSE)-res(k).n*log(res(k).n)+2*res(k).p;
    SBC(k)=res(k).n*log(res(k).SSE)-res(k).n*log(res(k).n)+log(res(k).n)*res(k).p;
    p(k)=res(k).p;
end

for k=1:size(criterion(:),1)
    if strcmp(criterion{k},'PRESS')
        PRESS=zeros(1,size(model,1));
        for j=1:1:length(res)
            for i=1:res(j).n
                PRESS(j)=PRESS(j)+(res(j).e(i)/(1-res(j).H(i,i)))^2;
            end
        end
        break
    end
end

results=struct('criterion',{},'models',{});
for i=1:size(criterion(:),1)
    results(i).criterion=criterion{i};
    for k=1:size(X,2)
        mod=find(p==k);
        if strcmp(criterion{i},'SSE')|strcmp(criterion{i},'R2')
            [sorted,IX]=sort(SSE(mod));
        elseif strcmp(criterion{i},'MSE')|strcmp(criterion{i},'aR2')
            [sorted,IX]=sort(MSE(mod));
        elseif strcmp(criterion{i},'C')
            [sorted,IX]=sort(C(mod));
        elseif strcmp(criterion{i},'AIC')
            [sorted,IX]=sort(AIC(mod));
        elseif strcmp(criterion{i},'SBC')
            [sorted,IX]=sort(SBC(mod));
        elseif strcmp(criterion{i},'PRESS')
            [sorted,IX]=sort(PRESS(mod));
        end
        if size(mod(:),1)<=num_models
            for m=1:size(mod(:),1)
                results(i).models(end+1).SSE=SSE(mod(m));
                results(i).models(end).R2=R2(mod(m));
                results(i).models(end).MSE=MSE(mod(m));
                results(i).models(end).aR2=aR2(mod(m));
                results(i).models(end).C=C(mod(m));
                results(i).models(end).AIC=AIC(mod(m));
                results(i).models(end).SBC=SBC(mod(m));
                if exist('PRESS','var')
                    results(i).models(end).PRESS=PRESS(mod(m));
                end
                results(i).models(end).p=p(mod(m));
                results(i).models(end).model=res(mod(m));
                results(i).models(end).parameters=find(model(mod(m),:));
            end
        else
            for m=1:min([num_models size(mod(:),1)])
                results(i).models(end+1).SSE=SSE(mod(IX(m)));
                results(i).models(end).R2=R2(mod(IX(m)));
                results(i).models(end).MSE=MSE(mod(IX(m)));
                results(i).models(end).aR2=aR2(mod(IX(m)));
                results(i).models(end).C=C(mod(IX(m)));
                results(i).models(end).AIC=AIC(mod(IX(m)));
                results(i).models(end).SBC=SBC(mod(IX(m)));
                if exist('PRESS','var')
                    results(i).models(end).PRESS=PRESS(mod(IX(m)));
                end
                results(i).models(end).p=p(mod(IX(m)));
                results(i).models(end).model=res(mod(IX(m)));
                results(i).models(end).parameters=find(model(mod(IX(m)),:));
            end
        end
    end
    for k=1:size(results(i).models(:),1)
        if strcmp(criterion{i},'SSE')|strcmp(criterion{i},'R2')
            indicator(k)=results(i).models(k).SSE;
        elseif strcmp(criterion{i},'MSE')|strcmp(criterion{i},'aR2')
            indicator(k)=results(i).models(k).MSE;
        elseif strcmp(criterion{i},'C')
            indicator(k)=results(i).models(k).C;
        elseif strcmp(criterion{i},'AIC')
            indicator(k)=results(i).models(k).AIC;
        elseif strcmp(criterion{i},'SBC')
            indicator(k)=results(i).models(k).SBC;
        elseif strcmp(criterion{i},'PRESS')
            indicator(k)=results(i).models(k).PRESS;
        end
    end
    [sorted,IX]=sort(indicator);
    results(i).models=results(i).models(IX);
end



disp(' ')
disp(' ')
for k=1:size(results(:),1)
    disp(['Criterion: ' criterion{k}])
    if exist('PRESS','var')
        [s,err]=sprintf(['p:' '\t' 'SSE:' '\t' 'R2:' '\t' 'MSE:' '\t' 'C:' '\t' 'AIC:' '\t' 'SBC:' '\t' 'PRESS:' '\t' 'Parameters:']);
        disp(s)
    else
        [s,err]=sprintf(['p:' '\t' 'SSE:' '\t' 'R2:' '\t' 'MSE:' '\t' 'C:' '\t' 'AIC:' '\t' 'SBC:' '\t' 'Parameters:']);
        disp(s)
    end
    disp(' ')
    if exist('PRESS','var')
        for j=1:size(results(k).models(:),1)
            [s,err]=sprintf([num2str(results(k).models(j).p) '\t' num2str(results(k).models(j).SSE) '\t' num2str(results(k).models(j).R2) '\t'...
                num2str(results(k).models(j).MSE) '\t' num2str(results(k).models(j).C) '\t' num2str(results(k).models(j).AIC) '\t'...
                num2str(results(k).models(j).SBC) '\t' num2str(results(k).models(j).PRESS) '\t' num2str(results(k).models(j).parameters)]);
            disp(s)
        end
    else
        for j=1:size(results(k).models(:),1)
            [s,err]=sprintf([num2str(results(k).models(j).p) '\t' num2str(results(k).models(j).SSE) '\t' num2str(results(k).models(j).R2) '\t'...
                num2str(results(k).models(j).MSE) '\t' num2str(results(k).models(j).C) '\t' num2str(results(k).models(j).AIC) '\t'...
                num2str(results(k).models(j).SBC) '\t' num2str(results(k).models(j).parameters)]);
            disp(s)
        end
    end
    disp(' ')
    disp(' ')
    disp(' ')
end

