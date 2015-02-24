function [results]=resid_analyze(X,Y,varargin)

% Usage:
%
% results=resid_analyze(X,Y,...)
%
% Calculates studentized deleted residuals, leverage values, DFFITS, Cook's
% Distance, DFBETAS, and the mean absolute percent difference (mpd) on all
% residuals for the multiple regression model specified by the design
% matrix X and the response variable Y. Other inputs may be included,
% matching the usage of mult_reg (see help mult_reg) to analyze
% non-standard regression models.
%
% In addition, the input 'P' will cause output to be printed to screen,
% while the input 'PF' followed by a valid filename will cause output to be
% printed to file. e.g. results=resid_analyze(X,Y,'PF','results.txt');
%
%

origin=0;
print=0;
for k=1:length(varargin)
    if ischar(varargin{k})
        if strcmp(varargin{k},'O')
            origin=1;
        elseif strcmp(varargin{k},'P')
            print=1;
            fid=1;
        elseif strcmp(varargin{k},'PF')
            print=1;
            fid=fopen(varargin{k+1},'w');
        end
    end
end

if (origin==0)&&~isequal(X(:,1),ones(size(X,1),1))
    X=[ones(size(X,1),1) X];
end

s='res=mult_reg(X,Y';
for k=1:length(varargin)
    s=[s ',varargin{' num2str(k) '}'];
end
s=[s ');'];
eval(s)

for k=1:size(X,1)
    results.student.t(k)=res.e(k)*((res.n-res.p-1)/(res.SSE*(1-res.H(k,k))-res.e(k)^2))^(1/2);
end

results.student.p=(1-tcdf(abs(results.student.t),res.n-res.p-1))*2;
results.student.outliers=find(results.student.p<=0.05);

results.leverage.values=diag(res.H);
results.leverage.outliers=find(results.leverage.values>=2*res.p/res.n);

results.DFFITS.values=results.student.t'.*(results.leverage.values./(1-results.leverage.values)).^(1/2);
results.DFFITS.outliers=find(results.DFFITS.values>=2*sqrt(res.p/res.n));

results.cooks.values=((res.e.^2)./res.p*res.MSE).*(results.leverage.values./((1-results.leverage.values).^2));
results.cooks.outliers=find(fcdf(results.cooks.values,res.p,res.n-res.p)>=.35);

for k=1:size(X,1)
    ares(k)=mult_reg(X([1:k-1 k+1:end],:),Y([1:k-1 k+1:end]));
    results.mpd(k)=sum(abs([ares(k).Yhat(1:k-1); X(k,:)*ares(k).b; ares(k).Yhat(k:end)]-res.Yhat)./res.Yhat)/res.n;
end

for i=1:size(X,2)
    for k=1:size(X,1)
        results.DFBETAS.values(k,i)=(res.b(i)-ares(k).b(i))/sqrt(ares(k).MSE*(res.var_cov_b(i,i)/res.MSE));
    end
end

for i=1:size(X,2)
    results.DFBETAS.beta(i).outliers=find(results.DFBETAS.values(:,i)>2/sqrt(res.n));
end

if print
    fprintf(fid,'Results of outlier analyses:\n\n');
    fprintf(fid,'Significant studentized deleted residuals:\n');
    if isempty(results.student.outliers)
        fprintf(fid,'None\n\n');
    else
        fprintf(fid,'Case #:\tp value:\n');
        for k=1:length(results.student.outliers)
            fprintf(fid,[num2str(results.student.outliers(k)) '\t' num2str(results.student.p(results.student.outliers(k))) '\n']);
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'Extreme leverage values:\n');
    if isempty(results.leverage.outliers)
        fprintf(fid,'None\n\n');
    else
        fprintf(fid,'Case #:\tleverage value:\n');
        for k=1:length(results.leverage.outliers)
            fprintf(fid,[num2str(results.leverage.outliers(k)) '\t' num2str(results.leverage.values(results.leverage.outliers(k))) '\n']);
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'Extreme DFFITS values:\n');
    if isempty(results.DFFITS.outliers)
        fprintf(fid,'None\n\n');
    else
        fprintf(fid,'Case #:\tDFFITS value:\n');
        for k=1:length(results.DFFITS.outliers)
            fprintf(fid,[num2str(results.DFFITS.outliers(k)) '\t' num2str(results.DFFITS.values(results.DFFITS.outliers(k))) '\n']);
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'Extreme Cook''s distance:\n');
    if isempty(results.cooks.outliers)
        fprintf(fid,'None\n\n');
    else
        fprintf(fid,'Case #:\tCook''s distance:\n');
        for k=1:length(results.cooks.outliers)
            fprintf(fid,[num2str(results.cooks.outliers(k)) '\t' num2str(results.cooks.values(results.cooks.outliers(k))) '\n']);
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'Extreme DFBETAS values:\n');
    for j=1:length(results.DFBETAS.beta)
        fprintf(fid,[res.param(j).name ':\n']);
        if isempty(results.DFBETAS.beta(j).outliers)
            fprintf(fid,'None\n\n');
        else
            fprintf(fid,'Case #:\tDFBETAS value:\n');
            for k=1:length(results.DFBETAS.beta(j).outliers)
                fprintf(fid,[num2str(results.DFBETAS.beta(j).outliers(k)) '\t' num2str(results.DFBETAS.values(results.DFBETAS.beta(j).outliers(k),j)) '\n']);
            end
            fprintf(fid,'\n');
        end
    end
end
if print&&fid~=1
    fclose(fid);
end