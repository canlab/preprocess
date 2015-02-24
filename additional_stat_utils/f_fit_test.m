function [results]=f_fit_test(x,y,varargin)

% 
% Usage = results=f_fit_test(X,Y,...)
% 
% x and y are vectors of equal length containing the regression variables.
% Results contains the F and p values for the test of whether any
% association between x and y is linear. Other inputs may be included,
% matching the usage of mult_reg (see help mult_reg) to analyze
% non-standard regression models.
% 

s='res=mult_reg(x,y';
for k=1:length(varargin)
    s=[s ',varargin{' num2str(k) '}'];
end
s=[s ');'];


if size(x,2)==1
    eval(s)
    results.data=struct('x',{},'y',{});
    while isempty(x)==0
        find(x(1)==x);
        results.data(end+1).x=x(ans);
        results.data(end).y=y(ans);
        find(x(1)~=x);
        x=x(ans);
        y=y(ans);
    end

    results.SSPE=0;
    for k=1:size(results.data(:),1)
        for j=1:size(results.data(k).x(:),1)
            results.SSPE=results.SSPE+(results.data(k).y(j)-mean(results.data(k).y(:)))^2;
        end
    end

    results.SSE=res.MSE*(res.n-2);

    results.df1=size(results.data(:),1)-2;
    results.df2=res.n-size(results.data(:),1);

    results.F=((results.SSE-results.SSPE)/results.df1)/(results.SSPE/results.df2);
    results.p=1-fcdf(results.F,results.df1,results.df2);
else
    string='';
    for k=1:size(x,2)
        if k~=1
            string=[string '&'];
        end
        string=[string 'x(1,' int2str(k) ')==x(:,' int2str(k) ')'];
    end
    
    eval(s)
    
    results.data=struct('x',{},'y',{});
    while isempty(x)==0
        eval(['find(' string ');'])
        results.data(end+1).x=x(ans,:);
        results.data(end).y=y(ans);
        eval(['find((' string ')==0);'])
        x=x(ans,:);
        y=y(ans);
    end
    
    results.SSPE=0;
    for k=1:size(results.data(:),1)
        for j=1:size(results.data(k).y(:),1)
            results.SSPE=results.SSPE+(results.data(k).y(j)-mean(results.data(k).y(:)))^2;
        end
    end
    
    results.SSE=res.SSE;
    results.SSLF=results.SSE-results.SSPE;
    
    results.c=size(results.data(:),1);
    results.df1=results.c-res.p;
    results.df2=res.n-results.c;
    
    results.F=(results.SSLF/results.df1)/(results.SSPE/results.df2);
    results.p=1-fcdf(results.F,results.df1,results.df2);
end

