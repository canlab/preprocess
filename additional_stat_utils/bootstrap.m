function results=bootstrap(data,it,out_type)

% Usage:
% results=bootstrap(data,it,out_type)
%
% Returns each calculated value of the statistic specified in 'out_type'
% based on 'it' iterations of a bootstrap applied to 'data'. The first
% dimension of data is assumed to be subjects. Consequently, if data is
% an n x p array, results will be an it x p array. Similarly, if you input
% an n x p x q array, results will be an it x p x q array (and so on to
% higher dimensionality).
%
% Currently supported out_types are (more will be added as they are
% needed):
%
% 'mean'
% 'reg' for regression--data is interpreted as the first column being a
% criterion variable and all other columns being predictors. Output is a
% structure corresponding to the results of a call to 'mult_reg' (see 'help
% mult_reg' for details).

s=size(data);
n=s(1);

if strcmp(out_type,'mean')
    results=zeros([it s(2:end)]);
elseif strcmp(out_type,'reg')
    results(1)=mult_reg(data(:,2:end),data(:,1));
    results(it)=results(1);
end

for i=1:it
    if strcmp(out_type,'mean')
        r=ceil(rand(n,1)*n);
        d=data(r,:);
        
        results(i,:)=mean(d);
    elseif strcmp(out_type,'reg')
        for k=1:size(data,2)
            r=ceil(rand(n,1)*n);
            d(:,k)=data(r,k);
        end
        results(i)=mult_reg(d(:,2:end),d(:,1));
    end
end