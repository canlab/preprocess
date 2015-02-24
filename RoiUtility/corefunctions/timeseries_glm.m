function [betas,t,p,resid] = timeseries_glm(X,y,varargin)
% [betas,t,p,resid] = timeseries_glm(X,y,[names])
%
% runs glmfit on a series of columns of y
% intercept is added as first covariate automatically

if length(varargin) > 0
    % print table if short (< 10 cols)
    names = varargin{1};
    %names = [{'Intercept'}, names];
end
    
for i = 1:size(y,2)
    
    [b,dev,stats] = glmfit(X,y(:,i));
    
    betas(:,i) = b;
    
    t(:,i) = stats.t;
    
    p(:,i) = stats.p;
    
    
    if length(varargin) > 0
        warning off
        fprintf(1,'\nColumn %3.0f\n------------------------------------\n', i);
        regression_table(y(:,i),X,names,0);
    end
    
    if nargout > 3
        
        resid(:,i) = stats.resid;
        
    end
    
end



return






function printtable(b,stat,y,X,names)



    if size(y,2) < 10   % small, print table
        
            fprintf(1,'Parameter\tb-hat\tt\tp\tConf. interval\t\tPartial Corr\tRob. Partial Corr.\n')

            for i = 1:length(b)

                if i == 1, %intercept
                        fprintf(1,'%s\t%3.2f\t%3.2f\t%3.4f\t%3.2f\t%3.2f\t\n', ...
                    names{i},b(i),stat.t(i),stat.p(i),BINT(i,1),BINT(i,2));

                else
                    figure; subplot(1,2,1);[r,str,sig]=prplot(y,X,i-1,0); xlabel([names{i}])
                    subplot(1,2,2); [r2,str,sig]=prplot(y,X,i-1,1); xlabel([names{i} ': Robust IRLS'])

                    fprintf(1,'%s\t%3.2f\t%3.2f\t%3.4f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t\n', ...
                    names{i},b(i),stat.t(i),stat.p(i),BINT(i,1),BINT(i,2),r,r2);
                end
            end


            fprintf(1,'\nOverall R2 = %3.2f, F = %3.2f, p = %3.4f\t\n\n',STATS(1),STATS(2),STATS(3));

    end
    
return