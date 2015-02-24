function glm_compare(y,X,names,wh)
% glm_compare(y,X,names,wh preds)
%
% Tests a full GLM model (X) against a reduced
% model, without whpreds columns.
% 
% A change in R-square test is used to assess
% significance of the set of predictors.



% full model

[b,dev,stat]=glmfit(X,y); 

[B,BINT,R,RINT,STATS] = regress(y,[ones(size(X,1),1) X]);

dffull = stat.dfe; r2full = STATS(1);



% reduced model

X(:,wh) = [];
[b,dev,stat]=glmfit(X,y); 
[B,BINT,R,RINT,STATS] = regress(y,[ones(size(X,1),1) X]);

dfred = stat.dfe; r2red = STATS(1);

fprintf(1,'\n')
fprintf(1,'Omnibus test of decrement in R^2 without the following predictors:\n');
fprintf(1,'%s\n',names{wh});
fprintf(1,'\n')

% compare

compare_rsquare(r2full,r2red,dffull,dfred);


return