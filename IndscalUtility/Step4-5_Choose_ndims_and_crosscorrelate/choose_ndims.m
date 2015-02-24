function [ndims,stats] = choose_ndims(resids)
% function [ndims,stats] = choose_ndims(DATA)
% 
% choose number of dimensions for data reduction

niter = 1000;    % iterations

numsub=length(resids);    %number of subjects

%m = celldat2matrix(DATA.r);
fprintf(1,'Choosing dimension for subject . ');

for i = 1:numsub  %size(m,3)     % for each subject
    
    fprintf(1,'%3.0f ',i);
    
    [pc,stats] = pca_npm(resids{i},niter,0,0,0);    %squeeze(m(:,:,i)),niter,0,0,0);
    drawnow
    
    nd(i) = size(pc,2);
    
end

fprintf(1,'\n');

ind_stats_ex = stats; 
%clear stats

stats.nd = nd;
%ndims = ceil(mean(nd));
ndims = round(mean(nd));

return

