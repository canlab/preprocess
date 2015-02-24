function [fitness,covmtx,ste_params,best_params] = prospect_ste_mc(x,y,p,ureport,iter)
%iter = 10;
noisestd = .2 * std(ureport);

best_params = zeros(iter,4);
noise = noisestd .* randn(size(ureport,1),iter);

%tic

for i = 1:iter
    
    repdat = ureport + noise(:,i);
    
    best_params(i,:) = prospect_fit_data(x,y,p,repdat,0);
end
ste_params = std(best_params);

fitness = 1 ./ sum(ste_params);

covmtx = [];
%toc
% % covmtx = cov(best_params);
% % p = size(best_params,2);
% % 
% % if rank(covmtx) < p
% %     % matrix is singular; ineligible design
% %     fitness = 0;
% % else
% %     fitness = 1./det(covmtx);
% % end



return