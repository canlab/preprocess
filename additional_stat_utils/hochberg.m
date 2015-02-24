function [out]=hochberg(in,fwe)

%
% Input must be a vector of p-values
% Output will be a vector of logicals where 1 indicates a significant
% result.
%


[sorted,ix]=sort(in,'descend');


t_out=zeros(size(in));
for k=1:size(sorted(:),1)
    if sorted(k)<=fwe/k
        t_out(k:end)=1;
        break
    end
end

out=zeros(size(in));
for k=1:size(sorted(:),1)
    out(ix(k))=t_out(k);
end


