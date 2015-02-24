function cl = contiguous_coordinates(xyzo,k)
% function cl = contiguous_coordinates(xyzo,k)
% spm_clusters.m kind of thing, except finds
% sets of clusters within k mm
% xyzo should be 3-column list of mm coordinates

% get contiguous sets
tmp = dist(xyzo') <= k;

% seed clusters

indx = 1; cl = zeros(size(xyzo,1),1);

while any(cl == 0)
    
    tmp3 = find(cl==0); first = tmp3(1);
    

    cl(first) = indx; cl(find(tmp(:,first))) = indx;

    msum = 0;
    while sum(cl==indx) > msum
        msum = sum(cl == indx);
        tmp2 = tmp(cl == indx,:); 
        wh = sum(tmp2,1) > 0;
        cl(find(wh)) = indx;
    end
    
    indx = indx + 1;
    %fprintf(1,'it %3.0f',indx)
    %cl'
end

return
