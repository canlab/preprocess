function [tmp,wh] = dominance_point_match(tmp,tmpm,tol)
% [tmp,wh] = dominance_point_match(tmp,tmpm,tol)
% Tor Wager
% 5/20/04
%
% Matches two lists of [n x 3] coordinates using a dominance metric
% With tolerance of tol 
%
% returns intersection points for first list only!
%
% wh is vector of which rows in tmp to save
%
%
% Used in cluster_mask.m for masking mm coordinates

if size(tmpm,1) > 1
    mmin = min(tmpm);
    mmax = max(tmpm);
else
    mmin = (tmpm);
    mmax = (tmpm);
end

    wh = zeros(size(tmp,1),1);    % which points to retain
    
    wh2 = ones(size(tmp,1),1);    % points to search thru point by point
    
    tmpm2 = [];                   % output
    
    % eliminate non-search-through pts
    % ------------------------------------
    
    % select cluster points that fall within min/max values in mask
    for j = 1:3
        wh2(tmp(:,j) < mmin(j) - tol | tmp(:,j) > mmax(j) + tol,:) = 0;
    end
    
    % select mask points that fall within min/max values of cluster
    for j = 1:3
        tmpm(tmpm(:,j) < min(tmp(:,j)) - tol | tmpm(:,j) > max(tmp(:,j)) + tol,:) = [];
    end
    
    % pt by pt search
    % ------------------------------------
    
    if ~isempty(tmpm) & ~sum(wh2) == 0
        for j = find(wh2')   % only subset of points that are w/i overall min.max
        
            % dominance metric
            a = [tmp(j,1) - tmpm(:,1) tmp(j,2) - tmpm(:,2) tmp(j,3) - tmpm(:,3)];
            a = max(abs(a)');   % max abs distance for cluster pt from each mask pt
            
            if any(a <= tol), 
                wh(j) = 1;, 
                %tmpm2 = [tmpm2; tmpm(find(a <= tol),:)];, % Oh, so slow!
            end
        
        end
    end  

    wh = find(wh);
    tmp = tmp(wh,:);
    
    return