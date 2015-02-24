function o =  match_rows(a,b)
% function o = match_rows(a,b)
% Tor Wager
% given two matrices of equal no. of columns
% and diff (or same) no. rows,
% find indices of b that match rows of a
% 
% i.e., find which elements in b
% match the rows of a, as in matching lists
% of coordinates.
% 
% o is vector of indices

for i = 1:size(a,1)
    
    a2 = repmat(a(i,:),size(b,1),1);
    o(i) = find(~any(a2-b,2));
    
end

return
