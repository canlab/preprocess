function [y,h] = hanning_taper(w,y)
%
% inputs:
% w, width of hanning window in elements OR full tapering function
% y, data vector (column) or data matrix (columns are data vectors)

if isempty(w), h = [];, return, end


if length(w) == 1

    h=hanning(w*2);

    h_start=h(1:w);h_end=flipud(h_start);

    h = [h_start; ones(size(y,1)-w*2,1); h_end];
    
    h = repmat(h,1,size(y,2));
    
else
    h = w;
    
end


y = y .* h;


return
