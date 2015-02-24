function [dx2,con2] = apply_dx_contrast(dx,numframes,con)
% [dx2,con2] = apply_dx_contrast(dx,numframes,con)
%
% for an FIR deconv matrix, a contrast set, and a specified number of conditions for each
% event type, builds and applies a matrix that collapses across event types
% by applying the contrast.
%
% e.g., con = [1 0 1 0; 0 1 0 1] and dx is an fir matrix
% numframes = [30 30 30 30]
% then this function will collapse the first and third sets of 30 columns
% into one set, and the 2nd and 4th sets of 30 into a second set.
%
% con2 returns an expanded contrast that gives a vector of which things to plot 
% in which series in plot_cluster_fir.  % this may need work to get right
% still.
%
% ONLY WORKS if numframes are all equal for contrasts to be combined.
%
% tor wager
% 

dx2 = [];
for i = 1:size(con,1)   % for each contrast entered
    dx2 = [dx2 expand_this_contrast(dx,numframes,con(i,:))];
end

dx2(:,end+1) = 1;   % add intercept

con2 = expand_contrast(ones(1,size(con,1)),repmat(numframes(1),1,size(con,1)),dx2);

return





function dx2 = expand_this_contrast(dx,numframes,con)
    
chk = numframes(find(con)); % where con values not 0
if ~all(chk == mean(chk)), error('numframes must be equal for basis functions with positive contrast weights to combine.'),end
chk = mean(chk);

c2 = []; % overall contrast matrix to apply
for i = 1:length(con),
    c = eye(numframes(i)) .* con(i);
    
    if con(i) == 0 & size(c,1)>chk, c = c(1:chk,:);,end   % if zeros w/diff elements, cut down
    if con(i) == 0 & size(c,1)<chk, c = [c; zeros(chk-size(c,1),size(c,2))];,end   % ...or pad
    
    c2 = [c2 c];
end
c2(:,end+1) = 0;    % add intercept
c2 = c2';   
dx2 = dx * c2;

return