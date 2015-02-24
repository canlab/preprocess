function [delta, wd, wb] = DX_find_delta(DX)
% [delta, wd] = DX_find_delta(DX)
%
% delta is matrix of stimulus onsets
% wd is indices of onset columns
% wb is columns between onsets
%
% Tor Wager
%
% for shape-free model matrix DX, finds columns
% coding onsets of new conditions.
% returns delta functions of all onsets and 
% vector of their column positions in DX.
%
% doesn't work for overlapping estimates

wd(1) = 1;
wb = [];

for i = 1:size(DX,2) - 2    % do not do this for intercept
    
    if sum(DX(:,i+1)) == 0
        % next column is empty; not an onset col, but maybe an offset col
        wd(i+1) = 0;
        if sum(DX(:,i)) ~= 0, wb(i) = 1;,end    % yes, an offset col
        
    elseif sum(DX(:,i)) == 0
        % this column is empty, next is not - an onset col
        wd(i+1) = 1;
    
    else
        a = shiftdown(DX(:,i));
        dif = a - DX(:,i+1);
        if any(dif)
            wd(i+1) = 1;    % next col is onset
            wb(i) = 1;      % this one is offset
        else
            wd(i+1) = 0;
            wb(i) = 0;
        end
        
        %elseif length(find(DX(:,i+1))) ~= length(find(DX(:,i)))
        %warning(['Col. ' num2str(i) ' and ' num2str(i+1) ' have different numbers of onsets.'])
        %wd(i+1) = 1;
        %wb(i) = 1;
        
        %else    
        % if next column is not just a shifted-over version of this one, then it's an onset
        %a = find(DX(:,i)) - find(DX(:,i+1));
        %if any(a ~= -1), wd(i+1) = 1;, wb(i) = 1;, else wd(i+1) = 0;, end
        %end
    end
    
end

wd = find(wd);
delta = DX(:,wd);
wb = find(wb);

%wb = [diff(wd) - 1 size(DX,2) - wd(end) - 1];
wb = [wb size(DX,2) - 1];
wb = wb - wd;
return


function out = shiftdown(in)

out(1) = 0;
out = [out; in];
out = out(1:length(in));

return