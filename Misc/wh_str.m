% wh = wh_str(regexp, cellstrs)
% Searches a cellstr array for all entries that match a regular expression

function wh = wh_str(str, strs)
    if(~ischar(str))
        error('Parameter str is not a string.');
    elseif(~iscellstr(strs) && ~ischar(strs))
        error('Parameter strs is not a cellstr or char matrix.');
    end
    
    strs = cellstr(strs);
    
    wh = ~cellfun(@isempty, regexp(strs, str, 'match'));
end