
% replace irrelevant chars

s = strrep(str, '<', ' ');
s = strrep(s, '>', ' ');
s = strrep(s, ',', ' ');

% Split into cells and trim
s = strsplit(s);

s = cellfun(@strtrim, s, 'UniformOutput', 0);

% ID valid (?) email addresses
indx = cellfun(@(x) any(x == '@'), s);

s = s(indx);

% remove duplicates
s = unique(s', 'stable');

% rejoin
s = strjoin(s', ',');

% split into cell - alt output
ss = strsplit(s, ',')

fprintf('Unique addresses: %3.0f\n', length(ss));

char(ss)
