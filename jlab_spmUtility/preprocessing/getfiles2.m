function P2 = getfiles2(instr)

P = dir(instr);

P = str2mat(P(:).name);

if isempty(P), error('No files match.'), end

% move sub10 to end
P(end+1,:) = P(1,:);

P(1,:) = [];

mydir = fileparts(instr);

for i = 1:size(P,1)
    P2(i,:) = [mydir filesep P(i,:)];
end

P2(P2 == '/') = '\';

if strcmp(P2(1),'\')
    P2 = P2(:,2:end);
end

return


