function subcluster_montage(subc)
% subcluster_montage(subc)
% 
% plots each subcluster in a different random color
% subc is an ordinary clusters vector of structures
%
% tor wager
%
% 

mycols = {'1 0 0' '0 0 1' '0 1 0' '1 1 0' '0 1 1' '1 0 1'};

str = ['montage_clusters([]'];
colstr = ['{'];

for i = 1:length(subc)
    
    str = [str ',subc(' num2str(i) ')'];
    
    if i > length(mycols)
        colstr = [colstr ' [' num2str(rand(1,3)) ']'];
    else
        colstr = [colstr ' [' mycols{i} ']'];
    end
end

str = [str ',' colstr '});'];
eval(str)

return
