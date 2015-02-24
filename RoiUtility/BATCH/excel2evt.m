function evt = excel2evt(exmat)
% function evt = excel2evt(exmat)
%
% input: matrix of event onset times from excel in TRs
% First scan (TR) is 0
% columns are conditions
% blank fields filled with NaN
%
% output: cell array of event onset times
% One cell per condition
% Column vectors of times

for i = 1:size(exmat,2)
    a = exmat(:,i);
    a(isnan(a)) = [];
    evt{i} = a;
end

return