function rt = trimmean3(rtatt)
% rt = trimmean3(rt)
%
% returns mean of vector trimmed to 3 sd
% ignores nans

rtatt(isnan(rtatt)) = []; 
rs = 3*std(rtatt); rtatt(rtatt>mean(rtatt)+rs | rtatt<mean(rtatt)-rs) = [];
rt = mean(rtatt);

return
