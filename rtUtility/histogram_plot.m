function histogram_plot(rt3,inc,com)
% function quintile_plot(rt3,inc,com)
%
% plots quintiles of the cumulative prob. distribution
% comparing two conditions, indicated by indicator vectors inc and com (one for
% each condition).
%
% rt3 is a subjects x trials matrix for all trials.
% inc indicates columns of rt3 belonging to trial type 1
% com is the same, for trial type 2
% names is for the plot legend



hold on;

rtz = rt3'; 
whom = find(any(isnan(rtz),2));
rtz(whom,:) = [];    % eliminate nan trials
rtz = scale(rtz)';

inc(whom) = [];
com(whom) = [];
inc = find(inc);
com = find(com);

rti = rtz(:,inc); rti = rti(:);
rtc = rtz(:,com); rtc = rtc(:);
[h,x] = hist(rti,100); 
han=bar(x,h./length(rti),'r'); 
set(han,'EdgeColor','none');
alpha .5; 
h2 = hist(rtc,x); han=bar(x,h2./length(rtc),'b');,
set(han,'EdgeColor','none');

xlabel('Z-score of RT'); ylabel('Probability density')

return
