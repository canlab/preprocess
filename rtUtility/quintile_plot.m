function quintile_plot(rt3,inc,com,names)
% function quintile_plot(rt3,inc,com,names)
%
% plots quintiles of the cumulative prob. distribution
% comparing two conditions, indicated by indicator vectors inc and com (one for
% each condition).
%
% rt3 is a subjects x trials matrix for all trials.
% inc indicates columns of rt3 belonging to trial type 1
% com is the same, for trial type 2
% names is for the plot legend

% eliminate trials with NaNs for all subjects
% not necessary - redundant w/prctile operation
%whom = find(all(isnan(rt3),1));
%rt3(:,whom) = [];    % eliminate nan trials
%inc(whom) = [];
%com(whom) = [];


for i = 1:size(rt3,1),
    tmp = rt3(i,inc);
    ii = prctile(tmp,[20 40 60 80 95]);
    idat(i,:) = ii;
    tmp = rt3(i,com);
    ii = prctile(tmp,[20 40 60 80 95]);
    cdat(i,:) = ii;
end

plot(mean(log(idat)),'r','LineWidth',2); hold on; 
plot(mean(log(cdat)),'b','LineWidth',2)

legend(names)

err = ste(log(idat)-log(cdat));
tor_bar_steplot(mean(log(idat)),err,{'r'})
tor_bar_steplot(mean(log(cdat)),err,{'b'})
tor_bar_steplot(mean(log(cdat)),-err,{'b'})
tor_bar_steplot(mean(log(idat)),-err,{'r'})

xlabel('Percentile'); set(gca,'XTick',[1:5],'XTickLabel',{'20%' '40%' '60%' '80%' '95%'})
ylabel('Log RT');