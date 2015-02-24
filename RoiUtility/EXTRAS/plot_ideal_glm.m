function plot_ideal_glm(meanest,min95est,max95est,ideal_data,meanfits,INFO,hrf,snr,TR)
%function plot_ideal_glm(meanest,min95est,max95est,ideal_data,meanfits,INFO,hrf,snr,TR)
% Tor Wager

%msdblim = max(max(abs(msdb)));
%bmlim = max(max(abs(biasmean))); bmlim = [-bmlim bmlim];

figure; subplot 221; hold on; set(gcf,'Color','w');
plot(meanfits','r-','LineWidth',.5)
plot(ideal_data,'k','LineWidth',2)
ylabel('Fitted values')
xlabel('Time')
title(['Est. vs. True response: ' INFO.condition_tested],'FontSize',14)
if max(get(gca,'XLim')) > 200, set(gca,'Xlim',[0 200]), end


subplot 222; hold on;
x = (1:length(hrf)) * TR;
plot(x, hrf,'k','LineWidth',2)
if ~isempty(INFO.autocorrelation)
    plot((1:length(INFO.autocorrelation))./(length(INFO.autocorrelation) / length(hrf)),INFO.autocorrelation,'k--','LineWidth',2)
    myleg = {'HRF','ACF'};
else
    myleg = {'HRF - white noise simulated'};
end

legend(myleg)
title('True response and ACF for simulation','FontSize',14)
xlabel('Time from trial onset')


subplot 223; hold on; 
ndesigns = size(min95est,1);
msdstd = [min95est(:,INFO.ttype) max95est(:,INFO.ttype)];
plot([(1:ndesigns); (1:ndesigns)], msdstd','k-','LineWidth',2)
    hold on;
    plot(meanest(:,INFO.ttype),'k.','LineWidth',4)
    xlabel('Design Realization')
ylabel('Mean parameter estimate')
title(['Mean and 95% CI for estimated betas for each design, ' INFO.condition_tested])


subplot 224; hold on; 
for j = 1:size(meanest,2)
    a = meanest(:,j);
    for i = 1:length(a)
        plot([j j],[a(i) a(i)],'k.','LineWidth',2)
    end
    myxl{j} = ['Beta' num2str(j)];
end
set(gca,'XTick',1:size(meanest,2))
ylabel('Mean Parameter Estimate')
title('Mean estimated betas for each trial type')


return

