function plot_logit_fits(y,fit,onsets)
%plot_logit_fits(y,fit,onset indicator matrix)
%
%t1 = clock; [m,fit]=fit_expfun(y,1.5,c.trigs'); etime(clock,t1)
% plot_logit_fits(y,fit,c.trigs');

% get fitted timeseries

for i = 1:size(fit,2)
    
    f = conv(fit(:,i),onsets(:,i));
    f = f(1:length(y));
    
    fmat(:,i) = f;
    
end

tor_fig(2,1)
subplot(2,1,1)

plot(fit,'LineWidth',2)

tmp = get(gca,'Position'); tmp(3) = tmp(3) .* .4;
set(gca,'Position',tmp)

title('HRF estimates for each condition')


subplot(2,1,2)

plot(y,'LineWidth',1); 
hold on;
plot(sum(fmat,2),'k--','LineWidth',2);

title('Fitted timeseries')

return
