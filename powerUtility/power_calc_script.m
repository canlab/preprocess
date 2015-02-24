%power_calc_script

clear a

figure;
set(gcf,'Color','w'); hold on;
sindx = 1;
trange = 3:15;

for SNR = 1:.5:4
 
    indx = 1;
    
    for tcrit = trange
        
        a(indx) = power_calc(SNR,tcrit, 20);
        indx = indx + 1;
        
    end
    
    plot(trange,a,'Color',rand(1,3),'LineWidth',2)
    myleg{sindx} = ['SNR = ' num2str(SNR)];
    sindx = sindx + 1;
    
end

legend(myleg)
grid on
xlabel('Critical t value for signficance','FontSize',14)
ylabel('Number of subjects required','FontSize',14)
title('Power as a function of threshold and SNR','FontSize',14)

plot([min(trange) max(trange)],[10 10],'k--')
plot([min(trange) max(trange)],[20 20],'k')