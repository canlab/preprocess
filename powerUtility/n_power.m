n   = 20;           % no. of subjects
O.snr = [.1 .2 .3 .5 .7 1 1.5 2 2.5 3];         % signal to noise ratio values to test
O.verbose = 1;        % verbose output
O.dofull = 1;         % do all-data analysis
O.dosplit = 0;        % do split-data analysis

%fig = figure;
index = 2;

for n = 20:10:60
    
    O.n = n;
    [p,ps] = snr_power(O,fig);
    
    leg{index} = ['n = ' num2str(n)];
    index = index + 1;
    
    legend(leg)
    
    saveas(gcf,'n_power','fig')
end

