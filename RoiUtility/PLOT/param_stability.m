function param_stability(OPT,clusters)
%  param_stability(OPT,clusters)
% Creates plots of parameter estimates given first session only, first two,
% etc. up to n sessions.  Created to assess stability of parameters as runs
% are added.
%
% timeseries field of clusters should be betas from beta imgs!
% OPT is input_stats_OPT struture for the study.
%
% Tor Wager, 2/13/02

for i = 1:length(clusters)
    
    figure; set(gcf,'Color','w')
    myname = clusters(i).name;
    myname(myname == '_') = ' ';
    
    subplot 221; hold on
    plot_scan_betas(clusters(i).timeseries,OPT,'params')
    title([myname ' avg betas'])
    
    subplot 222; hold on
    plot_scan_betas(clusters(i).all_data,OPT,'params')
    title(['Cl ' num2str(i) ' individual betas']),legend off
    
    subplot 223; hold on
    plot_scan_betas(clusters(i).timeseries,OPT,'diff')
    title(['Cl ' num2str(i) ' avg diffs']),legend off
        
    subplot 224; hold on
    plot_scan_betas(clusters(i).all_data,OPT,'diff')
    title(['Cl ' num2str(i) ' individual diffs']),legend off
    
end

return




% subfunctions
% ---------------------------------------------------------------

function plot_scan_betas(ts,OPT,method)

    %mycolors = {'r' 'b' 'g' 'y' 'k' 'c' 'm'};
    %mynames = unique(OPT.Cname); don't ever do this.
   
    mynames = OPT.Cname(1:(size(ts,1)-OPT.nruns) ./ OPT.nruns);
    
    hold on
    
    for i = 1:size(ts,2)
        
        myts    = ts(:,i);
        mytsk   = myts(1:end-OPT.nruns);
        myint   = myts(end-OPT.nruns + 1:end);
        
        mytsk   = reshape(mytsk,length(mynames),OPT.nruns)';
     
    
        switch method
        case 'params'
            % plot parameters by session
        
            plot(mytsk,'o-','LineWidth',2)
            legend(mynames)
            xlabel('Run number')
            ylabel('Parameter Estimate')
            set(gca,'XTick',1:OPT.nruns)
            set(gca,'XLim',[.5 OPT.nruns + .5])
            plot([get(gca,'XLim')],[0 0],'k')
        
        case 'diff'
            % plot running avg of difference between params for
            % runs 1:n-1 and runs 1:n
        
            myavg(1,:) = mytsk(1,:);
            for j = 2:size(mytsk,1)
                myavg(j,:) = mean(mytsk(1:j,:));
            end
            
            myavg = diff(myavg);
            myavg = [mytsk(1,:); myavg];
            
            plot(myavg,'o-','LineWidth',2)
            legend(mynames)
            xlabel('Number of Runs')
            ylabel('Parameter change due to adding this run')
            set(gca,'XTick',1:OPT.nruns)
            set(gca,'XLim',[.5 OPT.nruns + .5])
            plot([get(gca,'XLim')],[0 0],'k')
            
        end
            
    end
    
return