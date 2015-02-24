function EXPT = plot_dx_hrfs(EXPT,clusters)
% EXPT = plot_dx_hrfs(EXPT,clusters)
%
% uses EXPT.DX and clusters
%
% If not found, creates:
% EXPT.DX.regsofinterest = trial types
% EXPT.DX.mcol = colors
%

if ~isfield(EXPT.DX,'mcol')
    disp(EXPT.DX.dxnames)
    wh = input(['Enter vector of conditions to plot: ']);
    mcol = input('Enter vector of colors, e.g. {''ro-'' ''gd'' etc.}: ');
    EXPT.DX.regsofinterest = wh;
    EXPT.DX.mcol = mcol;
else
    wh = EXPT.DX.regsofinterest;
    mcol = EXPT.DX.mcol;
end
    
    
    for j = 1:length(clusters)
        if ~isempty(wh)
            figure('Color','w'); subplot(1,2,1)
            hold on;set(gca,'FontSize',18),xlabel(['Time (s)'])
            
            for k = 1:length(wh)
                %eval(['tmp = clusters(j).HRF.HRF' wh{k} ';'])
                str = (['tmp = 100*clusters(j).HRF.HRF{' num2str(wh(k)) '};']); 
                eval(str)
                
                x = 0:EXPT.DX.TR:length(tmp)*EXPT.DX.TR-1;
                h(k) = plot(x,tmp,mcol{k},'LineWidth',2);
                
                ww = find(tmp==max(tmp));
                tor_bar_steplot(tmp(ww),100*mean(clusters(j).HRF.STE{wh(k)}),mcol(wh(k)),x(ww)-1);
 
            end
            title(num2str(round(clusters(j).mm_center)))
            %legend(wh)
            legend(h,EXPT.DX.dxnames(wh));
            
            subplot(1,2,2)
            h = bar(clusters(j).HRF.BASE); set(h,'FaceColor',[.8 .8 .8]);
            tor_bar_steplot(clusters(j).HRF.BASE,clusters(j).HRF.BASESTE,{'k'});
            title('Baseline estimates')
            set(gcf,'Position',[109         215        1099         501])
        end
    end
    
    
    return
    
    