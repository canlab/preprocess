function plot_cl_summary(clusters,cl_summary,O,varargin)

if length(varargin) > 0, beta_clusters = varargin{1};,end

% ----------------------------------------------------------------
% process what to plot and set up plotting
% ----------------------------------------------------------------
O.plot = 1;
        nplots = 0;
        if O.plotSelAvg, nplots = nplots + 1;,end
        if O.plotCustomBeta, nplots = nplots + 1;,end
        if O.plotDeconv, nplots = nplots + 1;,end
        if O.plotFits, nplots = nplots + 1;,end
        a = 'generic_figure';

a1 = O.plotSelAvg;
a2 = O.plotCustomBeta;
a3 = O.plotDeconv;
a4 = O.plotFits;








for j = 1:length(clusters)
    

        
        figure; set(gcf,'Color','w')
        set(gcf,'Position',[36         563        1415         448])
        %set(gcf,'Position',[36.0000   97.0000  874.0000  914.0000])
        
        plotSelAvg = a1;
        plotCustomBeta = a2;
        plotDeconv = a3;
        plotindex = 1;
        
        O.title = ['Group average: con ' clusters(j).title ' cl ' num2str(j)];
        O.title(O.title == '_') = ' ';
        a = O.title(16:end); a(a == ' ') = '_';
        
        for i = 1:nplots

            if plotSelAvg
                plotSelAvg = 0;
		        subplot(1,nplots,plotindex);hold on
                plotindex = plotindex + 1;
        
		        tor_plot_avgs(cl_summary(j).sel_avg,cl_summary(j).sel_avg_ste,O)
 		
            elseif plotDeconv
                plotDeconv = 0;
		        subplot(1,nplots,plotindex);hold on;grid on
                plotindex = plotindex + 1;
                
		        index = 1;
		        for k = 1:O.tbefore+O.tp:length(clusters(j).DXb) - 1
		        	plot(cl_summary(j).DXb_avg(k:O.tbefore+k+O.tp-1),O.colors{index},'LineWidth',2)
                    DXstecol{1} = O.colors{index};
                    tor_bar_steplot(cl_summary(j).DXb_avg(k:O.tbefore+k+O.tp-1),cl_summary(j).DXb_ste(k:O.tbefore+k+O.tp-1),DXstecol)
		            index = index + 1;
                end
		        title(['Deconvolution grand average: ' O.title(16:end)])
                
            elseif plotCustomBeta
                plotCustomBeta = 0;

                if exist('beta_clusters') == 1
                    subplot(1,nplots,plotindex);hold on
                    plotindex = plotindex + 1;
                    mycbavg = cl_summary(j).custombeta_avg3; % - cl_summary(j).custombeta_avg3(1);
                    mycbavg = mycbavg(find(O.adjustOfInterest(1:length(mycbavg)))); % only plot those of interest
                    bar(mycbavg)
                    mycberr = cl_summary(j).custombeta_ste3(find(O.adjustOfInterest(1:length(cl_summary(j).custombeta_ste3))));
                    tor_bar_steplot(mycbavg,mycberr,[])
                    set(gca,'XTick',1:length(O.adjNames));
                    set(gca,'XTickLabel',(O.adjNames));
                    %set(gca,'XTickLabels',custombeta_titles)
                    ylabel('Beta values (arbitrary units)')
                    title('b from beta*img files')
                else
                    subplot(1,nplots,plotindex);hold on
                    plotindex = plotindex + 1;
                
                    bar(cl_summary(j).custombeta_avg)
                    tor_bar_steplot(cl_summary(j).custombeta_avg,cl_summary(j).custombeta_ste,[])
                    set(gca,'XTick',1:length(O.adjNames));
                    set(gca,'XTickLabel',(O.adjNames));
                    title('Betas estimated from adjusted only')
                    %set(gca,'XTickLabels',custombeta_titles)
                end
                    
            
            
            elseif O.plotFits
                plotFits = 0;
                subplot(1,nplots,plotindex);hold on
                plotindex = plotindex + 1;
                
                hrf = spm_hrf(O.TR); hrf = hrf ./ max(hrf);
                
                b = cl_summary(j).custombeta_avg3(find(O.adjustOfInterest(1:length(cl_summary(j).custombeta_avg3))));
                allse = cl_summary(j).custombeta_ste3(find(O.adjustOfInterest(1:length(cl_summary(j).custombeta_avg3))));
                
                clear fits,clear se
                for myb = 1:length(b), 
                    fits{myb} = (hrf .* b(myb))';,
                    se{myb} = (zeros(length(hrf)))';
                    se{myb}(hrf == max(hrf)) = allse(myb);
                end
                
                O.legend = O.adjNames;
                tor_plot_avgs(fits,se,O)
             end   
        end
        
        
        % equalize axes
        if a1 & a3
            subplot(1,nplots,1); r(1,:) = get(gca,'YLim');
            subplot(1,nplots,2); r(2,:) = get(gca,'YLim');
            %subplot(1,4,3); r(3,:) = get(gca,'YLim');
            ymin = min(r(:,1));
            ymax = max(r(:,2));
            %ymin = -.2;
            %ymax = .7;
            subplot(1,nplots,1); set(gca,'YLim',[ymin ymax]);
            ylabel('Percent signal change')
            subplot(1,nplots,2); set(gca,'YLim',[ymin ymax]);
            %subplot(1,4,3); set(gca,'YLim',[ymin ymax]);
        end    
               

		saveas(gcf,a,'fig')
		saveas(gcf,a,'jpg')
		%print -dpsc2 -Pptr-schnabel01
		%pause(5);
		%close
		
end



