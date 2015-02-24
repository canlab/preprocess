function plot_fitted_vs_xtracted_betas(cl_summary,beta_clusters,O)

for j = 1:length(cl_summary)
    
    figure; set(gcf,'Color','w')
    
    subplot 221; hold on; 
    mk_bar(cl_summary(j).custombeta_avg2,cl_summary(j).custombeta_ste2,O.adjNames)
    title('b fitted from extracted ts')
    
    subplot 222; hold on; 
    mk_bar(cl_summary(j).custombeta_avg3,cl_summary(j).custombeta_ste3,O.adjNames)
    title('b from beta*img files')
    
    subplot 223; hold on; 
    mk_bar(cl_summary(j).custombeta_avg4,cl_summary(j).custombeta_ste4,O.adjNames)
    title('b fitted from extracted ts')
    
    subplot 224; hold on; 
    b_sub = [];
    index = 1;
    for i = O.numsubs
        b = beta_clusters{i}(j).all_data(:,1);
        b = b(find(O.adjustfor));
        b_sub(index,:) = mean(reshape(b,length(b)/O.nruns,O.nruns)'); 
        index = index + 1;
    end
    bar(mean(b_sub))
    plot(b_sub','ko')
    title('1st voxel, beta*img files')
                    
                    
end
                
                    
  


% sub-functions
                    
function mk_bar(avg,err,names)
                    
                    
                    %mycbavg = avg - avg(1);
                    %mycbavg = mycbavg(find(O.adjustOfInterest(1:length(mycbavg)))); % only plot those of interest
                    bar(avg)
                    %mycberr = err(find(O.adjustOfInterest(1:length(cl_summary(j).custombeta_ste3))));
                    
                    tor_bar_steplot(avg,err,[])
                    set(gca,'XTick',1:length(names));
                    set(gca,'XTickLabel',(names));
                    %set(gca,'XTickLabels',custombeta_titles)
                    ylabel('Beta values (arbitrary units)')
                    
return
                    