function allrt = plot_rt_hists_stopsig(rt,ssoffset,stop_phat_ssrt,text)

% allrt = plot_rt_hists_stopsig(rt,ssoffset,stop_phat_ssrt,text)

n = length(rt);

allrt = cat(1,rt{:});

[h,x] = hist(allrt,50);

figure;

for i =1:n
    
    subplot(1,n,i); hold on;
    
    h = hist(rt{i},x);
    
    han = bar(x,h); set(han,'FaceColor',[.5 .5 .5],'EdgeColor',[.7 .7 .7]);
   
    myy = get(gca,'YLim');
    plot([ssoffset(i) ssoffset(i)] + stop_phat_ssrt,myy,'r','LineWidth',2);
    
    rt{i}(isnan(rt{i})) = [];
    stopperc = 100*sum(rt{i} > ssoffset(i) + stop_phat_ssrt) ./ length(rt{i});
    title([text num2str(stopperc) '%']);

end

return

