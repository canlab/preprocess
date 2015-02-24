function plot_cluster_hrfs2(cluster,nsess)
% function plot_cluster_hrfs2(cluster,nsess)
%

clname = cluster.name;
clname(clname == '_') = ' ';

set(gcf,'Color','w')
colors = {'r' 'b' 'g' 'y' 'm' 'c' 'k'};

subplot(1,2,1)
hold on,grid on

nInSess = floor(size(cluster.b,1) ./ nsess);
disp(['Found ' num2str(nInSess) ' conditions per session.'])

myleg = [];
for i = 1:nInSess
	plot(cluster.b(i:nInSess:length(cluster.b)),colors{i},'LineWidth',2)
	myleg{i} = ['Condition ' num2str(i)];
end
title(['Model betas across sessions: ' clname])

subplot(1,2,2)
hold on,grid on

for i = 1:length(cluster.hrf_est)
	plot(cluster.hrf_est{i},'LineWidth',2,'Color',colors{i})
end

xlabel('Seconds')
ylabel('Percent signal change')
legend(myleg)

center = round(mean(cluster.XYZmm,2)');
center2 = [num2str(center(1)) '_' num2str(center(2)) '_' num2str(center(3))];
cltitle = cluster.title;
cltitle(cltitle == ' ') = '_';

title(['Cluster center: ' num2str(center)])
saveas(gcf,['cluster_' cltitle '_' center2],'fig')


return