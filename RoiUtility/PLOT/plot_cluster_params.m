% ---------------------------------------------------------------------------
% * plot average of voxels by subject
% ---------------------------------------------------------------------------
clear col, clear leg, clear H
figure; hold on; set(gcf,'Color','w')

for i = 1:length(clusters)

	col = rand(1,3);
	
	plot(clusters(i).timeseries,'o','Color',col,'MarkerFaceColor',col)

	leg{i} = ['cluster ' num2str(i)];

end

ylabel('Parameter estimate','FontSize',14)
xlabel('Subject','FontSize',14)
legend(leg)
tit = clusters(1).title; tit(tit == '_') = ' ';
title(tit,'FontSize',14)

% ---------------------------------------------------------------------------
% * plot individual voxels by cluster
% ---------------------------------------------------------------------------
clear col, clear leg, clear H
figure; hold on; set(gcf,'Color','w')

for i = 1:length(clusters)

	col = rand(1,3);
	
	h = plot(clusters(i).all_data','.','Color',col,'MarkerFaceColor',col);

	H(i) = h(1);
	leg{i} = ['cluster ' num2str(i)];

end

ylabel('Parameter estimate','FontSize',14)
xlabel('Voxel','FontSize',14)
legend(H,leg)
tit = clusters(1).title; tit(tit == '_') = ' ';
title([tit ' all data'],'FontSize',14)


% ---------------------------------------------------------------------------
% * plot individual voxels by subject
% ---------------------------------------------------------------------------
clear col, clear leg, clear H
figure; hold on; set(gcf,'Color','w')

for i = 1:length(clusters)
	for j = 1:size(clusters(i).timeseries,1)
		if i == 1, col(j,:) = rand(1,3);, end			% define subject color if 1st cluster
	
		mycol = col(j,:);
		h = plot(clusters(i).all_data(j,:),'.','Color',mycol,'MarkerFaceColor',mycol);

		if i == 1, 
			H(j) = h(1);
			leg{j} = ['subject ' num2str(j)];			% save legend stuff if 1st cluster
		end
	end

end

ylabel('Parameter estimate','FontSize',14)
xlabel('Voxel','FontSize',14)
legend(H,leg)
tit = clusters(1).title; tit(tit == '_') = ' ';
title([tit ' all data'],'FontSize',14)



