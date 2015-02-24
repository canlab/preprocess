function cl_summary = batch_plot_clusters(Opt)
% cl_summary = batch_plot_clusters(Opt)
% compile output from biman_multisub_script running multisubject_cluster function
%
% Opt fields:
% plot		plot intermediate (individual subject) figures
% colors	cell array of plot colors
% legend	cell array of legend titles
% tp		required: time points in deconvolution
% Any fields for filterAdjust OK


Opt.plot = 0;
Opt.trimoverall = 0;
% Opt.trialbaseline = 1;
Opt.groups = [1 1 2 2 3 3 4 4];         
Opt.colors = {'r^-' 'bo-' 'gs-' 'mv-'};
Opt.legend = {'Incongruent' 'Congruent' 'Shape' 'Color'};
Opt.tp = 5;

tp = Opt.tp;

for mySub = [3 4 5 6 8 9 10]

	eval(['load sub' num2str(mySub) '_clusters;'])

	if Opt.plot
		figure;
		axisx = ceil(sqrt(length(clusters)));
		axisy = floor(sqrt(length(clusters)));
	end


	for j = 1:length(clusters)

		if Opt.plot,subplot(axisy,axisx,j),end

		% ----------------------------------------------------------------
		% get deconv hrfs for this cluster
		% ----------------------------------------------------------------
		% tp = 5;

		index = 1;
		for k = 1:tp:length(clusters(j).DXb) - 1
			cl_summary(j).hrf_ind{index}(mySub,:) = clusters(j).DXb(k:k+tp-1);
			index = index + 1;
		end
 	
		% ----------------------------------------------------------------
		% plot selective averages
		% ----------------------------------------------------------------
 		O.filtertype = 'none'; O.nscans = 6;O.scanadjust=1;O.y = clusters(j).timeseries;O.TR=2.5;O.percent=1;
 		y = filterAdjust(O);

		if j == 1
			Opt.title = ['Sub ' num2str(mySub) ' contrast ' clusters(j).title(end-3:end) ' cluster ' num2str(j) ' - tmax = ' num2str(max(clusters(j).Z))]; 
		else
			Opt.title = ['cluster' num2str(j)];
		end		

		myDelta = DX(:,[1:5:40]);
		
		[ROI,cluster(j).sel_avg_all_conds] = trialavg2(y,myDelta,[0 5],'options',Opt);
		cluster(j).sel_avg_ind = ROI.grpavg;		

		if Opt.plot
			legend off
			set(gcf,'Position',[ 224    75   870   871])
			set(gcf,'Color','w')
		end

		% ----------------------------------------------------------------
		% save group selective averages
		% ----------------------------------------------------------------
		for myCond = 1:length(cluster(j).sel_avg_ind)
			cl_summary(j).sel_avg_ind{myCond}(mySub,:) = cluster(j).sel_avg_ind{myCond}(1,:);
		end
		
		if Opt.plot,drawnow,end

	end

	if Opt.plot,close,end


	% ----------------------------------------------------------------
	% plot deconv hrfs for this cluster
	% ---------------------------------------------------------------
	if Opt.plot
		figure;
		axisx = ceil(sqrt(length(clusters)));
		axisy = floor(sqrt(length(clusters)));
	end

	for j = 1:length(clusters)

		if Opt.plot,subplot(axisy,axisx,j),end

		% ----------------------------------------------------------------
		% plot deconv hrfs for this cluster
		% ----------------------------------------------------------------

		if j == 1
			Opt.title = ['Deconvolution: Sub ' num2str(mySub) ' contrast ' clusters(j).title(end-3:end) ' cluster ' num2str(j) ' - tmax = ' num2str(max(clusters(j).Z))]; 
		else
			Opt.title = ['cluster' num2str(j)];
		end		

		if Opt.plot
			legend off
			set(gcf,'Position',[ 224    75   870   871])
			set(gcf,'Color','w')
			drawnow,end
		end

	end

	if Opt.plot,pause(3),close,end
end


for i = 1:length(cl_summary(1).sel_avg_ind)

	for j = 1:length(clusters)

		cl_summary(j).sel_avg_ind{i}( cl_summary(j).sel_avg_ind{i} == 0) = NaN;
		cl_summary(j).sel_avg{i} = nanmean(cl_summary(j).sel_avg_ind{i});
		cl_summary(j).sel_avg_ste{i} = nanstd(cl_summary(j).sel_avg_ind{i}) ./ sum(~isnan(cl_summary(j).sel_avg_ind{i}(:,1)));

		cl_summary(j).hrf_ind{i}( cl_summary(j).hrf_ind{i} == 0) = NaN;
		cl_summary(j).grp_hrf{i} = nanmean(cl_summary(j).hrf_ind{i});
		cl_summary(j).grp_hrf_ste{i} = nanstd(cl_summary(j).hrf_ind{i}) ./ sum(~isnan(cl_summary(j).hrf_ind{i}(:,1)));

	end

end

		
Opt.plot = 1;

for j = 1:length(clusters)
		figure
		Opt.title = ['Group average: contrast ' clusters(j).title(end-3:end) ' cluster' num2str(j)];
		tor_plot_avgs(cl_summary(j).sel_avg,cl_summary(j).sel_avg_ste,Opt)
 		
		a = Opt.title(16:end); a(a == ' ') = '_';
		%saveas(gcf,a,'fig')
		saveas(gcf,a,'jpg')
		%print -dpsc2 -Pptr-schnabel01
		%pause(5);
		close
		
end


return