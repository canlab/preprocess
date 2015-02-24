% biman_deconv_cluster

mytp = 5;
myTR = 2.5;
myeres = 35;
printstring = ['print -dpsc2 -append VOI_' xY.name '.ps'];

% -------------------------------------------------------------------
% * before running: get timeseries of cluster from SPM
% -------------------------------------------------------------------
% need: xY structure in workspace

% -------------------------------------------------------------------
% * filter and adjust data - nothing done with this right now.
% -------------------------------------------------------------------

O.y = xY.Y;
O.doHP = 0;
O.scanadjust = 1;
O.nscans = 6;
O.filtertYe  = 'none';
O.TR = myTR;
O.percent = 1;
% xY.Yp = filterAdjust(O);	% no point, don't have original data mean.

% -------------------------------------------------------------------
% * do deconvolution
% -------------------------------------------------------------------

if exist('DX') == 1
	[dummy,hrfs] = tor_cluster_deconv('data',xY.Y,'tp',mytp,'DX',DX.DX);
	xY.b = hrfs.b;
	xY.hrf = hrfs.hrf;
else

	[DX,hrfs] = tor_cluster_deconv('data',xY.Y,'tp',mytp,'samprate',myeres,'columns',1:8,'Sess',Sess);
	xY.b = hrfs.b;
	xY.hrf = hrfs.hrf;
end



%-------------------------------------------------------------------
% * get current axis handle
% -------------------------------------------------------------------
axish = gca;

% -------------------------------------------------------------------
% * plot
% -------------------------------------------------------------------
	hrf = xY.hrf;


	mycolors = {'ro-' 'r^--' 'bo-' 'b^--' 'k^--' 'ko-' 'g^--' 'go-' };
	cla, 
	hold on, 
	grid on
	title('Hemodynamic responses')
	for i = 1:length(hrf)
		plot(hrf{i},mycolors{i},'LineWidth',2),hold on
	end
    myleg = {'Incongruent Left' 'Incongruent Right' 'Congruent Left' 'Congruent Right' 'Shape Right' 'Shape Left' 'Color Right' 'Color Left'};
	legend(myleg)

eval(printstring)
pause(1)
    
% -------------------------------------------------------------------
% * trial averaging
% -------------------------------------------------------------------    

axes(axish)

xY = trialavg2(xY.Y,DX.sf,[0 4],'trialbaseline','trimoverall','append',xY);
xY.options.title = [xY.name ' Selective averages'];
xY.options.title(xY.options.title == '_') = ' ';
xY.options.legend = myleg;
xY.options.plot = 1;
xY.options.colors = mycolors;
cla
tor_plot_avgs(xY.avg,xY.ste,xY.options)


eval(printstring)
%pause(3)
axes(axish)



% -------------------------------------------------------------------
% * one more deconv collapsing across response hand
% ------------------------------------------------------------------- 
xY2 = xY;

% make collapsed stick function at TR resolution
a = DX.sf;
sf2{1} = a{1} + a{2}; sf2{2} = a{3} + a{4}; sf2{3} = a{5} + a{6}; sf2{4} = a{7} + a{8};

[DX2,dummy] = tor_make_deconv_mtx(sf2,mytp,1);
[dummy,dummy,xY2.b,xY2.hrf] = tor_cluster_deconv('data',xY2.Y,'DX',DX2,'tp',mytp);


mycolors = {'ro-' 'b^-' 'ks-' 'gd-'};
	cla,hold on, grid on,title('Hemodynamic responses')
	for i = 1:length(xY2.hrf)
		plot(xY2.hrf{i},mycolors{i},'LineWidth',2),hold on
	end
    myleg = {'Incongruent' 'Congruent' 'Shape only' 'Color only'};
	legend(myleg)

eval(printstring)
pause(1)



% -------------------------------------------------------------------
% * one more trial avg collapsing across response hand
% -------------------------------------------------------------------   

xY2.options.legend = {'Incongruent' 'Congruent' 'Shape only' 'Color only'};
xY2.options.nosteplot = 1;
xY2.options.groups = [1 1 2 2 3 3 4 4];
xY2.options.colors = {'ro-' 'b^-' 'ks-' 'gd-'};
xY2 = trialavg2(xY2.Y,DX.sf,[0 4],'options',xY2.options,'append',xY2);

eval(printstring)
pause(1)


% -------------------------------------------------------------------
% * make ideal data - does deconvolution work?
% -------------------------------------------------------------------
[ideal_hrf, ideal_data] = ideal_deconv(DX.DX,DX.hires_sf,mytp,myTR,myeres,0);   % DX must be structure with DX deconv matrix and sf, stick function
xY.ideal_hrf = ideal_hrf;
xY.ideal_data = ideal_data;

% -------------------------------------------------------------------
% * extra plots of ideal data
% ------------------------------------------------------------------- 
axes(axish)
ideal_opt = xY.options;
ideal_opt.title = [ideal_opt.title ' Ideal Data'];
trialavg2(ideal_data,DX.sf,[0 4],'options',ideal_opt);

eval(printstring)
pause(1)

% plot ideal deconv
[ideal_hrf, ideal_data] = ideal_deconv(DX.DX,DX.hires_sf,mytp,myTR,myeres,1); 

% finish with selective avg plot showing
%axes(axish)
%cla
%tor_plot_avgs(xY.avg,xY.ste,xY.options)




% -------------------------------------------------------------------
% * save the cluster data
% -------------------------------------------------------------------

str = ['save VOI_' xY.name ' xY xY2 DX DX2'];
eval(str)



