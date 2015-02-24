cd /data/vnl/plots

region = input('Enter area: ');
subj = input('Enter subject id: ');

title([region ' - ' subj])
saveas(gcf,[subj '_' region '_timeseries'],'fig')

% ============================
% * first plot 
% ============================
figure;subplot(2,2,1)
hold on
ROI.adjustedy = voistat('adjusty',ROI.ts.avg,HChoice,TR,numscans); 
[ROI.avgdata,ROI.avgste,ROI.ntrials] = voistat('trialavg',ROI.adjustedy,avgover,scolors,gca,vnl_condf,groups,basepoints,0,3);

title([subj '  ' region]);

 legend({'1','2','5','6','10','11'})





a = ROI.avgdata.group2 - ROI.avgdata.group1;
b = ROI.avgdata.group4 - ROI.avgdata.group3;
c = ROI.avgdata.group6 - ROI.avgdata.group5;

ROI.diff21 = a(3:end);
ROI.diff65 = b(11:end);
ROI.diff1110 = c(21:end);

% ============================
% * third plot 
% ============================
subplot(2,2,3)
hold on
plot(ROI.avgdata.group1,'k','LineWidth',2)
plot(a(3:end),'r','LineWidth',2)
plot(b(11:end),'g','LineWidth',2)
plot(c(21:end),'b','LineWidth',2)
legend({'first','2-1','6-5','11-10'})



% ============================
% filter these things.
% ============================
%cheby
%	Wn = [TR/HChoice .25];
% 	[B,A] = cheby2(12,20,Wn);

%for i = 1:6
%	eval(['d = ROI.avgdata.group' num2str(i) ';'])

	% moving average
	%for j = 2:length(d)-1		% avg of 3...
	%	d(j) = mean(d(j-1:j+1));
	%end
	
	%d = filter(B,A,d); 
	%eval(['d' num2str(i) ' = d;'])
%end


% ============================
% * second plot 
% ============================

subplot(2,2,2)

ROI.LPadjustedy = voistatHL('adjusty',ROI.ts.avg,HChoice,TR,numscans); 
[ROI.LPavgdata,ROI.LPavgste,ROI.ntrials] = voistat('trialavg',ROI.LPadjustedy,avgover,scolors,gca,vnl_condf,groups,basepoints,0,3);
a = ROI.LPavgdata.group2 - ROI.LPavgdata.group1;
b = ROI.LPavgdata.group4 - ROI.LPavgdata.group3;
c = ROI.LPavgdata.group6 - ROI.LPavgdata.group5;

%a = d2 - d1;
%b = d4 - d3;
%c = d6 - d5;

ROI.LPdiff21 = a(3:end);
ROI.LPdiff65 = b(11:end);
ROI.LPdiff1110 = a(21:end);



% ============================
% * second plot 
% ============================
%subplot(2,2,2)
%hold on
%plot(ROI.LPavgdata.group1,'r','LineWidth',2)
%plot(ROI.LPavgdata.group2,'r--','LineWidth',2)
%plot(ROI.LPavgdata.group3,'g','LineWidth',2)
%plot(ROI.LPavgdata.group4,'g--','LineWidth',2)
%plot(ROI.LPavgdata.group5,'b','LineWidth',2)
%plot(ROI.LPavgdata.group6,'b--','LineWidth',2)
 legend({'1','2','5','6','10','11'})


% ============================
% * fourth plot 
% ============================
subplot(2,2,4)
hold on
plot(ROI.LPavgdata.group1,'k','LineWidth',2)
plot(a(3:end),'r','LineWidth',2)
plot(b(11:end),'g','LineWidth',2)
plot(c(21:end),'b','LineWidth',2)
legend({'first','2-1','6-5','11-10'})

set(gcf,'Position',[298    80   598   579])		  % 360   137   737   797])

% ============================
% * saving and stuff 
% ============================
 saveas(gcf,[subj '_' region],'fig')
 saveas(gcf,[subj '_' region],'jpg')

name = ['roi_' subj '_' region];
eval([name ' = ROI;'])
eval(['save ' name ' ' name])
