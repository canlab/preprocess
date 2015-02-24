% New trial average function
% Tor Wager, 6/17/01
%
% [A,ROI] = trialavg(S,ROI)
%
% S is a structure with the following fields.
% required fields:
% 'data': the timeseries data to average
%
% 'condf' or 'events': the condition function or timing of the event onsets
%	...or just include a number in condf to average over every n timepoints.
%
% 'window'
%	the number of timepoints to include in the averaging window ([-5;60])
%	always [start;end], 0 is first timepoint in trial
%	or a vector whose elements are the # points for each group of conditions.
%	e.g., [-10 -10 -10; 60 70 82] for 3 condition groups
%
% optional fields:
%
% 'color': cell array of plot colors
%
% 'groups': group the conditions in condf or events into these groups
%	before averaging.
%	example: [1 0 1 2 2 0] groups conditions 1 & 3 into group 1,
%		and conds 4 & 5 into group 2, and excludes the others.
% 'method'
%	options: 
%	'deconv' 			deconvolution
%	'average' (default) trial averaging
%	'continuous' 		continuous, binned averages
%
%
% 'tstrim'
%	...followed by st. dev. threshold.  this trims the whole timeseries.
%
%
% 'trim'
%	...followed by outlier trimming threshold in st. deviations (default 0)
%	this trims averages for each timepoint within cells.
%
%
% 'basepoints' (default [])
%	the baseline points to use as the zero point within each trial 
%	subtracts each trial's values from basepoints for that trial
%	if empty, no subtraction of each trial from baseline.
%	First image in trial = 0 in basepoints vector
%	Example: [-5:1] uses the first image and the 5 from the last trial as baseline.
% 
% 'basepoints2' (default [])
%	baseline points defined relative to trial n + 1
%	use this for, e.g., defining baseline as points right before onset of trial n AND 
%	those right before n+1 - beginning and end of trial.
%	warning! for this to work right, window(:,2) must have the trial length for each
%		group.  won't work right if averaging over PART of a trial.
%
% 'basesubtract'
%	'eachtrial': subtract each trial from its own baseline
%	'overall':	subtract average response from average baseline
%
% 'detrend'
%   linear detrending of average vector for each trial condition group
%
% plot options: just include these with a 1 or a 0 to indicate yes or no
%
%	'trialbytrial' trial by trial interactive plots default: 	(0)
%	'summary' 	   summary plot of averages						(1)
%	'ste'		   plot standard error bars						(1)
%	'figH'		   existing figure handle should be in field.   []
%
% 	'print': if exists as a field, verbose printing to printer name in the field.
%		works when trialbytrial is on, and prints to printer instead of prompting.
%
% 	'label': followed by a text label for summary figure
%
% ROI field is for backward compatibility: will insert averages in there, too.
%
% example:
%	S.data = y;
%	S.condf = condfunct;
%	S.color = {'r' 'g' 'b'};
%	S.groups = [1 1 2 2 3 3];
%	S.method = 'average';

function [A,ROI] = trialavg(S,ROI)

summaryplot = 1;
conditions = 1;
groups = 1;
   
% data vector
% ----------------------------------------
try 
	y = S.data;
	avgover = S.window;	
	stop = size(y,1);
catch
	S
 	error('Either no data field or averaging window in S structure.')
	
end


% condition function
% ----------------------------------------

try 
	condfunct = S.condf;	
catch

	try
		events = S.events;
	catch

		% for single fixed trial length
		% --------------------------------
       		condfunct = zeros(stop,1);
       		for i = 1:avgover:stop-avgover
           		condfunct(i) = 1;
       		end
	end	
end



% groups
% ----------------------------------------
try 
	groups = S.groups;
catch
	if exist('condfunct'), groups = 1:max(condfunct);
	elseif exist('events'), groups = 1:length(events);
	else error('Group assignment stage: No events or condf found.')
	end
end




% window
% ----------------------------------------
if length(avgover) == 1, avgover(2,1) = avgover(1); avgover(1,1) = 0;, end

if size(avgover,2) < max(groups), 
    % replicate window if not entered
    avgover = avgover(:,ones(1,max(groups)));
	%avgover(1,:) = avgover(1) * ones(1,max(groups));
	%avgover(2,:) = avgover(2) * ones(1,max(groups));
end
S.window = avgover;

for g = 1:max(groups)
	avglength(g) = avgover(2,g) - avgover(1,g) + 1;
	triallength(g) = avgover(2,g)+1;
end 



% method
% ----------------------------------------
try
	avgmethod = S.method;
catch
	avgmethod = 'average';
	S.method = avgmethod;
end

if strcmp(avgmethod,'average') & ~exist('condfunct')
	error('Must use condition function for average method.')
end


% basepoints
% ----------------------------------------
try
	basepoints = S.basepoints;
catch
	basepoints = [];
end
try
	basesubtract = S.basesubtract;
catch
	basesubtract = 'overall';
	S.basesubtract = basesubtract;
end
try
	basepoints2 = S.basepoints2;
catch
	basepoints2 = [];
end

% trimming of overall ts
% ----------------------------------------

try 
	tstrim = S.tstrim;
	if tstrim
		meany = mean(y); stdy = std(y);
   		exclude = (y > meany + stdy * tstrim | y < meany - stdy * tstrim);
		excluded = y(exclude);
		y(exclude) = Inf./Inf;
		S.data = y;

   		disp(['		trialavg.m: removed ' num2str(sum(exclude)) ' images at z = ' num2str(tstrim)])
		disp(['			intensity of removed: ' num2str(excluded')]);
	else
		disp(['		trialavg.m: no timeseries trimming.'])
	end
catch
	disp(['		trialavg.m: no timeseries trimming.'])
   
end



% trimming of avgs across trials
% ----------------------------------------
try
	outthresh = S.trim;
catch
	outthresh = 0;
end



% color
% ----------------------------------------
try 
	plotcolor = S.color;
	for g = 1:max(groups)
		if g > size(plotcolor,2),plotcolor{g} = 'k',end
	end
	S.color = plotcolor;
catch
	plotcolor = {'r' 'r--' 'g' 'g--' 'b' 'b--'};
	S.color = plotcolor;
end


% plotting options
% ----------------------------------------

try
	steplot = S.ste;
catch
	steplot = 0;
end

try
	trialbytrial = S.trialbytrial;
catch
	trialbytrial = 0;
	S.trialbytrial = 0;
end

try
	summaryplot = S.summary;
catch
	summary = 1;
end

try
	steplot = S.ste;
catch
	steplot = 0;
end

try
	voifigH = S.figH;
	if ishandle(voifigH) & length(voifigH) > 1, voifigplot = 1;,else voifigplot = 0;,end
catch
	voifigplot = 0;
	figure;voifigH = gca;
end


% print options
% ----------------------------------------
if isfield(S,'print')
	printstr = ['print -P' S.print];
	printthis = 1;
else printthis = 0;
end


A.options = S;

% detrend
% ----------------------------------------
if isfield(S,'detrend')
    dodetrend = S.detrend;
    if dodetrend, disp(['   * Linear Detrending of Averages for Each Group Requested.']),end
else
    dodetrend = 0;
end
if ~dodetrend, disp(['   - No Detrending of Averages for Each Group Requested.']),end

% get averages
% ----------------------------------------

switch avgmethod


case 'average'
% =======================================================================================================

   conditions = 1:max(condfunct);       % index conditions
   trialsskipped = 0;
   % loop through groups.  
   % for each group, figure out which conditions it goes with, get those, and average.
   for g = 1:max(groups)    

	avgd = [];
	index = 1;
	zeromagnitude = [];
        
	whichconds = conditions(groups == g);					% select trial type
        
	for j = 1:stop-avglength(g)								% loop thru all images j.
            if sum(condfunct(j) == whichconds) > 0             	% if current element matches any of the conditions
	
	    % doesn't include partial last trial w/o full time in trial
	    % for 1st trial with no time-before-trial info, pads with NaNs.

	    if j-1 >= abs(avgover(1,g))
		trialy = y(j+avgover(1,g):j+avgover(2,g))';
            else
		pad = zeros(1,abs(avgover(1,g))-(j-1)) * Inf / Inf;
		trialy = [pad y(j:j+avgover(2,g))'];
	    end

		if ~isempty(basepoints) | ~isempty(basepoints2)

			% Baseline specification and removal
			% -----------------------------------------------------------------------------
			if ~isempty(basepoints)
				mybase = j+basepoints;	 						% mybase is index of timepoints
			else
				mybase = [];
			end
			if ~isempty(basepoints2)
				mybase = [mybase j+triallength(g)+basepoints];
			end
			mybase = mybase(mybase > 0);					% make sure not to average off end
			mybase = y(mybase);								% get baseline response

			% Trim outliers from trial baseline estimate
			% -----------------------------------------------------------------------------
			if outthresh									
				try
					mybase = mybase(mybase > nanmean(mybase)-outthresh*nanstd(mybase) ...
					& mybase < nanmean(mybase)+outthresh*nanstd(mybase));
				catch
					disp(['Group ' num2str(g) ', event ' num2str(j) ': no baseline available, using first timepoint in trial.'])
					mybase = y(j);
				end
			end

			if sum(isnan(mybase)) == length(mybase)
				BASELINE = mybase'
				disp(['Group ' num2str(g) ', event ' num2str(j) ': no baseline available due to bad data.  Skipping trial.'])
			%	mybase = y(j);
			end

			warning off
			zeromagnitude(index) = mean(mybase);	% save signal at baseline or 1st time for conversion to % change later
			warning on								% basepoints go from 0 = 1st el. to avg over; -1 is right before 1st el
													% zeromagnitude(index) = mean(y(j+basepoints-1));
			
			if strcmp(basesubtract,'eachtrial')
	        		avgd(index,:) = trialy - zeromagnitude(index);	% change from baseline or 1st timepoint
			elseif strcmp(basesubtract,'overall')
				avgd(index,:) = trialy;
			else
				error('Unknown value for basesubtract.  Use eachtrial or overall.')
			end	


		else
			
			% no trial baseline
			% -----------------------------------------------------------------------------
			avgd(index,:) = trialy;
	
		end
			
		% trial by trial plot
		% -----------------------------------------------------------------------------
		if trialbytrial
			figure;
			BASELINE = mybase
			TRIAL_LENGTH = triallength(g)
			%TRIALY = trialy
			plot(avgover(1,g):avgover(2,g),trialy,'b-.'); 
			hold on; 
			plot(avgover(1,g):avgover(2,g),avgd(index,:),'r--.'); 
			if ~isempty(zeromagnitude)
				plot([1 length(trialy)],[zeromagnitude(index) zeromagnitude(index)],'k')	
			end
			title(['Group ' num2str(g) ' condition ' num2str(condfunct(j)) ' trial ' num2str(index)])
			legend({'Original trial data' 'Final trial data in avgd matrix'})

			if printthis
				pause(3)
				trialbytrial = 0;
				disp([	'trialavg.m: skipping remaining trials due to print option.'])
			else
				if index == 1,
					docont = input('Return to continue or q to quit plotting. ','s');
				else
					docont = input('','s');
				end
				if strcmp(docont,'q'),trialbytrial = 0;, end
			end
			close
		end

		if sum(isnan(avgd(index,:))) == length(avgd(index,:)), trialsskipped = trialsskipped + 1;, end
                index = index + 1;

            end % if match condition in group
        end % loop thru tpoints

	A.avgd{g} = avgd;
	if ~isempty(basepoints) & strcmp(basesubtract,'overall') 
		A.zeropoint(g) = nanmean(zeromagnitude);
	end
	eval(['ROI.ntrials.group' num2str(g) ' = avgd;'])

    end % loop thru groups

    disp(['		...voistat trialavg: averaging done, skipped ' num2str(trialsskipped) ' trials due to outlier values.'])
    A.trialsskipped = trialsskipped;
    

case 'deconv'
% =======================================================================================================


case 'continuous'
% =======================================================================================================


end % end switch method (get averages)






% -----------------------------------------------------------------------------------
% remove outliers from timepoint estimates
% and save the averages and stds
% -----------------------------------------------------------------------------------

trialbytrial = S.trialbytrial;
for g = 1:max(groups)

	avgd = A.avgd{g};

	if outthresh

    		numout = 0;
    		avg = zeros(1,avglength(g));
		myste = zeros(1,avglength(g));
		exclude = zeros(size(avgd,1),size(avgd,2));

		% --------------------------------------
		% * compute which are outliers
		% --------------------------------------
    		for tpoint = 1:avglength(g)
			vals = avgd(:,tpoint); 
			meanval = nanmean(vals); 
			cutoff = outthresh * nanstd(vals);

			exclude(vals > meanval+cutoff | vals < meanval - cutoff,tpoint) = 1;
		end

		% --------------------------------------
		% * trial by trial plot
		% --------------------------------------

		if trialbytrial
			figure;
			plot(avgd','-x'); 
			hold on;
			for i = 1:size(exclude,1)
				for j = 1:size(exclude,2)
					if exclude(i,j) == 1, plot(j,avgd(i,j),'ko','LineWidth',2);,end
				end
			end
			title(['Group ' num2str(g) ': Outlier removal'])

			if printthis
				eval(printstr)
			else
				if index == 1,
					docont = input('Return to continue or q to quit plotting. ','s');
				else
					docont = input('','s');
				end
				if strcmp(docont,'q'),trialbytrial = 0;, end
			end
			close
		end 

		% --------------------------------------
		% * NaN out the values excluded and save
		% --------------------------------------
		avgd(exclude==1) = NaN;		
		A.avgd{g} = avgd;
		eval(['ROI.ntrials.group' num2str(g) ' = avgd;'])

		% --------------------------------------
		% * Get the means
		% --------------------------------------
		avg = nanmean(avgd);
		myste = nanstd(avgd) ./ sqrt(sum(~exclude));
		
		numout = sum(sum(exclude));
    		if numout > 0, disp(['		trialavg: removed ' num2str(numout) ' outliers at z = ' num2str(outthresh) ', group ' num2str(g)]), end
 
		
			

   	else
		% --------------------------------------
		% * get avg, std, avgd w/o outlier removal
		% --------------------------------------
 		avg = nanmean(avgd);
   		for tpoint = 1:avgover
			vals = avgd(:,tpoint);
			myste(tpoint) = std(vals(~isnan(vals)))/ sqrt(length(vals(~isnan(vals))));
		end
	
		% --------------------------------------
		% * trial by trial plot
		% --------------------------------------
		trialbytrial = S.trialbytrial;
		if trialbytrial
			figure;
			plot(avgd','-x'); 
			hold on;
			title(['Group ' num2str(g) ': All trials'])
			if index == 1,
				docont = input('Return to continue or q to quit plotting. ','s');
			else
				docont = input('','s');
			end
			if strcmp(docont,'q'),trialbytrial = 0;, end
			close
		end 


	end % end if outliers...else

    if dodetrend
        xx = avg;
        ym = mean(xx);
        yy = detrend(xx);yy = yy+ym;
        avg = yy;
    end

	% ----------------------------------------------------------------------------
	% overall baseline subtraction
	% ----------------------------------------------------------------------------
	if sum(basepoints) & strcmp(basesubtract,'overall')
		avg = avg - A.zeropoint(g);
	end


	A.avg{g} = avg;
	if exist('myste') == 1, A.ste{g} = myste;, else, A.ste{g} = ste(avgd);, end
	eval(['ROI.avgdata.group' num2str(g) ' = avg;'])
   	eval(['ROI.avgste.group' num2str(g) ' = avg;'])
   	

end % loop thru groups	




% plotting
% ----------------------------------------------

if summaryplot

   	if voifigplot & outthresh

		% ------ plot the timeseries with removed images --------
		try
			axes(voifigH(2));cla
			plot(y,'b'); 
			hold on;
			x = 1:length(y);
			plot(x(~include),y(~include),'rx','LineWidth',2)
		catch
			disp('Can''t find voifig window - try InitWindows. excluded points not plotted.')	
		end
	
	end

	if voifigplot

		try
			axes(voifigH(4))
		catch
			disp('Can''t find voifig window - try InitWindows. Plotting on current axis.')
		end

	else

		try,axes(voifigH),catch,end

	end

	
	% ------ plot the selective average --------
      	
	for g = 1:max(groups)
		
		avg = A.avg{g};
		myste = A.ste{g};
      		if g == 1, cla, end	% clear axis if first group.
      		plot(avg,plotcolor{g})
      		hold on
      		wid = get(gca,'XLim');
      		wid = (wid(2) - wid(1))/50;
      		if steplot
      			for i = 1:avglength(g)
         			plot([i i],[avg(i)-myste(i) avg(i)+myste(i)],plotcolor{g})
         			plot([i-wid i+wid],[avg(i)-myste(i) avg(i)-myste(i)],plotcolor{g})
         			plot([i-wid i+wid],[avg(i)+myste(i) avg(i)+myste(i)],plotcolor{g})
      			end
      		end
   
	end % loop thru groups

	if isfield(S,'label'), title(S.label), end
	if printthis, eval(printstr), end

end % summaryplot

return
