function [varargout] = voistat(action,varargin)
% X = voistat('getX',ref,irf,[samprate in s])
%		gets a design matrix from a condition function (.ref file)
%
% [voifigH,avgdata] = voistat('voifig',voifigH,y,y2,avgover,filter len. in s,TR(samprate))
% 		makes special figure window with 3 plots, returns handles to all in voifigH
%		plots timeseries,fft,and trial avg (all optional) in appropriate windows.
%		if voifigH = [], creates new window; otherwise activates old one.
%		if y2 = [], doesn't plot that - optional.
%
% [avg,ste,trials] = voistat('trialavg',y,avgover,color,voifigH,condfunction,condgroups)
% cond function: for selective averaging
%    condition number at start of trial,zeros between start points
%    e.g. [1;0;0;0;3;0;0;0;2;0;0;0;4;0;0;0;1;0;0;0; etc...]
% condition group: row vector as long as the number of conditions, specifying which to plot together
%    e.g. [1 1 2 3] creates three groups, with condition 1 or 2 plotted as one line, 3 and 4 as two separate lines
%    color can be a cell array the same size as the number of groups.
%
% H = voistat('plotall',gls)
%		makes plots for several diagnostics given a gls structure
%
% X = voistat('figure')
%		makes a specially sized figure window
%
% voistat('fitvsy',gls,whichp)
%		plots fits vs y given a gls structure and a plot number
%
% voistat('fitvse',gls,whichp)
%		plots fits vs e given a gls structure and a plot number
%
% voistat('fits',gls,whichp)
%		index plot of residuals, given a gls structure and a plot number
%
% voistat('adjusty',y,[high pass freq. cutoff],[sampling rate in s])
%		adjusts y to filter out low frequencies of no interest.
%		converts to % signal change from 0.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SWITCHYARD													  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch(action)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'getX',
   % load the user specified file as timeseries                                                                                                                                                                                        
	%EditHandle = findobj(gcbf, 'Tag', 'condAfile');                                                                                                                                                                                                                     
	%infile = get(EditHandle, 'String');                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
	%eval(['timeseries = load (''' infile ''');'])
   
   % set up HRF and variables
   if nargin < 2,error('X = voistat(''getX'',ref,irf,samprate in s'),end
   ref = varargin{1};
   irf = varargin{2};
   samprate = varargin{3};
   
   if not(strcmp(irf,'spm') | strcmp(irf,'SPM'))						% use IRF from input
   	disp('building HRF...using custom function')
   	HRF = irf / max(irf);
      HRF = interp(HRF,1000*TR);										
   else																				% use SPM's gamma function
  	   disp('building HRF...using spm_hrf gamma function')
   	HRF = spm_hrf(samprate);
  		HRF = HRF / max(HRF);
	end

	% make delta function
	for i = 1:max(ref(:,1))
      deltaF(:,i) = ref(:,1) == i;
      X(:,i) = conv(deltaF(:,i),HRF);
   end
   
   % cut off extra values created by convolution
   X = X(1:size(deltaF,1),:);
   
   disp('adding intercept...')
   X(:,i+1) = 1;
   varargout{1} = X;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
case 'voifig'
% function [voifigH,avgdata,y2] = voifig(voifigH,y,y2,avgover,[filt length, TR],numscans [opt for filtering], basepoints,sebars)
% basepoints is for selective avg, sebars 0 or 1.
% ********************************************
% voifigH is 4 handles: figure, tsH,fftH,avgH
% if voifigH is a figure handle, activate.  if voifigH = [], create figure.
%
% optional arguments -
% y = timeseries to plot
% y2 = fits, or any other timeseries

disp('*	voistat voifig')

% === if filter option is specified, make y2 the filtered data. ======
if nargin > 4
	if nargin > 8, basepoints = varargin{8};,else basepoints = 1;,end
	if nargin > 9, sebars = varargin{8};,else sebars = 0;,end

	if nargin > 7
		numscans = varargin{7};
	else
		numscans = 1;
	end
	y = varargin{2};
   	avgover = varargin{4};
	if nargin > 5,
	   	HChoice = varargin{5};, TR = varargin{6};
   		y2 = voistat('adjusty',y,HChoice,TR,numscans);
		disp('		...y2 (red) is filtered data.')
	end
end

if nargin > 7
	numscans = varargin{7};
else
	numscans = 1;
end


% === make the figure, if necessary, or activate it =====
voifigH = varargin{1};
if voifigH > 0
   figure(voifigH(1))
   varargout{1} = voifigH;
else
	voifigH(1) = figure;
	set(voifigH(1),'Position',[452   214   794   420])
    figure(voifigH(1))
   voifigH(2) = axes('position',[.05 .65 .9 .3]);
   voifigH(3) = axes('position',[.05 .15 .4 .4]);
   voifigH(4) = axes('position',[.55 .15 .4 .4]);
   set(gcf,'Name','VOIstat timeseries viewer')
   set(gcf,'NumberTitle','off')
   set(gcf,'PaperOrientation','landscape')
   set(gcf,'Color',[.7 .7 .7])
   varargout{1} = voifigH;
end

if nargin > 2
   y = varargin{2};
   axes(voifigH(2))
   plot(y,'b')
   title('Timeseries')
   axes(voifigH(3))
   power = abs(fft(y));
   if nargin > 5
	   TR = varargin{6};
	   x = 1:size(power,1); x = x/(TR*size(power,1));
	   nyquist = 1/(2*TR);
   else
	   x = 1:size(power,1);
	   nyquist = max(x)/2;
   end
   plot(x,(power),'bx')
   j = get(gca,'XLim');
   set(gca,'XLim',[j(1) nyquist])
	if ~(max(power(10:end)) == 0), 
   		set(gca,'YLim',[0 max(power(10:end))+.1*max(power(10:end))])
	end
   title('Power Spectrum')
end

% nargin == 4 ...user enters y2, but no filtering is specified
% nargin > 4 ...filtering of y to make y2 (plotted)
if nargin == 4,y2 = varargin{3};,disp('user entered y2 (red)'),end
if not(isempty(y2)) 
   axes(voifigH(2))
   if nargin == 4,hold on,elseif nargin > 4,hold off,end	% on if user entered, off if filtered.
   plot(y2,'r')
   if nargin == 4,title('User entered timeseries (blue and red)'),
   elseif nargin > 4,title('High-pass filtered timeseries (red)'),
   		ylabel('% Change')   
   end
   axes(voifigH(3))
   if nargin == 4,hold on,elseif nargin > 4,hold off,end	% on if user entered, off if filtered.
   plot(x,abs(fft(y2)),'ro')
   j = get(gca,'XLim');
   set(gca,'XLim',[j(1) nyquist])
   if nargin == 4,title('Power spectrum'),elseif nargin > 4,title('Filtered power spectrum'),end
end

if nargin > 4
   %figure;plot(y),title('unadjusted timeseries'),figure(voifigH(1))
   if (exist('y2'))
       if ~isempty(y2)    
			hold off, axes(voifigH(4)),cla,
			disp('		voistat ''voifig'': computing trial average.')
			[avg,ste,avgd]=voistat('trialavg',y2,avgover,'r',voifigH);
          varargout{2} = avgd;
      	    ylabel('% Change')  
      end
   end
   title('Trial Average')
end
if (exist('y2'))
   varargout{3} = y2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
case 'trialavg'
   % input: timeseries y, time pts to average over, [opt>>] color,axis handle,cond function, condition groups, ...
   % [opt] baseline timepoints, [opt] nosteplot, opt outlier threshold
   % cond function: for selective averaging
   %    condition number at start of trial,zeros between start points
   %    e.g. [1;0;0;0;3;0;0;0;2;0;0;0;4;0;0;0;1;0;0;0; etc...]
   %    condition group: row vector as long as the number of conditions, specifying which to plot together
   %    e.g. [1 1 2 3] creates three groups, with condition 1 or 2 plotted as one line, 3 and 4 as two separate lines
   %    COLOR can be a cell array the same size as the number of groups.

   disp('*	voistat trialavg')
   y = varargin{1};
   avgover = round(varargin{2});				% which time points to avg over
   stop = size(y,1);
   conditions = 1;
   groups = 1;
   
   if nargin > 9,	outthresh = varargin{9};,else outthresh = 0;,end
   if nargin > 8,	steplot = varargin{8}; , else steplot = 1; , end
   if nargin > 7,	basepoints = varargin{7}; , else basepoints = 1; , end
   if nargin > 4,	voifigH = varargin{4};,end
	try
		if ishandle(voifigH), voifigplot = 1;,else voifigplot = 0;,end
	catch
		voifigplot = 0;
	end

%basepoints
%if basepoints == 1,error('quitting!'),end


if outthresh
   % -------------------------------
   % * figure out images to exclude
   % -------------------------------
   meany = mean(y); stdy = std(y);
   include = (y < meany + stdy * outthresh & y > meany - stdy * outthresh);
   disp(['		voistat trialavg: removed ' num2str(sum(~include)) ' images at z = ' num2str(outthresh)])
else include = ones(1,length(y));
end


   % if condition function is specified, use it - otherwise build one to mark start of all trials
   if nargin > 5
       condfunct = varargin{5};
       groups = 1:max(condfunct);
   else
       condfunct = zeros(stop,1);
       for i = 1:avgover:stop-avgover
           condfunct(i) = 1;
       end
   end
   
   % if condition groups are specified, use them, otherwise each condition is a group.
   if nargin > 6
       groups = varargin{6};
   end
   
   conditions = 1:max(condfunct);       % index conditions
   trialsskipped = 0;
   % loop through groups.  
   % for each group, figure out which conditions it goes with, get those, and average.
   for g = 1:max(groups)    
        avg = zeros(1,avgover);
index = 1;
zeromagnitude = [];
        whichconds = conditions(groups == g);					% select trial type
        
	for j = 1:stop-avgover
									% loop thru all images j.
            if sum(condfunct(j) == whichconds) > 0             	% if current element matches any of the conditions
		trialy = y(j:j+avgover-1)';
		includey = include(j:j+avgover-1)';
		mybase = j+basepoints-1; 						% mybase is index of timepoints
		mybase = mybase(mybase > 0);					% make sure not to average off end
		mybaseinclude = include(mybase);				% points to include
		mybase = y(mybase);								% get baseline response
		mybase = mybase(mybaseinclude);					% only included points (images)
		
		if outthresh									
			try
				mybase = mybase(mybase > mean(mybase)-outthresh*std(mybase) ... 			% ----- remove outliers ---------
				& mybase < mean(mybase)+outthresh*std(mybase));
			catch
				disp(['Trial ' num2str(j) ': no baseline available, using first timepoint in trial.'])
				mybase = y(j);
			end
		end
		
		if ~isempty(mybase) & ~isempty(trialy) & ~(sum(includey == length(trialy))),

					%mybase = y(j);, % warning(['Group ' num2str(g) ' image ' num2str(j) ':mybase empty: using first timepoint in trial.']),end

			zeromagnitude(index) = mean(mybase);	% save signal at baseline or 1st time for conversion to % change later
													% basepoints go from 1 = 1st el. to avg over; 0 is right before 1st el
													% zeromagnitude(index) = mean(y(j+basepoints-1));

	        	avgd(index,:) = trialy - zeromagnitude(index);	% change from baseline or 1st timepoint
			avd(index,~includey) = Inf/Inf;					% NaN out the values for disincluded images.
			if sum(~includey) == length(trialy),trialsskipped = trialsskipped + 1;, end
                	index = index + 1;
		end
            end
        end

	disp(['		...voistat trialavg: averaging done, skipped ' num2str(trialsskipped) ' trials due to outlier values.'])

   % remove outliers
   % ----------------------------------------------
  if outthresh

    numout = 0;
    for tpoint = 1:avgover
	vals = avgd(:,tpoint); meanval = nanmean(vals); 
	cutoff = outthresh * std(vals);
			%avgd(vals > meanval+cutoff | vals < meanval - cutoff,tpoint) = 0/0; an alternative strategy.
	out = (vals > meanval+cutoff | vals < meanval - cutoff);
	avg(tpoint) = mean(vals(~out & ~isnan(vals)));
	ste(tpoint) = std(vals(~out & ~isnan(vals)))/ sqrt(size(vals(~out & ~isnan(vals)),1));
	numout = numout + sum(out);
    end
    if numout > 0, disp(['		voistat trialavg: removed ' num2str(numout) ' outliers at z = ' num2str(outthresh) ', group ' num2str(g)]), end
 
   else
 	avg = nanmean(avgd);
   	for tpoint = 1:avgover
		vals = avgd(:,tpoint);
		ste(tpoint) = std(vals(~isnan(vals)))/ sqrt(length(vals(~isnan(vals))));
	end
   end % end outliers
	


   % convert to percent change from first timepoint or baseline
	% bad idea if you've already mean centered your data (scan adjustment).
   % ----------------------------------------------
   %avg = 100*avg ./ mean(zeromagnitude);
   %ste = 100*ste ./ mean(zeromagnitude);


   if nargin > 4
      	color = varargin{3};
      	if iscell(varargin{3}),
          	if g > size(color,2),color{g} = 'k',end
          	plotcolor = color{g};,
	else plotcolor = color;
	end
   elseif nargin > 3,error('specify color and figure handle for plotting.')
   end
     

   if voifigplot
	if outthresh
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

	% ------ plot the selective average --------
      	try,axes(voifigH(4)),catch,disp('Can''t find voifig window - try InitWindows. Plotting on current axis.'),end
      	if g == 1, cla,end	% clear axis if first group.
      	plot(avg,plotcolor)
      	hold on
      	wid = get(gca,'XLim');
      	wid = (wid(2) - wid(1))/50;
      	if steplot
      	for i = 1:avgover-1
         	plot([i i],[avg(i)-ste(i) avg(i)+ste(i)],plotcolor)
         	plot([i-wid i+wid],[avg(i)-ste(i) avg(i)-ste(i)],plotcolor)
         	plot([i-wid i+wid],[avg(i)+ste(i) avg(i)+ste(i)],plotcolor)
      	end
      end
   end
   
   eval(['Avg.group' num2str(g) ' = avg;'])
   eval(['Ste.group' num2str(g) ' = avg;'])
   eval(['Avgtrials.group' num2str(g) ' = avgd;'])
   
end % loop thru groups

   varargout{1} = Avg;
   varargout{2} = Ste;
   varargout{3} = Avgtrials;
   
case 'plotall',
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
varargout{1} = voistat('figure');
gls = varargin{1};
voistat('fitvsy',gls,1);
voistat('fitvse',gls,2);
voistat('fits',gls,3);
if nargin > 2, 
   y = varargin{2};
   voistat('autocorr',gls,4,y);
else
   voistat('autocorr',gls,4);
end
   
case 'figure',
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
varargout{1} = figure;
%set(gcf,'Position',[360   123   619   811])
set(gcf,'Position',[232   258   560   420])

case 'fitvsy',
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
gls = varargin{1};
whichp = varargin{2};

% plot fits vs raw data
subplot(4,1,whichp)
plot(gls.y)
hold on
plot(gls.fits,'r--')
title('Raw data (blue) vs. Fits (red)')

case 'fitvse',
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
gls = varargin{1};
whichp = varargin{2};

% plot fits vs residuals
subplot(4,1,whichp)
plot(gls.fits,gls.e,'o')
hold on
plot([min(gls.fits) max(gls.fits)],[0 0],'r')
title('Residuals vs. fitted values')

case 'fits',
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
gls = varargin{1};
whichp = varargin{2};

% index plot of residuals
subplot(4,1,whichp)
plot(gls.e,'o')
hold on
plot([1 size(gls.e,1)],[0 0],'r')
title('Index plot of residuals')

case 'power',
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
gls = varargin{1};
whichp = varargin{2};
subplot(4,1,whichp)
if nargin > 3, y = varargin{3};, plot(abs(fft(y)),'gs'),
end
hold on
plot(abs(real(fft(gls.y))),'b+')
plot(abs(real(fft(gls.fits))),'rx')
set(gca,'XLim',[0 size(gls.y,1)/2])
title('Power spectra.  blue + = gls.y, red x = fits, green . = other var.')

case 'autocorr',
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
gls = varargin{1};
whichp = varargin{2};
subplot(4,1,whichp)	
e = gls.e;
xc = getv('get',e);
plot(xc)
title('Autocorrelation function of residuals.')

case 'adjusty',
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% out = voistat('adjusty',y,HChoice,TR,nscans, [opt] adjustmatrix); 
disp('*	voistat adjusty')
y = varargin{1};
HChoice = varargin{2};
TR = varargin{3};
npoints = size(y,1); 
if nargin > 4,nscans = varargin{4};,else nscans = 1;,end
if nargin > 5,adjustmatrix = varargin{5};,else adjustmatrix = [];,end

scanlen = npoints/nscans;
if ~(scanlen ~= round(scanlen)),warning(['Scan length is not an integer! scanlen = ' num2str(scanlen)]),end
scanlen = round(scanlen);


% MAKE FILTER - use_spm_filter.m
% --------------------------------------------------------
% Ran out of memory doing this with whole ts, so trying to make more efficient...do it scan by scan.
if ~isempty(HChoice)
	disp(	'		...voistat adjusty: making S filter')
	[S,KL,KH] = use_spm_filter(TR,scanlen,'none','specify',HChoice);
end

wholeyavg = mean(y);


% adjust scan means
% --------------------------------------------------------
disp('		...adjusting scan means to session mean.')	

index = 1;
for startimg = 1:scanlen:npoints
	X(startimg:startimg+scanlen-1,index) = 1;
	index = index + 1;
end

if ~isempty(adjustmatrix)
	disp('		...adjusting for user-entered covariates -e.g. movement params')
 	if sum(sum(isnan(adjustmatrix))) > 0,
		disp('		...WARNING: NaN values in user covariates!  Setting these to 0.')	
		adjustmatrix(isnan(adjustmatrix)) = 0;
	end
	try
		X = [X adjustmatrix];
	catch
		warning('		voistat ''adjusty'': X and adjustmatrix are different sizes.  no covariates entered.')
		whos X
		whos adjustmatrix
	end
end

%X = [X ones(size(X,1),1)];
%figure;imagesc(X);
%figure;plot(y);
%hold on; plot(X * (X \ y),'r');


y = y - (X * (X \ y));	% y-xB to adjust to mean.




filtmethod = 'spm';
	% options are: spm, cheby, fourier
switch filtmethod

case 'spm'
% filter each scan
% --------------------------------------------------------
if ~isempty(HChoice)
	filty = [];
	disp('		...applying to each session: high-pass filtering.')
	for startimg = 1:scanlen:npoints
		sessy = y(startimg:startimg+scanlen-1);
		sessy = S * sessy;
		filty = [filty; sessy];
	end
	y = filty;

       %Wn = .35;
       %[B,A] = cheby2(3,40,Wn);
       %y = filter(B,A,y); 

% Doug's custom fourier LP filter. 
len = length(y) /2;
ind = [1:len] ./ len; 
ff = 1./(1+exp((ind-.30)./0.05));
ff2 = [ff ff(len:-1:1)];
y = real(ifft(fft(y).*ff2'));

end


case 'cheby'
if ~isempty(HChoice)
	nyquist = TR/2;	% in seconds, not frequencies TR/nyquist = .5
	Wn = [TR/HChoice .35];
	disp([		'voistat adjusty: Chebyshev filter freqency is ' num2str(Wn)])
	for startimg = 1:scanlen:npoints
		sessy = y(startimg:startimg+scanlen-1);
 		[B,A] = cheby2(1,20,Wn);
		sessy = filter(B,A,y); 
	end
end

case 'fourier'
if ~isempty(HChoice)
	filty = []; nyquist = (1/TR)/2; lowerlim = 1/HChoice;
	hzfreqs = (1:scanlen) ./ (scanlen * TR);
	ftopass = hzfreqs >= lowerlim & hzfreqs <= nyquist;
	
	disp(['		...applying to each session: notch filtering: ' num2str(lowerlim) ' to ' num2str(nyquist) ' Hz.'])
	for startimg = 1:scanlen:npoints
		myfft = fft(y(startimg:startimg+scanlen-1));
		mypass = zeros(1,scanlen);
		mypass(ftopass) = myfft(ftopass);
		% plot(abs(myfft));hold on;plot(abs(mypass),'rx');	
		sessy = real(ifft(mypass))';
		filty = [filty; sessy];
	end
	y = filty;
end

end % end switch




disp('		...converting to % change from overall mean.')
y = (y-mean(y))*100 / wholeyavg;		% percent change from 0.
varargout{1} = y;


% figure; plot(abs(fft(y)),'ro');set(gca,'YLim',[0 1000])








case 'ortho',
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end switch   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
