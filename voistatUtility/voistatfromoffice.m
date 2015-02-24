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
% function [voifigH,avgdata] = voifig(voifigH,y,y2,avgover,[filt length, TR])
% ********************************************
% voifigH is 4 handles: figure, tsH,fftH,avgH
% if voifigH is a figure handle, activate.  if voifigH = [], create figure.
%
% optional arguments -
% y = timeseries to plot
% y2 = fits, or any other timeseries

% === if filter option is specified, make y2 the filtered data. ======
disp('y2 (red) is filtered data.')
if nargin > 4
	y = varargin{2};
   	avgover = varargin{4};
	if nargin > 5,
	   	HChoice = varargin{5};, TR = varargin{6};
   		y2 = voistat('adjusty',y,HChoice,TR);
	end
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

if nargin == 4,y2 = varargin{3};,disp('user entered y2 (red)'),end
if not(isempty(y2)) 
   axes(voifigH(2))
   if nargin == 4,hold on,elseif nargin > 4,hold off,end	% on if user entered, off if filtered.
   plot(y2,'r')
   if nargin == 4,title('User entered timeseries (blue and red)'),elseif nargin > 4,title('High-pass filtered timeseries (red)'),end
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
			hold off, axes(voifigH(4)),cla,[avg,ste,avgd]=voistat('trialavg',y2,avgover,'r',voifigH(4));
			varargout{2} = avgd;
		end
   end
   title('Trial Average')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
case 'trialavg'
   % input: timeseries y, time pts to average over, [opt>>] color,axis handle,cond function, condition groups
   % cond function: for selective averaging
   %    condition number at start of trial,zeros between start points
   %    e.g. [1;0;0;0;3;0;0;0;2;0;0;0;4;0;0;0;1;0;0;0; etc...]
   %    condition group: row vector as long as the number of conditions, specifying which to plot together
   %    e.g. [1 1 2 3] creates three groups, with condition 1 or 2 plotted as one line, 3 and 4 as two separate lines
   %    COLOR can be a cell array the same size as the number of groups.
   y = varargin{1};
   avgover = round(varargin{2});				% which time points to avg over
   stop = size(y,1);
   conditions = 1;
   groups = 1;
   
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
   
   % loop through groups.  
   % for each group, figure out which conditions it goes with, get those, and average.
   for g = 1:max(groups)    
        avg = zeros(1,avgover);
        index = 1;
        whichconds = conditions(groups == g)
        for j = 1:stop-avgover
            if sum(condfunct(j) == whichconds) > 0             % if current element matches any of the conditions
                avgd(index,:) = y(j:j+avgover-1)' - y(j);	   % change from 1st timepoint
                index = index + 1;
            end
        end
    
   %for i = 1:avgover:stop-avgover
	%   avgd(index,:) = y(i:i+avgover-1)' - y(i);	% change from 1st timepoint
	%   index = index+1;
    %end
   avg = mean(avgd);
   ste = std(avgd)/ sqrt(size(avgd,1));

   if nargin > 4
      color = varargin{3};
      if iscell(varargin{3}),
          if g > size(color,2),color{g} = 'k',end
          color = color{g};,end
   elseif nargin > 3,error('specify color and figure handle for plotting.')
   end
     
   if nargin > 4
      voifigH(4) = varargin{4};
      axes(voifigH(4))
      plot(avg,color)
      hold on
      wid = get(gca,'XLim');
      wid = (wid(2) - wid(1))/50;
      for i = 1:avgover-1
         plot([i i],[avg(i)-ste(i) avg(i)+ste(i)],color)
         plot([i-wid i+wid],[avg(i)-ste(i) avg(i)-ste(i)],color)
         plot([i-wid i+wid],[avg(i)+ste(i) avg(i)+ste(i)],color)
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
disp('Filtering and converting y to % change from zero.')
y = varargin{1};
HChoice = varargin{2};
TR = varargin{3};
npoints = size(y,1); 
% MAKE FILTER - use_spm_filter.m
[S,KL,KH] = use_spm_filter(TR,npoints,'none','specify',HChoice);
y = S * y;

y = (y-mean(y))*100 / mean(y);		% percent change from 0.

%freq = 1:npoints;
%freq = freq / (samprate * npoints);


varargout{1} = y;

case 'ortho',
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end switch   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
