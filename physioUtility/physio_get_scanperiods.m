function [scanner_onperiods,all_onperiods] = physio_get_scanperiods(scanner_pulse, task_length, samprate, varargin)
% [scanner_onperiods,all_onperiods] = phys_get_scanperiods(scanner_pulse, task_length, varargin)
%
% scanner_pulse         pulse data
% task_length           size of desired epoch in seconds (e.g., 500 if you have 250 vols with a TR of 2s)
% samprate              sample rate in Hz, (e.g. 1000 for .001 sec interval)
% varargin
%       'doplot'        whether or not to plot the scanner_pulse - defaults to 0
%       'tolerance'     percent difference of task length to include - defaults to 2
%       'nolength'      ignore task length and return all periods
%
% returns start and end indices for data
% scanner_pulse should be of the same length as the data (should be measured in 100ths or 1000ths of a s)

doplot = 0;
tol = .02;

if (length(varargin) > 0)
    for i = 1:length(varargin)
        if (strcmp(varargin{i},'doplot')); doplot = varargin{i+1}; end
        if (strcmp(varargin{i},'tolerance')); tol = varargin{i+1} / 100; end
        if (strcmp(varargin{i},'nolength')); tol = Inf; end
    end
end


kern = normpdf(-5:.1:5)./max(normpdf(-5:.1:5));

% locate most common value - assumes the most common val is when scanner is
% off, but this might not be true!
[h,x] = hist(scanner_pulse,100);
scanner_off_val = x(h == max(h)); if(~isscalar(scanner_off_val)), error('Multiple max vals determined... are you sure you passed in the scanner pulse?'); end

scanner_on = conv(abs(scanner_pulse - scanner_off_val),kern);
scanner_on = scanner_on(1:length(scanner_pulse));

scanner_on = scanner_on > std(scanner_on);
scanner_on = conv(scanner_on,kern);
scanner_on = scanner_on > 1;


if(doplot)
    figure; plot(scanner_pulse); hold on; plot(scanner_on,'r');
end


d = diff(scanner_on);
starts = find(d == 1);
ends = find(d == -1);
scanner_onperiods = [starts ends];
all_onperiods = scanner_onperiods;

len = round((scanner_onperiods(:,2) - scanner_onperiods(:,1)) ./ samprate);

if ~isinf(tol)
    scanner_onperiods = scanner_onperiods(find(len > ((1-tol)*task_length) & len < ((1+tol)*task_length)),:);
end

end


