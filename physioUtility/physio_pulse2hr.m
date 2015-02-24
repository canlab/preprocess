function [physio_data] = physio_pulse2hr(physio_data, doplot)
% [hr hrs] = physio_pulse2hr(pulse_data, [doplot])
%
% Analyses peripheral pulse data and returns heart rate information in raw
% and smoothed forms.
%
% hr = raw heart rate in BPM
% hrs = smoothed heart rate in BPM
% pulse_data = physio peripheral pulse data sampled at 100 Hz
% 
% E.g.: [physio.(subj)] = physio_pulse2hr(physio.(subj));
%
% THIS FUNCTION IS DEPRECATED. USE PHYSIO_PULSE2HR_NEW

disp('THIS FUNCTION IS DEPRECATED. USE PHYSIO_PULSE2HR_NEW')

samprate = 100; % sampling rate of pulse_data in Hz

pulse_data = physio_data.pulse;

hr = psd_welchlike_tor(physio_data.pulse,samprate*60,samprate,samprate,50/60,doplot);

% OLD way
% gp = gradient(pulse_data);
% iqrange = iqr(gp);
% thr = median(gp) + 10 * iqrange;
% wh = find(gp > thr);
% 
% % get earliest onsets in preceding 600 ms range
% whout = [];
% for i = 1:length(wh)
%     early = find(wh(1:i) - wh(i) > -60);
%     whout = [whout wh(early(1))];
% end
% 
% wh = unique(whout)';
% w2 = ([Inf; diff(wh)] ./ samprate);
% 
% omit = find(w2 < .6);  % omit > 120 BPM cause we're getting mid-beats
% w2(omit) = [];
% wh(omit) = [];
% 
% hr = 60 ./ w2;  % beat onset to onset / sampling -> seconds -> BPM
% 
% % interpolate
% hr = interp1(wh(1:length(hr)), hr, (1:length(pulse_data))');

hr = resample(hr,1,samprate);  % sampled in seconds
hr(isnan(hr)) = nanmean(hr);    % impute mean to NaNs

hr(1) = mean(hr(2:10));         % could be problem with extrapolation, first est.
hr = hr .* 60;                  % convert to beats ber min
hrs = smooth_timeseries(hr,12);

if(exist('doplot') && doplot)
    figure; plot(hr)
    hold on; plot(hrs,'r')
    legend({'original' 'smoothed'})
end

physio_data.hr = hr;
physio_data.hrs = hrs;

return