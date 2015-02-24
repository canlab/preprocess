function HRV = physio_hrv(onsets, samprate)
% HRV = physio_hrv(onsets, samprate)
% 
% Calculate heart-rate variability and related measures
% from pulse (R-wave) onsets (at some high time res) and sampling rate in
% Hz
% 
% See also physio_pulse2hr_new, which calls this function
% onsets = physio.HR.onsets;
% samprate = physio.raw.samprate;
%
% Tor Wager, Dec 2008
%
% NOTE: Time varying HRV is BROKEN and does not work.

% Get R-R and HRV
    RR = diff(find(onsets > 1000*eps)) ./ samprate;
    create_figure('R-R', 2, 2); 
    plot(RR, 'LineWidth', 2); title('R-R interval');
    subplot(2, 2, 2);
    [myfft, freq] = fft_plot_scnlab(scale(RR, 1), 1 , 'samefig');
    axis auto; axis tight
    
    wh_hfhrv = find(freq > .15 & freq <= .4);
    hfhrv = mean(myfft(wh_hfhrv));
    
    wh_lfhrv = find(freq > .04 & freq <= .15);
    lfhrv = mean(myfft(wh_lfhrv));
    
    subplot(2, 2, 2);
    hold on;
    yl = get(gca,'YLim');
    drawbox(freq(wh_hfhrv(1)), freq(wh_hfhrv(end)) - freq(wh_hfhrv(1)), yl(2)-yl(2) .* .1, yl(2), 'b');
    drawbox(freq(wh_lfhrv(1)), freq(wh_lfhrv(end)) - freq(wh_lfhrv(1)), yl(2)-yl(2) .* .1, yl(2), 'g');

    drawnow
    
    HRV = struct('RR', RR, 'power', myfft, 'freq', freq, 'wh_hfhrv', wh_hfhrv, ...
        'wh_lfhrv', wh_lfhrv, 'hfhrv', hfhrv, 'lfhrv', lfhrv);
  
    % time-varying estimate of HRV
% %     
% %     hhrvhan = @(onsets) get_hrv(onsets, samprate, .15, .4);
% %     lhrvhan = @(onsets) get_hrv(onsets, samprate, .04, .15);
% %     
% %     % set up 20 sec Tukey window. 1 sec time resolution
% %     HRV.hfhrv_over_time = time_varying_estimate('tukey',onsets, 20 * samprate, hhrvhan, samprate);
% %     HRV.lfhrv_over_time = time_varying_estimate('tukey',onsets, 20 * samprate, lhrvhan, samprate);
% %     
% %     HRV.hfhrv_over_time = downsample(HRV.hfhrv_over_time, samprate);
% %     HRV.lfhrv_over_time = downsample(HRV.lfhrv_over_time, samprate);
% % 
% %     subplot(2, 2, 3);
% %     plot(HRV.hfhrv_over_time, 'b', 'LineWidth', 2)
% %     plot(HRV.lfhrv_over_time, 'g', 'LineWidth', 2)
end

function hrvpower = get_hrv(onsets, samprate, freq1, freq2)

RR = diff(find(onsets > 1000*eps)) ./ samprate;

%rr_samprate = length(RR) ./ (length(onsets) ./ samprate);

n = length(RR);
nyq = rr_samprate ./ 2; %1 ./ (2 * TR);
timepts = floor(n ./ 2);
freq = (1:timepts)/timepts * nyq;
myfft = fft(RR); %real(abs(fft(dat)));
myfft = abs(myfft(1:timepts)) .^ 2;  % power
wh_hrv = find(freq > freq1 & freq <= freq2);
hrvpower = mean(myfft(wh_hrv));

end

    
