% [bpm_in_s, bpm, x_sec, onsets] = physio_pulse2hr_new(pulse)
%
% [bpm_in_s, bpm] = physio_pulse2hr_new(physio.cei1.pulse, ['start', 40000, 'clear'])
% 'start' : starts manual check at n points in timeseries
% 'clear' : clears old onsets from gui
% 'samprate' : followed by samples per second of pulse
% 'cutoff' : followed by cutoff_in_s for ID'ing onsets
% 
%
%
% tor wager, july 2007
% Matt Furstoss, June 2009

function [bpm_in_s, bpm, x_sec, onsets, pulse] = physio_pulse2hr_new(pulse, missing, varargin) 
    % defaults
    %  -------------------------------------------------------
    NPLOTS = 3;
    HR_AXES = 1;
    PULSE_AXES = 2;
    FIXED_AXES = 3;
    
    samprate = 100;
    cutoff_in_s = .5;   % this would be 150 BPM


    start_looking_at = 0;
    clearvals = 0;
    resume = 0;


    if ~isempty(varargin)
        for i = 1:length(varargin)
            if ischar(varargin{i})
                switch varargin{i}
                    case 'start', start_looking_at = varargin{i + 1};
                    case 'clear', clearvals = 1;
                    case 'resume', resume = 1;

                    case 'samprate', samprate = varargin{i + 1};
                    case 'cutoff', cutoff_in_s = varargin{i + 1};
                end
            end
        end
    end

    cutoff = samprate * cutoff_in_s;

    % init figure; first panel
    %  -------------------------------------------------------
    f1 = create_figure('Heart Rate Plot', NPLOTS, 1);
    axh(HR_AXES) = subplot(NPLOTS, 1, HR_AXES);

    if clearvals
        disp('Clearing old onsets from figure data.');
        data = [];
        data.last_start = start_looking_at;
        guidata(f1, data);
    end

    nobs = length(pulse);
    secs = linspace(0, nobs ./ samprate, nobs);

    plot(pulse);
    basefit = moving_average('gaussian', pulse, 3*samprate); % 3 sec FWHM moving average
    hold on; 
    plot(basefit, 'r')
    disp('Blue is pulse waveform, red is moving average that will be subtracted off.');
    pause(1) 
    cla
    pulse = pulse - basefit;
    runavg = basefit; % *MF* save moving_average in physio
    plot(pulse)

    axh(PULSE_AXES) = subplot(NPLOTS, 1, PULSE_AXES);
    
    % get vector of pulse onsets, pulse_indicator
    % --------------------------------------------
    % load, or store in figure
    data = guidata(f1);
    if ~isempty(data) && isfield(data, 'onsets') && isfield(data, 'pulse_indicator')
        disp('Found existing onsets in figure data.');
        if resume
            fprintf('Resuming data at position %d.\n', data.last_start);
            start_looking_at = data.last_start;
        end
    else
        pulse_indicator = pulse;
        pulse_indicator(pulse_indicator < 0) = 0;
        
        input('Adjust axes with zoom tools until you can see the HR peaks clearly (in blue) and press return.')
        disp('Click on the top axis so that just the r-spike peaks (in blue) are above the horizontal line.')
        [x, y] = ginput(1);  
        if y < 0, y = .01; end
        
        plot_horizontal_line(y);

        pulse_indicator = pulse_indicator > y; %.4;

        onsets = (diff(pulse_indicator) == 1);

        % get rid of onsets that occur too soon (false alarms)
        wh = find(onsets);
        samps_from_prev = [Inf; diff(wh)];
        bad = (samps_from_prev <= cutoff);

        onsets(wh(bad)) = 0;

        fprintf('Initial estimate: %3.0f heart-beat events.\n', sum(onsets));

        data.onsets = onsets;
        data.pulse_indicator = pulse_indicator;
        guidata(f1, data);
    end
    
    
    onsets = data.onsets;
    pulse_indicator = data.pulse_indicator;
    wh = find(onsets); %#ok

    axes(axh(PULSE_AXES))
    hold off;
    plot(pulse_indicator); 
    set(gca, 'YLim', [-.5 1.5]);
    hold on;
    han_onsets = plot_onsets(onsets, 'r', -.5, 2);
    set(axh(PULSE_AXES), 'YTick', []);
    
%     axes(axh(PULSE_AXES))
%     plot(zscore(pulse), 'k');
     
    drawnow
    %% Get missing beats manually
    % ------------------------------------------------
    disp('Showing a 20 s period of data.')
    fprintf('Click on middle panel to add onsets.\nShift+click to remove onsets.\nCtrl+click to move on to next series.\nApple-click when done.\n');

    viewsize = samprate * 20;
    xlim = [start_looking_at start_looking_at + viewsize];
    set(axh, 'XLim', xlim);
    linkaxes(axh, 'x');

    get_clicks_and_update_onsets();
    
    
    %% final plot of HR, sampled at samprate
    % 
    % onsets are 0, 1 for onsets of beats
    % beat2beat is b2b intervals for each heart-beat
    % x is the sample time of the midpoint between each beat and the next
    % bpm_at_x is the beats-per-min estimate for each heart-beat
    % bpm is at the original sampling resolution
    % ---------------------------------------------------

    data = guidata(f1);
    onsets = data.onsets;

    xlim = [0 nobs];
    set(axh, 'XLim', xlim);

    wh = find(onsets);
    x = wh(1:end-1) + diff(wh) ./ 2;    % x values: 1/2 the time to the next beat

    beat2beat = diff(wh ./ samprate);   % in s
    bpm_at_x = 60 ./ beat2beat;  % beats per min, given at x
    
    usable_pts = 1:nobs;
    usable_pts(missing) = NaN;
    
    % remove "bad" points from x and bpm_at_x, if any
    % round is used because some bad points fall between integers
    [dummy, whbadbeats] = intersect(round(x), missing);
    x(whbadbeats) = [];
    bpm_at_x(whbadbeats) = [];
 
     bpm = interp1(x, bpm_at_x, usable_pts, 'linear', 'extrap', NaN); % *MF*
    %bpm = moving_average('gaussian', bpm', 10 * samprate);   % 10-sec gaussian moving avg kernel

    axh(FIXED_AXES) = subplot(NPLOTS, 1, FIXED_AXES);
    
    linkaxes(axh, 'x');

    axes(axh(FIXED_AXES))
  %  plot(x, bpm_at_x, 'kx');  % black dots where we have actual samples
    hold on; 
    plot(bpm, 'r', 'LineWidth', 2);
%     plot(missing, 10*ones(length(missing)), 'ms', 'MarkerFaceColor', 'm')
    axis tight;

    disp('Red: beats per minute, magenta: missing data, black x''s: actual samples of each beat');
    
    % low-res, in seconds
    bpm_in_s = bpm(1 : samprate : end);

    x_sec = 0 : nobs ./ samprate;
    hold on; 
    plot(1 : samprate : nobs, bpm_in_s, 'g');

    disp('Green: beats per minute, sampled in seconds');
    
    % Get R-R and HRV
% see physio_hrv

    % --------------------------------------------------------
    % --------------------------------------------------------
    %% INLINE
    % --------------------------------------------------------
    % --------------------------------------------------------

    function get_clicks_and_update_onsets()
        X_TOL = 30;
        quitnow = 0;
        while ~quitnow
            data = guidata(f1);

            axes(axh(2));
            [xloc, yloc, button] = ginput(1);
            xloc = round(xloc);

            if ~isempty(button)
                switch(button)
                    case 1
                        data.onsets(xloc) = 1;
                        han_onsets(end+1) = plot_onsets(xloc, 'r', -.5, 2);
                    case 2
                        zero_range = max(1, xloc - X_TOL):min(length(data.onsets), xloc + X_TOL);
                        data.onsets(zero_range) = 0;
                        delete(han_onsets);
                        han_onsets = plot_onsets(data.onsets, 'r', -.5, 2);
                    case 3
                        xlim = get(axh(1), 'XLim');
                        stepsize = viewsize - 200;
                        data.last_start = xlim(1) + stepsize;
                        set(axh, 'XLim', xlim + stepsize);
                end
            else
                quitnow = 1;
                xlim = get(axh(1), 'XLim');
                fprintf('Last viewed range: %d %d\n', xlim(1), xlim(2));
            end

            guidata(f1, data);
        end
    end
end