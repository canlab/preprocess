function [GSRloc, GSRheight, humpdat, plothandles] = SCRdetect5(sig,thresh,fs, smoothw, timemin, humpmin, varargin)
    %
    % [GSRloc, GSRheight, humpdat, plothandles] = SCRdetect5(sig,thresh,fs, ... 
    % smoothw, timemin, humpmin, ['verbose'], ['noplot'])
    %
    % INPUTS:
    % ----------------------------------------------------------------------
    % sig, GSR signal, preprocessed (HP filtering and smoothing is
    % recommended)
    %
    % fs is sampling frequency in Hz, samples per second
    % smoothw is smoothing window in sec, or leave empty to use default of 1 s
    %       1st and 2nd derivatives will be smoothed with a Gaussian kernel
    %       of this FWHM
    %
    % thresh is the min threshold for smoothed gradient, or empty to use
    %       default of 90th percentile of derivative distribution
    %
    % timemin is min duration of time to count as a response (in sec),
    %       Empty to use default = 1.5
    %
    % humpmin is min height of inflection point to be counted as a response
    %
    % OUTPUTS:
    % ---------------------------------------------------------------------
    % GSRloc is a matrix: for each row (event), [onset peak offset]
    % plothandles is the handles to the final plots of time course, peaks,
    % and SCR event onsets/durations
    %
    % Examples:
    % ---------------------------------------------------------------------
    % % First, get GSR data with "bad" areas marked off from GUI output
    % (click save data, and then...)
    % y = cat(1, plotdata.clean_signal{1}{:});
    %
    % % Run SCRdetect with default values for all options, with 10 Hz
    % % signal, and 0 min SCR amplitude.
    %   [GSRloc, GSRheight, humpdat, plothandles] = SCRdetect5(y, [], 10, [], [], 0);
    %
    % % Run SCRdetect with arbitrary input threshold thresh, and 3 sec
    % % smoothing window, 2 sec min GSR duration
    %   [GSRloc, GSRheight, humpdat, plothandles] = SCRdetect5(y, thresh, 10, 3, 2, 0);
    %
    % % Run SCRdetect with .5 sec min GSR duration, plot/show all
    % output
    % [GSRloc, GSRheight, humpdat, plothandles] = SCRdetect5(y, [], 10, [], .5, 0, 'verbose');
    %
    % % Same as above, but suppress all plots
    % [GSRloc, GSRheight, humpdat, plothandles] = SCRdetect5(y, [], 10, [], .5, 0, 'noplot');

    % ---------------------------------------------------------------------
    % Defaults, and parse input arguments
    % ---------------------------------------------------------------------
    doverbose = 0;
    plothandles = [];
    plotfinal = 1;
    default_smoothw_s = 1;          % default smoothing window = 1 s
    default_thresh_prctile = 90;    % default gradient threshold for identifying possible events, in %ile of gradient distribution
    default_timemin = 1.5;
    
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}

                % functional commands
                case 'verbose', doverbose = 1;
                case 'noplot', plotfinal = 0;


                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

    % default smoothing window = 1 s
    if isempty(smoothw)

        smoothw = default_smoothw_s .* fs;
        if doverbose, fprintf('Empty smoothing window.  Using default of %3.0f samples = %3.2f sec\n', default_smoothw_s), end

    end

    if isempty(timemin)
        if doverbose, fprintf('Empty minimum time window.  Using default of %3.2f sec\n', default_timemin), end
        timemin = default_timemin;
    end
    
    if doverbose
        create_figure('Timeseries', 3, 1);
        plot(sig, 'k');
        ylim = get(gca, 'YLim');
    end

    n = length(sig);
    dsig = gradient(sig, 1/fs);          % derivative, fs is samp freq

    if isempty(thresh)
        if doverbose, fprintf('Empty Gradient threshold.  Using default of %3.0fth percentile, which is %3.2f.\n', default_thresh_prctile, thresh), end
        thresh = prctile(dsig, default_thresh_prctile);
    end

    if thresh < 0
        disp('Warning, threshhold set to < 0; changing threshhold to 0');
        thresh = 0;
    end

    % ---------------------------------------------------------------------
    % Get derivatives and possible SCR values
    % based on observations with large derivatives
    % ---------------------------------------------------------------------
    if doverbose
        create_figure('Timeseries', 3, 1, 1); subplot(3, 1, 2)
        plot(dsig);
    end

    dsig = smoothy(dsig, smoothw);      % smooth data, Gaussian

    ddsig = gradient(dsig,1/fs);        % ddsig is 2nd derivative
    ddsig = smoothy(ddsig, smoothw);

    poss_start = find(dsig > thresh);   % values with high slope, indicating events are occurring here

    if doverbose

        create_figure('Timeseries', 3, 1, 1); subplot(3, 1, 2)
        plot(dsig, 'r');
        title('Derivatives');
        legend({'Unsmoothed' 'Smoothed'});

        subplot(3, 1, 3)
        plot(ddsig);
        title('2nd Derivatives, smoothed');

        create_figure('Derivatives');
        hist(dsig, min(100, length(sig) ./ 10));
        plot_vertical_line(thresh);

        fprintf('Identified %3.0f GSR starts with large enough 1st derivatives.\n', length(poss_start))

    end


    % ---------------------------------------------------------------------
    % For each potential SCR, identify response
    % ---------------------------------------------------------------------

    min_start = 1;

    GSRloc = [];
    GSRheight = [];
    humpdat = [];
    i = 1;

    gsrindx = 1;

    while i <= length(poss_start)

        % Identify "zero points" or local maxima that may be SCR responses
        % ---------------------------------------------------------------------
        % for each high-slope point, go backwards to the end of the last
        % event and look for starting value
        % min_start is end of last response...     values w/decreasing
        % slope preceding "increasing slope" value

        % could put in limit?
        mask = dsig(min_start:poss_start(i)) > 0;                                   % only consider values with positive slope
        [start_val start_ind] = max(mask.*ddsig(min_start:poss_start(i)));          % max 2nd deriv : where slope is neg, but accel is most positive (inflection point)

        start_ind = start_ind + min_start - 1;  % adjusting...
        start_val = sig(start_ind);

        if doverbose, create_figure('Timeseries', 3, 1, 1); subplot(3, 1, 1), plot(start_ind, start_val, 'rx'); drawnow; end

        % Looking for max now ("hump"); these are candidates
        % We're looking for humps, so select only values with positive
        % acceleration as well.  Otherwise, we may identify dips, and this will
        % interfere with the identification of ending values at the
        % half-max point, causing failure to ID and end value and a skip to
        % the end of the timeseries.
        zero_points = find_intercept(dsig(start_ind:n), ddsig(start_ind:n));
        zero_points = zero_points + start_ind - 1;



        if doverbose
            fprintf('SCR #%3.0f: Found possible future peak values\n', gsrindx);
            gsrindx = gsrindx + 1;
            if exist('zerohandle', 'var') && all(ishandle(zerohandle)), delete(zerohandle), end
            zerohandle = plot(zero_points, sig(zero_points), 'kx', 'MarkerFaceColor', 'k'); drawnow;
        end


        % For each "zero point," identify possible ending values
        % based on return to 1/2 max response.
        % Look for later zero points within this window, which may
        % indicate "double humps" or additional SCRs built on last
        % ones, and stop when no higher humps can be identified.
        % ---------------------------------------------------------------------
        end_ind = 0;
        for j = 1:length(zero_points)

            point_height = sig(zero_points(j));

            % Need nanmean to avoid skipping to end
            end_val = nanmean([point_height start_val]);  % ending value to look for; 1/2 distance from starting (onset) value to peak

            poss_end = find(sig(zero_points(j):n) < end_val,1); % first point where signal drops below end value... or empty if it doesn't

            if isempty(poss_end)
                within_20_s_of_end = (n - start_ind) ./ fs < 20;
                if within_20_s_of_end
                    % signal never falls; use end of timeseries
                    end_ind = n;
                    break
                else
                    % signal never falls, but too far from end; maybe we're in an overall low point?
                    any_points_higher_than_hump = 0;  % this will skip this zero_point
                end
            else
                poss_end = poss_end + zero_points(j) - 1;  % little adjustment
                any_points_higher_than_hump = find(sig(zero_points(j)+2:poss_end) > sig(zero_points(j)), 1);
            end

            if doverbose && ~isempty(poss_end)
                create_figure('Timeseries', 3, 1, 1); subplot(3, 1, 1),
                if exist('possendhandle', 'var') && all(ishandle(possendhandle)), delete(possendhandle), end
                possendhandle = plot(poss_end, sig(poss_end), 'go', 'MarkerFaceColor', 'g');
                possendhandle(end+1) = plot(zero_points(j), sig(zero_points(j)), 'ks', 'MarkerFaceColor', 'k');
                possendhandle(end+1) = drawbox(start_ind, poss_end - start_ind, ylim(1), ylim(2) - ylim(1), [.5 .5 .5]);
                drawnow
                disp(sprintf('Found possible ending value for event starting at %3.0f, ending at %3.0f', poss_start(i), poss_end));

            end



            if isempty(any_points_higher_than_hump)
                % We can't find a later hump that is higher than this one,
                % but within the window of return to half-max;
                % We're done with this GSR event
                end_ind = poss_end;

                if doverbose, fprintf('Identified end of GSR event at %3.0f. Skipping rest of zero points\n', end_ind);
                    create_figure('Timeseries', 3, 1, 1); subplot(3, 1, 1), plot(end_ind, sig(end_ind), 'go', 'MarkerFaceColor', 'g');
                    if exist('zerohandle', 'var') && all(ishandle(zerohandle)), delete(zerohandle), end
                    if exist('possendhandle', 'var') && all(ishandle(possendhandle)), delete(possendhandle), end
                    drawnow
                end


                break
            end

        end

        if isempty(zero_points)
            end_ind = n;
        end

        % Test whether response is long enough to meet length
        % requirement and save it if it is.
        % ---------------------------------------------------------------------
        if (end_ind - start_ind)/fs > timemin && end_ind > poss_start(i)

            [peak_height peak_ind] = max(sig(start_ind:end_ind));
            peak_ind = peak_ind + start_ind - 1;
            peak_height = peak_height - start_val;

            if doverbose
                fprintf('GSR event starting at %3.0f and ending at %3.0f passed window length; saving.\n', start_ind, end_ind);
                drawbox(start_ind, end_ind - start_ind, ylim(1), ylim(2) - ylim(1), [.5 .7 .5]); drawnow
                input('Press return to continue.');
            end

            GSRloc = [GSRloc; [start_ind peak_ind end_ind]];
            GSRheight = [GSRheight; peak_height];
            humpdat = [humpdat; findhumps(sig(start_ind:end_ind), start_ind, humpmin)];

        elseif doverbose
            fprintf('Event does not meet window requirements.\n');
            if exist('possendhandle', 'var') && all(ishandle(possendhandle)), delete(possendhandle), end
        end
        if end_ind == 0
            min_start = poss_start(i) + 1;
        else
            min_start = end_ind + 1;
        end
        for k = 1:length(poss_start)
            if poss_start(k) > min_start;
                i = k;
                break
            elseif k==length(poss_start);
                i = k+1;
            end
        end

    end

    if plotfinal
        if doverbose, cla, end

        plothandles = plot_final_results(sig, GSRloc);

    end

end  % main function



% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Sub-functions
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------


function intercepts = find_intercept(dsig, ddsig)
    n = length(dsig);

    deriv_is_zero_btwn_points = sign(dsig(2:n)) ~= sign(dsig(1:n-1)) | abs(dsig(1:n-1)) < eps; % zero crossing, indicates local min or max

    % acceleration must be less than 0 (hump, not dip)
    possints = deriv_is_zero_btwn_points & ddsig(1:n-1) < 0;                   %%(dsig(2:n)~=0) .* (dsig(1:n-1)~=0)) + (dsig(1:n-1) == 0);

    possints = find(possints);

    % if zero-crossing is found between points but later point has smaller
    % derivative, then use later point
    for i = 1:length(possints)
        if possints(i) ~= n && (abs(dsig(possints(i))) > abs(dsig(possints(i)+1)));
            possints(i) = possints(i)+1;
        end
    end

    intercepts = possints;

end


function plothandles = plot_final_results(sig, GSRloc)

    hold on;
    
    plothandles = plot(GSRloc(:, 2), sig(GSRloc(:, 2)), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
    plothandles(end+1) = plot(sig, 'k', 'LineWidth', 2);
    
    nGSR = size(GSRloc, 1);

    ylim = get(gca, 'YLim');
    clear h1
    for i = 1:nGSR

        h1(i) = drawbox(GSRloc(i, 1), GSRloc(i, 3) - GSRloc(i, 1),ylim(1), ylim(2) - ylim(1), [.5 .5 .5]);

    end

    plothandles = [plothandles h1];
end




