% eyetrack_read - reads in, computes stats, draws and saves results
%
% contrastnames = {'Av - Neu Cue'; 'Av - Neu Picture'};
% contrasts = [.25 .25 -.25 -.25 0 0 .25 .25 -.25 -.25; ... %av - neu cue
%     0 0 0 0 1 -1 0 0 0 0];    % av - neu pic
% Example for iaps-load:
% eyetrack_read(nm, ...
%     {'Av cue 1s p' 'Av cue 4s p' 'Neu cue 1s p' 'Neu cue 4s p' 'Av pic' 'Neu pic' 'Av cue 1s load' 'Av cue 4s load' 'Neu cue 1s load' 'Neu cue 4s load' 'Av pic load' 'Neu pic load' 'Av probe' 'Neu probe' 'Letter encode'}, ...
%     [1 2 3 4 5 6 11 12 13 14 15 16 20 21 30], ...
%     'basepts', [60 60 60 60 0 0 60 60 60 60 0 0 0 0 0], ...
%     'groupnames', {'Av cue' 'Neu cue' 'Av pic' 'Neu pic'}, ...
%     'groups', {[1 2 7 8] [3 4 9 10] [5 11] [6 12]}, ...
%     'contrastnames', {'Av - Neu Cue'; 'Av - Neu Picture'}, ...
%     'contrasts', [.25 .25 -.25 -.25 0 0 .25 .25 -.25 -.25 0 0 0 0 0;
%         0 0 0 0 0.5 -0.5 0 0 0 0 0.5 -0.5 0 0 0]);

function eyetrack_read(filename, condnames, condcodes, varargin)
    % defaults
    using_baseline = 1;
    basepts = ones(1, length(condnames)) .* 3;       % how many samples of baseline data collected before event onset for the various epochs
    samprate = 120;     % in Hz
    colors = {'k' 'r' 'g' 'b' 'c' 'm' 'y'};
    max_epoch_length = 3 * 60 * samprate;       % 3 minutes at 120 Hz
    min_epoch_length = 10;
    doplot = 1;

    % defined strings
    PARAMETERS_RECORDED = 'Parameters Recorded:';
    TOTAL_PTS_RECORDED = 'Total Pts Recorded:';
    SAMPLE_NBR = 'Sample #';

    % constants
    DEFAULT_BASELINE_PTS = 30;
    SERIN0_COLUMN = 5;
    PUPILA1_COLUMN = 7;
    START_RECORDING_CODE = 132;
    STOP_RECORDING_CODE = 136;
    SMOOTHVAL = 10;             % smoothing applied to timeseries data; zero impact after SMOOTHVAL samples

    if(length(varargin) > 0)
        for i=1:length(varargin)
            if(strcmp(varargin{i}, 'samprate')), samprate = varargin{i+1}; end
            if(strcmp(varargin{i}, 'basepts')), basepts = varargin{i+1}; end
            if(strcmp(varargin{i}, 'groups')), groups = varargin{i+1}; end
            if(strcmp(varargin{i}, 'groupnames')), groupnames = varargin{i+1}; end
            if(strcmp(varargin{i}, 'contrasts')), contrasts = varargin{i+1}; end
            if(strcmp(varargin{i}, 'contrastnames')), contrastnames = varargin{i+1}; end
            if(strcmp(varargin{i}, 'evtypemapping')), evtypemapping = varargin{i+1}; end % override for evtypes..or in case you forgot to get them
            if(strcmp(varargin{i}, 'max_epoch_length')), max_epoch_length = varargin{i+1}; end
            if(strcmp(varargin{i}, 'min_epoch_length')), min_epoch_length = varargin{i+1}; end
            if(strcmp(varargin{i}, 'using_baseline')), using_baseline = varargin{i+1}; end
            if(strcmp(varargin{i}, 'doplot')), doplot = varargin{i+1}; end
        end
    end

    % Parameter checking
    if(length(basepts) ~= length(condcodes))
        error('basepts not same length as condcodes - check inputs');
    elseif(length(condnames) ~= length(condcodes))
        error('condnames not same length as condcodes - check inputs');
    end

%     if(~exist('groups', 'var') || isempty(groups))
%         error('groups parameter is empty');
%     elseif(~exist('groupnames', 'var') || isempty(groupnames))
%         error('groupnames parameter is empty');
%     elseif(length(groupnames) ~= length(groups))
%         error('Number of group names doesn''t match number of groups');
%     end
%     
%     if(~exist('contrasts', 'var') || isempty(contrasts))
%         error('contrasts parameter is empty');
%     elseif(~exist('contrastnames', 'var') || isempty(contrastnames))
%         error('contrastnames parameter is empty');
%     elseif(length(contrastnames) ~= size(contrasts, 1))
%         error('Number of contrast names doesn''t match number of contrasts');
%     end


    [pathstr, name, ext, versn] = fileparts(filename);

    fprintf('\n\nProcessing %s\n', filename);

    % locate start of raw data
    start_line_nbr = locate_line(filename, SAMPLE_NBR);
    if(start_line_nbr == -1)
        fprintf('ERROR: Unable to locate data in file %s.', filename);
        return;
    end
    fprintf(1, 'Raw data section starts at line %i\n', start_line_nbr);

    % determine number of rows/columns to read
    num_lines_to_read = num_after_string(filename, TOTAL_PTS_RECORDED);
    fprintf(1, 'Reading %i lines\n', num_lines_to_read);
    num_columns_recorded = num_after_string(filename, PARAMETERS_RECORDED);

    % read in data
    d = dlmread(filename, '\t', [start_line_nbr 0 (start_line_nbr+num_lines_to_read-1) (num_columns_recorded-1)]);

    % starts and ends of epochs
    ends = find(d(:, SERIN0_COLUMN) == STOP_RECORDING_CODE);   % find ends of epochs (events)
    starts = [1; ends(1:end-1)+1];

    % determine event types
    events = d(:, SERIN0_COLUMN);
    events(events==START_RECORDING_CODE) = 0;
    events(events==STOP_RECORDING_CODE) = 0;
    % event_types = unique(events);
    % event_types(event_types==0) = [];


    dat_pupil = d(:, PUPILA1_COLUMN);         % input data
    clear('d');
    outdat = cell(1, length(condcodes));    % one cell per event type
    basedat = cell(1, length(condcodes));   % for saving pre-trial baseline values
    means = cell(1, length(condcodes));     % means across trials within a condition
    sterrs = cell(1, length(condcodes));    % standard errors across trials within a condition

    % ---------------------------------------------------------
    % process each trial and add it to data cell array outdat
    % ---------------------------------------------------------

    for i = 1:length(starts)

        dat_epoch = dat_pupil(starts(i):ends(i))';
        if(length(dat_epoch) > max_epoch_length || length(dat_epoch) < min_epoch_length)
            disp('Insanely long/short trial... wtf?');
            %keyboard;
            continue;
        end

        events_epoch = events(starts(i):ends(i))';

        if(exist('evtypemapping', 'var') && ~isempty(evtypemapping))
            evtype = evtypemapping(i);
        else
            evtype = unique(events_epoch(events_epoch>0));
        end

        if(isempty(evtype))
            warning(['Event type not found for event ' num2str(i)]);
        else

            if(length(evtype) > 1)
                warning('WARNING: Multiple event types found for epoch. Using %d.\n', evtype(1));
            end
            
            evtype = evtype(1);
            ev_idx = find(condcodes == evtype);

            % pad, if necessary
            [dat_epoch outdat{ev_idx}] = padwithnan(dat_epoch, outdat{ev_idx}, 2);

            % a bit of smoothing
            tmp = smooth_timeseries(dat_epoch(~isnan(dat_epoch)), SMOOTHVAL);
            dat_epoch(1:length(tmp)) = tmp;

            % set low blink threshold based on median
            if(median(dat_epoch) < 400), lowthresh = 100; else lowthresh = 400; end

            % trimming and artifact/blink detection
            quarterlength = round(length(dat_epoch) / 4);
            trim_epoch = [fliplr(dat_epoch(1:quarterlength)) dat_epoch fliplr(dat_epoch(end-quarterlength:end))];
            [dat_epoch_tmp, ntr, whtrim] = splinetrim(trim_epoch(~isnan(trim_epoch)), 3, 3, lowthresh);
            dat_epoch_tmp = dat_epoch_tmp(quarterlength+1:(end-quarterlength-1));

            if(all(isnan(dat_epoch_tmp)))
                % overall signal too low! get rid of low thresh
                [dat_epoch_tmp, ntr, whtrim] = splinetrim(trim_epoch(~isnan(trim_epoch)), 3, 3);
            end

            dat_epoch = dat_epoch_tmp'; 

            if ntr > length(dat_epoch)./2
                % we've trimmed too many points; get rid of trial
                dat_epoch = dat_epoch .* NaN;
            end

            %blinks{evtype}(end+1) = blinks_from_spline(dat_epoch, whtrim);

            %         % Do not use interpolated values, just exclude them
            %dat_epoch(whtrim) = NaN;

            % subtract baseline values
            if(using_baseline)
                basedat{ev_idx}(end+1) = nanmean(dat_epoch(1:max(DEFAULT_BASELINE_PTS, basepts(ev_idx))));
               dat_epoch = dat_epoch - basedat{ev_idx}(end);
            end

            outdat{ev_idx}(end+1, 1:length(dat_epoch)) = dat_epoch;
        end
    end

    % ---------------------------------------------------------
    % eliminate empty cells
    % ---------------------------------------------------------
    whempty = [];
    for i = 1:length(outdat)
        if(isempty(outdat{i}))
            whempty(end+1) = i;
        end
    end
    outdat(whempty) = [];
    basedat(whempty) = [];
    means(whempty) = [];
    sterrs(whempty) = [];
    condcodes(whempty) = [];
    condnames(whempty) = [];

    % ---------------------------------------------------------
    % get means and standard errors for each event type
    % ---------------------------------------------------------
    warning off
    fprintf(1, 'Calculating robust means across events...');
    for i = 1:length(outdat)
        [means{i}, tr, pr, sterrs{i}] = robust_mean(outdat{i});
        fprintf(1, '%3.0f ', i);
    end
    fprintf(1, '\n');
    warning on

    % ---------------------------------------------------------
    % get means and standard errors for each condition group
    % this is better because we collapse across all trials in
    % a particular condition, so we won't average means later
    % with diff number of trials in them
    % ---------------------------------------------------------
    if(exist('groups', 'var') && ~isempty(groups))
        warning off
        fprintf(1, 'Calculating robust means across conditions...');

        % means_to_avg = cell(max(groups), 1);

        for i = 1:length(groups)
            % Collapse trials for this condition
            for j = 1:length(groups{i})
                if(j == 1)
                    mtoavg = outdat{groups{i}(j)};
                else
                    tmp = outdat{groups{i}(j)};

                    [mtoavg, tmp] = padwithnan(mtoavg, tmp, 2);
                    mtoavg = [mtoavg; tmp];
                end
            end

            [groupmeans{i}, tr, pr, groupsterrs{i}] = robust_mean(mtoavg);
            fprintf(1, '%3.0f ', i);
        end

        warning on
        fprintf(1, '\n');
    end


    %%%%%%%%%%%%%%%%%%%%%%%%
    % Contrasts
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    if(exist('contrasts', 'var') && ~isempty(contrasts))
        for i = 1:length(means)
            sz(i) = length(means{i});
        end;

        for i = 1:size(contrasts, 1)   % for each contrast
            indices = find(contrasts(i, :));
            minsz = min(sz(indices));

            conout{i} = zeros(1, minsz);
            for j = indices
                conout{i} = conout{i} + contrasts(i, j) * means{j}(1:minsz);
            end
        end
    end

    % ---------------------------------------------------------
    % plot
    % ---------------------------------------------------------

    % set x-axis (secs)
    secindex = ([-max(basepts):max(cellfun(@length, outdat))]) ./ samprate;

    if(doplot)

        %%%%%%%%%%%%%%%%%%%%%%%
        % Standard plot
        %%%%%%%%%%%%%%%%%%%%%%%

        if(exist('groupmeans', 'var') && ~isempty(groupmeans))
            tor_fig;
            h = [];
            for i = 1:length(groupmeans)

                color = colors{rem(i, length(colors)) + 1};
                h(i) = plot(secindex(1:length(groupmeans{i})), groupmeans{i}, color, 'LineWidth', 2);
                hold on; fill_around_line(groupmeans{i}, groupsterrs{i}, color, secindex(1:length(groupmeans{i})));

            end
            xlabel('Time from event onset (s)')
            ylabel('Pupil diameter');
            legend(h, groupnames);
            drawnow
            saveas(gcf, name, 'fig');
        end


        %%%%%%%%%%%%%%%%%%%%%%%
        % Bar plot
        %%%%%%%%%%%%%%%%%%%%%%%

        for i = 1:length(basedat)
            bm(i) = nanmean(basedat{i}');
            bste(i) = ste(basedat{i}');
        end
        title(name)

        tor_fig;
        bar(bm);
        tor_bar_steplot(bm, bste, {'k'})
        ylabel('Baseline pupil dilation')
        set(gca, 'XTick', 1:length(basedat), 'XTickLabel', condnames);
        title(name)
        saveas(gcf, [name '-bar'], 'fig');


        %%%%%%%%%%%%%%%%%%%%%%%
        % Contrast plots
        %%%%%%%%%%%%%%%%%%%%%%%

        if(exist('conout', 'var') && ~isempty(conout))
            tor_fig;
            for i = 1:length(conout)
                color = colors{rem(i, length(colors)) + 1};
                h(i) = plot(secindex(1:length(conout{i})), conout{i}, color, 'LineWidth', 2);
                %hold on;
                %fill_around_line(means{i}, sterrs{i}, color, secindex(1:length(means{i})));
            end
            legend(contrastnames);
            title(name);
            saveas(gcf, [name '-contrasts'], 'fig');
        end
    end

    % ---------------------------------------------------------
    % save
    % ---------------------------------------------------------
    [pathstr, name, ext, versn] = fileparts(filename);
    str = ['saving ' name]; disp(str);
    save_var_names = {};
    var_names = {'dat_pupil', 'outdat', 'means', 'sterrs', 'basedat', 'samprate', 'secindex', ...
        'condnames', 'contrastnames', 'contrasts', 'conout', 'groupmeans', 'groupsterrs', 'groups', 'groupnames'};
    for i=1:length(var_names)
        if(exist(var_names{i}, 'var'))
            save_var_names{end+1} = var_names{i};
        end
    end
    save(name, save_var_names{:});
end





% Returns the number after the first occurrence of str
function number = num_after_string(filename, str)
    [fid, message] = fopen(filename);
    if(fid == -1)
        error(message);
    end

    try
        while(1)
            tline = fgetl(fid);
            if(~ischar(tline))
                fclose(fid);
                error('Did not find data in file %s.', filename);
                return;
            end
            if(~isempty(strfind(tline, str)))
                number = str2num(tline(length(str)+1:end));
                break;
            end
        end
        fclose(fid);
    catch
        fclose(fid);
    end
end