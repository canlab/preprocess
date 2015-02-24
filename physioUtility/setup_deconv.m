function [physio, physio_onsets] = setup_deconv(physio, onsets, condition_names)

disdaqs = 7.5;      % disdaqs in s
samprate = 100;     % sampling rate of the scanner pulse in Hz
TR = 2;

if(~exist('condition_names') || isempty(condition_names))
    condition_names = {'Look_Neutral_Cue', 'Look_Neg_Cue', 'Reapp_Neg_Cue', 'Look_Neutral_Stim', 'Look_Neg_Stim', 'Reapp_Neg_Stim', 'Look_Neutral_Rating', 'Look_Neg_Rating', 'Reapp_Neg_Rating', 'Look_Neutral_Ant', 'Look_Neg_Ant', 'Reapp_Neg_Ant'};
end

subjs = {'s202' 's206' 's207' 's208' 's209' 's210' 's211'};
ni_codes = {'rea15' 'rea16' 'rea17' 'rea18' 'rea19' 'rea20' 'rea21'};
scanspersess = [192 196 196 184 190 192];     % images per session; length is num sessions
scanlengths = unique(scanspersess .* TR);

for h=1:length(subjs)
    subj = subjs{h};
    ni_code = ni_codes{h};

    % empty out and initialize onsets
    for i=1:length(condition_names)
        physio_onsets.(ni_code).(condition_names{i}) = [];
    end
    
    % empty out and initialize task data
    if(~isfield(physio.(ni_code), 'task_data'))
        physio.(ni_code).task_data = [];
    end
    for i=1:length(physio_fields)
        physio.(ni_code).task_data.(physio_fields{i}) = [];
    end


    % determine run periods
    scanner_periods = [];
    for i=1:length(scanlengths)
        scanner_periods = [scanner_periods; phys_get_scanperiods(physio.(ni_code).scanner_pulse, scanlengths(i))];
    end
    scanner_periods = unique(scanner_periods, 'rows');

    physio_fields = fieldnames(physio.(ni_code));
    physio_fields(strmatch('task_data', physio_fields)) = [];


    for i=1:size(scanner_periods, 1)
        st = scanner_periods(i,1) + (disdaqs * samprate);
        en = st + scanspersess(i) * samprate * TR - 1;
        for j=1:length(physio_fields)
            run_starts(i) = length(physio.(ni_code).task_data.(physio_fields{j})) / samprate; % what's the zero point of each run in the concatenated data?
            physio.(ni_code).task_data.(physio_fields{j}) = [physio.(ni_code).task_data.(physio_fields{j}); ...
                physio.(ni_code).(physio_fields{j})(st:en)];
        end

        for k=1:length(condition_names)
            physio_onsets.(ni_code).(condition_names{k}) = ...
                [physio_onsets.(ni_code).(condition_names{k}) (onsets.(subj).run{i}.(condition_names{k}) + run_starts(i))];
        end
    end
end

return