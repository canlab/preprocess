function [pupil_y, pupil_y_hf] = physio_make_pupil_matrix(flocation, subjlist)
% function [pupil_y, pupil_y_hf] = physio_make_pupil_matrix(y, y_hf)
%
% flocation = folder containing .mat files
% subjlist = list of subject numbers (give full path name)
%
% Matt Furstoss and Jason Buhle Sept 2009

%%
clear all

source_frequency = 1/60;

path_dir_temp = 'C:\Users\Jason\Documents\projects\research\Deception\TruthReversalRT\analyses\temp\';
file_behavioral= 'data_behavioral';
file_ID_list = 'ID_list_eye';
cd(path_dir_temp); %#ok<MFAMB>
load(file_behavioral);  %#ok<MFAMB>
load(file_ID_list);

flocation = 'C:\Users\Jason\Documents\projects\research\Deception\TruthReversalRT\analyses\temp\eye_data\pupil_processed\';
cd(flocation);

IDs_by_trial =data_behavioral.ID;
TrialOnsets = data_behavioral.TrialOnset;
TrialOffsets = data_behavioral.TrialOffset;
trial_durations = TrialOffsets - TrialOnsets;
Runs_by_trial = data_behavioral.Run;

% wh_ID = 12; wh_run = 1;

for wh_ID = 1:length(ID_list_eye)

    % Print Progress to the Command Window
    fprintf(1,'ID:  ');
    fprintf(1,'\b%d',wh_ID); pause(.1)
    fprintf('\n')

    wh = strcmp(IDs_by_trial, char(ID_list_eye(wh_ID)));

    if isempty(find(wh, 1)) %#ok<MFAMB>
        fprintf(1,'\tNo behavioral data'); pause(.1)
        fprintf('\n')
        continue
    end

    current_ID_TrialOnsets = TrialOnsets(wh);
    current_ID_TrialOffsets = TrialOffsets(wh);
    current_ID_trial_durations = trial_durations(wh);
    current_ID_runs = Runs_by_trial(wh);


    for wh_run = 1:4

        % Print Progress to the Command Window
        fprintf(1,'\trun:  ');
        fprintf(1,'\b%d',wh_run); pause(.1)
        fprintf('\n')

        wh = find(current_ID_runs==wh_run);

        if isempty(find(wh, 1)) %#ok<MFAMB>
            fprintf(1,'\tNo data for this run'); pause(.1)
            fprintf('\n')
            continue
        end

        current_ID_run_TrialOnsets = current_ID_TrialOnsets(wh);
        current_ID_run_TrialOffsets = current_ID_TrialOffsets(wh);
        current_ID_run_trial_durations = current_ID_trial_durations(wh)+1;
        current_max_duration = max(current_ID_run_trial_durations);

        % load file
        fname = strcat(flocation, (char(ID_list_eye(wh_ID))), '_', num2str(wh_run), '_ys'); %#ok<MFAMB>
        load(fname);

        %        pause(.1);

        % transpose y and y_hf
        y = y'; y_hf = y_hf';

%         resample from 60 to 1000 Hz
                y_1000hz = resample(y, 1000, 60);
                y_hf_1000hz = resample(y_hf, 1000, 60);

%         current_run_length = length(y);
% 
%         target_frequency = 1/1000;
% 
%         current_duration = source_frequency*current_run_length;
%         source_time_vector = source_frequency:source_frequency:current_duration;
%         target_time_vector = target_frequency:target_frequency:current_duration;
% 
%         y_1000hz = spline(source_time_vector, y, target_time_vector);
%         y_hf_1000hz = spline(source_time_vector, y_hf, target_time_vector);


        % create matrices
        pupil_y = NaN(48,current_max_duration);
        pupil_y_hf = NaN(48,current_max_duration);

        % divide the runs by trial and store in matrices, n rows for n
        % trials
        for wh_trial = 1:48
            start = current_ID_run_TrialOnsets(wh_trial);
            stop = current_ID_run_TrialOffsets(wh_trial);
            duration = (current_ID_run_trial_durations(wh_trial));
            pupil_y(wh_trial,1:duration)= y_1000hz(start:stop);
            pupil_y_hf(wh_trial,1:duration)=y_hf_1000hz(start:stop);
        end

        save_dir = 'C:\Users\Jason\Documents\projects\research\Deception\TruthReversalRT\analyses\temp\eye_data\pupil_processed\by_run_trial\';
        save_file = strcat(char(ID_list_eye(wh_ID)), '_', num2str(wh_run), '_pupil_y.mat');
        save_path = strcat(save_dir, save_file);
        save(save_path, 'pupil_y'); %#ok<MFAMB>
        save_file = strcat(char(ID_list_eye(wh_ID)), '_', num2str(wh_run), '_pupil_y_hf.mat');
        save_path = strcat(save_dir, save_file);
        save(save_path, 'pupil_y_hf');

    end
end