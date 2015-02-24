function [beat_onsets, bpm] = physio_make_hr_matrix(flocation, subjlist)
% function [beat_onsets, bpm] = physio_make_pupil_matrix(y, y_hf)
%
% flocation = folder containing .mat files
% subjlist = list of subject numbers (give full path name)
%
% Matt Furstoss and Jason Buhle Sept 2009

%%
clear all

path_dir_temp = '/Volumes/deception/analyses/temp/';
file_behavioral= 'data_behavioral';
cd(path_dir_temp); %#ok<MFAMB>
load(file_behavioral);  %#ok<MFAMB>
%%
flocation = '/Volumes/deception/analyses/temp/physio_processed/';
cd(flocation);

IDs_by_trial = data_behavioral.ID;
TrialOnsets = data_behavioral.TrialOnset;
TrialOffsets = data_behavioral.TrialOffset;
trial_durations = TrialOffsets - TrialOnsets + 1;
Runs_by_trial = data_behavioral.Run;
ID_list = unique(IDs_by_trial);

index_error_list = {};

%%

for wh_ID = 9:length(ID_list)

    % Print Progress to the Command Window
    fprintf(1,'ID:  ');
    fprintf(1,'\b%d',wh_ID); pause(.1)
    fprintf('\n')

    wh = strcmp(IDs_by_trial, char(ID_list(wh_ID)));
    
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
        fname = strcat(flocation, (char(ID_list(wh_ID))), '_TR_', num2str(wh_run), '_physio'); %#ok<MFAMB>
        % check if file exists *******
        chk_fname = strcat(fname, '.mat');
        temp = exist(chk_fname,'file');
        if temp == 0 
            fprintf(1,'\tFile does not exist'); pause(.1)
            fprintf('\n')
            continue
        end
        load(fname);

        %        pause(.1);

        % transpose y and y_hf **** modify for hr.mat files, bpm and onsets
        % instead of y and y_hf
        beat_onsets = physio.HR.beat_onsets; bpm = physio.HR.bpm;
        beat_onsets = beat_onsets'; 

        % resample from 60 to 1000 Hz
%         y_1000hz = resample(y, 1000, 60);
%         y_hf_1000hz = resample(y_hf, 1000, 60);

        % create matrices
        hr_beat_onsets = NaN(48,current_max_duration - 1);
        hr_bpm = NaN(48,current_max_duration);
        
        curr_length = length(bpm);
        
        % divide the runs by trial and store in matrices, n rows for n trials
        for wh_trial = 1:48
            start = current_ID_run_TrialOnsets(wh_trial);
            stop = current_ID_run_TrialOffsets(wh_trial);
            if stop < curr_length
                duration = (current_ID_run_trial_durations(wh_trial));
                hr_beat_onsets(wh_trial,1:duration-1) = beat_onsets(start:stop);
                hr_bpm(wh_trial,1:duration-1)=bpm(start:stop);
            else index_error_list{end+1} = strcat(char(ID_list(wh_ID)), '_', num2str(wh_run));     %#ok<AGROW>
            end
        end

        save_dir = '/Volumes/deception/analyses/temp/physio_processed/by_run_trial/';
        save_file = strcat(char(ID_list(wh_ID)), '_', num2str(wh_run), '_beat_onsets.mat');
        save_path = strcat(save_dir, save_file);
        save(save_path, 'beat_onsets'); %#ok<MFAMB>
        save_file = strcat(char(ID_list(wh_ID)), '_', num2str(wh_run), '_bpm.mat');
        save_path = strcat(save_dir, save_file);
        save(save_path, 'bpm');    
    end
    end
save '/Volumes/deception/analyses/temp/physio_processed/by_run_trial/index_error_list.mat' index_error_list;