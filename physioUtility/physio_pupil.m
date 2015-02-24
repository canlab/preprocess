subjcodes = unique(data_eyetracker.ID);

% 1 sample every 1000/60 = 16.6 ms


% to covert from ms to samples, x * 60 / 1000
% to convert samples to sec, / 60


eye_sample_rate = 60;




% presume area
a1 = data_eyetracker.Pupil_A1;

temp = data_eyetracker.Run;
run_var = cell2mat(temp);

%% LOAD DATA FOR EACH SUBJECT AND REMOVE BLINKS

%probs with
% subj = 8 : 4144
% subj = 17 : 4180

% subj = 6;
 for subj = 1:length(subjcodes)
   

    %


    wh = strcmp(data_eyetracker.ID, subjcodes{subj});
    


    a1subj = a1(wh);
    run_subj = run_var(wh);
    
    for wh_run = 1:4
       wh = find(run_subj==wh_run);
       a1_subj_run = a1subj(wh);
       fname = strcat('/Volumes/Truth_Reversal/tasks/TR_001/pupil_prepped/', subjcodes{subj}, '_', num2str(wh_run));
       save(fname, 'a1_subj_run');
    end
 end
%%
% for i = 1.5 %:4
    iqrmult = 1.5;
    
  for subj = 1:length(subjcodes)
    for wh_run = 1:4
       wh = find(run_subj==wh_run);
       a1_subj_run = a1subj(wh);
       fname = strcat('/Volumes/Truth_Reversal/tasks/TR_001/pupil_prepped/', subjcodes{subj}, '_', num2str(wh_run));
       load(fname);
    
    
    [y,ntrimmed,spikes, yfit, spikelen] = splinetrim_eyetrack(a1_subj_run,iqrmult,5,[],'p');


len_cutoff = 2 * 60;     % 2 sec * 60 Hz
     y(spikelen >= len_cutoff) = NaN;


    create_figure('pupil data', 1, 1, 1);
    ylim([0 100])
    plot(y, 'r')
    ylim([0 100])

  
    mov_av_width = 10 * 60; % 10 sec * 60 samples
    y_ma = moving_average('gaussian', y, mov_av_width);
%     yfit_ma = moving_average('gaussian', yfit, mov_av_width);
    hold on; plot(y_ma, 'k');
    ylim([0 100])
    y_hf = y - y_ma;
%     yfit_hf = yfit - y_ma;

    
    path_file_TEMP = strcat('/Volumes/Truth_Reversal/tasks/TR_001/pupil_processed/', subjcodes{subj}, '_', num2str(wh_run), '_ys');
%     path_file_TEMP = fullfile(path_dir_temp, 'eye_data', 'interpolated', strcat(subjcodes{subj}, '_', num2str(i)));
    save(path_file_TEMP, 'y', 'y_hf'); %, 'yfit_hf');

         
     path_file_TEMP = strcat('/Volumes/Truth_Reversal/tasks/TR_001/pupil_processed/', subjcodes{subj}, '_', num2str(wh_run), '_figure');
%      path_file_TEMP = fullfile(path_dir_output, 'eye_data', 'individual', 'raw_moving_average', strcat(subjcodes{subj}, '_', num2str(i)));
    saveas(gcf, path_file_TEMP, 'png');

      close all

%  plot([yfit_hf y_hf])
% plot([yfit_ma y_ma])
% plot([yfit_ma])
% end
  end
end