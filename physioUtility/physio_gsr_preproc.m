path_dir_temp = 'C:\Users\Jason\Documents\projects\research\Deception\TruthReversalRT\analyses\temp\';
file_behavioral= 'data_behavioral';

cd(path_dir_temp); %#ok<MFAMB>
load(file_behavioral);  %#ok<MFAMB>

subjcodes = unique(data_behavioral.ID);

% 1 sample every 1000/60 = 16.6 ms


% to covert from ms to samples, x * 60 / 1000
% to convert samples to sec, / 60


%% LOAD DATA FOR EACH SUBJECT AND REMOVE BLINKS

%probs with
% subj = 8 : 4144
% subj = 17 : 4180

% subj = 11; (maybe good)
 subj = 11; wh_run = 1;

    iqrmult = 1.5;
 for subj = 11:length(subjcodes)
    
    
    for wh_run = 1:4
       fname = strcat(path_dir_temp, '/physio_processed/', subjcodes{subj}, '_TR_', num2str(wh_run), '_', 'physio.mat');
       load(fname);

       curr_ID_run_gsr = physio.raw.data{3};
   % [y,ntrimmed,spikes, yfit, spikelen] = splinetrim_eyetrack(curr_ID_run_gsr,iqrmult,5,[],'p');


%len_cutoff = 2 * 60;     % 2 sec * 60 Hz
 %    y(spikelen >= len_cutoff) = NaN;

    y = curr_ID_run_gsr;
    time = (1:size(y, 1))';

    create_figure('gsr data', 1, 1, 1);
%     ylim([0 100])
    plot(y, 'r')
%     ylim([0 100])

  
    mov_av_width = 20000; % 10 sec * 60 samples
    y_ma = moving_average('gaussian', y, mov_av_width);
%     yfit_ma = moving_average('gaussian', yfit, mov_av_width);
    hold on; plot(y_ma, 'k');
%     ylim([0 100])
    y_hf = y - y_ma;
    plot(y_hf, 'b');
%     yfit_hf = yfit - y_ma;

    pause;
    
    path_file_TEMP = strcat(path_dir_temp, '/physio_processed/', subjcodes{subj}, '_', num2str(wh_run), '_ys');
%     path_file_TEMP = fullfile(path_dir_temp, 'eye_data', 'interpolated', strcat(subjcodes{subj}, '_', num2str(i)));
    save(path_file_TEMP, 'y', 'y_hf'); %, 'yfit_hf');

         
     path_file_TEMP = strcat(path_dir_temp, '/physio_processed/', subjcodes{subj}, '_', num2str(wh_run), '_figure');
%      path_file_TEMP = fullfile(path_dir_output, 'eye_data', 'individual', 'raw_moving_average', strcat(subjcodes{subj}, '_', num2str(i)));
    saveas(gcf, path_file_TEMP, 'png');

      close all

%  plot([yfit_hf y_hf])
% plot([yfit_ma y_ma])
% plot([yfit_ma])
% end
  end
end