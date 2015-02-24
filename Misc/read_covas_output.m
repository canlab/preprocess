function out = read_covas_output(fnames, num_trials)
% out = read_covas_output(fnames, num_trials)
%
% Reads a series of covas .txt files entered in char or cell array of file
% names, and parses them into a matlab structure (out)
%
% num_trials: vector with number of trials in each covas text file
%
% Tor Wager and Math Roy, Nov 2011
%
% Examples:
% cd('/Volumes/current/BMRK3/874/Functional')
% covasfiles = filenames('run*/covas/*TXT', 'char', 'absolute');
% covasfiles =
% 
% /Volumes/current/BMRK3/874/Functional/run1/covas/874_1.TXT
% /Volumes/current/BMRK3/874/Functional/run2/covas/874_2.TXT
% /Volumes/current/BMRK3/874/Functional/run3/covas/874_3.TXT
% /Volumes/current/BMRK3/874/Functional/run4/covas/874_4.TXT
% /Volumes/current/BMRK3/874/Functional/run5/covas/874_5.TXT
% /Volumes/current/BMRK3/874/Functional/run6/covas/874_6.TXT
% /Volumes/current/BMRK3/874/Functional/run7/covas/874_7.TXT
% /Volumes/current/BMRK3/874/Functional/run8/covas/874_8.TXT
% /Volumes/current/BMRK3/874/Functional/run9/covas/874_9.TXT
% out = read_covas_output(covasfiles, [12 12 11 12 12 12 11 12 12]);

if iscell(fnames), fnames = char(fnames{:}); end

if isempty(fnames)
    disp('No covas files to parse. Skipping.');
    out = [];
    return
end

for i = 1:size(fnames, 1)
    
    myfile = deblank(fnames(i, :));
    
    if ~exist(myfile, 'file')
        fprintf('File does not exist: %s\n', myfile);
        %keyboard
        continue
    end
    
    try
    out(i) = parse_file(myfile, num_trials(i));
    catch
        fprintf('Problem parsing: %s\n', myfile);
        
        % try running : [a,b,c,d,e,f,g] = textread(myfile, '%d%f%f%f%f%f%f', num_trials(i), 'delimiter', '\n', 'headerlines', 4);
        % [a,b,c,d] = textread(myfile, '%d%f%f%f', 'delimiter', '\n', 'headerlines', 6 + num_trials(i));
        
        keyboard
    end
       
end



end  % main function


function out = parse_file(fname, num_trials)


% Number of Stimuli:	Baseline (C):	Rate (C/sec):	Destination (C):	Duration (sec):	Return Rate (C/sec):	Time:

[a,b,c,d,e,f,g] = textread(fname, '%d%f%f%f%f%f%f', num_trials, 'delimiter', '\n', 'headerlines', 4);
dat1 = [a b c d e f g];

% skip first line
dat1 = dat1(2:end, :);

%%

% Number of Stimuli:	Time:	Temperature:	VAS

[a,b,c,d] = textread(fname, '%d%f%f%f', 'delimiter', '\n', 'headerlines', 6 + num_trials);

dat2 = [a b c d];

%%

out = struct('trial_num', dat1(:, 1), 'temp_by_trial', dat1(:, 4), 'trial_num_by_time', dat2(:, 1), ...
    'temp_across_time', dat2(:, 3), 'time_in_sec', dat2(:, 2), 'VAS_by_time', dat2(:, 4), 'all_trial_dat', dat1, 'all_time_dat', dat2, ...
    'trial_data_names', [], 'time_data_names', []);

out.trial_data_names = {'Number of Stimuli' 'Baseline Temp (C)'	'Rate of change (C/sec)' 'Destination (C)' 'Duration (sec)' 'Return Rate (C/sec)' 'Time'};

out.time_data_names = {'Number of Stimuli' 'Time' 'Temperature' 'VAS'};

%%

end




