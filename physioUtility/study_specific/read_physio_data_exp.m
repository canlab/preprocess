function [physio] = read_physio_data(filename,varargin)
%function [physio] = read_physio_data(filename,['physio', physio], ['nheaderrows', num header rows], etc.)
%
% OPTIONAL INPUTS:
% case 'nheaderrows', nheaderrows = varargin{i+1};
%  'mydelimiter', mydelimiter = varargin{i+1};
% 
%  'physio', physio = varargin{i+1};
%                     
%  'numc', followed by number of columns
%  'numrows', followed by number of rows
%                     
% 'plot', doplot = 1;
% 'noplot', doplot = 0;
% 'samprate', followed by sampling rate in Hz, e.g., 100 for samples every 0.01 sec
%                   
% 'nopulse'
% 'names', followed by variable names
% 'mydelimiter', followed by delimiter string for fields
%
% Examples:
%
% physio = read_physio_data(fname, 'numc', 6, 'numrows', 12200, 'nheaderrows', 5, 'plot');
%
% physio = read_physio_data(fname, 'numc', 6, 'nheaderrows', 5, 'plot', 'numrows', 12000, 'names', {'time' 'pulse' 'resp' 'gsr' 'gsr' 'scanner'});
%
% 3T Columbia fMRI Center Phillips
% -----------------------------
% First remove header rows, then check file for missing rows (#)
% physio = read_physio_data('/Users/tor/Desktop/SCANPHYS.log', 'numc', 10, ...
% 'nheaderrows', 0, 'plot', 'names', {'v1raw' 'v2raw' 'pulse' 'v2' 'ppu' 'resp' 'gx' 'gy' 'gz' 'mark'});
%
% physio = read_physio_data(fname, 'nheaderrows', 5, 'plot', 'samprate', ...
% 400, 'names', {'v1raw' 'v2raw' 'pulse' 'v2' 'ppu' 'resp' 'gx' 'gy' 'gz' 'mark'});
%
% 1.5T Columbia fMRI Center - NSF dataset
% ----------------------------------------
% physio2 = read_physio_data(fname, 'nheaderrows', 5, 'plot', 'samprate', ...
% 1/.01, 'names', {'time' 'pulse' 'resp1' 'resp2' 'gsr1' 'gsr2' 'scanner'}, 'numc', 7, 'mydelimiter', ' ');
%
% Location of HR files: /Volumes/NSF/physio_analysis/physio/

physio = [];

%% SET DEFAULTS and INPUTS
% --------------------------------------
% assume you have one header row
nheaderrows = 5;
mydelimiter = '\t';
numc = 0;
numrows = -1;
doplot = 1;
dopulse = 1;

samprate = [];
secs = [];
wh_scanner = [];
wh_time = [];
wh_pulse = [];

names = {};  %'time' 'pulse' 'resp' 'gsr' 'gsr' 'scanner'};

for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case 'nheaderrows', nheaderrows = varargin{i+1};
                case 'mydelimiter', mydelimiter = varargin{i+1};

                case 'physio', physio = varargin{i+1};
                    
                case 'numc', numc = varargin{i+1};
                case 'numrows', numrows = varargin{i + 1};
                    
                case 'plot', doplot = 1;
                case 'noplot', doplot = 0;
                    
                case {'nopulse', 'skippulse'}, dopulse = 0;
                    
                case 'names', names = varargin{i + 1};
                    
                case 'samprate', samprate = varargin{i + 1};
                    
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
end

if isempty(names)
    disp('You can enter descriptive names for each column as an input argument.')
    disp('Some of these are used to help process the data.')
    disp('Some used column names are : ''time'' ''pulse'' ''scanner''');
    
    for i = 1:numc, names{i} = ['col' num2str(i)]; end
end

%if ~iscell(files), files = {files}; end

wh_scanner = strmatch('scanner', names);
if isempty(wh_scanner)
    disp('No scanner pulse column name entered.  If these data were collected ')
    disp('during fMRI recording, enter ''scanner'' as one of the names ')
    disp('so I know which column contains scanner pulses.'); 
end

wh_time = strmatch('time', names);
if isempty(wh_time)
    disp('Enter ''time'' as one of the names so I know which column contains sample times in sec.'); 
end

wh_pulse = strmatch('pulse', names);
if isempty(wh_pulse)
    disp('No heart rate pulse column name entered.  ')
    disp('enter ''pulse'' as one of the names to run the heart-rate processing program here.')
end


%for i=1:length(files)
    disp(['Reading ' filename '...']);
    [tmp, basename, tmp] = fileparts(filename);
    
  
% READ THE FILE
% --------------------------------------
fprintf('Attempting to read your header. You entered %3.0f header rows.\n', nheaderrows);
hdr = textread(filename, '%s', 5, 'delimiter', mydelimiter, 'headerlines', 0);
disp(hdr)

disp('Reading the first line of your file:');

firstline = textread(filename, '%s', 1, 'delimiter', mydelimiter, 'headerlines', nheaderrows);
disp(firstline{:})
nums = nums_from_text(firstline{1});
if numc == 0
    fprintf('It looks like there are %3.0f columns.  Using this number to read data file.\n', length(nums));
    numc = length(nums);
    for i = 1:numc, names{i} = ['col' num2str(i)]; end
else
    fprintf('It looks like there are %3.0f columns.  You entered %3.0f.\n', length(nums), numc);
end

data = cell(1, numc);
[data{:}] = textread(filename, repmat('%f', 1, numc), numrows, 'delimiter', mydelimiter, 'headerlines', nheaderrows);

if isempty(data{1})
    disp('File may be empty!  No read errors, but no data read in.')
    return
end

% Get sampling rate
% Special format: sampling rate is first number in file
% Special format: time in seconds is first column in file, or whichever
% column has name 'time' in names input
% --------------------------------------
if isempty(samprate) && isempty(wh_time)
    fprintf('Trying to get sampling rate. ');

    try
        [fieldn, sec] = textread(filename, '%s%s%*[^\n]', 1, 'delimiter', '\t', 'headerlines', 0);
        sec = sec{1};
        wh = find(sec == ' ');
        sec(wh:end) = [];
        sec = str2num(sec);
        samprate = 1 ./ sec;
    catch
        disp('Error reading sampling rate. Using approx. from seconds in file.');
    end
end

if isempty(samprate) && isempty(wh_time)
    samprate = input('Enter sampling rate in Hz (e.g. 1000 for .001 sec interval): ');
elseif isempty(samprate) 
    % we have time column, but no samp rate
    samprate = round(1 ./ mean(diff(data{wh_time})));
end

if isempty(wh_time)
    % we have just samprate, no time column

    wh_time = length(data) + 1;
    data{wh_time} = (1:length(data{1}))' ./ samprate;
end
    
secs = data{wh_time};
len = length(secs);

if isempty(samprate)
    secs = data{wh_time}; %.* samprate;
end

if ~isempty(wh_scanner)
    
    run_length_TRs = input('Enter total run length, excluding disdaqs (discarded acquisitions), in images: ');
    disdaqs_TRs = input('Enter number of disdaqs: ');
    TR = input('Enter TR, repetition time for images (e.g, 2 for 2 sec TR): '); 
    
    total_run_secs = TR * (run_length_TRs + disdaqs_TRs);
    
    [scanner_onperiods, all_onperiods] = physio_get_scanperiods(data{wh_scanner}, total_run_secs, samprate); %max(secs), 'nolength');
    
    %**specific to this dataset**
    if isempty(scanner_onperiods)
        scanner_onperiods(1,1) = all_onperiods(1,1);
        scanner_onperiods(1,2) = all_onperiods(end);
    end
    
end

physio.raw.data = data; % save copy of raw data before removing stuff

% *MF*
% Finding skewed/missing spots in data
% Moving_average finds the spots of skewed/missing data & sets them to NaN
% missing, vector to store indices with value = NaN
% markers, store start and end indices of missing regions
% --------------------------------------
physio.movavg = moving_average('gaussian', data{wh_pulse}, 3*samprate); % 3 sec FWHM moving average

missing = [];
markers = [];
nreg = 0; % number of missing regions
stmrk = 0; %start marker, know whether to look for NaN or ~NaN
for i = 1:length(physio.movavg)     
    if isnan(physio.movavg(i))      
        missing = [missing i];
        
        if stmrk==0;
            markers(nreg+1,1) = i;
            stmrk = 1;
        end
    end
    if ~isnan(physio.movavg(i)) && stmrk==1
        markers(nreg+1,2) = i-1;
        nreg = nreg+1;
        stmrk = 0;
    end
end

% *MF*
% Make values of all other columns = NaN at appropriate indices *MF*
% for i=1:numc
%     for j=1:length(missing)
%         temp = missing(j);
%         data{i}(temp) = [NaN('double')]; 
%     end
% end
   data{wh_pulse}(missing) = [NaN('double')]; 


% *MF*
% PLOT
% -------------------------------------   
    % User decides whether or not to search/remove bad data in pulse
    yes = 'y';
    look_bad_data = input('Would you like to manually search for and remove skewed data?\n Enter y\\n: ','s');
    if strcmp(look_bad_data, yes)
        disp('You chose to search for bad data.')
        [new_missing, new_markers] = physio_remove_bad_data(data{wh_pulse}, missing, samprate);
        missing = new_missing;
        if ~isempty(new_markers)
            markers = [markers; new_markers];
        end
        % Make values of all other columns = NaN at appropriate indices
        % *MF*
%         if ~isempty(new_missing)
%             for i=1:numc
%                data{i}(missing) = NaN('double');
%             end
%         end
        data{wh_pulse}(missing) = NaN('double');
    else 
        disp('Moving on without searching for bad data.')
    end

% PLOT
% --------------------------------------
if doplot
    
    create_figure('physio file', numc, 1);
    
    for i = 1:numc
        subplot(numc, 1, i)
        
        hold on
        plot(secs, physio.raw.data{i}); %*MF*
        
        axis auto
        axis tight
    end
    
    if ~isempty(wh_scanner)
        disp('Detected scanner pulse, and saving times of scanner-on periods.');
        
        for j = 1:size(scanner_onperiods, 1)
           
            plot([scanner_onperiods(j, 1) scanner_onperiods(j, 2)] ./ samprate, [0 0], 'r', 'LineWidth', 6);
            
        end
        
    end
    
    subplot(numc, 1, 1)
    title(basename, 'FontSize', 18)
    
    ylabel(names{i})
    
    subplot(numc, 1, numc)
    xlabel('Time (sec)')
    drawnow
end

% save stuff

physio.filename = filename;
physio.basename = basename;

physio.names = names;
physio.wh_scanner = wh_scanner;
physio.wh_time = wh_time;

  
physio.raw.samprate = samprate;
physio.raw.secs = secs;
physio.HR.pulse = data{wh_pulse};
physio.missing = missing; % saves list to cell array


if ~isempty(wh_scanner)
    physio.raw.scanner_onperiods = scanner_onperiods;
    physio.raw.all_onperiods = all_onperiods;
    
    physio.TR = TR;
    physio.disdaqs_TRs = disdaqs_TRs;
    physio.vols_in_run = run_length_TRs;
    physio.total_run_secs = total_run_secs;
    physio.samples_in_scanning_run = scanner_onperiods(1, 1) + round(disdaqs_TRs * TR * samprate) : scanner_onperiods(1, 2);
    
    if doplot
        disp('Marking times I think the scanner starts and stops with vertical lines');
        disp('...and the run start time (after disdaqs) with a red vertical line');

        create_figure('physio file', numc, 1, 1);
        subplot(numc, 1, physio.wh_scanner);
        hold on;
        h1 = plot_vertical_line(scanner_onperiods(1, 1)./samprate); set(h1, 'Color', 'g', 'LineWidth', 3); 
        h1 = plot_vertical_line(scanner_onperiods(1, 2)./samprate); set(h1, 'Color', 'g', 'LineWidth', 3); 
        h1 = plot_vertical_line(physio.samples_in_scanning_run(1)./samprate); set(h1, 'Color', 'r', 'LineWidth', 3); 
        
        input('Press RETURN to continue');
        
    end
end 

% Downsample - to ~10 Hz, for easier data processing
% --------------------------------------
downsamplerate = round(samprate / 10);  % every n-th sample, target ~10 Hz

physio.downsampled.samprate = samprate ./ downsamplerate;
physio.downsampled.samprate_descrip = 'Sampling rate in Hz';
physio.downsampled.secs = downsample(data{wh_time}, downsamplerate);
physio.downsampled.secs_descrip = 'Times in sec';

for i = 1:numc
    physio.downsampled.data{i} = downsample(data{i}, downsamplerate);   
    
    % plot, if asked for
    if doplot
        if i ~= physio.wh_scanner
            subplot(numc, 1, i)
            
            plot(physio.downsampled.secs, physio.downsampled.data{i}, 'g');
        end
    end

end

if ~isempty(wh_scanner)
    physio.downsampled.scanner_onperiods = scanner_onperiods ./ downsamplerate;
    physio.downsampled.all_onperiods = all_onperiods ./ downsamplerate;
end 

disp('Check physio file figure: This shows you all data in the file.')
disp('Blue is raw data, and green is downsampled data.');

% HR
% --------------------------------------
if dopulse
    if ~isempty(wh_pulse)
        disp('Running heart rate processing')

        % *MF* set up to save (pulse-basefit) in physio,
        % added physio.missing to inputs
        [physio.HR.bpm_in_s, physio.HR.bpm, physio.HR.x_sec, physio.HR.beat_onsets, physio.HR.pulse] = physio_pulse2hr_new(data{wh_pulse}, ...
            physio.missing, 'clear', 'samprate', physio.raw.samprate);
        
        title(basename, 'FontSize', 18)
        
        if ~isempty(wh_scanner)
            physio.HR.scanner_onperiods = scanner_onperiods ./ physio.raw.samprate;
            physio.HR.all_onperiods = all_onperiods ./ physio.raw.samprate;
            
            physio.HR.HR_during_run = physio.HR.bpm(physio.samples_in_scanning_run);
            physio.HR.HR_during_run_in_TRs = downsample_scnlab(physio.HR.HR_during_run, physio.raw.samprate, 1./physio.TR);
            
            
        end

    else
        disp('Cannot do HR processing...I''m not sure which data has the pulse information.')
        disp('Enter ''pulse'' as a name in the names field.');
    end
else
    disp('HR pulse analysis option is skipped at user request')
end

end % main function