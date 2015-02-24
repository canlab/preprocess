function c = spm_config_gsranal(varargin)
% Configuration file for gsr analysis
%_______________________________________________________________________

addpath(fullfile(spm('dir'),'toolbox','GSRtoolbox'));


signals.type = 'files';
signals.name = 'GSR signal';
signals.tag  = 'signals';
signals.num  = [1 Inf];
signals.filter = '.mat';
signals.help   = {'Enter data for each session here. This should be a .mat file saved from the previous step "GSR preproc" '};

data.type = 'repeat';
data.name = 'Subject Data';
data.values = {signals};
data.num  = [1 Inf];
data.help = {'List of subjects.'};

normal.type = 'menu';
normal.name = 'Mode of normalization';
normal.tag  = 'normalisation';
normal.labels = {'by session', 'by subject', 'global', 'no normalisation'};
normal.values = {1, 2, 3, 0};
normal.help = { 'Range over which data is normalised.', ...
    'Choices are: by session, by subject, global, or none', ...
    'By session means that data will be transformed to z-scores within each session.', ...
    'By subject means that concatenated data will be transformed to z-scores across sessions.', ...
    'Global would (theoretically) z-score across all subjects/sessions (don''t do it.)', ...
};


smooth.type = 'entry';
smooth.name = 'Smoothing window';
smooth.strtype = 'r';
smooth.tag  = 'smoothingw';
smooth.num  = [1 1];
smooth.help = {'window size (FWHM in seconds) used for Gaussian data smoothing.'};


downs.type = 'entry';
downs.name = 'Output frequency';
downs.tag  = 'downsample';
downs.strtype = 'n';
downs.num  = [1 1];
downs.help = {'Desired output frequency for data (in Hz).', ...
        'If this is less than the acquisition frequency (should be the typical case), the data will be downsampled'...
        'If this is greater than the acquisition frequency, data will be interpolated'...
        'Otherwise, no adjustments are made.'};

    
minperiod.type = 'entry';
minperiod.name = 'Minimum SCR period';
minperiod.tag  = 'minperiod';
minperiod.strtype = 'r';
minperiod.num  = [1 1];
minperiod.help = {'The length of time in milliseconds that serves as the minimum acceptable SCR length.', ...
    'Very brief deflections should not be ID''d as SCRs.  We currently recommend 1.5 s.', ...
};

%minhump.type = 'entry';
%minhump.name = 'Minimum hump height';
%minhump.tag  = 'minhump';
%minhump.strtype = 'r';
%minhump.num  = [1 1];
%minhump.help = {'Minimum height for a hump to be thusly identified.', ...
%        'Value is a percentage of the total signal height. Default is .01'};


detect.type = 'menu';
detect.name = 'Mode of response detection';
detect.tag = 'detect';
detect.labels = {'Within each trial', 'Entire session', 'Entire subject', 'Global'}; 
detect.values = {1, 2, 3, 4};
detect.help = {'This is the mode that is used to produce SCR detection thresholds.', ...
    'Choices are : within each trial, entire session, entire subject, or global', ...
    'Within trial: Only consider periods between trial onset and offset', ...
    'Entire session: Find all SCRs for the session, irrespective of onsets', ...
    'Entire subject: Consider the whole time series at once when finding GSR events', ...
    '   ', ...
    'The algorithm works as follows: It uses the distribution of derivatives, and takes the ***th percentile  ' ...
    'of the distribution. Those are marked as likely events. It searches backwards for a point where the ' ...
    '2nd derivative is maximized to find the onset (inflection point), and forward to search for the peak ' ...
    '(zero in 1st derivative) and forward more to identify where the signal crosses half-max point ' ...
    ' defined by the onset-to-peak difference.' ...
    'So if you enter "within each trial", it does the entire analysis within each epoch.' ...
};

params.type = 'branch';
params.name = 'Perform analysis';
params.tag = 'parameters';
params.val = {detect, minperiod};
params.help = {'Set parameters for analysis.'};


nullb.type = 'branch';
nullb.name = 'Do nothing';
nullb.tag = 'nullb';
nullb.val = {};
nullb.help = {'Do nothing.'};


analchoice.type = 'choice';
analchoice.name = 'Autodetect SCRs?';
analchoice.tag = 'analchoice';
analchoice.values = {params, nullb};

c.type = 'branch';
c.name = 'Display data';
c.tag  = 'gsranal';
c.val  = {data, smooth, downs, normal, analchoice};
c.prog = @spm_gsranal;
c.help = {'Identifies and quantifies GSR data.',...
'Requires ascii text input.'};




end
%_______________________________________________________________________

%_______________________________________________________________________
function spm_gsranal(job)
    

    gradperc = .95;

    nsub = length(job.signals);
    
    %
    %load data
    %
    data = cell(nsub,1);
    trialdata = cell(nsub,1);
    nsess = zeros(nsub,1);
    fs = cell(nsub,1);
    
    total_loops = 0;
    for i=1:nsub
        nsess(i) = length(job.signals{i});
        subjdata = cell(nsess(i),1);
        subjtrialdata = cell(nsess(i),1);
        subjfs = cell(nsess(i),1);
        for j=1:nsess(i)
            sig = load(job.signals{i}{j});
            sig = sig.output;
            subjdata{j} = sig.signal;
            subjtrialdata{j} = sig.trialonsets;
            subjfs{j} = sig.fs;
        end;
        total_loops = total_loops + nsess(i);
   
        data{i} = subjdata;
        trialdata{i} = subjtrialdata;
        fs{i} = subjfs;
    end;
    
    %
    
    
    %
    % normalisation
    %
    
    switch job.normalisation
        case 1
            for i = 1:nsub
                for j = 1:nsess(i)
                    data{i}{j} = normalise(data{i}{j});
                end;
            end;
        case 2
            for i = 1:nsub
                data{i} = normalise(data{i});
            end;
        case 3
            data = normalise(data);
    end;
            
    %
    
    %
    % smoothing
    total_item = 0;
    spm_progress_bar('init');
    if job.smoothingw > 0;
        for i = 1:nsub
            for j = 1:nsess(i)
                  data{i}{j} = smoothy(data{i}{j}, job.smoothingw*fs{i}{j});
                total_item = total_item + 1;
                spm_progress_bar('set', total_item/total_loops);
            end;
        end;
        spm_progress_bar('clear');
    end;
    
    %
    
    
    %
    % downsampling
    %
    
    for i = 1:nsub
        for j = 1:nsess(i)
            if job.downsample < fs{i}{j}
                factor = fs{i}{j}/job.downsample;
                data{i}{j} = downsample(data{i}{j},factor);
                trialdata{i}{j}(:,1:2) = ceil(trialdata{i}{j}(:,1:2)/factor);
            end;
            fs{i}{j} = fs{i}{j} / factor;
        end;
    end;
        
    %
    
    
    %
    % SCR detection and quantification
    %
    
    actionloc = cell(nsub,1);
    heights = cell(nsub,1);
    humpdat = cell(nsub,1);
    for i = 1:nsub
         actionloc{i} = cell(nsess(i),1);
         heights{i} = cell(nsess(i),1);
         humpdat{i} = cell(nsess(i),1);
    end;
    %
%    if isfield(job.analchoice,'parameters') 
    %
    
   %     if job.analchoice.parameters.detect == 4
%         threshhold = pick_percentile(gradient(concatcell(trialdata),1/fs), gradperc);
%     end;
%     for i = 1:nsub
%        
%         if job.analchoice.parameters.detect == 3
%                threshhold = pick_percentile(gradient(concatcell(trialdata{i}),1/fs{i}), gradperc);
%         end;
%         for j = 1:nsess(i)            
%             if job.analchoice.parameters.detect == 1
%                 for k = 1:size(trialdata{i}{j},1)
%                     trialdat = data{i}{j}(trialdata{i}{j}(k,1):trialdata{i}{j}(k,2));
%                     threshhold = pick_percentile(gradient(trialdat,1/fs{i}{j}), gradperc);
%                     
%                     % Do the SCR detection
%                     % smoothingw may need to be a separate field for
%                     % derivative smoothing...
%                     [temptime, tempheight, temphump] = SCRdetect5(trialdat, threshhold, fs{i}{j}, ...
%                         job.smoothingw, job.analchoice.parameters.minperiod/1000, 0);
%                     
%                     if ~isempty(temptime)
%                         temptime = tempstart + trialdata{i}{j}(k,1);
%                         if ~isempty(temphump)
%                             temphump(:,1) = temphump(:,1) + trialdata{i}{j}(k,1);
%                             humpdat{i}{j} = [humpdat{i}{j}; temphump];
%                         end;
%                         actionloc{i}{j} = [actionloc{i}{j}; temptime];
%                         heights{i}{j} = [heights{i}{j}; tempheight];
%                     
%                     end;
%                 end;
%             else
%                 if job.analchoice.parameters.detect == 2
%                     threshhold = pick_percentile(gradient(data{i}{j},1/fs{i}{j}), gradperc);
%                 end
% 
%                 % Do the SCR detection
%                 % smoothingw may need to be a separate field for
%                 % derivative smoothing...
%                 [actionloc{i}{j}, heights{i}{j}, humpdat{i}{j}] = SCRdetect5(data{i}{j}, threshhold, fs{i}{j}, ...
%                     job.smoothingw, job.analchoice.parameters.minperiod/1000, 0);
% 
%             end
%         end
%     end
    
    %
%    end
    %
    
    
    
    %
    % set up display
    %
    
    artifacts = cell(nsub,1);
    for i = 1:nsub
        artifacts{i} = cell(nsess(i),1);
    end;
    
    
    % Initialize and save plotdata structure
    % This will store the data
    % ----------------------------------------
    plotdata.currentSession = 1;
    plotdata.currentSub = 1;
    plotdata.xloc = 1;
    plotdata.signal = data;
    plotdata.clean_signal = data; %until we've identified artifacts, so physio_find_SCRs works
    plotdata.trials = trialdata;
    plotdata.artifacts = artifacts;
    %     plotdata.SCR.location = actionloc;
    %     plotdata.SCR.height = heights;
    %     plotdata.SCR.humps = humpdat;
    plotdata.fs = fs;
    %plotdata.humpmin = job.analchoice.parameters.minhump;

    if isfield(job.analchoice, 'parameters'), plotdata.parameters = job.analchoice.parameters;
    else
        plotdata.parameters = [];
    end

    plotdata.parameters.smoothingw = job.smoothingw;
    plotdata.parameters.downsample = job.downsample;
    plotdata.parameters.normalisation = job.normalisation;

    create_figure('GSR Main Display');
    set(gcf, 'UserData', plotdata);

    
    if isfield(job.analchoice, 'parameters')
    
        plotdata.SCR = physio_find_SCRs(gcf, gradperc);
        set(gcf, 'UserData', plotdata);
        
    end
    
    
    %figure('Position',[scrsz(3)/3 1 scrsz(3)*2/3 scrsz(4)/2], 'UserData', plotdata, 'Name', 'Data display window');

    %clear prior plot windows
    oldaxes = findobj('UserData', 'mainPlot');
    if ~isempty(oldaxes)
        badwindow = get(oldaxes, 'Parent');
        close(badwindow);
    end;
    
    
    
    %top plot
    
    A = subplot(2,1,1);
    set(A,'Position', [.05, .85, .9, .1], 'UserData', 'topPlot');
    P = plot((1/fs{1}{1}):(1/fs{1}{1}):(length(data{1}{1})/fs{1}{1}), data{1}{1});
    
    
    
    % set up main plot
    
    B = subplot(2,1,2);
    set(B,'Position', [.1, .2, .8, .6], 'Userdata', 'mainPlot');
    
    setup_main_plot;
    
    set(P, 'ButtonDownFcn', 'zoom_in;');
    set(A, 'ButtonDownFcn', 'zoom_in;');
    
    
    % add buttons
    pos = get(gcf, 'Position');
    pos(1:2) = pos(3:4);
    
    uicontrol(gcf,'String','session ->',...
	'ToolTipString',...
	'Go to next session',...
	'Position',pos.*[.9 .10 .07 .05],...
	'CallBack','inc_session(1, gcbf);',...
	'ForegroundColor',[0 1 1])

    uicontrol(gcf,'String','<- session',...
	'ToolTipString',...
	'Go to previous session',...
	'Position',pos.*[.82 .10 .07 .05],...
	'CallBack','inc_session(-1, gcbf);',...
	'ForegroundColor',[0 1 1])
    

    uicontrol(gcf,'String','subject ->',...
	'ToolTipString',...
	'Go to next subject',...
	'Position',pos.*[.9 .03 .07 .05],...
	'CallBack','inc_subject(1, gcbf);',...
	'ForegroundColor',[0 1 1])

    uicontrol(gcf,'String','<- subject',...
	'ToolTipString',...
	'Go to previous session',...
	'Position',pos.*[.82 .03 .07 .05],...
	'CallBack','inc_subject(-1, gcbf);',...
	'ForegroundColor',[0 1 1])
    

    uicontrol(gcf,'String','Export Data',...
	'ToolTipString',...
	'Opens the data export Window',...
	'Position',pos.*[.1 .03 .12 .12],...
	'CallBack','export_data(gcbf);',...
	'ForegroundColor',[0 1 1])



    uicontrol(gcf,'String','Save Data',...
	'ToolTipString',...
	'Save your data',...
	'Position',pos.*[.25 .03 .12 .12],...
	'CallBack','save_data(gcbf);',...
	'ForegroundColor',[0 1 1])

    uicontrol(gcf,'String','Renormalize',...
	'ToolTipString',...
	'Renormalizes the current session, accounting for artifacts.',...
	'Position',pos.*[.40 .03 .12 .12],...
	'CallBack','renormalize(gcbf);',...
	'ForegroundColor',[0 1 1])
    
end







function output = normalise(input, m, s)
    
    if nargin == 1
        temp = concatcell(input);
        m = mean(temp);
        s = std(temp);
    end;
    
    if iscell(input)
        output = cell(length(input),1);
        for i = 1:length(input)
            output{i} = normalise(input{i},m,s);
        end;
    elseif isnumeric(input)
        output = (input - m)/s;
    else
        disp('data improperly formatted for normalisation');
    end;
end


function output = concatcell(incell)
    if isnumeric(incell)
        output = incell;
    elseif iscell(incell)
        if isnumeric(incell{1}) %is this the lowest level of the cell?
            
            output = [];
            for i = 1:length(incell)
                output = [output; incell{i}];
            end;
            
        else %not the lowest level, call self recursively
            output = [];
            for i = 1:length(incell)
                output = [output; concatcell(incell{i})];
            end;
        end;
    end
end

function outsignal = downsample(insignal, factor)
    minusside = ceil(-factor/2);
    plusside = ceil(factor/2)-1;
    
    outsignal = zeros(floor(length(insignal)/factor),1);
    j = 1;
    for i = (1-minusside):factor:(length(insignal)-plusside)
        outsignal(j) = mean(insignal(i+minusside:i+plusside));
        j = j + 1;
    end;
end

%funtion below is obsolete with current changes

% function value = pick_percentile(inarray, perc)
%     temp = sort(inarray);
%     pick_ind = perc * length(inarray);
%     
%     if pick_ind ~= floor(pick_ind)
%         floorco = pick_ind-floor(pick_ind);
%         value = temp(floor(pick_ind))*floorco + temp(ceil(pick_ind))*(1-floorco);
%     else
%         value = temp(pick_ind);
%     end;
% end
    


