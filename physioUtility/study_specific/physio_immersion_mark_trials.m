% function physio_immersion_mark_trials(trigger, samprate, seconds, onsets_fname)
%
% trigger        vector from physio file 
% seconds        length of time to look at in one window
% mark the start of each trial on a plot
% onsets_fname         path of onsets file
%
% run this after run_physio_data_imm.m
%
% Matt Furstoss Sept 2009
% --------------------------------------
function [trials] = physio_immersion_mark_trials(trigger, samprate, seconds, onsets_fname)
% MARK TRIALS
        disp('Mark the start of each trial')
        disp('Ctrl+click to advance to next window of data')
        disp('Left click to mark start of first trial (Green lines).  The end of this trial and the start/end of all other trials is automatically calculated')
        disp('Shift+click to remove most recent line.')
        disp('Left+click again to end.')

    rbd = create_figure('Trials', 1, 1);
    axh(1) = subplot(1, 1, 1);
     
    hold on
    plot(trigger);
    axis auto
    axis tight
    
    % -Defaults
    data = [];
    start_looking_at = 0;
    
    % CHANGE VIEW (ZOOM IN)
    disp('Showing a user defined period of data.')
    viewsize = samprate * seconds;
    xlim = [start_looking_at start_looking_at+viewsize];
    set(axh, 'XLim', xlim);
    linkaxes(axh, 'x');
   % set(gca, 'YLim', [-2 2]);
    
        
   quitnow = 0;
%   startloc = 0; % no open ended start line exists
   green = []; % for plotting green lines
%   red = [];  % for plotting red lines
   starts = []; % list of starts
%   stops = [];  % list of stops
   nreg = 1; % number of complete start/end sets
   scalef = max(trigger);
   miny = 0;
    while ~quitnow
        data = guidata(rbd);
       % axes(axh(2));
        [xloc, yloc, button] = ginput(1);
        xloc = round(xloc);

        % Actions to perform based on user's mouse clicks
        if ~isempty(button)
            switch(button)
                case 1
                        starts(nreg) = xloc;
                        green = plot_onsets(xloc, 'g', miny, scalef);
                        nreg = nreg+1;

                case 2
                            starts(end) = [];
                            delete(green);
                            green = plot_onsets(starts, 'g', miny, scalef);
                            if nreg>1
                               nreg = nreg-1;
                            end
                case 3
                    xlim = get(axh, 'XLim');
                    stepsize = viewsize - 200;
                    data.last_start = xlim(1) + stepsize;
                    set(axh, 'XLim', xlim + stepsize);
            end
        else
            quitnow = 1;
            xlim = get(axh, 'XLim');
            fprintf('Last viewed range: %d %d\n', xlim(1), xlim(2));
                
        end

        guidata(rbd, data);
    end

% Import Onsets Data to calculate starts of all other trials.
% This file should contain the distances between each trail start
%     onsets = fileread(onsets_fname);
%     onsets = str2num(onsets);
%     onsets = onsets';
%     onsets = [starts onsets];
    
% Using the time between trials and the location of the first trial, calculate the actual time point of each
% subsequent trial.  Add the location of trial 1 to the distance between T1
% and T2 to find the actual location of T2.  Find T3 by adding T2 to the
% distance between T2 and T3.  Etc...
%     for i=2:length(onsets)
%         onsets(i) = onsets(i)+onsets(i-1);
%     end
    stops = starts;
    for i=1:length(stops)
        stops(i) = stops(i)+15000;
    end
% The end of each trial can be calculated by adding 15 seconds to the start
% time of each trial.
    trials = [starts' stops'];
%     for i=1:length(onsets)
%         trials{i} = (onsets(i):(onsets(i)+15000));
%     end

 % replot
    plot(trigger);
    for i=1:length(stops)
        plot_onsets(starts(i), 'g');
        plot_onsets(stops(i), 'r');
    end
    axis auto
    axis tight
    hold off
    %physio.trials = trials;
    
    end