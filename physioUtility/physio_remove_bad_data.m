% [missing, markers] = physio_remove_bad_data(pulse, old_missing, samprate)
% 
% Plots a graph allowing you to mark the start and end points of "bad" sections of data
%
% pulse,             the heart pulse data
% old_missing        vector of indices of missing data points
% samprate              sample rate in Hz, (e.g. 1000 for .001 sec interval)
%
%
% Matt Furstoss, July 2009

function [missing, markers] = physio_remove_bad_data(pulse, old_missing, samprate)
    rbd = create_figure('Pulse', 1, 1);
    axh(1) = subplot(1, 1, 1);
     
    hold on
    plot(pulse);
    axis auto
    axis tight
    
    % -Defaults
    data = [];
    missing = old_missing;
    start_looking_at = 0;
    
    % CHANGE VIEW (ZOOM IN)
    disp('Showing a 20 s period of data.')
    fprintf('Find regions of bad data and remove it.\nLeft click to mark start\\end of bad data (Green\\Red lines).')
    fprintf('Shift+click to remove most recent red\\green line.\nCtrl+click to move on to next series.\nCmd+cilck when done')
    viewsize = samprate * 20;
    xlim = [start_looking_at start_looking_at+viewsize];
    set(axh, 'XLim', xlim);
    linkaxes(axh, 'x');
   % set(gca, 'YLim', [-2 2]);
    
        
   quitnow = 0;
   startloc = 0; % no open ended start line exists
   green = []; % for plotting green lines
   red = [];  % for plotting red lines
   starts = []; % list of starts
   stops = [];  % list of stops
   nreg = 0; % number of complete start/end sets
   scalef = max(pulse);
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
                    if startloc==0  
                        startloc = 1;
                        starts(nreg+1) = xloc;
                        green = plot_onsets(xloc, 'g', miny, scalef);
                    else
                        stops(nreg+1) = xloc;
                        red = plot_onsets(xloc, 'r', miny, scalef);
                        nreg = nreg+1;
                        startloc=0;
                    end
                case 2
                    if startloc==0
                        if ~isempty(stops)
                            stops(end) = [];
                            delete(red);
                            red = plot_onsets(stops, 'r', miny, scalef);
                            startloc = 1;
                            nreg = nreg-1;
                        end    
                    else
                        if ~isempty(starts)
                            starts(end) = [];
                            delete(green);
                            green = plot_onsets(starts, 'g', miny, scalef);
                            startloc = 0;
                        end  
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
            % make markers matrix from starts and stops
            if length(starts) ~= length(stops)
                starts(end)=[];
            end
            % make sure pulse length doesn't change by adding unnecessary NaN's past it's normal length
            if stops(end)>length(pulse)
                stops(end) = length(pulse);
            end
            markers = [starts' stops']; % keeps track of start and end points 
                
            % remove regions from pulse, & append missing list
            if ~isempty(markers)
                badlist = []; % defining for use later
                for i = 1:nreg
                    pulse(markers(i,1):markers(i,2)) = NaN('double');
                    missing = [missing markers(i,1):markers(i,2)];
                end
            end
        end

        guidata(rbd, data);
    end
    %replot
    plot(pulse);
    axis auto
    axis tight
    hold off
end