% [starts, stops] = physio_immersion_mark_blocks(trigger, samprate, seconds)
% 
% Plots a graph allowing you to mark notable sections in the data
%
% trigger             trigger channel 
% seconds             length of time to look at in one window
% name                name of plot
% samprate              sample rate in Hz, (e.g. 1000 for .001 sec interval)
% 
%
% Matt Furstoss, July 2009

function [blocks] = physio_immersion_mark_blocks(trigger, samprate, seconds)
    
% -------------------------------------   
% MARK BLOCKS
        disp('Mark the start\end of each block')
        fprintf('Left click to mark block start\\end (Green\\Red lines).  Your first of 2 clicks will mark start, second marks the end of a block')
        fprintf('Shift+click to remove most recent line.\nCtrl+click to move on to next series.\nCmd+cilck when done')

    rbd = create_figure('Blocks', 1, 1);
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
   startloc = 0; % no open ended start line exists
   green = []; % for plotting green lines
   red = [];  % for plotting red lines
   starts = []; % list of starts
   stops = [];  % list of stops
   nreg = 0; % number of complete start/end sets
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
                
        end

        guidata(rbd, data);
    end
    %replot
    plot(trigger);
    axis auto
    axis tight
    hold off
    
    if ~isempty(starts)
            if length(starts) ~= length(stops)
                stops(end+1) = data{wh_trigger}(end);
                fprintf('You did not mark an end point for the last block, so the end of the block will be the end of the trigger signal')
            end
        end
            
        blocks = [starts' stops']; % keeps track of start and end points 
%         blocks = {};
%         for i=1:length(stops)
%             blocks{i} = (markers(i,1):markers(i,2));
%         end
       %physio.blocks = blocks;
end