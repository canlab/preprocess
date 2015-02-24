% [starts, stops] = physio_immersion_mark_trigger2(trigger, color, samprate)
% 
% Plots a graph allowing you to mark notable sections in the data
%
% trigger             trigger channel 
% color               choose color of stop lines
% name                name of plot
% samprate              sample rate in Hz, (e.g. 1000 for .001 sec interval)
% 
%
% Matt Furstoss, July 2009

function [starts, stops] = physio_immersion_mark_trigger2(trigger, color, name, samprate)
   
   quitnow = 0;
   startloc = 0; % no open ended start line exists
   green = []; % for plotting green lines
   red = [];  % for plotting red lines
   starts = []; % list of starts
   stops = [];  % list of stops
   nreg = 0; % number of complete start/end sets
   scalef = max(trigger);
   miny = 0;
    
    rbd = create_figure(name, 1, 1);
    axh(1) = subplot(1, 1, 1);
   
   % -Defaults
    data = [];
    start_looking_at = 0;
    disp('Clearing old onsets from figure data.');
    guidata(rbd, data);
    
    hold on
    plot(trigger);
    axis auto
    axis tight
    stop_color = 'r';
    
    if strcmp('red', color)
        stop_color = 'r';
    if strcmp('green', color)
        stop_color = 'g';
    else
        disp('Do not recognize color name, so stop markers will be default, red')
        color = 'r';
    end
    
    
    % CHANGE VIEW (ZOOM IN)
    disp('Showing a 120 s period of data.')
    viewsize = samprate * 120;
    xlim = [start_looking_at start_looking_at+viewsize];
    set(axh, 'XLim', xlim);
    linkaxes(axh, 'x');
   % set(gca, 'YLim', [-2 2]);
    
        
   
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
                        red = plot_onsets(xloc, stop_color, miny, scalef);
                        nreg = nreg+1;
                        startloc=0;
                    end
                case 2
                    if startloc==0
                        if ~isempty(stops)
                            stops(end) = [];
                            delete(red);
                            red = plot_onsets(stops, stop_color, miny, scalef);
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
end