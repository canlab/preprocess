
function setup_main_plot()
    
    plotaxis = findobj('UserData', 'mainPlot');
    plotfigure = get(plotaxis, 'Parent');
    plotdata = get(plotfigure, 'UserData');
    
    % ---------------------------------------    
    % Overall session window, Tor added
    % ---------------------------------------    
    % plotdata has all the data in it
%     plotdata = get_clean_signal(plotdata);
%     
%     cleansig_allsess = cat(1, plotdata.clean_signal{plotdata.currentSub}{:});
%     rawsig_allsess = cat(1, plotdata.signal{plotdata.currentSub}{:});
%     
%     acrossfh = create_figure('Across Sessions'); 
%     %title('All sessions together');
%     
%     
%     % cat SCRs etc
%     % *Could put much of this in SETUP!
%     allscrs = [];
%     allartifacts = [];
%     
%     sess_starts = 1;
%      nsess = length(plotdata.signal{plotdata.currentSub});
%     
%     for i = 1:nsess - 1
%         sess_length = length(plotdata.signal{plotdata.currentSub}{i});
%         sess_starts = [sess_starts cumsum(sess_starts) + sess_length];
%     end
%     
%     for i = 1:nsess
%         %handle(i) = plot_vertical_line(sess_starts(i));
%         
%         allscrs = vertcat(allscrs, plotdata.SCR.location{plotdata.currentSub}{i} + (sess_starts(i) - 1) );
%         allartifacts = vertcat(allartifacts, plotdata.artifacts{plotdata.currentSub}{i} + (sess_starts(i) - 1) );
%         
%     end
%     
%     boxhandles = findobj(acrossfh, 'Type', 'Patch');
%     delete(boxhandles);
%     
%     plothandles = draw_boxes(rawsig_allsess, allscrs, [0 .7 0]);
%     plothandles = [plothandles draw_boxes(rawsig_allsess, allartifacts, [.7 0 0])];
%     set(plothandles, 'EdgeColor', 'none', 'FaceAlpha', .4);
%     
%     yl = get(gca, 'YLim');
%     for i = 1:nsess
%         handle(i) = plot_vertical_line(sess_starts(i));
%         text(sess_starts(i), yl(2) .* 1.02, ['Session ' num2str(i)], 'FontSize', 18, 'FontWeight', 'b');
%             
%     end
%     set(handle, 'LineWidth', 3);
%     set(gca, 'YLim', [yl(1) yl(2) .* 1.1]);
%     
%     lineh(1) = plot(rawsig_allsess, 'k:'); 
%     lineh(2) = plot(cleansig_allsess, 'r');
%     set(gca, 'YLim', [min(rawsig_allsess) max(rawsig_allsess)]);
%     
%     legend(lineh, {'Original' 'Cleaned'});
%     drawnow
    
    % sig = signal vector
    % boxlocations = GSRloc, or plotdata.SCR{i}{j}.location, or
    % plotdata.artifacts{i}{j}
    
    % ---------------------------------------    
    % Selective Average Plot
    % --------------------------------------- 
    % onsets  offsets trialtype
%     alltrials = [];
%     for i = 1:nsess
%         sesstrials = plotdata.trials{plotdata.currentSub}{i};
%         sesstrials(:, 1:2) = sesstrials(:, 1:2) + (sess_starts(i) - 1); % make onsets/offsets relative to entire data vector
%         alltrials = vertcat(alltrials, sesstrials);        
%     end
%     
%     trialtypes = unique(alltrials(:, 3));
%     ntypes = length(trialtypes);
%     [averages, stderrs, trialdata] = deal(cell(1, ntypes));
%     
%     for i = 1:ntypes
%         trials = alltrials(alltrials(:, 3) == trialtypes(i), 1:2);
%         duration = min(diff(trials'));
%         
%         [averages(i), stderrs(i), trialdata(i)] = selective_average(cleansig_allsess, {trials(:, 1)}, 't', duration);
%     end
%     
%     plotdata.averages.trialtypes = trialtypes;
%     plotdata.averages.averages = averages;
%     plotdata.averages.stderrs = stderrs;
%     plotdata.averages.trialdata = trialdata;
%     
%     set(plotfigure, 'UserData', plotdata);
%     
%     % Now plot it
%     fs = plotdata.fs{plotdata.currentSub}{1};
%     
%     colors = {'b' 'r' 'g' 'k' 'm' 'c' 'y'};
%     while length(colors) < length(trialtypes), colors = [colors colors]; end
% 
%     create_figure('GSR Selective Averages');
%     for i = 1:ntypes
%         avg = plotdata.averages.averages{i};
%         xvals = linspace(0, length(avg) ./ fs, length(avg)); 
%         plot(xvals, avg, colors{i}, 'LineWidth', 2);
%         fill_around_line(avg, plotdata.averages.stderrs{i}, colors{i}, xvals);
%     end
%     
%     xlabel('Time (sec)');
%     
    
    % ---------------------------------------    
    % Matt's original windows
    % ---------------------------------------       
    
    plotfigure = create_figure('GSR Main Display', 1, 1, 1);
    set(plotfigure, 'CurrentAxes', plotaxis);
    
    plotaxis = findobj('UserData', 'mainPlot');
    %plotfigure = get(plotaxis, 'Parent');
    plotdata = get(plotfigure, 'UserData');
    plotsize = get(plotaxis, 'Position');
    figsize = get(plotfigure, 'Position');
    
    
    subj = plotdata.currentSub;
    sess = plotdata.currentSession;
    fs = plotdata.fs{subj}{sess};
    

    
    %
    % set axis limits
    %
    
    n = length(plotdata.signal{subj}{sess});
    
    window = ceil(n/20/fs)*2;
    
    lowery = min(plotdata.signal{subj}{sess});
    uppery = max(plotdata.signal{subj}{sess});
    upperyplus = (uppery-lowery)*.08+uppery;
    
    
    lowerx = max(1/fs,plotdata.xloc-window/2+1);
    upperx = min(n/fs,lowerx + window);
    lowerx = upperx - window;
    
    
    textdiffy = (upperyplus - uppery)/2;
    textdiffx = (upperx - lowerx) / 100;
    
    
    dotsizey = (uppery - lowery)/40;
    dotsizex = (upperx - lowerx)/40 / plotsize(3) / figsize(3) * plotsize(4) * figsize(4);
    
    
    %
    %draw the plot
    %
    elplot = plot((1:n)/fs,plotdata.signal{subj}{sess});
    
    set(plotaxis, 'XLimMode', 'manual', 'XLim', [lowerx upperx], 'YLimMode', 'manual', 'YLim', [lowery upperyplus]);

    
    %
    % draw trial info
    %
    for i=1:size(plotdata.trials{subj}{sess},1)

        rectpos = [plotdata.trials{subj}{sess}(i,1)/fs, uppery, (plotdata.trials{subj}{sess}(i,2)-plotdata.trials{subj}{sess}(i,1)+1)/fs,upperyplus-uppery];
        if rectpos(3) <= 0 || rectpos(4) <= 0
            warning('Invalid rectangle! This should not happen.');
        else

            rectangle('Position', rectpos,'FaceColor','y');
            text(textdiffx + plotdata.trials{subj}{sess}(i,1)/fs, uppery+textdiffy, num2str(plotdata.trials{subj}{sess}(i,3)));
        end

    end

    %
    % draw SCRs
    %
    for i=1:size(plotdata.SCR.location{subj}{sess},1)
        xind = plotdata.SCR.location{subj}{sess}(i,1);

        rectpos = [xind/fs, lowery, (plotdata.SCR.location{subj}{sess}(i,3)-xind+1)/fs,uppery-lowery];
        if rectpos(3) <= 0 || rectpos(4) <= 0
            warning('Invalid rectangle! This should not happen.');
        else

            rectangle('Position', rectpos,'FaceColor','g');
            xind = plotdata.SCR.location{subj}{sess}(i,2);
            text(textdiffx + plotdata.SCR.location{subj}{sess}(i,2)/fs, plotdata.signal{subj}{sess}(xind)+textdiffy, num2str(plotdata.SCR.height{subj}{sess}(i)));
        end

    end

    %
    % draw humps
    %

    for i=1:size(plotdata.SCR.humps{subj}{sess},1)
        xind = plotdata.SCR.humps{subj}{sess}(i,1);

        rectpos = [xind/fs-dotsizex/2,plotdata.signal{subj}{sess}(xind)-dotsizey/2,dotsizex,dotsizey];
        if rectpos(3) <= 0 || rectpos(4) <= 0
            warning('Invalid rectangle! This should not happen.');
        else
            rectangle('Position',rectpos, 'Curvature', [1 1], 'FaceColor', 'b');
        end

    end
    %
    % draw artifacts
    %

    for i=1:size(plotdata.artifacts{subj}{sess},1)
        xind = plotdata.artifacts{subj}{sess}(i,1);

        rectpos = [xind/fs, lowery, (plotdata.artifacts{subj}{sess}(i,2)-xind+1)/fs,uppery-lowery];
        if rectpos(3) <= 0 || rectpos(4) <= 0
            warning('Invalid rectangle! This should not happen.');
        else
            rectangle('Position', rectpos, 'FaceColor','r');
        end
    end


    %rearrange


    set(plotaxis, 'Children', [elplot; findobj(plotaxis, 'Type', 'text');  findobj(plotaxis, 'Type', 'rectangle')]);
    
    
    
    
    
    % callbacks!
    
    set(plotaxis, 'UserData', 'mainPlot', 'ButtonDownFcn', 'mouse_down_main');
    
    
    children = get(plotaxis, 'Children');
    for i = 1:length(children)
        set(children(i), 'ButtonDownFcn', 'mouse_down_main();');
    end;
    rects = findobj(plotaxis, 'Type', 'rectangle');
    for i = 1:length(rects)
        set(rects(i), 'ButtonDownFcn', 'delete_rect();');
    end;    
    rects = findobj(plotaxis, 'Type', 'rectangle', 'FaceColor', 'y');
    for i = 1:length(rects)
        set(rects(i), 'ButtonDownFcn', ['Collapse_Trial(' num2str(plotdata.trials{subj}{sess}(i,3)) ');']);
    end;    


    
    
end



function [plothandles, peakhandles] = draw_boxes(sig, boxlocations, colorvector)
    % sig = signal vector
    % boxlocations = GSRloc, or plotdata.SCR{i}{j}.location, or plotdata.artifacts{i}{j}
    
    plothandles = [];
    
    if isempty(boxlocations)
        % nothing to do
        return
        
    elseif size(boxlocations, 2) == 3
        startvals = boxlocations(:, 1);
        endvals = boxlocations(:, 3);
        peakvals = boxlocations(:, 2);

    elseif size(boxlocations, 2) == 2
        startvals = boxlocations(:, 1);
        endvals = boxlocations(:, 2);
        
    else
        error('I don''t understand the boxlocations you gave me.');
    end
        
    hold on;
    
    if exist('peakvals', 'var')
        peakhandles = plot(peakvals, sig(peakvals), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10);
    end
    
    nGSR = size(boxlocations, 1);

    ylim = get(gca, 'YLim');
    clear h1
    for i = 1:nGSR

        plothandles(i) = drawbox(startvals(i), endvals(i) - startvals(i), ylim(1), ylim(2) - ylim(1), colorvector);

    end

    %plothandles = [plothandles h1];

end




function plotdata = get_clean_signal(plotdata)

    % Save vector of data will all bad times removed (marked NaN)
    for subj = 1:length(plotdata.signal)

        for sess = 1:length(plotdata.signal{subj})

            signal = plotdata.signal{subj}{sess};

            badtimes = round(plotdata.artifacts{subj}{sess});

            nbad = size(badtimes, 1);
            for j = 1:nbad, signal(badtimes(j, 1):badtimes(j, 2)) = NaN; end

            plotdata.clean_signal{subj}{sess} = signal;

        end

    end

end
