function zoom_in()

    fig = get(gca, 'Parent');
    xpos = get(gca, 'CurrentPoint');
    set(gca, 'UserData', 'topPlot');
    
    %xlims = get(gca, 'XLim')
    xpos = xpos(1,1);
    
    %window = ceil((xlims(2)-xlims(1))/20)*2;
    
    %lowerx = max(xlims(1),xpos-window/2+1);
    %upperx = min(xlims(2),lowerx + window);
    %lowerx = upperx - window;
    
    %set(plotaxis, 'XLim', [lowerx upperx])
    
    plotdata = get(fig, 'UserData');
    plotdata.xloc = xpos;
    set(fig, 'UserData', plotdata);
    
    setup_main_plot
    
    
    
end
    