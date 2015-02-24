function inc_subject(increment, mainfig)
    plotdata = get(mainfig, 'UserData');
    
    max_subject = length(plotdata.signal);
    
    plotdata.currentSub = plotdata.currentSub + increment;
    
    if plotdata.currentSub < 1;
        plotdata.currentSub = 1;
    elseif plotdata.currentSub > max_subject;
        plotdata.currentSub = max_subject;
    else
        plotdata.currentSession = 1;
    end;
    
    set(mainfig, 'UserData', plotdata);
    
    %top plot
    fs = plotdata.fs{plotdata.currentSub}{plotdata.currentSession};
    
    
    A = findobj(mainfig, 'UserData', 'topPlot');
    P = plot(A,[(1/fs):(1/fs):(length(plotdata.signal{plotdata.currentSub}{plotdata.currentSession})/fs)]', plotdata.signal{plotdata.currentSub}{plotdata.currentSession});
    set(P, 'ButtonDownFcn', 'zoom_in;');
    set(A, 'ButtonDownFcn', 'zoom_in;', 'UserData', 'topPlot');
    
    
    
    setup_main_plot;
    
end