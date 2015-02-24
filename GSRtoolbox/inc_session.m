function inc_session(increment, mainfig)
    plotdata = get(mainfig, 'UserData');
    
    max_session = length(plotdata.signal{plotdata.currentSub});
    
    plotdata.currentSession = plotdata.currentSession + increment;
    
    if plotdata.currentSession < 1;
        plotdata.currentSession = 1;
    elseif plotdata.currentSession > max_session;
        plotdata.currentSession = max_session;
    end;
    
    set(mainfig, 'UserData', plotdata);
    
    
    fs = plotdata.fs{plotdata.currentSub}{plotdata.currentSession};
    
    A = findobj(mainfig, 'UserData', 'topPlot');
    n = length(plotdata.signal{plotdata.currentSub}{plotdata.currentSession});
    
    P = plot(A,((1/fs):(1/fs):(n/fs))', plotdata.signal{plotdata.currentSub}{plotdata.currentSession});
    set(P, 'ButtonDownFcn', 'zoom_in;');
    set(A, 'ButtonDownFcn', 'zoom_in;', 'UserData', 'topPlot');
    
    
    setup_main_plot;
    
end