function mouse_down_main()
    
    
    firstx = get(gca, 'CurrentPoint');
    firstx = firstx(1,1);
    plotdata = get(gcf, 'UserData');
    fs = plotdata.fs;
    sess = plotdata.currentSession;
    subj = plotdata.currentSub;
    
    
    [secondx, y, but] = ginput(1);
    
    if secondx < firstx, tmp = firstx; firstx = secondx; secondx = tmp; end
    
    %kill boxes in the way
    %artifacts first
    boxestokeep = [];
    for i = 1:size(plotdata.artifacts{subj}{sess},1)
        if plotdata.artifacts{subj}{sess}(i,1) < firstx*fs{subj}{sess} || plotdata.artifacts{subj}{sess}(i,1) > secondx*fs{subj}{sess}
            boxestokeep = [boxestokeep; i];
        end;
    end
    plotdata.artifacts{subj}{sess} = plotdata.artifacts{subj}{sess}(boxestokeep,:);
    
    %then SCRs
    boxestokeep = [];
    for i = 1:size(plotdata.SCR.location{subj}{sess},1)
        if plotdata.SCR.location{subj}{sess}(i,1) < firstx*fs{subj}{sess} || plotdata.SCR.location{subj}{sess}(i,1) > secondx*fs{subj}{sess}
            boxestokeep = [boxestokeep; i];
        end;
    end;
    plotdata.SCR.location{subj}{sess} = plotdata.SCR.location{subj}{sess}(boxestokeep,:);
    plotdata.SCR.height{subj}{sess} = plotdata.SCR.height{subj}{sess}(boxestokeep,:);
    
    
    
    
    
    if but == 1;
        [region_max max_ind] = max(plotdata.signal{subj}{sess}(ceil(firstx*fs{subj}{sess}):ceil(secondx*fs{subj}{sess})));
        max_ind = max_ind+ceil(firstx*fs{subj}{sess})-1;
        plotdata.SCR.location{subj}{sess} = [plotdata.SCR.location{subj}{sess}; [ceil(firstx*fs{subj}{sess}) max_ind ceil(secondx*fs{subj}{sess})]];
        plotdata.SCR.height{subj}{sess} = [plotdata.SCR.height{subj}{sess}; region_max-plotdata.signal{subj}{sess}(ceil(firstx*fs{subj}{sess}))];
        %plotdata.SCR.humps = [plotdata.SCR.humps; findhumps(plotdata.signal{subj}{sess}(ceil(firstx*fs{subj}{sess}):ceil(secondx*fs{subj}{sess})), ceil(firstx*fs{subj}{sess}), plotdata.humpmin)];
        
    else
        plotdata.artifacts{subj}{sess} = [plotdata.artifacts{subj}{sess}; [firstx*fs{subj}{sess} secondx*fs{subj}{sess}]];        
    end;
    
    set(gcf, 'UserData', plotdata);
    setup_main_plot();
    
end