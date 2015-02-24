function delete_rect()

    x = get(gca, 'CurrentPoint');
    x = x(1,1);
    plotdata = get(gcf, 'UserData');
    sess = plotdata.currentSession;
    subj = plotdata.currentSub;
    fs = plotdata.fs{subj}{sess};
    
    
    xind = round(x*fs);
    %kill boxes in the way
    %artifacts first
    boxestokeep = [];
    for i = 1:size(plotdata.artifacts{subj}{sess},1)
        if plotdata.artifacts{subj}{sess}(i,1) > xind || plotdata.artifacts{subj}{sess}(i,2) < xind
            boxestokeep = [boxestokeep; i];
        end;
    end
    plotdata.artifacts{subj}{sess} = plotdata.artifacts{subj}{sess}(boxestokeep,:);
    
    %then SCRs
    boxestokeep = [];
    for i = 1:size(plotdata.SCR.location{subj}{sess},1)
        if plotdata.SCR.location{subj}{sess}(i,1) > xind || plotdata.SCR.location{subj}{sess}(i,3) < xind
            boxestokeep = [boxestokeep; i];
        end;
    end;
    plotdata.SCR.location{subj}{sess} = plotdata.SCR.location{subj}{sess}(boxestokeep,:);
    plotdata.SCR.height{subj}{sess} = plotdata.SCR.height{subj}{sess}(boxestokeep,:);
    
    
    
    
    
    set(gcf, 'UserData', plotdata);
    setup_main_plot();
    
end
    