function renormalize(datwin)

    gradperc = .95; %for now, hardcoded
    
    % plotdata has all the data in it, including 
    plotdata = get(datwin, 'UserData');
    plotdata = get_clean_signal(plotdata);
    
    subj = plotdata.currentSub;
    sess = plotdata.currentSession;
    
    rawsig = plotdata.signal{subj}{sess};
    
    
    arts = round(plotdata.artifacts{subj}{sess});
    
    if ~isempty(arts)
        arts = sort(arts,1);
        nart = size(arts);
        
        livesig = rawsig(1:(arts(1,1)-1));
        for i=2:nart(1)
            livesig = [livesig; rawsig((1+arts(i-1,2)):(arts(i,1)-1))];
        end
        
        livesig = [livesig; rawsig((arts(nart(1),2)+1):length(rawsig))];
        
    else
        livesig = rawsig;
    end
    
    rawsig = (rawsig - mean(livesig))/std(livesig);
        
    plotdata.signal{subj}{sess} = rawsig;

    
    
    %new stuff, april4 
    
    nsub = length(plotdata.signal);
    nsess = zeros(nsub,1);
    
    
    newSCRs = physio_find_SCRs(datwin, gradperc); 
    
    for i=1:nsub
        nsess(i) = length(plotdata.signal{i});
        for j=1:nsess(i)
            
    
            newSCRs.location{i}{j} = [plotdata.SCR.location{i}{j}; newSCRs.location{i}{j}];
            newSCRs.height{i}{j} = [plotdata.SCR.height{i}{j}; newSCRs.height{i}{j}];
            newSCRs.humps{i}{j} = [plotdata.SCR.humps{i}{j}; newSCRs.humps{i}{j}];
            
            keepers = clean_up(newSCRs.location{i}{j});          
            
            plotdata.SCR.location{i}{j} = newSCRs.location{i}{j}(keepers,:);
            plotdata.SCR.height{i}{j} = newSCRs.height{i}{j}(keepers,:);
            plotdata.SCR.humps{i}{j} = newSCRs.humps{i}{j}(keepers,:);
            
        end
    end
            
    
    % Make sure we don't duplicate SCRs or artifacts
    nsess = length(plotdata.signal{plotdata.currentSub});
    for i = 1:nsess
        scrs = plotdata.SCR.location{plotdata.currentSub}{i};
        scrs = unique(scrs, 'rows');
        plotdata.SCR.location{plotdata.currentSub}{i} = scrs;
        
        artifacts = plotdata.artifacts{plotdata.currentSub}{i};
        artifacts = unique(artifacts, 'rows');
        plotdata.artifacts{plotdata.currentSub}{i} = artifacts;
        
    end
    
    
    set(datwin, 'UserData', plotdata);
    setup_main_plot;
    
    
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

function keep_inds = clean_up(scr_locs)
    keep_inds = [];
    for i = size(scr_locs,1):-1:2 %throw out newests first
        if ~any((scr_locs(1:i-1,1)<scr_locs(i,1)&scr_locs(1:i-1,3)>scr_locs(i,1))|(scr_locs(1:i-1,1)<scr_locs(i,3)&scr_locs(1:i-1,3)>scr_locs(i,3))|(scr_locs(1:i-1,1)>scr_locs(i,1)&scr_locs(1:i-1,3)<scr_locs(i,3)))
            keep_inds = [keep_inds; i];
        end
    end;
end
        
        
        
    
    
    
    
    