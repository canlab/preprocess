function gsr_display(signal, overlay, tlabels, tspots)
    

    n = length(overlay);
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]); %, 'UserData', num2str(job.dispwin));
    A = subplot(2,1,1);
    
    P = plot(signal);
    drawnow
    
    
    set(A, 'XLimMode', 'manual', 'YLimMode', 'manual');

    xlimits = get(gca, 'XLim');   
    ylimits = get(gca, 'YLim');
    
    rectdims = find(overlay(1:n-1)~=overlay(2:n));
    if overlay(1) == 1
        rectdims = rectdims(2:length(rectdims));
    end;
    if overlay(length(overlay))==1
        rectdims = [rectdims; length(overlay)];
    end;
    
    rearrstr = '';
    for i = 1:2:length(rectdims)
        rectnum = num2str((i+1)/2);
        eval(['rect' rectnum ' = rectangle(''Position'',[rectdims(i), ylimits(1),rectdims(i+1)-rectdims(i),ylimits(2)-ylimits(1)],''FaceColor'',''g'', ''UserData'', tspots{str2num(rectnum)},  ''ButtonDownFcn'', ''rectZoom(' rectnum ');'');']);
        eval(['txt' rectnum ' = text(rectdims(i),ylimits(1)+(ylimits(2)-ylimits(1))/20,tlabels{str2num(rectnum)});']);
        rearrstr = [rearrstr ' txt' rectnum ' rect' rectnum];
    end;
    
    eval(['set(A, ''Children'', [P ' rearrstr ']);']);

end
        
