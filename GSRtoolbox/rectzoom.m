function rectzoom(rectnum)


    P = get(gca, 'children');
    linehandle = P(1);
    
    xvect = get(linehandle, 'XData');
    yvect = get(linehandle, 'YData');
   
    rectpos = get(P(1+rectnum*2), 'Position');
    
    xvect = xvect(rectpos(1):(rectpos(1)+rectpos(3)));
    yvect = yvect(rectpos(1):(rectpos(1)+rectpos(3)));
    
    subp = subplot(2,1,2);
    cla(subp);
    plot(xvect,yvect);
    
    items = get(P(1+rectnum*2), 'UserData');
    for i=1:length(items)
        text(items(i), 0, 'GSR');
    end;
    
    
end
    