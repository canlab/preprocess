
function zoomIn(varargin)
    set(gca, 'XLimMode', 'manual', 'YLimMode', 'manual');

    xlimits = get(gca, 'XLim');   
    ylimits = get(gca, 'YLim');
    
    
    gwin = 400;
    
    linehandle = findobj(gca, 'Type', 'line');
    
    xvect = get(linehandle, 'XData');
    yvect = get(linehandle, 'YData');  
    point = get(gca,'CurrentPoint');
    
    xpos = round(point(1,1));
    
    if xpos > xlimits(1) && xpos < xlimits(2)
    
        xlb = max(xlimits(1)+1,xpos-gwin);
        xub = min(xpos+gwin,xlimits(2)-1);
    
        recthandle = findobj(gca, 'Type', 'rectangle');
        if isempty(recthandle) 
            recthandle = rectangle('Position',[xlb, ylimits(1),xub-xlb+1,ylimits(2)-ylimits(1)],'FaceColor','g', 'ButtonDownFcn', 'pollMouse;');
        else
            set(recthandle, 'Position', [xlb, ylimits(1),xub-xlb+1,ylimits(2)-ylimits(1)]);
        end;
     
    
        set(gca, 'Children', [linehandle recthandle]);
        
        xvect = xvect(xlb:xub);
        yvect = yvect(xlb:xub);
    
        subp = subplot(2,1,2);
        cla(subp)
        linehandle = plot(xvect,yvect);
        
        set(subp, 'ButtonDownFcn', 'beginDrag;', 'ButtonUpFcn', 'endDrag;')
        set(linehandle, 'ButtonDownFcn', 'beginDrag;', 'ButtonUpFcn', 'endDrag;')
        
    end;
end