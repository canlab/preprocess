function histogram_tool(x, varargin)
%function histogram_tool(x, ['bins', num bins], ['notitle'])

wh = find(strcmp(varargin, 'bins'));
if isempty(wh)
    bins = input('Histogram tool: Enter # of bins, or 0 if data are integers: ');
else
    bins = varargin{wh(1)} + 1;
end


%figure('Color','w'); set(gca,'FontSize',16)
create_figure('histogram');

if ~bins
    range = [min(x) max(x)];
    h = histc(x,range(1):range(2));
    han = barh(range(1):range(2),h);
else
    hist(x, bins)
    [h, xv] = hist(x, bins);
    set(gca,'XTick', unique(round(xv*10)./10));
end

drawnow;

if ~any(strcmp(varargin, 'notitle'))
    
    title(input('Enter title: ','s'))
    
    xlabel(input('Enter x-axis label: ','s'))
end

ylabel('Frequency')

end

