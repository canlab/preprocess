function obs = groupErrorBarPlotObserver(varargin)
    colormapfun = @summer;
    names = {};

    for i=1:2:length(varargin)
        propName = varargin{i};
        val = varargin{i+1};
        switch propName
            case {'Colors', 'colors'}
                colors = val;
            case {'Names', 'names'}
                names = val;
            otherwise
                error([propName ' is not a valid property']);
        end
    end
    
    obs = IVObserver(@bar_plot_);




    function bar_plot_(iv, mmPos, voxPos)
        numDataSets = get(iv, 'NumDataSets');
        for i=1:numDataSets
            data = getDataVolProp(iv, i, 'Ts');
            meanData(i) = mean(data);
            errorData(i) = std(data) / sqrt(length(data));
        end
        bar_handle = bar(meanData);
        
        numBars = length(meanData);
        colormap(colormapfun(numBars));
        s = get(bar_handle, 'Children');
        set(s, 'CData', 1:numBars); 
        
        hold on;
        errorbar(meanData, errorData, 'k+');
        
        %     fnames = get(iv, 'AssociatedVolsFilenames');
        %     [dummy fnames] = cellfun(@fileparts, fnames, 'UniformOutput', 0);
        %     set(gca,'XTickLabel', fnames);

        %title(sprintf('[x,y,z] = %3.0f, %3.0f, %3.0f', mmPos(1), mmPos(2), mmPos(3)), 'FontSize', 14);
        %axis tight;
        hold off;
    end

end