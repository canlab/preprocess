% obs = groupTsErrorPlotObserver(varargin)
%
% groupTsErrorPlotObserver plots the mean timeseries across groups, with error bars around it

function obs = groupTsErrorPlotObserver(varargin)
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
            ts(i,:) = getDataVolProp(iv, i, 'Ts');
        end
        meanData = mean(ts);
        errorData = std(ts) / sqrt(length(ts));

        line_handle = plot(meanData);
        hold on;
        xVertices = [1:length(meanData) length(meanData):-1:1];
        yVertices = [meanData+errorData meanData(end:-1:1)-errorData(end:-1:1)];
        fill(xVertices, yVertices, get(line_handle, 'Color'), 'EdgeColor', 'none', 'FaceAlpha', .25);
        hold off;
    end

end