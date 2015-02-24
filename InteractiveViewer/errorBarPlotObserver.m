function obs = errorBarPlotObserver()
    obs = IVObserver(@bar_plot_);
end


function bar_plot_(iv, mmPos, voxPos)
    data = get(iv, 'CurrentVolsTs');
    errordata = (std(data) / sqrt(length(data))) * ones(size(data));
    bar_handle = bar(data);
    hold on;
    errorbar(data, errordata, '+');
    %set(bar_handle,'FaceColor', [.7 .7 .7]);

    %     fnames = get(iv, 'AssociatedVolsFilenames');
    %     [dummy fnames] = cellfun(@fileparts, fnames, 'UniformOutput', 0);
    %     set(gca,'XTickLabel', fnames);

    %title(sprintf('[x,y,z] = %3.0f, %3.0f, %3.0f', mmPos(1), mmPos(2), mmPos(3)), 'FontSize', 14);
    %axis tight;
    hold off;
end