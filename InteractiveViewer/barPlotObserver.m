function obs = barPlotObserver()
    obs = IVObserver(@bar_plot_);
end


function bar_plot_(iv, mmPos, voxPos)
    data = get(iv, 'CurrentVolsTs');

    bar_handle = bar(data);
    %set(bar_handle,'FaceColor', [.7 .7 .7]);

    %     fnames = get(iv, 'AssociatedVolsFilenames');
    %     [dummy fnames] = cellfun(@fileparts, fnames, 'UniformOutput', 0);
    %     set(gca,'XTickLabel', fnames);

    %title(sprintf('[x,y,z] = %3.0f, %3.0f, %3.0f', mmPos(1), mmPos(2), mmPos(3)), 'FontSize', 14);
    %axis tight;
end