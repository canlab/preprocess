% images = filenames('exp975/wtrial_height*.img', 'char');
% iv = InteractiveViewer('UseExistingGraphicsWindow', 1, 'AssociatedVols', images, 'IVObserver', tsPlotObserver());

function obs = tsPlotObserver()
    obs = IVObserver(@ts_plot_);
end

function ts_plot_(iv, mmPos, voxPos)
    ts = get(iv, 'CurrentVolsTs');
    plot(ts, 'Color', 'k');
    title(sprintf('[x,y,z] = %3.0f, %3.0f, %3.0f', mmPos(1), mmPos(2), mmPos(3)), 'FontSize', 14);
end