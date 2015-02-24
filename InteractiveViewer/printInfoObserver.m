function obs = printInfoObserver()
    obs = IVObserver(@print_info_);
end


function print_info_(iv, mmPos, voxPos)
    data = get(iv, 'CurrentVolsTs');

    fprintf('Values at %.1f %.1f %.1f:\n', mmPos(:));
    for i=1:length(data)
        fprintf('%.3f\n', data(i));
    end
    fprintf('\n');
end