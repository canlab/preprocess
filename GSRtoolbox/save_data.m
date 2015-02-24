function save_data(plotwin)

    plotdata = get(plotwin,'UserData');


    % Save vector of data will all bad times removed (marked NaN)
    for subj = 1:length(plotdata.signal)

        for sess = 1:length(plotdata.signal{subj})

            signal = plotdata.signal{subj}{sess};

            badtimes = round(plotdata.artifacts{subj}{sess});

            nbad = size(badtimes, 1);
            for j = 1:nbad, signal(badtimes(j, 1):badtimes(j, 2)) = NaN; end

            plotdata.clean_signal{subj}{sess} = signal;

        end

    end


    disp('Assigning variable plotdata in base workspace.');

    assignin('base', 'plotdata', plotdata);

    disp('Saving plotdata.mat');
    save('plotdata.mat', 'plotdata');

end


