function a_files = prepend_a_letter(func_files, num_vols_per_run, char_to_add, varargin)
    % a_files = prepend_a_letter(func_files, num_vols_per_run, char_to_add)
    %
    % This function adds a letter (or other char string) to the front of each
    % name in a list of images, keeping the paths the same
    %
    % func_files is a cell array of a single OR expanded 4-D image list, one
    % cell per run, OR a cell array of 3-D image names for each run
    %
    % Tor Wager, Apr 2010
    % Made to work with preproc_part1

    run_dir_base = '';
    if ~isempty(varargin)
        run_dir_base = varargin{1};
    end
    
    for i = 1:length(func_files)

        % EITHER cell array with full image list OR single 4-D image names

        if size(func_files{i}, 1) == num_vols_per_run(i)
            % we have the full list already % tor edit april 2010

            a_files{i} = [];
            for j = 1:num_vols_per_run(i)

                [pathstr filename ext] = fileparts(func_files{i}(j, :));
                
                if ~isempty(run_dir_base)
                    % use this instead of pathstr; assume last dir in path
                    % is run dir and should be kept
                    [tmp, rundir] = fileparts(pathstr);
                    pathstr = fullfile(run_dir_base, rundir);    
                end
                
                myfilename = fullfile(pathstr, [char_to_add filename ext]);
                a_files{i} = strvcat(a_files{i}, myfilename);

            end

        else

            % it's a single image name; just add string
            % run_dir_base not implemented yet!
            [pathstr filename ext] = fileparts(func_files{i}(1, :));
            
            if ~isempty(run_dir_base)
                % use this instead of pathstr; assume last dir in path
                % is run dir and should be kept
                [tmp, rundir] = fileparts(pathstr);
                pathstr = fullfile(run_dir_base, rundir);
            end

            a_files{i} = fullfile(pathstr, [char_to_add filename ext]);

        end
    end

    % checks
    if diff(size(a_files)) > 0, a_files = a_files'; end
    for i = 1:length(a_files)
        a_files{i} = deblank(a_files{i});
    end

end % function