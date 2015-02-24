% [preprocessed_files mean_image_name] = preproc_part1(PREPROC, ...
%   ['slice timing', 0|1], ['motion correction', 0|1], ['local', 0|1], ...
%   ['verbose', 0|1], ['clean_up', 0|1], ['save plots', 0|1], ['SPM2' | 'SPM5'], ...
%   ['movie', 0|1], ['mean func', 0|1], ['run_dir_base', subpath])
%
% PREPROC is a structure containing these fields:
%   .func_files: cellstr of Nifti .nii or Analyze .img files containing the data 
%       - use dicom2nifti() or dicom2analyze() to convert first, if needed
%   .TR: repeat time of scanner (essentially, time between each 3d volume)
%   .num_vols_per_run: vector of the number of 3d volumes in each run - NO DISDAQS
%
% Optional PREPROC fields:
%   .spm_ver: 'SPM2' or 'SPM5' (default: 'SPM5')
%   .run_dir_base - string of a subpath (or full path) under each subject's directory to place the functional files (default: '')
%
% Keyword parameters:
%   slice timing - boolean for whether or not to perform timing acquisition correction (default: 1)
%   motion correction - boolean for whether or not to perform motion correction (default: 1)
%   local - boolean for whether or not to copy files locally, operate, and then copy back (default: 0)
%   verbose - boolean for whether or not to display the shell commands (default: 0)
%   clean_up - boolean for saving plots of the results of various stages (default: 1)
%   save plots - boolean for saving plots of the results of various stages (default: 0)
%   SPM2 or SPM5 - flag to indicate which method to use (default: 'SPM5') - supercedes PREPROC.spm_ver
%   movie - boolean for whether or not to make an .avi movie of the middle slice (default: 0)
%   mean func - boolean for whether or not to make a mean image of all the functional images (default: 1)
%   acquisition - order of slice image acquisiiton - default interleaved, bottom up
%
% N.B.:
%   - Operates over network by default. If network is flaky, use 'local' option to somewhat alleviate.
%
% E.g.:
% PREPROC.func_files = {'run01/r01.nii' 'run02/r02.nii' 'run03/r03.nii' 'run04/r04.nii' 'run05/r05.nii' 'run06/r06.nii'}; % several 4d volumes, one per run - produced from stimulate.sdt files
%   or PREPROC.func_files = filenames('r*/vol*nii'); % list of all 3d volumes in experiment
%
% Or, for use with SPM2:
% PREPROC.func_files = {'run01/r01.img' 'run02/r02.img' 'run03/r03.img' 'run04/r04.img' 'run05/r05.img' 'run06/r06.img'}; % several 4d volumes, one per run - produced from stimulate.sdt files
%   or PREPROC.func_files = filenames('r*/vol*img'); % list of all 3d volumes in experiment
%
% PREPROC.TR = 2;
% PREPROC.num_vols_per_run = [120 120 120 120 120 120]; % DO NOT INCLUDE DISDAQS!!! REMOVE THEM YOURSELF OR WITH dicom2nifti / dicom2analyze!!!
% % To place run* dirs under a directory named 'Functional':
% PREPROC.run_dir_base = 'Functional'; % results in dirs like: 'subj03/Functional/run07', otherwise they take the form 'subj03/run07'
% save PREPROC PREPROC
%
% subjs = filenames('rea*');
% for i = 1:length(subjs)
%   cd(subjs{i});
%   preproc_part1(PREPROC);
%   cd('..');
% end
%
% % Or, if the denoised files are different for each subject, you could run something like this, 
% % to change the .func_files field for each subject:
% for i = 1:length(subjs)
%   cd(subjs{i});
%   PREPROC.func_files = filenames('*run*_denoised.img'); %or PREPROC.func_files = filenames('*run*_denoised.nii');
%   preproc_part1(PREPROC);
%   cd('..');
% end

% Notes:
% 8/2012  Tor changed default version to SPM8
%         Added minor improvements to output/error reporting
%         Changed default for making_mean_image to 0 (too slow; not needed)
%         Changed default for save_plots to 1
%         Added nuisance covs and higher-order transforms to PREPROC structure output

function PREPROC = preproc_part1(PREPROC, varargin)
    scn_setup();

    verbose = 0;
    save_plots = 1;
    clean_up = 1;
    spm_ver = 'SPM8';
    making_movie = 0;
    making_mean_image = 0;
    correcting_slice_timing = 1;
    correcting_motion = 1;
    using_local_copies = 0;
    run_dir_base = '';
    
    % Parse the inputs and update default behaviors with user options
    parse_inputs();

    if ~isempty(run_dir_base)
        fprintf('Saving preprocessed images in:\n%s\n', run_dir_base)
    end
    
    num_runs = length(PREPROC.num_vols_per_run);
    
    if(save_plots)
        if(length(PREPROC.func_files) <= num_runs)
            save_plot(PREPROC.func_files(:), 1);
        else
            save_plot(PREPROC.func_files{1}, 1);
        end
    end

    use_spm(spm_ver);
    spm_jobman('initcfg');
    
    switch(spm_ver)
        case{'SPM8'}
            PREPROC = preproc_part1_spm8(PREPROC, num_runs, correcting_slice_timing, correcting_motion, using_local_copies, ...
                verbose, clean_up, save_plots, run_dir_base, acq);
            PREPROC = review_motion_outliers(PREPROC, PREPROC.realigned_func_files, PREPROC.slice_time_corr_func_files, save_plots);
        case {'SPM5'}
            PREPROC = preproc_part1_spm5(PREPROC, num_runs, correcting_slice_timing, correcting_motion, using_local_copies, ...
                verbose, clean_up, save_plots, run_dir_base);
            PREPROC = review_motion_outliers(PREPROC, PREPROC.realigned_func_files, PREPROC.slice_time_corr_func_files, save_plots);
        case 'SPM2'
            PREPROC = preproc_part1_spm2(PREPROC, num_runs, correcting_slice_timing, correcting_motion, using_local_copies, run_dir_base, ...
                verbose, clean_up, save_plots);
    end
    
    if(making_movie)
        generate_movie(PREPROC.realigned_func_files);
    end

    if(making_mean_image)
        mean_image_name = generate_mean_func(PREPROC.realigned_func_files, run_dir_base);
        PREPROC.meanfilename_realigned = mean_image_name;
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inline functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function parse_inputs()
        if(isfield(PREPROC, 'spm_ver'))
            spm_ver = PREPROC.spm_ver;
        end
        if(isfield(PREPROC, 'run_dir_base'))
            run_dir_base = PREPROC.run_dir_base;
        end

        for i=1:length(varargin)
            if(ischar(varargin{i}))
                switch(varargin{i})
                    case {'verbose', 'display_commands', 'display commands'}
                        verbose = varargin{i+1};
                    case 'clean_up'
                        clean_up = varargin{i+1};
                    case {'save plots', 'save_plots'}
                        save_plots = varargin{i+1};
                    case 'SPM2'
                        spm_ver = 'SPM2';
                    case 'SPM5'
                        spm_ver = 'SPM5';
                    case 'SPM8'
                        spm_ver = 'SPM8';
                    case 'movie'
                        making_movie = varargin{i+1};
                    case {'local', 'copy local' 'use local copies'}
                        using_local_copies = varargin{i+1};
                    case {'slice_timing', 'slice timing'}
                        correcting_slice_timing = varargin{i+1};
                    case {'motion_correction', 'motion correction', 'correcting_motion', 'correcting motion'}
                        correcting_motion = varargin{i+1};
                    case {'mean', 'mean func', 'mean_func'}
                        making_mean_image = varargin{i+1};
                    case 'acquisition'
                        acq = varargin{i+1};
                end
            end
        end
        check_preproc(PREPROC);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------
% 
% * Main function for SPM5
%
% 
% -------------------------------------------

function PREPROC = preproc_part1_spm5(PREPROC, num_runs, correcting_slice_timing, correcting_motion, using_local_copies, verbose, clean_up, save_plots, run_dir_base)
    
    mvmt_param_files = [];
    if(using_local_copies)
        disp('Copying run files to local machine.');
        run_files = transfer_to_newdir(PREPROC.func_files, tempdir());
    else
        run_files = PREPROC.func_files;
    end
    if(save_plots), save_plot(run_files(:)); end


    % Correct for different times of acquisition in slices
    % add a to filename
    a_files = spm5_timing_correction(PREPROC, run_files, correcting_slice_timing, verbose, clean_up, save_plots, run_dir_base);
    if(using_local_copies), disp('Copying "a" files to remote location.'); transfer_to_newdir(a_files, pwd()); end

    % Motion correction
    % add r to filename
    [r_files, mvmt_param_files] = spm5_motion_correction(PREPROC, a_files, correcting_motion, using_local_copies, verbose, clean_up, save_plots, run_dir_base);
    if(using_local_copies), disp('Copying "ra" files to remote location.'); transfer_to_newdir(r_files, pwd()); end

    if(save_plots), save_plot(r_files{1}); end
    if(clean_up), disp('Do something for cleanup?'); end
    
    PREPROC.slice_time_corr_func_files = a_files;
    PREPROC.realigned_func_files = r_files;  
    PREPROC.mvmt_param_files = mvmt_param_files;
    
end

% -------------------------------------------
% 
% * Main function for SPM8
%
% 
% -------------------------------------------

function PREPROC = preproc_part1_spm8(PREPROC, num_runs, correcting_slice_timing, correcting_motion, using_local_copies, verbose, clean_up, save_plots, run_dir_base, acq)
    
    mvmt_param_files = [];
    if(using_local_copies)
        disp('Copying run files to local machine.');
        run_files = transfer_to_newdir(PREPROC.func_files, tempdir());
    else
        run_files = PREPROC.func_files;
    end
    if(save_plots), save_plot(run_files(:)); end


    % Correct for different times of acquisition in slices
    % add a to filename
    a_files = spm8_timing_correction(PREPROC, run_files, correcting_slice_timing, verbose, clean_up, save_plots, run_dir_base, acq);
    if(using_local_copies), disp('Copying "a" files to remote location.'); transfer_to_newdir(a_files, pwd()); end

    % Motion correction
    % add r to filename
    [r_files, mvmt_param_files] = spm5_motion_correction(PREPROC, a_files, correcting_motion, using_local_copies, verbose, clean_up, save_plots, run_dir_base);
    if(using_local_copies), disp('Copying "ra" files to remote location.'); transfer_to_newdir(r_files, pwd()); end

    if(save_plots), save_plot(r_files{1}); end
    if(clean_up), disp('Do something for cleanup?'); end
    
    PREPROC.slice_time_corr_func_files = a_files;
    PREPROC.realigned_func_files = r_files;  
    PREPROC.mvmt_param_files = mvmt_param_files;
    
end




% -------------------------------------------

function files_to_xfer = transfer_to_newdir(files_to_xfer, newdir)
    for i=1:length(files_to_xfer)
        [path image_file ext] = fileparts(files_to_xfer{i});
        switch(ext)
            case '.nii'
                cmd_string = sprintf('cp "%s".{nii,mat} "%s"', fullfile(path, image_file), newdir);
            case {'.hdr' '.img'}
                cmd_string = sprintf('cp "%s".{hdr,img,mat} "%s"', fullfile(path, image_file), newdir);
            otherwise
                cmd_string = sprintf('cp "%s" "%s"', files_to_xfer{i}, newdir);
        end
        
        [status, result] = system(cmd_string);
        if(status ~= 0)
            error(result);
        end
        files_to_xfer{i} = fullfile(newdir, files_to_xfer{i});
    end
end






% -------------------------------------------
% 
% * Main function for SPM2
%
% 
% -------------------------------------------


function PREPROC = preproc_part1_spm2(PREPROC, num_runs, correcting_slice_timing, correcting_motion, using_local_copies, run_dir_base, verbose, clean_up, save_plots)
    % Make into run files if needed
    if(length(PREPROC.func_files) ~= num_runs)
        run_files = create_run_files(PREPROC.func_files, PREPROC.num_vols_per_run, verbose, save_plots);
    else
        run_files = PREPROC.func_files;
    end

    if(save_plots), save_plot(run_files(:)); end


    % Correct for different times of acquisition in slices
    % add a to filename
    a_files = fsl_timing_correction(PREPROC, run_files, correcting_slice_timing, verbose, clean_up, save_plots);

    % Motion correction
    % add r to filename
    r_files = fsl_motion_correction(a_files, correcting_motion, verbose, clean_up, save_plots);

    %%% *** the line below seems like it might contain a bug; splits first
    %%% run only??
    r_files = split_4d_expt_analyze_img(r_files{1}, PREPROC, num_runs, run_dir_base, verbose);

    PREPROC.slice_time_corr_func_files = a_files;
    PREPROC.realigned_func_files = r_files;  
    PREPROC.mvmt_param_files = mvmt_param_files;
    
    if(save_plots), save_plot(filenames('run*/ravol0001.img')); end
    if(clean_up), delete_ana_imgs('ravols.img'); end
end

% -------------------------------------------

function check_preproc(PREPROC)
    errors = {};
    if(~isstruct(PREPROC))
        errors{end+1} = 'Variable PREPROC does not appear to be a structure';
    else
        if(~isfield(PREPROC, 'func_files'))
            errors{end+1} = 'Missing field "func_files"';
        end
        if(~isfield(PREPROC, 'TR'))
            errors{end+1} = 'Missing field "TR"';
        end
        if(~isfield(PREPROC, 'num_vols_per_run'))
            errors{end+1} = 'Missing field "num_vols_per_run"';
        end

        if(length(PREPROC.func_files) ~= 1 && length(PREPROC.func_files) ~= length(PREPROC.num_vols_per_run) && length(PREPROC.func_files) ~= sum(PREPROC.num_vols_per_run))
            errors{end+1} = 'Field "func_files" must be a cell array containing a) a single 4D volume for the whole experiment, b) a list of 4D volumes, one per run, or c) a list of 3D volumes for the whole experiment';
        end
    end

    if(~isempty(errors))
        error(['Errors encountered in PREPROC object:\n' implode(errors, '\n')], []);
    end
end

% -------------------------------------------

function run_files = create_run_files(func_files, num_vols_per_run, verbose, save_plots)
    disp('Creating run files.');

%     [func_path dummy dummy] = fileparts(func_files{1});
%     if(isempty(func_path))
        func_path = '.';
%     end

    if(length(func_files) == 1)
        split_command = sprintf('%s %s', ANA4DTO3D, func_files{1});
        if(verbose), disp(split_command); end
        unix(split_command);
        [func_path func_base func_ext] = fileparts(func_files{1}); %#ok
        func_files = filenames(fullfile(func_path, [func_base '_*.img']));
        if(save_plots), save_plot(func_files{1}); end
    end

    num_vols_per_run = num_vols_per_run(:);
    st = [1; cumsum(num_vols_per_run(1:end-1)) + 1];
    en = cumsum(num_vols_per_run);
    for i = 1:length(num_vols_per_run)
        image_list = func_files(st(i):en(i));
        
        run_files{i} = merge_vols(image_list, sprintf('%s/run%02d', func_path, i), 'verbose', verbose);
    end
    if(save_plots), save_plot(run_files{1}); end

    disp('Run file creation finished.');
end

% -------------------------------------------

function save_plot(img_files, new_step_num)
    persistent step_num

    img_files = cellstr(img_files);
    if(exist('new_step_num', 'var') && ~isempty(new_step_num))
        step_num = new_step_num;
    end

    for i = 1:length(img_files)
        spm_image('init', img_files{i});
        spm_print();
        %scn_export_papersetup(1000)

        [path, basename] = fileparts(img_files{i}(1,:));
        png_file = fullfile(path, sprintf('%02d-%s.png', step_num, basename));
        saveas(gcf, png_file);
    end

    step_num = step_num + 1;
    
    try
        snapnow
    catch
        disp('Error executing snapnow.  Old version of Matlab??');
    end

end

% -------------------------------------------

function a_files = spm5_timing_correction(PREPROC, func_files, correcting_slice_timing, verbose, clean_up, save_plots, run_dir_base)
    
    if(correcting_slice_timing)
        disp('Acquisition timing correction started.');

        load('slice_timing_spm5_job.mat');

        Vfirst_vol = spm_vol(func_files{1});
        num_slices = Vfirst_vol(1).dim(3);

        slice_timing_job.temporal{1}.st.nslices = num_slices;
        slice_timing_job.temporal{1}.st.tr = PREPROC.TR;
        slice_timing_job.temporal{1}.st.ta = PREPROC.TR - (PREPROC.TR / num_slices);
        slice_timing_job.temporal{1}.st.so = [1:2:num_slices 2:2:num_slices]; %might have to change
        slice_timing_job.temporal{1}.st.refslice = round(num_slices/2);

        slice_timing_job.temporal{1}.st.scans = spm5_image_list(PREPROC, func_files);

        spm_jobman('run', {slice_timing_job});

        disp('Acquisition timing correction finished.');
    else
        disp('Skipping timing correction.');
    end
   
    % Tor created this, April 2010
  	a_files = prepend_a_letter(func_files, PREPROC.num_vols_per_run, 'a', run_dir_base);
    
    if(save_plots), save_plot(a_files{1}); end
end

function a_files = spm8_timing_correction(PREPROC, func_files, correcting_slice_timing, verbose, clean_up, save_plots, run_dir_base, acq)
    if correcting_slice_timing
        disp('Acquisition timing correction started.');
        
        Vfirst_vol = spm_vol(func_files{1});
        num_slices = Vfirst_vol(1).dim(3);
        
        slice_timing_job{1}.spm.temporal.st.scans = spm5_image_list(PREPROC, func_files);
        slice_timing_job{1}.spm.temporal.st.nslices = num_slices;
        slice_timing_job{1}.spm.temporal.st.tr = PREPROC.TR;
        slice_timing_job{1}.spm.temporal.st.ta = PREPROC.TR - (PREPROC.TR / num_slices);
        switch acq
            case 'interleaved_BU'
                slice_timing_job{1}.spm.temporal.st.so = [1:2:num_slices 2:2:num_slices]; %interleaved acquisition bottom up
            case 'interleaved_TD'
                slice_timing_job{1}.spm.temporal.st.so = [num_slices:-2:1, num_slices-1:-2:1]; %interleaved acquisition top down
            case 'ascending'
                slice_timing_job{1}.spm.temporal.st.so = (1:1:num_slices);
            case 'descending'
                slice_timing_job{1}.spm.temporal.st.so = (num_slices:-1:1);
            case 'interleaved_MT'
                for k = 1:num_slices
                    slice_timing_job{1}.spm.temporal.st.so(k) = (round((num_slices-k)/2 + (rem((num_slices-k),2) * (num_slices - 1)/2)) + 1);
                end
            otherwise
                error('STOP! Unrecognized image acquisition method');
        end
        slice_timing_job{1}.spm.temporal.st.refslice = round(num_slices/2);
        slice_timing_job{1}.spm.temporal.st.prefix = 'a';
        
        spm_jobman('run', {slice_timing_job});
        
        disp('Acquisition timing correction finished.');
    else
        disp('Skipping timing correction.');
    end
    
    % SMS modified from spm5_timing_correction, June 2010
    a_files = prepend_a_letter(func_files, PREPROC.num_vols_per_run, 'a', run_dir_base);
    
    if(save_plots), save_plot(a_files{1}); end
end

% -------------------------------------------

function [r_files, mvmt_param_files] = spm5_motion_correction(PREPROC, a_files, correcting_motion, using_local_copies, verbose, clean_up, save_plots, run_dir_base)
    if(correcting_motion)
        disp('Motion correction started.');

        load('motion_correction_spm5_job.mat');

        motion_correction_job.spatial{1}.realign{1}.estwrite.eoptions.rtm = 0;
        motion_correction_job.spatial{1}.realign{1}.estwrite.data = spm5_image_list(PREPROC, a_files);

        spm_jobman('run', {motion_correction_job});

        disp('Motion correction finished.');
    else
        disp('Skipping motion correction.');
    end

    r_files = prepend_a_letter(a_files, PREPROC.num_vols_per_run, 'r', run_dir_base);
    
    for i = 1:length(a_files)
         [pathstr filename ext] = fileparts(a_files{i}(1, :));
%         r_files{i} = fullfile(pathstr, ['r' filename ext]);
        mvmt_param_files{i} = fullfile(pathstr, ['rp_' filename '.txt']);
    end
    r_files = r_files(:);
    
    if(using_local_copies), transfer_to_newdir(mvmt_param_files, pwd()); end

    if(save_plots), save_plot(r_files{1}); end
    
end

% -------------------------------------------

function a_file = fsl_timing_correction(PREPROC, run_files, correcting_slice_timing, verbose, clean_up, save_plots)
    global FSLDIR;
    
    if(correcting_slice_timing)
        disp('Acquisition timing correction started.');
        for i = 1:length(run_files)
            [run_path run_base run_ext] = fileparts(run_files{i});
            arun_files{i} = fullfile(run_path, ['a' run_base run_ext]);
            slicetimer_command = sprintf('export FSLDIR=%s && . %s/etc/fslconf/fsl.sh && FSLOUTPUTTYPE=ANALYZE && %s/bin/slicetimer --odd -r %s -i %s -o %s', ...
                FSLDIR, FSLDIR, FSLDIR, num2str(PREPROC.TR), run_files{i}, arun_files{i});
            if(verbose), disp(slicetimer_command); end
            unix(slicetimer_command);
        end
        
        if(save_plots), save_plot(arun_files); end
        merge_vols(arun_files, 'avols', 1);
        if(save_plots), save_plot('avols.img'); end
        if(clean_up), delete_ana_imgs(arun_files); end

        disp('Acquisition timing correction finished.');
    else
        disp('Skipping timing correction.');
    end

    a_file = 'avols.img';
end

% -------------------------------------------

function r_files = fsl_motion_correction(a_file, correcting_motion, verbose, clean_up, save_plots)
    global FSLDIR;

    if(correcting_motion)
        disp('Motion correction started.');
        mcflirt_command = sprintf('export FSLDIR=%s && . %s/etc/fslconf/fsl.sh && FSLOUTPUTTYPE=ANALYZE && %s/bin/mcflirt -in %s -out ravols.img -refvol 1 -rmsrel -rmsabs -plots', ...
            FSLDIR, FSLDIR, FSLDIR, a_file);
        if(verbose), disp(mcflirt_command); end
        unix(mcflirt_command);
        disp('Motion correction finished.');
    else
        disp('Skipping motion correction.');
    end
    
    r_files = {'ravols.img'};
    
    if(save_plots), save_plot('ravols.img'); end
    if(clean_up), delete_ana_imgs('avols.img'); end

    flip_if_neg_dim('ravols.img');
    if(save_plots), save_plot('ravols.img'); end
end

% -------------------------------------------

function split_files = split_4d_expt_analyze_img(img, PREPROC, num_runs, run_dir_base, verbose)
    global ANA4DTO3D;
    warning('off', 'MATLAB:MKDIR:DirectoryExists');

    split_files = {};
    mkdir('temp');
    split_command = sprintf('%s %s temp', ANA4DTO3D, img);
    if(verbose), disp(split_command); end
    unix(split_command);
    img_files = dir('temp/ravols*.img');
    files_copied = 0;
    for i = 1:num_runs
        run_dir = fullfile(run_dir_base, sprintf('run%02d', i));
        mkdir(run_dir);
        for j = 1:PREPROC.num_vols_per_run(i)
            files_copied = files_copied + 1;

            [d f e] = fileparts(img_files(files_copied).name); %#ok
            % NB: movefile() is insanely slow, do not use it for multiple files
            unix(sprintf('mv temp/%s %s/ravol%04d.hdr', [f '.hdr'], run_dir, j));
            unix(sprintf('mv temp/%s %s/ravol%04d.img', [f '.img'], run_dir, j));
            if(exist(fullfile('temp', [f '.mat']), 'file'))
                unix(sprintf('mv temp/%s %s/ravol%04d.mat', [f '.mat'], run_dir, j));
            end
            
            flip_if_neg_dim(sprintf('%s/ravol%04d.img', run_dir, j));
            split_files{end+1} = sprintf('%s/ravol%04d.img', run_dir, j);
        end
        fprintf('Moved %d files for run %d\n', PREPROC.num_vols_per_run(i), i);
    end
    rmdir('temp','s');
end

% -------------------------------------------

function run_image_list = spm5_image_list(PREPROC, file_list)
    num_runs = length(PREPROC.num_vols_per_run);

    if(length(file_list) == num_runs)
        % EITHER cell array with full image list OR single 4-D image names

        for i = 1:num_runs
            if size(file_list{i}, 1) == PREPROC.num_vols_per_run(i)
                % we have the full list already % tor edit april 2010
                for j = 1:PREPROC.num_vols_per_run(i)
                    run_image_list{i}{j} = deblank(file_list{i}(j, :));
                end
            else
                % it's a single image name; expand it and create names
                % (could use expand_4d_images as well)
                printf_str = ['%' int2str(size(int2str(max(PREPROC.num_vols_per_run)), 2)) 'd'];
                for j = 1:PREPROC.num_vols_per_run(i)
                    run_image_list{i}{j} = [file_list{i} ',' sprintf(printf_str, j)];
                end
            end
        end
        
    else
        % another format; not cell array of length(nruns)
        st = [0 cumsum(PREPROC.num_vols_per_run(1:end-1))];
        for i = 1:num_runs
            for j = 1:PREPROC.num_vols_per_run(i)
                run_image_list{i}{j} = [file_list{st(i) + j} ',1'];
            end
        end
    end
end % function

% -------------------------------------------

function generate_movie(img_files)
    Vfirst_func = spm_vol(img_files{1});
    ref_slice_num = round(Vfirst_func(1).dim(1) / 2);

    fprintf('Generating movie of slice %d.\n', ref_slice_num);
    movie_of_slice_timeseries(char(img_files), ref_slice_num, sprintf('expt-slice%d.avi', ref_slice_num), 'sagittal');
end

function mean_image_name = generate_mean_func(func_files, run_base_dir)
    [d f ext] = fileparts(func_files{1}(1, :));
    mean_image_name = ['mean_ra_func' ext];
    mean_image(char(func_files), mean_image_name, []);
end

function PREPROC = review_motion_outliers(PREPROC, r_files, a_files, save_plots)
    disp('Reviewing motion outliers');

    if(length(a_files) == sum(PREPROC.num_vols_per_run))
        st = cumsum([0 PREPROC.num_vols_per_run(1:end-1)]) + 1;
    else
        st = 1:length(a_files);
    end

    mvmt = [];
    for i = 1:length(st)
        [d f e] = fileparts(a_files{st(i)}(1,:)); %#ok
        mvmt_file = fullfile(d, ['rp_' spm_str_manip(f, 'r') '.txt']);
        current_mvmt = load(mvmt_file);
        mvmt = vertcat(mvmt, current_mvmt);
    end

    %     Vfunc = spm_vol(r_files{1});
    %     y = iimg_get_data(round(Vfunc(1).dim(1:3)/2), r_files);
    y = tor_global(char(r_files));
    PREPROC.global_realigned_values = y;

    OPT = scnlab_outlier_id('setup', 'tr', PREPROC.TR, 'spersess', PREPROC.num_vols_per_run, 'hp', 100, 'mad', 4, 'niter', 3, 'mvmt', mvmt);
    scnlab_outlier_id('data', y, 'options', OPT);
    subplot(4, 2, 3); title('Global mean of realigned imgs');
    
    if save_plots
        savedirname = fullfile(pwd, 'qc_images');
        if ~exist(savedirname, 'dir'), mkdir(savedirname); end
        scn_export_papersetup(600);
        saveas(gcf, fullfile(savedirname, 'movement_analysis.png'));
    end
end
