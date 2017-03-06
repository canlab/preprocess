% [preprocessed_files mean_image_name] = preproc_part1_2012(PREPROC, ...
%   ['slice timing', 0|1], ['motion correction', 0|1], ['local', 0|1], ...
%   ['verbose', 0|1], ['clean_up', 0|1], ['save plots', 0|1], ['SPM2' | 'SPM5'], ...
%   ['movie', 0|1], ['mean func', 0|1], ['run_dir_base', subpath])
%
% PREPROC is a structure containing these fields:
%   .func_files: cellstr of Nifti .nii or Analyze .img files containing the data
%       - use dicom2nifti() or dicom2analyze() to convert first, if needed
%   .TR: repeat time of scanner (essentially, time between each 3d volume)
%   .images_per_session: vector of the number of 3d volumes in each run - NO DISDAQS
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
% PREPROC.images_per_session = [120 120 120 120 120 120]; % DO NOT INCLUDE DISDAQS!!! REMOVE THEM YOURSELF OR WITH dicom2nifti / dicom2analyze!!!
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
%         More readable headers
%         Changed default for making_mean_image to 0 (too slow; not needed)
%         Changed default for save_plots to 1
%         Avoid using old/obsolete code
%         Added non-interactive mode option
%         Build SPM batch jobs from scratch, rather than relying on loading from .mat.
%         Save batch jobs to files
%         Show data plots for images post-realignment
%         Save additional output in PREPROC.Timeseries and PREPROC.Nuisance
% %           Added nuisance covs and higher-order transforms to PREPROC structure output
%         Eliminate extraneous plots for unused procedures (motion outliers...)
%         Movement parameters are collected and plotted in main canlab_preproc function
%             Using new subfunction, canlab_preproc_motion_covariates.m
%         Slice timing now set to 1st slice
%             - Was set to the middle slice, but most users probably not
%             adjusting onsets accordingly in model (e.g., microtime onset in model).
%             Therefore: default is now 1.
%         Improved realignment quality setting (from 0.9 to 1)

function PREPROC = preproc_part1_2012(PREPROC, varargin)
  spm('defaults', 'fmri');

  verbose = 0;
  save_plots = 1;
  clean_up = 1;               % this now moves files to Preprocessed directly; will not work right w/o this.
  spm_ver = 'SPM8';
  making_mean_image = 1;      % required for later steps - they will not work right w/o this
  correcting_slice_timing = 1;
  correcting_motion = 1;
  using_local_copies = 0;
  run_dir_base = '';

     % Parse the inputs and update default behaviors with user options
  parse_inputs();

  if ~isempty(run_dir_base)
    fprintf('Saving preprocessed images in:\n%s\n', run_dir_base)
  end

  spm_jobman('initcfg');

  switch(spm_ver)
    case{'SPM8', 'SPM12'}
      PREPROC = preproc_part1_spm8(PREPROC, correcting_slice_timing, correcting_motion, using_local_copies, ...
				   verbose, clean_up, save_plots, run_dir_base, acq);
      
    case {'SPM5', 'SPM2', 'SPM99'}
      print_header1('Warning! You are using an old version of SPM.', 'Use the original canlab_preproc for legacy pre-SPM8.');
      error('quitting now.')
  end

				% Mean image and diagnostic plots
  print_header1('Preproc_part1_2012', 'Diagnostic plots post_realignment');

  dat = fmri_data(char(PREPROC.ra_func_files{:}), PREPROC.implicit_mask_file);

			       % Write mean and save global timeseries
% -------------------------------------------------------------------------
  mm = mean(dat);
  mm.fullpath = PREPROC.meanfilename_realigned;
  write(mm);

  PREPROC.Timeseries.global_realigned = mean(dat.dat)';

				% Show/save plots of mean and data
% -------------------------------------------------------------------------
  savefilename = fullfile(PREPROC.qcdir, 'montage_mean_ra_func.png');
  canlab_preproc_show_montage(PREPROC.meanfilename_realigned, savefilename);

  plot(dat)
  fh = findobj('Type', 'Figure', 'Tag', 'fmri data matrix');
  if ~isempty(fh)
    sz = get(0, 'screensize');
    set(fh, 'Position', [50 sz(4)-50 sz(3)./2 sz(4)./2])
    delete(subplot(2, 3, 6));
    figure(fh);
    axh = subplot(2, 3, 5);
    set(axh, 'Position', [0.4108    0.1100    0.5134    0.3394]);
    axes(axh);
    title('Global values. circle size = spatial STD');
    scn_export_papersetup(500)
  end
  
  disp('First realigned functional of each run');
  savefilename = fullfile(PREPROC.qcdir, 'ra_func_first_volumes.png');
  canlab_preproc_show_montage(PREPROC.ra_func_files, savefilename)

  try_snapnow_for_publish;  % Data plot and mean heatmap

  if clean_up
    PREPROC = canlab_preproc_move_part1(PREPROC);
  end

				% End of main function


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
	  case 'SPM12'
            spm_ver = 'SPM12';
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
          case {'nomean'}
            making_mean_image = 0;
            
        end
      end
    end
    check_preproc(PREPROC);
  end % PARSE INPUTS


end % MAIN FUNCTION






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------------
%
% * Main function for SPM8
%
%
% --------------------------------------------------------------------------

function PREPROC = preproc_part1_spm8(PREPROC, correcting_slice_timing, correcting_motion, using_local_copies, verbose, clean_up, save_plots, run_dir_base, acq)

  mvmt_param_files = [];

  if(using_local_copies)
    disp('Copying run files to local machine.');
    run_files = transfer_to_newdir(PREPROC.func_files, tempdir());
  else
    run_files = PREPROC.func_files;
  end

		% Correct for different times of acquisition in slices
		% add a to filename
  a_files = spm8_timing_correction(PREPROC, run_files, correcting_slice_timing, verbose, clean_up, save_plots, run_dir_base, acq);
  if(using_local_copies), disp('Copying "a" files to remote location.'); transfer_to_newdir(a_files, pwd()); end

				% Motion correction
				% add r to filename
  [r_files, mvmt_param_files] = spm8_motion_correction(PREPROC, a_files, correcting_motion, using_local_copies, verbose, clean_up, save_plots, run_dir_base);
  if(using_local_copies), disp('Copying "ra" files to remote location.'); transfer_to_newdir(r_files, pwd()); end

  if(clean_up)
				% Remove a_files
    print_header2('ra_files created successfully. Removing a_files');
    for i = 1:length(a_files)
      if exist(a_files{i}, 'file')
        [dirname, fname] = fileparts(a_files{i});
        eval(['!rm ' fullfile(dirname, fname) '.*']);
      end
    end
  else
    PREPROC.a_func_files = a_files;
  end

  PREPROC.ra_func_files = r_files;
  PREPROC.mvmt_param_files = mvmt_param_files;

end


% -------------------------------------------------------------------------
% SLICE TIMING
% -------------------------------------------------------------------------

function a_files = spm8_timing_correction(PREPROC, func_files, correcting_slice_timing, verbose, clean_up, save_plots, run_dir_base, acq)

  if correcting_slice_timing
    print_header2('Slice-timing correction started.', acq);
    
    Vfirst_vol = spm_vol(func_files{1});
    num_slices = Vfirst_vol(1).dim(3);
    
    slice_timing_job{1}.spm.temporal.st.scans = spm8_image_list(PREPROC, func_files);
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
    slice_timing_job{1}.spm.temporal.st.refslice = 1; % round(num_slices/2);
    slice_timing_job{1}.spm.temporal.st.prefix = 'a';
    
				% save slice_timing_job
    savefile = fullfile(PREPROC.basedir, 'Functional', 'Preprocessed', 'canlab_preproc_slice_timing_job.mat');
    save(savefile, 'slice_timing_job');
    
    spm_jobman('run', {slice_timing_job});
    
    disp('Acquisition timing correction finished.');
  else
    disp('Skipping timing correction.');
  end

  a_files = prepend_a_letter(func_files, PREPROC.images_per_session, 'a', run_dir_base);

end

% -------------------------------------------------------------------------
% REALIGNMENT / MOTION CORRECTION
% -------------------------------------------------------------------------

function [r_files, mvmt_param_files] = spm8_motion_correction(PREPROC, a_files, correcting_motion, using_local_copies, verbose, clean_up, save_plots, run_dir_base)
  if(correcting_motion)
    print_header2('Motion correction started.');
    
			     % Build job ourself.
			     % load('motion_correction_spm5_job.mat');
    
    def = spm_get_defaults('realign');
    matlabbatch = {};
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions = def.estimate;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions = def.write;
    
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % do not register to mean (twice as long)
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 0; % do not mask (will set data to zero at edges!)
    
    matlabbatch{1}.spm.spatial.realign.estwrite.data = spm8_image_list(PREPROC, a_files);
    
    savefile = fullfile(PREPROC.basedir, 'Functional', 'Preprocessed', 'canlab_preproc_realignment_job.mat');
    save(savefile, 'matlabbatch');
    
    spm_jobman('run', {matlabbatch});
    
    disp('Motion correction finished.');
  else
    disp('Skipping motion correction.');
  end

  r_files = prepend_a_letter(a_files, PREPROC.images_per_session, 'r', run_dir_base);

  for i = 1:length(a_files)
    [pathstr filename ext] = fileparts(a_files{i}(1, :));
    mvmt_param_files{i} = fullfile(pathstr, ['rp_' filename '.txt']);
  end
  r_files = r_files(:);

  if(using_local_copies), transfer_to_newdir(mvmt_param_files, pwd()); end

  if(save_plots), save_plot(r_files{1}); end

end


% -------------------------------------------------------------------------
% FILE COPYING AND CHECKING UTILITY FUNCTIONS
% -------------------------------------------------------------------------

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
    if(~isfield(PREPROC, 'images_per_session'))
      errors{end+1} = 'Missing field "images_per_session"';
    end
    
    if(length(PREPROC.func_files) ~= 1 && length(PREPROC.func_files) ~= length(PREPROC.images_per_session) && length(PREPROC.func_files) ~= sum(PREPROC.images_per_session))
      errors{end+1} = 'Field "func_files" must be a cell array containing a) a single 4D volume for the whole experiment, b) a list of 4D volumes, one per run, or c) a list of 3D volumes for the whole experiment';
    end
  end

  if(~isempty(errors))
    error(['Errors encountered in PREPROC object:\n' implode(errors, '\n')], []);
  end
end


% -------------------------------------------------------------------------
% LIST IMAGES
% -------------------------------------------------------------------------

function run_image_list = spm8_image_list(PREPROC, file_list)
  num_runs = length(PREPROC.images_per_session);

  if(length(file_list) == num_runs)
    % EITHER cell array with full image list OR single 4-D image names
    
    for i = 1:num_runs
      if size(file_list{i}, 1) == PREPROC.images_per_session(i)
		 % we have the full list already % tor edit april 2010
        for j = 1:PREPROC.images_per_session(i)
          run_image_list{i}{j} = deblank(file_list{i}(j, :));
        end
      else
		% it's a single image name; expand it and create names
		% (could use expand_4d_images as well)
        printf_str = ['%' int2str(size(int2str(max(PREPROC.images_per_session)), 2)) 'd'];
        for j = 1:PREPROC.images_per_session(i)
          run_image_list{i}{j} = [file_list{i} ',' sprintf(printf_str, j)];
        end
      end
    end
    
  else
		     % another format; not cell array of length(nruns)
    st = [0 cumsum(PREPROC.images_per_session(1:end-1))];
    for i = 1:num_runs
      for j = 1:PREPROC.images_per_session(i)
        run_image_list{i}{j} = [file_list{st(i) + j} ',1'];
      end
    end
  end
end % function


% -------------------------------------------------------------------------
% MAKE MOVIE
% -------------------------------------------------------------------------

function generate_movie(img_files)

  Vfirst_func = spm_vol(img_files{1});
  ref_slice_num = round(Vfirst_func(1).dim(1) / 2);

  fprintf('Generating movie of slice %d.\n', ref_slice_num);
  movie_of_slice_timeseries(char(img_files), ref_slice_num, sprintf('expt-slice%d.avi', ref_slice_num), 'sagittal');

end

% -------------------------------------------------------------------------
% HEADER AND SAVE SUPPORT FUNCTIONS
% -------------------------------------------------------------------------

function print_header1(str, str2)

  s = '======================================================================================';
  len = length(s);

  disp('======================================================================================');
  disp('=                                                                                    =');
  fprintf('= %s%s=\n', str, repmat(' ', 1, max([1 length(s) - length(str) - 3])));
  if nargin > 1
    fprintf('= %s%s=\n', str2, repmat(' ', 1, max([1 length(s) - length(str2) - 3])));
  end
  disp('=                                                                                    =');
  disp('======================================================================================');


end


function print_header2(str, str2)

  s =  ' ______________________________________________________________________________________';
  s2 = '|                                                                                      |';
  s3 = '|______________________________________________________________________________________|';

  len = length(s);

  fprintf('%s\n%s\n', s, s2);

  fprintf('|%s%s|\n', str, repmat(' ', 1, max([1 length(s) - length(str) - 1])));

  if nargin > 1
    fprintf('|%s%s|\n', str2, repmat(' ', 1, max([1 length(s) - length(str2) - 1])));
  end

  fprintf('%s\n%s\n', s2, s3);


end



function try_snapnow_for_publish

  fhi = findobj('Type', 'Figure', 'Tag', 'Interactive');
  if ishandle(fhi), set(fhi, 'Visible', 'off'); end

  try
    snapnow
  catch
    warning('Error executing snapnow.  Old version of Matlab??');
  end


end


