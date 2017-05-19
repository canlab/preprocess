function canlab_preproc_2012(basedir, run_wildcard, image_wildcard, TR, struct_wildcard, varargin)
% canlab_preproc_2012(basedir, run_wildcard, image_wildcard, TR, struct_wildcard, [options])
%
% This is the basic preprocessing function for fmri images.  It uses
% preproc_part1 and preproc_part2.  It will output files into a
% directory based on the child directory from basedir.  If images are
% stored as 4-D image files, make sure that only one image file is in each
% run directory, otherwise the code will assume that the images are
% actually 3-D.
%
%
% INPUT:
%   basedir:
%      the working directory that is common to all the necessary
%      files.  With the suggested set-up, this should be the
%      the subject's base directory directory
%   run_wildcard:
%      a string that indicates the set [*] of directories that
%      hold all the raw data to be processed. Must be located
%      in basedir/Functional/Raw/[*].  If there are multiple wilcards
%      in your run directory listing, make a char array where each row
%      is one of the required wildcards.  This feature only works with
%      up to two rows at the moment.
%      For example, ['run*'; 'APP*'] Note that both elements must be
%      the same length!  You can pad with spaces if necessary (i.e.
%      ['run*', 'A*  '], but do not use '*    ' -- you must have at
%      least one non-wildcard character or it will take the . and ..
%      directories as well, which is bad.
%   image_wildcard:
%      a string that indicates the raw files that should be processed
%   TR:
%      the TR of the scan, in seconds.
%   struct_wildcard:
%      a string that indicates the the set of structural
%      images.  can have two wildcards, i.e. '/APP*/*.nii'
%
% OPTIONS:
%   'nodatachecks'  Skip initial data checks/plots
%   'nomovie':      no movie file is created
%   'nomean':       no mean image mask for spike analysis is generated
%   'nospikeid':	the spike analysis is not run
%   'nopart1':      preproc_part1 is not run. This assumes that it has
%                      been run in the past and that the *.mat file is
%                      located in the proper location so that it may be
%                      loaded properly.
%   'nopart2':      do not run part 2
%
%   'nosliceorder': ?????
%   'disdaqs':      the number of functional images at the beginning of
%                      a scan to discard.  The string should be followed by
%                      either a single integer that represents the number
%                      of disdaqs to remove from each run, or a column
%                      array that specifies how many disdaqs to remove from
%                      each individual run.  Removing disdaqs for 4D images
%                      requires FSL to run. Default value of disdaqs is 0.
%
%         **WARNING:** if disdaqs already happened, and you
%                      are re-running the data from a midway point, do NOT
%                      do disdaqs again!
%
%   'acquisition':  used to specify the acquisition type in the
%                      following string.  Possible choices are:
%
%      'interleaved_BU': Image acquisition was interleaved, bottom-up (Default)
%      'interleaved_TD': Image acquisition was interleaved, top-down
%      'ascending':      Image acquisition was in order, bottom to top
%      'descending':     Image acquisition was in order, top to bottom
%      'interleavedMT':  Image acquisition was interleaved, middle - top
%
% Part 2 options (these will be passed into preproc_part2_2012
% {'verbose' 'display_commands' 'display commands'}, verbose = varargin{i+1};
% {'clean_up' 'clean up'}, clean_up = varargin{i+1};
% {'save_plots' 'save plots'}, save_plots = varargin{i+1};
% {'set_origins' 'set origins'}
% setting_origins = varargin{i+1};
% {'coregistering_anat_to_func' 'coregistering anat to func' 'coreg anat to func'}
% coregistering_anat_to_func = varargin{i+1};
% {'warp' 'warping'}
% warping = varargin{i+1};
% {'generating_mean' 'generating mean' 'generate mean'}
% generating_mean = varargin{i+1};
% {'smooth' 'smoothing'}
% smoothing = varargin{i+1};
% {'check_norms' 'check norms'}
% check_norms = varargin{i+1};
% 'nointeractive'
% dointeractive = 0;
% 'noorigin'
% setting_origins = 0;
% 'nocoreg'
% coregistering_anat_to_func = 0;
% 'nowarp'
% warping = 0;
% 'nomean'
% generating_mean = 0;
% 'nosmooth'
% smoothing = 0;

% -------------------------------------------------------------------------
% OUTPUT:
%
% Diagnostic files are output to qc_dir.
% Preprocessed files are output
%    to basedir/../Preprocessed.
%
% Images include:
%  * Implicit mask histogram:  intensity values of the image and the
%  implicit mask.  The values of the first 20 images are averaged and then
%  displayed in this histogram.  This is created in the function
%  fmri_mask_thresh_canlab.m.  Voxels with low values are likely out-of-brain
%  and will be excluded from the implicit mask.
%  * Implicit mask montage:  a picture of the implicit mask
%  * Mean mask before preproc:   This is the (fairly broad) mask used for
%  preprocessing.  It is a binary mask produced by taking all voxels whose
%  mean signal is in the top 80% of mean signal values (i.e., applying a
%  threshold of 80% to mean_func_before_preproc.img).  you can use
%  spm_check_registration(<filename>) to view this image.
%  * Motion params:  in the top left subplot, each line represents movement
%  in one direction (X, Y, or Z).
%  * Scn_session_spike:  this figure is generated by
%  scn_session_spike_id.m.  In the color bars, each point represents the
%  mean intensity or STD for a given slice at a given time/TR.  Values are
%  z-scored within slice.  A spike may be identifiable if a given image
%  (i.e. vertical column) is several standard deviations from other images.
%  Spikes are usually identified at the level of image rather than slice.
%  The bar graphs are generated by the function trimts.m, see the help
%  there for more information.  Each red vertical line representes an
%  identified spike
%  * Sources of variance:  This is generated in the inline function
%  summarize_multisession_output in scn_session_spike_id.m.  The total
%  variance in global signal is divided into spikes, linear trends, and
%  within-session error.
%  * FFT of unexplained global signal:  This is generated in the inline
%  function summarize_multisession_output in scn_session_spike_id.m.  It is
%  the FFT of the residual error, i.e. variance not due to spikes or linear
%  drift
%
% Examples:
% -------------------------------------------------------------------------
% canlab_preproc('/data/projects/wagerlab/labdata/current/PAL/Imaging/PALP2324', 'r*', '*run*nii', 'IM*dcm', 2);
%
% Notes:
% To clean up a partially run processing job, try:
% 1) load the PREPROC_SETUP file
% 2) PREPROC = canlab_preproc_clean_up_and_move_files(PREPROC);
% 3) allerrors = canlab_preproc_check_PREPROC(PREPROC);

% Changes in 2012 scripts:
% DIAGNOSTICS AND PLOTS
% More diagnostic plots, including plots from fmri_display object
% one implict mask and one mean image only during spike artifact detection
% New movie of successive difference images
% New additional outlier detection based on RMSSD (root mean sq. successive
% diffs) of images.
% Saves motion parameter estimates and higher-order transformations in
% PREPROC.Nuisance and ***** file compatible with SPM batch mode.
%
% PROCESS CONTROL
% Now make non-interactive mode, compatible with publish.m
%  - also allows re-running of multiple subjects without interactive input
%  - particularly OK if origins and coreg (main interactive steps) are
%  already done
% Builds and saves SPM8-compatible SPM batch jobs
% Returns errors for SPM2/5
%
% READABILITY
% More consistent use of names in PREPROC structure and various
% subfunctions
% Uses existing object-oriented functions to reduce amount of code
%
% TO-DO
% -- interpolate bad outliers
% -- hide SPM windows before snapnow - progress bar
% ventricle timeseries
%
% Fix (omit) mean writing in part 2
% Publish all after initial, within script?  Do interactive first...
% I think is flipping the images - fix that.
% Fix reliance on saved jobs in normalize/smooth
% Set header on VOLS so we can run in non-interactive mode after that.
%
% Temporal SNR images - montage plot
% Add successive diffs in images to band plot of globals, etc.
%
% * Move the origin and coreg parts to separate function.
% * put setup and checks for whether dirs and files exist in
% canlab_preproc_list_files.  that should be separate, and set up -- then
% next part runs.

%% ------------------------------------------------------------------------
% Defaults
% -------------------------------------------------------------------------

dodatachecks = 1;
writemovie = 1;
writemean = 1;
dospikeid = 1;
dopart1 = 1;
dopart2 = 1;
disdaqs = 0;
acq = 'interleaved_BU';
dointeractive = 1;
clean_up = 1;

meanfilename_realigned = []; % not called otherwise, so we need to here.
qcdir = [];
meanfilename = [];

setup_inputs;

PREPROC = struct();

if exist(preprocmatname, 'file')
    print_header2('Loading existing PREPROC structure and adding/overwriting fields.');
    load(preprocmatname, 'PREPROC');
end

% Replace fields in PREPROC based on inputs
myvar = [];
for i = {'basedir', 'logfilename', 'run_wildcard', 'image_wildcard', 'qcdir', 'TR', 'meanfilename', 'meanfilename_realigned'} %, 'meanmaskname'}
    
    eval(['myvar = ' i{1} ';']);
    PREPROC.(i{1}) = myvar;
    clear myvar
    
end

PREPROC.run_dir_base = []; % preprocdir;  % could do this to move files, but need to debug in preproc_part1

%% -------------------------------------------------------------------------
% List files and print
% Get the anatomical ready for part 2 if requested
%
% Expand 4-D image filenames for SPM analysis if needed
% Check that files exist
% Move discarded acquisitions to subfolder, if requested
% -------------------------------------------------------------------------

% Check to see what is already done
print_header1('Checking existing PREPROC structure', 'and setting up file names')

[errors, ISDONE] = canlab_preproc_check_PREPROC(PREPROC, 'erroroff');


if ~ISDONE.func_files
    
    [PREPROC.func_files, PREPROC.rundirs, PREPROC.images_per_session, PREPROC.num_disdaqs] = canlab_preproc_list_files(basedir, run_wildcard, image_wildcard, disdaqs);
    diary off
    
end

if dopart2 && ~ISDONE.anat_files
    
    PREPROC.anat_files = prep_structural_files(basedir, struct_wildcard);
    
    if ~exist(PREPROC.anat_files,'file')
        fprintf('Anatomical file specified in PREPROC.anat_files does not exist...Looking for:\n%s\n. Quitting.\n', PREPROC.anat_files);
        return
    end
    
end

try_snapnow_for_publish;

%
% % Add to fields in PREPROC based on file search
% for i = {'rundirs', 'func_files', 'images_per_session', 'num_disdaqs', 'anat_files'};
%
%     eval(['myvar = ' i{1} ';']);
%     PREPROC.(i{1}) = myvar;
%     clear myvar
%
% end
%
% if dopart2
%     if ~exist(PREPROC.anat_files,'file')
%         fprintf('Anatomical file specified in PREPROC.anat_files does not exist...Looking for:\n%s\n. Quitting.\n', PREPROC.anat_files);
%         return
%     end
% end

% Run this part always
[dummy, dummy, dummy, dummy, outputname] = fmri_mask_thresh_canlab(char(PREPROC.func_files{:}), 'implicit_mask.img');
PREPROC.implicit_mask_file = fullfile(pwd, outputname);

try_snapnow_for_publish;  % Mask figures

%% -------------------------------------------------------------------------
% Preliminary data check
% Create and show implicit mask
% Vizualize in-brain data
% Spike (transient) detection and global signal
% Save in PREPROC.Nuisance and PREPROC.Timeseries
% -------------------------------------------------------------------------
if dodatachecks
    
    %     [dummy, dummy, dummy, dummy, outputname] = fmri_mask_thresh_canlab(char(PREPROC.func_files{:}), 'implicit_mask.img');
    %     PREPROC.implicit_mask_file = fullfile(pwd, outputname);
    %     % used later for mean image, etc. - moved later by canlab_preproc_move_part1
    %
    %     try_snapnow_for_publish;  % Mask figures
    
    print_header1('Showing all fMRI data before preprocessing.')
    dat = fmri_data(char(PREPROC.func_files{:}), 'implicit_mask.img');
    
    plot(dat)
    
    disp('First raw functional of each run');
    savefilename = fullfile(PREPROC.qcdir, 'func_first_volumes.png');
    canlab_preproc_show_montage(PREPROC.func_files, savefilename)
    
    
    try_snapnow_for_publish;  % Data plot and mean heatmap
    
    diary(logfilename)
    
    print_header2('Outliers based on mahalanobis dist of global vals/spatial STD/RMSSD.')
    dat.images_per_session = PREPROC.images_per_session';
    dat = preprocess(dat, 'outliers', 'plot');  % Spike detect and globals by slice
    %     drawnow
    
    subplot(5, 1, 5);
    dat = preprocess(dat, 'outliers_rmssd', 'plot');  % RMSSD Spike detect
    
    sz = get(0, 'screensize'); % Wani added two lines to make this visible (but it depends on the size of the monitor)
    set(gcf, 'Position', [sz(3)*.02 sz(4)*.05 sz(3) *.45 sz(4)*.85]);
    
    qcspikefilename = fullfile(PREPROC.qcdir, 'qc_spike_plot.png'); % Scott added some lines to actually save the spike images
    saveas(gcf,qcspikefilename);
    
    try_snapnow_for_publish;  % Data plot and mean heatmap
    diary off
    
    % Show results
    fh = findobj('Type', 'Figure', 'Tag', 'fmri data matrix');
    if ~isempty(fh)
        sz = get(0, 'screensize');
        % set(fh, 'Position', [50 sz(4)-50 sz(3)./2 sz(4)./2])
        set(fh, 'Position', [50 sz(4)*.2 sz(3)*.65 sz(4)*.75])
	set(fh, 'Visible', 'off');
        figure(fh);
        delete(subplot(2, 3, 6));
        figure(fh);
	set(fh, 'Visible', 'off');
        axh = subplot(2, 3, 5);
        set(axh, 'Position', [0.4108    0.1100    0.5134    0.3394]);
        axes(axh);
        title('Global values. circle size = spatial STD');
        scn_export_papersetup(500)
        
        fmrimatfilename = fullfile(PREPROC.qcdir, 'fmri_designmatrix.png'); %Scott added code to save images
        saveas(gcf,fmrimatfilename)
    end
    
    try_snapnow_for_publish;
    
    %clean up spikes and split by runs
    dat = split_and_clean_nuisance_covs(dat,1);
    
    % Save results
    PREPROC.Timeseries.global_raw = mean(dat.dat)';
    PREPROC.Nuisance.spikes = dat.covariates;
    
end

savePREPROC(PREPROC, preprocmatname);


%% -------------------------------------------------------------------------
% Movie of montages
% -------------------------------------------------------------------------

if writemovie
    
    % Movie of successive differences (sagittal slice)
    % ------------------------------------------------------
    
    movieoutfile = fullfile(PREPROC.qcdir, 'slice_movie.tiff');
    
    sagg_slice_movie(dat, movieoutfile);
    
    try_snapnow_for_publish
end


%% Interpolate bad points here

% *************************

%% -------------------------------------------------------------------------
% Create and threshold mean image; plot
% NOTE: SPIKE ID NO LONGER NEEDS/USES THIS
% -------------------------------------------------------------------------

if writemean
    print_header2('Mean functional image')
    
    % montage of mean image
    savefilename = fullfile(PREPROC.qcdir, 'mean_image_montage.png');
    
    mm = mean(dat);
    mm.fullpath = PREPROC.meanfilename;
    write(mm);
    
    canlab_preproc_show_montage(PREPROC.meanfilename, savefilename);
    title('Mean functional before Part 1');
    
    if dointeractive
        spm_image('init', PREPROC.meanfilename);
        saveimagenameo2 = fullfile(PREPROC.qcdir, 'mean_image_atmarker.png');
        input('Click on vitamin e capsule marker and press return: ','s');
        saveas(gcf, saveimagenameo2);
    end
    
end

try_snapnow_for_publish

%% -------------------------------------------------------------------------
% Preproc Part 1
% -------------------------------------------------------------------------

if dopart1
    diary(logfilename)
    print_header1('Running preproc_part1: Slice timing, realignment', ['Using' acq ' acquisition'])
    
    
    PREPROC = preproc_part1_2012(PREPROC, ...
        'slice timing', 1, 'motion correction', 1, 'local', 0, ...
        'verbose', 1, 'clean_up', 1, 'save plots', 0, spm_ver, ...
        'movie', 0, 'mean func', 1, 'run_dir_base', PREPROC.run_dir_base, 'acquisition', acq);
    
    fprintf('Saving PREPROC setup structure with UPDATED avol and ravol image names in: \n%s\n', preprocmatname)
    save(preprocmatname, '-append', 'PREPROC')
    
    diary off
    
    % Move Part1 files to Processed dir, and save names in PREPROC
    PREPROC = canlab_preproc_clean_up_and_move_files(PREPROC);
    
    savePREPROC(PREPROC, preprocmatname);
    
end

%% -------------------------------------------------------------------------
% Plot and save movement parameters
% Save nuisance covariate set with higher-order transformations
% Add spikes
% -------------------------------------------------------------------------
print_header2('Building and saving nuisance regressors in PREPROC.Nuisance');

if isfield(PREPROC, 'mvmt_param_files')
    
    [dummy, dummy, PREPROC.Nuisance.motion_cov_set] = canlab_preproc_motion_covariates(PREPROC.mvmt_param_files, PREPROC.images_per_session, 'savedir', PREPROC.qcdir);
    
    %generate temporary dat to adjust movement parameters
    mdat.covariates = PREPROC.Nuisance.motion_cov_set;
    mdat.images_per_session = PREPROC.images_per_session;
    %split parameters into runs without removing any of them
    mdat = split_and_clean_nuisance_covs(mdat,0);
    %save output values to PREPROC
    PREPROC.Nuisance.motion_cov_set = mdat.covariates;
    
    try_snapnow_for_publish
    
    savePREPROC(PREPROC, preprocmatname);
    
else
    print_header1('Movement parameter filenames not saved in PREPROC. Part1 not run?');
end

canlab_preproc_save_nuisance(PREPROC);


%% ------------------------------------------------------------------------
% Preproc Part 2
% -------------------------------------------------------------------------

if dopart2
    
    % Check some required files
    allerrors = canlab_preproc_check_PREPROC(PREPROC, 'anat_files', 'Anatomy files for Part 2', 1);
    allerrors = canlab_preproc_check_PREPROC(PREPROC, 'ra_func_files', 'Realigned images for Part 2', 1);
    
    PREPROC = preproc_part2_2012(PREPROC, varargin{:});
    
    try_snapnow_for_publish
    savePREPROC(PREPROC, preprocmatname);
    
    
end

%% ------------------------------------------------------------------------
% Clean up and finish
% -------------------------------------------------------------------------

if clean_up
    % Move wra, swra files (and any remaining ra_ files), delete a_files if
    % not done before, save filenames with new locations in PREPROC
    PREPROC = canlab_preproc_clean_up_and_move_files(PREPROC);
end

savePREPROC(PREPROC, preprocmatname);

disp('======================================================================================');
disp('=                                                                                    =');
disp('= Saved names of ra, wra, swra image files                                           =');
disp('= in PREPROC_SETUP.mat in PREPROC variable.                                          =');
disp('=                                                                                    =');
disp('= CANLAB_PREPROC done.                                                               =');
disp('=                                                                                    =');
disp('======================================================================================');

try_snapnow_for_publish

% html showing output files in directory and labels that describe what they are
% canlab_html_creator(basedir,'preproc')

% End of main function.

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% INLINE FUNCTIONS

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


    function setup_inputs()
        
        i=1;
        while i<=length(varargin)
            if(ischar(varargin{i}))
                switch(varargin{i})
                    case 'nodatachecks'
                        dodatachecks = 0; writemean = 0; writemovie = 0;
                    case {'nomovie'}
                        writemovie = 0;
                    case {'nomean'}
                        writemean = 0;
                    case {'nospikeid'}
                        dospikeid = 0;
                    case {'nopart1'}
                        dopart1 = 0;
                    case 'nopart2'
                        dopart2 = 0;
                    case {'disdaqs','disdaq'}
                        i=i+1;
                        disdaqs = varargin{i};
                    case 'nointeractive'
                        dointeractive = 0;
                        
                    case {'acquisition'}
                        i=i+1;
                        switch varargin{i}
                            case {'interleaved_BU','interleaved_TD','ascending','descending','interleavedMT'}
                                acq = varargin{i};
                            otherwise
                                error(['Unrecognized option to ''acquisition'': ' varargin{i}]);
                        end
                        
                        % Part 2 options
                    case {'verbose' 'display_commands' 'display commands'}, verbose = varargin{i+1};
                    case {'clean_up' 'clean up'}, clean_up = varargin{i+1};
                    case {'save_plots' 'save plots'}, save_plots = varargin{i+1};
                    case {'set_origins' 'set origins'}
                        setting_origins = varargin{i+1};
                    case {'coregistering_anat_to_func' 'coregistering anat to func' 'coreg anat to func'}
                        coregistering_anat_to_func = varargin{i+1};
                    case {'warp' 'warping'}
                        warping = varargin{i+1};
                    case {'generating_mean' 'generating mean' 'generate mean'}
                        generating_mean = varargin{i+1};
                    case {'smooth' 'smoothing'}
                        smoothing = varargin{i+1};
                    case {'check_norms' 'check norms'}
                        check_norms = varargin{i+1};
                    case 'noorigin'
                        setting_origins = 0;
                    case 'nocoreg'
                        coregistering_anat_to_func = 0;
                    case 'nowarp'
                        warping = 0;
                    case 'nosmooth', smoothing = 0;
                        
                    otherwise
                        error(['Unknown input string option: ' varargin{i}]);
                end
            end
            i=i+1;
        end
        
        %check if the base directory is valid
        if ~exist(basedir, 'dir'), error([basedir ' does not exist!']); end
        cd(basedir)
        
        if ~exist('TR', 'var') || isempty(TR)
            TR = input('Enter TR in sec: ');
        end
        
        % Define directories and file names for output
        % ----------------------------------------------------------------
        %create a directory for quality control images
        qcdir = fullfile(basedir, 'qc_images');
        if ~exist(qcdir, 'dir'), mkdir(qcdir); end
        
        logfilename = fullfile(qcdir, 'preprocessing_log.txt');
        meanfilename = fullfile(basedir, 'Functional', 'Raw', 'mean_func_before_preproc.img');
        
        preprocbase = fullfile(basedir, 'Functional', 'Preprocessed');
        if ~exist(preprocbase, 'dir'), mkdir(preprocbase); end
        meanfilename_realigned = fullfile(preprocbase, 'mean_ra_func.img');
        
        %meanmaskname = fullfile(qcdir, 'mean_mask_before_preproc.img');
        
        diary(logfilename)
        spm_ver = spm('ver');
        
        fprintf('Running for basedir = %s\n', basedir)
        fprintf('Saving in file: %s\n', logfilename)
        fprintf('Running the following version of SPM: %s\n',spm_ver);
        fprintf('\n');
        
        preprocdir = fullfile(basedir, 'Functional', 'Preprocessed');
        if ~exist(preprocdir, 'dir')
            mkdir(preprocdir);
        end
        preprocmatname = fullfile(preprocdir, 'PREPROC_SETUP.mat');
        
    end




end % main function





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


function savePREPROC(PREPROC, preprocmatname)
fprintf('Saving PREPROC setup structure with image names, etc., in file: \n%s\n', preprocmatname)

if exist(preprocmatname, 'file')
    try
        save(preprocmatname, '-append', 'PREPROC')
    catch
        % may not be able to find; perhaps doesn't exist, but file does?
        disp('Warning: Overwriting PREPROC because PREPROC variable not found in .mat file.')
        save(preprocmatname, 'PREPROC')
    end
else
    save(preprocmatname, 'PREPROC')
end
end


function anat = prep_structural_files(basedir, struct_wildcard)

struct_dir_path = fullfile(basedir, 'Structural', 'SPGR', struct_wildcard);
anat = filenames(struct_dir_path, 'char', 'absolute');

% There should only be one anatomical file here.
% If there are multiple structurals, average them before running this.
% You will have to deal with multiple structurals separately until we edit
% the code.
% For ex: for multiple structs on diff days, 1) coreg all structs to one
% another.  2) Make a mean (if same image type, e.g., T1) and use that mean
% in canlab_preproc.m.  3) After preproc part 2, the mean T1 will be
% coregistered, so you need to coreg (mutual info!!) all the other structurals to that one.

%the anat file
if (size(anat, 1) > 1)
    warning('THIS CODE ONLY SUPPORTS ONE SINGLE ANATOMICAL FILE NOW.\nTo handle multiple structurals, For ex: for multiple structs on diff days\n');
    warning('\t1) coreg all structs to one\n');
    warning('\t2) Make a mean (if same image type, e.g., T1) and use that mean in canlab_preproc.m.\n');
    warning('\t3) After preproc part 2, the mean T1 will be coregistered, so you need to coreg (mutual info!!) all the other structurals to that one.\n\n');
    warning('\tOnly the first anatomical file will be used, THE EXTRAS WILL BE IGNORED!\n')
elseif (size(anat, 1) < 1)
    fprintf('Looking for anatomical in:\n%s\n]', struct_dir_path);
    disp(['No anatomical files found using wildcard: ' struct_wildcard])
    error('Directory structure does not conform to CANlab standards or files are missing .')
end
end % setup anat


function try_snapnow_for_publish

fhi = findobj('Type', 'Figure', 'Tag', 'Interactive');
if ishandle(fhi), set(fhi, 'Visible', 'off'); end

try
    snapnow
catch
    warning('Error executing snapnow.  Old version of Matlab??');
end

end

function dat = split_and_clean_nuisance_covs(dat,doclean)

if nargin < 2
    %This eliminates doubled spikes.  The default is 'on'
    doclean = 1;
end
%create a copy of the covariates
spikes = dat.covariates;
%figure out how many images in each run
n_imgs = dat.images_per_session;
%get the number of runs
n_runs = numel(n_imgs);
%initialize the new covariate structure.  The first cells will be run
%specific and the final cell will contain covariates for use with
%concatenated data.
dat.covariates = cell(n_runs+1,1);
for i = 1:n_runs
    %choose the rows that contain data for this run and get those covs
    wh_rows = (1:n_imgs(i))+sum(n_imgs(1:i-1));
    covs = spikes(wh_rows,:);
    %if cleaning is requested,
    %check for copies in these spike regressors.
    %This may occur since spikes are defined both by Mahalanobis
    %distance and standard deviation
    covs = covs(:,var(covs)~=0);
    n_covs = size(covs,2);
    if doclean
        cov_pos = zeros(n_covs,1);
        for j = 1:n_covs
            cov_pos(j,1) = find(covs(:,j));
        end
        %only keep unique covs
        [dummy wh_covs] = unique(cov_pos);
    else
        wh_covs = 1:n_covs;
    end
    dat.covariates{i} = covs(:,wh_covs);
end

%repeat above loop, but use the entire concatenated set of covariates
n_spikes = size(spikes,2);
if doclean
    spike_pos = zeros(n_spikes,1);
    for i=1:n_spikes
        spike_pos(i) = find(spikes(:,i));
    end
    [dummy wh_spikes] = unique(spike_pos);
else
    wh_spikes = 1:n_spikes;
end
dat.covariates{n_runs+1} = spikes(:,wh_spikes);

%do a quick check to make sure that all the rows sum together properly
all_scans = size(dat.covariates{n_runs+1},1);
run_scans = 0;
for i = 1:n_runs
    run_scans = run_scans + size(dat.covariates{i},1);
end
if all_scans ~= run_scans
    error('The nuisance covariate processing has encountered an error: The number of scans per run do not match.')
end

end

