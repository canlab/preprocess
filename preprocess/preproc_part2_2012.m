% PREPROC = preproc_part2(PREPROC, ['set_origins', 0|1], ['warp', 0|1], ['smooth', 0|1], ['check norms', 0|1], ['SPM2' | 'SPM5' | 'SPM8'] ...
%   ['coreg anat to func', 0|1], ['generate mean', 0|1], ['verbose', 0|1], ['clean up', 0|1], ['save plots', 0|1])
%
% Coregistration, normalization and smoothing.
% Takes images from realigned (ravols) to smoothed, normalized functional images (swravols)
%
% PREPROC is a structure containing these fields:
%   .anat_files - cellstr of subjects' anatomical images
%   .files_to_warp - cell array of cellstrs of subjects' functional images or result images (also falls back to .ra_func_files)
%   .meanfilename_realigned - cellstr of subjects' mean prewarped functional images
%
% Optional PREPROC fields:
%   .inplane_files - cellstr of subjects' T1 images with the same number of slices as a functional image,
%        but much higher resolution within the slice - useful as an intermediate image, but not required
%   .template_image - canonical image to warp to - (default: avg152T1.nii);
%   .spm_ver: 'SPM2' or 'SPM5' - (default: SPM5)
%
% Other parameters
%   set_origins - set origins of images (default: 1)
%   warp - set normalization of functional images (default: 1)
%   smooth - set smoothing of functional images (default: 1)
%   check_norms - check the fitness of the normalization (default: 0)
%       - this is better run separately, after all preprocessing is done
%   coreg anat to func - coregister the anatomical image too the mean func (default: 1)
%       - don't turn this off unless you know for a fact they're already coregistered; if in doubt, leave on
%   generate mean - generate mean wra img (default: 1)
%       - meaningless if warping result images (e.g., con*nii, spmT*nii, beta*nii, etc)
%   SPM2 or SPM5 - flag to indicate which method to use (default: 'SPM5') - supercedes PREPROC.spm_ver
%
% E.g.:
% %load('subjs'); %{'rea1' 'rea2' 'rea3' ...}
% PREPROC.anat_files = filenames('rea*/structural/T1.nii');
%
% PREPROC.anat_files = filenames('rea*/structural/T1.nii');
% for i=1:length(subjs)
%   % or if the above anat line didn't work for you, maybe:
%   % PREPROC.anat_files{i} = filenames('rea*/structural/T1.nii');
%   PREPROC.files_to_warp{i} = filenames(sprintf('%s/run[0-9][0-9]/ravol*nii', subjs{i}));
% end
% preproc_part2(PREPROC);
%
% % with an inplane .img:
% PREPROC.inplane_files = filenames('rea*/structural/T1inplane.img');
% preproc_part2(PREPROC);
%
% % without smoothing:
% preproc_part2(PREPROC, 'smooth', 0);
% preproc_part2(PREPROC, 'nosmooth');
%
% without warping (smoothing not done either):
% preproc_part2(PREPROC, 'nowarp');
%
% Notes:
% 8/2012  Tor changed default version to SPM8
%         Changed default for making_mean_image to 0 (too slow)
%         Changed default for saving files to 1
%         Saves more plots
%         Takes single PREPROC structure as input, like preproc_part1,
%         which is automatically loaded and updated in canlab_preproc.m
%         Non-interactive mode
% ****need to show final files from each run

function PREPROC = preproc_part2_2012(PREPROC, varargin)

% -------------------------------------------------------------------------------------------------------------------------------------------
% Default values.
% These are updated during parse_inputs()
% -------------------------------------------------------------------------------------------------------------------------------------------

verbose = 0;
save_plots = 1;
clean_up = 1;

dointeractive = 1;
setting_origins = 1;
coregistering_anat_to_func = 1;
warping = 1;
generating_mean = 1;
smoothing = 1;
template_image = which('avg152T1.nii');
spm_ver = 'SPM8';  % Tor changed default version, 8/9/12

parse_inputs();

use_spm(spm_ver);
spm('defaults', 'fmri')

PREPROC = check_preproc(PREPROC); % Sets PREPROC.files_to_warp based on ra_func_files if needed.

% -------------------------------------------------------------------------------------------------------------------------------------------
% Do all the work here.
% -------------------------------------------------------------------------------------------------------------------------------------------

if(setting_origins)
    if ~dointeractive
        print_header2('Running in non-interactive mode: Skipping origin setting.');
    else
        set_origins(PREPROC, dointeractive);
    end
end

if coregistering_anat_to_func
    coreg_anat_to_func(PREPROC, dointeractive);
end

if warping
    PREPROC = segment_and_normalize(PREPROC);
end

if generating_mean  % Generate means for each run 
  generate_mean_wra_funcs(PREPROC.wra_func_files);
end

if smoothing
    PREPROC = smooth_funcs(PREPROC);
end

% End of main function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inline functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function parse_inputs()
        if(isfield(PREPROC, 'spm_ver'))
            spm_ver = PREPROC.spm_ver;
        end
        
        
        for i=1:length(varargin)
            if(ischar(varargin{i}))
                switch(varargin{i})
                    case {'verbose' 'display_commands' 'display commands'}
                        verbose = varargin{i+1};
                    case {'clean_up' 'clean up'}
                        clean_up = varargin{i+1};
                    case {'save_plots' 'save plots'}
                        save_plots = varargin{i+1};
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
                    case 'nointeractive'
                        dointeractive = 0;
                    case 'noorigin'
                        setting_origins = 0;
                    case 'nocoreg'
                        coregistering_anat_to_func = 0;
                    case 'nowarp'
                        warping = 0;
                        smoothing = 0;
                    case 'nomean'
                        generating_mean = 0;
                    case 'nosmooth'
                        smoothing = 0;
                        
                    case {'SPM99', 'SPM2', 'SPM5'}
                        print_header1('Warning! You are using an old version of SPM.', 'Use the original canlab_preproc for legacy pre-SPM8.');
                        error('quitting now.')
                end
            end
        end
        
    end % parse_inputs

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PREPROC = check_preproc(PREPROC)
errors = {};


if(~isfield(PREPROC, 'files_to_warp') && ~isfield(PREPROC, 'ra_func_files'))
    errors{end+1} = 'Must include either field "files_to_warp" or "ra_func_files"';
end
if(~isfield(PREPROC, 'meanfilename_realigned'))
    errors{end+1} = 'Missing field "meanfilename_realigned"';
end
if(~isfield(PREPROC, 'anat_files'))
    errors{end+1} = 'Missing field "anat_files"';
end

for f = {'inplane_files' 'ra_func_files' 'meanfilename_realigned' 'files_to_warp'}
    if ~isfield(PREPROC, f{1})
        PREPROC.(f{1}) = [];
    end
end

for f = {'ra_func_files' 'files_to_warp'}
    
    if ischar(PREPROC.(f{1}))
        PREPROC.(f{1}) = cellstr(PREPROC.(f{1}));
    end
end

for f = {'meanfilename_realigned' 'anat_files'}
    
    if iscell(PREPROC.(f{1}))
        PREPROC.(f{1}) = char(PREPROC.(f{1}));
    end
end

if isempty(PREPROC.files_to_warp)
    PREPROC.files_to_warp = PREPROC.ra_func_files;
end

% Expand to make sure OK for setting origins, etc.
for i = 1:length(PREPROC.files_to_warp)
    PREPROC.files_to_warp{i} = expand_4d_filenames(PREPROC.files_to_warp{i});
end

% THE FOLLOWING TEST IS PERFORMED IN THE List files and print SECTION
% WHEN RUNNING PART1. IT IS REPEATED HERE AS A REMINDER TO THE USER
if (size(PREPROC.anat_files, 1) > 1)
    warning('THIS CODE ONLY SUPPORTS ONE SINGLE ANATOMICAL FILE NOW.\nTo handle multiple structurals, For ex: for multiple structs on diff days\n');
    warning('\t1) coreg all structs to one\n');
    warning('\t2) Make a mean (if same image type, e.g., T1) and use that mean in canlab_preproc.m.\n');
    warning('\t3) After preproc part 2, the mean T1 will be coregistered, so you need to coreg (mutual info!!) all the other structurals to that one.\n\n');
    warning('\tOnly the first anatomical file will be used, THE EXTRAS WILL BE IGNORED!\n')
    PREPROC.anat_files = deblank(PREPROC.anat_files(1, :));
end

if ~exist(PREPROC.anat_files, 'file')
    print_header1('Anatomy file in PREPROC does not exist!', PREPROC.anat_files);
    errors{end+1} = 'Anatomy file in PREPROC does not exist!';
end

if length(PREPROC.files_to_warp) < 1
    errors{end+1} = 'No files_to_warp: should be functional images in cell array';
end

if any([ size(PREPROC.meanfilename_realigned, 1) size(PREPROC.anat_files, 1) ] > 1)
    errors{end+1} = 'meanfilename_realigned and anat_files cell arrays must have only one image each.';
end

if(~isempty(errors))
    error(['Errors encountered in PREPROC object:\n' implode(errors, '\n')], []);
end
end

% -------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
% set_origins
% -------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
%
% Make sure origins for structs and functionals are APPROXIMATELY
% on the anterior commissure.  Images will be coregistered after this, so
% this is just a starting point.
%
% Set the header of this:
% PREPROC.anat_files
%
% Also set the header of this (depending on whether we have an in-plane struct):
%
% PREPROC.inplane_files       % in-plane structural
% PREPROC.meanfilename_realigned  % mean realigned, slice-timing corr'd func across all runs
%
% Apply the second set of origins to these (functional) images,
% which are assumed to be in register at this point
%
% PREPROC.files_to_warp       % realigned, slice-timing corr'd funcs

function set_origins(PREPROC, dointeractive)
fprintf('Setting origins for: \n%s\n', PREPROC.anat_files);
display_key_subj_files(PREPROC, 'before');

try_snapnow_for_publish

if ~dointeractive
    return
end

% Set them
% --------------------------------------------------------------------------------------------------------------------------------------------------
disp('Setting origin of structural(s) and applying to functionals.');
set_hdr_current_coords(PREPROC.anat_files);

if(~isempty(PREPROC.inplane_files))
    set_hdr_current_coords(PREPROC.inplane_files, PREPROC.meanfilename_realigned, char(PREPROC.files_to_warp{:}));
else
    set_hdr_current_coords(PREPROC.meanfilename_realigned, char(PREPROC.files_to_warp{:}));
end

display_key_subj_files(PREPROC, 'after');
disp('Check images: Cursor is at [0 0 0]. ')
disp('Crosshairs for the images should be the same!!!');
disp('After setting, they should all be at the origin.');
input('Press return to continue');

% saveas(gcf, '
% try_snapnow_for_publish;

end



function display_key_subj_files(PREPROC, status)
imgs_to_display = char(PREPROC.files_to_warp{:});

if( ~isempty(PREPROC.anat_files) && ~isempty(PREPROC.meanfilename_realigned) && ~isempty(imgs_to_display))
    num_imgs = size(imgs_to_display, 1);
    wh_vols = [1 round(num_imgs/2) num_imgs];
    display_imgs = strvcat(PREPROC.anat_files, PREPROC.meanfilename_realigned, expand_4d_filenames(imgs_to_display(wh_vols,:), 1));
    
    spm_check_registration(display_imgs);
    if(strcmp(status, 'before')),
        fig_name = 'Before setting origin';
    else
        fig_name = 'After setting origin';
    end
    set(gcf, 'Name', fig_name);
    
    for j = 1:size(display_imgs, 1)
        [d, f] = fileparts(deblank(display_imgs(j,:)));
        spm_orthviews_name_axis(f, j);
    end
    spm_orthviews('Reposition', [0 0 0]);
else
    error('Something is wrong with image inputs. You need a mean func and individual func images.');
end

end % function

% -------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
% coreg_anat_to_func
% -------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
%
% Line the (single) structural image up with the (single) mean functional image.
% anat = PREPROC.anat_files;
% mean_func = PREPROC.meanfilename_realigned;
%
function coreg_anat_to_func(PREPROC, dointeractive)
%load('coreg_spm5_job');

print_header2('Coregisterng anatomical to mean functional');

% get images: reference is functional, source is anatomical
coregreference = deblank(PREPROC.meanfilename_realigned);    %PREPROC.ra_func_files{1}(1,:));
coregsource = PREPROC.anat_files;

coregdef = spm_get_defaults('coreg');
coreg_job = {};
%matlabbatch{1}.spm.spatial.coreg.estimate

coreg_job{1}.spatial{1}.coreg{1}.estimate.eoptions = coregdef.estimate;
coreg_job{1}.spatial{1}.coreg{1}.estimate.ref{1} = coregreference; % make robust file get %PREPROC.meanfilename_realigned;
coreg_job{1}.spatial{1}.coreg{1}.estimate.source{1} = coregsource;

spm_jobman('run', {coreg_job});

if dointeractive
    status = '';
else
    status = 'done';
end

fh = findobj('Type', 'Figure', 'Tag', 'Graphics');
if ishandle(fh), set(fh, 'Visible', 'on'), end

while(~strcmp(status, 'done'))
    %     anat = PREPROC.anat_files;
    %     mean_func = PREPROC.meanfilename_realigned;
    %
    %     coreg_job.spatial{1}.coreg{1}.estimate.source = {[anat ',1']};
    %     coreg_job.spatial{1}.coreg{1}.estimate.ref = {[mean_func ',1']};
    
    spm_check_registration(strvcat(coregsource, coregreference));
    spm_orthviews('Reposition', [0 0 0]);
    status = input('Type "done" if finished, or press return to adjust anatomical: ', 's');
    
    if(~strcmp(status, 'done'))
        fprintf('In dir: %s\n', pwd());
        fprintf('Use spm_image to shift anatomical image to better match the mean functional.\n');
        fprintf('Adjust the parameters (right, forward, up, pitch, roll, yaw) to alter the anatomical image.\n');
        fprintf('When you think you have it, click on "Reorient images..." and choose "%s" to apply the changes.\n', anat);
        fprintf('To check your changes, type "spm_check_registration(strvcat(anat, mean_func))"\n');
        fprintf('When you''re all done, type "return" and it''ll go through another pass using your changes.\n');
        
        spm_image('init', coregsource);
        spm_orthviews('Reposition', [0 0 0]);
        
        keyboard;
        
        spm_jobman('run', {coreg_job});
    end
end

if ~dointeractive
    spm_check_registration(strvcat(coregsource, coregreference, PREPROC.meanfilename_realigned));
    spm_orthviews('Reposition', [0 0 0]);
end

scn_export_papersetup(400);
saveas(gcf, fullfile(PREPROC.qcdir, 'coregistration.png'));

try_snapnow_for_publish
end% function


% -------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
% segment_and_normalize
% -------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
%
% Segment+warp the structural to the atlas template image
%
% Seg+warp this:   PREPROC.anat_files
%
% Apply the warps to these:
%                 PREPROC.anat_files   % single anatomical image
%                 PREPROC.meanfilename_realigned  % single func image
%                 PREPROC.files_to_warp   % list of ALL func images for the subject

function PREPROC = segment_and_normalize(PREPROC)

% Define the output names so we save them for later
% Saved in: wra_func_files and wanat_files fields
% -------------------------------------------------------------------------------------------------

PREPROC.wra_func_files = prepend_a_letter(PREPROC.files_to_warp, PREPROC.images_per_session, 'w');

wanat_files = prepend_a_letter({PREPROC.anat_files}, 1, 'w');
wanat_files = wanat_files{1};
PREPROC.wanat_files = wanat_files;

wmean_ra_func_files = prepend_a_letter({PREPROC.meanfilename_realigned}, 1, 'w');
wmean_ra_func_files = wmean_ra_func_files{1};
PREPROC.wmean_ra_func_files = wmean_ra_func_files;

[d f] = fileparts(PREPROC.anat_files);
PREPROC.norm_parameter_file = fullfile(d, [f '_seg_sn.mat']);

% Load (define) and customize the spm jobs for segment+norm and
% application to functional images (writenorm)
% -------------------------------------------------------------------------------------------------

load('segment_spm5_job');
load('writenorm_spm5_job');

% NOTE: THIS WILL NOT WORK for multiple subjects!!!
segment_job.spatial{1}.preproc.data = spm_image_list(PREPROC.anat_files);
segment_job.spatial{1}.preproc.opts.tpm = {which('grey.nii'); which('white.nii'); which('csf.nii')};

cnt = 1;

% define the name of the normalization parameter file for later
writenorm_job.spatial{1}.normalise{1}.write.subj(cnt).matname = {PREPROC.norm_parameter_file};
writenorm_job.spatial{1}.normalise{1}.write.subj(cnt).resample = spm_image_list([PREPROC.anat_files; PREPROC.meanfilename_realigned; PREPROC.files_to_warp]);
writenorm_job.spatial{1}.normalise{1}.write.subj(cnt).resample = cellstr(strvcat(writenorm_job.spatial{1}.normalise{1}.write.subj(cnt).resample{:}));

% Save the job
% -------------------------------------------------------------------------------------------------
savefile = fullfile(PREPROC.basedir, 'Functional', 'Preprocessed', 'canlab_preproc_norm_job.mat');
save(savefile, 'segment_job', 'writenorm_job');

% Now run it
% -------------------------------------------------------------------------------------------------

disp('Segmenting and normalizing all subjects');
spm_jobman('run', {segment_job writenorm_job});

% Now run it
% -------------------------------------------------------------------------------------------------
spm_check_registration(strvcat(wanat_files, wmean_ra_func_files, which('avg152T1.nii')));
spm_orthviews('Reposition', [0 0 0]);

scn_export_papersetup(400);
saveas(gcf, fullfile(PREPROC.qcdir, 'normalization.png'));

canlab_preproc_show_montage(PREPROC.wanat_files, fullfile(PREPROC.qcdir, 'warped_anatomical.png'));
canlab_preproc_show_montage(PREPROC.wra_func_files, fullfile(PREPROC.qcdir, 'wra_func_files.png'));

try_snapnow_for_publish

end % function

% -------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
% generate session-specific means
% -------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------

function mean_wra_func_name = generate_mean_wra_funcs(wra_func_files)

print_header2('Generating mean wra images for each run');

% Make sure these are in expanded format
wra_func_files = spm_image_list(wra_func_files, 0); % do not expand to individual cells

for i = 1:length(wra_func_files)
    
    [d f ext] = fileparts(wra_func_files{i}(1, :));
    d = fileparts(d);
    mean_wra_func_name{i, 1} = fullfile(d, sprintf('mean_wra_func_run_%02d%s', i, ext));
    
    dat = fmri_data(wra_func_files{i});
    m = mean(dat);
    clear dat
    m.fullpath = mean_wra_func_name{i, 1};
    write(m);
    
end

% display - but now done later
% canlab_preproc_montage_first_volumes(mean_wra_func_name);
% try_snapnow_for_publish;

end % function


% -------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
% smoothing
% -------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------

function PREPROC = smooth_funcs(PREPROC)
%load('smooth_spm5_job');

print_header2('Smoothing functional images');

%smooth_job.spatial{1}.smooth.data = spm_image_list(PREPROC.wra_func_files);

matlabbatch = {};
matlabbatch{1}.spm.spatial.smooth = spm_get_defaults('smooth');
matlabbatch{1}.spm.spatial.smooth.dtype = 0; % data type; 0 = same as before
matlabbatch{1}.spm.spatial.smooth.im = 0; % implicit mask; 0 = no
matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8]; % override whatever the defaults were with this
matlabbatch{1}.spm.spatial.smooth.data = spm_image_list(PREPROC.wra_func_files, 1); % individual cells for each volume

% Save the job
% -------------------------------------------------------------------------------------------------
savefile = fullfile(PREPROC.basedir, 'Functional', 'Preprocessed', 'canlab_preproc_smooth_job.mat');
save(savefile, 'matlabbatch');

% Run the job
% -------------------------------------------------------------------------------------------------
spm_jobman('run', matlabbatch);

% Display and save output
% -------------------------------------------------------------------------------------------------
PREPROC.swra_func_files = prepend_a_letter(PREPROC.wra_func_files, PREPROC.images_per_session, 's');
canlab_preproc_show_montage(PREPROC.swra_func_files, fullfile(PREPROC.qcdir, 'swra_func_files.png'));

end % function


% -------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
% utility functions
% -------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------

function image_list = spm_image_list(image_list, do_indiv_cellstrs)

% spm's batch handling appears to be changing in some versions of SPM8
% the optional flag do_indiv_cellstrs will make the files into a single,
% individual list of cells, one volume per cell.
% this may be appropriate for some uses in canlab_preproc_2012, but not
% all calls to this function.

if nargin < 2, do_indiv_cellstrs = 0; end

if ~iscell(image_list)
    image_list = cellstr(image_list); % should already be cells; just in case
end

for i = 1:length(image_list)
    image_list{i} = expand_4d_filenames(image_list{i});
end

doflatten = 0;
if do_indiv_cellstrs
   for i = 1:length(image_list)
    if size(image_list{i}, 1) > 1
        doflatten = 1;
    end
   end
end

if doflatten
    image_list = char(image_list{:});
    image_list = cellstr(image_list);
end


end % function



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
