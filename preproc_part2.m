%preproc_part2(PREPROC2, ['set_origins', 0|1], ['warp', 0|1], ['smooth', 0|1], ['check norms', 0|1], ['SPM2' | 'SPM5'], ['wh_subjs', n] ...
%   ['coreg anat to func', 0|1], ['generate mean', 0|1], ['verbose', 0|1], ['clean up', 0|1], ['save plots', 0|1])
%
% Coregistration, normalization and smoothing.
% Takes images from realigned (ravols) to smoothed, normalized functional images (swravols)
%
% PREPROC2 is a structure containing these fields:
%   .anat_files - cellstr of subjects' anatomical images
%   .files_to_warp - cell array of cellstrs of subjects' functional images or result images (also falls back to .ra_func_files)
%   .mean_ra_func_files - cellstr of subjects' mean prewarped functional images
%
% Optional PREPROC2 fields:
%   .inplane_files - cellstr of subjects' T1 images with the same number of slices as a functional image,
%        but much higher resolution within the slice - useful as an intermediate image, but not required
%   .templ - canonical image to warp to - (default: avg152T1.nii);
%   .spm_ver: 'SPM2' or 'SPM5' - (default: SPM5)
%
% Other parameters
%   wh_subjs - vector of numbers indicating which subjects to preprocess (default: all subjects)
%   set_origins - set origins of images (default: 1)
%   warp - set normalization of functional images (default: 1)
%   smooth - set smoothing of functional images (default: 1)
%   check_norms - check the fitness of the normalization (default: 0) 
%       - this is better run separately, after all preprocessing is done
%   coreg anat to func - coregister the anatomical image too the mean func (default: 1)
%       - don't turn this off unless you know for a fact they're already coregistered; if in doubt, leave on
%   generate mean - generate mean wra img (default: 1)
%       - meaningless if warping result images (e.g., con*nii, spmT*nii, beta*nii, etc)
%   SPM2 or SPM5 - flag to indicate which method to use (default: 'SPM5') - supercedes PREPROC2.spm_ver
%
% E.g.:
% %load('subjs'); %{'rea1' 'rea2' 'rea3' ...}
% PREPROC2.anat_files = filenames('rea*/structural/T1.nii');
%
% PREPROC2.anat_files = filenames('rea*/structural/T1.nii');
% for i=1:length(subjs)
%   % or if the above anat line didn't work for you, maybe:
%   % PREPROC2.anat_files{i} = filenames('rea*/structural/T1.nii');
%   PREPROC2.files_to_warp{i} = filenames(sprintf('%s/run[0-9][0-9]/ravol*nii', subjs{i}));
% end
% preproc_part2(PREPROC2);
%
% % with an inplane .img:
% PREPROC2.inplane_files = filenames('rea*/structural/T1inplane.img');
% preproc_part2(PREPROC2);
%
% % only the first 8 subjs
% preproc_part2(PREPROC2, 'wh_subjs', 1:8);
%
% % without smoothing:
% preproc_part2(PREPROC2, 'smooth', 0);

% Notes:
% 8/2012  Tor changed default version to SPM8
%         Changed default for making_mean_image to 0 (too slow)
%         Changed default for saving files to 1
%         Saves more plots

function preproc_part2(PREPROC2, varargin)

    % scn_setup refers to old FSL, etc. stuff
    % just call spm defaults
    % scn_setup();
    spm('defaults', 'fmri')

    verbose = 0;
    save_plots = 1;  
    clean_up = 1;
    wh_subjs = 1:size(PREPROC2.anat_files, 1);
    
    if length(wh_subjs) > 1 || wh_subjs ~= 1
        warning('STOP!!! ONLY ONE ANATOMICAL SHOULD BE ENTERED, AND WH_SUBJS SHOULD ALWAYS BE 1.');
        keyboard;
    end
    
    % ----------------------------------------------
    % Default values. 
    % These are updated during parse_inputs()
    % ----------------------------------------------
    
    dointeractive = 1;
    setting_origins = 1;
    coregistering_anat_to_func = 1;
    warping = 1;
    generating_mean = 0;
    smoothing = 1;
    check_norms = 0;
    templ = which('avg152T1.nii');
    spm_ver = 'SPM8';  % Tor changed default version, 8/9/12
    
    use_spm(spm_ver);
    
    parse_inputs();
    
    check_preproc(PREPROC2);

    % ----------------------------------------------
    % Do all the work here.
    % ----------------------------------------------
    
    if(setting_origins)
        set_origins(PREPROC2, wh_subjs, dointeractive);
    end
    
    if(coregistering_anat_to_func)
        coreg_anat_to_func(PREPROC2, wh_subjs, dointeractive);
    end
    
    [wra_func_files wanat_files] = segment_and_normalize(PREPROC2, warping, wh_subjs);

    wra_func_files = deblank(wra_func_files);

    mean_wra_funcs = generate_mean_wra_funcs(wra_func_files, generating_mean, wh_subjs);
    
    if(smoothing)
        smooth_funcs(wra_func_files, wh_subjs);
    end
    
    if(check_norms)
        scnlab_norm_check(templ, wanat_files, mean_wra_funcs, subjects);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inline functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function parse_inputs()
        if(isfield(PREPROC2, 'spm_ver'))
            spm_ver = PREPROC2.spm_ver;
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
                    case {'wh_subjs' 'wh subjs'}
                        wh_subjs = varargin{i+1};
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
                    case 'nointeractive'
                        dointeractive = 0;
                    case 'SPM2'
                        spm_ver = 'SPM2';
                    case 'SPM5'
                        spm_ver = 'SPM5';
                end
            end
        end
        if(~isfield(PREPROC2, 'inplane_files') || isempty(PREPROC2.inplane_files))
            PREPROC2.inplane_files = cell(size(PREPROC2.anat_files));
        end
        if(isfield(PREPROC2, 'anat_files') && ischar(PREPROC2.anat_files))
            PREPROC2.anat_files = cellstr(PREPROC2.anat_files);
        end
        if(isfield(PREPROC2, 'mean_ra_func_files') && ischar(PREPROC2.mean_ra_func_files))
            PREPROC2.mean_ra_func_files = cellstr(PREPROC2.mean_ra_func_files);
        end
        if(isfield(PREPROC2, 'files_to_warp'))
            for i=1:length(PREPROC2.files_to_warp)
                PREPROC2.files_to_warp{i} = cellstr(PREPROC2.files_to_warp{i});
            end
        elseif(isfield(PREPROC2, 'ra_func_files'))
            for i=1:length(PREPROC2.ra_func_files)
                PREPROC2.files_to_warp{i} = cellstr(PREPROC2.ra_func_files{i});
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_preproc(PREPROC2)
    errors = {};
    if(~isfield(PREPROC2, 'files_to_warp') && ~isfield(PREPROC2, 'ra_func_files'))
        errors{end+1} = 'Missing either field "files_to_warp" or "ra_func_files"';
    end
    if(~isfield(PREPROC2, 'mean_ra_func_files'))
        errors{end+1} = 'Missing field "mean_ra_func_files"';
    end
    if(~isfield(PREPROC2, 'anat_files'))
        errors{end+1} = 'Missing field "anat_files"';
    end

    if length(PREPROC2.files_to_warp) < 1
        errors{end+1} = 'No files_to_warp: should be functional images in cell array';
    end
    
    if any([ length(PREPROC2.mean_ra_func_files) length(PREPROC2.anat_files) ] > 1)
        errors{end+1} = 'mean_ra_func_files and anat_files cell arrays must have only one image each.';
    end
    
    if(~isempty(errors))
        error(['Errors encountered in PREPROC2 object:\n' implode(errors, '\n')], []);
    end
end

% -----------------------------------
% -----------------------------------
% set_origins
% -----------------------------------
% -----------------------------------
%
% Make sure origins for structs and functionals are APPROXIMATELY 
% on the anterior commissure.  Images will be coregistered after this, so
% this is just a starting point.
%
% Set the header of this:
% PREPROC2.anat_files
%
% Also set the header of this (depending on whether we have an in-plane struct):
%
% PREPROC2.inplane_files       % in-plane structural
% PREPROC2.mean_ra_func_files  % mean realigned, slice-timing corr'd func across all runs
%
% Apply the second set of origins to these (functional) images,
% which are assumed to be in register at this point
%
% PREPROC2.files_to_warp       % realigned, slice-timing corr'd funcs

function set_origins(PREPROC2, wh_subjs, dointeractive)
    for i = wh_subjs
        fprintf('Setting origins for subject %d\n', i);
        display_key_subj_files(PREPROC2, i, 'before');

        % Set them
        % -----------------------------------------------------
        disp('Setting origin of structural(s) and applying to functionals.');
        set_hdr_current_coords(PREPROC2.anat_files{i});

        if(~isempty(PREPROC2.inplane_files{i}))
            set_hdr_current_coords(PREPROC2.inplane_files{i}, PREPROC2.mean_ra_func_files{i}, char(PREPROC2.files_to_warp{i}{1:end}));
        else
            set_hdr_current_coords(PREPROC2.mean_ra_func_files{i}, char(PREPROC2.files_to_warp{i}{1:end}));
        end

        display_key_subj_files(PREPROC2, i, 'after');
        disp('Check images: Cursor is at [0 0 0]. ')
        disp('Crosshairs for the images should be the same!!!'); 
        disp('After setting, they should all be at the origin.');
        if dointeractive
            input('Press return to continue');
        end
        
        %disp(displayimgs);
        keyboard
        % SAVE IMAGE HERE
        
    end
end

function display_key_subj_files(PREPROC2, i, status)
    subj_ravols = char(PREPROC2.files_to_warp{i}{1:end});
    
    if( ~isempty(PREPROC2.anat_files{i}) && ~isempty(PREPROC2.mean_ra_func_files{i}) && ~isempty(subj_ravols))
        num_imgs = size(subj_ravols, 1);
        wh_vols = [1 round(num_imgs/2) num_imgs];
        display_imgs = strvcat(PREPROC2.anat_files{i}, PREPROC2.mean_ra_func_files{i}, expand_4d_filenames(subj_ravols(wh_vols,:), 1));

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
end

% -----------------------------------
% -----------------------------------
% coreg_anat_to_func
% -----------------------------------
% -----------------------------------
%
% Line the (single) structural image up with the (single) mean functional image.
% anat = PREPROC2.anat_files;
% mean_func = PREPROC2.mean_ra_func_files;
% 
function coreg_anat_to_func(PREPROC2, wh_subjs, dointeractive)
    load('coreg_spm5_job');

    for i = wh_subjs
        fprintf('Coregisterng anat to func for subject %d\n', i);
        
        if dointeractive
            status = '';
        else
            status = 'done';
        end
        
        while(~strcmp(status, 'done'))
            anat = PREPROC2.anat_files{i};
            mean_func = PREPROC2.mean_ra_func_files{i};
            coreg_job.spatial{1}.coreg{1}.estimate.source = {[anat ',1']};
            coreg_job.spatial{1}.coreg{1}.estimate.ref = {[mean_func ',1']};

            spm_jobman('run', {coreg_job});

            spm_check_registration(strvcat(anat, mean_func));
            spm_orthviews('Reposition', [0 0 0]);
            spm_print();
            status = input('Type "done" if finished, or press return to adjust anatomical: ', 's');
            if(~strcmp(status, 'done'))
                fprintf('In dir: %s\n', pwd());
                fprintf('Use spm_image to shift anatomical image to better match the mean functional.\n');
                fprintf('Adjust the parameters (right, forward, up, pitch, roll, yaw) to alter the anatomical image.\n');
                fprintf('When you think you have it, click on "Reorient images..." and choose "%s" to apply the changes.\n', anat);
                fprintf('To check your changes, type "spm_check_registration(strvcat(anat, mean_func))"\n');
                fprintf('When you''re all done, type "return" and it''ll go through another pass using your changes.\n');

                spm_image('init', anat);
                spm_orthviews('Reposition', [0 0 0]);

                keyboard;
            end
        end
    end
end


% -----------------------------------
% -----------------------------------
% segment_and_normalize
% -----------------------------------
% -----------------------------------
%
% Segment+warp the structural to the atlas template image
%
% Seg+warp this:   PREPROC2.anat_files
%
% Apply the warps to these:
%                 PREPROC2.anat_files   % single anatomical image
%                 PREPROC2.mean_ra_func_files  % single func image
%                 PREPROC2.files_to_warp   % list of ALL func images for the subject
function [wra_func_files wanat_files] = segment_and_normalize(PREPROC2, warping, wh_subjs)

    for i = wh_subjs
        
        % Define the output names so we save them for later
        % Saved in: wra_func_files and wanat_files fields
        % -----------------------------------
        
        for j = 1:length(PREPROC2.files_to_warp{i})
            [d f e] = fileparts(PREPROC2.files_to_warp{i}{j});
            wra_func_files{i}{j} = fullfile(d, ['w' f e]);
        end

        [d f e] = fileparts(PREPROC2.anat_files{i});
        wanat_files{i} = fullfile(d, ['w' f e]);
        
    end

    if(warping)
        
        % Load (define) and customize the spm jobs for segment+norm and
        % application to functional images (writenorm)
        % -----------------------------------
        
        load('segment_spm5_job');
        load('writenorm_spm5_job');

        % NOTE: THIS WILL NOT WORK for multiple subjects!!! wh_subjs should
        % ALWAYS be 1.
        segment_job.spatial{1}.preproc.data = spm_image_list(PREPROC2.anat_files(wh_subjs));
        segment_job.spatial{1}.preproc.opts.tpm = {which('grey.nii'); which('white.nii'); which('csf.nii')};

        cnt = 1;
        for i = wh_subjs
            
            [d f] = fileparts(PREPROC2.anat_files{i});
            
            % define the name of the normalization parameter file for later
            writenorm_job.spatial{1}.normalise{1}.write.subj(cnt).matname = {fullfile(d, [f '_seg_sn.mat'])};
            
            % define names of images to apply warps to later
            % PREPROC2.anat_files
            % PREPROC2.mean_ra_func_files
            % PREPROC2.files_to_warp
            
            writenorm_job.spatial{1}.normalise{1}.write.subj(cnt).resample = spm_image_list([PREPROC2.anat_files{i}; PREPROC2.mean_ra_func_files{i}; PREPROC2.files_to_warp{i}]);
            writenorm_job.spatial{1}.normalise{1}.write.subj(cnt).resample = cellstr(strvcat(writenorm_job.spatial{1}.normalise{1}.write.subj(cnt).resample{:}));

            [d f e] = fileparts(PREPROC2.mean_ra_func_files{i});
            wmean_ra_func_files{i} = fullfile(d, ['w' f e]);

            cnt = cnt + 1;
        end

        % Now run it
        % -----------------------------------
        
        disp('Segmenting and normalizing all subjects');
        spm_jobman('run', {segment_job writenorm_job});

        for i = wh_subjs
            spm_check_registration(strvcat(wanat_files{i}, wmean_ra_func_files{i}, which('avg152T1.nii')));
            spm_orthviews('Reposition', [0 0 0]);

            spm_print();
        end
    end
end

function mean_wra_funcs = generate_mean_wra_funcs(wra_func_files, generating_mean, wh_subjs)
    if(generating_mean), fprintf('Generating mean wra images.\n'); end
    for i = wh_subjs
        
        [d f ext] = fileparts(wra_func_files{i}{1});
        d = fileparts(d);
        mean_wra_funcs{i} = fullfile(d, ['mean_wra_func' ext]);
        
        if(generating_mean)
            %mean_image(char(wra_func_files{i}), mean_wra_funcs{i}, []);
            dat = fmri_data(char(wra_func_files{i}));
            m = mean(dat);
            m.fullpath = mean_wra_funcs{i};
            write(m);
            
        end
    end
end


function smooth_funcs(wra_func_files, wh_subjs)
    load('smooth_spm5_job');
    
    for i = wh_subjs
        fprintf('Smoothing for subject %d\n', i);
        
        smooth_job.spatial{1}.smooth.data = spm_image_list(wra_func_files{i});
        spm_jobman('run', {smooth_job});
    end
end


function image_list = spm_image_list(image_list)
    for i = 1:length(image_list)
        image_list{i} = expand_4d_filenames(image_list{i});
    end
end
