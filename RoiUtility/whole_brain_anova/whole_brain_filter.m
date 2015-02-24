% function Pw = whole_brain_filter(P, DX, TR, HP, nruns, numframes, [var args])
%
% single-subject timeseries extraction and model fitting
% saves filtered images f* and full-model betas (model fits)
%
%
% Inputs:
% P = list of image names, raw functionals or beta images (see "type"
%       below)
% DX = deconvolution (FIR) model matrix, see tor_make_deconv_mtx3.m
% TR = repetition time of your study (sampling rate)
% HP = high-pass filter cutoff in seconds to apply (SPM style)
% nruns = number of sessions (intercepts) OR num. of images in each sess,
%           e.g., [169 169 172]
% numframes = number of beta images per event type, e.g., [20 20 20] for FIR
%           model
% [dographics] = optional graphics
% [mask] = optional name of mask image to limit processing to brain, or []
% [smoothlen] = exponential smoothing filter length for timeseries and
%           betas; influence is 0 after smoothlen seconds; default = 6
% [startslice] = slice to start at, default = 1
% [type] = optional type of analysis to run.  Options:
%        'full'             Everything, from filtering to height/time/width, enter
%                           raw image names and model matrix DX
%
%        'filt', 'data'     Filter data only and save images in f*img
%                           Enter raw image names and empty model matrix []
%
%        'nofilt'           Everything, but no filtering. (Hint: Place filtered image in EXPT.FILES.im_files for whole_brain_fir)
%                           WARNING: Experimental
%
%        'htw'              Smooths beta series and computes height/time/width only
%                           Enter beta images from FIR model within condition,
%                            e.g., beta1_event1 beta2_event1 ... beta1_event2 b2_e2
%                            ...
%                           One possibility is to use images from SPM2 at this
%                           stage.
%         'glm'             Gen linear model settings; DX should be design
%                           matrix; DO NOT include intercept, this is
%                           automatically added.
%                           Pw = whole_brain_filter(fp, canon, TR, 120, 1, 10, 1, [], 0, 1, 'glm');
%
% Var. args: optional
% String, followed by input (param / value pair format)
% 'dographics'      followed by flag to make plots (1/0)
% 'mask'     followed by mask image name
% 'smoothlen'     followed by # images to 0 weight in exponential smoothing
% 'startslice'     followed by start slice #
% 'trialX'     followed by structure of trial-level design matrix stuff
%               (see single_trial_analysis)
%
%
% For a shell to run this function, see whole_brain_fir.m
%
% Tor Wager, 2/20/05, Modified 10/15/2006
% Modified 3/15/2007 single-trial update with full condition effects
% analysis
% More help in programmers' notes inside the function
%
% %
% Examples:
% -------------------------------------------------------------------------
% The analysis sequence for deconvolution (DX):
% Run whole_brain_fir, which mainly just calls whole_brain_filter
% whole_brain_filter gives you:
%       1) trimmed, filtered images, saved in f*img
%           these can be passed to SPM2, e.g., if desired
%       2) extraction of beta images for FIR model, saved in dx_beta*img
%           in individual subject folders; also smooths betas if desired
%       3) calculation of height, delay, and width for each FIR extracted
%           saved on images in individual subject folders.
%
% After running this, you'll need to collect image names for use in random
% effects analyses, and most likely create contrasts across conditions for
% differences in height, delay, and width.
% to do this, first go to main model directory, with subject subfolders
% to save a list of height/delay/width images:
%EXPT = get_htw_image_names(EXPT)
% to create contrasts across those images and save in EXPT.SNPM for rfx
% analysis:
%EXPT = make_htw_contrast_images(EXPT);
% Now you can run robust regression across the images in EXPT.SNPM.
% The robust reg uses image names in EXPT.SNPM.P, so you may have to save
% images from another field (e.g., EXPT.SNPM.heightP) in the .P field.
% you may also want to test individual contrasts against zero, in which
% case you can use EXPT.NLCON.height images, for example (and others in
% NLCON) instead.
% EXPT = robfit(EXPT);
%
% Example:
% -------------------------------------------------------------------------
% Single-trial estimation
% OLD NEEDS UPDATING as of 2007
% whole_brain_filter(images, glm_model, TR, HPfilter, nruns, not_used, 'glm'
% ..., 'trialX', trialX, 'bf', bf);
% images = list of img names
% glm_model = any glm matrix to fit: not single trial
%           this ends up being PSDXS in the function right now
%           Filtering matrix is added to this des. matrix before fitting
% trialX = trial-level design matrix; fit exactly as input; add any
%           filtering matrix, etc. to design matrix trialX first.
% bf = basis functions for trial-level model
% see trial_level_beta3 for explanations of trialX and bf


% Explanations of some more variables
% vb    = [optional] verbose output level: 0 none, 1 some, 2 lots
% mask  = [optional] mask 3-D volume to apply before extracting
%
% P     = image file names (str matrix)
% S     = filtering matrix (e.g., high-pass)
% DX    = full model to fit, unfiltered
% vb    = [optional] verbose output level: 0 none, 1 some, 2 lots
% mask  = [optional] mask 3-D volume to apply before extracting
% nsess = number of sessions (intercepts): assumes last nsess columns
%         of DX are run intercepts to be removed before trimming
%
% dims  = dimensions of image files in data
% cols  = no. of columns of DX, also no. of beta images to write
% SDX   = smoothed (filtered) full model to fit
% PSDX  = pseudoinverse of filtered full model
% PSDXS = PSDX * S, ready to multiply with data y for pinv(SX) * Sy
% betas = 4-D array of x, y, z, column
%
% ntrimmed = number of outliers trimmed from timeseries at each voxel
% Pw    = string matrix of output file names
%
% Pw = whole_brain_filter(d, c.model(1:166, :), 2, 120, 5, c.numframes, 1, 10);
function Pw = whole_brain_filter(P, DX, TR, HP, nruns, numframes, varargin)

    % -------------------------------------------------------------------
    % * configure input arguments
    % -------------------------------------------------------------------

    if isempty(P), warning('Empty image names!  No processing done.'), Pw = []; end

    % defaults
    % ---------------
    % for data filtering
    global dofilt
    global writefilt

    % for extracting and parameterizing hrfs
    global doysmooth
    global dosmooth
    global dowritebetas
    global dohtw

    % for plotting
    global f1
    global f2
    global f3

    % for knowing what to do
    global type
    global dographics

    % for single-trial model %%now:pass these in
    global trialmodels
    % %     global bf
    global trialpx
    % %     global ntrials

    global dims
    global cols

    type = 'full';    % predefined set of functions to do -- 'data' or 'full'
    mask = [];
    smoothlen = 0;
    startslice = 1;

    % glm
    X = [];

    % trial-level estimation
    trialmodels = [];    % single-trial design matrix
    %bf = [];        % basis functions

    % -------------------------------------------------------------------
    % * set inputs and flags
    % -------------------------------------------------------------------
    setup_inputs();

    mask;

    if dographics, drawnow; end

    % verbose output
    vb = 2; %if length(varargin) > 0, vb = varargin{1};  end


    setup_images();

    %if dographics, close, end

    % -------------------------------------------------------------------
    % * set filter
    % -------------------------------------------------------------------
    setup_filter();


    % -------------------------------------------------------------------
    % * adjust model DX based on filter, if not empty
    % -------------------------------------------------------------------
    % main outputs: PSDXS, F
    % if filter matrix DX is entered and type is not glm, adds filter to model
    % matrix, takes pinv, and returns in PSDXS
    % Otherwise, if glm, sets PSDXS to design matrix with smoothing covs

    setup_model();

    % finish setup of all components
    if vb > 0, fprintf(1, '\n\t\tFinished in %3.0f s\n', etime(clock, t1)), end

    % need to return some values so we can use them in another inline
    dotriallevel;
    dostandardglm;

    report_options();

    % -------------------------------------------------------------------
    % * for each slice...
    % -------------------------------------------------------------------

    for slicei = firstslice:dims(3)

        if vb > 0, t1 = clock;  fprintf(1, '\nSlice %3.0f \n------------>\n ', slicei), end

        [a, b] = process_slice(slicei, P, rfS, X, PSDXS, TR, V, Vall, F, numframes, mask(:, :, slicei), smoothlen);

        ntrimmed(:, :, slicei) = b;

        %dographics = 0;
        if vb > 0, fprintf(1, '\t%6.0f pts. trimmed, %6.0f s total', sum(sum(sum(~isnan(b)))), etime(clock, t1)), end
    end

    % -------------------------------------------------------------------
    % * write final images -- number of trimmed outliers image
    % -------------------------------------------------------------------
    if dofilt
        V.fname = 'ntrimmed.img';
        V.descrip = 'Number of outliers trimmed from timeseries at each voxel.';
        spm_write_vol(V, ntrimmed);

    end


    if vb > 0, fprintf(1, '%6.0f s Total, to dir: %s ', etime(clock, t1), pwd), end

    fprintf(1,'\nAll done.\n')

    Pw = [];




    % -------------------------------------------------------------------
    %
    %
    %
    %
    % * Inline functions
    %
    %
    %
    %
    % -------------------------------------------------------------------

    function setup_inputs

        % filtering
        % ---------------
        dofilt = 1;     % filter data
        writefilt = 1;  % write filtered data image set f*

        % model fitting
        % ---------------
        dohtw = 1;      % extract height, delay, and width
        dosmooth = 1;   % smooth beta series before hdw
        doysmooth = 1;  % smooth data before model fitting
        dowritebetas = 1;  % write beta images to disk as dx_*img

        % high-pass smoothing (defaults)
        smoothlen = 6;      % smoothing for beta series, in seconds (time to no effect)
        firstslice = 1;  % first slice -- option to start at something besides 1
        dographics = 1;

        mask = [];

        for i = 1:length(varargin)
            if ischar(varargin{i})
                switch varargin{i}
                    % reserved keywords
                    case 'full', type = 'full';
                    case 'htw', type = 'htw';
                    case 'glm', type = 'glm';
                    case 'nofilt', type = 'nofilt';

                        % functional commands
                    case 'smoothlen', smoothlen = varargin{i+1};
                    case {'startslice', 'firstslice'}, firstslice = varargin{i+1};


                    case 'dofilt', dofilt = varargin{i+1};
                    case 'writefilt', writefilt = varargin{i+1};
                    case 'dohtw', dohtw = varargin{i+1};
                    case 'dosmooth', dosmooth = varargin{i+1};
                    case 'dowritebetas', dowritebetas = varargin{i+1};
                    case {'dographics', 'doplot', 'graphics'}, dographics = varargin{i+1};
                    case 'mask', mask = varargin{i+1}; varargin{i+1} = [];

                        % GLM
                    case 'X', X = varargin{i+1};

                        % single-trial
                    case 'trialX', trialmodels = varargin{i+1}; %trialX = varargin{i+1}; trialpx = pinv(trialX);
                        %case 'bf', bf = varargin{i+1};
                        %case 'ntrials', ntrials = varargin{i+1};


                    otherwise, warning(['Unknown input string option:' varargin{i}]);
                end
            end
        end


        switch type
            case 'full', % do nothing
            case 'nofilt'
                dofilt = 0; writefilt = 0;
            case {'data', 'filt'} % data filtering only
                dohtw = 0; dosmooth = 0; doysmooth = 0; DX = [];
                fprintf(1, '\nWriting filtered data f*img only in:\n  %s\n', P(1, :))
            case 'htw', % enter beta images, do dosmooth and height/time/width only
                dofilt = 0; writefilt = 0;  doysmooth = 0;
                dohtw = 1; dosmooth = 1; dowritebetas = 0; DX = [];

            case 'glm'
                dohtw = 0; dosmooth = 0; doysmooth = 0; writefilt = 0; dofilt = 1;
                dowritebetas = 0;
                % but dofilt is turned off later

            otherwise
                error('Unknown instruction type')
        end

    end



    function setup_images

        fprintf(1,'\tWhole_brain_filter\n\t-----------------------------------------');

        if vb > 0, t1 = clock;  fprintf(1, '\n\t\tSubject setup (mapping images): %s', pwd), end

        % get rid of empty image names
        wh = find(all(P == ' ', 2));
        if ~isempty(wh), disp('WARNING: Some empty rows in image names.  Removing these.');
            P(find(wh), :) = [];
        end

        % defaults
        switch spm('Ver')
            case 'SPM2'
                % image dimensions, etc.
                V = spm_vol(P(1, :));
                Vall = spm_vol(char(expand_4d_filenames(P)));

            case 'SPM5'
                Vall = spm_vol(P);
                V = Vall(1);

            case 'SPM8'
                Vall = spm_vol(P);
                V = Vall(1);

            otherwise
                % unknown SPM
                disp('Unknown version of SPM!');

        end

        dims = V(1).dim(1:3);
        cols = size(DX, 2);

        % initialize beta matrix
        betas = NaN * zeros([dims cols]);

        if vb > 0, fprintf(1, '\n\t\tNew image dims: %3.0f %3.0f %3.0f %3.0f ', dims(1), dims(2), dims(3), cols), end

        % brain masking
        if ~exist('mask', 'var')
            mask = ones(dims);
            fprintf(1,'\nNo mask found.  Analyzing all voxels.\n');

        elseif ischar(mask)
            fprintf(1,'\nUsing mask: %s\n',mask);
            Vm = spm_vol(mask);
            %mask = spm_read_vols(Vm);

            % Read mask data in space of image data
            mask = scn_map_image(mask, P(1, :), 'write', 'mask.img');

        end
        if isempty(mask)
            mask = ones(dims);
            fprintf(1,'\nNo mask found.  Analyzing all voxels.\n');
        end

        ntrimmed = NaN * zeros(dims);

        fprintf(1,'\n\t\tFirst image: %s\n',Vall(1).fname);
        fprintf(1,'\t\tLast image: %s\n',Vall(end).fname);


    end




    function setup_filter

        npoints = size(P, 1);

        S = [];
        rfS = [];
        F = [];

        if dofilt

            fprintf(1, '\n\t\tSetting up filter for %3.0f images in %3.0f runs', npoints, length(nruns))

            % Set up rfS, residual-forming matrix based on nuisance params in
            % HP and intercepts
            % and F, LP filter matrix
            % -------------------------------------------------------------

            if isempty(HP)
                fprintf(1, '\n\t\tHP filter parameter empty: Windsorizing only.\n');

                % dofilt is windsorize only
                S = [];
                rfS = [];

            else
                % Get Residual-forming matrix hat such that hat * y  removes HP filter and
                % intercepts for each session; save in S

                [dummy, KL, KH] = use_spm_filter(TR, npoints, 'none', 'specify', HP);		% HP filter mtx
                clear dummy KL
                IX = intercept_model(nruns); % nruns has number of images in each run

                if size(IX, 1) ~= size(KH, 1), error('Filter and intercept len do not match.  Wrong number of runs?'), end

                S = [KH IX];                             % filtering matrix: nuisance covariates

                PS = pinv(S);
                hat = S * PS;                           % get residual forming matrix
                rfS = eye(length(hat)) - hat;             % DX-hat*DX to get residuals = (I-hat)DX
                clear hat PS IX

            end


        else
            % no filtering
            fprintf(1, '\n\t\tNo Windsorizing or pre-filtering. \n\t\t%3.0f images in %3.0f runs', npoints, nruns)


        end  % if dofilt

        if doysmooth
            [tmp, F] = smooth_timeseries(ones(npoints, 1), smoothlen); % F is LP filter matrix for data
        end

    end




    % -------------------------------------------------------------------
    % * adjust model DX based on filter, if not empty
    % -------------------------------------------------------------------
    % main outputs: PSDXS, F
    % if filter matrix DX is entered and type is not glm, adds filter to model
    % matrix, takes pinv, and returns in PSDXS
    % Otherwise, if glm, sets PSDXS to design matrix with smoothing covs

    function setup_model

        fprintf(1,'\n')

        if ~isempty(DX) && ~strcmp(type, 'glm')
            fprintf(1,'\t\tSetting up deconvolution design.\n');

            % For GLM, S is the actual predictor matrix, PSDXS
            % is the actual model matrix. this is to get valid single subj
            % stats
            fprintf(1, '\t\t...adjusting model matrix\n ')
            if size(S, 2) ~= size(DX, 1), error('Images and model size do not match!');  end

            SDX = rfS * DX;                           % model with intercept, HP terms removed, LP smooth
            PSDX = pinv(SDX);                       % pseudoinverse of filtered design matrix
            PSDXS = PSDX * rfS;                       % beta-forming matrix, post-multiplied by filter
            % last S is
            % applying filter S to
            % y
            % so that PSDXS * y = betas
            if dographics
                figure('Color', 'w'); subplot(1, 2, 1); set(gca, 'FontSize', 16)
                imagesc(S), title('Filtering matrix'), subplot(1, 2, 2); imagesc(SDX), title('Design matrix'), colormap hot, drawnow,
                scn_export_papersetup(500);
                saveas(gcf, 'scnlab_filters', 'png')
            end

        elseif strcmp(type, 'glm')
            SDX = []; PSDX = []; PSDXS = [];

            if ~isempty(DX)
                dostandardglm = 1;
                fprintf(1,'\t\tSetting up GLM design.\n');
                X = [DX S]; % X is the orig matrix, with intercepts and smoothing added based on your inputs
            else
                dostandardglm = 0;
                fprintf(1,'\t\tNo GLM design: assuming this is a single-trial design only.\n');
            end

            dotriallevel = 0;


            % Set up and check single-trial model
            % -------------------------------------------------------------

            if ~isempty(trialmodels) && isfield(trialmodels,'trialX')

                dotriallevel = 1;
                fprintf(1,'\t\tFound trial-level design matrix trialX.\n');
                if isfield(trialmodels,'bf')  %~isempty(bf)
                    fprintf(1,'\t\tFound trial-level basis set bf.\n');
                else
                    fprintf(1,'\t\tProblem: Basis functions bf not entered.\n');
                end

                if isfield(trialmodels,'ntrials') %~isempty(ntrials)
                    fprintf(1,'\t\tFound ntrials.\n');
                else
                    fprintf(1,'\t\tProblem: ntrials not entered.\n');
                end

                if ~exist('trialpx','var') || isempty(trialpx)
                    trialpx = pinv(trialmodels.trialX);
                end

                % check for ratings and mean-center if they exist
                if isfield(trialmodels,'ratings')

                    if size(trialmodels.ratings, 2) > 1
                        error('You cannot have more than one set of ratings (right now).')

                        fprintf(1,'\t\tFound ratings (behavioral data): Centering.\n');
                        trialmodels.ratings = scale(trialmodels.ratings, 1);
                    end
                end

            else
                fprintf(1,'\t\tTrial-level design structure trialmodels not entered: no single-trial model.\n');
            end

            clear DX

            if dographics
                figure('Color', 'w'); set(gca, 'FontSize', 16)
                imagesc(PSDXS);
                title('Design matrix');
                colormap gray;
                drawnow
                scn_export_papersetup(500);
                saveas(gcf, 'scnlab_design', 'png')
            end
        end

    end




    function report_options
        ynstr = {'No' 'Yes'};
        fprintf(1,'\tOptions:\n\t-----------------------------------------\n');
        fprintf(1,'\tAnalysis type: %s\n',type)

        switch type
            case 'glm'
                fprintf(1,'\t\tPre-fitting processing (dofilt): %s\n',ynstr{dofilt + 1})
                if dofilt, fprintf(1,'\t\tApply filter matrix: %s\n',ynstr{~isempty(rfS) + 1}); end
                fprintf(1,'\t\tStandard GLM option: %s\n',ynstr{dostandardglm + 1})
                fprintf(1,'\t\tTrial-level option: %s\n',ynstr{dotriallevel + 1})

            case {'full', 'htw', 'nofilt'}
                fprintf(1,'\t\tDeconvolution design\n')
                fprintf(1,'\t\tPre-fitting filtering (dofilt): %s\n',ynstr{dofilt + 1})
                if dofilt, fprintf(1,'\t\tApply filter matrix: %s\n',ynstr{~isempty(rfS) + 1}); end
                if dofilt, fprintf(1,'\t\tWrite filtered images: %s\n',ynstr{writefilt + 1}); end

                fprintf(1,'\t\tHeight, time to peak, width est.: %s\n',ynstr{dohtw + 1})
                fprintf(1,'\t\tSmooth betas before HTW est.: %s\n',ynstr{dosmooth + 1})
                fprintf(1,'\t\tWrite DX beta images: %s\n',ynstr{dowritebetas + 1})

        end

        fprintf(1,'\t-----------------------------------------\n');

    end



end % END MAIN FUNCTION






% -------------------------------------------------------------------
%
%
%
%
% * SUB-FUNCTIONS
%
%
%
%
% -------------------------------------------------------------------








function [fbetas, ntrimmed] = process_slice(slicei, P, rfS, X,  PSDXS, TR, V, Vall, F, numframes, varargin)
    %[a, b] = process_slice(slicei, P, S, PSDXS, TR, V, numframes, mask(:, :, slicei), smoothlen);

    % Notes: TR only needed for beta smoothing and HTW estimation

    % graphics
    global dographics
    global f1
    global f2
    global f3

    % options
    global type
    global dosmooth
    global doysmooth
    global dohtw
    global dofilt
    global writefilt
    global dowritebetas

    % trial level analysis
    global trialmodels
    % %     global bf
    global trialpx
    % %     global ntrials

    % general
    global dims
    global cols

    % initialize outputs
    fbetas = NaN * zeros([dims(1:2) cols]);
    ntrimmed = NaN * zeros([dims(1:2) 1]);

    % initialize inputs
    mask = []; if length(varargin) > 0, mask = varargin{1};  end
    smoothlen = 6; if length(varargin) > 1, smoothlen = varargin{2};  end

    % -------------------------------------------------------------------
    % * load the slice
    % -------------------------------------------------------------------
    if ~isempty(mask) && ~(nansum(mask(:)))  %sum(sum(mask))) > 0)
        % skip it
        fprintf(1, '...Empty slice...')
        return
    else
        sl = load_slice_data(mask,Vall,slicei);
    end

    % -------------------------------------------------------------------
    % * setup stuff we need in general about the slice data
    % -------------------------------------------------------------------

    nsize = size(sl);
    fsl = zeros(nsize);     % filtered output slices -- initialize
    nsize = nsize(1:2);

    whzero = all(sl == 0, 3);
    whnan = any(isnan(sl), 3);

    fprintf('Selecting voxels. %3.0f are zero, %3.0f are NaN, ', sum(whzero(:)), sum(whnan(:)));

    wvox = find(~whzero & ~whnan);
    fprintf(' %3.0f are valid.\n', length(wvox));

    [i, j] = ind2sub(size(sl(:, :, 1)), wvox);

    % -------------------------------------------------------------------
    % * filter data and store in fsl
    % Windsorize to 3 st. deviations
    % Apply rfS, if rfS is not empty
    % -------------------------------------------------------------------

    if dofilt

        % Windsorize to 3 st. deviations
        % Apply rfS, if rfS is not empty

        [fsl,ntrimmed] = windsorize_and_filter(sl,i,j,rfS);

    else
        fsl = sl;   % skip filtering
    end % if dofilt

    clear sl


    % -------------------------------------------------------------------
    % * write filtered images
    % -------------------------------------------------------------------

    emptyimg = zeros(dims);     % in case we need to create a new volume, used later as well

    if writefilt
        warning off
        fprintf(1, '\n\tWriting f (filtered) plane > ')
        et = clock;

        Vout = scn_write_plane(P, fsl, slicei, V(1));

        % %         for i = 1:size(P, 1)
        % %
        % %             V(i) = V(1);
        % %             [d, f, e]=fileparts(deblank(P(i, :)));
        % %             V(i).fname = fullfile(d, ['f' f e]);
        % %
        % %             if i == 1, fprintf(1, '%s .', V(i).fname), end
        % %
        % %             % create volume, if necessary
        % %             if ~(exist(V(i).fname) == 2), spm_write_vol(V(i), emptyimg); end
        % %
        % %             if ~(exist(V(i).fname) == 2)
        % %             Vout = struct('fname', V.fname, 'mat', V.mat, 'dim', V.dim);
        % %             if isfield(V, 'dt'), Vout.dt = V.dt; end      % SPM5 only
        % %             if isfield(V, 'n'), Vout.n = V.n;  end        % SPM5 only
        % %             Vout = spm_create_vol(Vout);
        % %             else
        % %                 Vout = spm_vol(V(i).fname);
        % %             end
        % %     %********
        % %             % Write data...in SPM5, by accessing file_array object in
        % %             % V.private directly.  The file_array object points to data in
        % %             % the actual file, so when values are assigned to the array
        % %             % object, they are written directly in the file.
        % %             % These values depend on the offset and slope values (spm
        % %             % scaling factors!) in V.private.dat as well, so care must be
        % %             % taken to assign data to a file_array object with the correct
        % %             % scaling factors.  This is why it is better to load the
        % %             % structure one wants to write to with spm_vol first, so that
        % %             % the name in V.private.dat.fname will be correct, and
        % %             % V.private.dat.scl_slope and ...inter will be correct as well.
        % %             Vout = spm_vol(V(i));  % loads correct .private info from V.fname
        % %
        % %             % in SPM5, this simply assigns data in Vout.private.dat
        % %             % in SPM2, it does something different, but should be
        % %             % compatible, since spm_vol was used above...
        % %             spm_write_plane(Vout, fsl(:, :, i), slicei);
        % %
        % %         end
        fprintf(1, 'Done in %3.0f s\n', etime(clock, et))
        warning on
    end


    % -------------------------------------------------------------------
    % * Model fitting: GLM/trial GLM or FIR
    % -------------------------------------------------------------------

    % Model-fitting stuff (and TS smoothing if req.)
    switch type
        % can be 'glm' 'nofilt' 'htw' 'full'

        case 'glm'
            % GLM fitting and/or trial-level model
            % do not do deconvolution

            if ~isempty(X)
                % Fit GLM model with AR(nphis)
                % PSDXS is original design matrix for GLM (not pinv of design)
                %function [betas,tvals,stes,phis,dfs,sigmas] = glm_slice(fsl,i,j,PSDXS,cols,nphis,slicei,V,dims)
                % no need to end values, actually, because this writes them to disk
                nphis = 2;  % AR process order
                glm_slice(fsl,i,j,X,cols,nphis,slicei,V,dims);

            end

            if ~isempty(trialmodels) % modified so trialX is structure, 3/15/07 && ~isempty(bf)
                % Fit trial-level GLM model and then AR(1) at 2nd level
                trial_glm_slice(fsl,i,j,trialmodels,cols,slicei,V,dims);
            end



        case {'nofilt' 'htw' 'full'}

            if ~isempty(PDSXS)
                % to fit FIR model
                betas = fir_slice(fsl,i,j,PSDXS,F);

            elseif (dosmooth || dohtw || dowritebetas)
                % If we have no model, but we want to do stuff with betas
                % Assume that we've entered Beta images in timeseries as slice
                betas = fsl;

            else   %  & isempty(PSDXS)
                % -------------------------------------------------------------------
                % * skip timeseries smooth and model fit if no model matrix is entered
                % -------------------------------------------------------------------
                % We're done, nothing left to do.
                betas = NaN * zeros([dims(1:2) cols]);
                return
            end

            if dosmooth
                fbetas = smooth_beta_series(betas,numframes,smoothlen,TR,dims);
            else
                fbetas = betas;
            end

            % -------------------------------------------------------------------
            % * get height, width, delay parameters
            % -------------------------------------------------------------------
            if dohtw
                get_htw(fbetas,TR,numframes,st,en,emptyimg);

            end

        otherwise
            error('Unknown type keyword.')
    end


end % END SLICE FUNCTION






% -------------------------------------------------------------------
% * load the slice
% -------------------------------------------------------------------
function sl = load_slice_data(mask,Vall,slicei)

    global dographics
    global f1

    fprintf(1, '\tLoading data >')
    et = clock;

    sl = timeseries_extract_slice(Vall, slicei);

    % load it, check values for middle slice only
    %if slicei == round(dims(3)./2),
    %    [ts, sl] = timeseries4(round(dims./2), P);
    %else,
    %    [ts, sl] = timeseries4(round(dims./2), P, [], 1);
    %end

    mask = abs(mask) > 0;  % returns zero for 0 or nan values

    fprintf(1, ' loaded in %3.2f s.', etime(clock, et))
    %if ~isempty(mask), sl(:, :, 1) = sl(:, :, 1) .* mask;  end

    if ~isempty(mask), sl = sl .* repmat(mask, [1 1 size(sl, 3)]); end


    if dographics
        if ishandle(f1), figure(f1); else f1 = figure; set(gcf, 'Color', 'w');  end
        subplot 221; imagesc(mask); colormap(gray); title('Mask'), colorbar('horiz'), drawnow;
        subplot 222; imagesc(sl(:, :, 2)); colormap(gray); title('Slice timepoint 2'), colorbar('horiz'), drawnow;
        subplot 223; cla
        subplot 224; cla
    end

end





% -------------------------------------------------------------------
% * filter data and store in fsl
% Windsorize to 3 st. deviations
% Apply rfS, if rfS is not empty
% -------------------------------------------------------------------
function [fsl,ntrimmed] = windsorize_and_filter(sl,i,j,rfS)

    global dographics
    global f1
    global f2
    global dims

    ntrimmed = NaN * zeros([dims(1:2) 1]);

    fsl = NaN * zeros(size(sl));

    if isempty(rfS), dofilt = 0; filtstr = ''; else, dofilt = 1; filtstr = 'and filtering'; end

    fprintf(1, '\n\tPreprocessing (Windsorizing to 3SD %s) %3.0f voxels > ', filtstr,length(i))

    et = clock;
    for k = 1:length(i)

        % first image usually off -- replace w/mean of first 10 images
        % (now this is handled by separate betas in GLM for first
        % images; optional)
        y = squeeze(sl(i(k), j(k), :));
        %y(1) = mean(y(min(length(y), 10)));

        % filter
        if dofilt
            y = rfS * y;
        end

        % trim
        [y, ntrimmed(i(k), j(k))] = trimts(y, 3, []);
        fsl(i(k), j(k), :) = y;   % fsl is filtered slices

        if k == 1000, fprintf(1, '%3.0f s per 1000 vox.', etime(clock, et)), end
    end

    if dographics
        if ishandle(f2), figure(f2);subplot(2, 1, 1); else f2 = tor_fig(2, 1);  end, cla
        plot(scale(squeeze(sl(round(dims(1)./2), round(dims(2)./2), :))), 'k', 'LineWidth', 2)
        plot(scale(squeeze(fsl(round(dims(1)./2), round(dims(2)./2), :))), 'r', 'LineWidth', 2),
        legend({'Raw' 'Filtered'}); title('Standardized data from center voxel')
        subplot(2, 1, 2); cla;
    end

    clear sl

    if dographics
        figure(f1); subplot 223; imagesc(ntrimmed); colormap(gray); title('ntrimmed'), colorbar('horiz'), drawnow;
    end

end





% -------------------------------------------------------------------
% * write images
% -------------------------------------------------------------------
function Pw = write_beta(bi, V, betas)
    % Volumetric version
    if bi < 10, myz = '000';  elseif bi < 100, myz = '00';  else myz = '000'; end
    V.fname = ['dx_beta_' myz num2str(bi) '.img'];
    V.descrip = ['Beta image for column ' num2str(bi)];

    spm_write_vol(V, betas);
    Pw = which(V.fname);
end



function Pw = write_beta_slice(slicei, V, betas, emptyimg, varargin)
    % Pw = write_beta_slice(sliceindex, V, data, empty, prefix, suppressnum)
    % Slice-a-metric version
    warning off % due to NaN to int16 zero conversions

    % alter data type if needed
    switch(spm('Ver'))
        case 'SPM2'
            V.dim(4) = 16; % set to float to preserve decimals

        case 'SPM5'
            V.dt(1) = 16;
    end


    prefix = 'dx_beta_';
    suppress = 0;
    if length(varargin) > 0, prefix = varargin{1}; end
    if length(varargin) > 1, suppress = varargin{2}; end

    nbetas = size(betas, 3);

    P = {};
    for voli = 1:nbetas
        if suppress
            toappend = '';
        else
            if voli < 10, myz = '000';  elseif voli < 100, myz = '00';  else myz = '000'; end
            toappend = [myz num2str(voli)];
        end

        P{voli} = [prefix toappend '.img'];
    end


    Vout = scn_write_plane(P, betas, slicei, V(1));


    % % %         V.fname = [prefix toappend '.img'];
    % % %         V.descrip = ['Image for column ' num2str(voli)];
    % % %
    % % %         % create volume, if necessary
    % % %         if ~(exist(V.fname) == 2), spm_write_vol(V, emptyimg); end
    % % %
    % % %         % SPM5 uses info from private; can't just replace filename
    % % %         % alter data type if needed
    % % %         switch(spm('Ver'))
    % % %             case 'SPM2'
    % % %                 Vout = V;
    % % %
    % % %             case 'SPM5'
    % % %                 Vout = spm_vol(V.fname);
    % % %
    % % %             case 'SPM8'
    % % %                 error('Fix or verify image writing for SPM8. Not implemented yet.');
    % % %
    % % %             otherwise
    % % %                 error('Unknown SPM version.');
    % % %         end
    % % %
    % % %         spm_write_plane(Vout, betas(:, :, voli), slicei);
    % % %     end

    Pw = which(Vout(1).fname);

    warning on
end






function [betas,tvals,stes,phis,dfs,sigmas] = glm_slice(fsl,i,j,PSDXS,cols,nphis,slicei,V,dims)
    % no need to end values, actually, because this writes them to disk

    % * GLM estimation!
    % -------------------------------------------------------------------
    fprintf(1, '\tGLM estimation for %3.0f voxels > ', length(i))

    % setup inputs
    % -------------------------------------------------------------

    nextra = size(PSDXS, 2) - cols;
    c = blkdiag(eye(cols), zeros(nextra));
    % indicator for each regressor of interest; could be changed to do contrasts!
    % c should be column vectors of contrasts

    ncons = cols;
    px = pinv(PSDXS);

    % setup outputs
    % -------------------------------------------------------------

    % one per contrast of interest
    betas = NaN * zeros([dims(1:2) ncons]);
    tvals = NaN * zeros([dims(1:2) ncons]);
    stes = NaN * zeros([dims(1:2) ncons]);

    % other (ar# or 1 param)
    phis = NaN * zeros([dims(1:2) nphis]);
    dfs = NaN * zeros([dims(1:2) 1]);
    sigmas = NaN * zeros([dims(1:2) 1]);

    % -------------------------------------------------------------------

    et = clock;
    for k = 1:length(i)

        voxdat = squeeze(fsl(i(k), j(k), :));

        if all(voxdat == mean(voxdat))
            % nothing to do
        else

            %[t, df, beta, Phi, sigma, stebeta] = fit_gls(voxdat, PSDXS, c, nphis, px);
            [beta, t, pvals, convals, con_t, con_pvals, sigma, Phi, df, stebeta, conste, F] = fit_gls(voxdat, PSDXS, c, nphis, px);

            % one per contrast of interest
            betas(i(k), j(k), :) = beta(1:ncons);
            tvals(i(k), j(k), :) = t(1:ncons);
            stes(i(k), j(k), :) = stebeta(1:ncons);

            % other (ar# or 1 param)
            phis(i(k), j(k), :) = Phi;
            dfs(i(k), j(k), :) = df;
            sigmas(i(k), j(k), :) = sigma;
            if k == 1000, fprintf(1, '%3.0f s per 1000 vox.', etime(clock, et)), end

        end
        if dographics && (mod(k, 1000) == 0)
            if ishandle(f1), figure(f1); else f1 = figure; set(gcf, 'Color', 'w');  end
            subplot 223; imagesc(squeeze(tvals(:, :, 1))); colormap(hot); title('t-vals(1)'), colorbar('horiz'), drawnow;

        end


    end  % loop voxels

    if dographics
        scn_export_papersetup(500);
        saveas(f1, ['slice_' num2str(slicei)], 'png');
    end


    % * write betas
    % -------------------------------------------------------------
    fprintf(1, '\n\tWriting GLM images for this plane > ')
    et = clock;
    emptyimg = zeros(dims);     % in case we need to create a new volume

    write_beta_slice(slicei, V(1), betas, emptyimg, 'beta_');
    write_beta_slice(slicei, V(1), tvals, emptyimg, 't_');
    write_beta_slice(slicei, V(1), stes, emptyimg, 'ste_');
    write_beta_slice(slicei, V(1), phis, emptyimg, 'phi_');
    write_beta_slice(slicei, V(1), dfs, emptyimg, 'df', 1);
    write_beta_slice(slicei, V(1), sigmas, emptyimg, 'sigma', 1);

end






function trial_glm_slice(fsl,i,j,trialmodels,cols,slicei,V,dims)

    % setup inputs
    % -------------------------------------------------------------

    % ***to-do: update to be nice if empty/missing.
    trialX = trialmodels.trialX;
    bf =   trialmodels.bf ;
    ntrials =  trialmodels.ntrials;
    trialpx = trialmodels.trialpx;

    % already entered/but not used!  need this!!: trialpx =  trialmodels.trialpx;
    trial2ndX = trialmodels.subjlevelX;
    trial2ndpx = trialmodels.subjlevelpx;

    nconds = size(trial2ndX,2);

    c = trialmodels.conditioncontrasts;
    cnames = trialmodels.conditionnames;

    if isfield(trialmodels, 'contrastnames')
        contrastnames = trialmodels.contrastnames;
    end

    if isempty(trial2ndX)
        trial2ndX = ones(ntrials, 1);
        trial2ndpx = pinv(trial2ndX);
    end

    if isfield(trialmodels,'ratings')
        ratings = trialmodels.ratings;
        ratings(:,end+1) = 1; % intercept
        ratingspx = pinv(ratings);
    end

    % set up number of trials (should enter this, actually!)
    if isempty(ntrials)
        ntrials = (size(trialX, 2) - mod(size(trialX, 2), size(bf, 2))) ./ size(bf, 2);
        if round(ntrials) ~= ntrials, error('Check trial-level model: doesn''t match bf.'), end
    end

    fprintf(1, '\n\tTrial-level estimation: %3.0f trials, %3.0f basis functions) > ', ntrials, size(bf, 2))

    % setup outputs
    % -------------------------------------------------------------

    trialheight = NaN * zeros([dims(1:2) ntrials]);
    trialdelay = NaN * zeros([dims(1:2) ntrials]);
    trialwidth = NaN * zeros([dims(1:2) ntrials]);
    trialauc = NaN * zeros([dims(1:2) ntrials]);

    trialphi = NaN * zeros([dims(1:2) 1]);  % AR(1) model
    sigma_level1 = NaN * zeros([dims(1:2) 1]);  % error

    conditiont = NaN * zeros([dims(1:2) nconds]);
    conditionmean = NaN * zeros([dims(1:2) nconds]);
    conditionp = NaN * zeros([dims(1:2) nconds]);
    conditionste = NaN * zeros([dims(1:2) nconds]);


    if ~isempty(c)
        if size(c,1) ~= size(trial2ndX,2)
            disp('Contrasts are wrong size.  Probably transposed...transposing now.');
            c = c';
        end

        ncontrasts = size(c,2); % contrasts are column vectors
        contrastt =  NaN * zeros([dims(1:2) ncontrasts]);
        contrastste = NaN * zeros([dims(1:2) ncontrasts]);
        contrastp = NaN * zeros([dims(1:2) ncontrasts]);
        contrastest = NaN * zeros([dims(1:2) ncontrasts]);
    end

    if exist('ratings','var')
        ratingsest =  NaN * zeros([dims(1:2) 1]);
        ratingst =  NaN * zeros([dims(1:2) 1]);
        ratingsste = NaN * zeros([dims(1:2) 1]);
        ratingsp = NaN * zeros([dims(1:2) 1]);

        intcptest =  NaN * zeros([dims(1:2) 1]);
        intcptt =  NaN * zeros([dims(1:2) 1]);
        intcptste = NaN * zeros([dims(1:2) 1]);
        intcptp = NaN * zeros([dims(1:2) 1]);
    end

    et = clock;
    for k = 1:length(i)

        voxdat = squeeze(fsl(i(k), j(k), :));

        if all(voxdat == mean(voxdat))
            % nothing to do
        else

            % trial-level betas and heights for one voxel
            % -------------------------------------------------------------
            [b1tmp, f, X2, px2, bf, h, t, w, auc] = trial_level_beta3('X', trialX, 'pinvx', trialpx, 'bf', bf, ...
                'output', 'betas', 'plot', 0, 'data', voxdat, 'ntrials', ntrials);

            % If two basis sets, do it twice, and then make sure the h, t,
            % w associated with the correct basis set for each trial is
            % returned.
            if length(trialmodels.basistype) > 1 && isfield(trialmodels, 'bf2')
                if ~isfield(trialmodels, 'wh_bf2')
                    error('Enter trialmodels.wh_bf2, which is the list of index values belonging to 2nd basis set.');
                end

                [dummy, dummy, dummy, dummy, dummy, h2, t2, w2, auc2] = trial_level_beta3('X', trialX, 'pinvx', trialpx, 'bf', trialmodels.bf2, ...
                    'output', 'betas', 'plot', 0, 'data', voxdat, 'ntrials', ntrials);

                % replace est. for events modeled with bf2
                h(trialmodels.wh_bf2) = h2(trialmodels.wh_bf2);
                t(trialmodels.wh_bf2) = t2(trialmodels.wh_bf2);
                w(trialmodels.wh_bf2) = w2(trialmodels.wh_bf2);
                auc(trialmodels.wh_bf2) = auc2(trialmodels.wh_bf2);

            end

            resid = voxdat - trialX * b1tmp;
            ny = length(resid);
            nparams = size(trialX, 2);
            sigma = sqrt((1 / (ny - nparams) * sum(resid.^2)));  % Estimate of Sigma

            sigma_level1(i(k), j(k), 1) = sigma;

            trialheight(i(k), j(k), :) = h;
            trialdelay(i(k), j(k), :) = t;
            trialwidth(i(k), j(k), :) = w;
            trialauc(i(k), j(k), :) = auc;

            % 2nd-level estimation of t-value with AR(1) model
            % -------------------------------------------------------------
            %[t, df, beta, Phi, sigma, stebeta] = fit_gls(h, trial2ndX, [], 1, trial2ndpx);
            [beta, t, pvals, convals, con_t, con_pvals, sigma, Phi, df, stebeta, conste, F] = fit_gls(h, trial2ndX, [], 1, trial2ndpx);

            conditiont(i(k), j(k), :) = t;
            conditionste(i(k), j(k), :) = stebeta;
            conditionmean(i(k), j(k), :) = beta;
            trialphi(i(k), j(k), :) = Phi;

            conditionp(i(k), j(k), :) = 2 .* (1 - tcdf(abs(t), df));

            if ~isempty(c)
                stecons = diag(c' * diag(stebeta).^2 * c).^.5;
                con_est = c' * beta;
                tcons = con_est ./ stecons;

                contrastest(i(k), j(k), :) = con_est;
                contrastt(i(k), j(k), :) = tcons;
                contrastste(i(k), j(k), :) = stecons;
                contrastp(i(k), j(k), :) = 2 .* (1 - tcdf(abs(tcons), df));
            end


            % 2nd-level estimation of ratings model with AR(1) model
            % use first param, ignore intercept
            % -------------------------------------------------------------
            if exist('ratings','var')
                %[t, df, beta, Phi, sigma, stebeta] = fit_gls(h, ratings, [], 1, ratingspx);
                [beta, t, pvals, convals, con_t, con_pvals, sigma, Phi, df, stebeta, conste, F] = fit_gls(h, ratings, [], 1, ratingspx);

                % 1st col. is rating effect
                ratingsest(i(k), j(k), 1) = beta(1);
                ratingst(i(k), j(k), 1) = t(1);
                ratingsste(i(k), j(k), 1) = stebeta(1);
                ratingsp(i(k), j(k), 1) = 2 .* (1 - tcdf(abs(t(1)), df));

                % 2nd column is intercept
                intcptest(i(k), j(k), 1) = beta(2);
                intcptt(i(k), j(k), 1) = t(2);
                intcptste(i(k), j(k), 1) = stebeta(2);
                intcptp(i(k), j(k), 1) = 2 .* (1 - tcdf(abs(t(2)), df));

            end
        end

        if k == 1000, fprintf(1, '%3.0f s per 1000 vox.', etime(clock, et)), end


    end % loop through voxels

    % * write images
    % -------------------------------------------------------------
    fprintf(1, '\n\tWriting trial GLM images for this plane > ')
    et = clock;
    emptyimg = zeros(dims);     % in case we need to create a new volume


    write_beta_slice(slicei, V(1), sigma_level1, emptyimg, 'sigma_level1_');

    write_beta_slice(slicei, V(1), trialheight, emptyimg, 'trial_height_');
    write_beta_slice(slicei, V(1), trialdelay, emptyimg, 'trial_delay_');
    write_beta_slice(slicei, V(1), trialwidth, emptyimg, 'trial_width_');
    write_beta_slice(slicei, V(1), trialauc, emptyimg, 'trial_AUC_');

    write_beta_slice(slicei, V(1), conditiont, emptyimg, 'condition_t_');
    write_beta_slice(slicei, V(1), conditionste, emptyimg, 'condition_ste_');
    write_beta_slice(slicei, V(1), conditionmean, emptyimg, 'condition_height_');
    write_beta_slice(slicei, V(1), trialphi, emptyimg, 'trial_ar1', 1);
    write_beta_slice(slicei, V(1), conditionp, emptyimg, 'condition_p');

    if ~isempty(c)
        if ~exist('contrastnames','var')
            write_beta_slice(slicei, V(1), contrastest, emptyimg, 'contrast_est_');
            write_beta_slice(slicei, V(1), contrastt, emptyimg, 'contrast_t_');
            write_beta_slice(slicei, V(1), contrastste, emptyimg, 'contrast_ste_');
            write_beta_slice(slicei, V(1), contrastp, emptyimg, 'contrast_p_');
        else
            for mycon = 1:ncontrasts
                write_beta_slice(slicei, V(1), squeeze(contrastt(:,:,mycon)), emptyimg, [contrastnames{mycon} '_t'],1);
                write_beta_slice(slicei, V(1), squeeze(contrastste(:,:,mycon)), emptyimg, [contrastnames{mycon} '_ste'],1);
                write_beta_slice(slicei, V(1), squeeze(contrastp(:,:,mycon)), emptyimg, [contrastnames{mycon} '_p'],1);
            end
        end

    end

    if exist('ratings','var')
        write_beta_slice(slicei, V(1), ratingsest, emptyimg, 'ratings_est_');
        write_beta_slice(slicei, V(1), ratingst, emptyimg, 'ratings_t_');
        write_beta_slice(slicei, V(1), ratingsste, emptyimg, 'ratings_ste_');
        write_beta_slice(slicei, V(1), ratingsp, emptyimg, 'ratings_p_');

        write_beta_slice(slicei, V(1), intcptest, emptyimg, 'avg_trial_est_');
        write_beta_slice(slicei, V(1), intcptt, emptyimg, 'avg_trial_t_');
        write_beta_slice(slicei, V(1), intcptste, emptyimg, 'avg_trial_ste_');
        write_beta_slice(slicei, V(1), intcptp, emptyimg, 'avg_trial_p_');
    end

    fprintf(1, 'Done in %3.0f s\n', etime(clock, et))
end




function betas = fir_slice(fsl,i,j,PSDXS,F)

    % needed info
    global slicei
    global dims
    global cols
    global V

    % graphics
    global dographics
    global f1
    global f2

    % processing options
    global doysmooth
    global dosmooth
    global dohtw
    global dowritebetas


    % -------------------------------------------------------------------
    % * smooth timeseries data before fitting
    % -------------------------------------------------------------------
    if doysmooth

        fprintf(1, '\tSmoothing timeseries for %3.0f voxels > ', length(i))

        et = clock;
        for k = 1:length(i)

            fsl(i(k), j(k), :) = F * squeeze(fsl(i(k), j(k), :));
            if k == 1000, fprintf(1, '%3.0f s per 1000 vox.', etime(clock, et)), end

        end

        if dographics
            if ishandle(f2), figure(f2);subplot(2, 1, 1); else f2 = tor_fig(2, 1);  end
            plot(scale(squeeze(fsl(round(dims(1)./2), round(dims(2)./2), :))), 'b', 'LineWidth', 2),
            legend({'Raw' 'Filtered' 'Smoothed'}); title('Standardized data from center voxel')
        end

    end  % if doysmooth

    % -------------------------------------------------------------------
    % * fit the model to each nonzero voxel
    % -------------------------------------------------------------------
    betas = NaN * zeros([dims(1:2) cols]);
    fprintf(1, '\n\tFitting model matrix > ')

    et = clock;
    for k = 1:length(i)
        % b = pinv(rfS * X) * (rfS * y) do PSDX * rfS first
        % equivalent to: betas(i(k), j(k), :) = PSDX * (rfS * squeeze(sl(i(k), j(k), :)));

        betas(i(k), j(k), :) = PSDXS * squeeze(fsl(i(k), j(k), :));   %for NaN replacement: PSDXS(:, ~isnan(y)) * y(~isnan(y));

        if any(isnan(betas(i(k), j(k), :))), disp('Some betas are NaN for non-empty voxels.  This should not happen.'); keyboard, end

        if k == 1000, fprintf(1, '%3.0f s per 1000 vox.', etime(clock, et)), end
    end

    if dographics
        figure(f1); subplot 224; imagesc(betas(:, :, end)); colormap(gray); title('Last beta (intercept)'), colorbar('horiz'), drawnow;
        scn_export_papersetup(500);
        saveas(gcf, ['slices_and_trimming_' num2str(slicei)], 'png')
        figure(f2); subplot(2, 1, 2); plot(squeeze(betas(round(dims(1)./2), round(dims(2)./2), :)), 'k', 'LineWidth', 2)
    end
end




function [fbetas,bname] = smooth_beta_series(betas,numframes,smoothlen,TR,dims)
    global f1
    global f2
    global dographics

    global dowritebetas

    global slicei
    global V


    bname = [];
    % -------------------------------------------------------------------
    % * smooth beta series within condition to regularize
    % -------------------------------------------------------------------
    st = cumsum([1 numframes]);
    en = st(2:end) - 1;         % ending values
    st = st(1:end-1);           % starting values

    flen = max(1, round(smoothlen/TR));    % smoothing filter length: 6 s to zero weight in exponential LP smooth

    fprintf(1, '\n\tRegularizing betas > ')
    et = clock;
    fbetas = NaN*zeros(size(betas));    %initialize to NaN to avoid dim problems

    for bet = 1:length(st)      % for each beta series

        % get LP filter
        [x, F] = smooth_timeseries(squeeze(betas(round(dims(1)./2), round(dims(1)./2), st(bet):en(bet))), flen);

        % loop through voxels, save in fbetas
        for k = 1:length(i)
            fbetas(i(k), j(k), st(bet):en(bet)) = F * squeeze(betas(i(k), j(k), st(bet):en(bet)));   %for NaN replacement: PSDXS(:, ~isnan(y)) * y(~isnan(y));
            if k == 1000 & bet == 1, fprintf(1, '%3.0f s per 1000 vox.', etime(clock, et)), end

            if any(isnan(fbetas(i(k), j(k), st(bet):en(bet)))), disp('Smoothing has introduced NaN beta values.  This should not happen.'); keyboard, end
        end
    end

    if dographics,
        if ishandle(f2), figure(f2);subplot(2, 1, 1); else f2 = tor_fig(2, 1);  end
        figure(f2); subplot(2, 1, 2); plot(squeeze(fbetas(round(dims(1)./2), round(dims(2)./2), :)), 'r', 'LineWidth', 2)
        legend({'Beta series' 'Smoothed'})
        scn_export_papersetup(500);
        saveas(gcf, ['timeseries_center_voxel_' num2str(slicei)], 'png')
    end

    fprintf(1, 'Done in %3.0f s\n', etime(clock, et))

    % -------------------------------------------------------------------
    % * write smoothed betas
    % -------------------------------------------------------------------
    if dowritebetas
        fprintf(1, '\n\tWriting beta plane > ')
        et = clock;
        emptyimg = zeros(dims);     % in case we need to create a new volume

        bname = write_beta_slice(slicei, V(1), fbetas, emptyimg);

        fprintf(1, '%s. Done in %3.0f s\n', bname, etime(clock, et))
    end
end






function get_htw(fbetas,TR,numframes,st,en,emptyimg)
    global dographics
    global dosmooth
    global doysmooth
    global dohtw
    global dofilt
    global writefilt
    global dowritebetas
    global f1
    global f2
    global f3

    global dims
    global cols

    % -------------------------------------------------------------------
    % * get height, width, delay parameters
    % -------------------------------------------------------------------
    hconst = round(30 ./ TR);   % constrain peak height to within first 30 s
    hconst = min([repmat(hconst, 1, length(numframes)); (numframes)]);  %min([hconst min(numframes)]);  % min or we get an error
    fprintf(1, '\n\tExtracting height, delay, width using first %3.0f betas > ', hconst)

    sz = size(fbetas); sz = [sz(1:2) length(st)];
    h = NaN*zeros(sz);    %initialize to NaN to avoid dim problems
    t = NaN*zeros(sz);    %initialize to NaN to avoid dim problems
    w = NaN*zeros(sz);    %initialize to NaN to avoid dim problems

    warning off
    et = clock;
    for bet = 1:length(st)      % for each beta series

        % loop through voxels, save in fbetas
        for k = 1:length(i)

            % use version 2, which parameterizes from first beta value
            [h(i(k), j(k), bet), t(i(k), j(k), bet), w(i(k), j(k), bet)] = fir2htw2(squeeze(fbetas(i(k), j(k), st(bet):en(bet))), hconst(bet), 0, {'k'});

            if k == 1000 & bet == 1, fprintf(1, '%3.0f s per 1000 vox.', etime(clock, et)), end

        end
    end
    warning on

    if dographics
        if ishandle(f3)
            figure(f3);
        else
            f3 = figure;
            set(gcf, 'Color', 'w');
            tmp = get(gcf, 'Position');
            tmp(3) = tmp(3)+.5*tmp(3);
            set(gcf, 'Position', tmp)
        end
        subplot(1, 4, 1);
        imagesc(h(:, :, 1)); colormap(hot); title('Height', 'FontSize', 24), colorbar('horiz'), axis image; drawnow;
        subplot(1, 4, 2); imagesc(t(:, :, 1)); title('Time to peak', 'FontSize', 24), colorbar('horiz'), axis image, drawnow;
        subplot(1, 4, 3); imagesc(w(:, :, 1)); title('Width', 'FontSize', 24), colorbar('horiz'), axis image, drawnow;
        drawnow

        k = find(h == max(h(:)));    % max height
        [ii, jj, kk] = ind2sub(size(h), k);

        subplot(1, 4, 4), cla
        y = squeeze(fbetas(ii, jj, st(kk):en(kk)));
        plot(y, 'ro-', 'LineWidth', 2)
        [h2, t2, w2] = fir2htw2(y, hconst(kk), 1, {'k'});

        scn_export_papersetup(500);
        saveas(gcf, ['height_time_wid_' num2str(slicei)], 'png')
    end

    % write images for all beta series
    write_beta_slice(slicei, V(1), h, emptyimg, 'height_event_');

    write_beta_slice(slicei, V(1), t, emptyimg, 'delay_event_');

    write_beta_slice(slicei, V(1), w, emptyimg, 'width_event_');


    fprintf(1, '\tDone in %3.0f s\n', etime(clock, et))
end
