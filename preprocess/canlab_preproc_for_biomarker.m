function dat = canlab_preproc_for_biomarker(dat, TR, varargin)
% dat = canlab_preproc_for_biomarker(dat, TR, [optional inputs])
%
% Take images (typically swravols or wravols) and perform additional
% preprocessing steps found to be helpful in pain biomarker analyses.
%
% Saves a single 4-D output image with the whole dataset.
% (This can be omitted.)
%
% Processing steps include:
% 1. 'session_outliers': Identify session-wise (run-wise) outliers with significant
% based on mahalanobis distance with FDR-corrected P-values in chi-square test.
% Impute session grand mean outliers.
%
% 2. 'windsorize': Windsorize entire data multisession matrix to 3 STD,
% replacing values > +-3 standard deviations with the 3 SD value.
%
% 3. 'percentchange': Rescale each voxel to percent signal change, with a mean of 100,
% based on a locally smoothed (16 mm) baseline signal estimate across the
% entire multisession dataset. Note: SPM does it's own rescaling, I believe.
%
% 4. 'hpfilter': High-pass filter with 180 sec cutoff by default (a
% different cutoff can be entered as an optional input argument).
%
% Optional inputs
% ------------------------------------------------------------------------
% 'hpfilter', followed by HP filter cutoff in sec
%
% 'images_per_session', followed by vector of num images in each scanning
%   run, e.g., [150 150 150] for 3 runs of 150 images each.
%
% Examples:
% =========================================================================
% imgs = filenames('dpsp043/r*/swra*.img', 'char', 'absolute'); imgs
% outpath = fileparts(imgs(1, :)); outpath = fileparts(outpath); % get subject dir
% [~, outputname] = fileparts(outpath); % get subject name
% outputname = ['swra_canlabpreproc_' outputname '.img']; % define output image name
% dat = canlab_preproc_for_biomarker(imgs, 2, 'images_per_session', repmat(180, 1, 8));
%
% Example script (define inputs/outputs and run):
% =========================================================================
% subjname = 'dpsp043';
%
% % Define image names
% % -------------------------------------------------------------------------
% rundirs = dir(fullfile(subjname, 'r*'));
% imgs = cell(1, length(rundirs));
% for i = 1:length(rundirs)
%     imgs{i} = filenames(fullfile(subjname, rundirs(i).name, 'swra*.img'), 'char', 'absolute');
%     imgs_per_session(1, i) = size(imgs{i}, 1);
% end
%
% imgs = char(imgs{:});
%
% % Define output image name (a single 4-D image)
% % -------------------------------------------------------------------------
% outpath = fileparts(imgs(1, :)); outpath = fileparts(outpath); % get subject dir
%
% outputname = ['swra_canlabpreproc_' subjname '.img']; % define output image name
%
% % Run
% % -----------------------------------------------------------------------
% dat = canlab_preproc_for_biomarker(imgs, 2, 'images_per_session', repmat(180, 1, 8), 'path', outpath, 'outputname', outputname);

% defaults
% ---------------------------------------------------------------
outputname = 'canlab_bmrk_preproc.img';
images_per_session = [];
outputpath = [];
dowrite = 1;

% Process optional inputs
% ---------------------------------------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            % functional commands
            case 'images_per_session', images_per_session = varargin{i+1};
            case {'HP', 'hp' 'hpfilter'}, y = varargin{i+1};
                
            case 'outputname', outputname = varargin{i + 1};
            case {'path', 'outputpath'}, outputpath = varargin{i + 1};
                
            case 'nowrite', dowrite = 0;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if ischar(dat)
    % Load into object
    fprintf('Loading images entered in string matrix.\n');
    dat = fmri_data(dat);
    
    if isempty(images_per_session)
        error('Warning! Images per session is empty. For a multi-session design, enter as optional input.');
    else
        dat.images_per_session = images_per_session;
    end
    
elseif isa(dat, 'fmri_data')
    % do nothing - but check imgs per sess
    
    if isempty(images_per_session)
        disp('Warning! Images per session is empty. Define this before running.');
        images_per_session = size(dat.dat, 2);
    end
    
else
    error('Enter dat in char array of image names or fmri_data object.');
end

if isempty(outputpath)
    outputpath = fileparts(dat.fullpath(1, :));
end


% Do the work
% ---------------------------------------------------------------
fprintf('Starting preprocessing.\n');

dat = preprocess(dat, 'session_outliers');

dat = preprocess(dat, 'windsorize');

dat = rescale(dat, 'percentchange');

dat = preprocess(dat, 'hpfilter', 180, TR);

fprintf('Completed preprocessing. Summary:\n');
history(dat)

% write output images

dat.fullpath = fullfile(outputpath, outputname);

if dowrite
    write(dat);
end

end % function