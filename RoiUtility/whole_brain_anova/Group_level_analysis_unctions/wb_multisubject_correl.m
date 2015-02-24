function EXPT = wb_multisubject_correl(EXPT,varargin)
% EXPT = wb_multisubject_correl(EXPT,[type],[subject indices],[start slice])
%
% Critical input fields:
% EXPT.FILES.im_files       all image names, in cell array, one cell per subject
% EXPT.seeds                seeds to correlate, separately, one cell per subject
% EXPT.TR                   TR in sec
% EXPT.FIR.nruns            vector of images per run, redundant with numframes; vector if same for all subjects, cell for each subject if not
% EXPT.FIR.numframes        vector of images per run, vector if same for all subjects, cell for each subject if not
% EXPT.HP.len               high-pass filter length (sec); if not included, no filtering
% EXPT.FIR.smoothlen        % not used... enter 0
%
% start slice reverts to 1 after the first subject
%
% Examples:
% EXPT = wb_multisubject_correl(EXPT,'full',29:length(EXPT.subjects));
%
%
% betas:
% EXPT = getfunctnames2(EXPT,'dx*img','FIR.dxbetas');
% EXPT = wb_multisubject_correl(EXPT,'htw');
%
%
% The analysis sequence:
% Run wb_multisubject_correl, which mainly just calls whole_brain_correl
% whole_brain_correl gives you:
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

% start at subject and slice (2nd subject starts at slice 1 auto)
if length(varargin) > 0, type = varargin{1};, else, type = 'full';, end
if length(varargin) > 1, subjidxs = varargin{2};, else, subjidxs = 1:length(EXPT.subjects);, end
if length(varargin) > 2, firstslice = varargin{3};, else, firstslice = 1;, end

imgnames = EXPT.FILES.im_files;       % for full model
if strcmp(type,'htw'),
    imgnames = EXPT.FIR.dxbetas;    % for htw ONLY
end


if isfield(EXPT.FIR,'mask'), 
    mask = spm_read_vols(spm_vol(EXPT.FIR.mask));
else
    mask = [];
end

if ~isfield(EXPT, 'HP') || isempty(EXPT.HP) || ~isfield(EXPT.HP, 'len') || isempty(EXPT.HP.len)
    EXPT.HP.len = Inf;
end

if ~isfield(EXPT, 'FIR') || ~isfield(EXPT.FIR, 'smoothlen') || isempty(EXPT.FIR.smoothlen)
    EXPT.FIR.smoothlen = 0;
end

for i = subjidxs
    
    %spmname = fullfile(EXPT.studydir{i},EXPT.subjects{i},'SPM.mat');
    %warning off; DX = spm2dx(spmname); warning on
    
    if iscell(EXPT.FIR.nruns)
        myruns = EXPT.FIR.nruns{i};
    else
        myruns = EXPT.FIR.nruns;
    end
    
    if iscell(EXPT.FIR.numframes)
        myframes = EXPT.FIR.numframes{i};
    else
        myframes = EXPT.FIR.numframes;
    end
    
        
    mkdir(EXPT.subjects{i})
    
    cd(EXPT.subjects{i})
    
    % get sample image name to see if it exists
    % take off ,### from end if spm-style volume reference is used
    sampleimage = deblank(imgnames{i}(1,:));
    whcomma = find(sampleimage == ',');
    if ~isempty(whcomma), sampleimage = sampleimage(1:whcomma-1); end
        
    if exist(sampleimage, 'file')
        % only if we can find images
        
        %Pw = whole_brain_correl(P,DX,TR,HP,nruns,numframes,[doplot],[mask],[smoothlen])
        Pw = whole_brain_correl( ...
        imgnames{i}, ...
        EXPT.seeds{i}, ...
        [], ...      %(:,1:end-1), ...
        EXPT.TR, ...
        EXPT.HP.len, ...
        myruns, ...
        myframes, ...
        0, ...                      % do graphics
        mask, ...
        EXPT.FIR.smoothlen, ...
        firstslice, ...
        type ...
        );
    
        
    else
        disp(['Missing images for ' EXPT.subjects{i}])
    end
            
    firstslice = 1;     % re-set first slice to 1
    
    cd ..
    
end

return
