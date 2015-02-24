function EXPT = wb_multisubject_parammod(EXPT,varargin)
% EXPT = wb_multisubject_parammod(EXPT,[type],[startat subject #],[start slice])
%
% start slice reverts to 1 after the first subject
%
% Examples:
% EXPT = wb_multisubject_parammod(EXPT,'full',29);
%
%
%
% The analysis sequence:
% Run whole_brain_fir, which mainly just calls whole_brain_parammod
% whole_brain_filter gives you:
%       1) trimmed, filtered images, saved in f*img
%           these can be passed to SPM2, e.g., if desired
%       2) extraction of beta images for FIR model, saved in dx_beta*img
%           in individual subject folders; also smooths betas if desired
%       3) calculation of height, delay, and width for each FIR extracted
%           saved on images in individual subject folders.
%
% After running this, you'll need to collect image names for use in random
% effects analyses.  Contrasts are created in this function.
%
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
if length(varargin) > 1, startat = varargin{2};, else, startat = 1;, end
if length(varargin) > 2, firstslice = varargin{3};, else, firstslice = 1;, end
if length(varargin) > 3, subjectlist = varargin{4};, else, subjectlist = startat:length(EXPT.subjects);, end
try
    imgnames = EXPT.FILES.im_files;       % for full model
    disp('Found image names in: EXPT.FILES.im_files');
catch
    imgnames = EXPT.im_files;       % for full model
    disp('Found image names in: EXPT.im_files');
end
    
if strcmp(type,'htw'),
    imgnames = EXPT.RT.dxbetas;    % for htw ONLY
end

% remove empty trailing images
TMP.imgnames = imgnames;
TMP = remove_empty_image_names(TMP,'imgnames');
imgnames = TMP.imgnames; clear TMP;

if isfield(EXPT.FIR,'mask'), 
    mask = spm_read_vols(spm_vol(EXPT.RT.mask));
else
    mask = [];
end


for i = subjectlist  % [6 7 8 9 27] %length(EXPT.subjects)
    
    %spmname = fullfile(EXPT.studydir{i},EXPT.subjects{i},'SPM.mat');
    %warning off; DX = spm2dx(spmname); warning on

    % adjust number of runs if necesary
    % will return an error (likely) if dropping a session won't work
    % or if the n imgs per session varies and a middle session is missing.
    nruns = EXPT.RT.nruns;
    mydif = cumsum(nruns) - size(imgnames{i},1);
    if mydif(end) ~= 0,
        newend = find(mydif == 0);
        nruns = nruns(1:newend);
    end

    mkdir(EXPT.subjects{i})
    
    cd(EXPT.subjects{i})
    
    if exist(EXPT.im_files{1}(1,:)) == 2
        % only if we can find images
        
        %Pw = whole_brain_filter(P,DX,TR,HP,nruns,numframes,[doplot],[mask],[smoothlen])
        Pw = whole_brain_parammod( ...
        imgnames{i}, ...
        EXPT.RT.model{i}, ...      %(:,1:end-1), ...
        EXPT.TR, ...
        EXPT.RT.HP, ...
        nruns, ...                  % vector of images in each run
        EXPT.RT.whpredictors, ...   % which predictors to save betas of
        EXPT.RT.contrasts, ...      % contrasts across betas(whpredictors)
        EXPT.RT.connames, ...       % names of contrasts, no spaces (for filenames)
        0, ...                      % do graphics
        mask, ...
        EXPT.RT.smoothlen, ...
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
