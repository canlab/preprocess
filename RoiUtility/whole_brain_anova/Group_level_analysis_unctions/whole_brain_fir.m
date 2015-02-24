function EXPT = whole_brain_fir(EXPT,varargin)
% EXPT = whole_brain_fir(EXPT,[type],[startat subject #],[start slice])
%
% start slice reverts to 1 after the first subject
%
% Examples:
% EXPT = whole_brain_fir(EXPT,'full',29);
%
%
% betas: 
% EXPT = getfunctnames2(EXPT,'dx*img','FIR.dxbetas');
% EXPT = whole_brain_fir(EXPT,'htw');
%
%
% The analysis sequence:
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

% start at subject and slice (2nd subject starts at slice 1 auto)
if length(varargin) > 0, type = varargin{1};, else, type = 'full';, end
if length(varargin) > 1
    if(length(varargin{2}) == 1)
        subjidxs = varargin{2}:length(EXPT.subjects);
    else
        subjidxs = varargin{2};
    end
else
    subjidxs = 1:length(EXPT.subjects);
end
if length(varargin) > 2, firstslice = varargin{3};, else, firstslice = 1;, end

imgnames = EXPT.FILES.im_files;       % for full model
if strcmp(type,'htw'),
    imgnames = EXPT.FIR.dxbetas;    % for htw ONLY
end


if isfield(EXPT,'mask')
    mask = spm_read_vols(spm_vol(EXPT.mask));
elseif isfield(EXPT.FIR,'mask')
    mask = spm_read_vols(spm_vol(EXPT.FIR.mask));
else
    mask = [];
end

fprintf(1,'Processing subject indices: %s\n', num2str(subjidxs));

for i = subjidxs
    
    %spmname = fullfile(EXPT.studydir{i},EXPT.subjects{i},'SPM.mat');
    %warning off; DX = spm2dx(spmname); warning on

    
    mkdir(EXPT.subjects{i})
    
    cd(EXPT.subjects{i})
    
    if exist(EXPT.FILES.im_files{1}(1,:)) == 2
        % only if we can find images
        
        %Pw = whole_brain_filter(P,DX,TR,HP,nruns,numframes,[doplot],[mask],[smoothlen])
        Pw = whole_brain_filter( ...
        imgnames{i}, ...
        EXPT.FIR.model{i}, ...      %(:,1:end-1), ...
        EXPT.TR, ...
        EXPT.FIR.HP, ...
        EXPT.FIR.nruns, ...
        EXPT.FIR.numframes, ...
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
