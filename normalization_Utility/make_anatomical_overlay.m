function outname = make_anatomical_overlay(my_mean_anat, varargin)
    % Segments a mean anatomical image and makes a brain mask from it.
    % outname = function make_anatomical_overlay(my_mean_anat, [option:enter 1 to mask with canonical brain as well])
    %
    % tor wager, aug 2009
    %
    % Example:
    % outname = make_anatomical_overlay(my_mean_anat, 1);
    
    smoothval = 2;
    
load('/Users/tor/Documents/matlab_code/normalization_Utility/warp_mean_anatomical_job.mat')

jobs{1}.spatial{1}.preproc.data{1} = my_mean_anat;

jobs{1}.spatial{1}.preproc.opts.tpm{1} = which('grey.nii');
jobs{1}.spatial{1}.preproc.opts.tpm{2} = which('white.nii');
jobs{1}.spatial{1}.preproc.opts.tpm{3} = which('csf.nii');


% threshold each partition (no, don't)
% img = 'c1mean_wT1.nii';
% mask_image(img, img, img, 'minmask', .05);
% img = 'c2mean_wT1.nii';
% mask_image(img, img, img, 'minmask', .05);

spm_jobman('run', jobs);

[dd, ff, ee] = fileparts(my_mean_anat);

%
gray = fullfile(dd, ['c1' ff ee]);
white = fullfile(dd, ['c2' ff ee]);
smoothed = fullfile(dd, ['sc1' ff ee]);

%
% smooth a *tiny* bit
spm_smooth(gray, smoothed, smoothval);

% maybe add this? 'erode' a bit
spm_imcalc_ui(smoothed, smoothed,'i1 > .7')
spm_smooth(smoothed, smoothed, round(smoothval./2));
spm_imcalc_ui(smoothed, smoothed,'i1 > 0')

% create brain mask
maskname = 'brainmask.img';
mask_union([], maskname, smoothed, white);

% mask the actual image
outname = fullfile(dd, 'mean_wT1_brain.img');

mask_image(my_mean_anat, 'brainmask.img', outname);

spm_image('init', outname)

if nargin > 1
    disp('Masking with brainmask.nii')
    mask = which('brainmask.nii');
    scn_map_image(mask, outname, 'write', 'canonical_brainmask.img');
    mask_image(outname, 'canonical_brainmask.img', outname, 'minmask', .2);
    spm_image('init', outname);
end

end