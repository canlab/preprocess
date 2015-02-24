% function flip_X_pixdim(imgs, ['flip_negative_only', 0|1], ['verbose', 0|1])
%   Takes a list of .img files (as a cellstr or char array) and flips the x dimension in the affine matrix
%
%   Parameter Options:
%      flip_negative_only - flip only images with negative X pixdims (default: 0)
%      verbose - print status (default: 0)
%
%   Use this when you have weird cases where the X dimension/voxel size is negative and you want to make it positive, or vice versa.
%   E.g., turning:
%
%   -3.5000    0.0004   -0.0003  111.9878
%    0.0004    3.5000   -0.0016 -111.9920
%   -0.0003    0.0013    4.5000  -67.4930
%         0         0         0    1.0000
%
%   into:
%
%    3.5000    0.0004   -0.0003 -111.9878
%   -0.0004    3.5000   -0.0016 -111.9920
%    0.0003    0.0013    4.5000  -67.4930
%         0         0         0    1.0000

function flip_X_pixdim(imgs, varargin)
flip_negative_only = 0;
verbose = 0;

for i=1:length(varargin)
    if(ischar(varargin{i}))
        switch(varargin{i})
            case 'flip_negative_only'
                flip_negative_only = varargin{i+1};
            case 'verbose'
                verbose = varargin{i+1};
        end
    end
end

imgs = cellstr(imgs);
if(verbose), fprintf('Flipping %d files: %3d%%', length(imgs), 0); end

for i = 1:length(imgs)
    skipped = invert_file(imgs{i}, flip_negative_only, verbose);
    
    if(verbose)
        if(skipped)
            fprintf('Flipping %d files: %3d%%', length(imgs), round(i * 100.0 / length(imgs)));
        else
            fprintf('\b\b\b\b%3d%%', round(i * 100.0 / length(imgs)));
        end
    end
end

fprintf('\n');
end


function skipped = invert_file(img, flip_negative_only, verbose)
skipped = invert_file_SPM(img, flip_negative_only, verbose);
end

function skipped = invert_file_SPM(img, flip_negative_only, verbose)
skipped = 0;
[pathstr, imgbase] = fileparts(img);
mat_file = fullfile(pathstr, [imgbase '.mat']);
if(~exist(mat_file, 'file'))
    Vimg = spm_vol(img);
    M = Vimg.mat;
    mat = M; %#ok
    save(mat_file, 'M', 'mat');
end

load(mat_file);
if ~exist('M', 'var'), M = mat; end

if(flip_negative_only && M(1,1) >= 0)  %#ok
    if(verbose), fprintf('\nSkipping file %s: already positive\n', img); end
    skipped = 1;
else
    % tor updated for spm8...M was not a variable, mat is. use with
    % caution.
    %if exist('M', 'var')
    M = invert_affine_matrix(M);
    mat = M; %#ok
    %else
    %    mat = invert_affine_matrix(mat);
    %end
    save(mat_file, 'M', 'mat');
end
end

% Seriously broken, since FSL doesn't alter Analyze hdrs, only Nifti hdrs
% function invert_file_FSL(img, display_commands)
%     global FSLDIR;
%     scn_setup();
%
%     Vimg = spm_vol(img);
%     flip_command = sprintf('export FSLDIR=%s && . %s/etc/fslconf/fsl.sh && FSLOUTPUTTYPE=ANALYZE && %s/bin/avwchpixdim %s %f %f %f', ...
%         FSLDIR, FSLDIR, FSLDIR, Vimg(1).fname, -1 * Vimg(1).mat(1,1), Vimg(1).mat(2,2), Vimg(1).mat(3,3));
%     if(display_commands), disp(flip_command); end
%     unix(flip_command);
% end

function M = invert_affine_matrix(M)
for i = 1:size(M, 3)
    % tor updated for 4-D files, 8/2012
    params = spm_imatrix(M(:, :, i));
    params(7) = -params(7);
    params(1) = -params(1);
    M(:, :, i) = spm_matrix(params);
end
end