function canlab_combine_structurals_luka(meanstruct, structs)
% canlab_combine_structurals_luka(meanstruct, structs)
%
% DESCRIPTION
%   This function takes multiple structural images, aligns them using SPM,
%   and creates an average image.
%
% ARGUMENTS
%   meanstruct
%     filename for created image
%   structs
%     a cell array of structural image filenames
%
% EXAMPLE
%   structs = filenames('*/Structural/SPGR/mprage_rms.nii','absolute');
%   canlab_combine_structurals_luka('mprage_rms_study_average', structs);
%

% add .nii image file extension if necessary
if ~regexp(meanstruct,'\.nii$'), meanstruct = [meanstruct '.nii']; end

% create temporary working directory
wd = fullfile(pwd, ['tempdir_' num2str(sprintf('%05d',randi(10000)))]);
mkdir(wd);

% combine structurals with SPM realign job
fprintf('... STARTED structural realignment\n')
nimgs = numel(structs);
matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = cell(nimgs,1);
% make temporary copies of target files and add to SPM realign job
for n = 1:nimgs
    [d f e] = fileparts(structs{n});
    wf{n} = fullfile(wd, [f '_' num2str(n) e]);
    copyfile(structs{n},wf{n});
    matlabbatch{1}.spm.spatial.realign.estwrite.data{1}{n} = wf{n};
end
% set options for SPM realign
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {''};
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
% run SPM realign estimate+write
spm('defaults', 'FMRI');
spm_jobman('serial', matlabbatch, '', cell(0,1));
fprintf('... FINISHED structural realignment\n')

% save created file
[d f e] = fileparts(matlabbatch{1}.spm.spatial.realign.estwrite.data{1}{1});
movefile(fullfile(d, ['mean' f e]), meanstruct);

% delete temporary working directory
rmdir(wd,'s');

fprintf('... WROTE %s\n', meanstruct)

end
